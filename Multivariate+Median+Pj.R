################################################################
# In same order as our plots
#-----------------------------------------------
library("MASS")
library("Gmedian") 
library("OjaNP")
library("RColorBrewer")
library("robust")
library("SpatialNP")
library("lattice")
#-----------------------------------------------
# self implemented function
ellipse <- function(mu, sigma, c2, output=FALSE) {
  # mu: vector
  # sigma: matrix
  # c2 <- qchisq(p,df)
  # output: print the ellipse if TRUE
  # return dot series of ellipse
  
  es <- eigen(sigma)
  c1 <- sqrt(c2)
  e1 <- es$vec%*%diag(sqrt(es$val))
  theta <- seq(0,2*pi,len=1000)
  v1 <- cbind(c1*cos(theta),c1*sin(theta))
  pts=t(mu+(e1%*%t(v1)))
  # graph
  if(output==TRUE){plot(pts, type="l",xlab="x", ylab="y")}
  return(pts)
}
#-----------------------------------------------
# marginal median
mmed <- function(data){
  # data[,i]: all data from ith dimension
  # di: number of dimensions
  di <- dim(sample)[2]
  result <- rep(0,di)
  for (i in 1:di){ result[i] <- median(data[,i]) }
  return(result)
}
# marginal sign
marsign <- function(mm, sample){
  n <- length(mm)
  sign <- rep(0,n)
  for (i in 1:n){
    if (sample[i] < mm[i]){sign[i]<- (-1)}
    if (sample[i] == mm[i]){sign[i]<- 0}
    if (sample[i] > mm[i]){sign[i]<- 1}    
  }
  return(sign)
}
# marginal sign covariance matrix
MSCM <- function(mm, sample){
  # mm: marginal median
  # sample : original sample
  n <- length(mm)
  scm <- matrix(0, nrow=n, ncol=n)
  for (i in 1:dim(sample)[1]) {
    tmp <- marsign(mm,sample[i,])
    scm <- scm + tmp %*% t(tmp)
  }
  return(scm/(dim(sample)[1]))
}
# calculate covariance matrix based on marginal sign
MSCM_covariance <- function(sample){
  m_med <- mmed(sample) # center: spatial median
  scm <- MSCM(m_med, sample) # spatial sign covariance matrix
  U <- eigen(scm)$vector # step 1
  # step 2
  trans <- cbind(rep(0,10),rep(0,10))
  for (i in 1:10){
    trans[i,] <- t(U%*%sample[i,])
  }
  # marginal variance
  lambda <- matrix(c(mad(trans[,1]),0,0,mad(trans[,2])),2,2)
  return(t(U)%*%lambda%*%U)
}
#-----------------------------------------------
# spatial sign covariance matrix
# calculate spatial median: Gmedian(sample)
SSCM <- function(mm, sample){
  # mm: spatial median
  # sample: original sample
  scm <- matrix(c(0,0,0,0),2,2)
  for (i in 1:dim(sample)[1]){
    tmp <- sample[i,]-mm # direction
    long <-sqrt((tmp[1])^2+(tmp[2])^2) # length
    tmp <- tmp/long # unit vector
    scmi <- t(tmp)%*%tmp
    scm <- scm + scmi
  }
  return(scm/(dim(sample)[1]))
}
# calculate covariance matrix based on spatial/oja sign
SCM_covariance <- function(sample,type){
  if (type=="OJA") {scm <- ojaSCM(sample)}
  if (type=="Spatial") {scm <- SCov(sample)}
  U <- eigen(scm)$vector # step 1
  # step 2
  trans <- matrix(0,ncol=dim(sample)[1], nrow=dim(sample)[2])
  for (i in 1:dim(sample)[2]){ # loop every eigenvector
    ui <- U[,i]
    trans[i,] <- ui%*%t(sample)
  }
  # marginal variance
  lambda <- rep(0, dim(sample)[2])
  for (j in 1:dim(sample)[2]) {lambda[j] <- mad(trans[j,])^2}
  # diagonal matrix
  H <- matrix(0, ncol=dim(sample)[2], nrow=dim(sample)[2])
  for (k in 1:dim(sample)[2]) {H[k,k] <- lambda[k]}
  
  return(U%*%H%*%t(U)) # scatter matrix
}
#-----------------------------------------------
# OJA
OSCM <- function(sample){
  ojaSCM(sample)
}
#-----------------------------------------------
# draw sign vector
drawsign <- function(mm, sample, m="marginal", color){
  if (m=="oja") {sign <- ojaSign(sample)}
  for (i in 1:dim(sample)[1]){
    samplei <- sample[i,]
    if (m=="spatial") {direction <- (samplei-mm)/sqrt((mm[1] - samplei[1])^2+(mm[2] - samplei[2])^2)}
    if (m=="marginal") {direction <- marsign(mm, samplei)/sqrt(2)}
    if (m=="oja") {direction <- sign[i,]}
    end <- direction + mm
    arrows(mm[1],mm[2],end[1],end[2],length = 0.1, angle = 15,col=color)
  }
}
################################################
# Figure 1: draw sign and sign covariance
set.seed(0)
sample <- mvrnorm(10, c(0,0), matrix(c(1,0,0,1),2,2))
c2 <- qchisq(0.8,2)

# marginal
mar <- mmed(sample)
color = "brown2"
png(filename="Marginal_Median_and_Sign.png")
plot(sample, xlab="x", ylab="y", asp=1,xlim=c(-2,3))
points(mar[1],mar[2],pch=20, col=color)
abline(h=mar[2],lty=6,col="grey")
abline(v=mar[1],lty=6,col="grey")
drawsign(mar, sample,"marginal",color)
sigma <- MSCM(mar, sample)
mu <- c(mar[1],mar[2])
ellipse_mar <- ellipse(mu, sigma, c2)
lines(ellipse_mar, col=color, lty=1)
dev.off()

# spatial
color = "darkorange"
l1 <- Gmedian(sample)
png(filename="Spatial_Median_and_Sign.png")
plot(sample, xlab="x", ylab="y", asp=1)
points(l1[1],l1[2],pch=20, col=color)
abline(h=l1[2],lty=6,col="grey")
abline(v=l1[1],lty=6,col="grey")
drawsign(l1, sample, "spatial", col=color)
sigma <- SSCM(l1,sample)
mu <- c(l1[1],l1[2])
ellipse_sp <- ellipse(mu, sigma, c2)
lines(ellipse_sp, col=color, lty=1)
dev.off()

# OJA
color = "dimgray"
oja <- ojaMedian(sample)
png(filename="OJA_Median_and_Sign.png")
plot(sample, xlab="x", ylab="y", asp=1)
points(oja[1],oja[2],pch=20, col=color)
abline(h=oja[2],lty=6,col="grey")
abline(v=oja[1],lty=6,col="grey")
drawsign(oja, sample, "oja", color)
mu <- c(oja[1],oja[2])
sigma <- OSCM(sample)
ellipse_oja <- ellipse(mu, sigma, c2)
lines(ellipse_oja, col=color, lty=1)
dev.off()
################################################
# Figure 2
# condition number (shape)
cond <- function(matrix){
  value <- eigen(matrix)$values
  return(max(value)/min(value))
}
# change dimension
set.seed(321)
SSCM_ratio <- matrix(0, nrow = 19, ncol = 19)
MSCM_ratio <- matrix(0, nrow = 19, ncol = 19)
for (i in 2:20){
  print(i)
  sigma <- diag(x=1,i,i)
  mu <- rep(0, i)
  sample <- rmvnorm(n = 10, mean = mu, sigma = sigma)
  x <- seq(-1000,1000,110)
  SSCM_cond <- cond(SCov(sample))
  m_med <- mmed(sample)
  MSCM_cond <- cond(MSCM(m_med,sample))
  SSCM_IF <- rep(0, length(x))
  MSCM_IF <- rep(0, length(x))
  for (j in 1:length(x)){
    tmp <- sample[10,]
    tmp[1] <- j
    samplei <- rbind(sample[1:9,], tmp)
    m_medi <- mmed(samplei)
    SSCM_IF[j] <- cond(SCov(samplei))/SSCM_cond
    MSCM_IF[j] <- cond(MSCM(m_medi,samplei))/MSCM_cond
  }
  SSCM_ratio[i-1,] <- SSCM_IF
  MSCM_ratio[i-1,] <- MSCM_IF
}

ratio <- colorRampPalette(brewer.pal(9, "Blues"))(100)
png(filename="SSCM_ratio.png",width = 400, height = 400, units = "px")
levelplot(SSCM_ratio, xlab = list(label="dimension",cex=1.3), ylab = list(label="outlier",cex=1.3),col.regions=ratio, row.values=seq(2,20,1), column.values = seq(-1000,1000,110),asp=1)
dev.off()
png(filename="MSCM_ratio.png",width = 400, height = 400, units = "px")
levelplot(MSCM_ratio, xlab = list(label="dimension",cex=1.3), ylab = list(label="outlier",cex=1.3),col.regions=ratio, row.values=seq(2,20,1), column.values = seq(-1000,1000,110),asp=1)
dev.off()

# generate a symmetric sample
set.seed(0)
mu <- c(0,0,0)
sigma <- matrix(c(1,0,0,0,1,0,0,0,1),3,3)
sample <- rmvnorm(n = 10, mean = mu, sigma = sigma)
x <- seq(-500,500)

OSCM_cond <- cond(ojaSCM(sample))
SSCM_cond <- cond(SCov(sample))
m_med <- mmed(sample)
MSCM_cond <- cond(MSCM(m_med,sample))
OSCM_IF <- rep(0, length(x))
SSCM_IF <- rep(0, length(x))
MSCM_IF <- rep(0, length(x))
for (i in 1:length(x)){
  print(i)
  tmp <- sample[10,]
  tmp[1] <- i
  samplei <- rbind(sample[1:9,], tmp)
  m_medi <- mmed(samplei)
  OSCM_IF[i] <- cond(ojaSCM(samplei))/OSCM_cond
  SSCM_IF[i] <- cond(SCov(samplei))/SSCM_cond
  MSCM_IF[i] <- cond(MSCM(m_medi,samplei))/MSCM_cond
}
png(filename="OJA_IF.png",width = 400, height = 400, units = "px")
plot(x,OSCM_IF, xlim=c(-500,500),type="l",xlab="Outlier",ylab="IF",cex.lab=1.3)
legend("topleft", legend=c("Oja median"),bty="n",lty=1,cex=2)
dev.off()
png(filename="mar_spa_IF.png",width = 400, height = 400, units = "px")
plot(x,SSCM_IF, xlim=c(-500,500),type="l",ylim=c(0.5,2),xlab="Outlier",ylab="IF",col="blue3",cex.lab=1.3)
lines(x, MSCM_IF,col="darkred")
legend("topright", legend=c("marginal median", "spatial median"),col=c("darkred","blue3"), bty="n", lty=1, ncol=1, cex=2)
dev.off()

################################################
# Figure 3: Affine
set.seed(0)
mu <- rbind(c(-40,-40), c(10,0), c(-20,20),c(30,30),c(20,-40))
Sigma <- rbind(matrix(c(40,0,0,70),2,2),matrix(c(10,0,0,10),2,2),matrix(c(50,0,0,10),2,2),matrix(c(5,3,2,8),2,2),matrix(c(10,0,0,8),2,2))
pop <- c(100,200,300,400,500)
sample <- cbind(rep(0, sum(pop)),rep(0, sum(pop)))
flag <- 1
for (i in 1:5){
  mui <- mu[i,]
  sig <- Sigma[(2*i-1):(2*i),]
  p <- pop[i]
  sample[flag:(flag+p-1),]<-mvrnorm(p, mui, sig)
  flag <- flag + p
  print(flag)
}

l1 <- Gmedian(sample)
mm <- mmed(sample)

theta <- seq(0,2*pi,len=1000)
t.spatial.rotate <- matrix(nrow = length(theta), ncol = 2)
t.marginal.rotate <- matrix(nrow = length(theta), ncol = 2)
for (i in 1:length(theta)){
  angle <- theta[i]
  A <- rbind(cbind(cos(angle),-sin(angle)),cbind(sin(angle),cos(angle)))
  t.spatial.rotate[i,] <- Gmedian(sample %*% A)
  t.marginal.rotate[i,] <- mmed(sample %*% A)
  # m <- test %*% A
  # t.oja.rotate[i,] <- ojaMedian(m, alg="exact")
}
color=c("blue3", "blueviolet")
png(filename="rotate.png",width = 400, height = 400, units = "px")
plot(sample,asp=1,col="grey",ylab="y",xlab="x")
lines(t.spatial.rotate,col=color[1])
lines(t.marginal.rotate, col=color[2])
dev.off()

library("pracma")
set.seed(0)
A.rand <- rand(2,2)
for (i in 1:length(theta)){
  angle <- theta[i]
  A <- rbind(cbind(cos(angle),-sin(angle)),cbind(sin(angle),cos(angle)))
  t.spatial.rotate[i,] <- Gmedian(sample %*% A %*% A.rand)
}
test <- sample %*% A.rand
png(filename="affine.png",width = 400, height = 400, units = "px")
plot(test,asp=1,col="grey",ylab="y",xlab="x")
lines(t.spatial.rotate, col=color[1])
legend("topleft", legend=c("Spatial Median"), lty=c(1),col=color[1], bty="n", ncol=1, cex=2)
legend("bottomright", legend=c("Marginal Median"), lty=c(1),col=color[2], bty="n", ncol=1, cex=2)
dev.off()
################################################
# Weiszfeld's algorithm in Fig.4
set.seed(0)
test <- mvrnorm(500, c(1,1), matrix(c(3,4,5,7),2,2))
test <- rbind(test,mvrnorm(100, c(-5,5), matrix(c(3,0,0,7),2,2)))
png(filename="Weiszfeld.png",width = 600, height = 400, units = "px")
plot(test[501:600,],col="cornflowerblue",xlab="x", ylab="y",pch=2,xlim=c(-10,7),ylim=c(-8,14))
points(test[1:500,],col="grey")
legend("topright", legend=c("random start point", "spatial median"),col=c("black","blue"), bty="n", ncol=2, cex=1.8,pch=c(1,2,8,20))
for (k in 1:10){
  u <- matrix(nrow = 5000, ncol = 2)
  u[1,] <- c(runif(1,-10,6),runif(1,-6,12))
  points(u[1,1], u[1,2],pch=8)
  for (i in 2:5000){ # in one iteration
    nu <-c(0,0)
    dn <-0
    for (j in 1:600){
      di <-u[i-1,]-test[j,]
      wi <- sqrt((di[1])^2+(di[2])^2)
      dn <- dn + (1/wi)
      nu <- nu + test[j,]/wi
    }
    ui <- nu/dn
    u[i,]<-ui
  }
  lines(u, lty=6)
  points(u[5000,1], u[5000,2],pch=20,col="blue")
}
G <- Gmedian(test)
# title(sub=paste("estimated spatial median:(",round(u[5000,1],digits = 3),",",round(u[5000,2],digits = 3),")", ",result from Gmedian funciton:(",round(G[1],digits = 3),",",round(G[2],digits = 3),")"), cex=1.8)
dev.off()
################################################
# Figure 5 outlier in two dimension
# Simulation
set.seed(0)
# generate 10 sample and 5 outlier
sample <- mvrnorm(25, c(0,0), matrix(c(4,0,0,1),2,2))
out <- mvrnorm(5, c(0,-10), matrix(c(1,0,0,5),2,2))
# with 1 outlier
out_1 <- rbind(sample, out[1,])
# with 5 outlier
out_5 <- rbind(sample, out)

# estimation
c2 <- qchisq(0.95,2) # confidence region
# threotical
scatter <- ellipse(c(0,0),matrix(c(4,0,0,1),2,2),c2)
# no outlier
S <- SCM_covariance(sample,"Spatial") # covariance based on spatial sign
O <- SCM_covariance(sample,"OJA")
cov_0 <- ellipse(c(0,0),cov(sample),c2)
s_0 <- ellipse(c(0,0),S,c2)
o_0 <- ellipse(c(0,0),O,c2)
# 1 outlier
S_1 <- SCM_covariance(out_1,"Spatial")
O_1 <- SCM_covariance(out_1,"OJA")
cov_1 <- ellipse(c(0,0),cov(out_1),c2)
s_1 <- ellipse(c(0,0),S_1,c2)
o_1 <- ellipse(c(0,0),O_1,c2)
# 5 outlier
S_5 <- SCM_covariance(out_5,"Spatial")
O_5 <- SCM_covariance(out_5,"OJA")
cov_5 <- ellipse(c(0,0),cov(out_5),c2)
s_5 <- ellipse(c(0,0),S_5,c2)
o_5 <- ellipse(c(0,0),O_5,c2)
################################################
# visualization and save png
png(filename="scm1.png",width = 400, height = 400, units = "px")
color = c("black","red","chartreuse3","blue3")
plot(sample, asp=1, xlim=c(-8,8), pch=20, xlab="x", ylab="y")
lines(cov_0)
lines(s_0, col=color[2])
lines(o_0, col=color[3])
lines(scatter, lty=6, col=color[4])
dev.off()

png(filename="scm2.png",width = 400, height = 400, units = "px")
plot(out_1, asp=1, xlim=c(-10,10), ylim=c(-10,6), pch=20, xlab="x", ylab="y")
lines(cov_1)
lines(s_1, col=color[2])
lines(o_1, col=color[3])
lines(scatter, lty=6, col=color[4])
dev.off()

png(filename="scm3.png",width = 400, height = 400, units = "px")
plot(out_5, asp=1, xlim=c(-10,25), ylim=c(-15,12), pch=20, xlab="x", ylab="y")
lines(cov_5)
lines(s_5, col=color[2])
lines(o_5, col=color[3])
lines(scatter, lty=6, col=color[4])
legend("bottomright", legend=c("sample","spatial","oja","theoretical"),col=color, bty="n",lty=c(1,1,1,6),ncol=1,cex=1.8)
dev.off()
################################################
# # Figure 6
# Use statistical distance to detect outlier
# based on spatial & Oja
data("woodmod.dat")
outlier <- c(4, 6, 8, 19)
XX <- as.matrix(woodmod.dat) # 5 Variables
mean <- c(mean(XX[,1]),mean(XX[,2]),mean(XX[,3]),mean(XX[,4]),mean(XX[,5]))
CovS <- cov(XX)
cov_inverse <- solve(CovS)

l1 <- as.vector(Gmedian(XX))
S_spatial <- SCM_covariance(XX,"Spatial")
S_inverse_spatial <- solve(S_spatial)

d_classical <- rep(0,20)
d_robust <- rep(0,20)

for(i in 1:20){
  d_classical[i] <- sqrt(t(XX[i,]-mean)%*%cov_inverse%*%(XX[i,]-mean))
  d_robust[i] <- sqrt(t(XX[i,]-l1)%*%S_inverse_spatial%*%(XX[i,]-l1))
}

png(filename="outlier1.png",width = 400, height = 400, units = "px")
plot(d_classical, d_robust, xlab="Classical distance", ylab="Robust distance (Spatial)", pch=20, col="blue3")
points(d_classical[outlier], d_robust[outlier], pch=19, col="red")
lines(0:7,0:7, lty=6, col="dimgray")
legend("bottomright", legend=c("outliers"),col=c("red"),pch=20,ncol=1)
dev.off()

oja <- as.vector(ojaMedian(XX))
S_oja <- SCM_covariance(XX,"OJA")
S_inverse_oja <- solve(S_oja)

d_classical <- rep(0,20)
d_robust <- rep(0,20)
for(i in 1:20){
  d_classical[i] <- sqrt(t(XX[i,]-mean)%*%cov_inverse%*%(XX[i,]-mean))
  d_robust[i] <- sqrt(t(XX[i,]-oja)%*%S_inverse_oja%*%(XX[i,]-oja))
}

png(filename="outlier2.png",width = 400, height = 400, units = "px")
plot(d_classical, d_robust, xlab="Classical distance", ylab="Robust distance (Oja)", pch=20, col="blue3")
points(d_classical[outlier], d_robust[outlier], pch=19, col="red")
lines(0:7,0:7, lty=6, col="dimgray")
legend("bottomright", legend=c("outliers"),col=c("red"),pch=20,ncol=1)
dev.off()

png(filename="outlier3.png",width = 400, height = 400, units = "px")
woodm.fm=fit.models(list(Robust="covRob", Classical="covClassic"),
                    data=woodmod.dat)
ddPlot.covfm(woodm.fm, pch=4, col="purple", xlab="Classical distance",
             ylab="Robust distance", id.n=5)
dev.off()
################################################################

