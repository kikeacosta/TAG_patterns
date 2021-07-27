## simple R-code for estimating log-mortality age-pattern 
## in a populations by sex in two different years
## in which for a given sex,the first is simply a smooth function over age 
## and the second is the log-mortality of the first + 
## scaling factor +
## smooth age-factor
## for a given year, males mortality is female mortality 
## + smooth sex factor over age

## issues: 
## grouping structures which could be eventually be different
## scaling factor different by sex?
## smooth-age factor different by sex?
## !! underdetermined linear system of equation
## 1) simulated data
## by C.G. Camarda 2021.07.22

## OPTION 2: 
# sex specific scaling and age-factor

## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## R-studio to get the same dir as for the .R file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(MortalitySmooth)
library(magic)

x <- 0:100
m <- length(x)

setwd("~/WORK/TAG_patterns/")
## loading exposures
EE0 <- read.csv("Output/offsets.csv", header=TRUE)
## remove only totals
EE1 <- subset(EE0, Sex!="t")
## remove 2021
EE1 <- subset(EE1, Year<=2020)
## take a given country: 
cou <- "France"
EE2 <- subset(EE1, Country=="France")

## females
EE2F <- subset(EE2, Sex=="f")
ee1F <- subset(EE2F, Year<2020)
e1F <- rowSums(matrix(ee1F$Population, m))
e2F <- subset(EE2F, Year==2020)$Population
plot(x, e1F)
lines(x, e2F)
## males
EE2M <- subset(EE2, Sex=="m")
ee1M <- subset(EE2M, Year<2020)
e1M <- rowSums(matrix(ee1M$Population, m))
e2M <- subset(EE2M, Year==2020)$Population
plot(x, e1M)
lines(x, e2M)

## plotting exposures
rany <- range(e1F,e2F,e1M,e2M)
plot(x, e1F, t="l", col=2, lwd=2, ylim=rany)
lines(x, e2F, col=2, lwd=2, lty=2)
lines(x, e1M, t="l", col=4, lwd=2, ylim=rany)
lines(x, e2M, col=4, lwd=2, lty=2)

## known mortality age patterns for females in year 1

## for illustrative purposes: 
## assuming as combination of two components (in a log-scale)
## component 1: simple exponential

## population 1 (pre-C19): lower mortality and more gompertzian
alpha1 <- -6.8
alpha2 <- -0.7
alphas <- c(alpha1, alpha2)
X1 <- cbind(1, x)
eta1 <- X1%*%alphas
## component 2: exponential + cos
ome <- pi/45
beta1 <- -13
beta2 <- 0.12
beta3 <- 0.8
betas <- c(beta1, beta2, beta3)
X2 <- cbind(1, x, cos(x*ome))
eta2 <- X2%*%betas
## overall (log-)mortality
muT1F <- exp(eta1) + exp(eta2)
etaT1F <- log(muT1F)

## plotting true log-mortality
## with components that we will disregard afterward
rany <- c(-12, 0)
plot(x, eta2, t="l", lwd=2, col=2, ylim=rany, axes=FALSE)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy)
axis(1);box()
lines(x, eta1, t="l", lwd=2, col=1)
lines(x, etaT1F, col=2, lwd=2)

## males in year 1: etaT1M(x) = etaT1F(x) + s(x)
## constructing true s(x) vector
sT <-  dgamma(x, shape=20, scale=1.5)*30 + 0.5
plot(x, sT)
## true log-mortality for the males in year 1
etaT1M <- etaT1F + sT
muT1M <- exp(etaT1M)

## plot males+females in year 1
rany <- range(etaT1F, etaT1M)
plot(x, etaT1F, t="l", lwd=2, col=2, ylim=rany, axes=FALSE)
yy <- 10^seq(-20, 2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
lines(x, etaT1M, col=4, lwd=2)


## females in year 2: etaT2F(x) = etaT1F(x) + cF + deltaF(x)
## true scaling factor
cTF <- 0.4
## constructing true delta vector
ome <- pi/40
beta1 <- 0.3
beta2 <- 0.02
beta3 <- 0.5
betas <- c(beta1, beta2, beta3)
X2 <- cbind(1, x, cos(x*ome))
deltaTF <- X2%*%betas
deltaTF <-  deltaTF - mean(deltaTF)
plot(x, deltaTF)
## true log-mortality for the second population (C19-2020)
etaT2F <- etaT1F + cTF + deltaTF
muT2F <- exp(etaT2F)


## males in year 2: etaTM2(x) = etaT1M(x) + c + delta(x)
## true scaling factor
cTM <- 0.5
## constructing true delta vector
ome <- pi/40
beta1 <- 0.3
beta2 <- 0.02
beta3 <- -0.5
betas <- c(beta1, beta2, beta3)
X2 <- cbind(1, x, cos(x*ome))
deltaTM <- X2%*%betas
deltaTM <-  deltaTM - mean(deltaTM)
plot(x, deltaTM)
lines(x, deltaTF)
## true log-mortality for the second population (C19-2020)
etaT2M <- etaT1M + cTM + deltaTM
muT2M <- exp(etaT2M)


rany <- range(etaT1F, etaT1M, etaT2F, etaT2M)
plot(x, etaT1F, t="l", lwd=2, col=2, ylim=rany, axes=FALSE)
yy <- 10^seq(-20, 2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
lines(x, etaT1M, col=4, lwd=2)
lines(x, etaT2F, col=2, lwd=2, lty=2)
lines(x, etaT2M, col=4, lwd=2, lty=2)


## simulating deaths for both populations
## as realization from a Poisson distribution
d1TF <- e1F*muT1F
d1TM <- e1M*muT1M
d2TF <- e2F*muT2F
d2TM <- e2M*muT2M

d1F <- rpois(m, d1TF)
d1M <- rpois(m, d1TM)
d2F <- rpois(m, d2TF)
d2M <- rpois(m, d2TM)

lmx1F <- log(d1F/e1F)
lmx1M <- log(d1M/e1M)
lmx2F <- log(d2F/e2F)
lmx2M <- log(d2M/e2M)


## plotting simulated and true log-mortality
rany <- range(etaT1F, etaT1M, etaT2F, etaT2M)
plot(x, etaT1F, t="l", lwd=2, col=2, ylim=rany, axes=FALSE)
yy <- 10^seq(-20, 2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
lines(x, etaT1M, col=4, lwd=2)
lines(x, etaT2F, col=2, lwd=2, lty=2)
lines(x, etaT2M, col=4, lwd=2, lty=2)
points(x, lmx1F, col=2, pch=16)
points(x, lmx1M, col=4, pch=16)
points(x, lmx2F, col=2, pch=17)
points(x, lmx2M, col=4, pch=17)

## assuming different grouping structure

## grouping for year 1
low1 <- c(0, 1, seq(5, 90,5))
up1  <- c(0, seq(4, 89,5), 100)
n1 <- length(low1)
age.gr1 <- paste(low1,up1,sep="-")
len1 <- up1-low1+1
## matrix to group both deaths and exposures
## (grouped expoures will be used only for plotting purposes)
G1 <- matrix(0, n1, m)
rownames(G1) <- age.gr1
colnames(G1) <- x
for(i in 1:n1){
  age.low.i <- low1[i]
  age.up.i <- up1[i]
  whi <- which(x>=age.low.i & x<=age.up.i)
  G1[i,whi] <- 1
}
all(colSums(G1)==1)

## what we would observed in an actual world for pop 1
d1Fg <- G1%*%d1F
e1Fg <- G1%*%e1F
lmx1Fg <- log(d1Fg/e1Fg)
d1Mg <- G1%*%d1M
e1Mg <- G1%*%e1M
lmx1Mg <- log(d1Mg/e1Mg)



## grouping for pop 2
low2 <- seq(0, 80, 20)
up2  <- c(seq(19, 79, 20), 100)
# low2 <- low1
# up2  <- up1
n2 <- length(low2)
age.gr2 <- paste(low2,up2,sep="-")
len2 <- up2-low2+1
G2 <- matrix(0, n2, m)
rownames(G2) <- age.gr2
colnames(G2) <- x
for(i in 1:n2){
  age.low.i <- low2[i]
  age.up.i <- up2[i]
  whi <- which(x>=age.low.i & x<=age.up.i)
  G2[i,whi] <- 1
}
all(colSums(G2)==1)
## what we would observed in an actual world for pop 2
d2Fg <- G2%*%d2F
e2Fg <- G2%*%e2F
lmx2Fg <- log(d2Fg/e2Fg)
d2Mg <- G2%*%d2M
e2Mg <- G2%*%e2M
lmx2Mg <- log(d2Mg/e2Mg)

## plotting what we actually observed
rany <- range(etaT1F, etaT1M, etaT2F, etaT2M)
plot(x, etaT1F, t="n", lwd=2, col=2, ylim=rany, axes=FALSE)
yy <- 10^seq(-20, 2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
lines(x, etaT1M, col=4, lwd=2)
lines(x, etaT2F, col=2, lwd=2, lty=2)
lines(x, etaT2M, col=4, lwd=2, lty=2)
# points(x, lmx1F, col=2, pch=16)
# points(x, lmx1M, col=4, pch=16)
# points(x, lmx2F, col=2, pch=17)
# points(x, lmx2M, col=4, pch=17)
for(i in 1:n1){
  segments(x0=low1[i], x1=up1[i],
           y0=lmx1Fg[i], y1=lmx1Fg[i], col=2, lwd=3)
  segments(x0=low1[i], x1=up1[i],
           y0=lmx1Mg[i], y1=lmx1Mg[i], col=4, lwd=3)
}
for(i in 1:n2){
  segments(x0=low2[i], x1=up2[i],
           y0=lmx2Fg[i], y1=lmx2Fg[i], col=2, lwd=3)
  segments(x0=low2[i], x1=up2[i],
           y0=lmx2Mg[i], y1=lmx2Mg[i], col=4, lwd=3)
}


## build composite matrices for the two population
C1F <- C1M <- G1
for(i in 1:n1){
  age.low.i <- low1[i]
  age.up.i <- up1[i]
  whi <- which(x>=age.low.i & x<=age.up.i)
  C1F[i,whi] <- e1F[whi]
  C1M[i,whi] <- e1M[whi]
}
C2F <- C2M <- G2
for(i in 1:n2){
  age.low.i <- low2[i]
  age.up.i <- up2[i]
  whi <- which(x>=age.low.i & x<=age.up.i)
  C2F[i,whi] <- e2F[whi]
  C2M[i,whi] <- e2M[whi]
}


## estimation
y <- c(d1Fg, d1Mg, d2Fg, d2Mg)
C <- adiag(C1F, C1M, C2F, C2M)
y1F.st0 <- rep(d1Fg/len1, len1)
fit1F.0 <- Mort1Dsmooth(x=x, y=y1F.st0,
                        offset=log(e1F),
                        method=3, lambda=10^4)
plot(fit1F.0)
y1M.st0 <- rep(d1Mg/len1, len1)
fit1M.0 <- Mort1Dsmooth(x=x, y=y1M.st0,
                        offset=log(e1M),
                        method=3, lambda=10^4)
plot(fit1M.0)
y2F.st0 <- rep(d2Fg/len2, len2)
fit2F.0 <- Mort1Dsmooth(x=x, y=y2F.st0,
                        offset=log(e2F),
                        method=3, lambda=10^4)
plot(fit2F.0)
y2M.st0 <- rep(d2Mg/len2, len2)
fit2M.0 <- Mort1Dsmooth(x=x, y=y2M.st0,
                        offset=log(e2M),
                        method=3, lambda=10^4)
plot(fit2M.0)



eta.st <- c(fit1F.0$logmortality, fit1M.0$logmortality,
            fit2F.0$logmortality, fit2M.0$logmortality)
## penalty stuff
Deta1F <- diff(diag(m), diff=2)
tDDeta1F <- t(Deta1F)%*%Deta1F
Ds <- diff(diag(m), diff=2)
tDDs <- t(Ds)%*%Ds
DdeltaF <- diff(diag(m), diff=2)
tDDdeltaF <- t(DdeltaF)%*%DdeltaF
DdeltaM <- diff(diag(m), diff=2)
tDDdeltaM <- t(DdeltaM)%*%DdeltaM

lambda.eta1F <- 10^3
lambda.s <- 10^2.5
lambda.deltaF <- 10^3.5
lambda.deltaM <- 10^3.5

P <- adiag(lambda.eta1F*tDDeta1F,
           lambda.s*tDDs,
           0,
           lambda.deltaF*tDDdeltaF,
           0,
           lambda.deltaM*tDDdeltaM)
## ridge penalty
Pr <- 1e-4*diag(ncol(P))

## model matrix
Ueta1F <- kronecker(rep(1,4), diag(m))
Us <- rbind(0*diag(m), diag(m), 0*diag(m), diag(m))
UcF <- matrix(c(rep(0,2*m), rep(1,m), rep(0,m)), 4*m)
UdeltaF <- rbind(0*diag(m), 0*diag(m), diag(m), 0*diag(m))
UcM <- matrix(c(rep(0,3*m), rep(1,m)), 4*m)
UdeltaM <- rbind(0*diag(m), 0*diag(m), 0*diag(m), diag(m))
U <- cbind(Ueta1F,Us,UcF,UdeltaF,UcM,UdeltaM)

## constraining delta to sum up to 0
H <- rbind(c(rep(0,2*m+1), rep(1,m), rep(0,m+1)),
           c(rep(0,2*m+1), rep(0,m), 0, rep(1,m)))
## check
H%*%c(etaT1F, sT, cTF, deltaTF, cTM, deltaTM)
kappa <- c(0,0)

eta <- eta.st
max.it <- 100
zeros <- matrix(0,2,2)
for(it in 1:max.it){
  gamma   <- exp(eta)
  mu      <- c(C %*% gamma)
  X       <- ( C * ((1 / mu) %*% t(gamma)) ) %*% U
  w       <- as.vector(mu)
  r       <- y - mu + C %*% (gamma * eta)
  G       <- t(X) %*% (w * X) 
  GpP     <- G + P + Pr
  tXr     <- t(X) %*% r
  ## adding constraints
  LHS <- rbind(cbind(GpP, t(H)),
               cbind(H, zeros))
  RHS <- matrix(c(tXr, kappa), ncol=1)
  coeff   <- solve(LHS, RHS)
  
  betas   <- coeff[1:(m*4+2)]
  eta.old <- eta
  eta     <- U%*%betas
  dif.eta <- max(abs((eta - eta.old)/eta.old) )
  if(dif.eta < 1e-04 & it > 4) break
  cat(it, dif.eta, "\n")
}
eta1F.hat <- betas[1:m]
s.hat <- betas[1:m+m]
cF.hat <- betas[2*m+1]
deltaF.hat <- betas[1:m+2*m+1]
cM.hat <- betas[3*m+2]
deltaM.hat <- betas[1:m+3*m+2]


eta1M.hat <- eta1F.hat + s.hat
eta2F.hat <- eta1F.hat + cF.hat + deltaF.hat
eta2M.hat <- eta1F.hat + s.hat + cM.hat + deltaM.hat # eta1M.hat + cM.hat + deltaM.hat

# par(mfrow=c(1,3))
# ## eta1
# rany <-range(etaT1, etaT2, lmx1g, lmx2g, finite=TRUE)
# ranx <- range(x)
# plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE, 
#      xlab="age", ylab="mortality, log-scale",
#      main="mortality age-pattern")
# yy <- 10^seq(-7, 2)
# axis(2, at=log(yy), labels=yy)
# axis(1);box()
# lines(x, eta1.hat, col=2, lwd=4)
# lines(x, etaT1, col=4, lwd=2, lty=2)
# for(i in 1:n1){
#   segments(x0=low1[i], x1=up1[i],
#            y0=lmx1g[i], y1=lmx1g[i], col=4, lwd=3)
# }
# ## eta2
# plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE, 
#      xlab="age", ylab="mortality, log-scale",
#      main="mortality age-pattern")
# yy <- 10^seq(-7, 2)
# axis(2, at=log(yy), labels=yy)
# axis(1);box()
# lines(x, eta2.hat, col=2, lwd=4)
# lines(x, etaT2, col=4, lwd=2, lty=2)
# for(i in 1:n1){
#   segments(x0=low2[i], x1=up2[i],
#            y0=lmx2g[i], y1=lmx2g[i], col=4, lwd=3)
# }
# 
# plot(x, delta.hat, col=4, lwd=4, t="l",
#      main=paste("c = ", signif(c.hat,4)))
# abline(h=0, col=8, lty=3, lwd=2)
# lines(x, deltaT, col=4, lty=2, lwd=2)
# par(mfrow=c(1,1))

H0 <- solve(GpP)
Vbetas <- H0 %*% G %*% H0
se.betas <- sqrt(diag(Vbetas))
se.eta1F <- se.betas[1:m]
se.s <- se.betas[1:m+m]
se.cF <- se.betas[2*m+1]
se.deltaF <- se.betas[1:m+1+2*m]
se.cM <- se.betas[3*m+2]
se.deltaM <- se.betas[1:m+3*m+2]

V.eta <- U %*% Vbetas %*% t(U)
se.eta <- sqrt(diag(V.eta))
se.eta1F <- se.eta[1:m]
se.eta1M <- se.eta[1:m+m]
se.eta2F <- se.eta[1:m+2*m]
se.eta2M <- se.eta[1:m+3*m]

eta1F.hatL <- eta1F.hat - 2*se.eta1F
eta1F.hatU <- eta1F.hat + 2*se.eta1F
eta1M.hatL <- eta1M.hat - 2*se.eta1M
eta1M.hatU <- eta1M.hat + 2*se.eta1M

eta2F.hatL <- eta2F.hat - 2*se.eta2F
eta2F.hatU <- eta2F.hat + 2*se.eta2F
eta2M.hatL <- eta2M.hat - 2*se.eta2M
eta2M.hatU <- eta2M.hat + 2*se.eta2M

s.hatL <- s.hat - 2*se.s
s.hatU <- s.hat + 2*se.s
cF.hatL <- cF.hat - 2*se.cF
cF.hatU <- cF.hat + 2*se.cF
deltaF.hatL <- deltaF.hat - 2*se.deltaF
deltaF.hatU <- deltaF.hat + 2*se.deltaF
cM.hatL <- cM.hat - 2*se.cM
cM.hatU <- cM.hat + 2*se.cM
deltaM.hatL <- deltaM.hat - 2*se.deltaM
deltaM.hatU <- deltaM.hat + 2*se.deltaM


par(mfrow=c(3,2))
## year1
rany <-range(etaT1F, etaT1M, etaT2F, etaT2M, lmx1Fg, lmx1Mg, lmx2Fg, lmx2Mg, finite=TRUE)
ranx <- range(x)
plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE, 
     xlab="age", ylab="mortality, log-scale",
     main="mortality age-pattern")
yy <- 10^seq(-7, 2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
xx <- c(x, rev(x))
yy <- c(eta1F.hatL, rev(eta1F.hatU))
polygon(xx, yy, col=adjustcolor(2, 0.5), border=adjustcolor(2, 0.5))
lines(x, eta1F.hat, col=2, lwd=4)
lines(x, etaT1F, col=4, lwd=2, lty=2)
for(i in 1:n1){
  segments(x0=low1[i], x1=up1[i],
           y0=lmx1Fg[i], y1=lmx1Fg[i], col=3, lwd=3)
}
yy <- c(eta1M.hatL, rev(eta1M.hatU))
polygon(xx, yy, col=adjustcolor(2, 0.5), border=adjustcolor(2, 0.5))
lines(x, eta1M.hat, col=2, lwd=4)
lines(x, etaT1M, col=4, lwd=2, lty=2)
for(i in 1:n1){
  segments(x0=low1[i], x1=up1[i],
           y0=lmx1Mg[i], y1=lmx1Mg[i], col=3, lwd=3)
}
## year2
rany <-range(etaT1F, etaT1M, etaT2F, etaT2M, lmx1Fg, lmx1Mg, lmx2Fg, lmx2Mg, finite=TRUE)
ranx <- range(x)
plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE, 
     xlab="age", ylab="mortality, log-scale",
     main="mortality age-pattern")
yy <- 10^seq(-7, 2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
xx <- c(x, rev(x))
yy <- c(eta2F.hatL, rev(eta2F.hatU))
polygon(xx, yy, col=adjustcolor(2, 0.5), border=adjustcolor(2, 0.5))
lines(x, eta2F.hat, col=2, lwd=4)
lines(x, etaT2F, col=4, lwd=2, lty=2)
for(i in 1:n2){
  segments(x0=low2[i], x1=up2[i],
           y0=lmx2Fg[i], y1=lmx2Fg[i], col=3, lwd=3)
}
yy <- c(eta2M.hatL, rev(eta2M.hatU))
polygon(xx, yy, col=adjustcolor(2, 0.5), border=adjustcolor(2, 0.5))
lines(x, eta2M.hat, col=2, lwd=4)
lines(x, etaT2M, col=4, lwd=2, lty=2)
for(i in 1:n2){
  segments(x0=low2[i], x1=up2[i],
           y0=lmx2Mg[i], y1=lmx2Mg[i], col=3, lwd=3)
}
## s(x)
ranx <- range(x)
rany <- range(s.hatL,s.hatU)
plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE,
     xlab="ages", ylab="s")
axis(2)
axis(1);box()
abline(h=0, col=8, lty=3, lwd=2)
xx <- c(x, rev(x))
yy <- c(s.hatL, rev(s.hatU))
polygon(xx, yy, col=adjustcolor(2, 0.5), border=adjustcolor(4, 0.5))
lines(x, s.hat, col=2, lwd=4)
lines(x, sT, col=4, lwd=2, lty=2)
## deltaF(x)
ranx <- range(x)
rany <- range(deltaF.hatL,deltaF.hatU)
plot(1, 1, t="n", xlim=ranx, ylim=rany,
     main=paste("c = ", signif(cF.hat,4), "+ -", signif(2*se.cF,4)), axes=FALSE,
     xlab="ages", ylab="delta")
axis(2)
axis(1);box()
abline(h=0, col=8, lty=3, lwd=2)
xx <- c(x, rev(x))
yy <- c(deltaF.hatL, rev(deltaF.hatU))
polygon(xx, yy, col=adjustcolor(2, 0.5), border=adjustcolor(4, 0.5))
lines(x, deltaF.hat, col=2, lwd=4)
lines(x, deltaTF, col=4, lwd=2, lty=2)
## deltaM(x)
ranx <- range(x)
rany <- range(deltaM.hatL,deltaM.hatU)
plot(1, 1, t="n", xlim=ranx, ylim=rany,
     main=paste("c = ", signif(cM.hat,4), "+ -", signif(2*se.cM,4)), axes=FALSE,
     xlab="ages", ylab="delta")
axis(2)
axis(1);box()
abline(h=0, col=8, lty=3, lwd=2)
xx <- c(x, rev(x))
yy <- c(deltaM.hatL, rev(deltaM.hatU))
polygon(xx, yy, col=adjustcolor(2, 0.5), border=adjustcolor(4, 0.5))
lines(x, deltaM.hat, col=2, lwd=4)
lines(x, deltaTM, col=4, lwd=2, lty=2)
legend("topleft", inset=0.1,
       legend=c("True", "Simulated+Grouped", "Fit + 95% CI"),
       col=c(4, 3, adjustcolor(2, 0.5)), pch=c(NA, NA, 15), 
       pt.cex=2, lty=c(2,1,1),
       lwd=c(2,3,4), cex=1.5)


par(mfrow=c(1,1))




## colors
## colors
col1F <- "red"
col1FT <- adjustcolor(col1F, 0.5)
col2F <- "orange"
col2FT <- adjustcolor(col2F, 0.5)

col1M <- "green2"
col1MT <- adjustcolor(col1M, 0.5)
col2M <- "cyan3"
col2MT <- adjustcolor(col2M, 0.5)

coldF <- "orange3"
coldFT <- adjustcolor(coldF, 0.5)
colcF <- "orange3"
colcFT <- adjustcolor(colcF, 0.1)

coldM <- "deepskyblue2"
coldMT <- adjustcolor(coldM, 0.5)
colcM <- "deepskyblue2"
colcMT <- adjustcolor(colcM, 0.1)

cols <- "coral"
colsT <- adjustcolor(cols, 0.5)



pdf("/home/gccamarda/WORK/TAG_patterns/Slides/Figures/SimulEta3FSexSat.pdf",
    width = 8, height = 8)
par(mar=c(5,5,0.5, 0.5))
rany <-range(etaT1F, etaT2F, 
             etaT1M, etaT2M, finite=TRUE)
ranx <- range(x)
plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE, 
     xlab="", ylab="")
yy <- 10^seq(-7, 2)
axis(2, at=log(yy), labels=yy, las=2, cex.lab=1.3)
axis(1, cex.lab=1.3);box()
abline(h=log(yy), lty=2, col="grey85")
abline(v=seq(-100, 100, 20), lty=2, col="grey85")
mtext("age", 1, cex=1.8, line=3)
mtext(expression(paste("log-mortality, ", eta)), 2, cex=1.8, line=3)
lines(x, etaT1F, col=col1F, lwd=6)
lines(x, etaT1F+cTF+deltaTF, col=col2F, lwd=6, lty=1)
text(80, -4.5, expression(paste(eta^F1)), col=col1F, cex=3)
text(25, -7.5, expression(paste(eta^F2, "=", eta^F1, "+", c^F, "+", delta^F)), col=col2F, cex=3)
dev.off()

pdf("/home/gccamarda/WORK/TAG_patterns/Slides/Figures/SimulEta3MSexSat.pdf",
    width = 8, height = 8)
par(mar=c(5,5,0.5, 0.5))
rany <-range(etaT1F, etaT2F, 
             etaT1M, etaT2M, finite=TRUE)
ranx <- range(x)
plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE, 
     xlab="", ylab="")
yy <- 10^seq(-7, 2)
axis(2, at=log(yy), labels=yy, las=2, cex.lab=1.3)
axis(1, cex.lab=1.3);box()
abline(h=log(yy), lty=2, col="grey85")
abline(v=seq(-100, 100, 20), lty=2, col="grey85")
mtext("age", 1, cex=1.8, line=3)
mtext(expression(paste("log-mortality, ", eta)), 2, cex=1.8, line=3)
lines(x, etaT1F+ sT, col=col1M, lwd=6, lty=1)
lines(x, etaT1F+ sT + cTM+deltaTM, col=col2M, lwd=6, lty=1)
text(80, -6,  expression(paste(eta^M1, "=", eta^F1, "+", s)), col=col1M, cex=3)
text(25, -5, expression(paste(eta^M2, "=", eta^M1, "+", c^M, "+", delta^M)), col=col2M, cex=3)
dev.off()


pdf("/home/gccamarda/WORK/TAG_patterns/Slides/Figures/SimulDelta3SexSat.pdf",
    width = 8, height = 8)
par(mar=c(5,5,0.5, 0.5))
ranx <- range(-10, x)
rany <- range(cTM, cTF, deltaTF, deltaTM)
plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE,
     xlab="", ylab="")
axis(2, cex.lab=1.3, las=2)
axis(1, cex.lab=1.3);box()
abline(v=seq(-100, 100, 20), lty=2, col="grey85")
abline(h=seq(-100, 100, 0.5), lty=2, col="grey85")
abline(h=0, col="grey40", lty=3, lwd=3)
mtext("age", 1, cex=1.8, line=3, at=50)
mtext(expression(paste(c, " , ", delta)), 2, cex=1.8, line=3)
abline(v=-1, col=8, lwd=4, lty=2)
rect(xleft=-15, xright=-1, ybottom=-20, ytop=20, col=colcFT, border = colcFT)
points(-9, cTF, col=colcF, pch=16, cex=3, lwd=3)
lines(x, deltaTF, col=coldF, lwd=6)
points(-7, cTM, col=colcM, pch=16, cex=3, lwd=3)
lines(x, deltaTM, col=coldM, lwd=6)

mtext("c", 1, at=-8, cex=2, line=3, col=1)

text(40, -1, expression(delta^F), col=coldF, cex=3)
text(40, 0.7, expression(delta^M), col=coldM, cex=3)

text(-9, 0.25, expression(c^F), col=coldF, cex=3)
text(-6, 0.65, expression(c^M), col=coldM, cex=3)

dev.off()

pdf("/home/gccamarda/WORK/TAG_patterns/Slides/Figures/Simuls3SexSat.pdf",
    width = 8, height = 8)
par(mar=c(5,5,0.5, 0.5))
ranx <- range(x)
rany <- range(0, sT)
plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE,
     xlab="", ylab="")
axis(2, cex.lab=1.3, las=2)
axis(1, cex.lab=1.3);box()
abline(v=seq(-100, 100, 20), lty=2, col="grey85")
abline(h=seq(-100, 100, 0.5), lty=2, col="grey85")
abline(h=0, col="grey40", lty=3, lwd=3)
mtext("age", 1, cex=1.8, line=3, at=50)
mtext(expression(paste(s)), 2, cex=1.8, line=3)
lines(x, sT, col=cols, lwd=6)
text(92, 0.6, expression(s), col=cols, cex=3)
dev.off()