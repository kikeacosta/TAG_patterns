## simple R-code for estimating log-mortality age-pattern 
## in two populations, in which the first 
## is simply a smooth function over age 
## and the second is the log-mortality of the first + 
## scaling factor +
## smooth age-factor
## issues: 
## grouping structures which could be eventually be different
## !! underdetermined linear system of equation
## 1) simulated data
## by C.G. Camarda 2021.06.11


## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## R-studio to get the same dir as for the .R file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(MortalitySmooth)

## assuming same exposures
## exposures population 1
e1 <- c(29588,29295,29002,28707,28410,28111,27809,27504,27196,26883,26566,26245,25919,25587,25251,24909,24561,24208,23849,23484,23113,22736,22353,21965,21570,21170,20800,20489,20227,20006,19819,19657,19514,19383,19258,19131,18996,18848,18681,18489,18266,18008,17711,17370,16984,16549,16065,15531,14948,14320,13648,12991,12394,11848,11340,10863,10410,9973,9547,9128,8713,8297,7878,7456,7029,6598,6164,5726,5289,4853,4422,4000,3589,3193,2815,2458,2131,1839,1579,1349,1147,969,815,681,566,468,384,313,254,204,163,129,102,80,62,48,36,28,21,15,11)
## exposures population 2
e2 <- e1#c(38233,38253,38274,38295,38316,38337,38357,38375,38392,38406,38418,38427,38431,38432,38429,38420,38406,38386,38359,38326,38286,38238,38182,38118,38045,37962,37935,38020,38205,38476,38823,39234,39696,40198,40726,41267,41806,42329,42820,43261,43637,43930,44123,44199,44142,43937,43571,43033,42314,41410,40318,39212,38243,37382,36603,35881,35197,34528,33856,33165,32437,31660,30820,29907,28914,27836,26670,25417,24080,22669,21192,19663,18098,16516,14935,13375,11892,10526,9273,8129,7090,6152,5309,4556,3887,3296,2779,2327,1937,1601,1314,1072,867,697,556,440,345,269,208,159,121)

## increasing sample size?
e1 <- e1*10
e2 <- e2*10
## dimension & age
m <- length(e1)
x <- 1:m-1

## plotting exposures
rany <- range(e1,e2)
plot(x, e1, t="l", lwd=2, ylim=rany)
lines(x, e2, col=2, lwd=2)

## known mortality age patterns

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
muT1 <- exp(eta1) + exp(eta2)
etaT1 <- log(muT1)

## plotting true log-mortality
## with components that we will disregard afterward
rany <- c(-12, 0)
plot(x, eta2, t="l", lwd=2, col=2, ylim=rany, axes=FALSE)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy)
axis(1);box()
lines(x, eta1, t="l", lwd=2, col=1)
lines(x, etaT1, col=4, lwd=2)


## population 2: etaT2 = etaT1 + c + delta
## true scaling factor
cT <- 0.5
## constructing true delta vector
ome <- pi/40
beta1 <- 0
beta2 <- 0.02
beta3 <- 1
betas <- c(beta1, beta2, beta3)
X2 <- cbind(1, x, cos(x*ome))
deltaT <- X2%*%betas
deltaT <-  deltaT - mean(deltaT)
plot(x, deltaT)
## true log-mortality for the second population (C19-2020)
etaT2 <- etaT1 + cT + deltaT
muT2 <- exp(etaT2)

# ## population 2 (C19-2020): higher mortality and more "sinousoid"
# ## (higher mortality at higher ages)
# alpha1 <- -5.8
# alpha2 <- -0.7
# alphas <- c(alpha1, alpha2)
# X1 <- cbind(1, x)
# eta1 <- X1%*%alphas
# ## component 2: exponential + cos
# ome <- pi/45
# beta1 <- -9.5
# beta2 <- 0.105
# beta3 <- 0.6
# betas <- c(beta1, beta2, beta3)
# X2 <- cbind(1, x, cos(x*ome))
# eta2 <- X2%*%betas
# ## overall (log-)mortality
# muT2 <- exp(eta1) + exp(eta2)
# etaT2 <- log(muT2)

# ## plotting true log-mortality
# ## with components that we will disregard afterward
# rany <- c(-10, 1)
# plot(x, eta2, t="l", lwd=2, col=2, ylim=rany, axes=FALSE)
# yy <- 10^seq(-20, 2)
# axis(2, at=log(yy), labels=yy)
# axis(1);box()
# lines(x, eta1, t="l", lwd=2, col=1)
# lines(x, etaT2, col=4, lwd=2)


## simulating deaths for both populations
## as realization from a Poisson distribution
d1T <- e1*muT1
d1 <- rpois(m, d1T)
lmx1 <- log(d1/e1)
d2T <- e2*muT2
d2 <- rpois(m, d2T)
lmx2 <- log(d2/e2)

## plotting simulated and true log-mortality
rany <- range(etaT1, etaT2)
plot(x, etaT1, t="l", lwd=2, col=2, ylim=rany, axes=FALSE)
yy <- 10^seq(-20, 2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
points(x, lmx1, col=2)
points(x, lmx2, col=3, pch=2)
lines(x, etaT2, col=3, lwd=2)
## assuming different grouping structure

## grouping for pop 1
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
d1g <- G1%*%d1
e1g <- G1%*%e1
lmx1g <- log(d1g/e1g)
## grouping for pop 2
# low2 <- seq(0, 80, 20)
# up2  <- c(seq(19, 79, 20), 100)
low2 <- low1
up2  <- up1
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
d2g <- G2%*%d2
e2g <- G2%*%e2
lmx2g <- log(d2g/e2g)

## plotting what we actually observed
rany <-range(etaT1, etaT2, lmx1, lmx2, lmx1g, lmx2g, finite=TRUE)
plot(x, etaT1, t="l", lwd=2, ylim=rany, axes=FALSE, col=2,
     xlab="age", ylab="mortality, log-scale",
     main="mortality age-pattern")
yy <- 10^seq(-7, 2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
points(x, lmx1, col=2)
for(i in 1:n1){
  segments(x0=low1[i], x1=up1[i],
           y0=lmx1g[i], y1=lmx1g[i], col=2, lwd=3)
}
lines(x, etaT2, lwd=2, col=3)
points(x, lmx2, col=3, pch=2)
for(i in 1:n2){
  segments(x0=low2[i], x1=up2[i],
           y0=lmx2g[i], y1=lmx2g[i], col=3, lwd=3)
}


## build composite matrices for the two population
C1 <- G1
for(i in 1:n1){
  age.low.i <- low1[i]
  age.up.i <- up1[i]
  whi <- which(x>=age.low.i & x<=age.up.i)
  C1[i,whi] <- e1[whi]
}
C2 <- G2
for(i in 1:n2){
  age.low.i <- low2[i]
  age.up.i <- up2[i]
  whi <- which(x>=age.low.i & x<=age.up.i)
  C2[i,whi] <- e2[whi]
}


## estimation
y <- c(d1g, d2g)
C <- adiag(C1, C2)
y1.st0 <- rep(d1g/len1, len1)
fit1.0 <- Mort1Dsmooth(x=x, y=y1.st0,
                       offset=log(e1),
                       method=3, lambda=10^3)
plot(fit1.0)
y2.st0 <- rep(d2g/len2, len2)
fit2.0 <- Mort1Dsmooth(x=x, y=y2.st0,
                       offset=log(e2),
                       method=3, lambda=10^3)
plot(fit2.0)

eta.st <- c(fit1.0$logmortality, fit2.0$logmortality)
## penalty stuff
Deta1 <- diff(diag(m), diff=2)
tDDeta1 <- t(Deta1)%*%Deta1
Ddelta <- diff(diag(m), diff=2)
tDDdelta <- t(Ddelta)%*%Ddelta

lambda.eta1 <- 10^3
lambda.delta <- 10^3.5

P <- adiag(lambda.eta1*tDDeta1,
           0,
           lambda.delta*tDDdelta)
## ridge penalty
Pr <- 1e-4*diag(ncol(P))

## model matrix
U0 <- rbind(diag(m), diag(m))
U1 <- matrix(c(rep(0,m), rep(1,m)), 2*m)
U2 <- rbind(0*diag(m), diag(m))
U <- cbind(U0,U1,U2)

eta <- eta.st
max.it <- 100
for(it in 1:max.it){
  gamma   <- exp(eta)
  mu      <- c(C %*% gamma)
  X       <- ( C * ((1 / mu) %*% t(gamma)) ) %*% U
  w       <- as.vector(mu)
  r       <- y - mu + C %*% (gamma * eta)
  G       <- t(X) %*% (w * X) 
  GpP     <- G + P + Pr
  tXr     <- t(X) %*% r
  betas   <- solve(GpP, tXr)
  eta.old <- eta
  eta     <- U%*%betas
  dif.eta <- max(abs((eta - eta.old)/eta.old) )
  if(dif.eta < 1e-04 & it > 4) break
  cat(it, dif.eta, "\n")
}
eta1.hat <- betas[1:m]
c.hat <- betas[m+1]
delta.hat <- betas[1:m+1+m]
eta2.hat <- eta1.hat + c.hat + delta.hat


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
se.eta1 <- se.betas[1:m]
se.c <- se.betas[m+1]
se.delta <- se.betas[1:m+1+m]
V.eta12 <- U %*% Vbetas %*% t(U)
se.eta2 <- sqrt(diag(V.eta12))[1:m+m]


eta1.hatL <- eta1.hat - 2*se.eta1
eta1.hatU <- eta1.hat + 2*se.eta1
eta2.hatL <- eta2.hat - 2*se.eta2
eta2.hatU <- eta2.hat + 2*se.eta2
c.hatL <- c.hat - 2*se.c
c.hatU <- c.hat + 2*se.c
delta.hatL <- delta.hat - 2*se.delta
delta.hatU <- delta.hat + 2*se.delta

par(mfrow=c(1,3))
## eta1
rany <-range(etaT1, etaT2, lmx1g, lmx2g, finite=TRUE)
ranx <- range(x)
plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE, 
     xlab="age", ylab="mortality, log-scale",
     main="mortality age-pattern")
yy <- 10^seq(-7, 2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
xx <- c(x, rev(x))
yy <- c(eta1.hatL, rev(eta1.hatU))
polygon(xx, yy, col=adjustcolor(2, 0.5), border=adjustcolor(2, 0.5))
lines(x, eta1.hat, col=2, lwd=4)
lines(x, etaT1, col=4, lwd=2, lty=2)
for(i in 1:n1){
  segments(x0=low1[i], x1=up1[i],
           y0=lmx1g[i], y1=lmx1g[i], col=3, lwd=3)
}
## eta2
plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE, 
     xlab="age", ylab="mortality, log-scale",
     main="mortality age-pattern")
yy <- 10^seq(-7, 2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
xx <- c(x, rev(x))
yy <- c(eta2.hatL, rev(eta2.hatU))
polygon(xx, yy, col=adjustcolor(2, 0.5), border=adjustcolor(4, 0.5))
lines(x, eta2.hat, col=2, lwd=4)
lines(x, etaT2, col=4, lwd=2, lty=2)
for(i in 1:n1){
  segments(x0=low2[i], x1=up2[i],
           y0=lmx2g[i], y1=lmx2g[i], col=3, lwd=3)
}

ranx <- range(x)
rany <- range(delta.hatL,delta.hatU)
plot(1, 1, t="n", xlim=ranx, ylim=rany,
     main=paste("c = ", signif(c.hat,4), "+ -", signif(2*se.c,4)), axes=FALSE,
     xlab="ages", ylab="delta")
axis(2)
axis(1);box()
abline(h=0, col=8, lty=3, lwd=2)
xx <- c(x, rev(x))
yy <- c(delta.hatL, rev(delta.hatU))
polygon(xx, yy, col=adjustcolor(2, 0.5), border=adjustcolor(4, 0.5))
lines(x, delta.hat, col=2, lwd=4)
lines(x, deltaT, col=4, lwd=2, lty=2)
legend("topleft", inset=0.1,
       legend=c("True", "Simulated+Grouped", "Fit + 95% CI"),
       col=c(4, 3, adjustcolor(2, 0.5)), pch=c(NA, NA, 15), 
       pt.cex=2, lty=c(2,1,1),
       lwd=c(2,3,4), cex=1.5)

par(mfrow=c(1,1))







