## simple R-code for estimating log-mortality age-pattern 
## in a population, as sum of the 
## log-mortality of another population + 
## scaling factor +
## smooth age-factor
## issue: grouping structures which could be eventually be different
## in the two populations
## 1) simulated data
## by C.G. Camarda 2021.06.11


## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## R-studio to get the same dir as for the .R file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(MortalitySmooth)

## function to compute matrix for computing 1st derivative for equally-spaced x
Dfun <- function(m){
  D0 <- diff(diag(m), diff=1)
  D1 <- diff(diag(m), diff=1, lag=2)*0.5
  D <- rbind(D0[1,],
             D1,
             D0[nrow(D0),])
  return(D)
}


## assuming same exposures
## exposures population 1
e1 <- c(29588,29295,29002,28707,28410,28111,27809,27504,27196,26883,26566,26245,25919,25587,25251,24909,24561,24208,23849,23484,23113,22736,22353,21965,21570,21170,20800,20489,20227,20006,19819,19657,19514,19383,19258,19131,18996,18848,18681,18489,18266,18008,17711,17370,16984,16549,16065,15531,14948,14320,13648,12991,12394,11848,11340,10863,10410,9973,9547,9128,8713,8297,7878,7456,7029,6598,6164,5726,5289,4853,4422,4000,3589,3193,2815,2458,2131,1839,1579,1349,1147,969,815,681,566,468,384,313,254,204,163,129,102,80,62,48,36,28,21,15,11)
## exposures population 2
e2 <- e1#c(38233,38253,38274,38295,38316,38337,38357,38375,38392,38406,38418,38427,38431,38432,38429,38420,38406,38386,38359,38326,38286,38238,38182,38118,38045,37962,37935,38020,38205,38476,38823,39234,39696,40198,40726,41267,41806,42329,42820,43261,43637,43930,44123,44199,44142,43937,43571,43033,42314,41410,40318,39212,38243,37382,36603,35881,35197,34528,33856,33165,32437,31660,30820,29907,28914,27836,26670,25417,24080,22669,21192,19663,18098,16516,14935,13375,11892,10526,9273,8129,7090,6152,5309,4556,3887,3296,2779,2327,1937,1601,1314,1072,867,697,556,440,345,269,208,159,121)

## increasing sample size?
e1 <- e1
e2 <- e2
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
alpha1 <- -6
alpha2 <- -0.7
alphas <- c(alpha1, alpha2)
X1 <- cbind(1, x)
eta1 <- X1%*%alphas
## component 2: exponential + cos
ome <- pi/85
beta1 <- -11
beta2 <- 0.11
beta3 <- 0.7
betas <- c(beta1, beta2, beta3)
X2 <- cbind(1, x, sin(x*ome))
eta2 <- X2%*%betas
## overall (log-)mortality
muT1 <- exp(eta1) + exp(eta2)
etaT1 <- log(muT1)

## plotting true log-mortality
## with components that we will disregard afterward
rany <- c(-16, -4)
plot(x, eta2, t="l", lwd=2, col=2, ylim=rany, axes=FALSE)
yy <- 10^seq(-7, -2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
lines(x, eta1, t="l", lwd=2, col=1)
lines(x, etaT1, col=4, lwd=2)

## population 2 (C19-2020): higher mortality and more "sinousoid"
## (higher mortality at higher ages)
alpha1 <- -5.8
alpha2 <- -0.7
alphas <- c(alpha1, alpha2)
X1 <- cbind(1, x)
eta1 <- X1%*%alphas
## component 2: exponential + cos
ome <- pi/45
beta1 <- -9.5
beta2 <- 0.105
beta3 <- 0.6
betas <- c(beta1, beta2, beta3)
X2 <- cbind(1, x, cos(x*ome))
eta2 <- X2%*%betas
## overall (log-)mortality
muT2 <- exp(eta1) + exp(eta2)
etaT2 <- log(muT2)

## plotting true log-mortality
## with components that we will disregard afterward
rany <- c(-16, -4)
plot(x, eta2, t="l", lwd=2, col=2, ylim=rany, axes=FALSE)
yy <- 10^seq(-7, -2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
lines(x, eta1, t="l", lwd=2, col=1)
lines(x, etaT2, col=4, lwd=2)



rany <- range(etaT1, etaT2)
plot(x, etaT1, t="l", lwd=2, col=2, ylim=rany, axes=TRUE)
yy <- 10^seq(-7, -2)
# axis(2, at=log(yy), labels=log(yy))
# axis(1);box()
lines(x, etaT2, col=4, lwd=2)




## simulating deaths for both populations
## as realization from a Poisson distribution
d1T <- e1*muT1
d1 <- rpois(m, d1T)
lmx1 <- log(d1/e1)
d2T <- e2*muT2
d2 <- rpois(m, d2T)
lmx2 <- log(d2/e2)

## plotting simulated and true log-mortality
rany <- c(-16, -2)
plot(x, etaT1, t="l", lwd=2, ylim=rany, axes=FALSE)
points(x, lmx1, col=2)
lines(x, etaT2, col=3, lwd=2)
points(x, lmx2, col=3, pch=2)
yy <- 10^seq(-7, -2)
axis(2, at=log(yy), labels=yy)
axis(1);box()

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

# pdf("C19agepattern.pdf", width = 10, height = 10)
par(mfrow=c(1,2))
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
plot(m0)#x, eta1.hat-eta2.hat)
# dev.off()


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

## penalty stuff
D <- diff(diag(m), diff=2)
tDD <- t(D)%*%D
lambda <- 10^2.5
P <- lambda*tDD

## PCLM independent to each population

## population 1
y <- d1g
e <- e1
C <- C1
y.st0 <- rep(y/len1, len1)
fit0 <- Mort1Dsmooth(x=x, y=y.st0,
                     offset=log(e),
                     method=3, lambda=10^3)
plot(fit0)
eta.st <- fit0$logmortality
eta <- eta.st
max.it <- 100
for(it in 1:max.it){
  gamma   <- exp(eta)
  mu      <- c(C %*% gamma)
  X       <- C * ((1 / mu) %*% t(gamma)) 
  w       <- as.vector(mu)
  r       <- y - mu + C %*% (gamma * eta)
  G       <- t(X) %*% (w * X) 
  GpP     <- G + P
  tXr     <- t(X) %*% r
  eta.old <- eta
  eta    <- solve(GpP, tXr)
  dif.eta <- max(abs((eta - eta.old)/eta.old) )
  if(dif.eta < 1e-04 & it > 4) break
  cat(it, dif.eta, "\n")
}
H0 <- solve(GpP)
Veta1 <- H0 %*% G %*% H0
se.eta1 <- sqrt(diag(Veta1))
eta1.hat <- eta
eta1.hatL <- eta1.hat - 2*se.eta1
eta1.hatU <- eta1.hat + 2*se.eta1
## rate-of-aging
Droa <- Dfun(m)
Vroa1 <- Droa %*% Veta1 %*% t(Droa)
se.roa1 <- sqrt(diag(Vroa1))
roa1.hat <- Droa%*%eta1.hat
roa1.hatL <- roa1.hat - 2*se.roa1
roa1.hatU <- roa1.hat + 2*se.roa1


## population 2
y <- d2g
e <- e2
C <- C2
y.st0 <- rep(y/len2, len2)
fit0 <- Mort1Dsmooth(x=x, y=y.st0,
                     offset=log(e),
                     method=3, lambda=10^3)
plot(fit0)
eta.st <- fit0$logmortality
eta <- eta.st
max.it <- 100
for(it in 1:max.it){
  gamma   <- exp(eta)
  mu      <- c(C %*% gamma)
  X       <- C * ((1 / mu) %*% t(gamma)) 
  w       <- as.vector(mu)
  r       <- y - mu + C %*% (gamma * eta)
  G       <- t(X) %*% (w * X) 
  GpP     <- G + P
  tXr     <- t(X) %*% r
  eta.old <- eta
  eta    <- solve(GpP, tXr)
  dif.eta <- max(abs((eta - eta.old)/eta.old) )
  if(dif.eta < 1e-04 & it > 4) break
  cat(it, dif.eta, "\n")
}
H0 <- solve(GpP)
Veta2 <- H0 %*% G %*% H0
se.eta2 <- sqrt(diag(Veta2))
eta2.hat <- eta
eta2.hatL <- eta2.hat - 2*se.eta2
eta2.hatU <- eta2.hat + 2*se.eta2
## rate-of-aging
Droa <- Dfun(m)
Vroa2 <- Droa %*% Veta2 %*% t(Droa)
se.roa2 <- sqrt(diag(Vroa2))
roa2.hat <- Droa%*%eta2.hat
roa2.hatL <- roa2.hat - 2*se.roa2
roa2.hatU <- roa2.hat + 2*se.roa2

## true roa
roaT1 <- Droa%*%etaT1
roaT2 <- Droa%*%etaT2




par(mfrow=c(2,2))
## eta1
rany <-range(etaT1, etaT2, lmx1g, lmx2g, finite=TRUE)
ranx <- range(x)
plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE, 
     xlab="age", ylab="mortality, log-scale",
     main="mortality age-pattern")
yy <- 10^seq(-7, 2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
lines(x, eta1.hat, col=2, lwd=4)
lines(x, etaT1, col=4, lwd=2, lty=2)
for(i in 1:n1){
  segments(x0=low1[i], x1=up1[i],
           y0=lmx1g[i], y1=lmx1g[i], col=4, lwd=3)
}
## eta2
plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE, 
     xlab="age", ylab="mortality, log-scale",
     main="mortality age-pattern")
yy <- 10^seq(-7, 2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
lines(x, eta2.hat, col=2, lwd=4)
lines(x, etaT2, col=4, lwd=2, lty=2)
for(i in 1:n1){
  segments(x0=low2[i], x1=up2[i],
           y0=lmx2g[i], y1=lmx2g[i], col=4, lwd=3)
}
## roa1
rany <-range(roa1.hat, roa2.hat, finite=TRUE)
ranx <- range(x)
plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE, 
     xlab="age", ylab="rate-of-aging",
     main="mortality age-pattern")
axis(2)
axis(1);box()
abline(h=0, col=8, lwd=2, lty=3)
lines(x, roa1.hat, col=2, lwd=4)
lines(x, roaT1, col=4, lwd=2, lty=2)
## roa2
plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE, 
     xlab="age", ylab="rate-of-aging",
     main="mortality age-pattern")
axis(2)
axis(1);box()
abline(h=0, col=8, lwd=2, lty=3)
lines(x, roa2.hat, col=2, lwd=4)
lines(x, roaT2, col=4, lwd=2, lty=2)
par(mfrow=c(1,1))

## both roa
rany <-range(roa1.hatL, roa1.hatU, roa2.hatL, roa2.hatU, finite=TRUE)
plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE, 
     xlab="age", ylab="rate-of-aging",
     main="mortality age-pattern")
axis(2)
axis(1);box()
abline(h=0, col=8, lwd=2, lty=3)
xx <- c(x, rev(x))
yy <- c(roa1.hatL, rev(roa1.hatU))
polygon(xx, yy, col=adjustcolor(2, 0.5), border=adjustcolor(2, 0.5))
lines(x, roa1.hat, col=2, lwd=4)

yy <- c(roa2.hatL, rev(roa2.hatU))
polygon(xx, yy, col=adjustcolor(4, 0.5), border=adjustcolor(2, 0.5))
lines(x, roa2.hat, col=4, lwd=4)


## working on the differences in roa
eta12.hat <- c(eta1.hat, eta2.hat)
U <- adiag(Droa, Droa)
H <- cbind(diag(m), -diag(m))
HU <- H%*%U
delta <- HU%*%eta12.hat
# delta1 <- roa1.hat-roa2.hat
# range(delta-delta1)

Veta12 <- adiag(Veta1, Veta2)
dim(Veta12)
Vdelta <- HU %*% Veta12 %*% t(HU)
se.delta <- sqrt(diag(Vdelta))
deltaL <- delta - 2*se.delta
deltaU <- delta + 2*se.delta

ranx <- range(x)
rany <- range(deltaL, deltaU)
plot(1,1, xlim=ranx, ylim=rany, t="n",
     xlab="ages", ylab="differences in rate-of-aging")
abline(h=0, col=8, lwd=2, lty=3)
xx <- c(x, rev(x))
yy <- c(deltaL, rev(deltaU))
polygon(xx, yy, col=adjustcolor(4, 0.5), border=adjustcolor(2, 0.5))
lines(x, delta, col=4, lwd=4)



library(mgcv)
m0 <- gam(d2 ~ s(x), offset=log(e1)+eta1.hat, family=poisson())
summary(m0)
plot(m0)




