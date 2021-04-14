## simple R-code for estimating mortality age-pattern from
## two populations with different grouping structures
## assuming a common age-pattern
## useful(?) for estimating Covid-19 mortality age-pattern
## 1) simulated data
## 2) England&Wales and USA
## by C.G. Camarda 2021.04.09


## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## R-studio to get the same dir as for the .R file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


## exposures population 1
e1 <- c(29588,29295,29002,28707,28410,28111,27809,27504,27196,26883,26566,26245,25919,25587,25251,24909,24561,24208,23849,23484,23113,22736,22353,21965,21570,21170,20800,20489,20227,20006,19819,19657,19514,19383,19258,19131,18996,18848,18681,18489,18266,18008,17711,17370,16984,16549,16065,15531,14948,14320,13648,12991,12394,11848,11340,10863,10410,9973,9547,9128,8713,8297,7878,7456,7029,6598,6164,5726,5289,4853,4422,4000,3589,3193,2815,2458,2131,1839,1579,1349,1147,969,815,681,566,468,384,313,254,204,163,129,102,80,62,48,36,28,21,15,11)
## exposures population 2
e2 <- c(38233,38253,38274,38295,38316,38337,38357,38375,38392,38406,38418,38427,38431,38432,38429,38420,38406,38386,38359,38326,38286,38238,38182,38118,38045,37962,37935,38020,38205,38476,38823,39234,39696,40198,40726,41267,41806,42329,42820,43261,43637,43930,44123,44199,44142,43937,43571,43033,42314,41410,40318,39212,38243,37382,36603,35881,35197,34528,33856,33165,32437,31660,30820,29907,28914,27836,26670,25417,24080,22669,21192,19663,18098,16516,14935,13375,11892,10526,9273,8129,7090,6152,5309,4556,3887,3296,2779,2327,1937,1601,1314,1072,867,697,556,440,345,269,208,159,121)

## increasing sample size?
e1 <- e1*100
e2 <- e2*100
## dimension & age
m <- length(e1)
x <- 1:m-1

## plotting exposures
rany <- range(e1,e2)
plot(x, e1, t="l", lwd=2, ylim=rany)
lines(x, e2, col=2, lwd=2)

## known common Covid-19 age pattern

## for illustrative purposes: 
## assuming as combination of two components (in a log-scale)
## component 1: simple exponential
alpha1 <- -11
alpha2 <- -0.7
alphas <- c(alpha1, alpha2)
X1 <- cbind(1, x)
eta1 <- X1%*%alphas
## component 2: exponential + cos
ome <- pi/45
beta1 <- -16
beta2 <- 0.11
beta3 <- 1.1
betas <- c(beta1, beta2, beta3)
X2 <- cbind(1, x, cos(x*ome))
eta2 <- X2%*%betas
## overall (log-)mortality
muT <- exp(eta1) + exp(eta2)
etaT <- log(muT)

## plotting true log-mortality
## with components that we will disregard afterward
rany <- c(-16, -4)
plot(x, eta2, t="l", lwd=2, col=2, ylim=rany, axes=FALSE)
yy <- 10^seq(-7, -2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
lines(x, eta1, t="l", lwd=2, col=1)
lines(x, etaT, col=4, lwd=2)


## simulating deaths for both populations
## as realization from a Poisson distribution
d1T <- e1*muT
d1 <- rpois(m, d1T)
lmx1 <- log(d1/e1)
d2T <- e2*muT
d2 <- rpois(m, d2T)
lmx2 <- log(d2/e2)

## plotting simulated and true log-mortality
rany <- c(-16, -2)
plot(x, etaT, t="l", lwd=2, ylim=rany, axes=FALSE)
points(x, lmx1, col=2)
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
low2 <- seq(0, 80, 20)
up2  <- c(seq(19, 79, 20), 100)
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
rany <- c(-15, -3)
plot(x, etaT, t="l", lwd=2, ylim=rany, axes=FALSE,
     xlab="age", ylab="mortality, log-scale",
     main="'typical' C19 mortality age-pattern")
points(x, lmx1, col=2)
points(x, lmx2, col=3, pch=2)
yy <- 10^seq(-7, -2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
for(i in 1:n1){
  segments(x0=low1[i], x1=up1[i],
           y0=lmx1g[i], y1=lmx1g[i], col=2, lwd=3)
}
for(i in 1:n2){
  segments(x0=low2[i], x1=up2[i],
           y0=lmx2g[i], y1=lmx2g[i], col=3, lwd=3)
}
# dev.off()


## what one observe in a single response vector
y <- c(d1g, d2g)

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
## composite matrix for the overall model
C <- rbind(C1,C2)

## penalty stuff
D <- diff(diag(m), diff=2)
tDD <- t(D)%*%D
lambda <- 10^2.5
P <- lambda*tDD

## starting eta
library(MortalitySmooth)
d1.st0 <- rep(d1g/len1, len1)
d2.st0 <- rep(d2g/len2, len2)
y.st0 <- d1.st0+d2.st0
fit0 <- Mort1Dsmooth(x=x, y=y.st0,
                     offset=log(e1+e2),
                     method=3, lambda=10^3)
plot(fit0)
eta.st <- fit0$logmortality
## PCLM regression
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
## fitted common log-mortality
eta.hat <- eta

## plotting outcomes
rany <- c(-15, -3)
plot(x, etaT, t="l", lwd=2, ylim=rany, axes=FALSE,
     xlab="age", ylab="mortality, log-scale",
     main="'typical' C19 mortality age-pattern")
yy <- 10^seq(-7, -1)
axis(2, at=log(yy), labels=yy)
axis(1);box()
for(i in 1:n1){
  segments(x0=low1[i], x1=up1[i],
           y0=lmx1g[i], y1=lmx1g[i], col=2, lwd=3)
}
for(i in 1:n2){
  segments(x0=low2[i], x1=up2[i],
           y0=lmx2g[i], y1=lmx2g[i], col=3, lwd=3)
}
lines(x, eta.hat, col=4, lwd=4)
legend("top", c("True", "Obs 1", "Obs 2", "Fitted"),
       col=1:4, lwd=c(2,2,2,4))




## actual example
## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## R-studio to get the same dir as for the .R file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## ages and dimensions
x <- 0:100
m <- length(x)


## England&Wales
low1 <- c(0, 1, seq(5,90,5))
up1  <- c(0, 4, seq(9,89,5), max(x))
## Covid-19 deaths, England & Wales (ONS), both sexes, up to 2021.04.08
d1g <- c(2,1,3,9,21,56,121,223,399,672,1324,2441,4013,6068,8193,13033,18173,24694,27148,29685)
## HMD England&Wales rounded exposures, Both sexes, 2018
e1 <- c(669211,691861,712005,715740,724496,740452,754532,749898,734785,724574,726745,708642,688934,669227,652195,638155,630687,646509,665833,690974,712259,736908,738219,754322,774901,780075,806704,818993,808466,803668,808064,795247,794789,798010,780421,782176,782238,786424,781490,758139,709801,694318,706222,721686,736031,764986,791021,819060,801836,820680,819767,830890,830177,830331,820612,807172,787402,758392,724263,713066,690461,665084,642302,619906,616814,606549,591775,594740,603224,619066,644019,674298,552941,519674,496972,467864,414459,374726,372418,364513,349290,323838,299739,274434,248360,226145,206984,186341,162814,139348,116900,97594,80572,64176,50282,38456,29222,22483,14335,7535,4528)#,3024,1970,1191,688,355)

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
## only for plotting
e1g <- G1%*%e1
lmx1g <- log(d1g/e1g)

## USA
low2 <- c(0, 1, seq(5,85,10))
up2  <- c(0, 4, seq(14,84,10), max(x))
## Covid-19 deaths, USA (CDC), both sexes, up to 2021.04.07
d2g <- c(59,31,89,804,3543,9259,25840,65781,118850,149813,165653)
## HMD USA rounded exposures, Both sexes, 2018
e2 <- c(3836098,3902806,3979332,4022696,4021466,4014217,4012505,4047202,4040051,4081791,4167920,4200339,4161575,4159962,4153749,4135745,4146522,4237521,4278578,4267685,4263407,4288529,4326594,4422873,4512135,4607788,4711052,4779568,4764408,4623439,4495740,4428088,4436002,4423408,4337589,4365936,4354958,4343036,4298810,4159753,4015540,3985238,3892901,3929940,3864760,3932324,4060009,4305776,4266404,4156279,4016044,4045375,4121559,4293401,4383064,4391719,4384033,4439260,4384204,4328658,4249479,4183061,4061669,3963409,3825247,3661756,3510499,3409261,3289099,3208093,3140174,3177042,2451204,2343223,2261360,2256172,1973215,1798484,1658210,1549383,1443026,1318675,1226448,1134483,1007730,936276,857456,784358,705371,637990,569323,476733,389846,316895,249272,188930,141920,102322,68008,44763,29692)#,18491,11272,6796,3945,2182)
n2 <- length(low2)
age.gr2 <- paste(low2,up2,sep="-")
len2 <- up2-low2+1
## matrix to group both deaths and exposures
## (grouped expoures will be used only for plotting purposes)
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
## only for plotting
e2g <- G2%*%e2
lmx2g <- log(d2g/e2g)


rany <- c(-15, -3)
plot(1,1, t="n", ylim=rany, xlim=range(x), axes=FALSE,
     xlab="age", ylab="mortality, log-scale",
     main="C19 mortality age-pattern")
yy <- 10^seq(-7, -2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
for(i in 1:n1){
  segments(x0=low1[i], x1=up1[i],
           y0=lmx1g[i], y1=lmx1g[i], col=2, lwd=3)
}
for(i in 1:n2){
  segments(x0=low2[i], x1=up2[i],
           y0=lmx2g[i], y1=lmx2g[i], col=3, lwd=3)
}
lege <- c("England & Wales","USA")
legend("top", legend=lege, col=2:3, lwd=2)


## what one observe in a single response vector
y <- c(d1g, d2g)

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
## composite matrix for the overall model
C <- rbind(C1,C2)

## penalty stuff
D <- diff(diag(m), diff=2)
tDD <- t(D)%*%D
lambda <- 10^4
P <- lambda*tDD

## starting eta
library(MortalitySmooth)
d1.st0 <- rep(d1g/len1, len1)
d2.st0 <- rep(d2g/len2, len2)
y.st0 <- d1.st0+d2.st0
fit0 <- Mort1Dsmooth(x=x, y=y.st0,
                     offset=log(e1+e2),
                     method=3, lambda=10^4)
plot(fit0)
eta.st <- fit0$logmortality
## PCLM regression
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
## fitted common log-mortality
eta.hat <- eta

## plotting outcomes
rany <- c(-15, -3)
plot(1,1, t="n", ylim=rany, xlim=range(x), axes=FALSE,
     xlab="age", ylab="mortality, log-scale",
     main="C19 mortality age-pattern")
yy <- 10^seq(-7, -2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
for(i in 1:n1){
  segments(x0=low1[i], x1=up1[i],
           y0=lmx1g[i], y1=lmx1g[i], col=2, lwd=3)
}
for(i in 1:n2){
  segments(x0=low2[i], x1=up2[i],
           y0=lmx2g[i], y1=lmx2g[i], col=3, lwd=3)
}
lines(x, eta.hat, col=4, lwd=4)
lege <- c("England & Wales","USA","Fitted common pattern")
legend("top", legend=lege, col=2:4, lwd=2)























# 
# ## actual patterns from England&Wales and USA
# rm(list = ls())
# options(device="X11")
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# 
# ## ENG
# pdf("C19agepatternENWandUSA.pdf", width = 10, height = 10)
# low <- c(0, 1, seq(5,90,5))
# up  <- c(1, 4, seq(9,89,5), 100)
# cbind(low, up)
# d <- c(2,1,3,9,21,56,121,223,399,672,1324,2441,4013,6068,8193,13033,18173,24694,27148,29685)
# e <- c(669797,2845633,3708320,3450782,3270795,3717960,4022272,3976030,3900901,3559955,4005397,4137131,3785564,3234026,3006776,2900152,1985125,1491797,918437,528959)
# lmx <- log(d/e)
# rany <- c(-16, -2)
# plot(1,1,t="n", ylim=rany, xlim=c(0,100), axes=FALSE,
#      xlab="age", ylab="mortality, log-scale",
#      main="C19 mortality (C19 deaths / Population)")
# yy <- 10^seq(-7, -1)
# axis(2, at=log(yy), labels=yy)
# axis(1);box()
# segments(x0=0, x1=1,
#          y0=lmx[1], y1=lmx[1], col=1, lwd=3)
# for(i in 2:length(low)){
#   segments(x0=low[i], x1=up[i]+1,
#            y0=lmx[i], y1=lmx[i], col=1, lwd=3)
# }
# 
# 
# ## USA  
# low <- c(0, 1, seq(5,85,10))
# up  <- c(1, 4, seq(14,84,10), 100)
# cbind(low, up)
# d <- c(59,31,89,804,3543,9259,25840,65781,118850,149813,165653)
# e <- c(3745000,15713000,42255000,43166000,45866000,42932000,41356000,43259000,32456000,15786000,6063000)
# lmx <- log(d/e)
# segments(x0=0, x1=1,
#          y0=lmx[1], y1=lmx[1], col=2, lwd=3)
# for(i in 2:length(low)){
#   segments(x0=low[i], x1=up[i]+1,
#            y0=lmx[i], y1=lmx[i], col=2, lwd=3)
# }
# 
# lege <- c("England & Wales (ONS), up to 2021.04.08",
#           "USA (CDC), up to 2021.04.07")
# legend("top", legend=lege, col=1:2, lwd=2)
# dev.off()
## END