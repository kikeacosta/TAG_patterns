## simple R-code for estimating log-mortality age-pattern 
## in a population for each sex in two different periods, 
## e.g. pre-pandemic years and 2020 (or 2021)

## notation:
## pre-pandemic year(s) is(are) denoted by 1
## 2020 (or 2021) is denoted by 2
## Males and Females are denoted by M and F
## x is the vector of ages
## eta denoted log-mortality

## eta.F1 : log-mortality for females in pre-pandemic years
## eta.F2 :  log-mortality for females in 2020
## eta.M1 : log-mortality for males in pre-pandemic years
## eta.M2 :  log-mortality for males in 2020

## assumptions:

## eta.F1 : smooth function over age/x

## eta.F2 = eta.F1 + c.F + delta.F
## where 
## c.F is a "pandemic" age-independent scaling factor for females
## delta.F is a smooth perturbation function over age/x for females

## eta.M1 = eta.F1 + s
## where 
## s is a sex age-factor, smooth function over age

## eta.M2 = eta.M1 + c.M + delta.M
## where 
## c.M is a "pandemic" age-independent scaling factor for males
## delta.M is a smooth perturbation function over age/x for males

## eta.M2 = eta.F1 + s + c.M + delta.M

## SIMULATED DATA

## ISSUES:
## 1) underdetermined linear system of equation
## 2) smoothing parameters selection
## by C.G. Camarda 2022.04.01


## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## R-studio to get the same dir as for the .R file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(MortalitySmooth)
library(magic)

## ages and dimensions
x <- 0:100
m <- length(x)

## exposures, given information over age
## assuming a much larger population in pre-pandemic years
## which is typical if more years are merged
e.F1 <- c(1789868,1826814,1859435,1885465,1904620,1918796,1930258,1939072,1945262,1948767,1949568,1947748,1943060,1935261,1924527,1911532,1896861,1881408,1866384,1852991,1841847,1833282,1829217,1831976,1842624,1860692,1886240,1916500,1947086,1974740,1999588,2021911,2039110,2047957,2047197,2036804,2014789,1990193,1977174,1984198,2010024,2056672,2114331,2166637,2203154,2225496,2232360,2228487,2222699,2220467,2220223,2222260,2223887,2220424,2209278,2191546,2167114,2139396,2113708,2093292,2076664,2062442,2051143,2043461,2038760,2037559,2043348,2042208,2012871,1943753,1837646,1691219,1527352,1382142,1277859,1210271,1181200,1178486,1180563,1173168,1159113,1138281,1109330,1072315,1027341,973119,908202,835426,759641,683533,606385,532710,458909,382843,304263,228977,166113,117157,80185,47820,83710)
e.M1 <- c(1876826,1912217,1944204,1970433,1990359,2005474,2018156,2028349,2035941,2040722,2042461,2041013,2036171,2027809,2016059,2001316,1983998,1964543,1943638,1922281,1901326,1881264,1864381,1853546,1850360,1854687,1866923,1884846,1904609,1923279,1941058,1958316,1972416,1979947,1979440,1970695,1951498,1930647,1922036,1933833,1964556,2016055,2077832,2132653,2169615,2190410,2193676,2184652,2172960,2164564,2157807,2153029,2147483,2136320,2116872,2090339,2056639,2019352,1984045,1954132,1928245,1905228,1885363,1868867,1854844,1843638,1838076,1825899,1788623,1716292,1611579,1472005,1316942,1176631,1069529,991845,944904,918627,895409,863932,827016,784968,737488,685716,630747,572201,509575,445283,383198,325435,271778,223983,180741,141045,104891,73991,50339,33010,20717,10756,14351)
e.F2 <- c(337152,346477,354847,362282,368806,374441,379210,383136,386240,388546,390156,391172,391219,390159,388255,386100,383810,380940,377405,373588,370014,366648,364281,363384,363740,364297,364828,367634,373614,381564,389306,397188,403819,408259,411012,414209,418243,418970,414804,407939,401362,393521,391735,399644,413495,426186,439484,448645,450701,448115,446203,443769,442019,442045,442934,442933,442689,440701,436144,429980,423837,417194,411993,409300,408133,405608,401303,398814,399318,400644,401182,403347,394101,367464,330326,294614,256807,229513,219851,221498,220856,220457,218094,211244,201329,192670,184856,174274,159702,142581,124143,109772,97460,80441,58716,45370,40645,33603,24243,12565,16450)
e.M2 <- c(354172,362805,370771,378050,384625,390480,395596,399956,403543,406339,408423,409874,410192,409167,407071,404592,401838,398160,393385,387971,382687,377547,372963,369273,366482,363933,361366,361305,364818,370651,376376,382345,387623,391452,394153,397431,401698,402942,399576,393705,388175,381408,380637,389425,403995,417331,431228,440480,441913,438187,435031,431175,427925,426465,425835,424156,422127,418209,411550,403167,394765,385810,378355,373536,370325,365853,359839,355303,353185,351671,349436,348455,337882,312728,278695,245784,211174,184696,171991,168160,162671,157641,151029,140764,128242,117170,107091,96187,84016,71267,57512,46539,39580,30734,20002,13679,11798,9394,6468,3020,2993)

## plotting exposures
rany <- range(e.F1,e.F2,e.M1,e.M2)
plot(x, e.F1, t="l", col=2, lwd=2, ylim=rany)
lines(x, e.F2, col=2, lwd=2, lty=2)
lines(x, e.M1, t="l", col=4, lwd=2, ylim=rany)
lines(x, e.M2, col=4, lwd=2, lty=2)

## hypothesized true log-mortality for each sex and period

## eta.F1 

## combination of two components (in a log-scale)
## component 1: simple exponential (~Gompertz)
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
## true (log-)mortality for females in year(s) 1
mu.F1T <- exp(eta1) + exp(eta2)
eta.F1T <- log(mu.F1T)

## plotting true log-mortality
## with components that we will disregard afterward
rany <- c(-12, 0)
plot(x, eta2, t="l", lwd=2, col=2, ylim=rany, axes=FALSE)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy)
axis(1);box()
lines(x, eta1, t="l", lwd=2, col=1)
lines(x, eta.F1T, col=2, lwd=2)

## males in year(s) 1: eta.M1 = eta.F1 + s
## constructing true s vector
s.T <-  dgamma(x, shape=20, scale=1.5)*30 + 0.5
plot(x, s.T)
## true (log-)mortality for the males in year 1
eta.M1T <- eta.F1T + s.T
mu.M1T <- exp(eta.M1T)

## plot males+females in year(s) 1
rany <- range(eta.F1T, eta.M1T)
plot(x, eta.F1T, t="l", lwd=2, col=2, ylim=rany, axes=FALSE)
yy <- 10^seq(-20, 2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
lines(x, eta.M1T, col=4, lwd=2)


## females in year 2: ## eta.F2 = eta.F1 + c.F + delta.F
## true scaling factor
c.FT <- 0.4
## constructing true delta vector
ome <- pi/40
beta1 <- 0.3
beta2 <- 0.02
beta3 <- 0.5
betas <- c(beta1, beta2, beta3)
X2 <- cbind(1, x, cos(x*ome))
delta.FT <- X2%*%betas
delta.FT <-  delta.FT - mean(delta.FT)
plot(x, delta.FT)
## true log-mortality for the females in year(s) 2
eta.F2T <- eta.F1T + c.FT + delta.FT
mu.F2T <- exp(eta.F2T)

## males in year 2: eta.M2 = eta.M1 + c.M + delta.M
## true scaling factor
c.MT <- 0.5
## constructing true delta vector
ome <- pi/40
beta1 <- 0.3
beta2 <- 0.02
beta3 <- -0.5
betas <- c(beta1, beta2, beta3)
X2 <- cbind(1, x, cos(x*ome))
delta.MT <- X2%*%betas
delta.MT <-  delta.MT - mean(delta.MT)
plot(x, delta.MT)
lines(x, delta.MT)
## true log-mortality for the second population (C19-2020)
eta.M2T <- eta.M1T + c.MT + delta.MT
mu.M2T <- exp(eta.M2T)
## check if eta.M2 = eta.M1 + c.M + delta.M
## is equal to eta.M2 = eta.F1 + s + c.M + delta.M
range(eta.M2T  - (eta.F1T + s.T + c.MT + delta.MT))


## plotting all 4 true log-mortality patterns
rany <- range(eta.F1T, eta.M1T, eta.F2T, eta.M2T)
plot(x, eta.F1T, t="l", lwd=2, col=2, ylim=rany, axes=FALSE)
yy <- 10^seq(-20, 2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
lines(x, eta.M1T, col=4, lwd=2)
lines(x, eta.F2T, col=2, lwd=2, lty=2)
lines(x, eta.M2T, col=4, lwd=2, lty=2)


## simulating deaths for both populations
## as realization from a Poisson distribution
d.F1T <- e.F1*mu.F1T
d.M1T <- e.M1*mu.M1T
d.F2T <- e.F2*mu.F2T
d.M2T <- e.M2*mu.M2T
## simulating deaths
d.F1 <- rpois(m, d.F1T)
d.M1 <- rpois(m, d.M1T)
d.F2 <- rpois(m, d.F2T)
d.M2 <- rpois(m, d.M2T)
## log-rates
lmx.F1 <- log(d.F1/e.F1)
lmx.M1 <- log(d.M1/e.M1)
lmx.F2 <- log(d.F2/e.F2)
lmx.M2 <- log(d.M2/e.M2)


## plotting simulated and true log-mortality
rany <- range(eta.F1T, eta.M1T, eta.F2T, eta.M2T,
              lmx.F1, lmx.M1, lmx.F2, lmx.M2)
plot(x, eta.F1T, t="l", lwd=2, col=2, ylim=rany, axes=FALSE)
yy <- 10^seq(-20, 2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
lines(x, eta.M1T, col=4, lwd=2)
lines(x, eta.F2T, col=2, lwd=2, lty=2)
lines(x, eta.M2T, col=4, lwd=2, lty=2)
points(x, lmx.F1, col=2, pch=16)
points(x, lmx.M1, col=4, pch=16)
points(x, lmx.F2, col=2, pch=17)
points(x, lmx.M2, col=4, pch=17)

## assuming different grouping structure

## grouping for year(s) 1
low1 <- c(0, 1, seq(5, 90,5))
up1  <- c(0, seq(4, 89,5), 100)
n1 <- length(low1)
age.gr1 <- paste(low1,up1,sep="-")
len1 <- up1-low1+1
## matrix to group both deaths and exposures
## (grouped exposures will be used only for plotting purposes)
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

## what we would observed in an actual world for year(s) 1
d.F1g <- G1%*%d.F1
e.F1g <- G1%*%e.F1
lmx.F1g <- log(d.F1g/e.F1g)
d.M1g <- G1%*%d.M1
e.M1g <- G1%*%e.M1
lmx.M1g <- log(d.M1g/e.M1g)



## grouping for year 2
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
## what we would observed in an actual world for year 2
d.F2g <- G2%*%d.F2
e.F2g <- G2%*%e.F2
lmx.F2g <- log(d.F2g/e.F2g)
d.M2g <- G2%*%d.M2
e.M2g <- G2%*%e.M2
lmx.M2g <- log(d.M2g/e.M2g)

## plotting what we actually observed
rany <- range(eta.F1T, eta.M1T, eta.F2T, eta.M2T,
              lmx.F1, lmx.M1, lmx.F2, lmx.M2)
plot(x, eta.F1T, t="l", lwd=2, col=2, ylim=rany, axes=FALSE)
yy <- 10^seq(-20, 2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
lines(x, eta.M1T, col=4, lwd=2)
lines(x, eta.F2T, col=2, lwd=2, lty=2)
lines(x, eta.M2T, col=4, lwd=2, lty=2)
# points(x, lmx.F1, col=2, pch=16)
# points(x, lmx.M1, col=4, pch=16)
# points(x, lmx.F2, col=2, pch=17)
# points(x, lmx.M2, col=4, pch=17)
for(i in 1:n1){
  segments(x0=low1[i], x1=up1[i],
           y0=lmx.F1g[i], y1=lmx.F1g[i], col=2, lwd=3)
  segments(x0=low1[i], x1=up1[i],
           y0=lmx.M1g[i], y1=lmx.M1g[i], col=4, lwd=3)
}
for(i in 1:n2){
  segments(x0=low2[i], x1=up2[i],
           y0=lmx.F2g[i], y1=lmx.F2g[i], col=2, lwd=3)
  segments(x0=low2[i], x1=up2[i],
           y0=lmx.M2g[i], y1=lmx.M2g[i], col=4, lwd=3)
}


## build composite matrices for both periods
C.F1 <- C.M1 <- G1
for(i in 1:n1){
  age.low.i <- low1[i]
  age.up.i <- up1[i]
  whi <- which(x>=age.low.i & x<=age.up.i)
  C.F1[i,whi] <- e.F1[whi]
  C.M1[i,whi] <- e.M1[whi]
}
C.F2 <- C.M2 <- G2
for(i in 1:n2){
  age.low.i <- low2[i]
  age.up.i <- up2[i]
  whi <- which(x>=age.low.i & x<=age.up.i)
  C.F2[i,whi] <- e.F2[whi]
  C.M2[i,whi] <- e.M2[whi]
}


## estimation
y <- c(d.F1g, d.M1g, d.F2g, d.M2g)
C <- adiag(C.F1, C.M1, C.F2, C.M2)
y.F1.st0 <- rep(d.F1g/len1, len1)
fit.F1.0 <- Mort1Dsmooth(x=x, y=y.F1.st0,
                         offset=log(e.F1),
                         method=3, lambda=10^4)
plot(fit.F1.0)
y.M1.st0 <- rep(d.M1g/len1, len1)
fit.M1.0 <- Mort1Dsmooth(x=x, y=y.M1.st0,
                        offset=log(e.M1),
                        method=3, lambda=10^4)
plot(fit.M1.0)
y.F2.st0 <- rep(d.F2g/len2, len2)
fit.F2.0 <- Mort1Dsmooth(x=x, y=y.F2.st0,
                        offset=log(e.F2),
                        method=3, lambda=10^4)
plot(fit.F2.0)
y.M2.st0 <- rep(d.M2g/len2, len2)
fit.M2.0 <- Mort1Dsmooth(x=x, y=y.M2.st0,
                         offset=log(e.M2),
                         method=3, lambda=10^4)
plot(fit.M2.0)



eta.st <- c(fit.F1.0$logmortality, fit.M1.0$logmortality,
            fit.F1.0$logmortality, fit.M2.0$logmortality)
## penalty stuff
Deta.F1 <- diff(diag(m), diff=2)
tDDeta.F1 <- t(Deta.F1)%*%Deta.F1
Ds <- diff(diag(m), diff=2)
tDDs <- t(Ds)%*%Ds
Ddelta.F <- diff(diag(m), diff=2)
tDDdelta.F <- t(Ddelta.F)%*%Ddelta.F
Ddelta.M <- diff(diag(m), diff=2)
tDDdelta.M <- t(Ddelta.M)%*%Ddelta.M

## model matrix
Ueta.F1 <- kronecker(rep(1,4), diag(m))
Us <- rbind(0*diag(m), diag(m), 0*diag(m), diag(m))
Uc.F <- matrix(c(rep(0,2*m), rep(1,m), rep(0,m)), 4*m)
Udelta.F <- rbind(0*diag(m), 0*diag(m), diag(m), 0*diag(m))
Uc.M <- matrix(c(rep(0,3*m), rep(1,m)), 4*m)
Udelta.M <- rbind(0*diag(m), 0*diag(m), 0*diag(m), diag(m))
U <- cbind(Ueta.F1, Us, Uc.F, Udelta.F, Uc.M, Udelta.M)

## ridge penalty
Pr <- 1e-4*diag(ncol(U))


## constraining delta to sum up to 0
H <- rbind(c(rep(0,2*m+1), rep(1,m), rep(0,m+1)),
           c(rep(0,2*m+1), rep(0,m), 0, rep(1,m)))
## check
H%*%c(eta.F1T, s.T, c.FT, delta.FT, c.MT, delta.MT)
kappa <- c(0,0)

## try three lambdas:
## common lambdas for both lambda.delta.F and lambda.delta.M
lambdas.eta.F1 <- 10^seq(1, 5, 1)
nl.eta <- length(lambdas.eta.F1)
lambdas.s <- 10^seq(1, 5, 1)
nl.s <- length(lambdas.s)
lambdas.delta <- 10^seq(2, 5, 1)
nl.delta <- length(lambdas.delta)

BICs <- array(NA, dim=c(nl.eta, nl.s, nl.delta), 
              dimnames = list(log10(lambdas.eta.F1), 
                              log10(lambdas.s), 
                              log10(lambdas.delta)))
tim1 <- Sys.time()
l.eta <- l.s <- l.delta <- 1
for(l.eta in 1:nl.eta){
  for(l.s in 1:nl.s){
    for(l.delta in 1:nl.delta){
      lambda.eta.F1 <- lambdas.eta.F1[l.eta]
      lambda.s <- lambdas.s[l.s]
      lambda.delta.F <- lambda.delta.M <- lambdas.delta[l.delta]
  
      ## penalty term
      P <- adiag(lambda.eta.F1*tDDeta.F1,
                 lambda.s*tDDs,
                 0,
                 lambda.delta.F*tDDdelta.F,
                 0,
                 lambda.delta.M*tDDdelta.M)
      
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
        ## cat(it, dif.eta, "\n")
      }
      ## diagnostics
      Hat <- solve(GpP, G)
      ED <- sum(diag(Hat))
      y1 <- y
      y1[y == 0] <- 10^(-4)
      DEV <- 2 * sum( y1 * log(y1/mu) )
      AIC <- DEV + 2 * ED
      BIC <- DEV + log(length(y)) * ED
      psi <- DEV/(length(y)-ED)
      QIC <- length(y) + ED + length(y)*log(psi)
      
      if(it<100){
        BICs[l.eta, l.s, l.delta] <- BIC
      }
      cat(log10(lambda.eta.F1), 
          log10(lambda.s),
          log10(lambda.delta.F), it, dif.eta, "\n")
    }
  }
}
tim2 <- Sys.time()
tim2-tim1
(whimin <- which(BICs==min(BICs, na.rm=TRUE), arr.ind=TRUE))
nl.eta;nl.s;nl.delta

lambda.eta.F1 <- lambdas.eta.F1[whimin[1]]
lambda.s <- lambdas.s[whimin[2]]
lambda.delta.F <- lambda.delta.M <- lambdas.delta[whimin[3]]


## penalty term
P <- adiag(lambda.eta.F1*tDDeta.F1,
           lambda.s*tDDs,
           0,
           lambda.delta.F*tDDdelta.F,
           0,
           lambda.delta.M*tDDdelta.M)


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

eta.F1.hat <- betas[1:m]
s.hat <- betas[1:m+m]
c.F.hat <- betas[2*m+1]
delta.F.hat <- betas[1:m+2*m+1]
c.M.hat <- betas[3*m+2]
delta.M.hat <- betas[1:m+3*m+2]


eta.M1.hat <- eta.F1.hat + s.hat
eta.F2.hat <- eta.F1.hat + c.F.hat + delta.F.hat
eta.M2.hat <- eta.F1.hat + s.hat + c.M.hat + delta.M.hat # 
## check
range(eta.M2.hat - (eta.M1.hat + c.M.hat + delta.M.hat))

## check constraints
H%*%betas
sum(delta.F.hat);sum(delta.M.hat)

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
se.c.F <- se.betas[2*m+1]
se.delta.F <- se.betas[1:m+1+2*m]
se.c.M <- se.betas[3*m+2]
se.delta.M <- se.betas[1:m+3*m+2]

V.eta <- U %*% Vbetas %*% t(U)
se.eta <- sqrt(diag(V.eta))
se.eta.F1 <- se.eta[1:m]
se.eta.M1 <- se.eta[1:m+m]
se.eta.F2 <- se.eta[1:m+2*m]
se.eta.M2 <- se.eta[1:m+3*m]

eta.F1.hatL <- eta.F1.hat - 2*se.eta.F1
eta.F1.hatU <- eta.F1.hat + 2*se.eta.F1
eta.M1.hatL <- eta.M1.hat - 2*se.eta.M1
eta.M1.hatU <- eta.M1.hat + 2*se.eta.M1

eta.F2.hatL <- eta.F2.hat - 2*se.eta.F2
eta.F2.hatU <- eta.F2.hat + 2*se.eta.F2
eta.M2.hatL <- eta.M2.hat - 2*se.eta.M2
eta.M2.hatU <- eta.M2.hat + 2*se.eta.M2

s.hatL <- s.hat - 2*se.s
s.hatU <- s.hat + 2*se.s
c.F.hatL <- c.F.hat - 2*se.c.F
c.F.hatU <- c.F.hat + 2*se.c.F
delta.F.hatL <- delta.F.hat - 2*se.delta.F
delta.F.hatU <- delta.F.hat + 2*se.delta.F
c.M.hatL <- c.M.hat - 2*se.c.M
c.M.hatU <- c.M.hat + 2*se.c.M
delta.M.hatL <- delta.M.hat - 2*se.delta.M
delta.M.hatU <- delta.M.hat + 2*se.delta.M


par(mfrow=c(3,2))
## year(s) 1
rany <-range(eta.F1T, eta.M1T, eta.F2T, eta.M2T, 
             lmx.F1g, lmx.M1g, lmx.F2g, lmx.M2g, finite=TRUE)
ranx <- range(x)
plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE, 
     xlab="age", ylab="mortality, log-scale",
     main="mortality age-pattern")
yy <- 10^seq(-7, 2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
xx <- c(x, rev(x))
yy <- c(eta.F1.hatL, rev(eta.F1.hatU))
polygon(xx, yy, col=adjustcolor(2, 0.5), border=adjustcolor(2, 0.5))
lines(x, eta.F1.hat, col=2, lwd=4)
lines(x, eta.F1T, col=4, lwd=2, lty=2)
for(i in 1:n1){
  segments(x0=low1[i], x1=up1[i],
           y0=lmx.F1g[i], y1=lmx.F1g[i], col=3, lwd=3)
}
yy <- c(eta.M1.hatL, rev(eta.M1.hatU))
polygon(xx, yy, col=adjustcolor(2, 0.5), border=adjustcolor(2, 0.5))
lines(x, eta.M1.hat, col=2, lwd=4)
lines(x, eta.M1T, col=4, lwd=2, lty=2)
for(i in 1:n1){
  segments(x0=low1[i], x1=up1[i],
           y0=lmx.M1g[i], y1=lmx.M1g[i], col=3, lwd=3)
}
## year 2
rany <-range(eta.F1T, eta.M1T, eta.F2T, eta.M2T, 
             lmx.F1g, lmx.M1g, lmx.F2g, lmx.M2g, finite=TRUE)
ranx <- range(x)
plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE, 
     xlab="age", ylab="mortality, log-scale",
     main="mortality age-pattern")
yy <- 10^seq(-7, 2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
xx <- c(x, rev(x))
yy <- c(eta.F2.hatL, rev(eta.F2.hatU))
polygon(xx, yy, col=adjustcolor(2, 0.5), border=adjustcolor(2, 0.5))
lines(x, eta.F2.hat, col=2, lwd=4)
lines(x, eta.F2T, col=4, lwd=2, lty=2)
for(i in 1:n1){
  segments(x0=low2[i], x1=up2[i],
           y0=lmx.F2g[i], y1=lmx.F2g[i], col=3, lwd=3)
}
yy <- c(eta.M2.hatL, rev(eta.M2.hatU))
polygon(xx, yy, col=adjustcolor(2, 0.5), border=adjustcolor(2, 0.5))
lines(x, eta.M2.hat, col=2, lwd=4)
lines(x, eta.M2T, col=4, lwd=2, lty=2)
for(i in 1:n1){
  segments(x0=low2[i], x1=up2[i],
           y0=lmx.M2g[i], y1=lmx.M2g[i], col=3, lwd=3)
}
## s
ranx <- range(x)
rany <- range(s.hatL, s.hatU)
plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE,
     xlab="ages", ylab="s")
axis(2)
axis(1);box()
abline(h=0, col=8, lty=3, lwd=2)
xx <- c(x, rev(x))
yy <- c(s.hatL, rev(s.hatU))
polygon(xx, yy, col=adjustcolor(2, 0.5), border=adjustcolor(4, 0.5))
lines(x, s.hat, col=2, lwd=4)
lines(x, s.T, col=4, lwd=2, lty=2)
## delta.F
ranx <- range(x)
rany <- range(delta.F.hatL, delta.F.hatU)
plot(1, 1, t="n", xlim=ranx, ylim=rany,
     main=paste("c = ", signif(c.F.hat,4), "+ -", signif(2*se.c.F,4)), axes=FALSE,
     xlab="ages", ylab="delta")
axis(2)
axis(1);box()
abline(h=0, col=8, lty=3, lwd=2)
xx <- c(x, rev(x))
yy <- c(delta.F.hatL, rev(delta.F.hatU))
polygon(xx, yy, col=adjustcolor(2, 0.5), border=adjustcolor(4, 0.5))
lines(x, delta.F.hat, col=2, lwd=4)
lines(x, delta.FT, col=4, lwd=2, lty=2)
## delta.M
ranx <- range(x)
rany <- range(delta.M.hatL, delta.M.hatU)
plot(1, 1, t="n", xlim=ranx, ylim=rany,
     main=paste("c = ", signif(c.M.hat,4), "+ -", signif(2*se.c.M,4)), axes=FALSE,
     xlab="ages", ylab="delta")
axis(2)
axis(1);box()
abline(h=0, col=8, lty=3, lwd=2)
xx <- c(x, rev(x))
yy <- c(delta.M.hatL, rev(delta.M.hatU))
polygon(xx, yy, col=adjustcolor(2, 0.5), border=adjustcolor(4, 0.5))
lines(x, delta.M.hat, col=2, lwd=4)
lines(x, delta.MT, col=4, lwd=2, lty=2)
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