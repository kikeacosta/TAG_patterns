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

## by C.G. Camarda 2022.06.01

## ACTUAL DATA


## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## !!!! to be changed
setwd("~/WORK/TAG_patterns/")
library(MortalitySmooth)
library(magic)
library(colorspace)
library(ggplot2)
library(plotly)


## DATA ISSUES 
# - no deaths above 85 in Tunisia in pre-pandemic years
# - need to remove deaths above age 100 for certain (HMD) countries
# - Bolivia at oldest ages
# - age 100 in Spain and Japan
# - single death in Aruba at age 15-19
# - Oman males 70-74, especially in 2020
## MODEL ISSUES
# - some issues with CI for deltas in US (practically solved)
# - I may think to impose monotonicity in eta after a certain age (implemented)
# - maybe it would be better to start after childhood (10?)
# - undersmooth for eta in certain population at oldest ages?


## loading population
offset <- read.csv("Output/offsets.csv", header=TRUE)
## loading deaths
deaths <- read.csv("Output/annual_deaths_countries_selected_sources.csv", header=TRUE)

## remove deaths above 100 (mismatch with exposures)
deaths <- subset(deaths, Age<=100)

## ages and dimensions
x <- 0:100
m <- length(x)
## pre-pandemic years
t1 <- 2015:2019
n1 <- length(t1)
## pandemic years, independently analyzed
t2 <- 2020 ## 2021

## common countries between deaths and offset
codes.offset <- offset$Code
codes.deaths <- deaths$Code
codes <- intersect(codes.offset,codes.deaths)
nc <- length(codes)
nc

OUTPUT <- list()
j=1

ds <- rep(0,nc)
ds[c(1,50)] <- 0.8
tim1 <- Sys.time()
for(j in 1:nc){
  ## for a given population
  code <- codes[j]
  ## exposures in pre-pandemic and pandemic years
  offset.i <- subset(offset, Code==code)
  offset.i.1 <- subset(offset.i, Year<=2019)
  offset.i.2 <- subset(offset.i, Year==2020)
  ## summing up over periods and sex
  offset.i.1F <- subset(offset.i.1, Sex=="f")
  offset.i.1M <- subset(offset.i.1, Sex=="m")
  offset.i.2F <- subset(offset.i.2, Sex=="f")
  offset.i.2M <- subset(offset.i.2, Sex=="m")
  e.F1 <- c(tapply(offset.i.1F$Population, offset.i.1F$Age, sum))
  e.M1 <- c(tapply(offset.i.1M$Population, offset.i.1M$Age, sum))
  e.F2 <- c(tapply(offset.i.2F$Population, offset.i.2F$Age, sum))
  e.M2 <- c(tapply(offset.i.2M$Population, offset.i.2M$Age, sum))
  
  # ## plotting exposures
  # rany <- range(e.F1,e.F2,e.M1,e.M2)
  # plot(x, e.F1, t="l", col=2, lwd=2, ylim=rany)
  # lines(x, e.F2, col=2, lwd=2, lty=2)
  # lines(x, e.M1, t="l", col=4, lwd=2, ylim=rany)
  # lines(x, e.M2, col=4, lwd=2, lty=2)
  
  
  ## deaths in pre-pandemic and pandemic years
  deaths.i <- subset(deaths, Code==code)
  deaths.i.1 <- subset(deaths.i, Year<=2019)
  deaths.i.2 <- subset(deaths.i, Year==2020)
  ## summing up over periods and sex
  deaths.i.1F <- subset(deaths.i.1, Sex=="f")
  deaths.i.1M <- subset(deaths.i.1, Sex=="m")
  deaths.i.2F <- subset(deaths.i.2, Sex=="f")
  deaths.i.2M <- subset(deaths.i.2, Sex=="m")
  d.F1g <- tapply(deaths.i.1F$Deaths, deaths.i.1F$Age, sum)
  d.M1g <- tapply(deaths.i.1M$Deaths, deaths.i.1M$Age, sum)
  d.F2g <- tapply(deaths.i.2F$Deaths, deaths.i.2F$Age, sum)
  d.M2g <- tapply(deaths.i.2M$Deaths, deaths.i.2M$Age, sum)
  
  ## create composite, eventually for both periods
  ## grouping for pre-pandemic years
  low1 <- unique(deaths.i.1F$Age)
  up1 <- c(low1[-1]-1, max(x))
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
  ## aggregating exposures only for plotting 
  e.F1g <- c(G1%*%e.F1)
  e.M1g <- c(G1%*%e.M1)
  lmx.F1g <- log(d.F1g/e.F1g)
  lmx.M1g <- log(d.M1g/e.M1g)
  
  ## grouping for pandemic year 
  low2 <- unique(deaths.i.2F$Age)
  up2 <- c(low2[-1]-1, max(x))
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
  ## aggregating exposures only for plotting 
  e.F2g <- c(G2%*%e.F2)
  e.M2g <- c(G2%*%e.M2)
  lmx.F2g <- log(d.F2g/e.F2g)
  lmx.M2g <- log(d.M2g/e.M2g)
  
  
  
  ## plotting what we actually observed
  rany <- range(lmx.F1g, lmx.M1g, lmx.F2g, lmx.M2g, finite=TRUE)
  ranx <- range(x)
  plot(1, 1, t="l", lwd=2, col=2, ylim=rany, xlim=ranx, axes=FALSE)
  title(main=paste(code))
  yy <- 10^seq(-20, 2)
  axis(2, at=log(yy), labels=yy)
  axis(1);box()
  for(i in 1:n1){
    segments(x0=low1[i], x1=up1[i],
             y0=lmx.F1g[i], y1=lmx.F1g[i], col=2, lwd=3)
    segments(x0=low1[i], x1=up1[i],
             y0=lmx.M1g[i], y1=lmx.M1g[i], col=4, lwd=3)
  }
  for(i in 1:n2){
    segments(x0=low2[i], x1=up2[i],
             y0=lmx.F2g[i], y1=lmx.F2g[i], col=2, lwd=3, lty=2)
    segments(x0=low2[i], x1=up2[i],
             y0=lmx.M2g[i], y1=lmx.M2g[i], col=4, lwd=3, lty=2)
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
  wei.F1 <- rep(1,m)
  wei.F1[e.F1==0] <- 0
  fit.F1.0 <- Mort1Dsmooth(x=x, y=y.F1.st0,
                           offset=log(e.F1),
                           w=wei.F1,
                           method=3, lambda=10^4)
  # plot(fit.F1.0)
  y.M1.st0 <- rep(d.M1g/len1, len1)
  wei.M1 <- rep(1,m)
  wei.M1[e.M1==0] <- 0
  fit.M1.0 <- Mort1Dsmooth(x=x, y=y.M1.st0,
                           offset=log(e.M1),
                           w=wei.M1,
                           method=3, lambda=10^4)
  # plot(fit.M1.0)
  y.F2.st0 <- rep(d.F2g/len2, len2)
  wei.F2 <- rep(1,m)
  wei.F2[e.F2==0] <- 0
  fit.F2.0 <- Mort1Dsmooth(x=x, y=y.F2.st0,
                           offset=log(e.F2),
                           w=wei.F2,
                           method=3, lambda=10^4)
  # plot(fit.F2.0)
  y.M2.st0 <- rep(d.M2g/len2, len2)
  wei.M2 <- rep(1,m)
  wei.M2[e.M2==0] <- 0
  fit.M2.0 <- Mort1Dsmooth(x=x, y=y.M2.st0,
                           offset=log(e.M2),
                           w=wei.M2,
                           method=3, lambda=10^4)
  # plot(fit.M2.0)
  
  eta.st <- c(fit.F1.0$logmortality, fit.M1.0$logmortality,
              fit.F1.0$logmortality, fit.M2.0$logmortality)
  ## penalty stuff
  Deta.F1 <- diff(diag(m), diff=2)
  veta <- rep(1, nrow(Deta.F1))
  if(up1[1]==0){
    veta[1] <- 0 ## kink at age 0 only 
                 ## when age 0 is explicitly given in pre-pandemic years
  }
  Veta <- diag(veta)
  tDDeta.F1 <- t(Deta.F1)%*%Veta%*%Deta.F1
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
  ## including monotonicity over age (after a certain age)
  kappa.mon <- 10^6
  Dmon <- diff(diag(m), diff=1)
  wmon <- rep(0, nrow(Dmon))
  Wmon <- diag(wmon)
  eps <- 0.01
  epsvec <- rep(eps, m - 1)
  
  ## constraining delta to sum up to 0
  H <- rbind(c(rep(0,2*m+1), rep(1,m), rep(0,m+1)),
             c(rep(0,2*m+1), rep(0,m), 0, rep(1,m)))
  kappa <- c(0,0)

  ## try three lambdas:
  ## common lambdas for both lambda.delta.F and lambda.delta.M
  lambdas.eta.F1 <- 10^seq(3, 7, 1)
  nl.eta <- length(lambdas.eta.F1)
  lambdas.s <- 10^seq(3, 7, 1)
  nl.s <- length(lambdas.s)
  lambdas.delta <- 10^seq(2, 6, 1)
  nl.delta <- length(lambdas.delta)

  BICs <- array(NA, dim=c(nl.eta, nl.s, nl.delta),
                dimnames = list(log10(lambdas.eta.F1),
                                log10(lambdas.s),
                                log10(lambdas.delta)))
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
          d <- ifelse(it<=20, ds[j], 0)
          ## penalty for the monotonicity (only for etaF1)
          Pmon.eta <- kappa.mon * t(Dmon) %*% Wmon %*% Dmon
          Pmon <- matrix(0,ncol(P), ncol(P))
          Pmon[1:m,1:m] <- Pmon.eta
          vmon.eta <- kappa.mon * t(Dmon) %*% (wmon * epsvec)
          vmon <- rep(0, ncol(P))
          vmon[1:m] <- vmon.eta

          gamma   <- exp(eta)
          mu      <- c(C %*% gamma)
          mu[mu==0] <- 1e-8
          X       <- ( C * ((1 / mu) %*% t(gamma)) ) %*% U
          w       <- as.vector(mu)
          r       <- y - mu + C %*% (gamma * eta)
          G       <- t(X) %*% (w * X)
          GpP     <- G + P + Pr + Pmon
          tXr     <- t(X) %*% r + vmon
          ## adding constraints
          LHS <- rbind(cbind(GpP, t(H)),
                       cbind(H, zeros))
          RHS <- matrix(c(tXr, kappa), ncol=1)

          ## check eventual problem with convergence?
          tryInv <- try(solve(LHS, RHS), silent=TRUE)
          ## if so, set "conv" equal to FALSE
          ## and break the loop w/o interrupting
          ## the overall lambdas optimization
          if(!is.matrix(tryInv) |
             any(is.na(tryInv)) |
             any(is.nan(tryInv))){
            it <- 100
            break
          }
          coeff <- tryInv

          betas   <- coeff[1:(m*4+2)]
          eta.old <- eta
          eta <- U%*%betas
          eta <- d*eta.old + (1-d)*eta
          # eta.old <- eta
          # eta     <- U%*%betas

          ## update wmon
          betas.eta.F1 <- betas[1:m]
          diff.betas.eta.F1 <- diff(betas.eta.F1) >= eps
          Wmon.old <- Wmon
          wmon <- rep(1, nrow(Dmon))
          wmon[diff.betas.eta.F1] <- 0
          ## only from age 50
          wmon[1:50] <- 0
          Wmon <- diag(wmon)
          ## convergence check
          dif.eta <- max(abs((eta - eta.old)/eta.old) )
          if(dif.eta < 1e-04 & it > 4) break
          ## cat(it, dif.eta, "\n")
        }
        ## diagnostics
        if(it<100){
          Hat <- solve(GpP, G)
          ED <- sum(diag(Hat))
          y1 <- y
          y1[y == 0] <- 10^(-4)
          DEV <- 2 * sum( y1 * log(y1/mu) )
          AIC <- DEV + 2 * ED
          BIC <- DEV + log(length(y)) * ED
          psi <- DEV/(length(y)-ED)
          QIC <- length(y) + ED + length(y)*log(psi)
          BICs[l.eta, l.s, l.delta] <- BIC
        }
        # cat(log10(lambda.eta.F1),
        #     log10(lambda.s),
        #     log10(lambda.delta.F), it, dif.eta, "\n")
      }
    }
  }
  (whimin <- which(BICs==min(BICs, na.rm=TRUE), arr.ind=TRUE))
  nl.eta;nl.s;nl.delta

  lambda.eta.F1 <- lambdas.eta.F1[whimin[1]]
  lambda.s <- lambdas.s[whimin[2]]
  lambda.delta.F <- lambda.delta.M <- lambdas.delta[whimin[3]]
  # # 
  # lambda.eta.F1 <- 10^3
  # lambda.s <- 10^5
  # lambda.delta.F <- lambda.delta.M <- 10^6

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
    d <- ifelse(it<=20, ds[j], 0)
    ## penalty for the monotonicity (only for etaF1)
    Pmon.eta <- kappa.mon * t(Dmon) %*% Wmon %*% Dmon
    Pmon <- matrix(0,ncol(P), ncol(P))
    Pmon[1:m,1:m] <- Pmon.eta
    vmon.eta <- kappa.mon * t(Dmon) %*% (wmon * epsvec)
    vmon <- rep(0, ncol(P))
    vmon[1:m] <- vmon.eta
    
    gamma   <- exp(eta)
    mu      <- c(C %*% gamma)
    mu[mu==0] <- 1e-8
    X       <- ( C * ((1 / mu) %*% t(gamma)) ) %*% U
    w       <- as.vector(mu)
    r       <- y - mu + C %*% (gamma * eta)
    G       <- t(X) %*% (w * X) 
    GpP     <- G + P + Pr + Pmon
    tXr     <- t(X) %*% r + vmon
    ## adding constraints
    LHS <- rbind(cbind(GpP, t(H)),
                 cbind(H, zeros))
    RHS <- matrix(c(tXr, kappa), ncol=1)
    coeff   <- solve(LHS, RHS)

    betas   <- coeff[1:(m*4+2)]
    eta.old <- eta
    eta <- U%*%betas
    eta <- d*eta.old + (1-d)*eta
    
    # eta.old <- eta
    # eta     <- U%*%betas
    # 
    ## update wmon
    betas.eta.F1 <- betas[1:m]
    diff.betas.eta.F1 <- diff(betas.eta.F1) >= eps
    Wmon.old <- Wmon
    wmon <- rep(1, nrow(Dmon))
    wmon[diff.betas.eta.F1] <- 0
    ## only from age 50
    wmon[1:50] <- 0
    Wmon <- diag(wmon)
    ## convergence check
    dif.eta <- max(abs((eta - eta.old)/eta.old) )
    if(dif.eta < 1e-04 & it > 4) break
    cat(j, code, it, dif.eta, "\n")
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
  
  ## about CI
  H0 <- solve(GpP)
  Vbetas <- H0 %*% G %*% H0
  se.betas <- sqrt(diag(Vbetas))
  se.eta1F <- se.betas[1:m]
  se.s <- se.betas[1:m+m]
  se.c.F <- se.betas[2*m+1]
  se.delta.F <- se.betas[1:m+1+2*m]
  se.c.M <- se.betas[3*m+2]
  se.delta.M <- se.betas[1:m+3*m+2]
  
  whina <- which(is.nan(se.delta.F))
  se.delta.F[whina] <- mean(se.delta.F[-whina])
  whina <- which(is.nan(se.delta.M))
  se.delta.M[whina] <- mean(se.delta.M[-whina])
  
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
  
  ## preparing output 
  out <- list(## observed data, intervals and log-rates
              low1=low1, up1=up1, low2=low2, up2=up2, 
              lmx.F1g=lmx.F1g, lmx.M1g=lmx.M1g,
              lmx.F2g=lmx.F2g, lmx.M2g=lmx.M2g,
              ## pre-pandemic years, log-mortality
              eta.F1.hatL=eta.F1.hatL, eta.F1.hatU=eta.F1.hatU, 
              eta.F1.hat=eta.F1.hat,
              eta.M1.hatL=eta.M1.hatL, eta.M1.hatU=eta.M1.hatU, 
              eta.M1.hat=eta.M1.hat,
              ## pandemic year, log-mortality
              eta.F2.hatL=eta.F2.hatL, eta.F2.hatU=eta.F2.hatU, 
              eta.F2.hat=eta.F2.hat,
              eta.M2.hatL=eta.M2.hatL, eta.M2.hatU=eta.M2.hatU, 
              eta.M2.hat=eta.M2.hat,
              ## model parameters
              ## sex-factor over age
              s.hatL=s.hatL, s.hatU=s.hatU, s.hat=s.hat,
              ## scaling factor
              c.F.hatL=c.F.hatL, c.F.hatU=c.F.hatU, 
              c.F.hat=c.F.hat,
              c.M.hatL=c.M.hatL, c.M.hatU=c.M.hatU, 
              c.M.hat=c.M.hat,
              ## perturbation function over age
              delta.F.hatL=delta.F.hatL, delta.F.hatU=delta.F.hatU,
              delta.F.hat=delta.F.hat,
              delta.M.hatL=delta.M.hatL, delta.M.hatU=delta.M.hatU,
              delta.M.hat=delta.M.hat)
  
  OUTPUT[[j]] <- out
  
  
  par(mfrow=c(3,2))
  ## year(s) 1
  rany <-range(lmx.F1g, lmx.M1g, lmx.F2g, lmx.M2g, finite=TRUE)
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
  for(i in 1:n1){
    segments(x0=low1[i], x1=up1[i],
             y0=lmx.F1g[i], y1=lmx.F1g[i], col=3, lwd=3)
  }
  yy <- c(eta.M1.hatL, rev(eta.M1.hatU))
  polygon(xx, yy, col=adjustcolor(2, 0.5), border=adjustcolor(2, 0.5))
  lines(x, eta.M1.hat, col=2, lwd=4)
  for(i in 1:n1){
    segments(x0=low1[i], x1=up1[i],
             y0=lmx.M1g[i], y1=lmx.M1g[i], col=3, lwd=3)
  }
  ## year 2
  rany <-range(lmx.F1g, lmx.M1g, lmx.F2g, lmx.M2g, finite=TRUE)
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
  for(i in 1:n2){
    segments(x0=low2[i], x1=up2[i],
             y0=lmx.F2g[i], y1=lmx.F2g[i], col=3, lwd=3)
  }
  yy <- c(eta.M2.hatL, rev(eta.M2.hatU))
  polygon(xx, yy, col=adjustcolor(2, 0.5), border=adjustcolor(2, 0.5))
  lines(x, eta.M2.hat, col=2, lwd=4)
  for(i in 1:n2){
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
  legend("topleft", inset=0.1,
         legend=c("Actual", "Fit + 95% CI"),
         col=c(3, adjustcolor(2, 0.5)), pch=c(NA, NA, 15), 
         pt.cex=2, lty=c(1,1),
         lwd=c(3,4), cex=1.5)
  par(mfrow=c(1,1))
}
tim2 <- Sys.time()
tim2-tim1


#save.image("AllDataNoOptimization.RData")



## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## !!!! to be changed
setwd("~/WORK/TAG_patterns/R")
library(MortalitySmooth)
library(magic)
library(colorspace)
library(ggplot2)
library(plotly)
load("AllDataOptimization.RData")


countries.offset <- offset$Country
countries.deaths <- deaths$Country
countries <- intersect(countries.offset,countries.offset)


names(OUTPUT) <- countries


## spaghetti plots with only estimated log-mortality

## 2x2 by sex and pre- and pandemic years

DF <- expand.grid(list(ages=x, 
                       country=countries, 
                       sex=factor(c("Females", "Males")),
                       period=factor(c("pre-pandemic", "pandemic"),
                                     levels=c("pre-pandemic", "pandemic"))
                       ))


DF$eta <- c(unlist(lapply(OUTPUT, "[[", "eta.F1.hat")),
            unlist(lapply(OUTPUT, "[[", "eta.M1.hat")),
            unlist(lapply(OUTPUT, "[[", "eta.F2.hat")),
            unlist(lapply(OUTPUT, "[[", "eta.M2.hat")))




p1 <- DF %>%
  #filter(period=="pre-pandemic" & sex=="Females") %>%
  ggplot(aes(x=ages, y=eta)) + 
  geom_line() +
  xlab("age") +
  ylab("log-mortality") + 
  aes(colour = country)+
  facet_grid(vars(sex), vars(period))+
  scale_colour_viridis_d(option = "plasma")+
  scale_y_continuous(limits = c(-11, 1))+
  theme(legend.position='none',
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        text = element_text(size=20))+
  ggtitle("Fitted log-mortality")

# ggplotly(p1)
pdf("/home/gccamarda/WORK/TAG_patterns/R/LogMortalityAll.pdf", 
    width = 16, height = 10)
p1
dev.off()

## plotting a actual and fitted log-mortality 
## for a specific population
j=63

pdf("/home/gccamarda/WORK/TAG_patterns/R/LogMortalityObsHatEach.pdf", 
    width = 16, height = 10)

for(j in 1:nc){
country <- countries[j]
OUTPUT.i <- OUTPUT[[j]]

## observed
df.obs1 <- expand.grid(list(x=OUTPUT.i$low1, 
                            sex=factor(c("Females", "Males")),
                            period="pre-pandemic"))
df.obs1$xend <- OUTPUT.i$up1+1
df.obs1$y    <- c(OUTPUT.i$lmx.F1g, OUTPUT.i$lmx.M1g)
df.obs1$yend <- df.obs1$y
df.obs2 <- expand.grid(list(x=OUTPUT.i$low2, 
                            sex=factor(c("Females", "Males")),
                            period="pandemic"))
df.obs2$xend <- OUTPUT.i$up2+1
df.obs2$y    <- c(OUTPUT.i$lmx.F2g, OUTPUT.i$lmx.M2g)
df.obs2$yend <- df.obs2$y
df.obs <- rbind(df.obs1,df.obs2)
df.obs$up <- NA
df.obs$low <- NA
df.obs$data <- "Observed"

## fitted
df.hat1 <- expand.grid(list(x=x+0.5, 
                            sex=factor(c("Females", "Males")),
                            period="pre-pandemic"))
df.hat1$xend <- NA
df.hat1$y    <- c(OUTPUT.i$eta.F1.hat, OUTPUT.i$eta.M1.hat)
df.hat1$yend <- NA
df.hat1$up   <- c(OUTPUT.i$eta.F1.hatU, OUTPUT.i$eta.M1.hatU)
df.hat1$low  <- c(OUTPUT.i$eta.F1.hatL, OUTPUT.i$eta.M1.hatL)
df.hat2 <- expand.grid(list(x=x+0.5, 
                            sex=factor(c("Females", "Males")),
                            period="pandemic"))
df.hat2$xend <- NA
df.hat2$y    <- c(OUTPUT.i$eta.F2.hat, OUTPUT.i$eta.M2.hat)
df.hat2$yend <- NA
df.hat2$up   <- c(OUTPUT.i$eta.F2.hatU, OUTPUT.i$eta.M2.hatU)
df.hat2$low  <- c(OUTPUT.i$eta.F2.hatL, OUTPUT.i$eta.M2.hatL)

df.hat <- rbind(df.hat1, df.hat2)
df.hat$data <- "Fitted"

df <- rbind(df.obs, df.hat)

p <- ggplot(df, aes(x=x, y=y, color=data)) +
  geom_segment(data=filter(df, data=="Observed"),
               aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_line(data=filter(df, data=="Fitted"),
            aes(y=y))+
  geom_ribbon(data=filter(df, data=="Fitted"),
    aes(ymin=low, ymax=up), alpha=.3)+
  facet_grid(vars(sex), vars(period))+
  xlab("age") +
  ylab("log-mortality")+
  scale_color_brewer(palette="Dark2")+
  theme(legend.position = c(0.65, 0.9),
        legend.background = element_rect(fill = "white", 
                                         color = "black"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text = element_markdown(size=16),
        legend.title = element_blank(),
        text = element_text(size=20))+
  ggtitle(paste(country))
print(p)
}
dev.off()
## plotting deltas all together, only fitted

DF <- expand.grid(list(ages=x, 
                       country=countries, 
                       sex=factor(c("Females", "Males"))
                       ))
DF$delta <- c(unlist(lapply(OUTPUT, "[[", "delta.F.hat")),
              unlist(lapply(OUTPUT, "[[", "delta.M.hat")))


## issues with Tunisia, Bolivia

p1 <- DF %>%
  filter(country!="Tunisia" & country!="Bolivia") %>%
  ggplot(aes(x=ages, y=delta)) + 
  geom_line() +
  xlab("age") +
  ylab("delta") + 
  aes(colour = country)+
  facet_wrap(~sex)+
  scale_colour_viridis_d(option = "plasma")+
  theme(legend.position='none',
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        text = element_text(size=20))+
  ggtitle("Delta, w/o Tunisia and Bolivia")

# ggplotly(p1)
pdf("/home/gccamarda/WORK/TAG_patterns/R/DeltaAll.pdf", 
    width = 16, height = 10)
p1
dev.off()



## s and deltas for each population


j=34

pdf("/home/gccamarda/WORK/TAG_patterns/R/ModelParametersEach.pdf", 
    width = 16, height = 10)
for(j in 1:nc){
(country <- countries[j])
OUTPUT.i <- OUTPUT[[j]]

df.par <- expand.grid(list(x=x, 
                           parameter=factor(c("s", "deltaF", "deltaM"))
                           ))
df.par$hat <- c(OUTPUT.i$s.hat, 
                OUTPUT.i$delta.F.hat,
                OUTPUT.i$delta.M.hat)
df.par$low <- c(OUTPUT.i$s.hatU,
                 OUTPUT.i$delta.F.hatU,
                 OUTPUT.i$delta.M.hatU)
df.par$up <- c(OUTPUT.i$s.hatL,
                 OUTPUT.i$delta.F.hatL,
                 OUTPUT.i$delta.M.hatL)


p <- df.par %>% 
  #filter(parameter=="s") %>%
  ggplot(aes(x = x, y = hat, color=parameter)) + 
  geom_line()+
  geom_ribbon(aes(ymin=low, ymax=up), alpha=.3)+
  facet_wrap(~parameter, scales = "free")+
  xlab("age") +
  ylab("model parameters") + 
  theme(legend.position='none',
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        text = element_text(size=20))+
  geom_hline(yintercept=0, linetype="dashed", color = "grey20")+
  ggtitle(paste(country))
print(p)
}
dev.off()

## scaling factors

pdf("/home/gccamarda/WORK/TAG_patterns/R/ScalingFactor.pdf", 
    width = 16, height = 10)


DFm <- data.frame(countries=countries)
DFm$hat <- unlist(lapply(OUTPUT, "[[", "c.M.hat"))
DFm$low <- unlist(lapply(OUTPUT, "[[", "c.M.hatL"))
DFm$up <- unlist(lapply(OUTPUT, "[[", "c.M.hatU"))

aa <- DFm %>%
  filter(countries!="Tunisia") %>%
  mutate(countries = fct_reorder(countries, hat, .desc = TRUE))%>%
  ggplot(aes(x=hat, y=countries, xmin=low, xmax=up, 
             label=countries))+
  geom_pointrange()+
  geom_point()+
  annotate("rect", ymin=0, ymax=nc,xmin=0,xmax=-Inf, alpha=0.1,
           fill="#4DAF4A")+
  annotate("rect", ymin=-Inf,ymax=Inf,xmin=0,xmax=Inf, alpha=0.1,
           fill="red")+
  geom_text(aes(fontface=1, x=low), 
            hjust=1.1, vjust=0.5, size=3)+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_line(colour="white", size=0.5),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))+
  labs(x = "", y = "")+
  ggtitle("Age-independent scaling factor, Males (w/o Tunisia)")
print(aa)


DFf <- data.frame(countries=countries)
DFf$hat <- unlist(lapply(OUTPUT, "[[", "c.F.hat"))
DFf$low <- unlist(lapply(OUTPUT, "[[", "c.F.hatL"))
DFf$up <- unlist(lapply(OUTPUT, "[[", "c.F.hatU"))

aa <- DFf %>%
  #filter(countries!="Tunisia") %>%
  mutate(countries = fct_reorder(countries, hat, .desc = TRUE))%>%
  ggplot(aes(x=hat, y=countries, xmin=low, xmax=up, 
             label=countries))+
  geom_pointrange()+
  geom_point()+
  annotate("rect", ymin=0, ymax=nc+1,xmin=0,xmax=-Inf, alpha=0.1,
           fill="#4DAF4A")+
  annotate("rect", ymin=-Inf,ymax=Inf,xmin=0,xmax=Inf, alpha=0.1,
           fill="red")+
  geom_text(aes(fontface=1, x=low), 
            hjust=1.1, vjust=0.5, size=3)+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_line(colour="white", size=0.5),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16))+
  labs(x = "", y = "")+
  ggtitle("Age-independent scaling factor, Females")
print(aa)
dev.off()

















## END




