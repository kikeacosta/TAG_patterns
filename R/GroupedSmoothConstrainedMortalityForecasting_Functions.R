################################################################
################################################################
## R-code with a set of functions useful 
## for modelling and forecasting 
## based on CP-splines
## for data provided by age-groups

## Associated to a future paper 

## function to build up B-splines and associated bases for derivatives
BsplineGrad <- function(x, xl, xr, ndx=NULL, deg, knots=NULL){
  if(is.null(knots)){
    dx <- (xr - xl)/ndx
    knots <- seq(xl - deg * dx, xr + deg * dx, by=dx)
    knots <- round(knots, 8)
  }else{
    knots <- knots
    dx <- diff(knots)[1]
  }
  P <- outer(x, knots, MortSmooth_tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff=deg+1)/(gamma(deg+1)*dx^deg)
  B <- (-1)^(deg + 1) * P %*% t(D)
  ##
  knots1 <- knots[-c(1,length(knots))]
  P <- outer(x, knots1, MortSmooth_tpower, deg-1)
  n <- dim(P)[2]
  D <- diff(diag(n),diff=deg)/(gamma(deg)*dx^(deg-1))
  BB <- ((-1)^(deg) * P %*% t(D))/dx
  D <- diff(diag(ncol(BB) + 1))
  C <- BB %*% D
  ##
  out <- list(dx=dx, knots=knots, B=B, C=C)
}

PSinfant <- function(Y, E, lambdas, WEI, infant=TRUE, verbose=FALSE){
  ## dimensions
  m <- nrow(Y)
  n <- ncol(Y)
  a <- 1:m
  t <- 1:n
  ## w/o age 0
  a0 <- a[-1]
  m0 <- m-1
  ## original offset
  OFF <- log(E)
  ## B-splines basis
  ## with infant-specilized coeff
  if(infant){
    ## over ages w/o age 0
    a0min <- min(a0)
    a0max <- max(a0)
    nda0 <- floor(m0/5)
    dega <- 3
    BCa0 <- BsplineGrad(a0, a0min, a0max, nda0, dega)
    Ba0 <- BCa0$B
    nba0 <- ncol(Ba0)
    ## adding infant-specific basis
    Ba <- cbind(0, Ba0)
    Ba <- rbind(c(1, rep(0,nba0)), Ba)
    nba <- ncol(Ba)
    ## basis for the derivatives
    ## over ages w/o age 0
    Ca0 <- BCa0$C
    ## adding infant-specific basis
    Ca <- cbind(0, Ca0)
    Ca <- rbind(c(-1,Ba[2,1:dega+1],rep(0,nba0-dega)),
                Ca)
  }else{
    amin <- min(a)
    amax <- max(a)
    nda <- floor(m/5)
    dega <- 3
    BCa <- BsplineGrad(a, amin, amax, nda, dega)
    Ba <- BCa$B
    nba <- ncol(Ba)
    ## basis for the derivatives
    Ca <- BCa$C
  }
  ## over years
  tmin <- min(t)
  tmax <- max(t)
  ndt <- floor(n/5)
  degt <- 3
  BCt <- BsplineGrad(t, tmin, tmax, ndt, degt)
  Bt <- BCt$B
  nbt <- ncol(Bt)
  ## basis for the derivatives
  Ct <- BCt$C
  ## weights for exposures=0
  WEI[E==0] <- 0
  
  ## tensor product of B-splines for the GLAM
  Ba1 <- kronecker(matrix(1, ncol=nba, nrow=1), Ba)
  Ba2 <- kronecker(Ba, matrix(1, ncol=nba, nrow=1))
  RTBa <- Ba1 * Ba2
  Bt1 <- kronecker(matrix(1, ncol=nbt, nrow=1), Bt)
  Bt2 <- kronecker(Bt, matrix(1, ncol=nbt, nrow=1))
  RTBt <- Bt1 * Bt2
  
  ## penalty terms
  ## over ages
  Da <- diff(diag(nba), diff=2)
  ## no penalization over age for age 0,
  ## with infant=TRUE
  if(infant){
    Da[1,1] <- 0
  }
  tDDa <- t(Da) %*% Da
  ## over years
  Dt <- diff(diag(nbt), diff=2)
  tDDt <- t(Dt) %*% Dt
  ## kronecker product of difference matrices
  Pa <- kronecker(diag(nbt), tDDa)
  Pt <- kronecker(tDDt, diag(nba))
  ## smoothing parameters
  lambda.a <- lambdas[1]
  lambda.t <- lambdas[2]
  P <- lambda.a * Pa + lambda.t * Pt
  ## data in vector for the starting values
  y <- c(Y)
  e <- c(E)
  wei <- c(WEI)
  off <- log(e)
  off0 <- off
  off0[wei==0] <- 100
  ## "other" offset in matrices
  OFF0 <- matrix(off0, m, n)
  ## starting coeff
  aa <- rep(a, n)
  tt <- rep(t, each = m)
  fit0 <- glm(round(y) ~ aa + tt + offset(off0), 
              family = poisson, weights = wei)
  etaGLM <- matrix(log(fit0$fitted) - c(off0), m, n)
  eta0 <- log((Y + 1)) - OFF0
  eta0[WEI == 0] <- etaGLM[WEI == 0]
  BBa <- solve(t(Ba) %*% Ba + 1e-06*diag(nba), t(Ba))
  BBt <- solve(t(Bt) %*% Bt + 1e-06*diag(nbt), t(Bt))
  alphas <- MortSmooth_BcoefB(BBa, BBt, eta0)
  ## Poisson iteration within a GLAM framework
  for(it in 1:20){
    eta <- MortSmooth_BcoefB(Ba, Bt, alphas)
    mu <- exp(OFF0 + eta)
    W <- mu
    z <- eta + (1/mu) * (Y - mu)
    z[which(WEI == 0)] <- 0
    WW <- WEI * W
    BWB <- MortSmooth_BWB(RTBa, RTBt, nba, nbt, WW)
    BWBpP <- BWB + P
    BWz <- MortSmooth_BcoefB(t(Ba), t(Bt), (WW * z))
    alphas0 <- solve(BWBpP, c(BWz))
    alphas.old <- alphas
    alphas <- matrix(alphas0, nrow = nba)
    dalphas <- max(abs(alphas-alphas.old)/abs(alphas))
    if(verbose) cat(it, dalphas, "\n")
    if(dalphas<=10^-6 & it>=4) break
  }
  ## fitted values
  ALPHAS.hat <- alphas
  ETA.hat <- eta
  Y.hat <- exp(OFF0 + ETA.hat)
  ## fitted derivatives over ages
  ETA.hat1a <- MortSmooth_BcoefB(Ca, Bt, ALPHAS.hat)
  ## fitted derivatives over years
  ETA.hat1t <- MortSmooth_BcoefB(Ba, Ct, ALPHAS.hat)
  ## diagnostics
  BWBpP1 <- solve(BWBpP)
  H <- BWBpP1%*%BWB
  h <- diag(H)
  ed <- sum(h)
  y0 <- y
  y0[y==0] <- 10^-8
  dev <- 2 * sum(wei * (y * log(y0/mu) - (y0-mu)))
  bic <- dev + log(sum(wei))*ed
  # ## deviance residuals
  # res0 <- sign(Y - Y.hat)
  # res1 <- log(Y/Y.hat)
  # res2 <- Y - Y.hat
  # res <- res0 * sqrt(2 * (Y * res1 - res2))
  # res[which(is.nan(res))] <- NA
  # res <- matrix(res,m,n)
  # 
  # ## standard errors for the penalized coefficients
  # se <- matrix(Mort2Dsmooth_se(RTBa, RTBt,
  #                              nbx=nba, nby=nbt,
  #                              BWB.P1=BWBpP1),
  #              m, n,
  #              dimnames=list(a, t))
  
  ## return objects
  out <- list(
    ## original data
    Y=Y, E=E, lambdas=lambdas,
    ## diagnostic
    dev=dev, ed=ed, bic=bic, #res=res,
    ## fitted values
    ALPHAS=ALPHAS.hat, ETA=ETA.hat, MU=Y.hat,
    ETA1a=ETA.hat1a, ETA1t=ETA.hat1t, #SE=se,
    ## basis and associated objects
    Ba=Ba, Bt=Bt, Ca=Ca, Ct=Ct,
    ## convergence objects
    dalphas=dalphas, it=it
  )
  return(out)
}




## function for estimating two-dimensional P-CLM
## with B-spline bases penalized over both dimentions
## and addressing infant mortality (when infant=TRUE)
## for a given set of smoothing parameters lambdas
## Whereas Yg (Deaths) should be given by age-groups and single years
## E (Exposures) should be given by single age and year 
## - length of "a" (single ages) MUST be equal to ncol(E)
## - max(a) MUST be smaller than max(a.low)
## - min(a.low) MUST be EQUAL to min(a)
## (previously ungrouped by pclm(ungroup) with an independent estimate
## for each year)


# a.low=ag.low
# Yg=Yg1
# a=a
# E=E1hat
# lambdas=c(10^-1,10^1)
# WEIg=WEIg1
# infant=TRUE
# verbose=FALSE
# kappa.shape=10^5
# age.mon=40
# upper.der=0.15

# a.low=a.low
# Yg1=Yg1
# WEIg1=WEIg1
# a=ages
# E1=E1
# obs.yrs=years1
# for.hor=for.hor
# ETA1hat=ETA1hat
# deltas=deltasi
# S=S
# lambdas = OPT$par
# verbose = FALSE
# kappa.shape=10^5
# age.mon=40
# upper.der=0.15
# a.low = ag.low
# Yg=Yg1
# a=a
# E=E1
# WEIg=WEIg1
# lambdas=c(10^5, 10^5)
# verbose=TRUE
# kappa.shape=0
# age.mon=40
# upper.der=0.15
# infant=FALSE

# a.low = ag.low
# Yg=Yg1
# a=a
# E=E1
# WEIg=WEIg1
# lambdas=c(10^4,10^4)
# verbose=TRUE
# kappa.shape=0
# infant=infantj
# age.mon=40
# upper.der=0.15

PSinfantGrouped <- function(a.low, Yg, a, E, WEIg, lambdas,
                            infant=TRUE, verbose=FALSE, 
                            kappa.shape=10^5, age.mon=40, 
                            upper.der=0.15){
  ## dimensions
  mg <- nrow(Yg) ## length(a.low)
  m <- nrow(E) # length(a)
  n <- ncol(Yg) # == ncol(E)
  t <- 1:n
  
  ## w/o age 0
  a0 <- a[-1]
  m0 <- m-1
  ## B-splines basis
  ## with infant-specialized coeff
  if(infant){
    ## over ages w/o age 0
    a0min <- min(a0)
    a0max <- max(a0)
    nda0 <- floor(m0/5)
    dega <- 3
    BCa0 <- BsplineGrad(a0, a0min, a0max, nda0, dega)
    Ba0 <- BCa0$B
    nba0 <- ncol(Ba0)
    ## adding infant-specific basis
    Ba <- cbind(0, Ba0)
    Ba <- rbind(c(1, rep(0,nba0)), Ba)
    nba <- ncol(Ba)
    ## basis for the derivatives
    ## over ages w/o age 0
    Ca0 <- BCa0$C
    ## adding infant-specific basis
    Ca <- cbind(0, Ca0)
    Ca <- rbind(c(-1,Ba[2,1:dega+1],rep(0,nba0-dega)),
                Ca)
  }else{
    amin <- min(a)
    amax <- max(a)
    nda <- floor(m/5)
    dega <- 3
    BCa <- BsplineGrad(a, amin, amax, nda, dega)
    Ba <- BCa$B
    nba <- ncol(Ba)
    ## basis for the derivatives
    Ca <- BCa$C
  }
  ## over years
  tmin <- min(t)
  tmax <- max(t)
  ndt <- floor(n/5)
  degt <- 3
  BCt <- BsplineGrad(t, tmin, tmax, ndt, degt)
  Bt <- BCt$B
  nbt <- ncol(Bt)
  ## basis for the derivatives
  Ct <- BCt$C
  
  ## complete composite matrix, only for observed period
  ## with fitted ungrouped exposures
  a.up  <- c(a.low[-1]-1, max(a))
  ## length of each age-group
  len <- a.up-a.low+1 
  ## Boolean matrix mapping age-group and single ages
  G <- matrix(0, mg, m)
  for(i in 1:mg){
    a.low.i <- a.low[i]
    a.up.i <- a.up[i]
    wc <- which(a>=a.low.i & a<=a.up.i)
    G[i,wc] <- 1
  }
  In <- diag(n)
  C <- kronecker(In, G)
  ## including exposures that will work as offset
  C[which(C==1)] <- c(E)
  
  ## complete model matrix (no GLAM here since we are in a CLM framework)
  B <- kronecker(Bt, Ba)
  ## penalty
  Da <- diff(diag(nba), diff=2)
  ## no penalization over age for age 0,
  ## with infant=TRUE
  if(infant){
    Da[1,1] <- 0
  }
  tDDa <- t(Da) %*% Da
  ## over years
  Dt <- diff(diag(nbt), diff=2)
  tDDt <- t(Dt) %*% Dt
  ## kronecker product of difference matrices
  Pa <- kronecker(diag(nbt), tDDa)
  Pt <- kronecker(tDDt, diag(nba))
  ## smoothing parameters
  lambda.a <- lambdas[1]
  lambda.t <- lambdas[2]
  P <- lambda.a * Pa + lambda.t * Pt
  
  ## monotonicity after age a0
  G.a.low <- matrix(0, m, n)
  G.a.up <- matrix(upper.der, m, n)
  Sshape <- matrix(1, m, n)
  Sshape[which(a<age.mon),] <- 0
  ## tensor product for the derivatives
  Bt1 <- kronecker(matrix(1, ncol=nbt, nrow=1), Bt)
  Bt2 <- kronecker(Bt, matrix(1, ncol=nbt, nrow=1))
  RTBt <- Bt1 * Bt2
  Ca1 <- kronecker(matrix(1, ncol=nba, nrow=1), Ca)
  Ca2 <- kronecker(Ca, matrix(1, ncol=nba, nrow=1))
  RTCa <- Ca1 * Ca2
  ## data in vector since we are not in a GLAM 
  yg <- c(Yg)
  weig <- c(WEIg)
  
  ## preparing starting values
  ## first approximation: uniform
  y0 <- rep(yg/rep(len,n), rep(len,n))
  e0 <- c(E)
  Y0 <- matrix(y0, m, n)
  E0 <- matrix(e0, m, n)
  # fit0 <- Mort2Dsmooth(x=a, y=t, Z=Y0, offset=log(E0),
  #                      method=3, lambdas=c(100,10000),
  #                      ndx = c(nda, ndt))
  # ALPHAS <- fit0$coef
  
  fit0 <- PSinfant(Y=Y0, E=E, WEI=matrix(1,m,n),
                   lambdas = lambdas*10, infant=infant)
  ALPHAS <- fit0$ALPHAS
  # ETA0 <- log(Y0/E0)
  # BBa <- solve(t(Ba) %*% Ba + 1e-06*diag(nba), t(Ba))
  # BBt <- solve(t(Bt) %*% Bt + 1e-06*diag(nbt), t(Bt))
  # ALPHAS <- MortSmooth_BcoefB(BBa, BBt, ETA0)
  if(verbose) cat(lambdas, "\n")
  for(it in 1:20){
    ## penalty for the monotonicity
    ## over ages
    CaBt.alphas <- MortSmooth_BcoefB(Ca, Bt, ALPHAS)
    ## low
    v.a.low <- CaBt.alphas < G.a.low
    v.a.low <- v.a.low * Sshape
    P.a.low <- MortSmooth_BWB(RTCa, RTBt,
                              nba, nbt, v.a.low)
    P.a.low <- kappa.shape * P.a.low
    ## up
    v.a.up <- CaBt.alphas > G.a.up
    v.a.up <- v.a.up * Sshape
    P.a.up <- MortSmooth_BWB(RTCa, RTBt,
                             nba, nbt, v.a.up)
    P.a.up <- kappa.shape * P.a.up
    p.a.up <- MortSmooth_BcoefB(t(Ca), t(Bt),
                                v.a.up*G.a.up)
    p.a.up <- kappa.shape * c(p.a.up)
    ##
    P.a <- P.a.low + P.a.up
    
    eta <- c(MortSmooth_BcoefB(Ba, Bt, ALPHAS))
    gamma <- exp(eta)
    mu <- C %*% gamma
    mu[mu==0] <- 10^-8
    X <- (C * ((1 / mu) %*% t(gamma)) ) %*% B
    w <- as.vector(weig*mu)
    r <- weig*(yg - mu + C %*% (gamma * eta))
    tXWX <- t(X) %*% (w * X)
    tXWXpP <- tXWX + P + P.a
    tXr <- t(X) %*% r + p.a.up
    alphas.old <- c(ALPHAS)
    alphas <- solve(tXWXpP, tXr)
    ALPHAS <- matrix(alphas, nba, nbt)
    dalphas <- max(abs(alphas-alphas.old)/abs(alphas))
    if(verbose) cat(it, dalphas, "\n")
    if(dalphas < 1e-04 & it > 4) break
  }
  ## fitted values
  ALPHAS.hat <- matrix(alphas, nba, nbt)
  ETA.hat <- matrix(eta, m, n)
  Y.hat <- E * matrix(gamma, m, n)
  ## fitted derivatives over ages
  ETA.hat1a <- MortSmooth_BcoefB(Ca, Bt, ALPHAS.hat)
  ## fitted derivatives over years
  ETA.hat1t <- MortSmooth_BcoefB(Ba, Ct, ALPHAS.hat)
  ## diagnostics
  H <- solve(tXWXpP, tXWX)
  h <- diag(H)
  ed <- sum(h)
  y0 <- yg
  y0[yg==0] <- 10^-8
  dev <- 2 * sum(weig * (yg * log(y0/mu) - (y0-mu)))
  bic <- dev + log(sum(weig))*ed
  ## return objects
  out <- list(
    ## original data
    Yg=Yg, E=E, lambdas=lambdas,
    ## diagnostic
    dev=dev, ed=ed, bic=bic,
    ## fitted values
    ALPHAS=ALPHAS.hat, ETA=ETA.hat, MU=Y.hat,
    ETA1a=ETA.hat1a, ETA1t=ETA.hat1t, 
    ## basis and associated objects
    Ba=Ba, Bt=Bt, Ca=Ca, Ct=Ct, G=G,
    ## convergence objects
    dalphas=dalphas, it=it
  )
  return(out)
}



## function to extract deltas from a PSinfant object
## for a given confidence level, default 95% and 50% 
## over age and year as reccommened/used in the manuscript
deltasFUN <- function(object, levels=c(95,50)){
  ETA1a <- object$ETA1a
  ETA1t <- object$ETA1t
  ## compute levels and deltas over ages
  p.a.up <- (100 - (100-levels[1])/2)/100
  p.a.low <- ((100-levels[1])/2)/100
  delta.a.up <- apply(ETA1a, 1, quantile,
                      probs=p.a.up)
  delta.a.low <- apply(ETA1a, 1, quantile,
                       probs=p.a.low)
  ## compute levels and deltas over years
  p.t.up <- (100 - (100-levels[2])/2)/100
  p.t.low <- ((100-levels[2])/2)/100
  delta.t.up <- apply(ETA1t, 1, quantile,
                      probs=p.t.up)
  delta.t.low <- apply(ETA1t, 1, quantile,
                       probs=p.t.low)
  ## return objects
  out <- list(delta.a.up=delta.a.up,
              delta.a.low=delta.a.low,
              delta.t.up=delta.t.up,
              delta.t.low=delta.t.low)
  return(out)
}


## function for estimating CP-splines with grouped data
## for a given set of smoothing parameters lambdas
## commonly taken from the estimated object, PSinfant()

## Like in PSinfantGrouped:
## Yg (Deaths) should be given by age-groups and single years
## E (Exposures) should be given by single age and year 
## - length of "a" (single ages) MUST be equal to ncol(E)
## - max(a) MUST be smaller than max(a.low)
## - min(a.low) MUST be EQUAL to min(a)
## (previously ungrouped by pclm(ungroup) with an independent estimate
## for each year)
## Additionally:
## - obs.yrs  (vector of observed years) 
##   must have the same # of columns in Yg and E : length(obs.yrs)==ncol(E)==ncol(Yg)
## - for.hor (last year to forecast) must be larger than max(obs.yrs)


# rm(list=ls(all=TRUE))
# options(device="X11")
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# library(MortalitySmooth)
# library(HMDdata)
# library(ggplot2)
# library(hrbrthemes)
# library(tidyverse)
# 
# load("Trial.RData")
# 
# a.low = ag.low
# Yg1=Yg1
# WEIg1=WEIg1
# 
# a=a
# E1=E1hat
# obs.yrs=years1
# for.hor=2050
# 
# ETA1hat <- ETA1hat
# deltas=deltasi
# S=S
# lambdas=OPT$par
# kappas=c(10^4, 10^4)
# infant=TRUE
# verbose=TRUE
# st.alphas=NULL

# a.low = a.low
# Yg=Yg1
# a=ages
# E=E1
# WEIg=WEIg1
# lambdas=OPT$par
# verbose=TRUE
# kappas=c(10^4, 10^4)
# infant=TRUE
# st.alphas=NULL

# a.low=ag.low
# Yg1=Yg1
# WEIg1=WEIg1
# a=a
# E1=E1hat
# obs.yrs=years1
# for.hor=2050
# ETA1hat=ETA1hat
# deltas=deltasi
# S=S
# lambdas = OPT$par
# verbose = TRUE
# kappas=c(10^4, 10^4)
# infant=TRUE
# a.low=ag.low
# Yg1=Yg1
# WEIg1=WEIg1
# a=a
# E1=E1hat
# obs.yrs=years1
# for.hor=max(years)
# ETA1hat=ETA1hat
# deltas=deltasi
# S=S
# lambdas = OPT$par
# kappas=c(10^4, 10^4)
# kappa.shape=10^5
# age.mon=40
# infant=TRUE
# verbose = TRUE
# st.alphas=NULL
# upper.der=0.15

#a.low=ag.low
#Yg1=Yg1
# WEIg1=WEIg1
# a=a
# E1=E1
# obs.yrs=t1
# for.hor=max(t)
# ETA1hat=ETA1hat
# deltas=deltasi
# S=S
# lambdas = OPT$par
# verbose = TRUE
# infant = infantj
# kappa.shape=0
# kappas=c(10^4, 10^4)
# age.mon=40
# upper.der=0.15
CPSfunctionGrouped <- function(a.low, Yg1, WEIg1,
                               a, E1,
                               obs.yrs, for.hor,
                               ETA1hat, deltas,
                               S,lambdas,
                               kappas=c(10^4, 10^4),
                               kappa.shape=10^5, age.mon=40, upper.der=0.15,
                               infant=TRUE,
                               verbose=FALSE,
                               st.alphas=NULL){
  ## dimensions
  mg <- nrow(Yg1) ## length(a.low)
  m <- nrow(E1) # length(a)
  n1 <- ncol(Yg1) # ncol(E1), length(obs.yrs)
  n <- length(min(obs.yrs):for.hor)
  nF <- n - n1  
  t1 <- 1:n1
  t <- 1:n
  tF <- (n1+1):n
  
  ## arbitrary values for forecasting
  Yg <- matrix(10, mg, n)
  Yg[1:mg,1:n1] <- Yg1
  ## 0/1 weights for the arbitrary values
  WEIg <- cbind(WEIg1, matrix(0, mg, n-n1))
  WEI1 <- matrix(1, m, n1)
  WEI1[E1==0] <- 0
  WEI <- matrix(0, m, n)
  WEI[,1:n1] <- WEI1
  
  
  
  ## w/o age 0
  a0 <- a[-1]
  m0 <- m-1
  ## B-splines basis
  ## with infant-specialized coeff
  if(infant){
    ## over ages w/o age 0
    a0min <- min(a0)
    a0max <- max(a0)
    nda0 <- floor(m0/5)
    dega <- 3
    BCa0 <- BsplineGrad(a0, a0min, a0max, nda0, dega)
    Ba0 <- BCa0$B
    nba0 <- ncol(Ba0)
    ## adding infant-specific basis
    Ba <- cbind(0, Ba0)
    Ba <- rbind(c(1, rep(0,nba0)), Ba)
    nba <- ncol(Ba)
    ## basis for the derivatives
    ## over ages w/o age 0
    Ca0 <- BCa0$C
    ## adding infant-specific basis
    Ca <- cbind(0, Ca0)
    Ca <- rbind(c(-1,Ba[2,1:dega+1],rep(0,nba0-dega)),
                Ca)
  }else{
    amin <- min(a)
    amax <- max(a)
    nda <- floor(m/5)
    dega <- 3
    BCa <- BsplineGrad(a, amin, amax, nda, dega)
    Ba <- BCa$B
    nba <- ncol(Ba)
    ## basis for the derivatives
    Ca <- BCa$C
  }
  
  
  ## over years
  tmin <- min(t)
  tmax <- max(t)
  ndt <- floor(n/5)
  BCt <- BsplineGrad(t, tmin, tmax, ndt, deg=3)
  Bt <- BCt$B
  nbt <- ncol(Bt)
  ## basis for the derivatives
  Ct <- BCt$C
  
  ## complete composite matrix
  ## with fitted ungrouped exposures
  a.up  <- c(a.low[-1]-1, max(a))
  ## length of each age-group
  len <- a.up-a.low+1 
  ## Boolean matrix mapping age-group and single ages
  G <- matrix(0, mg, m)
  for(i in 1:mg){
    a.low.i <- a.low[i]
    a.up.i <- a.up[i]
    wc <- which(a>=a.low.i & a<=a.up.i)
    G[i,wc] <- 1
  }
  In1 <- diag(n1)
  InF <- diag(nF)
  C1 <- kronecker(In1, G)
  ## including exposures that will work as offset
  C1[which(C1==1)] <- c(E1)
  CF <- kronecker(InF, G)
  ## including exposures that will work as offset
  CF[which(CF==1)] <- 999
  # complete C matrix
  C <- bdiag(C1, CF)
  
  ## complete model matrix (no GLAM here)
  B <- kronecker(Bt, Ba)
  
  ## penalty
  Da <- diff(diag(nba), diff=2)
  ## no penalization over age for age 0,
  ## with infant=TRUE
  if(infant){
    Da[1,1] <- 0
  }
  tDDa <- t(Da) %*% Da
  
  ## over years
  Dt <- diff(diag(nbt), diff=2)
  tDDt <- t(Dt) %*% Dt
  ## kronecker product of difference matrices
  Pa  <- kronecker(diag(nbt), tDDa)
  Pt <- kronecker(tDDt, diag(nba))
  ## smoothing parameters
  lambda.a <- lambdas[1]
  lambda.t <- lambdas[2]
  P <- lambda.a * Pa + lambda.t * Pt
  
  ## data in vector, no GLAM here
  yg <- c(Yg)
  weig <- c(WEIg)
  wei <- c(WEI)
  
  ## penalty terms for constraints
  ## extract deltas
  delta.a.up <- deltas$delta.a.up
  delta.a.low <- deltas$delta.a.low
  delta.t.up <- deltas$delta.t.up
  delta.t.low <- deltas$delta.t.low
  ## construct G matrices
  ones12 <- matrix(1, n, 1)
  g.a.up <- kronecker(ones12, delta.a.up)
  G.a.up <- matrix(g.a.up, m, n)
  g.a.low <- kronecker(ones12, delta.a.low)
  G.a.low <- matrix(g.a.low, m, n)
  g.t.up <- kronecker(ones12, delta.t.up)
  G.t.up <- matrix(g.t.up, m, n)
  g.t.low <- kronecker(ones12, delta.t.low)
  G.t.low <- matrix(g.t.low, m, n)
  ## monotonicity after age.mon
  G.a.mon <- matrix(0, m, n)
  Sshape <- matrix(1, m, n)
  Sshape[which(a<age.mon),] <- 0
  ## upper limit to the derivatives over ages, as in the fitted section
  G.a.upper.der <- matrix(upper.der, m, n)
  ## tensor product of B-splines for the GLAM
  Ba1 <- kronecker(matrix(1, ncol=nba, nrow=1), Ba)
  Ba2 <- kronecker(Ba, matrix(1, ncol=nba, nrow=1))
  RTBa <- Ba1 * Ba2
  Bt1 <- kronecker(matrix(1, ncol=nbt, nrow=1), Bt)
  Bt2 <- kronecker(Bt, matrix(1, ncol=nbt, nrow=1))
  RTBt <- Bt1 * Bt2
  ## tensor product for the derivatives
  Ca1 <- kronecker(matrix(1, ncol=nba, nrow=1), Ca)
  Ca2 <- kronecker(Ca, matrix(1, ncol=nba, nrow=1))
  RTCa <- Ca1 * Ca2
  Ct1 <- kronecker(matrix(1, ncol=nbt, nrow=1), Ct)
  Ct2 <- kronecker(Ct, matrix(1, ncol=nbt, nrow=1))
  RTCt <- Ct1 * Ct2
  ##
  
  ## starting values
  if(is.null(st.alphas)){
    Y0 <- cbind(E1*exp(ETA1hat), matrix(10,m,nF))
    E0 <- cbind(E1, matrix(999,m,nF))
    ETA0 <- log((Y0+1)/(E0+1))
    BBa <- solve(t(Ba) %*% Ba + 1e-06*diag(nba), t(Ba))
    BBt <- solve(t(Bt) %*% Bt + 1e-06*diag(nbt), t(Bt))
    ALPHAS <- MortSmooth_BcoefB(BBa, BBt, ETA0)
  }else{
    alphas=st.alphas
  }
  for(it in 1:200){
    ## over ages
    CaBt.alphas <- MortSmooth_BcoefB(Ca, Bt, ALPHAS)
    ## up
    v.a.up <- CaBt.alphas > G.a.up
    v.a.up <- v.a.up * S
    P.a.up <- MortSmooth_BWB(RTCa, RTBt,
                             nba, nbt, v.a.up)
    P.a.up <- kappas[1] * P.a.up
    p.a.up <- MortSmooth_BcoefB(t(Ca), t(Bt),
                                v.a.up*G.a.up)
    p.a.up <- kappas[1] * p.a.up
    ## low
    v.a.low <- CaBt.alphas < G.a.low
    v.a.low <- v.a.low * S
    P.a.low <- MortSmooth_BWB(RTCa, RTBt,
                              nba, nbt, v.a.low)
    P.a.low <- kappas[1] * P.a.low
    p.a.low <- MortSmooth_BcoefB(t(Ca), t(Bt),
                                 v.a.low*G.a.low)
    p.a.low <- kappas[1] * p.a.low
    
    ## add monotonicity
    v.a.mon <- CaBt.alphas < G.a.mon
    v.a.mon <- v.a.mon * Sshape
    P.a.mon <- MortSmooth_BWB(RTCa, RTBt,
                              nba, nbt, v.a.mon)
    P.a.mon <- kappa.shape * P.a.mon
    ## adding upper limit to the derivatives
    v.a.upper.der <- CaBt.alphas > G.a.upper.der
    v.a.upper.der <- v.a.upper.der * Sshape
    P.a.upper.der <- MortSmooth_BWB(RTCa, RTBt,
                                    nba, nbt, v.a.upper.der)
    P.a.upper.der <- kappa.shape * P.a.upper.der
    p.a.upper.der <- MortSmooth_BcoefB(t(Ca), t(Bt),
                                       v.a.upper.der*G.a.upper.der)
    p.a.upper.der <- kappa.shape * c(p.a.upper.der)
    
    P.a <- P.a.up + P.a.low + P.a.mon + P.a.upper.der
    p.a <- p.a.up + p.a.low + p.a.upper.der
    ## over years
    BaCt.alphas <- MortSmooth_BcoefB(Ba, Ct, ALPHAS)
    ## up
    v.t.up <- BaCt.alphas > G.t.up
    v.t.up <- v.t.up * S  
    P.t.up <- MortSmooth_BWB(RTBa, RTCt, nba, nbt,
                             v.t.up)
    P.t.up <- kappas[2] * P.t.up
    p.t.up <- MortSmooth_BcoefB(t(Ba), t(Ct),
                                v.t.up*G.t.up)
    p.t.up <- kappas[2] * p.t.up
    ## low
    v.t.low <- BaCt.alphas < G.t.low
    v.t.low <- v.t.low * S  
    P.t.low <- MortSmooth_BWB(RTBa, RTCt, nba, nbt,
                              v.t.low)
    P.t.low <- kappas[2] * P.t.low
    p.t.low <- MortSmooth_BcoefB(t(Ba), t(Ct),
                                 v.t.low*G.t.low)
    p.t.low <- kappas[2] * p.t.low
    
    P.t <- P.t.up + P.t.low
    p.t <- p.t.up + p.t.low
    
    ######
    eta <- c(MortSmooth_BcoefB(Ba, Bt, ALPHAS))
    gamma <- exp(eta)
    mu <- C %*% gamma
    mu[mu==0] <- 10^-8
    X <- (C * ((1 / mu) %*% t(gamma)) ) %*% B
    w <- as.vector(weig*mu)
    r <- weig*(yg - mu + C %*% (gamma * eta))
    tXWX <- t(X) %*% (w * X) 
    tXWXpP <- tXWX + P + P.a + P.t
    tXr <- t(X) %*% r
    tXrpP <- tXr + c(p.a) + c(p.t)
    alphas.old <- c(ALPHAS)
    alphas <- solve(tXWXpP, tXrpP)
    ALPHAS <- matrix(alphas, nba, nbt)
    dalphas <- max(abs(alphas-alphas.old)/abs(alphas))
    if(verbose) cat(it, dalphas, "\n")
    if(dalphas < 1e-04 & it > 4) break
  }
  ## fitted values
  ALPHAS.hat <- ALPHAS
  ETA.hat <- matrix(eta, m, n)
  Y.hat <- cbind(E1, matrix(999, m, nF)) *  exp(ETA.hat)
  ## fitted derivatives over ages
  ETA.hat1a <- MortSmooth_BcoefB(Ca, Bt, ALPHAS.hat)
  ## fitted derivatives over years
  ETA.hat1t <- MortSmooth_BcoefB(Ba, Ct, ALPHAS.hat)
  
  ## variance-covariance matrix for the coefficients alpha
  tXWX <- t(X) %*% (w * X) 
  tXWXpP1 <- solve(tXWX + P)#solve(tXWX + P + P.a + P.t)
  Valpha <- tXWXpP1 %*% tXWX %*% tXWXpP1
  
  H <- solve(tXWX + P, tXWX)
  h <- diag(H)
  ed <- sum(h)
  
  y0 <- yg
  y0[yg==0] <- 10^-8
  dev <- 2 * sum(weig * (yg * log(y0/mu) - (y0-mu)))
  
  ## return objects
  out <- list(
    ## original data
    a.low = a.low, Yg1=Yg1, WEIg1=WEIg1,
    a=a,E1=E1,obs.yrs=obs.yrs,for.hor=for.hor,
    lambdas=lambdas,
    kappas=kappas, deltas=deltas, S=S,
    ## fitted values
    ALPHAS=ALPHAS.hat, ETA=ETA.hat, MU=Y.hat,
    ETA1a=ETA.hat1a, ETA1t=ETA.hat1t, 
    ## basis and associated objects
    Ba=Ba, Bt=Bt, Ca=Ca, Ct=Ct, G=G,
    ## convergence
    dalphas=dalphas, it=it,
    ## uncertainty
    Valpha=Valpha,
    ## diagnostics
    dev=dev, ed=ed
  )
  return(out)
}

## function for estimating CP-splines
## for a given set of smoothing parameters lambdas
## commonly taken from the estimated object, PSinfant()
CPSfunction <- function(Y, E, lambdas, WEI,
                        kappas=c(10^4, 10^4),
                        deltas, S,
                        infant=TRUE, verbose=FALSE,
                        st.alphas=NULL){
  ## dimensions
  m <- nrow(Y)
  n <- ncol(Y)
  a <- 1:m
  t <- 1:n
  ## w/o age 0
  a0 <- a[-1]
  m0 <- m-1
  ## original offset
  OFF <- log(E)
  ## B-splines basis
  ## with infant-specilized coeff
  if(infant){
    ## over ages w/o age 0
    a0min <- min(a0)
    a0max <- max(a0)
    nda0 <- floor(m0/5)
    dega <- 3
    BCa0 <- BsplineGrad(a0, a0min, a0max, nda0, dega)
    Ba0 <- BCa0$B
    nba0 <- ncol(Ba0)
    ## adding infant-specific basis
    Ba <- cbind(0, Ba0)
    Ba <- rbind(c(1, rep(0,nba0)), Ba)
    nba <- ncol(Ba)
    ## basis for the derivatives
    ## over ages w/o age 0
    Ca0 <- BCa0$C
    ## adding infant-specific basis
    Ca <- cbind(0, Ca0)
    Ca <- rbind(c(-1,Ba[2,1:dega+1],rep(0,nba0-dega)),
                Ca)
  }else{
    amin <- min(a)
    amax <- max(a)
    nda <- floor(m/5)
    dega <- 3
    BCa <- BsplineGrad(a, amin, amax, nda, dega)
    Ba <- BCa$B
    nba <- ncol(Ba)
    ## basis for the derivatives
    Ca <- BCa$C
  }
  ## over years
  tmin <- min(t)
  tmax <- max(t)
  ndt <- floor(n/5)
  degt <- 3
  BCt <- BsplineGrad(t, tmin, tmax, ndt, degt)
  Bt <- BCt$B
  nbt <- ncol(Bt)
  ## basis for the derivatives
  Ct <- BCt$C
  ## weights for exposures=0
  WEI[E==0] <- 0
  
  ## tensor product of B-splines for the GLAM
  Ba1 <- kronecker(matrix(1, ncol=nba, nrow=1), Ba)
  Ba2 <- kronecker(Ba, matrix(1, ncol=nba, nrow=1))
  RTBa <- Ba1 * Ba2
  Bt1 <- kronecker(matrix(1, ncol=nbt, nrow=1), Bt)
  Bt2 <- kronecker(Bt, matrix(1, ncol=nbt, nrow=1))
  RTBt <- Bt1 * Bt2
  ## tensor product for the derivatives
  Ca1 <- kronecker(matrix(1, ncol=nba, nrow=1), Ca)
  Ca2 <- kronecker(Ca, matrix(1, ncol=nba, nrow=1))
  RTCa <- Ca1 * Ca2
  Ct1 <- kronecker(matrix(1, ncol=nbt, nrow=1), Ct)
  Ct2 <- kronecker(Ct, matrix(1, ncol=nbt, nrow=1))
  RTCt <- Ct1 * Ct2
  
  
  ## penalty terms for smoothing
  ## over ages
  Da <- diff(diag(nba), diff=2)
  ## no penalization over age for age 0,
  ## with infant=TRUE
  if(infant){
    Da[1,1] <- 0
  }
  tDDa <- t(Da) %*% Da
  ## over years
  Dt <- diff(diag(nbt), diff=2)
  tDDt <- t(Dt) %*% Dt
  ## kronecker product of difference matrices
  Pa <- kronecker(diag(nbt), tDDa)
  Pt <- kronecker(tDDt, diag(nba))
  ## smoothing parameters
  lambda.a <- lambdas[1]
  lambda.t <- lambdas[2]
  P <- lambda.a * Pa + lambda.t * Pt
  
  ## penalty terms for constraints
  ## extract deltas
  delta.a.up <- deltas$delta.a.up
  delta.a.low <- deltas$delta.a.low
  delta.t.up <- deltas$delta.t.up
  delta.t.low <- deltas$delta.t.low
  ## construct G matrices
  ones12 <- matrix(1, n, 1)
  g.a.up <- kronecker(ones12, delta.a.up)
  G.a.up <- matrix(g.a.up, m, n)
  g.a.low <- kronecker(ones12, delta.a.low)
  G.a.low <- matrix(g.a.low, m, n)
  g.t.up <- kronecker(ones12, delta.t.up)
  G.t.up <- matrix(g.t.up, m, n)
  g.t.low <- kronecker(ones12, delta.t.low)
  G.t.low <- matrix(g.t.low, m, n)
  
  ## data in vector for the starting values
  y <- c(Y)
  e <- c(E)
  wei <- c(WEI)
  off <- log(e)
  off0 <- off
  off0[wei==0] <- 100
  ## "other" offset in matrices
  OFF0 <- matrix(off0, m, n)
  if(is.null(st.alphas)){
    ## starting coeff
    aa <- rep(a, n)
    tt <- rep(t, each = m)
    fit0 <- glm(round(y) ~ aa + tt + offset(off0), 
                family = poisson, weights = wei)
    etaGLM <- matrix(log(fit0$fitted) - c(off0), m, n)
    eta0 <- log((Y + 1)) - OFF0
    eta0[WEI == 0] <- etaGLM[WEI == 0]
    BBa <- solve(t(Ba) %*% Ba + 1e-06*diag(nba), t(Ba))
    BBt <- solve(t(Bt) %*% Bt + 1e-06*diag(nbt), t(Bt))
    alphas <- MortSmooth_BcoefB(BBa, BBt, eta0)
  }else{
    alphas=st.alphas
  }
  ## Poisson iteration with asymmetric penalty in a GLAM framework
  for(it in 1:100){
    ## over ages
    CaBt.alphas <- MortSmooth_BcoefB(Ca, Bt, alphas)
    ## up
    v.a.up <- CaBt.alphas > G.a.up
    v.a.up <- v.a.up * S
    P.a.up <- MortSmooth_BWB(RTCa, RTBt,
                             nba, nbt, v.a.up)
    P.a.up <- kappas[1] * P.a.up
    p.a.up <- MortSmooth_BcoefB(t(Ca), t(Bt),
                                v.a.up*G.a.up)
    p.a.up <- kappas[1] * p.a.up
    ## low
    v.a.low <- CaBt.alphas < G.a.low
    v.a.low <- v.a.low * S
    P.a.low <- MortSmooth_BWB(RTCa, RTBt,
                              nba, nbt, v.a.low)
    P.a.low <- kappas[1] * P.a.low
    p.a.low <- MortSmooth_BcoefB(t(Ca), t(Bt),
                                 v.a.low*G.a.low)
    p.a.low <- kappas[1] * p.a.low
    
    P.a <- P.a.up + P.a.low
    p.a <- p.a.up + p.a.low
    ## over years
    BaCt.alphas <- MortSmooth_BcoefB(Ba, Ct, alphas)
    ## up
    v.t.up <- BaCt.alphas > G.t.up
    v.t.up <- v.t.up * S  
    P.t.up <- MortSmooth_BWB(RTBa, RTCt, nba, nbt,
                             v.t.up)
    P.t.up <- kappas[2] * P.t.up
    p.t.up <- MortSmooth_BcoefB(t(Ba), t(Ct),
                                v.t.up*G.t.up)
    p.t.up <- kappas[2] * p.t.up
    ## low
    v.t.low <- BaCt.alphas < G.t.low
    v.t.low <- v.t.low * S  
    P.t.low <- MortSmooth_BWB(RTBa, RTCt, nba, nbt,
                              v.t.low)
    P.t.low <- kappas[2] * P.t.low
    p.t.low <- MortSmooth_BcoefB(t(Ba), t(Ct),
                                 v.t.low*G.t.low)
    p.t.low <- kappas[2] * p.t.low
    
    P.t <- P.t.up + P.t.low
    p.t <- p.t.up + p.t.low
    
    eta <- MortSmooth_BcoefB(Ba, Bt, alphas)
    mu <- exp(OFF0 + eta)
    W <- mu
    z <- eta + (1/mu) * (Y - mu)
    z[which(WEI == 0)] <- 0
    WW <- WEI * W
    BWB <- MortSmooth_BWB(RTBa, RTBt, nba, nbt, WW)
    BWBpP <- BWB + P + P.a + P.t
    BWz <- MortSmooth_BcoefB(t(Ba), t(Bt), (WW * z))
    BWzpP <- BWz + p.a + p.t
    alphas0 <- solve(BWBpP, c(BWzpP))
    alphas.old <- alphas
    alphas <- matrix(alphas0, nrow = nba)
    dalphas <- max(abs(alphas-alphas.old)/abs(alphas))
    if(verbose) cat(it, dalphas, "\n")
    if(dalphas<=10^-4 & it>=4) break   
  }
  ## fitted values
  ALPHAS.hat <- alphas
  ETA.hat <- eta
  Y.hat <- exp(OFF0 + ETA.hat)
  ## fitted derivatives over ages
  ETA.hat1a <- MortSmooth_BcoefB(Ca, Bt, ALPHAS.hat)
  ## fitted derivatives over years
  ETA.hat1t <- MortSmooth_BcoefB(Ba, Ct, ALPHAS.hat)
  ## diagnostics
  H <- solve(BWBpP, BWB)
  h <- diag(H)
  ed <- sum(h)
  y0 <- y
  y0[y==0] <- 10^-8
  dev <- 2 * sum(wei * (y * log(y0/mu) - (y0-mu)))
  bic <- dev + log(sum(wei))*ed
  
  ## variance-covariance matrix for the coefficients alpha
  BWB <- MortSmooth_BWB(RTBa, RTBt, nba, nbt, WEI *Y.hat)
  BWBpP1 <- solve(BWB + P)#solve(BWB + P + P.a + P.t)
  Valpha <- BWBpP1 %*% BWB %*% BWBpP1
  #SE.ETA.hat <- Mort2Dsmooth_se(RTBa, RTBt, nba, nbt, BWB.P1=BWBpP1)
  
  ## return objects
  out <- list(
    ## original data
    Y=Y, E=E, lambdas=lambdas,
    kappas=kappas, deltas=deltas, S=S,
    ## diagnostic
    dev=dev, ed=ed, bic=bic,
    ## fitted values
    ALPHAS=ALPHAS.hat, ETA=ETA.hat, MU=Y.hat,
    ETA1a=ETA.hat1a, ETA1t=ETA.hat1t, 
    ## basis and associated objects
    Ba=Ba, Bt=Bt, Ca=Ca, Ct=Ct,
    ## convergence
    dalphas=dalphas, it=it,
    ## uncertainty
    Valpha=Valpha#, SE.ETA.hat=SE.ETA.hat
  )
  return(out)
}





## function to estimate \delta and c
PertubFUN <- function(ag.low, obs.y, obs.e, a, eta.for,
                      lambda){
  m <- length(a)
  ag.up1  <- ag.low-1
  if(ag.low[2]==1){
    ag.up <- c(0, ag.up1[ag.up1>0], max(a))
  }else{
    ag.up <- c(ag.up1[ag.up1>0], max(a))
  }
  ag.mid <- (ag.low+ag.up)/2 
  ag.lab <- paste(ag.low, ag.up,sep="-")
  
  # lg <- ag.up-ag.low+1
  mg <- length(ag.low)
  # G20 <- matrix(0, mg, m)
  # rownames(G20) <- ag.mid
  # colnames(G20) <- a
  # for(ii in 1:mg){
  #   ag.low.i <- ag.low[ii]
  #   ag.up.i <- ag.up[ii]
  #   wc <- which(a>=ag.low.i & a<=ag.up.i)
  #   G20[ii,wc] <- 1
  # }
  # ag.low <- seq(0,95,by=5)
  # mg <- length(ag.low)
  A1 <- tibble(age = a, exposure = obs.e)
  B  <- tibble(age = ag.low, value = rep(1, mg))
  C20  <- left_join(A1, B, by = "age") |> 
    mutate(ag = age * value) |> 
    fill(ag,.direction = "down") |> 
    select(-value) |> 
    pivot_wider(names_from=age, values_from = exposure, values_fill = 0) |> 
    column_to_rownames("ag") |> 
    as.matrix()

  # 
  # TR: note to self, create G20 more efficiently
  # rows age groups
  # cols single ages
  # cells exposures
  
  
  #all(colSums(G20)==1)
  ## create C matrix with e20
  # C20 <- G20
  # C20[G20==1] <- c(obs.e)
  
  ## design matrix
  Ba <- MortSmooth_bbase(x=a, min(a), max(a), floor(m/7), 3)
  Ba <- zapsmall(Ba, 8)
  
  nb <- ncol(Ba)
  U <- cbind(1, Ba)
  kappa <- 0
  ## penalty stuff
  D <- diff(diag(nb), diff=2)
  tDD <- t(D)%*%D
  
  ## constraining \delta to sum up to 0
  H0 <- matrix(c(0, rep(1,m)), 1, m+1)
  H1 <- adiag(0, Ba)
  H <- H0 %*% H1
  
  P0 <- lambda*tDD
  P <- matrix(0,nb+1,nb+1)
  P[-1,-1] <- P0
  
  ## starting eta (log-mortality minus forecast log-mortality in 2020)
  eta <- rep(0.01,m)
  ## PCLM regression with eta20 as offset
  max.it <- 100
  for(it in 1:max.it){
    gamma   <- exp(eta + eta.for)
    mu      <- c(C20 %*% gamma)
    X       <- (C20 * ((1 / mu) %*% t(gamma)) ) %*% U
    w       <- as.vector(mu)
    r       <- obs.y - mu + C20 %*% (gamma * eta)
    tXWX    <- t(X) %*% (w * X) 
    tXWXpP  <- tXWX + P
    tXr     <- t(X) %*% r
    ## adding constraints
    LHS     <- rbind(cbind(tXWXpP, t(H)),
                     cbind(H, 0))
    RHS     <- matrix(c(tXr, kappa), ncol=1)
    coeff   <- solve(LHS, RHS)
    ##
    betas   <- coeff[1:(nb+1)]
    eta.old <- eta
    eta     <- U%*%betas
    dif.eta <- max(abs((eta - eta.old)/eta.old) )
    if(dif.eta < 1e-04 & it > 4) break
    #cat(it, dif.eta, "\n")
  }
  gamma   <- exp(eta + eta.for)
  mu      <- c(C20 %*% gamma)
  ## computing BIC
  Pr <- 10^-6 * diag(nb+1)
  Hat <- solve(tXWXpP+Pr, tXWX)
  ED <- sum(diag(Hat))
  y1 <- obs.y
  y1[obs.y == 0] <- 10^(-4)
  DEV <- 2 * sum( (obs.y * log(y1/mu) - (y1-mu)))
  #BICs[l] <- DEV + log(mg)*ED
  ## try QIC
  psi <- DEV/(mg - ED)
  QIC <- mg + ED + mg*log(psi)
  
  ## variance-covariance matrix
  Pr <- 10^-3 * diag(nb+1)
  H0 <- solve(tXWXpP+Pr)
  Vbetas <- H0 %*% tXWX %*% H0
  
  ## returning object   
  out <- list(QIC=QIC, Vbetas=Vbetas,
              coeff=coeff, Ba=Ba, U=U)
  return(out)
}
