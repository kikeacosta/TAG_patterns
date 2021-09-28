## R-code for estimating log-mortality age-pattern 
## in two populations, in which the first 
## is simply a smooth function over age 
## and the second is the log-mortality of the first + 
## scaling factor +
## smooth age-factor
## issues: 
## only age-grouped data for input and output
## open age 85+
## !! underdetermined linear system of equation
## 1) actual data from different WM sources
## by C.G. Camarda 2021.09.27


## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## !!!! to be changed
setwd("~/WORK/TAG_patterns/")
library(MortalitySmooth)
library(magic)
library(colorspace)

## loading baseline
data1 <- read.csv("Output/GHE2019_baseline.csv", header=TRUE)
## loading 2020
data2 <- read.csv("Output/WM2020_observed.csv", header=TRUE)
## removed doubles Colombia in data2
data2 <- data2[-c(2161:2268), ]

## select pop which are in both datasets
## populations
pop <- as.character(sort(unique(data2$iso3)))
p <- length(pop)


## ages
x <- c(seq(2.5, 82.5, 5), 90)
m <- length(x)
x.low <- unique(data2$age)
x.up <- c(unique(data2$age)[-1], 105)
x.lab <- c(paste(x.low[-m], x.up[-m]-1,sep="-"), "85+")

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

colpop <- rainbow_hcl(p)
colpopT <- adjustcolor(colpop, 0.3)


OUT <- list()
PLOT <- FALSE
j=59
for(j in 1:p){
  OUT.j <- list()
  
  data1.j <- subset(data1, iso3==pop[j])
  data2.j <- subset(data2, iso3==pop[j])
  
  data1F.j <- subset(data1.j, sex=="Female")
  data2F.j <- subset(data2.j, sex=="Female")
  data1M.j <- subset(data1.j, sex=="Male")
  data2M.j <- subset(data2.j, sex=="Male")
  
  y1F <- data1F.j$deaths
  e1F <- data1F.j$Nx
  lmx1F <- log(y1F/e1F)
  y2F <- data2F.j$deaths
  e2F <- data2F.j$Nx
  lmx2F <- log(y2F/e2F)
  
  y1M <- data1M.j$deaths
  e1M <- data1M.j$Nx
  lmx1M <- log(y1M/e1M)
  y2M <- data2M.j$deaths
  e2M <- data2M.j$Nx
  lmx2M <- log(y2M/e2M)
  
  OUT.j$lmx1F <- lmx1F
  OUT.j$lmx2F <- lmx2F
  OUT.j$lmx1M <- lmx1M
  OUT.j$lmx2M <- lmx2M
  
  ## plotting
  if(PLOT){
    par(mfrow=c(1,2))
    rany <- range(lmx1F, lmx2F, lmx1M, lmx2M, finite=TRUE)
    ranx <- range(x.low, x.up)
    
    plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE,
         xlab="age", ylab="mortality, log-scale",
         main="mortality age-pattern, Females")
    yy <- 10^seq(-7, 2)
    axis(2, at=log(yy), labels=yy)
    axis(1);box()
    for(i in 1:m){
      segments(x0=x.low[i], x1=x.up[i]+1, y0=lmx1F[i], y1=lmx1F[i], col=col1F, lwd=3)
      segments(x0=x.low[i], x1=x.up[i]+1, y0=lmx2F[i], y1=lmx2F[i], col=col2F, lwd=3)
    }
    legend("topleft", inset=0.1,
           legend=c("2019", "2020"),
           col=c(col1F,col2F), lwd=3)

    plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE,
         xlab="age", ylab="mortality, log-scale",
         main="mortality age-pattern, Males")
    yy <- 10^seq(-7, 2)
    axis(2, at=log(yy), labels=yy)
    axis(1);box()
    for(i in 1:m){
      segments(x0=x.low[i], x1=x.up[i]+1, y0=lmx1M[i], y1=lmx1M[i], col=col1M, lwd=3)
      segments(x0=x.low[i], x1=x.up[i]+1, y0=lmx2M[i], y1=lmx2M[i], col=col2M, lwd=3)
    }
    legend("topleft", inset=0.1,
           legend=c("2019", "2020"),
           col=c(col1M,col2M), lwd=3)
  }
  
  ## FEMALES
  ## estimation
  ## response
  y <- c(y1F, y2F)
  ## offset
  e <- c(e1F, e2F)
  ## finding starting values
  eta.st <- c(lmx1F, lmx2F)
  
  ## B-splines
  B <- MortSmooth_bbase(x=x, xl=min(x), xr=max(x), ndx=8, deg=3)
  nb <- ncol(B)
  
  ## penalty stuff
  Deta1 <- diff(diag(nb), diff=2)
  tDDeta1 <- t(Deta1)%*%Deta1
  Ddelta <- diff(diag(nb), diff=2)
  tDDdelta <- t(Ddelta)%*%Ddelta
  
  ## penalty stuff
  lambda.eta1 <- 10^-1
  lambda.delta <- 10^2
  P <- adiag(lambda.eta1*tDDeta1,
             0,
             lambda.delta*tDDdelta)
  ## ridge penalty
  Pr <- 1e-4*diag(ncol(P))
  
  ## model matrix
  X0 <- rbind(B, B)
  X1 <- matrix(c(rep(0,m), rep(1,m)), 2*m)
  X2 <- rbind(0*B, B)
  X <- cbind(X0,X1,X2)
  ## constraining delta to sum up to 0
  H <- matrix(c(rep(0,nb+1), rep(1,nb)), nrow=1)
  kappa <- 0
  
  eta <- eta.st
  max.it <- 100
  for(it in 1:max.it){
    mu <- exp(eta)
    z <- eta + (y - e*mu)/(e*mu)
    w <- as.vector(e*mu)
    tXWX <- t(X)%*%(w*X)
    tXWXpP <- tXWX + P + Pr
    tXWz <- t(X)%*%(w*z)
    
    LHS <- rbind(cbind(tXWXpP, t(H)),
                 cbind(H, 0))
    RHS <- matrix(c(tXWz, kappa), ncol=1)
    coeff   <- solve(LHS, RHS)
    betas   <- coeff[1:(nb*2+1)]
    eta.old <- eta
    eta     <- X%*%betas
    dif.eta <- max(abs((eta - eta.old)/eta.old))
    if(dif.eta < 1e-04 & it > 4) break
    #cat(it, dif.eta, "\n")
  }
  betasF.hat <- betas
  eta1F.hat <- B%*%betas[1:nb]
  cF.hat <- betas[nb+1]
  deltaF.hat <- B%*%betas[1:nb+1+nb]
  eta2F.hat <- eta1F.hat + cF.hat + deltaF.hat
  
  ## standard errors
  H0 <- solve(tXWXpP)
  Vbetas <- H0 %*% tXWX %*% H0
  diagVbetas <- diag(Vbetas)
  se.betas <- sqrt(diagVbetas)
  se.betas[is.na(se.betas)] <- mean(se.betas[!is.na(se.betas)])
  se.c <- se.betas[nb+1]
  V.eta12 <- X %*% Vbetas %*% t(X)
  se.eta1 <- sqrt(diag(V.eta12))[1:m]
  se.eta2 <- sqrt(diag(V.eta12))[1:m+m]
  V.delta <- B %*% Vbetas[1:nb+nb+1, 1:nb+nb+1] %*% t(B)
  se.delta <- sqrt(diag(V.delta))
  se.delta[is.na(se.delta)] <- mean(se.delta[!is.na(se.delta)])
  
  ## confidence intervals
  eta1F.hatL <- eta1F.hat - 2*se.eta1
  eta1F.hatU <- eta1F.hat + 2*se.eta1
  eta2F.hatL <- eta2F.hat - 2*se.eta2
  eta2F.hatU <- eta2F.hat + 2*se.eta2
  cF.hatL <- cF.hat - 2*se.c
  cF.hatU <- cF.hat + 2*se.c
  deltaF.hatL <- deltaF.hat - 2*se.delta
  deltaF.hatU <- deltaF.hat + 2*se.delta
  
  
  
  ## MALES
  ## estimation
  ## response
  y <- c(y1M, y2M)
  ## offset
  e <- c(e1M, e2M)
  ## finding starting values
  eta.st <- c(lmx1M, lmx2M)
  
  ## B-splines
  B <- MortSmooth_bbase(x=x, xl=min(x), xr=max(x), ndx=8, deg=3)
  nb <- ncol(B)
  
  ## penalty stuff
  Deta1 <- diff(diag(nb), diff=2)
  tDDeta1 <- t(Deta1)%*%Deta1
  Ddelta <- diff(diag(nb), diff=2)
  tDDdelta <- t(Ddelta)%*%Ddelta
  
  ## penalty stuff
  lambda.eta1 <- 10^-1
  lambda.delta <- 10^2
  P <- adiag(lambda.eta1*tDDeta1,
             0,
             lambda.delta*tDDdelta)
  ## ridge penalty
  Pr <- 1e-4*diag(ncol(P))
  
  ## model matrix
  X0 <- rbind(B, B)
  X1 <- matrix(c(rep(0,m), rep(1,m)), 2*m)
  X2 <- rbind(0*B, B)
  X <- cbind(X0,X1,X2)
  ## constraining delta to sum up to 0
  H <- matrix(c(rep(0,nb+1), rep(1,nb)), nrow=1)
  kappa <- 0
  
  eta <- eta.st
  max.it <- 100
  for(it in 1:max.it){
    mu <- exp(eta)
    z <- eta + (y - e*mu)/(e*mu)
    w <- as.vector(e*mu)
    tXWX <- t(X)%*%(w*X)
    tXWXpP <- tXWX + P + Pr
    tXWz <- t(X)%*%(w*z)
    
    LHS <- rbind(cbind(tXWXpP, t(H)),
                 cbind(H, 0))
    RHS <- matrix(c(tXWz, kappa), ncol=1)
    coeff   <- solve(LHS, RHS)
    betas   <- coeff[1:(nb*2+1)]
    eta.old <- eta
    eta     <- X%*%betas
    dif.eta <- max(abs((eta - eta.old)/eta.old))
    if(dif.eta < 1e-04 & it > 4) break
    #cat(it, dif.eta, "\n")
  }
  betasM.hat <- betas
  eta1M.hat <- B%*%betas[1:nb]
  cM.hat <- betas[nb+1]
  deltaM.hat <- B%*%betas[1:nb+1+nb]
  eta2M.hat <- eta1M.hat + cM.hat + deltaM.hat
  
  ## standard errors
  H0 <- solve(tXWXpP)
  Vbetas <- H0 %*% tXWX %*% H0
  diagVbetas <- diag(Vbetas)
  se.betas <- sqrt(diagVbetas)
  se.betas[is.na(se.betas)] <- mean(se.betas[!is.na(se.betas)])
  se.c <- se.betas[nb+1]
  V.eta12 <- X %*% Vbetas %*% t(X)
  se.eta1 <- sqrt(diag(V.eta12))[1:m]
  se.eta2 <- sqrt(diag(V.eta12))[1:m+m]
  V.delta <- B %*% Vbetas[1:nb+nb+1, 1:nb+nb+1] %*% t(B)
  se.delta <- sqrt(diag(V.delta))
  se.delta[is.na(se.delta)] <- mean(se.delta[!is.na(se.delta)])
  ## confidence intervals
  eta1M.hatL <- eta1M.hat - 2*se.eta1
  eta1M.hatU <- eta1M.hat + 2*se.eta1
  eta2M.hatL <- eta2M.hat - 2*se.eta2
  eta2M.hatU <- eta2M.hat + 2*se.eta2
  cM.hatL <- cM.hat - 2*se.c
  cM.hatU <- cM.hat + 2*se.c
  deltaM.hatL <- deltaM.hat - 2*se.delta
  deltaM.hatU <- deltaM.hat + 2*se.delta
  
  
  OUT.j$betasF.hat <- betasF.hat
  OUT.j$betasM.hat <- betasM.hat
  
  OUT.j$eta1F.hat <- eta1F.hat
  OUT.j$eta1F.hatL <- eta1F.hatL
  OUT.j$eta1F.hatU <- eta1F.hatU
  OUT.j$eta2F.hat <- eta2F.hat
  OUT.j$eta2F.hatL <- eta2F.hatL
  OUT.j$eta2F.hatU <- eta2F.hatU
  
  OUT.j$eta1M.hat <- eta1M.hat
  OUT.j$eta1M.hatL <- eta1M.hatL
  OUT.j$eta1M.hatU <- eta1M.hatU
  OUT.j$eta2M.hat <- eta2M.hat
  OUT.j$eta2M.hatL <- eta2M.hatL
  OUT.j$eta2M.hatU <- eta2M.hatU
  
  OUT.j$deltaF.hat <- deltaF.hat
  OUT.j$deltaF.hatL <- deltaF.hatL
  OUT.j$deltaF.hatU <- deltaF.hatU
  OUT.j$cF.hat <- cF.hat
  OUT.j$cF.hatL <- cF.hatL
  OUT.j$cF.hatU <- cF.hatU
  
  OUT.j$deltaM.hat <- deltaM.hat
  OUT.j$deltaM.hatL <- deltaM.hatL
  OUT.j$deltaM.hatU <- deltaM.hatU
  OUT.j$cM.hat <- cM.hat
  OUT.j$cM.hatL <- cM.hatL
  OUT.j$cM.hatU <- cM.hatU
  
  
  OUT[[j]] <- OUT.j
  
  if(PLOT){
    par(mfrow=c(2,3))
    
    ## females
    rany <- range(lmx1F, lmx2F, lmx1M, lmx2M, 
                  lmx1F, lmx2F, lmx1M, lmx2M, finite=TRUE)
    ranx <- range(x)
    plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE,
         xlab="age", ylab="mortality, log-scale",
         main="mortality age-pattern, Females, 2019")
    legend("topleft", legend=paste(pop[j]), cex=2.5, text.col=8, bty="n")
    yy <- 10^seq(-7, 2)
    axis(2, at=log(yy), labels=yy)
    axis(1);box()
    points(x, lmx1F, col=col1FT, pch=16)
    points(x, eta1F.hat, col=col1F, pch=3)
    arrows(x0=x, x1=x, y0=eta1F.hatL, y1=eta1F.hatU, col=col1F, length=0.1, angle=90, code=3)
    
    plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE,
         xlab="age", ylab="mortality, log-scale",
         main="mortality age-pattern, Females, 2020")
    yy <- 10^seq(-7, 2)
    axis(2, at=log(yy), labels=yy)
    axis(1);box()
    points(x, lmx2F, col=col2FT, pch=16)
    points(x, eta2F.hat, col=col2F, pch=3)
    arrows(x0=x, x1=x, y0=eta2F.hatL, y1=eta2F.hatU, col=col2F, length=0.1, angle=90, code=3)
    
    ranx <- range(-10, x)
    rany <- range(deltaF.hatL,deltaF.hatU, cF.hatU, cF.hatL, 
                  na.rm = TRUE)
    plot(1, 1, t="n", xlim=ranx, ylim=rany,
         main="", axes=FALSE,
         xlab="ages", ylab="delta")
    axis(2)
    axis(1);box()
    abline(h=0, col=8, lty=3, lwd=2)
    points(x, deltaF.hat, col=coldF, pch=3)
    arrows(x0=x, x1=x, y0=deltaF.hatL, y1=deltaF.hatU, col=coldF, length=0.1, angle=90, code=3)
    abline(v=-1, col=8, lwd=4, lty=2)
    points(-8, cF.hat, col=colcF, pch=3, lwd=3)
    arrows(x0=-8, x1=-8, y0=cF.hatL, y1=cF.hatU, col=colcF, lwd=3, 
           angle = 90, code=3, length=0.15)
    
    ## Males
    rany <- range(lmx1F, lmx2F, lmx1M, lmx2M, 
                  lmx1F, lmx2F, lmx1M, lmx2M, finite=TRUE)
    ranx <- range(x)
    plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE,
         xlab="age", ylab="mortality, log-scale",
         main="mortality age-pattern, Males, 2019")
    yy <- 10^seq(-7, 2)
    axis(2, at=log(yy), labels=yy)
    axis(1);box()
    points(x, lmx1M, col=col1MT, pch=16)
    points(x, eta1M.hat, col=col1M, pch=3)
    arrows(x0=x, x1=x, y0=eta1M.hatL, y1=eta1M.hatU, col=col1M, length=0.1, angle=90, code=3)
    
    plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE,
         xlab="age", ylab="mortality, log-scale",
         main="mortality age-pattern, Males, 2020")
    yy <- 10^seq(-7, 2)
    axis(2, at=log(yy), labels=yy)
    axis(1);box()
    points(x, lmx2M, col=col2MT, pch=16)
    points(x, eta2M.hat, col=col2M, pch=3)
    arrows(x0=x, x1=x, y0=eta2M.hatL, y1=eta2M.hatU, col=col2M, length=0.1, angle=90, code=3)
    
    ranx <- range(-10, x)
    rany <- range(deltaM.hat, deltaF.hat, deltaM.hatL,deltaM.hatU, deltaF.hatL,deltaF.hatU, cM.hatU, cM.hatL, 0,
                  na.rm = TRUE)
    plot(1, 1, t="n", xlim=ranx, ylim=rany,
         main="", axes=FALSE,
         xlab="ages", ylab="delta")
    axis(2)
    axis(1);box()
    abline(h=0, col=8, lty=3, lwd=2)
    points(x, deltaM.hat, col=coldM, pch=3)
    arrows(x0=x, x1=x, y0=deltaM.hatL, y1=deltaM.hatU, col=coldM, length=0.1, angle=90, code=3)
    
    abline(v=-1, col=8, lwd=4, lty=2)
    points(-8, cM.hat, col=colcM, pch=3, lwd=3)
    arrows(x0=-8, x1=-8, y0=cM.hatL, y1=cM.hatU, col=colcM, lwd=3, 
           angle = 90, code=3, length=0.15)
    locator(1)
  }
} 
names(OUT) <- pop

save.image("Output/OutPrepandemic2020sex_StratifiedBySex_WMdata.Rdata")


DELTAS_WM <- expand.grid(age=seq(0,85, 5), sex=c("Female", "Male"), iso3=pop)
DELTAS_WM$deltas <- NA
j=1
for(j in 1:p){
  OUT.j <- OUT[[j]]
  deltas <- c(OUT.j$deltaF.hat, OUT.j$deltaM.hat)
  wr <- 1:(m*2)+(j-1)*(m*2)
  DELTAS_WM$deltas[wr] <- deltas
}
write.table(DELTAS_WM, "DELTAS_WM.txt")




