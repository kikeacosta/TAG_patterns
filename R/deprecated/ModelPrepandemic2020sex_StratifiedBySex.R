## R-code for estimating log-mortality age-pattern 
## in two populations, in which the first 
## is simply a smooth function over age 
## and the second is the log-mortality of the first + 
## scaling factor +
## smooth age-factor
## issues: 
## grouping structures which could be eventually be different
## !! underdetermined linear system of equation
## 1) actual data from different sources
## by C.G. Camarda 2021.07.07


## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## !!!! to be changed
setwd("~/WORK/TAG_patterns/")
library(MortalitySmooth)
library(magic)
library(colorspace)

## loading deaths
YY0 <- read.csv("Output/annual_deaths_countries_selected_sources.csv", header=TRUE)
## remove totals
YY1 <- subset(YY0, Sex!="t")

## loading exposures
EE0 <- read.csv("Output/offsets.csv", header=TRUE)
## remove only totals
EE1 <- subset(EE0, Sex!="t")
## remove 2021
EE1 <- subset(EE1, Year<=2020)

## remove uncommon countries
## "Liechtenstein" in YY1
## "Ecuador" in EE1
YY1 <- subset(YY1, Country!="Liechtenstein")
EE1 <- subset(EE1, Country!="Ecuador")
rowSums(table(YY1$Source, YY1$Country)!=0)
## manual checking
cbind(as.character(sort(unique(YY1$Country))),
      as.character(sort(unique(EE1$Country))))

## sorting by name
YY <- YY1[order(YY1$Country),]
YYF <- subset(YY, Sex=="f")
YYM <- subset(YY, Sex=="m")
EE <- EE1[order(EE1$Country),]
EEF <- subset(EE, Sex=="f")
EEM <- subset(EE, Sex=="m")

# YY2020 <- subset(YY, Year==2020)
# tapply(YY2020$Deaths, YY2020$Source, length)
# colSums(table(YY2020$Country, YY2020$Source))


## populations
pop <- as.character(sort(unique(YY1$Country)))
p <- length(pop)
##
x <- unique(EE1$Age)
m <- length(x)
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


colcou <- rainbow_hcl(p)
colcouT <- adjustcolor(colcou, 0.3)


OUT <- list()
PLOT <- FALSE
j=21
for(j in 1:p){
  OUT.j <- list()
  ## subset for a given country
  YYF.i <- subset(YYF, Country==pop[j])
  EEF.i <- subset(EEF, Country==pop[j])
  YYM.i <- subset(YYM, Country==pop[j])
  EEM.i <- subset(EEM, Country==pop[j])
  
  ## pre-pandemic years
  ## deaths
  Y1F <- subset(YYF.i, Year<2020)
  t1 <- unique(Y1F$Year)
  nt1 <- length(t1)
  y1F <- rowSums(matrix(Y1F$Deaths, ncol=nt1))
  Y1M <- subset(YYM.i, Year<2020)
  y1M <- rowSums(matrix(Y1M$Deaths, ncol=nt1))
  ## exposures
  E1F <- subset(EEF.i, Year%in%t1)
  e1F <- rowSums(matrix(E1F$Population, ncol=nt1))
  E1M <- subset(EEM.i, Year%in%t1)
  e1M <- rowSums(matrix(E1M$Population, ncol=nt1))
  
  ## grouping structure
  low1 <- unique(Y1F$Age)
  up1  <- c(unique(Y1F$Age)[-1]-1, 100)
  n1 <- length(low1)
  age.gr1 <- paste(low1,up1,sep="-")
  len1 <- up1-low1+1
  ## matrix to group exposures (for plotting purposes)
  ## and composite matrix
  G1 <- matrix(0, n1, m)
  rownames(G1) <- age.gr1
  colnames(G1) <- x
  C1F <- C1M <- G1
  for(i in 1:n1){
    age.low.i <- low1[i]
    age.up.i <- up1[i]
    whi <- which(x>=age.low.i & x<=age.up.i)
    G1[i,whi] <- 1
    C1F[i,whi] <- e1F[whi]
    C1M[i,whi] <- e1M[whi]
  }
  e1Fg <- G1%*%e1F
  lmx1Fg <- log(y1F/e1Fg)
  e1Mg <- G1%*%e1M
  lmx1Mg <- log(y1M/e1Mg)
  
  ## 2020
  ## deaths
  Y2F <- subset(YYF.i, Year==2020)
  y2F <- Y2F$Deaths
  Y2M <- subset(YYM.i, Year==2020)
  y2M <- Y2M$Deaths
  ## exposures
  e2F <- subset(EEF.i, Year==2020)$Population
  e2M <- subset(EEM.i, Year==2020)$Population
  ## grouping structure
  low2 <- unique(Y2F$Age)
  up2  <- c(unique(Y2F$Age)[-1]-1, 100)
  n2 <- length(low2)
  age.gr2 <- paste(low2,up2,sep="-")
  len2 <- up2-low2+1
  ## matrix to group exposures (for plotting purposes)
  ## and composite matrix
  G2 <- matrix(0, n2, m)
  rownames(G2) <- age.gr2
  colnames(G2) <- x
  C2F <- C2M <- G2
  for(i in 1:n2){
    age.low.i <- low2[i]
    age.up.i <- up2[i]
    whi <- which(x>=age.low.i & x<=age.up.i)
    G2[i,whi] <- 1
    C2F[i,whi] <- e2F[whi]
    C2M[i,whi] <- e2M[whi]
  }
  e2Fg <- G2%*%e2F
  lmx2Fg <- log(y2F/e2Fg)
  e2Mg <- G2%*%e2M
  lmx2Mg <- log(y2M/e2Mg)
  
  
  OUT.j$t1 <- t1
  OUT.j$low1 <- low1
  OUT.j$up1 <- up1
  OUT.j$low2 <- low2
  OUT.j$up2 <- up2
  OUT.j$lmx1Fg <- lmx1Fg
  OUT.j$lmx1Mg <- lmx1Mg
  OUT.j$lmx2Fg <- lmx2Fg
  OUT.j$lmx2Mg <- lmx2Mg
  
  # ## plotting
  # if(PLOT){
  #   rany <- range(lmx1g, lmx2g, finite=TRUE)
  #   ranx <- range(x)
  #   plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE,
  #        xlab="age", ylab="mortality, log-scale",
  #        main="mortality age-pattern")
  #   yy <- 10^seq(-7, 2)
  #   axis(2, at=log(yy), labels=yy)
  #   axis(1);box()
  #   for(i in 1:n1){
  #     segments(x0=low1[i], x1=up1[i], y0=lmx1g[i], y1=lmx1g[i], col=col1, lwd=3)
  #   }
  #   for(i in 1:n2){
  #     segments(x0=low2[i], x1=up2[i], y0=lmx2g[i], y1=lmx2g[i], col=col2, lwd=3)
  #   }
  #   legend("topleft", inset=0.1,
  #          legend=c(paste(min(t1), max(t1), sep="-"), "2020"),
  #          col=c(col1,col2), lwd=3)
  # }
  
  ## FEMALES
  ## estimation
  ## response
  y <- c(y1F, y2F)
  ## overall composite matrix
  C <- adiag(C1F, C2F)
  
  ## finding starting values
  y1F.st0 <- rep(y1F/len1, len1)
  wei <- rep(1,m)
  wei[e1F==0] <- 0
  fit1F.0 <- Mort1Dsmooth(x=x, y=y1F.st0,
                          offset=log(e1F),
                          w=wei,
                          method=3, lambda=10^4)
  y2F.st0 <- rep(y2F/len2, len2)
  wei <- rep(1,m)
  wei[e2F==0] <- 0
  fit2F.0 <- Mort1Dsmooth(x=x, y=y2F.st0,
                          offset=log(e2F),
                          w=wei,
                          method=3, lambda=10^4)
  
  
  eta.st <- c(fit1F.0$logmortality, 
              fit2F.0$logmortality)
  
  ## penalty stuff
  Deta1 <- diff(diag(m), diff=2)
  tDDeta1 <- t(Deta1)%*%Deta1
  Ddelta <- diff(diag(m), diff=2)
  tDDdelta <- t(Ddelta)%*%Ddelta
  
  lambda.eta1 <- 10^3
  lambda.delta <- 10^5
  
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
  ## constraining delta to sum up to 0
  H <- matrix(c(rep(0,m+1), rep(1,m)), nrow=1)
  kappa <- 0
  
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
    ## adding constraints
    LHS <- rbind(cbind(GpP, t(H)),
                 cbind(H, 0))
    RHS <- matrix(c(tXr, kappa), ncol=1)
    coeff   <- solve(LHS, RHS)
    betas   <- coeff[1:(m*2+1)]
    #betas   <- solve(GpP, tXr)
    eta.old <- eta
    eta     <- U%*%betas
    dif.eta <- max(abs((eta - eta.old)/eta.old))
    if(dif.eta < 1e-04 & it > 4) break
    #cat(it, dif.eta, "\n")
  }
  eta1F.hat <- betas[1:m]
  cF.hat <- betas[m+1]
  deltaF.hat <- betas[1:m+1+m]
  eta2F.hat <- eta1F.hat + cF.hat + deltaF.hat
  ## standard errors
  H0 <- solve(GpP)
  Vbetas <- H0 %*% G %*% H0
  diagVbetas <- diag(Vbetas)
  se.betas <- sqrt(diagVbetas)
  se.eta1 <- se.betas[1:m]
  se.c <- se.betas[m+1]
  se.delta <- se.betas[1:m+1+m]
  se.delta[is.na(se.delta)] <- mean(se.delta[!is.na(se.delta)])
  V.eta12 <- U %*% Vbetas %*% t(U)
  se.eta2 <- sqrt(diag(V.eta12))[1:m+m]
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
  ## overall composite matrix
  C <- adiag(C1M, C2M)
  
  ## finding starting values
  y1M.st0 <- rep(y1M/len1, len1)
  wei <- rep(1,m)
  wei[e1M==0] <- 0
  fit1M.0 <- Mort1Dsmooth(x=x, y=y1M.st0,
                          offset=log(e1M),
                          w=wei,
                          method=3, lambda=10^4)
  y2M.st0 <- rep(y2M/len2, len2)
  wei <- rep(1,m)
  wei[e2M==0] <- 0
  fit2M.0 <- Mort1Dsmooth(x=x, y=y2M.st0,
                          offset=log(e2M),
                          w=wei,
                          method=3, lambda=10^4)
  
  
  eta.st <- c(fit1M.0$logmortality, 
              fit2M.0$logmortality)
  
  ## penalty stuff
  Deta1 <- diff(diag(m), diff=2)
  tDDeta1 <- t(Deta1)%*%Deta1
  Ddelta <- diff(diag(m), diff=2)
  tDDdelta <- t(Ddelta)%*%Ddelta
  
  lambda.eta1 <- 10^3
  lambda.delta <- 10^5
  
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
  ## constraining delta to sum up to 0
  H <- matrix(c(rep(0,m+1), rep(1,m)), nrow=1)
  kappa <- 0
  
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
    ## adding constraints
    LHS <- rbind(cbind(GpP, t(H)),
                 cbind(H, 0))
    RHS <- matrix(c(tXr, kappa), ncol=1)
    coeff   <- solve(LHS, RHS)
    betas   <- coeff[1:(m*2+1)]
    #betas   <- solve(GpP, tXr)
    eta.old <- eta
    eta     <- U%*%betas
    dif.eta <- max(abs((eta - eta.old)/eta.old))
    if(dif.eta < 1e-04 & it > 4) break
    #cat(it, dif.eta, "\n")
  }
  eta1M.hat <- betas[1:m]
  cM.hat <- betas[m+1]
  deltaM.hat <- betas[1:m+1+m]
  eta2M.hat <- eta1M.hat + cM.hat + deltaM.hat
  ## standard errors
  H0 <- solve(GpP)
  Vbetas <- H0 %*% G %*% H0
  diagVbetas <- diag(Vbetas)
  se.betas <- sqrt(diagVbetas)
  se.eta1 <- se.betas[1:m]
  se.c <- se.betas[m+1]
  se.delta <- se.betas[1:m+1+m]
  se.delta[is.na(se.delta)] <- mean(se.delta[!is.na(se.delta)])
  V.eta12 <- U %*% Vbetas %*% t(U)
  se.eta2 <- sqrt(diag(V.eta12))[1:m+m]
  ## confidence intervals
  eta1M.hatL <- eta1M.hat - 2*se.eta1
  eta1M.hatU <- eta1M.hat + 2*se.eta1
  eta2M.hatL <- eta2M.hat - 2*se.eta2
  eta2M.hatU <- eta2M.hat + 2*se.eta2
  cM.hatL <- cM.hat - 2*se.c
  cM.hatU <- cM.hat + 2*se.c
  deltaM.hatL <- deltaM.hat - 2*se.delta
  deltaM.hatU <- deltaM.hat + 2*se.delta
  
  
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
    par(mfrow=c(1,3))
    ## eta1 + eta2
    rany <- range(lmx1Fg, lmx2Fg, lmx1Mg, lmx2Mg,
                  eta1F.hatL, eta1F.hatU, eta2F.hatL, eta2F.hatU, 
                  eta1M.hatL, eta1M.hatU, eta2M.hatL, eta2M.hatU, 
                  finite=TRUE)
    ranx <- range(x)
    plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE, 
         xlab="age", ylab="mortality, log-scale",
         main="mortality age-pattern, females")
    yy <- 10^seq(-7, 2)
    axis(2, at=log(yy), labels=yy)
    axis(1);box()
    ## Females
    for(i in 1:n1){
      segments(x0=low1[i], x1=up1[i], y0=lmx1Fg[i], y1=lmx1Fg[i], col=col1F, lwd=1)
    }
    for(i in 1:n2){
      segments(x0=low1[i], x1=up1[i], y0=lmx2Fg[i], y1=lmx2Fg[i], col=col2F, lwd=1)
    }
    xx <- c(x, rev(x))
    yy <- c(eta1F.hatL, rev(eta1F.hatU))
    polygon(xx, yy, col=col1FT, border=col1FT)
    lines(x, eta1F.hat, col=col1F, lwd=2)
    yy <- c(eta2F.hatL, rev(eta2F.hatU))
    polygon(xx, yy, col=col2FT, border=col2FT)
    lines(x, eta2F.hat, col=col2F, lwd=2)
    legend("topleft", inset=0.1,
           legend=c(paste(min(t1), max(t1), sep="-"), "2020"),
           col=c(col1F,col2F), lwd=3)
    
    plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE, 
         xlab="age", ylab="mortality, log-scale",
         main="mortality age-pattern, males")
    yy <- 10^seq(-7, 2)
    axis(2, at=log(yy), labels=yy)
    axis(1);box()
    ## Males
    for(i in 1:n1){
      segments(x0=low1[i], x1=up1[i], y0=lmx1Mg[i], y1=lmx1Mg[i], col=col1M, lwd=1)
    }
    for(i in 1:n2){
      segments(x0=low1[i], x1=up1[i], y0=lmx2Mg[i], y1=lmx2Mg[i], col=col2M, lwd=1)
    }
    xx <- c(x, rev(x))
    yy <- c(eta1M.hatL, rev(eta1M.hatU))
    polygon(xx, yy, col=col1MT, border=col1MT)
    lines(x, eta1M.hat, col=col1M, lwd=2)
    yy <- c(eta2M.hatL, rev(eta2M.hatU))
    polygon(xx, yy, col=col2MT, border=col2MT)
    lines(x, eta2M.hat, col=col2M, lwd=2)
    legend("topleft", inset=0.1,
           legend=c(paste(min(t1), max(t1), sep="-"), "2020"),
           col=c(col1M,col2M), lwd=3)

    ## deltas
    ranx <- range(-10, x)
    rany <- range(deltaF.hatL,deltaF.hatU, cF.hatU, cF.hatL, 
                  deltaM.hatL,deltaM.hatU, cM.hatU, cM.hatL,
                  na.rm = TRUE)
    plot(1, 1, t="n", xlim=ranx, ylim=rany,
         main="", axes=FALSE,
         xlab="ages", ylab="delta")
    axis(2)
    axis(1);box()
    abline(h=0, col=8, lty=3, lwd=2)
    xx <- c(x, rev(x))
    yy <- c(deltaF.hatL, rev(deltaF.hatU))
    polygon(xx, yy, col=coldFT, border=coldFT)
    lines(x, deltaF.hat, col=coldF, lwd=4)
    yy <- c(deltaM.hatL, rev(deltaM.hatU))
    polygon(xx, yy, col=coldMT, border=coldMT)
    lines(x, deltaM.hat, col=coldM, lwd=4)
    
    abline(v=-1, col=8, lwd=4, lty=2)
    points(-8, cF.hat, col=colcF, pch=3, lwd=3)
    arrows(x0=-8, x1=-8, y0=cF.hatL, y1=cF.hatU, col=colcF, lwd=3, 
           angle = 90, code=3, length=0.15)
    points(-6, cM.hat, col=colcM, pch=3, lwd=3)
    arrows(x0=-6, x1=-6, y0=cM.hatL, y1=cM.hatU, col=colcM, lwd=3, 
           angle = 90, code=3, length=0.15)
    mtext("c", 1, at=-7, cex=2, line=1.5)
    legend("top", inset=0.1,
           legend=c("Females", "Males"),
           col=c(colcF,colcM), lwd=3)
    locator(1)
  }
}
names(OUT) <- pop

save.image("Output/OutPrepandemic2020sex_StratifiedBySex.Rdata")
