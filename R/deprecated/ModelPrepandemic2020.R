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
## take only totals
YY1 <- subset(YY0, Sex=="t")

## loading exposures
EE0 <- read.csv("Output/offsets.csv", header=TRUE)
## take only totals
EE1 <- subset(EE0, Sex=="t")
## remove 2021
EE1 <- subset(EE1, Year<=2020)

## remove uncommon countries
## "Liechtenstein" in Y1
YY1 <- subset(YY1, Country!="Liechtenstein")
rowSums(table(YY1$Source, YY1$Country)!=0)
## manual checking
cbind(as.character(sort(unique(YY1$Country))),
      as.character(sort(unique(EE1$Country))))

## sorting by name
YY <- YY1[order(YY1$Country),]
EE <- EE1[order(EE1$Country),]

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
col1 <- "red"
col1T <- adjustcolor(col1, 0.5)
col2 <- "blue"
col2T <- adjustcolor(col2, 0.5)
cold <- "orange"
coldT <- adjustcolor(cold, 0.5)
colc <- "violet"
colcT <- adjustcolor(colc, 0.1)
colcou <- rainbow_hcl(p)
colcouT <- adjustcolor(colcou, 0.3)


OUT <- list()
PLOT <- FALSE
for(j in 1:p){
  OUT.j <- list()
  ## subset for a given country
  YY.i <- subset(YY, Country==pop[j])
  EE.i <- subset(EE, Country==pop[j])
  
  ## pre-pandemic years
  ## deaths
  Y1 <- subset(YY.i, Year<2020)
  t1 <- unique(Y1$Year)
  nt1 <- length(t1)
  y1 <- rowSums(matrix(Y1$Deaths, ncol=nt1))
  ## exposures
  E1 <- subset(EE.i, Year%in%t1)
  e1 <- rowSums(matrix(E1$Population, ncol=nt1))
  ## grouping structure
  low1 <- unique(Y1$Age)
  up1  <- c(unique(Y1$Age)[-1]-1, 100)
  n1 <- length(low1)
  age.gr1 <- paste(low1,up1,sep="-")
  len1 <- up1-low1+1
  ## matrix to group exposures (for plotting purposes)
  ## and composite matrix
  G1 <- matrix(0, n1, m)
  rownames(G1) <- age.gr1
  colnames(G1) <- x
  C1 <- G1
  for(i in 1:n1){
    age.low.i <- low1[i]
    age.up.i <- up1[i]
    whi <- which(x>=age.low.i & x<=age.up.i)
    G1[i,whi] <- 1
    C1[i,whi] <- e1[whi]
  }
  e1g <- G1%*%e1
  lmx1g <- log(y1/e1g)
  ## 2020
  ## deaths
  Y2 <- subset(YY.i, Year==2020)
  y2 <- Y2$Deaths
  ## exposures
  e2 <- subset(EE.i, Year==2020)$Population
  ## grouping structure
  low2 <- unique(Y2$Age)
  up2  <- c(unique(Y2$Age)[-1]-1, 100)
  n2 <- length(low2)
  age.gr2 <- paste(low2,up2,sep="-")
  len2 <- up2-low2+1
  ## matrix to group exposures (for plotting purposes)
  ## and composite matrix
  G2 <- matrix(0, n2, m)
  rownames(G2) <- age.gr2
  colnames(G2) <- x
  C2 <- G2
  for(i in 1:n2){
    age.low.i <- low2[i]
    age.up.i <- up2[i]
    whi <- which(x>=age.low.i & x<=age.up.i)
    G2[i,whi] <- 1
    C2[i,whi] <- e2[whi]
  }
  e2g <- G2%*%e2
  lmx2g <- log(y2/e2g)
  
  
  OUT.j$t1 <- t1
  OUT.j$low1 <- low1
  OUT.j$up1 <- up1
  OUT.j$low2 <- low2
  OUT.j$up2 <- up2
  OUT.j$lmx1g <- lmx1g
  OUT.j$lmx2g <- lmx2g
  
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
  ## estimation
  ## response
  y <- c(y1, y2)
  ## overall composite matrix
  C <- adiag(C1, C2)
  
  ## finding starting values
  y1.st0 <- rep(y1/len1, len1)
  fit1.0 <- Mort1Dsmooth(x=x, y=y1.st0,
                         offset=log(e1),
                         method=3, lambda=10^3)
  y2.st0 <- rep(y2/len2, len2)
  fit2.0 <- Mort1Dsmooth(x=x, y=y2.st0,
                         offset=log(e2),
                         method=3, lambda=10^3)
  eta.st <- c(fit1.0$logmortality, fit2.0$logmortality)
  
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
  eta1.hat <- betas[1:m]
  c.hat <- betas[m+1]
  delta.hat <- betas[1:m+1+m]
  eta2.hat <- eta1.hat + c.hat + delta.hat
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
  eta1.hatL <- eta1.hat - 2*se.eta1
  eta1.hatU <- eta1.hat + 2*se.eta1
  eta2.hatL <- eta2.hat - 2*se.eta2
  eta2.hatU <- eta2.hat + 2*se.eta2
  c.hatL <- c.hat - 2*se.c
  c.hatU <- c.hat + 2*se.c
  delta.hatL <- delta.hat - 2*se.delta
  delta.hatU <- delta.hat + 2*se.delta
  
  OUT.j$eta1.hat <- eta1.hat
  OUT.j$eta1.hatL <- eta1.hatL
  OUT.j$eta1.hatU <- eta1.hatU
  OUT.j$eta2.hat <- eta2.hat
  OUT.j$eta2.hatL <- eta2.hatL
  OUT.j$eta2.hatU <- eta2.hatU
  OUT.j$delta.hat <- delta.hat
  OUT.j$delta.hatL <- delta.hatL
  OUT.j$delta.hatU <- delta.hatU
  OUT.j$c.hat <- c.hat
  OUT.j$c.hatL <- c.hatL
  OUT.j$c.hatU <- c.hatU
  
  OUT[[j]] <- OUT.j
  
  if(PLOT){
    par(mfrow=c(1,2))
    ## eta1 + eta2
    rany <- range(lmx1g, lmx2g, eta1.hatL, eta1.hatU, eta2.hatL, eta2.hatU, finite=TRUE)
    ranx <- range(x)
    plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE, 
         xlab="age", ylab="mortality, log-scale",
         main="mortality age-pattern")
    yy <- 10^seq(-7, 2)
    axis(2, at=log(yy), labels=yy)
    axis(1);box()
    ## 1
    for(i in 1:n1){
      segments(x0=low1[i], x1=up1[i], y0=lmx1g[i], y1=lmx1g[i], col=col1, lwd=1)
    }
    xx <- c(x, rev(x))
    yy <- c(eta1.hatL, rev(eta1.hatU))
    polygon(xx, yy, col=col1T, border=col1T)
    lines(x, eta1.hat, col=col1, lwd=2)
    ## 2
    for(i in 1:n1){
      segments(x0=low2[i], x1=up2[i], y0=lmx2g[i], y1=lmx2g[i], col=col2, lwd=1)
    }
    xx <- c(x, rev(x))
    yy <- c(eta2.hatL, rev(eta2.hatU))
    polygon(xx, yy, col=col2T, border=col2T)
    lines(x, eta2.hat, col=col2, lwd=2)
    legend("topleft", inset=0.1,
           legend=c(paste(min(t1), max(t1), sep="-"), "2020"),
           col=c(col1,col2), lwd=3)
    ## deltas
    ranx <- range(-10, x)
    rany <- range(delta.hatL,delta.hatU, c.hatU, c.hatL, na.rm = TRUE)
    plot(1, 1, t="n", xlim=ranx, ylim=rany,
         main=paste("c = ", signif(c.hat,4), "+ -", signif(2*se.c,4)), axes=FALSE,
         xlab="ages", ylab="delta")
    axis(2)
    axis(1);box()
    abline(h=0, col=8, lty=3, lwd=2)
    xx <- c(x, rev(x))
    yy <- c(delta.hatL, rev(delta.hatU))
    polygon(xx, yy, col=coldT, border=coldT)
    lines(x, delta.hat, col=cold, lwd=4)
    abline(v=-1, col=8, lwd=4, lty=2)
    points(-7, c.hat, col=colc, pch=3, lwd=3)
    arrows(x0=-7, x1=-7, y0=c.hatL, y1=c.hatU, col=colc, lwd=3, 
           angle = 90, code=3, length=0.15)
    mtext("c", 1, at=-7, cex=2, line=1.5)
    par(mfrow=c(1,1))
    #locator(1)
  }
}
names(OUT) <- pop
save.image("Output/OutPrepandemic2020.Rdata")
