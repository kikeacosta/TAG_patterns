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
col2F <- "blue"
col2FT <- adjustcolor(col2F, 0.5)

col1M <- "brown1"
col1MT <- adjustcolor(col1M, 0.5)
col2M <- "deepskyblue2"
col2MT <- adjustcolor(col2M, 0.5)

cold <- "orange"
coldT <- adjustcolor(cold, 0.5)
colc <- "violet"
colcT <- adjustcolor(colc, 0.1)

cols <- "coral"
colsT <- adjustcolor(cols, 0.5)


colcou <- rainbow_hcl(p)
colcouT <- adjustcolor(colcou, 0.3)


OUT <- list()
PLOT <- FALSE
j=48
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
  ## estimation
  ## response
  y <- c(y1F, y1M, y2F, y2M)
  ## overall composite matrix
  C <- adiag(C1F, C1M, C2F, C2M)
  
  ## finding starting values
  y1F.st0 <- rep(y1F/len1, len1)
  wei <- rep(1,m)
  wei[e1F==0] <- 0
  fit1F.0 <- Mort1Dsmooth(x=x, y=y1F.st0,
                          offset=log(e1F),
                          w=wei,
                          method=3, lambda=10^4)
  y1M.st0 <- rep(y1M/len1, len1)
  wei <- rep(1,m)
  wei[e1M==0] <- 0
  fit1M.0 <- Mort1Dsmooth(x=x, y=y1M.st0,
                          offset=log(e1M),
                          w=wei,
                          method=3, lambda=10^4)
  y2F.st0 <- rep(y2F/len2, len2)
  wei <- rep(1,m)
  wei[e2F==0] <- 0
  fit2F.0 <- Mort1Dsmooth(x=x, y=y2F.st0,
                          offset=log(e2F),
                          w=wei,
                          method=3, lambda=10^4)
  y2M.st0 <- rep(y2M/len2, len2)
  wei <- rep(1,m)
  wei[e2M==0] <- 0
  fit2M.0 <- Mort1Dsmooth(x=x, y=y2M.st0,
                          offset=log(e2M),
                          w=wei,
                          method=3, lambda=10^4)
  
  
  eta.st <- c(fit1F.0$logmortality, fit1M.0$logmortality,
              fit2F.0$logmortality, fit2M.0$logmortality)
  
  ## penalty stuff
  Deta1F <- diff(diag(m), diff=2)
  tDDeta1F <- t(Deta1F)%*%Deta1F
  Ds <- diff(diag(m), diff=2)
  tDDs <- t(Ds)%*%Ds
  Ddelta <- diff(diag(m), diff=2)
  tDDdelta <- t(Ddelta)%*%Ddelta
  
  lambda.eta1F <- 10^3
  lambda.s <- 10^3.5
  lambda.delta <- 10^5
  
  P <- adiag(lambda.eta1F*tDDeta1F,
             lambda.s*tDDs,
             0,
             lambda.delta*tDDdelta)
  ## ridge penalty
  Pr <- 1e-4*diag(ncol(P))
  
  ## model matrix
  Ueta1F <- kronecker(rep(1,4), diag(m))
  Us <- rbind(0*diag(m), diag(m), 0*diag(m), diag(m))
  Uc <- matrix(c(rep(0,2*m), rep(1,2*m)), 4*m)
  Udelta <- rbind(0*diag(m), 0*diag(m), diag(m), diag(m))
  U <- cbind(Ueta1F,Us,Uc,Udelta)
  
  ## constraining delta to sum up to 0
  H <- matrix(c(rep(0,2*m+1), rep(1,m)), nrow=1)
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
    betas   <- coeff[1:(m*3+1)]
    eta.old <- eta
    eta     <- U%*%betas
    dif.eta <- max(abs((eta - eta.old)/eta.old))
    if(dif.eta < 1e-04 & it > 4) break
    #cat(it, dif.eta, "\n")
  }
  
  eta1F.hat <- betas[1:m]
  s.hat <- betas[1:m+m]
  c.hat <- betas[2*m+1]
  delta.hat <- betas[1:m+2*m+1]
  eta1M.hat <- eta1F.hat + s.hat
  eta2F.hat <- eta1F.hat + c.hat + delta.hat
  eta2M.hat <- eta1F.hat + s.hat + c.hat + delta.hat # eta1M.hat + c.hat + delta.hat
  
  ## standard errors
  H0 <- solve(GpP)
  Vbetas <- H0 %*% G %*% H0
  se.betas <- sqrt(diag(Vbetas))
  se.eta1F <- se.betas[1:m]
  se.s <- se.betas[1:m+m]
  se.c <- se.betas[2*m+1]
  se.delta <- se.betas[1:m+1+2*m]
  
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
  c.hatL <- c.hat - 2*se.c
  c.hatU <- c.hat + 2*se.c
  delta.hatL <- delta.hat - 2*se.delta
  delta.hatU <- delta.hat + 2*se.delta
  
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
  
  OUT.j$s.hat <- s.hat
  OUT.j$s.hatL <- s.hatL
  OUT.j$s.hatU <- s.hatU
  
  OUT.j$delta.hat <- delta.hat
  OUT.j$delta.hatL <- delta.hatL
  OUT.j$delta.hatU <- delta.hatU
  OUT.j$c.hat <- c.hat
  OUT.j$c.hatL <- c.hatL
  OUT.j$c.hatU <- c.hatU
  
  OUT[[j]] <- OUT.j
  
  if(PLOT){
    par(mfrow=c(2,2))
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
    ## sex-smooth-factor
    ranx <- range(x)
    rany <- range(s.hatL,s.hatU)
    plot(1, 1, t="n", xlim=ranx, ylim=rany, axes=FALSE,
         xlab="ages", ylab="s", main="Sex-factor")
    axis(2)
    axis(1);box()
    abline(h=0, col=8, lty=3, lwd=2)
    xx <- c(x, rev(x))
    yy <- c(s.hatL, rev(s.hatU))
    polygon(xx, yy, col=colsT, border=colsT)
    lines(x, s.hat, col=cols, lwd=4)
    
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
    locator(1)
  }
}
names(OUT) <- pop

save.image("Output/OutPrepandemic2020sex_CommonC19factor.Rdata")
