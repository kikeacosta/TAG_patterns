## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## R-studio to get the same dir as for the .R file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(magic)
library(colorspace)
library(MortalitySmooth)

Y <- readRDS("C19_use.rds", refhook = NULL)
## delete last row for Argentina
Y <- Y[-115,]
E0 <- readRDS("offsets.rds", refhook = NULL)
## select E where we have data in D
E <- subset(E0, Country%in%unique(Y$Country))

names(table(Y$Country))
names(table(E$Country))
## divide data by regions
EUR <- c("Austria", "Belgium", "Czechia", "Denmark", 
         "Finland","France", "Germany", "Greece", "Hungary",
         "Ireland", "Israel", "Italy", "Netherlands", "Northern Ireland", "Norway",
         "Portugal","Romania", "Scotland", "Slovenia",
         "Spain", "Sweden","Switzerland", "Turkey","Ukraine")
AFR <- c("Chad", "Eswatini", "Kenya", "Malawi","Nigeria","Somalia", "South Africa")
ASI <- c("Afghanistan", "Bangladesh", "China", "India", "Indonesia", "Iraq", "Jordan",
         "Japan", "Nepal","Pakistan","Philippines", "South Korea")
LAM <- c("Argentina", "Brazil", "Chile", "Colombia", "Costa Rica", "Cuba", "Ecuador",
         "Mexico",
         "Panama", "Paraguay", "Peru", "Suriname", "Uruguay")
NAM <- c("Canada", "Jamaica", "USA")
OCE <- c("Australia")
length(EUR)+length(AFR)+length(ASI)+length(SAM)+length(NAM)+length(OCE)
length(unique(Y$Country))

## only for a given region
PLOT=FALSE

REGIONS <- list(EUR=EUR, AFR=AFR, ASI=ASI, LAM=LAM, NAM=NAM, OCE=OCE)
k=1
for(k in 1:length(REGIONS)){
  ## only for a given region
  REG <- REGIONS[[k]]
  Y.REG <- subset(Y, Country%in%REG)
  E.REG <- subset(E, Country%in%REG)
  ## number of age-group within each pop
  n <- table(Y.REG$Country)
  ## number of pop
  nc <- length(n)
  ## ages
  x <- 0:104
  m <- length(x)
  
  ## building country-specific composite matrix
  ## and place them in a list
  Clist <- Glist <- LMX <- Egr <- list()
  for(i in 1:nc){
    E.i <- subset(E.REG, Country==REG[i])
    Y.i <- subset(Y.REG, Country==REG[i])
    e.i <- E.i$Population
    age.low.i <- Y.i$Age
    age.up.i <- Y.i$Age+Y.i$AgeInt-1
    age.gr.i <- paste(age.low.i, age.up.i, sep="-")
    C.i <- matrix(0, n[i], m)
    rownames(C.i) <- age.gr.i
    colnames(C.i) <- x 
    G.i <- C.i
    for(j in 1:n[i]){
      age.low.ij <- age.low.i[j]
      age.up.ij <- age.up.i[j]
      whi <- which(x>=age.low.ij & x<=age.up.ij)
      C.i[j,whi] <- e.i[whi]
      G.i[j,whi] <- 1
    }
    Clist[[i]] <- C.i
    Glist[[i]] <- G.i
    egr.i <- G.i%*%e.i
    Egr[[i]] <- egr.i
    LMX[[i]] <- log(Y.i$Deaths/egr.i)
  }
  ## composite matrix for the overall model
  C <- do.call("adiag",Clist)
  
  ## response
  y <- Y.REG$Deaths
  
  ## penalty stuff
  D <- diff(diag(m), diff=2)
  tDD <- t(D)%*%D
  lambda <- 10^4
  if(k==3){lambda<-10^5}
  if(k==5){lambda<-10^5}
  P0 <- lambda*tDD
  P <- adiag(P0,0*diag(nc-1))
  
  ## model matrix
  U0 <- kronecker(rep(1,nc), diag(m))
  U1 <- kronecker(diag(nc), rep(1,m))
  U1 <- U1[,-1]
  U <- cbind(U0, U1)
  
  ## starting eta
  eta.st.list <- list()
  for(i in 1:nc){
    E.i <- subset(E.REG, Country==REG[i])
    Y.i <- subset(Y.REG, Country==REG[i])
    d.i <- Y.i$Deaths
    e.i <- E.i$Population
    e.i[e.i==0] <- 0.5
    len.i <- Y.i$AgeInt
    dst.i <- rep(d.i/len.i, len.i)
    fit0.i <- Mort1Dsmooth(x=x, y=dst.i,
                           offset=log(e.i),
                           method=3, lambda=10^5)
    eta.st.i <- fit0.i$logmortality
    eta.st.list[[i]] <- eta.st.i
  }
  eta.st <- unlist(eta.st.list)
  
  ## regression
  eta <- eta.st
  max.it <- 100
  for(it in 1:max.it){
    gamma   <- exp(eta)
    mu      <- c(C %*% gamma)
    X       <- (C * ((1 / mu) %*% t(gamma)) ) %*% U
    w       <- as.vector(mu)
    r       <- y - mu + C %*% (gamma * eta)
    G       <- t(X) %*% (w * X) 
    GpP     <- G + P
    tXr     <- t(X) %*% r
    beta    <- solve(GpP, tXr) 
    eta.old <- eta
    eta     <- U%*%beta
    dif.eta <- max(abs((eta - eta.old)/eta.old) )
    if(dif.eta < 1e-04 & it > 4) break
    #cat(it, dif.eta, "\n")
  }
  
  ETA <- ETA.LOW <- ETA.UP <- matrix(NA, m, nc)
  ## fitted log-mortality
  eta.hat <- list()
  for(i in 1:nc){
    eta.hat[[i]] <- eta[1:m+(i-1)*m]
    ETA[,i] <- eta[1:m+(i-1)*m]
  }
  
  theta.hat <- beta[-c(1:m)]
  
  ## standard errors
  H0 <- solve(GpP)
  H1 <- H0 %*% G %*% H0
  se.theta <- sqrt(diag(H1)[-c(1:m)])
  se <- sqrt(diag(U %*% H1 %*% t(U)))
  
  ## upper and lower bound for eta.hat and theta.hat
  eta.hatL <- eta.hatU <- list()
  se.hat <- list()
  for(i in 1:nc){
    se.hat[[i]] <- se[1:m+(i-1)*m]
    eta.hatL[[i]] <- eta.hat[[i]] - 2*se.hat[[i]]
    eta.hatU[[i]] <- eta.hat[[i]] + 2*se.hat[[i]]
    ETA.LOW[,i] <- eta.hatL[[i]]
    ETA.UP[,i] <- eta.hatU[[i]]
  }
  theta.hatL <- theta.hat - 2*se.theta
  theta.hatU <- theta.hat + 2*se.theta
  
  
  if(PLOT){
    colocou <- rainbow_hcl(nc)
    par(mfrow=c(1,2))
    ## plotting outcomes
    rany <- c(-15, -3)
    ranx <- range(x)
    plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
         axes=FALSE,
         xlab="age", ylab="mortality, log-scale",
         main="OCEANIA C19 mortality age-pattern")
    yy <- 10^seq(-7, -2)
    axis(2, at=log(yy), labels=yy)
    axis(1);box()
    xx <- c(x, rev(x))
    for(i in 1:nc){
      Y.i <- subset(Y.REG, Country==REG[i])
      age.low.i <- Y.i$Age
      age.up.i <- Y.i$Age+Y.i$AgeInt
      lmx.i <- LMX[[i]]
      for(j in 1:n[i]){
        segments(x0=age.low.i[j], x1=age.up.i[j],
                 y0=lmx.i[j], y1=lmx.i[j], col=colocou[i], lwd=3)
      }
      yy <- c(eta.hatL[[i]], rev(eta.hatU[[i]]))
      polygon(xx, yy, col=adjustcolor(colocou[i], 0.5), border=adjustcolor(colocou[i], 0.5))
      lines(x, eta.hat[[i]], col=colocou[i], lwd=4)
    }
    legend("topleft", inset=0.01,
           legend=REG, col=colocou, pch=15, lwd=2)
    
    plot(c(0,theta.hat), 1:nc, pch=15, cex=1.2, col=colocou,
         ylim=c(1,nc+0.5), xlim=range(0, theta.hatL,theta.hatU), axes=FALSE,
         xlab="", ylab="")
    axis(1);box()
    abline(v=0, col=8, lwd=2, lty=3)
    grid(ny=NA)
    for(i in 2:nc){
      arrows(x0=theta.hatL[i-1], y0=i, x1 = theta.hatU[i-1], col=colocou[i], 
             lwd=2, code=3, length=0.1, angle=90)
    }
    for(i in 1:nc){
      text(c(0,theta.hat)[i], y=i+0.5, REG[i], col=colocou[i])
    }
    par(mfrow=c(1,1))
  }
  
  ## saving outcomes
  nome <- paste("OUTCOMES/", names(REGIONS)[k], "region.RData", sep="")
  save.image(nome)
}
