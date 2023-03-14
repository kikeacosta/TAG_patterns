## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## R-studio to get the same dir as for the .R file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(MortalitySmooth)
library(colorspace)

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
length(EUR)+length(AFR)+length(ASI)+length(LAM)+length(NAM)+length(OCE)
length(unique(Y$Country))

## ages for the latent pattern
x <- 0:104
m <- length(x)
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
  ## penalty stuff
  D <- diff(diag(m), diff=2)
  tDD <- t(D)%*%D
  lambda <- 10^4
  P <- lambda*tDD

  ETA <- ETA.LOW <- ETA.UP <- matrix(NA, m, nc)
  ## PCLM iterations
  for(i in 1:nc){
    C.i <- Clist[[i]]
    E.i <- subset(E.REG, Country==REG[i])
    Y.i <- subset(Y.REG, Country==REG[i])
    y.i <- Y.i$Deaths
    e.i <- E.i$Population
    e.i[e.i==0] <- 0.5
    len.i <- Y.i$AgeInt
    dst.i <- rep(y.i/len.i, len.i)
    fit0.i <- Mort1Dsmooth(x=x, y=dst.i,
                           offset=log(e.i),
                           method=3, lambda=10^4)
    eta.st.i <- fit0.i$logmortality
    eta <- eta.st.i
    max.it <- 100
    for(it in 1:max.it){
      gamma   <- exp(eta)
      mu      <- c(C.i %*% gamma)
      X       <- C.i * ((1 / mu) %*% t(gamma))
      w       <- as.vector(mu)
      r       <- y.i - mu + C.i %*% (gamma * eta)
      G       <- t(X) %*% (w * X) 
      GpP     <- G + P
      tXr     <- t(X) %*% r
      eta.old <- eta
      eta    <- solve(GpP, tXr) 
      dif.eta <- max(abs((eta - eta.old)/eta.old) )
      if(dif.eta < 1e-04 & it > 4) break
      #cat(it, dif.eta, "\n")
    }
    eta.hat <- eta
    H0 <- solve(GpP)
    H1 <- H0 %*% G %*% H0
    se <- sqrt(diag(H1))
    eta.hatL <- eta.hat - 2*se
    eta.hatU <- eta.hat + 2*se
    ETA[,i] <- eta.hat
    ETA.LOW[,i] <- eta.hatL
    ETA.UP[,i] <- eta.hatU
    
    if(PLOT){
      ## population specific colors
      colocou <- rainbow_hcl(nc)
      rany <- c(-15, -3)
      ranx <- range(x)
      plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
           axes=FALSE,
           xlab="age", ylab="mortality, log-scale",
           main=REG[i])
      Y.i <- subset(Y.REG, Country==REG[i])
      age.low.i <- Y.i$Age
      age.up.i <- Y.i$Age+Y.i$AgeInt
      lmx.i <- LMX[[i]]
      for(j in 1:n[i]){
        segments(x0=age.low.i[j], x1=age.up.i[j],
                 y0=lmx.i[j], y1=lmx.i[j], col=colocou[i], lwd=3)
      }
      xx <- c(x, rev(x))
      yy <- c(eta.hatL, rev(eta.hatU))
      polygon(xx, yy, col=adjustcolor(colocou[i], 0.5), border=adjustcolor(colocou[i], 0.5))
      lines(x, eta.hat, col=colocou[i], lwd=4)
      locator(1)
    }
  }

  ## saving outcomes
  nome <- paste("OUTCOMES/", names(REGIONS)[k], "indiv.RData", sep="")
  save.image(nome)
}

##


















## END