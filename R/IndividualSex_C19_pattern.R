## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## R-studio to get the same dir as for the .R file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(MortalitySmooth)
library(colorspace)

## reading deaths
YFM0 <- readRDS("C19_use_sex.rds", refhook = NULL)
## delete last row for Argentina
YFM0 <- YFM0[-106,]

## deaths
YFM <- YFM0[,c(1,3,4,5,6)]

## exposures
E0 <- readRDS("offsets.rds", refhook = NULL)
unique(YFM$Country)
unique(E0$Country)
## delete if no exposures available by sex

EFM <- subset(E0, Country%in%unique(YFM$Country))

## delete "England and Wales" "Moldova" "United Kingdom"
## from deaths since no exposures available
whiNA <- c("England and Wales", "Moldova", "United Kingdom")
whi <- !YFM$Country%in%whiNA
YFM <- YFM[whi,]

#cbind(sort(unique(YFM$Country)),sort(unique(E$Country)))

## delete "Ecuador" since no exposures by sex is available
YFM <- subset(YFM, Country!="Ecuador")
EFM <- subset(EFM, Country!="Ecuador")

# cbind(sort(unique(YFM$Country)),sort(unique(E$Country)))

## !!!!!! select sex
sex <- "M"
## !!!!!! select sex

if(sex=="F"){
  whi <- c(1,2,3,4)
  Y <- YFM[,whi]
  E <- EFM[,whi]
}
if(sex=="M"){
  whi <- c(1,2,3,5)
  Y <- YFM[,whi]
  E <- EFM[,whi]
}
if(sex=="B"){
  Y <- data.frame(YFM[,1:3], YFM[,4]+YFM[,5])
  E <- data.frame(EFM[,1:3], EFM[,4]+EFM[,5])
}


## replace name
names(Y)[4] <- "Deaths"
names(E)[4] <- "Exposures"



## replace NA with zero in the deaths
Y$Deaths[which(is.na(Y$Deaths))] <- 0


Y <- Y[order(Y$Country),]
E <- E[order(E$Country),]


## ages for the latent pattern
x <- 0:104
m <- length(x)
PLOT=FALSE

## populations
pop <- unique(E$Country)
## number of pop
p <- length(pop)
## number of age-group within each pop
n <- table(Y$Country)

Clist <- Glist <- LMX <- Egr <- list()
ETA <- ETA.LOW <- ETA.UP <- matrix(NA, m, p)
colnames(ETA) <- colnames(ETA.LOW) <- colnames(ETA.UP) <- pop
i=1
## penalty stuff
D <- diff(diag(m), diff=2)
tDD <- t(D)%*%D
lambda <- 10^4
P <- lambda*tDD

for(i in 1:p){

  ## building country-specific composite matrix
  ## and place them in a list
  E.i <- subset(E, Country==pop[i])
  Y.i <- subset(Y, Country==pop[i])
  y.i <- Y.i$Deaths
  e.i <- E.i$Exposures
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
  
  ## PCLM iterations
  
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
    colocou <- rainbow_hcl(p)
    rany <- c(-15, -3)
    ranx <- range(x)
    plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
         axes=TRUE,
         xlab="age", ylab="mortality, log-scale",
         main=pop[i])
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
nome <- paste("OUTCOMES/ALLIndiv", sex, ".RData", sep="")
save.image(nome)
##


















## END