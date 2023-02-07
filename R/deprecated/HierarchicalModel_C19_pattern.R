## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## R-studio to get the same dir as for the .R file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(magic)
library(MortalitySmooth)
library(colorspace)
library(spam)
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

E <- subset(E0, Country%in%unique(YFM$Country))

## delete "England and Wales" "Moldova" "United Kingdom"
## from deaths since no exposures available
whiNA <- c("England and Wales", "Moldova", "United Kingdom")
whi <- !YFM$Country%in%whiNA
YFM <- YFM[whi,]

#cbind(sort(unique(YFM$Country)),sort(unique(E$Country)))

## delete "Ecuador" since no exposures by sex is available
YFM <- subset(YFM, Country!="Ecuador")
EFM <- subset(E, Country!="Ecuador")

# cbind(sort(unique(YFM$Country)),sort(unique(E$Country)))

## replace NA with zero in the deaths
YFM[which(is.na(YFM[,4])),4] <- 0
YFM[which(is.na(YFM[,5])),5] <- 0

## sorting in the same way: dunno if useful :-)
YFM <- YFM[order(YFM$Country),]
EFM <- EFM[order(EFM$Country),]

## take only 3 pop: 
whi <- c("Argentina", "Australia", "Belgium")
Y <- YFM#subset(YFM, Country%in%whi)
E <- EFM#subset(EFM, Country%in%whi)

## number of age-groups for each pop
n <- table(Y$Country)


## ages for the latent pattern
x <- 0:104
m <- length(x)

## factors
sex <- c("F", "M")
pop <- unique(E$Country)
p <- length(pop)
s <- length(sex)

## building model matrix
X0 <- expand.grid(list(pop=factor(pop),
                       sex=factor(sex)))
X1 <- model.matrix(~-1 + X0$pop + X0$sex)
X11 <- cbind(1, X1)
X2 <- matrix(c(X11), nrow=p*s)

## expand X2 using B-splines
B <- MortSmooth_bbase(x=x, xl=min(x), xr=max(x), ndx=20, deg=3)
B[B<1e-8] <- 0
B <- spam(x=c(B), nrow(B), ncol(B))
nb <- ncol(B)

U <- spam(0, m*nrow(X2), nb*ncol(X2))
whi <- which(X2==1, arr.ind=TRUE)
for(i in 1:nrow(whi)){
  bla <- whi[i,]
  wr <- 1:m+(bla[1]-1)*m
  wc <- 1:nb+(bla[2]-1)*nb
  U[wr,wc] <- B
}

U <- as.matrix(U)
image(X2)## image(U)


## 
CFlist <- CMlist <- Glist <- LMXM <- LMXF <- EFgr <- EMgr <- yFlist <- yMlist <- list()
i=1
for(i in 1:p){
  E.i <- subset(E, Country==pop[i])
  Y.i <- subset(Y, Country==pop[i])
  yF.i <- Y.i$C19f
  yM.i <- Y.i$C19m
  eF.i <- E.i$f
  eM.i <- E.i$m
  age.low.i <- Y.i$Age
  age.up.i <- Y.i$Age+Y.i$AgeInt-1
  age.gr.i <- paste(age.low.i, age.up.i, sep="-")
  CF.i <- matrix(0, n[i], m)
  rownames(CF.i) <- age.gr.i
  colnames(CF.i) <- x 
  CM.i <- CF.i
  G.i <- CF.i
  for(j in 1:n[i]){
    age.low.ij <- age.low.i[j]
    age.up.ij <- age.up.i[j]
    whi <- which(x>=age.low.ij & x<=age.up.ij)
    CF.i[j,whi] <- eF.i[whi]
    CM.i[j,whi] <- eM.i[whi]
    G.i[j,whi] <- 1
  }
  CFlist[[i]] <- CF.i
  CMlist[[i]] <- CM.i
  Glist[[i]] <- G.i
  eFgr.i <- G.i%*%eF.i
  eMgr.i <- G.i%*%eM.i
  EFgr[[i]] <- eFgr.i
  EMgr[[i]] <- eMgr.i
  LMXF[[i]] <- log(yF.i/eFgr.i)
  LMXM[[i]] <- log(yM.i/eMgr.i)
  yFlist[[i]] <- yF.i
  yMlist[[i]] <- yM.i
}

## composite matrix
C <- do.call("adiag",c(CFlist,CMlist))

## starting eta
etaF.st.list <- list()
etaM.st.list <- list()
i=1
for(i in 1:p){
  E.i <- subset(E, Country==pop[i])
  Y.i <- subset(Y, Country==pop[i])
  yF.i <- Y.i$C19f
  yM.i <- Y.i$C19m
  eF.i <- E.i$f
  eM.i <- E.i$m
  eF.i[eF.i==0] <- 0.5
  eM.i[eM.i==0] <- 0.5
  len.i <- Y.i$AgeInt
  yFst.i <- rep(yF.i/len.i, len.i)
  yMst.i <- rep(yM.i/len.i, len.i)
  
  fitF0.i <- Mort1Dsmooth(x=x, y=yFst.i,
                          offset=log(eF.i),
                          method=3, lambda=10^5)
  etaF.st.i <- fitF0.i$logmortality
  etaF.st.list[[i]] <- etaF.st.i
  fitM0.i <- Mort1Dsmooth(x=x, y=yMst.i,
                          offset=log(eM.i),
                          method=3, lambda=10^5)
  etaM.st.i <- fitM0.i$logmortality
  etaM.st.list[[i]] <- etaM.st.i
}
etaF.st <- unlist(etaF.st.list)
etaM.st <- unlist(etaM.st.list)

## penalty stuff
D <- diff(diag(nb), diff=2)
tDD <- t(D)%*%D
lambda <- 10^4
P0 <- lambda*tDD
P <- kronecker(diag(ncol(X2)), P0)
P[1:nb,1:nb] <- 10^8*tDD
## ridge penalty
kappa <- 1e-6
Pr <- kappa*diag.spam(ncol(P))

## regression
eta <- c(etaF.st, etaM.st)
y <- c(unlist(yFlist), unlist(yMlist))
max.it <- 100
for(it in 1:max.it){
  gamma   <- exp(eta)
  mu      <- c(C %*% gamma)
  X       <- (C * ((1 / mu) %*% t(gamma)) ) %*% U
  w       <- as.vector(mu)
  r       <- y - mu + C %*% (gamma * eta)
  G       <- t(X) %*% (w * X) 
  GpP     <- G + P + Pr
  tXr     <- t(X) %*% r
  beta    <- solve(GpP, tXr) 
  eta.old <- eta
  eta     <- U%*%beta
  dif.eta <- max(abs((eta - eta.old)/eta.old) )
  if(dif.eta < 1e-04 & it > 4) break
  cat(it, dif.eta, "\n")
}
save.image("HierarOutGomp.RDdata")
rm(list = ls())
## plotting in a difference device
options(device="X11")
## R-studio to get the same dir as for the .R file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(magic)
library(MortalitySmooth)
library(colorspace)
library(spam)
library(factoextra)
load("HierarOut2.RDdata")


ETAF.hat <- matrix(eta[1:(p*m)], m, p)
ETAM.hat <- matrix(eta[1:(p*m)+p*m], m, p)
beta0 <- beta[1:nb]
eta0 <- B%*%beta0
beta.pop <- beta[1:(p*nb)+nb]
Bpop <- kronecker(diag(p), B)
pi <- Bpop%*%beta.pop
Pi <- matrix(pi, m, p)
beta.sex <- beta[1:nb+(p*nb)+nb]
si <- B%*%beta.sex

for(i in 1:p){
  par(mfrow=c(1,2))
  ## population specific colors
  colocou <- rainbow_hcl(p)
  rany <- c(-15, -3)
  ranx <- range(x)
  Y.i <- subset(Y, Country==pop[i])
  age.low.i <- Y.i$Age
  age.up.i <- Y.i$Age+Y.i$AgeInt
  
  plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
       axes=TRUE,
       xlab="age", ylab="mortality, log-scale",
       main=pop[i])
  lmxF.i <- LMXF[[i]]
  for(j in 1:n[i]){
    segments(x0=age.low.i[j], x1=age.up.i[j],
             y0=lmxF.i[j], y1=lmxF.i[j], col=colocou[i], lwd=3)
  }
  lines(x, ETAF.hat[,i], col=colocou[i], lwd=4)
  legend("topleft", inset = 0.01, legend="Females", bty="n")
  
  plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
       axes=TRUE,
       xlab="age", ylab="mortality, log-scale",
       main=pop[i])
  lmxM.i <- LMXM[[i]]
  for(j in 1:n[i]){
    segments(x0=age.low.i[j], x1=age.up.i[j],
             y0=lmxM.i[j], y1=lmxM.i[j], col=colocou[i], lwd=3)
  }
  lines(x, ETAM.hat[,i], col=colocou[i], lwd=4)
  legend("topleft", inset = 0.01, legend="Males", bty="n")
  locator(1)
}

## country effects
matplot(x, PSI.pop, t="l")

plot(x, psi.sex)

colnames(PSI.pop) <- pop[-1]

## what to cluster
dati <- t(PSI.pop) 

res.dist <- get_dist(dati, stand = TRUE, method = "pearson")
fviz_dist(res.dist, 
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
fviz_nbclust(dati, kmeans, method = "wss")

## decide the number clusters
nc <- 4
set.seed(123)
km.res <- kmeans(dati, nc, nstart = 25)
# Visualize
pdf("ClusterDeviance.pdf")
fviz_cluster(km.res, data = dati,
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal())
dev.off()
res.hc <- hclust(dist(dati,
                      method = "euclidean"),
                 method = "ward.D2")

mycol <- c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07", "darkblue", "darkred")
# plot(res.hc)
# par(mfrow=c(1,2))
fviz_dend(res.hc, k = nc, # Cut in nk groups
          cex = 1, # label size
          k_colors = mycol[1:nc],
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE, # Add rectangle around groups
          horiz=TRUE
)




## !!!!!! select sex
sex <- "F"
## !!!!!! select sex

if(sex=="F"){
  whi <- c(1,2,3,4)
}
if(sex=="M"){
whi <- c(1,2,3,5)
}
Y <- YFM[,whi]
E <- E[,whi]

## replace name
names(Y)[4] <- "Deaths"
names(E)[4] <- "Exposures"



## replace NA with zero in the deaths
Y$Deaths[which(is.na(Y$Deaths))] <- 0


## divide data by regions
EUR <- c("Belgium", "Czechia", "France", "Germany", "Greece", "Hungary",
         "Israel", "Italy", "Netherlands", "Norway",
         "Portugal","Romania", "Scotland", "Slovenia",
         "Spain", "Switzerland", "Turkey","Ukraine")

AFR <- c("Chad", "Eswatini", "Kenya", "Malawi","Nigeria")

ASI <- c("India", "Iraq", "Japan", "Jordan", "Nepal","Pakistan","Philippines")
LAM <- c("Argentina", "Brazil", "Chile", "Colombia", "Cuba", "Mexico",
         "Panama", "Paraguay", "Peru",  "Uruguay")

NAM <- c("Canada", "USA")

OCE <- c("Australia")
length(EUR)+length(AFR)+length(ASI)+length(LAM)+length(NAM)+length(OCE)
length(unique(Y$Country))

## ages for the latent pattern
x <- 0:104
m <- length(x)
PLOT=TRUE


ALL <- c(EUR, AFR, ASI, LAM, NAM, OCE)
## number of pop
nc <- length(ALL)
## number of age-group within each pop
n <- table(Y$Country)
n <- n[match(ALL,names(n))]

Clist <- Glist <- LMX <- Egr <- list()
ETA <- ETA.LOW <- ETA.UP <- matrix(NA, m, nc)
colnames(ETA) <- colnames(ETA.LOW) <- colnames(ETA.UP) <- ALL
i=1
## penalty stuff
D <- diff(diag(m), diff=2)
tDD <- t(D)%*%D
lambda <- 10^4
P <- lambda*tDD

for(i in 1:nc){

  ## building country-specific composite matrix
  ## and place them in a list
  E.i <- subset(E, Country==ALL[i])
  Y.i <- subset(Y, Country==ALL[i])
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
    colocou <- rainbow_hcl(nc)
    rany <- c(-15, -3)
    ranx <- range(x)
    plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
         axes=TRUE,
         xlab="age", ylab="mortality, log-scale",
         main=ALL[i])
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