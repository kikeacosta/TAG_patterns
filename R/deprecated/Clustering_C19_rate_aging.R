## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## R-studio to get the same dir as for the .R file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(magic)
library(colorspace)
library(MortalitySmooth)
library(colorspace)
library(gridExtra)

source("EmpiricDeriv.R")


Y <- readRDS("C19_use.rds", refhook = NULL)
## delete last row for Argentina
Y <- Y[-115,]
## delete England and England and Wales
Y <- subset(Y, Country!="England")
Y <- subset(Y, Country!="England and Wales")
E0 <- readRDS("Offsets.rds", refhook = NULL)
## only both sex
E1 <- subset(E0, Sex=="b")
# delete specific Regions
E1 <- subset(E1, Region=="All")
## select E where we have data in D
E <- subset(E1, Country%in%unique(Y$Country))


cbind(names(table(Y$Country)), names(table(E$Country)))

## sorting by name
Y <- Y[order(Y$Country),]
E <- E[order(E$Country),]


## define the age range on which we aim to do our analysis
x1 <- 30
x2 <- 70
xs <- x1:x2


## whether to plot outcomes 
PLOT <- FALSE

## populations
pop <- unique(E$Country)
## number of pop
p <- length(pop)
## country colors
colocou <- rainbow_hcl(p)

Clist <- Glist <- LMX <- Egr <- list()
ETA <- ETA.LOW <- ETA.UP <- ETA1 <- list()
ETAs <- ETA1s <- matrix(0, length(x1:x2), p)
colnames(ETAs) <- colnames(ETA1s) <- pop
i=37
for(i in 1:p){
  
  ## building country-specific composite matrix
  ## and place them in a list
  
  ## given selected age-range
  
  ## select age-grouped deaths
  Y.i0 <- subset(Y, Country==pop[i])
  pmin <- max(which(Y.i0$Age<=x1))
  if(max(Y.i0$Age)>=x2){
    pmax <- max(which(Y.i0$Age<=x2))
  }else{
    pmax <- which.max(Y.i0$Age)
  }
  
  Y.i <- Y.i0[pmin:pmax, ]
  y.i <- Y.i$Deaths
  ## number of age-group for the ith pop
  n <- length(y.i)
  
  age.low.i <- Y.i$Age
  age.up.i <- Y.i$Age+Y.i$AgeInt-1
  age.min <- min(age.low.i)
  age.max <- max(age.up.i)
  ## select exposures
  E.i <- subset(E, Country==pop[i])
  e.i0 <- E.i$Population
  whi.ar <- which(E.i$Age>=age.min & E.i$Age<=age.max)
  e.i <- e.i0[whi.ar]
  ## ages for the latent pattern
  x <- E.i$Age[whi.ar]
  m <- length(x)
  
  
  C.i <- matrix(0, n, m)
  age.gr.i <- paste(age.low.i, age.up.i, sep="-")
  rownames(C.i) <- age.gr.i
  colnames(C.i) <- x 
  G.i <- C.i
  for(j in 1:n){
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
  
  ## PCLM procedure
  
  ## starting values
  e.i[e.i==0] <- 0.5
  len.i <- Y.i$AgeInt
  dst.i <- rep(y.i/len.i, len.i)
  fit0.i <- Mort1Dsmooth(x=x, y=dst.i,
                         offset=log(e.i),
                         method=3, lambda=10^4)
  eta.st.i <- fit0.i$logmortality
  eta <- eta.st.i
  max.it <- 100
  
  ## penalty stuff
  D <- diff(diag(m), diff=2)
  tDD <- t(D)%*%D
  
  ## optimizing lambda??
  lambdas <- 10^seq(4, 7, 1)
  if(pop[i]=="India"){
    lambdas <- 10^seq(5, 7, 1)
  }
  nl <- length(lambdas)
  BICs <- numeric(nl)
  k=1
  for(k in 1:nl){
    P <- lambdas[k]*tDD
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
    }
    ## diagnostic
    ok <- y.i > 0 & mu > 0
    Dev <- 2 * sum(y.i[ok] * log(y.i[ok] / mu[ok]) )
    H <- solve(GpP, G)
    Ed <- sum(diag(H))
    BICs[k] <- Dev + log(n)*Ed
  }
  lambda.hat <- lambdas[which.min(BICs)]
  P <- lambda.hat*tDD
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
  ETA[[i]] <- eta.hat
  ETA.LOW[[i]] <- eta.hatL
  ETA.UP[[i]] <- eta.hatU
  ## compute rate-of-aging
  eta1 <- emp.der(x=x, y=eta)
  ETA1[[i]] <- eta1
  ## take only the selected ages
  sel <- which(c(age.min:age.max)%in%xs)
  ETAs[,i] <- eta.hat[sel]
  ETA1s[,i] <- eta1[sel]
  cat(pop[i], lambda.hat, "\n")
  
  if(PLOT){
    ## population specific colors
    rany <- c(-15, -3)
    ranx <- range(x)
    plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
         axes=TRUE,
         xlab="age", ylab="mortality, log-scale",
         main=pop[i])
    age.low.i <- Y.i$Age
    age.up.i <- Y.i$Age+Y.i$AgeInt
    lmx.i <- LMX[[i]]
    for(j in 1:n){
      segments(x0=age.low.i[j], x1=age.up.i[j],
               y0=lmx.i[j], y1=lmx.i[j], col=colocou[i], lwd=3)
    }
    xx <- c(x, rev(x))
    yy <- c(eta.hatL, rev(eta.hatU))
    polygon(xx, yy, col=adjustcolor(colocou[i], 0.5), border=adjustcolor(colocou[i], 0.5))
    lines(x, eta.hat, col=colocou[i], lwd=4)
    abline(v=c(x1,x2), col=2, lwd=2, lty=2)
    #locator(1)
    Sys.sleep(0.5)
  }
}
names(ETA) <- names(ETA.LOW) <- names(ETA.UP) <- names(ETA1) <- pop

par(mfrow=c(2,1))
matplot(xs, ETAs, col=colocou, t="l", lty=1)
matplot(xs, ETA1s, col=colocou, t="l", lty=1)
par(mfrow=c(1,1))

## clustering
library(factoextra)
library(cluster)
library(ppclust)
## data to cluster
tETA1s <- t(ETA1s)
## number of clusters
nc <- 5
coloclu <- rainbow_hcl(nc)
## Run FCM 
fcm.res <- fcm(tETA1s, centers=nc, numseed=2021, m=3)
## fuzzy memberships degrees of the data objects
par(las=1) # make label text perpendicular to axis
par(mar=c(5,8,4,2)) # increase y-axis margin.
barplot(t(fcm.res$u), horiz=TRUE, col=coloclu,
        names.arg=pop, cex.names=0.8)
## plot cluster prototypes
matplot(xs, t(fcm.res$v), t="l", lwd=5, lty=1)
## add pop-specific data
matplot(xs, ETA1s, col=coloclu[fcm.res$cluster], t="l", lty=3)
matlines(xs, t(fcm.res$v), col=coloclu, t="l", lwd=8, lty=1)
## clusplot
fcm.res2 <- ppclust2(fcm.res, "kmeans")
fviz_cluster(fcm.res2, data = tETA1s, 
             ellipse.type = "convex",
             palette = coloclu,
             repel = TRUE)

fcm.res3 <- ppclust2(fcm.res, "fclust")
# Fuzzy Silhouette Index:
idxsf <- SIL.F(fcm.res3$Xca, fcm.res3$U, alpha=1)
paste("Fuzzy Silhouette Index: ",idxsf)

ETA1s.hat <- ETA1s*0
for(i in 1:p){
  ETA1s.hat[,i] <- fcm.res$u[i,1]*fcm.res$v[1,]
  j=2
  for(j in 2:nc){
    ETA1s.hat[,i] <- ETA1s.hat[,i] + fcm.res$u[i,j]*fcm.res$v[j,]
  }
}


write.table(fcm.res$u, "./OutCluster/Membership5cluster3m.txt")
write.table(fcm.res$v, "./OutCluster/Centers5cluster3m.txt")
write.table(ETA1s.hat, "./OutCluster/FittedRateAging5cluster3m.txt")




write.table(ETAs, "./OutCluster/LogMortality.txt")
write.table(ETA1s, "./OutCluster/RateAgingClusteredData.txt")












SVD <- svd(ETA1s)

U <- SVD$u
V <- SVD$v
dim(V)

matplot(V[,1:2])

matplot(xs, U[,1:2])
SVD$d

BLA <- matrix(0, p, 3)

for(i in 1:p){
  BLA[i,] <- V[i,1:3]*SVD$d[1:3]
}

i=
bla <- SVD$d[1]*U[,1]*V[i,1] + SVD$d[2]*U[,2]*V[i,2] + SVD$d[3]*U[,3]*V[i,3] 
matplot(xs, ETA1s, t="n")
lines(xs, tETA1s[i,], col=colocou[i], lwd=5, t="l")
lines(xs, bla, col=2, lwd=5)

lines(xs, bla, col=3, lwd=5)



require(fclust)
## data to cluster
tETA1s <- (ETA1s)

CM <- cor(tETA1s)
ED <- eigen(CM)

U <- ED$vectors
V <- ED$values

## number of clusters
nc <- 3
coloclu <- rainbow_hcl(nc)
## Run FCM 
fcm.res <- fcm(tETA1s, centers=nc, numseed=2021, m=1)
## fuzzy memberships degrees of the data objects
par(las=1) # make label text perpendicular to axis
par(mar=c(5,8,4,2)) # increase y-axis margin.
barplot(t(fcm.res$u), horiz=TRUE, col=coloclu,
        names.arg=pop, cex.names=0.8)
## plot cluster prototypes
matplot(xs, t(fcm.res$v), t="l", lwd=5, lty=1)
## add pop-specific data
matplot(xs, ETA1s, col=coloclu[fcm.res$cluster], t="l", lty=3)
matlines(xs, t(fcm.res$v), col=coloclu, t="l", lwd=8, lty=1)
## clusplot
fcm.res2 <- ppclust2(fcm.res, "kmeans")
fviz_cluster(fcm.res2, data = tETA1s, 
             ellipse.type = "convex",
             palette = coloclu,
             repel = TRUE)

fcm.res3 <- ppclust2(fcm.res, "fclust")
# Fuzzy Silhouette Index:
idxsf <- SIL.F(fcm.res3$Xca, fcm.res3$U, alpha=1)
paste("Fuzzy Silhouette Index: ",idxsf)

## trials 


matplot(xs, ETA1s, t="n")
fcm.res$clus
i=1
lines(xs, tETA1s[i,], col=colocou[i], lwd=5, t="l")
lines(xs, fcm.res$u[i,1]*fcm.res$v[1,]+fcm.res$u[i,2]*fcm.res$v[2,]+fcm.res$u[i,3]*fcm.res$v[3,],
      col=colocou[i], lwd=5, lty=2)
i=13
lines(xs, tETA1s[i,], col=colocou[i], lwd=5, t="l")
lines(xs, fcm.res$u[i,1]*fcm.res$v[1,]+fcm.res$u[i,2]*fcm.res$v[2,]+fcm.res$u[i,3]*fcm.res$v[3,],
      col=colocou[i], lwd=5, lty=2)

res.fcm <- fcm(tETA1s, centers=3, memberships=0.5)


## final  matrices
res.fcm$v










library(factoextra)
library(cluster)
library(e1071)


## trial
tETA1s <- t(ETA1s)
FC <- fanny(tETA1s, 3, memb.exp = 1.1)
p1 <- fviz_cluster(FC, data = dati,
             ellipse.type = "convex",
             palette = "jco")
fviz_silhouette(FC, palette = "jco",
                ggtheme = theme_minimal())

plot(FC)
## reconstructing Argentina
str(FC)

library(fclust)


FC <- FKM(X = tETA1s, k=3, m = 1.2, RS = 50, stand = 1, index = "SIL.F")
summary(FC)
matlines(xs, t(Hraw(X = tETA1s, H = FC$H)), t="l", lwd=2)

library(fclust)
library(ppclust)
FC1 <- fcm(tETA1s, centers=3)
## plotcluster(FC1, cp=1, trans=TRUE)

res.fcm2 <- ppclust2(FC1, "kmeans")

p2 <- fviz_cluster(res.fcm2, data = tETA1s, 
             ellipse.type = "convex",
             palette = "jco",
             repel = TRUE)
grid.arrange(p1, p2, ncol=2)


fviz_cluster(FC, data = dati,
             ellipse.type = "convex",
             palette = "jco")

plot(FC, pca=TRUE)
## possible number of clusters
NC <- 3:5






## number of cluster
nc <- 3
## age range to consider
ar <- which(x%in%30:70)
## what to cluster
dati <- t(ETA1[ar,]) 
## hierarchical clustering
HC <- hkmeans(dati, nc)


HCout[,i] <- 



## membership exponent for the fuzzy clustering
me <- 1.2
## fuzzy clustering
FC <- fanny(dati, nc, memb.exp = 1.3)
head(round(FC$membership*100,1), 10)
p2 <- fviz_cluster(FC, data = dati,
                   ellipse.type = "convex",
                   palette = "jco",
                   ggtheme = theme_minimal(),
                   main="Fuzzy")
p2



write.table("~OutCluster/bla.txt", )




matplot(x, ETA, t="l", col=colocou, lty=1)
lines(x, ETA[,which(pop=="Denmark")], col=2, lwd=4)



matplot(x, t(dati), t="l", col=colocou, lty=1)
lines(x, dati[which(pop=="Denmark"),], col=2, lwd=4)
abline(h=0, col="grey", lty=2)


## decide the number clusters
nc <- 3
mycol <- c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07", "darkblue", "darkred")[1:nc]





p1 <- fviz_cluster(HC, data = dati,
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal(),
             main="Hierarchical, all ages")
p2 <- fviz_cluster(FC, data = dati,
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal(),
             main="Fuzzy, all ages")



dev.off()



res.fanny <- fanny(df, 2)
res.fanny$membership

fviz_cluster(res.fanny, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(),
             legend = "right")
fviz_silhouette(res.fanny, palette = "jco",
                ggtheme = theme_minimal())

fviz_dend(res.fanny, # Cut in nk groups
          cex = 1, # label size
          k_colors = mycol,
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE, # Add rectangle around groups
          horiz=TRUE,
          main="C-PCLM"
)
fviz_silhouette(res.fanny, palette = "jco",
                ggtheme = theme_minimal())
fviz_cluster(res.hk, data = dati,
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal())

fviz_nbclust(dati, hkmeans, method = "wss")

matplot(x, t(res.hk$centers), col=mycol[c(2,1,3,4)])

par(mfrow=c(2,1))
rany <- c(-19, -1)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main="Log-mortality", col.main=1, cex.main=1.5)
matlines(x, ETA, col=colocou, t="l", lty=1)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy, las=2)
axis(1);box()
abline(h=log(yy), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")

rany <- c(-0.2, 0.4)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main="Rate of aging", col.main=1, cex.main=1.5)
yy <- seq(-1,1,0.1)
axis(2, at=yy, labels=yy, las=2)
axis(1);box()
abline(h=seq(-1,1,0.1), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
abline(h=0, col = "gray", lty = 2, lwd=2)
matlines(x, t(dati), col=colocou, t="l", lty=1)

par(mfrow=c(2,1))
rany <- c(-19, -1)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main="Log-mortality", col.main=1, cex.main=1.5)

yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy, las=2)
axis(1);box()
abline(h=log(yy), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
mycol1 <- mycol[c(2,1,4,3)]
matlines(x, ETA, col=mycol1[res.hk$cluster], t="l", lty=1)

rany <- c(-0.2, 0.4)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main="Clustering Rate of aging", col.main=1, cex.main=1.5)
yy <- seq(-1,1,0.1)
axis(2, at=yy, labels=yy, las=2)
axis(1);box()
abline(h=seq(-1,1,0.1), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
abline(h=0, col = "gray", lty = 2, lwd=2)

mycol1 <- mycol[c(2,1,4,3)]
matlines(x, t(dati), col=mycol1[res.hk$cluster], t="l", lty=3)
matlines(x, t(res.hk$centers), col=mycol1, t="l", lwd=7, lty=1)
dev.off()

