## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## R-studio to get the same dir as for the .R file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(magic)
library(colorspace)
library(MortalitySmooth)
library(factoextra)

path <- "/home/gccamarda/WORK/Applications/TAG_COVID19_UN/2021_04_09_TAG/ModelPattern/TAG_patterns/Slides/Figures/"

source("EmpiricDeriv.R")
library(dendextend)

## both
load("OUTCOMES/ALLIndivB.RData")
ETAB <- ETA
x <- x
## compute rate-of-aging
ETAB1 <- ETAB*0
for(i in 1:p){
  ETAB1[,i] <- emp.der(x=x, y=ETAB[,i])
}

## what to cluster
dati <- t(ETAB1) 

## decide the number clusters
nc <- 4
mycol <- c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07", "darkblue", "darkred")[1:4]
res.hk <-hkmeans(dati, nc)

pdf(paste(path, "CombDendo.pdf", sep=""), width = 6, height = 9)
fviz_dend(res.hk, k = nc, # Cut in nk groups
          cex = 1, # label size
          k_colors = mycol,
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE, # Add rectangle around groups
          horiz=TRUE,
          main="C-PCLM"
)
dev.off()
fviz_cluster(res.hk, data = dati,
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal())

fviz_nbclust(dati, hkmeans, method = "wss")

matplot(x, t(res.hk$centers), col=mycol[c(2,1,3,4)])

pdf(paste(path, "CombData.pdf", sep=""), width = 8, height = 9)
par(mfrow=c(2,1))
rany <- c(-19, -1)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main="Log-mortality", col.main=1, cex.main=1.5)
matlines(x, ETAB, col=colocou, t="l", lty=1)
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
dev.off()

pdf(paste(path, "CombClust.pdf", sep=""), width = 8, height = 12)
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
matlines(x, ETAB, col=mycol1[res.hk$cluster], t="l", lty=1)

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


## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## R-studio to get the same dir as for the .R file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(magic)
library(colorspace)
library(MortalitySmooth)
library(factoextra)

path <- "/home/gccamarda/WORK/Applications/TAG_COVID19_UN/2021_04_09_TAG/ModelPattern/TAG_patterns/Slides/Figures/"

source("EmpiricDeriv.R")
library(dendextend)

## males
load("OUTCOMES/ALLIndivB.RData")
load("OUTCOMES/ALLIndivM.RData")
ETAM <- ETA
x <- x
## compute rate-of-aging
ETAM1 <- ETAM*0
for(i in 1:p){
  ETAM1[,i] <- emp.der(x=x, y=ETAM[,i])
}

## females
load("OUTCOMES/ALLIndivF.RData")
ETAF <- ETA
x <- x
## compute rate-of-aging
ETAF1 <- ETAF*0
for(i in 1:p){
  ETAF1[,i] <- emp.der(x=x, y=ETAF[,i])
}

ETAFM1 <- rbind(ETAF1,ETAM1)

## what to cluster
dati <- t(ETAFM1) 


## decide the number clusters
nc <- 4
mycol <- c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07", "darkblue", "darkred")[1:4]
res.hk <-hkmeans(dati, nc)

pdf(paste(path, "SexDendo.pdf", sep=""), width = 6, height = 9)
fviz_dend(res.hk, k = nc, # Cut in nk groups
          cex = 1, # label size
          k_colors = mycol,
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE, # Add rectangle around groups
          horiz=TRUE,
          main="S-PCLM"
)
dev.off()
fviz_cluster(res.hk, data = dati,
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal())

fviz_nbclust(dati, hkmeans, method = "wss")


pdf(paste(path, "SexData.pdf", sep=""), width = 14, height = 9)
par(mfrow=c(2,2))
rany <- c(-19, -1)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main="Log-mortality, females", col.main=1, cex.main=1.5)
matlines(x, ETAF, col=colocou, t="l", lty=1)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy, las=2)
axis(1);box()
abline(h=log(yy), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")

rany <- c(-19, -1)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main="Log-mortality, males", col.main=1, cex.main=1.5)
matlines(x, ETAM, col=colocou, t="l", lty=1)
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
title(main="Rate of aging, females", col.main=1, cex.main=1.5)
yy <- seq(-1,1,0.1)
axis(2, at=yy, labels=yy, las=2)
axis(1);box()
abline(h=seq(-1,1,0.1), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
abline(h=0, col = "gray", lty = 2, lwd=2)

matlines(x, ETAF1, col=colocou, t="l", lty=1)

rany <- c(-0.2, 0.4)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main="Rate of aging, males", col.main=1, cex.main=1.5)
yy <- seq(-1,1,0.1)
axis(2, at=yy, labels=yy, las=2)
axis(1);box()
abline(h=seq(-1,1,0.1), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
abline(h=0, col = "gray", lty = 2, lwd=2)

matlines(x, ETAM1, col=colocou, t="l", lty=1)

dev.off()

pdf(paste(path, "SexClust.pdf", sep=""), width = 8, height = 12)
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
matlines(x, ETAB, col=mycol1[res.hk$cluster], t="l", lty=1)

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












## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## R-studio to get the same dir as for the .R file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(magic)
library(colorspace)
library(MortalitySmooth)
library(factoextra)

path <- "/home/gccamarda/WORK/Applications/TAG_COVID19_UN/2021_04_09_TAG/ModelPattern/TAG_patterns/Slides/Figures/"
load("OUTCOMES/ALLIndivB.RData")
source("EmpiricDeriv.R")
library(dendextend)

## males
load("HierarOut.RDdata")

ETAF.hat <- matrix(eta[1:(p*m)], m, p)
ETAM.hat <- matrix(eta[1:(p*m)+p*m], m, p)
beta0 <- beta[1:nb]
eta0 <- B%*%beta0
beta.pop <- beta[1:(p*nb)+nb]
Bpop <- kronecker(diag(p), B)
pi <- Bpop%*%beta.pop
Pi <- matrix(pi, m, p)
beta.sex <- beta[1:nb+(p*nb)+nb]
si <- as.matrix(B)%*%beta.sex

colnames(Pi) <- pop
Pinorm <- Pi
for(i in 1:p){
  Pinorm[,i] <- Pi[,i]-mean(Pi[,i])
}

pdf(paste(path, "HierOutGomp.pdf", sep=""), width = 14, height = 7)
par(mfrow=c(1,3))
rany <- c(-19, -1)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main="Gompertz reference log-mortality", col.main=1, cex.main=1.5)

yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy, las=2)
axis(1);box()
abline(h=log(yy), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
lines(x, eta0, col="darkred", t="l", lty=1, lwd=4)


rany <- c(-0, 0.8)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main="Sex factor over age", col.main=1, cex.main=1.5)
yy <- seq(-1,1,0.1)
axis(2, at=yy, labels=yy, las=2)
axis(1);box()
abline(h=seq(-1,1,0.1), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
abline(h=0, col = "gray", lty = 2, lwd=2)

lines(x, si, col="darkblue", lwd=4, lty=1)

rany <- c(-8, 4)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main="Population-specific deviance", col.main=1, cex.main=1.5)
yy <- seq(-10,10,1)
axis(2, at=yy, labels=yy, las=2)
axis(1);box()
abline(h=seq(-10,10,1), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
abline(h=0, col = "gray", lty = 2, lwd=2)

matlines(x, Pi, col=colocou, lwd=1, lty=1)
dev.off()


colnames(Pi) <- pop



## what to cluster
dati <- t(Pinorm) 


## decide the number clusters
nc <- 4
mycol <- c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07", "darkblue", "darkred")[1:4]
res.hk <-hkmeans(dati, nc)

pdf(paste(path, "HierDendoGomp.pdf", sep=""), width = 6, height = 9)
fviz_dend(res.hk, k = nc, # Cut in nk groups
          cex = 1, # label size
          k_colors = mycol,
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE, # Add rectangle around groups
          horiz=TRUE,
          main="A-PCLM with Gompertz reference"
)
dev.off()
fviz_cluster(res.hk, data = dati,
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal())

fviz_nbclust(dati, hkmeans, method = "wss")


pdf(paste(path, "HierClust.pdf", sep=""), width = 9, height = 9)
rany <- c(-4, 4)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main="Population-specific deviance", col.main=1, cex.main=1.5)
yy <- seq(-10,10,1)
axis(2, at=yy, labels=yy, las=2)
axis(1);box()
abline(h=seq(-10,10,1), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
abline(h=0, col = "gray", lty = 2, lwd=2)
mycol1 <- mycol[c(1,2,3,4)]
matlines(x, Pinorm, col=mycol1[res.hk$cluster], t="l",lty=3)
matlines(x, t(res.hk$centers), col=mycol1, lwd=7, lty=1)
#lines(x, Pi[,3], col=2, lwd=4, lty=2)

dev.off()































pdf(paste(path, "SexData.pdf", sep=""), width = 14, height = 9)
par(mfrow=c(2,2))
rany <- c(-19, -1)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main="Log-mortality, females", col.main=1, cex.main=1.5)
matlines(x, ETAF, col=colocou, t="l", lty=1)
yy <- 10^seq(-10, 0)
axis(2, at=log(yy), labels=yy, las=2)
axis(1);box()
abline(h=log(yy), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")

rany <- c(-19, -1)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main="Log-mortality, males", col.main=1, cex.main=1.5)
matlines(x, ETAM, col=colocou, t="l", lty=1)
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
title(main="Rate of aging, females", col.main=1, cex.main=1.5)
yy <- seq(-1,1,0.1)
axis(2, at=yy, labels=yy, las=2)
axis(1);box()
abline(h=seq(-1,1,0.1), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
abline(h=0, col = "gray", lty = 2, lwd=2)

matlines(x, ETAF1, col=colocou, t="l", lty=1)

rany <- c(-0.2, 0.4)
ranx <- range(x)
par(mar=c(3,4,3,1))
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,xlab="", ylab="")
title(main="Rate of aging, males", col.main=1, cex.main=1.5)
yy <- seq(-1,1,0.1)
axis(2, at=yy, labels=yy, las=2)
axis(1);box()
abline(h=seq(-1,1,0.1), col = "lightgray", lty = "dotted")
abline(v=seq(0,120,20), col = "lightgray", lty = "dotted")
abline(h=0, col = "gray", lty = 2, lwd=2)

matlines(x, ETAM1, col=colocou, t="l", lty=1)

dev.off()

pdf(paste(path, "SexClust.pdf", sep=""), width = 8, height = 12)
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
matlines(x, ETAB, col=mycol1[res.hk$cluster], t="l", lty=1)

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
matlines(x, t(dati), col=mycol1[res.hk$cluster], t="l", lwd=0.3, lty=3)
matlines(x, t(res.hk$centers), col=mycol1, t="l", lwd=7, lty=1)
dev.off()


























res.hc <- hclust(dist(dati,
                      method = "euclidean"),
                 method = "ward.D2")

mycol <- c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07", "darkblue", "darkred")
# plot(res.hc)
# par(mfrow=c(1,2))
pdf(paste(path, "CombDendo.pdf", sep=""), width = 6, height = 9)
fviz_dend(res.hc, k = nc, # Cut in nk groups
          cex = 1, # label size
          k_colors = mycol[1:nc],
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE, # Add rectangle around groups
          horiz=TRUE
)
dev.off()

## 
matplot(x, t(km.res$centers), col=mycol[1:nc])



## plotting rate-of-aging
Centers <- matrix(NA, m, nc)
for(i in 1:nc){
   Centers[,i] <- fit.i$B %*% km.res$centers[i,]
}



## males
load("OUTCOMES/ALLIndivM.RData")
ETAM <- ETA
x <- x
## compute rate-of-aging
ETAM1 <- ETAM*0
for(i in 1:nc){
  ETAM1[,i] <- emp.der(x=x, y=ETAM[,i])
}

## females
load("OUTCOMES/ALLIndivF.RData")
ETAF <- ETA
x <- x
## compute rate-of-aging
ETAF1 <- ETAF*0
for(i in 1:nc){
  ETAF1[,i] <- emp.der(x=x, y=ETAF[,i])
}

par(mfrow=c(2,2))
matplot(x, ETAM, t="l", lty=1)
matplot(x, ETAF, t="l", lty=1)

matplot(x, ETAM1, t="l", lty=1)
matplot(x, ETAF1, t="l", lty=1)
par(mfrow=c(1,1))

ETAFM1 <- rbind(ETAF1,ETAM1)

## what to cluster
dati <- t(ETAFM1) 

res.dist <- get_dist(dati, stand = TRUE, method = "pearson")
fviz_dist(res.dist, 
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
fviz_nbclust(dati, kmeans, method = "wss")

## decide the number clusters
nc <- 3
set.seed(123)
km.res <- kmeans(dati, nc, nstart = 25)
# Visualize
pdf("ClusterRateAging.pdf")
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

E <- subset(E0, Country%in%unique(YFM$Country))

## delete "England and Wales" "Moldova" "United Kingdom"
## from deaths since no exposures available
whiNA <- c("England and Wales", "Moldova", "United Kingdom")
whi <- !YFM$Country%in%whiNA
YFM <- YFM[whi,]

#cbind(sort(unique(YFM$Country)),sort(unique(E$Country)))

## delete "Ecuador" since no exposures by sex is available
YFM <- subset(YFM, Country!="Ecuador")
E <- subset(E, Country!="Ecuador")

# cbind(sort(unique(YFM$Country)),sort(unique(E$Country)))

## !!!!!! select sex
sex <- "M"
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
length(unique(YF$Country))

## ages for the latent pattern
x <- 0:104
m <- length(x)
PLOT=TRUE

REGIONS <- list(EUR=EUR, AFR=AFR, ASI=ASI, LAM=LAM, NAM=NAM, OCE=OCE)
k=4
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
    e.i <- E.i$Exposures
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
  nome <- paste("OUTCOMES/", names(REGIONS)[k], "indivM.RData", sep="")
  save.image(nome)
}

##


















## END