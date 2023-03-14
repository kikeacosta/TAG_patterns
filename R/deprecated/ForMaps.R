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
library(dendextend)

source("EmpiricDeriv.R")

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

## clustering
res.hk3 <-hkmeans(dati, 3)
res.hk4 <-hkmeans(dati, 4)
res.hk5 <-hkmeans(dati, 5)


clusterB <- cbind(res.hk3$cluster, res.hk4$cluster, res.hk5$cluster)

## only adult ages

## what to cluster
dati <- t(ETAB1[31:71,]) 

## clustering
res.hk3 <-hkmeans(dati, 3)
res.hk4 <-hkmeans(dati, 4)
res.hk5 <-hkmeans(dati, 5)


clusterBadult <- cbind(res.hk3$cluster, res.hk4$cluster, res.hk5$cluster)


## males
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

## clustering
res.hk3 <-hkmeans(dati, 3)
res.hk4 <-hkmeans(dati, 4)
res.hk5 <-hkmeans(dati, 5)

clusterFM <- cbind(res.hk3$cluster, res.hk4$cluster, res.hk5$cluster)

## only adult ages

## what to cluster
dati <- t(ETAFM1[31:71,]) 

## clustering
res.hk3 <-hkmeans(dati, 3)
res.hk4 <-hkmeans(dati, 4)
res.hk5 <-hkmeans(dati, 5)


clusterFMadult <- cbind(res.hk3$cluster, res.hk4$cluster, res.hk5$cluster)


## hierarchical model
load("HierarOut.RDdata")

ETAF.hat <- matrix(eta[1:(p*m)], m, p)
ETAM.hat <- matrix(eta[1:(p*m)+p*m], m, p)
beta0 <- beta[1:nb]
eta0 <- as.matrix(B)%*%beta0
beta.pop <- beta[1:(p*nb)+nb]
Bpop <- kronecker(diag(p), B)
pi <- as.matrix(Bpop)%*%beta.pop
Pi <- matrix(pi, m, p)
beta.sex <- beta[1:nb+(p*nb)+nb]
si <- as.matrix(B)%*%beta.sex

## what to cluster
colnames(Pi) <- pop
Pinorm <- Pi
for(i in 1:p){
  Pinorm[,i] <- Pi[,i]-mean(Pi[,i])
}



dati <- t(Pinorm)

## clustering
res.hk3 <-hkmeans(dati, 3)
res.hk4 <-hkmeans(dati, 4)
res.hk5 <-hkmeans(dati, 5)


clusterH <- cbind(res.hk3$cluster, res.hk4$cluster, res.hk5$cluster)


## only adult
dati <- t(Pinorm[31:71,])

## clustering
res.hk3 <-hkmeans(dati, 3)
res.hk4 <-hkmeans(dati, 4)
res.hk5 <-hkmeans(dati, 5)


clusterHadult <- cbind(res.hk3$cluster, res.hk4$cluster, res.hk5$cluster)



CLUSTERS <- data.frame(clusterB, clusterBadult,
                       clusterFM, clusterFMadult,
                       clusterH, clusterHadult)
colnames(CLUSTERS) <- c(paste("C", 3:5, sep=""),
                        paste("Cadult", 3:5, sep=""),
                        paste("S", 3:5, sep=""),
                        paste("Sadult", 3:5, sep=""),
                        paste("A", 3:5, sep=""),
                        paste("Aadult", 3:5, sep=""))
                    
write.table(CLUSTERS, "CLUSTERS.txt")



