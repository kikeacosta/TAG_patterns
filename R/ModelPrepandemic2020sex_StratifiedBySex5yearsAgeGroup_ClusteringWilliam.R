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
## by C.G. Camarda 2021.09.15


## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## !!!! to be changed
setwd("~/WORK/TAG_patterns/")
library(MortalitySmooth)
library(magic)
library(colorspace)

load("Output/OutPrepandemic2020sex_StratifiedBySexAgeGroup.Rdata")


WMout <- read.csv("Output/WM_constraints_and_clusters.csv", header=TRUE)


nc <- length(table(WMout$Cluster))


cbind(names(OUT), 
      sort(WMout$iso3[WMout$iso3%in%names(OUT)]))

## appearing in OUT, not in WMout, at least with the same name
## "MAC"
## "PYF"
## "TWN"
OUT <- OUT[-c(40,54,61)]
p <- length(OUT)
pop <- pop[-c(40,54,61)]

cbind(sort(names(OUT)), 
      sort(WMout$iso3[WMout$iso3%in%names(OUT)]))

## take for each cluster, population in which we have data
## and average the coefficients betas associated to deltas
betas.delta.cluF <- betas.delta.cluM <- matrix(NA, nb, nc)

i=1
for(i in 1:nc){
  ## pop in a given cluster
  cou.i <- WMout$iso3[WMout$Cluster==i]
  ## take those with data
  cou.data.i <- names(OUT)[names(OUT)%in%cou.i]
  cat(i, cou.data.i, "\n")
  ## take the betas associated to deltas from cou.data.i
  out.i <- OUT[cou.data.i]
  betasF.i <- betasM.i <- matrix(NA, nb, length(out.i))
  for(j in 1:length(out.i)){
    betasF.i[,j] <- out.i[[j]]$betasF.hat[1:nb+nb+1]
    betasM.i[,j] <- out.i[[j]]$betasM.hat[1:nb+nb+1]
  }
  betas.delta.cluF[,i] <- rowMeans(betasF.i)
  betas.delta.cluM[,i] <- rowMeans(betasM.i)
}

## ages at which we wanna evaluate outcomes
xs <- c(0, 2, seq(7.5,97.5,5))
ms <- length(xs)
Bs <- MortSmooth_bbase(x=xs, xl=min(x), xr=max(x), ndx=15, deg=3)
deltasF.clu <- matrix(NA, ms, nc)
rownames(deltasF.clu) <- c(0,1,seq(5,95,5))
colnames(deltasF.clu) <- 1:8
deltasM.clu <- deltasF.clu
for(i in 1:nc){
  deltasF.clu[,i] <- Bs%*%betas.delta.cluF[,i]
  deltasM.clu[,i] <- Bs%*%betas.delta.cluM[,i]
}

par(mfrow=c(1,2))
matplot(xs, deltasF.clu, t="l")
matplot(xs, deltasM.clu, t="l")

write.table(deltasF.clu, "deltasF.txt")
write.table(deltasM.clu, "deltasM.txt")

deltasF.obs <- matrix(NA, ms, p)
rownames(deltasF.obs) <- c(0,1,seq(5,95,5))
colnames(deltasF.obs) <- pop
deltasM.obs <- deltasF.obs
for(i in 1:p){
  betasF.i <- OUT[[i]]$betasF.hat[1:nb+nb+1]
  betasM.i <- OUT[[i]]$betasM.hat[1:nb+nb+1]
  deltasF.obs[,i] <- Bs%*%betasF.i
  deltasM.obs[,i] <- Bs%*%betasM.i
}

write.table(deltasF.obs, "deltasFobs.txt")
write.table(deltasM.obs, "deltasMobs.txt")


