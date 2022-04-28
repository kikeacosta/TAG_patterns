## R-code for estimating log-mortality age-pattern 
## in two periods for a group of populations belonging to the same cluster
## (clusters provided by WM)
## for each population within a cluster we assume that 2019 log-mortality 
## is a smooth function over age  
## and log-mortality in 2020 is 
## estimated log-mortality in 2019 + 
## scaling factor +
## smooth cluster-specific age-factor
## by C.G. Camarda 2021.09.29


## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## !!!! to be changed
setwd("~/WORK/TAG_patterns/")
library(MortalitySmooth)
library(magic)
library(colorspace)
library(ggplot2)
library(plotly)


## loading baseline
data1 <- read.csv("Output/GHE2019_baseline.csv", header=TRUE)
## loading 2020
data2 <- read.csv("Output/WM2020_observed.csv", header=TRUE)
## removed doubles Colombia in data2
data2 <- data2[-c(2161:2268), ]

# whi.data1 <- which(data1$lmx==data2$lmx[37])
# whi.data2 <- which(data2$lmx==data2$lmx[37])
# 
# data1$deaths[whi.data1] <- 0
# data2$deaths[whi.data2] <- 0
# data1$lmx <- log(data1$deaths/data1$Nx)
# data2$lmx <- log(data2$deaths/data2$Nx)

## select pop which are in both datasets
## populations
pop <- as.character(sort(unique(data2$iso3)))
p <- length(pop)

## loading clusters
clu <- read.csv("Output/iso.clusters.csv", header=TRUE)
data1$cluster <- NA
data2$cluster <- NA
for(c in 1:8){
  clu.c <- subset(clu, Cluster==c)
  whi.c1 <- which(data1$iso3 %in% clu.c$iso3)
  data1$cluster[whi.c1] <- c
  whi.c2 <- which(data2$iso3 %in% clu.c$iso3)
  data2$cluster[whi.c2] <- c
}

data1$lmx <- log(data1$deaths/data1$Nx)
data2$lmx <- log(data2$deaths/data2$Nx)

data1$lmx.hat <- NA
data1$lmx.hatU <- NA
data1$lmx.hatL <- NA
data2$lmx.hat <- NA
data2$lmx.hatU <- NA
data2$lmx.hatL <- NA

## objects for scaling factor and delta
DELTAs <- expand.grid(age=unique(data1$age), sex=unique(data1$sex), cluster=as.factor(1:8))
DELTAs$delta.hat <- NA
DELTAs$delta.hatU <- NA
DELTAs$delta.hatL <- NA
Cs <- expand.grid(pop=pop, sex=unique(data1$sex))
Cs$c.hat <- NA
Cs$c.hatU <- NA
Cs$c.hatL <- NA

## ages
x <- c(seq(2.5, 82.5, 5), 90)
m <- length(x)

## over cluster
c <- 8
for(c in 1:8){
  
  whi.pop.c <- which(clu$Cluster==c)
  pop.c <- clu$iso3[whi.pop.c][which(clu$iso3[whi.pop.c]%in%pop)]
  p.c <- length(pop.c)
  ## data1
  whi.pop1F <- data1$iso3%in%pop.c & data1$sex=="Female"
  whi.pop1M <- data1$iso3%in%pop.c & data1$sex=="Male"
  data1F.c <- data1[whi.pop1F,]
  data1M.c <- data1[whi.pop1M,]
  
  ## data2
  whi.pop2F <- data2$iso3%in%pop.c & data2$sex=="Female"
  whi.pop2M <- data2$iso3%in%pop.c & data2$sex=="Male"
  data2F.c <- data2[whi.pop2F,]
  data2M.c <- data2[whi.pop2M,]
  
  ## B-splines
  B <- MortSmooth_bbase(x=x, xl=min(x), xr=max(x), ndx=8, deg=3)
  nb <- ncol(B)
  
  ## building model matrix
  ## for pop-specific log-mortality in 2019
  Xeta0 <- kronecker(diag(p.c), B)
  Xeta1 <- kronecker(rep(1,2), Xeta0)
  ## for pop-specific scaling factor
  Xc0 <- kronecker(diag(p.c), matrix(1,m))
  Xc <- rbind(0*Xc0, Xc0)
  ## for the cluster=specific delta
  Xdelta0 <- kronecker(rep(1,p.c), B)
  Xdelta <- rbind(0*Xdelta0, Xdelta0)
  ## final X
  X <- cbind(Xeta1, Xc, Xdelta)
  
  ## penalty stuff
  Deta1 <- diff(diag(nb), diff=2)
  tDDeta1 <- t(Deta1)%*%Deta1
  Ddelta <- diff(diag(nb), diff=2)
  tDDdelta <- t(Ddelta)%*%Ddelta
  lambda.eta1 <- 10^-1
  lambda.delta <- 10^3
  Peta1 <- kronecker(lambda.eta1*diag(p.c), tDDeta1)
  Pdelta <- lambda.delta*tDDdelta
  P <- adiag(Peta1, 0*diag(p.c), Pdelta)
  ## ridge penalty
  Pr <- 1e-4*diag(ncol(P))
  ## constraining delta to sum up to 0
  H <- matrix(c(rep(0,(nb+1)*p.c), rep(1,nb)), nrow=1)
  kappa <- 0
  
  ## FEMALES
  
  ## data in vector
  y <- c(data1F.c$deaths, data2F.c$deaths)
  e <- c(data1F.c$Nx,     data2F.c$Nx)
  
  ## starting values
  eta <- log((y+1)/(e+1))
  max.it <- 100
  for(it in 1:max.it){
    mu <- exp(eta)
    z <- eta + (y - e*mu)/(e*mu)
    w <- as.vector(e*mu)
    tXWX <- t(X)%*%(w*X)
    tXWXpP <- tXWX + P + Pr
    tXWz <- t(X)%*%(w*z)
    LHS <- rbind(cbind(tXWXpP, t(H)),
                 cbind(H, 0))
    RHS <- matrix(c(tXWz, kappa), ncol=1)
    coeff   <- solve(LHS, RHS)
    betas   <- coeff[1:((nb+1)*p.c+nb)]
    eta.old <- eta
    eta     <- X%*%betas
    dif.eta <- max(abs((eta - eta.old)/eta.old))
    if(dif.eta < 1e-04 & it > 4) break
    cat(it, dif.eta, "\n")
  }
  betas.hat <- betas
  eta1.hat <- eta[1:(m*p.c)]
  c.hat <- betas[1:p.c+nb*p.c]
  delta.hat <- B%*%betas[1:nb+p.c+nb*p.c]
  eta2.hat <- eta[1:(m*p.c)+m*p.c]
  
  ## standard errors
  H0 <- solve(tXWXpP)
  Vbetas <- H0 %*% tXWX %*% H0
  diagVbetas <- diag(Vbetas)
  se.betas <- sqrt(diagVbetas)
  se.betas[is.na(se.betas)] <- mean(se.betas[!is.na(se.betas)])
  se.c <- se.betas[1:p.c+nb*p.c]
  V.eta12 <- X %*% Vbetas %*% t(X)
  se.eta1 <- sqrt(diag(V.eta12))[1:(m*p.c)]
  se.eta2 <- sqrt(diag(V.eta12))[1:(m*p.c)+m*p.c]
  V.delta <- B %*% Vbetas[1:nb+p.c+nb*p.c, 1:nb+p.c+nb*p.c] %*% t(B)
  se.delta <- sqrt(diag(V.delta))
  se.delta[is.na(se.delta)] <- mean(se.delta[!is.na(se.delta)])
  
  ## confidence intervals
  eta1.hatL <- eta1.hat - 2*se.eta1
  eta1.hatU <- eta1.hat + 2*se.eta1
  eta2.hatL <- eta2.hat - 2*se.eta2
  eta2.hatU <- eta2.hat + 2*se.eta2
  c.hatL <- c.hat - 2*se.c
  c.hatU <- c.hat + 2*se.c
  delta.hatL <- delta.hat - 2*se.delta
  delta.hatU <- delta.hat + 2*se.delta
  
  data1$lmx.hat[whi.pop1F] <- eta1.hat
  data1$lmx.hatL[whi.pop1F] <- eta1.hatL
  data1$lmx.hatU[whi.pop1F] <- eta1.hatU
  data2$lmx.hat[whi.pop2F] <- eta2.hat
  data2$lmx.hatL[whi.pop2F] <- eta2.hatL
  data2$lmx.hatU[whi.pop2F] <- eta2.hatU
  
  whi.cpF <- which(Cs$pop%in%pop.c & Cs$sex=="Female")
  Cs[whi.cpF, ]$c.hat <- c.hat
  Cs[whi.cpF, ]$c.hatU <- c.hatU
  Cs[whi.cpF, ]$c.hatL <- c.hatL
  whi.deltaF <- which(DELTAs$cluster%in%c & DELTAs$sex=="Female")
  DELTAs[whi.deltaF, ]$delta.hat <- delta.hat
  DELTAs[whi.deltaF, ]$delta.hatU <- delta.hatU
  DELTAs[whi.deltaF, ]$delta.hatL <- delta.hatL
  
  ## MALES
  y <- c(data1M.c$deaths, data2M.c$deaths)
  e <- c(data1M.c$Nx, data2M.c$Nx)
  
  ## starting values
  eta <- log((y+1)/(e+1))
  max.it <- 100
  for(it in 1:max.it){
    mu <- exp(eta)
    z <- eta + (y - e*mu)/(e*mu)
    w <- as.vector(e*mu)
    tXWX <- t(X)%*%(w*X)
    tXWXpP <- tXWX + P + Pr
    tXWz <- t(X)%*%(w*z)
    LHS <- rbind(cbind(tXWXpP, t(H)),
                 cbind(H, 0))
    RHS <- matrix(c(tXWz, kappa), ncol=1)
    coeff   <- solve(LHS, RHS)
    betas   <- coeff[1:((nb+1)*p.c+nb)]
    eta.old <- eta
    eta     <- X%*%betas
    dif.eta <- max(abs((eta - eta.old)/eta.old))
    if(dif.eta < 1e-04 & it > 4) break
    cat(it, dif.eta, "\n")
  }
  betas.hat <- betas
  eta1.hat <- eta[1:(m*p.c)]
  c.hat <- betas[1:p.c+nb*p.c]
  delta.hat <- B%*%betas[1:nb+p.c+nb*p.c]
  eta2.hat <- eta[1:(m*p.c)+m*p.c]
  
  ## standard errors
  H0 <- solve(tXWXpP)
  Vbetas <- H0 %*% tXWX %*% H0
  diagVbetas <- diag(Vbetas)
  se.betas <- sqrt(diagVbetas)
  se.betas[is.na(se.betas)] <- mean(se.betas[!is.na(se.betas)])
  se.c <- se.betas[1:p.c+nb*p.c]
  V.eta12 <- X %*% Vbetas %*% t(X)
  se.eta1 <- sqrt(diag(V.eta12))[1:(m*p.c)]
  se.eta2 <- sqrt(diag(V.eta12))[1:(m*p.c)+m*p.c]
  V.delta <- B %*% Vbetas[1:nb+p.c+nb*p.c, 1:nb+p.c+nb*p.c] %*% t(B)
  se.delta <- sqrt(diag(V.delta))
  se.delta[is.na(se.delta)] <- mean(se.delta[!is.na(se.delta)])
  
  ## confidence intervals
  eta1.hatL <- eta1.hat - 2*se.eta1
  eta1.hatU <- eta1.hat + 2*se.eta1
  eta2.hatL <- eta2.hat - 2*se.eta2
  eta2.hatU <- eta2.hat + 2*se.eta2
  c.hatL <- c.hat - 2*se.c
  c.hatU <- c.hat + 2*se.c
  delta.hatL <- delta.hat - 2*se.delta
  delta.hatU <- delta.hat + 2*se.delta
  
  data1$lmx.hat[whi.pop1M] <- eta1.hat
  data1$lmx.hatL[whi.pop1M] <- eta1.hatL
  data1$lmx.hatU[whi.pop1M] <- eta1.hatU
  data2$lmx.hat[whi.pop2M] <- eta2.hat
  data2$lmx.hatL[whi.pop2M] <- eta2.hatL
  data2$lmx.hatU[whi.pop2M] <- eta2.hatU
  
  whi.cpM <- which(Cs$pop%in%pop.c & Cs$sex=="Male")
  Cs[whi.cpM, ]$c.hat <- c.hat
  Cs[whi.cpM, ]$c.hatU <- c.hatU
  Cs[whi.cpM, ]$c.hatL <- c.hatL
  whi.deltaM <- which(DELTAs$cluster%in%c & DELTAs$sex=="Male")
  DELTAs[whi.deltaM, ]$delta.hat <- delta.hat
  DELTAs[whi.deltaM, ]$delta.hatU <- delta.hatU
  DELTAs[whi.deltaM, ]$delta.hatL <- delta.hatL
}


save.image("Output/OutPrepandemic2020sex_StratifiedBySex_WMdata.Rdata")



## 2019
pdf("ActualFittedLogMortalityClusterDelta.pdf", width = 12, height = 10)
data1.obs <- subset(data1, !is.na(lmx.hat))
aa <- ggplot(data=data1.obs, aes(x=age, y=lmx, col=sex))+
  geom_point()+
  geom_line(aes(y=lmx.hat))+
  facet_wrap(~cluster+iso3, 6, 10, scales="free_y")+
  theme(axis.text.y = element_blank())+
labs(x = "age", y = "", title = "log-mortality, 2019")
plot(aa)

## 2020
aa <- ggplot(data=data2, aes(x=age, y=lmx, col=sex))+
  geom_point()+
  geom_line(aes(y=lmx.hat))+
  facet_wrap(~cluster+iso3, 6, 10, scales="free_y")+
  theme(axis.text.y = element_blank())+
  labs(x = "age", y = "", title = "log-mortality, 2020")
plot(aa)
# dev.off()
# ## deltas
# pdf("/home/gccamarda/WORK/TAG_patterns/Slides/Figures/ClusterDeltas.pdf", 
#    width=12, height = 8)
ggplot(DELTAs, aes(x = age, y = delta.hat)) +
  geom_smooth( aes(ymin = delta.hatL, ymax = delta.hatU, fill = sex, colour = sex),
    stat = "identity")+
  facet_wrap(~cluster, 2, 4)+
  geom_hline(yintercept=0,linetype=2)+
  theme()+
  labs(x = "age", y = expression(hat(delta)))
dev.off()





# 


DELTAS_WM <- expand.grid(age=seq(0,85, 5), sex=c("Female", "Male"), iso3=pop)
DELTAS_WM$deltas <- NA
j=1
for(j in 1:p){
  OUT.j <- OUT[[j]]
  deltas <- c(OUT.j$deltaF.hat, OUT.j$deltaM.hat)
  wr <- 1:(m*2)+(j-1)*(m*2)
  DELTAS_WM$deltas[wr] <- deltas
}
write.table(DELTAS_WM, "DELTAS_WM.txt")

## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## !!!! to be changed
setwd("~/WORK/TAG_patterns/")
library(MortalitySmooth)
library(magic)
library(colorspace)

## loading outcomes from model
load("Output/OutPrepandemic2020sex_StratifiedBySex_WMdata.Rdata")

## loading baseline
compare <- read.csv("Output/age_sex_compare.csv", header=TRUE)






