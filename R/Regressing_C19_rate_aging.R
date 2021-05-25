## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## R-studio to get the same dir as for the .R file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(magic)
library(MortalitySmooth)
library(colorspace)
library(gridExtra)
library(scales)
library(PerformanceAnalytics)
library(mgcv)
## loading fitted data by PCLM
load("FittedData.RData")


ETAs <- ETAs[,-c(40,49)]
ETA1s <- ETA1s[,-c(40,49)]
pop <- pop[-c(40,49)]
colocou <- colocou[-c(40,49)]
p <- length(pop)

## loading covariates
X0 <- readRDS("WPP_Covariates.rds", refhook = NULL)
## select for countries with available C19 profile
X1 <- subset(X0, Country%in%pop)
## only both sex data
X2 <- subset(X1, Sex=="b")

## remove Northern Ireland and Scotland
which(unique(X0$Country)=="Northern Ireland")
which(unique(X0$Country)=="Scotland")



## sort both datasets 
## sorting by name
ETAs <- ETAs[,order(pop)]
ETA1s <- ETA1s[,order(pop)]
X2 <- X2[order(X2$Country),]
pop <- pop[order(pop)]
colocou <- colocou[order(pop)]


par(mfrow=c(2,1))
matplot(xs, ETAs, col=colocou, t="l", lty=1)
matplot(xs, ETA1s, col=colocou, t="l", lty=1)
par(mfrow=c(1,1))


SVD <- svd(ETA1s)

U <- SVD$u
V <- SVD$v
s <- SVD$d

## take the first k singular values/vectors
k <- 1

## plot first k components
colok <- rainbow_hcl(k)
Ck <- U[,1]
Vk <- V[,1]
names(Vk) <- pop
sk <- s[1]

## plotting first k components
mylwd <- rescale(sk, c(3,15))
matplot(xs, Ck, t="l", lty=1, lwd=mylwd, col=colok)

## fitted values by svd with k components
ETA1s.hat1 <- (Ck*sk)%*%t(Vk)


## take the first k singular values/vectors
k <- 2

## plot first k components
colok <- rainbow_hcl(k)
Ck <- U[,1:k]
Vk <- V[,1:k]
rownames(Vk) <- pop
sk <- s[1:k]

## plotting first k components
mylwd <- rescale(sk, c(3,15))
matplot(xs, Ck, t="l", lty=1, lwd=mylwd, col=colok)

## fitted values by svd with k components
ETA1s.hat2 <- Ck%*%diag(sk)%*%t(Vk)

## take the first k singular values/vectors
k <- 3

## plot first k components
colok <- rainbow_hcl(k)
Ck <- U[,1:k]
Vk <- V[,1:k]
rownames(Vk) <- pop
sk <- s[1:k]

## plotting first k components
#png("SingularVectors.png", width = 960, height = 540)
mylwd <- rescale(sk^2, c(3,15))
matplot(xs, Ck, t="l", lty=1, lwd=mylwd, col=colok,
        xlab="ages", ylab="singular vectors",
        cex.lab=1.5, cex.axis=1.3)
title(main="First 3 singular vectors from Rate of Aging",
      cex.main=2)
grid()

legend("bottomright", inset=0.05,
       legend=c("Vector 1 ~ 97% (Gompertzian)", 
                "Vector 2 ~ 1.8% (Deceleration)",
                "Vector 3 ~ 1% (Sinusoid)"),
       col=colok, lwd=mylwd, bg="grey90", cex=1.3)
#dev.off()

## fitted values by svd with k components
ETA1s.hat3 <- Ck%*%diag(sk)%*%t(Vk)

#write.table(Vk, "SingularValues.txt")
((sk^2)/sum(s^2))*100



## plotting PCLM-fitted and svd-reduced rate-of-aging
rany <- range(ETA1s)
whi <- c(7, 13, 20, 35, 58, 60)
#png("FittedSVD.png", width = 960, height = 540)
par(mfrow=c(2,3))
for(i in whi){
  #rany <- range(ETA1s[,i], ETA1s.hat2[,i], ETA1s.hat3[,i])
  plot(xs, ETA1s[,i], col="grey60", pch=16, cex=1.2, ylim=rany,
       xlab="ages", ylab="Rate of Aging",
       cex.lab=1.5, cex.axis=1.3)
  grid()
  lines(xs, ETA1s.hat2[,i], col="darkred", lwd=4)
  lines(xs, ETA1s.hat3[,i], col="darkblue", lwd=4)
  title(main=pop[i], col.main=colocou[i], cex.main=1.5)
  legend("bottomleft", inset=0.01,
         legend=c("Rate-of-aging", "Fitted with 2 PC", "Fitted with 3 PC"),
         col=c("grey60","darkred","darkblue"), lty=c(NA,1,1), 
         lwd=c(1,4,4), pch=c(16,NA,NA),
         pt.cex=1.2, cex=1.2, bg="grey90")
#  locator(1)
}
par(mfrow=c(1,1))
#dev.off()

## plotting country-specific k singular values
mycex <- rescale(sk, c(1,2))
matplot(Vk, cex=mycex, pch=16:18)



## data for regression
dati <- data.frame(sv1=Vk[,1],
                   sv2=Vk[,2],
                   sv3=Vk[,3],
                   e0=X2$e0,
                   ed=X2$edagger,
                   b=X2$B,
                   roa=X2[,-c(1:7)])
chart.Correlation(dati)
library(Hmisc)


# png("LifeExp.png", width = 960, height = 540)
# par(mfrow=c(1,3))

## for TR
write.table(dati, "ResponsesCovariates.txt")


#pdf("RegressionSVD.pdf", width = 14, height = 10)
par(mfrow=c(3,3))

plot(dati$e0, dati$sv1, xlab="life expectancy", ylab="1st singular values", t="n",
     cex.lab=1.5, cex.axis=1.3)
grid()
text(dati$e0, dati$sv1, pop, col=colocou)
lines(lowess(dati$sv1~dati$e0), col="darkred", lwd=4)

plot(dati$e0, dati$sv2, xlab="life expectancy", ylab="2nd singular values", t="n",
     cex.lab=1.5, cex.axis=1.3)
grid()
text(dati$e0, dati$sv2, pop, col=colocou)
lines(lowess(dati$sv2~dati$e0), col="darkred", lwd=4)

plot(dati$e0, dati$sv3, xlab="life expectancy", ylab="3rd singular values", t="n",
     cex.lab=1.5, cex.axis=1.3)
grid()
text(dati$e0, dati$sv3, pop, col=colocou)
lines(lowess(dati$sv3~dati$e0), col="darkred", lwd=4)

# dev.off()



# png("Edagger.png", width = 960, height = 540)
# par(mfrow=c(1,3))

plot(dati$ed, dati$sv1, xlab="e-dagger", ylab="1st singular values", t="n",
     cex.lab=1.5, cex.axis=1.3)
grid()
text(dati$ed, dati$sv1, pop, col=colocou)
lines(lowess(dati$sv1~dati$ed), col="darkred", lwd=4)

plot(dati$ed, dati$sv2, xlab="e-dagger", ylab="2nd singular values", t="n",
     cex.lab=1.5, cex.axis=1.3)
grid()
text(dati$ed, dati$sv2, pop, col=colocou)
lines(lowess(dati$sv2~dati$ed), col="darkred", lwd=4)

plot(dati$ed, dati$sv3, xlab="e-dagger", ylab="3rd singular values", t="n",
     cex.lab=1.5, cex.axis=1.3)
grid()
text(dati$ed, dati$sv3, pop, col=colocou)
lines(lowess(dati$sv3~dati$ed), col="darkred", lwd=4)


#dev.off()

#png("bGompertz.png", width = 960, height = 540)
#par(mfrow=c(1,3))

plot(dati$b, dati$sv1, xlab="Gompertz b", ylab="1st singular values", t="n",
     cex.lab=1.5, cex.axis=1.3)
grid()
text(dati$b, dati$sv1, pop, col=colocou)
lines(lowess(dati$sv1~dati$b), col="darkred", lwd=4)

plot(dati$b, dati$sv2, xlab="Gompertz b", ylab="2nd singular values", t="n",
     cex.lab=1.5, cex.axis=1.3)
grid()
text(dati$b, dati$sv2, pop, col=colocou)
lines(lowess(dati$sv2~dati$b), col="darkred", lwd=4)

plot(dati$b, dati$sv3, xlab="Gompertz b", ylab="3rd singular values", t="n",
     cex.lab=1.5, cex.axis=1.3)
grid()
text(dati$b, dati$sv3, pop, col=colocou)
lines(lowess(dati$sv3~dati$b), col="darkred", lwd=4)

#dev.off()


## independent regression for each singular values

## simple modelling using covariates


#par(mfrow=c(3,3))
m1 <- gam(sv1 ~ s(e0) + s(ed) + s(b), data=dati)
#plot(m1)
m2 <- gam(sv2 ~ s(e0) + s(ed) + b, data=dati)
#plot(m2)
m3 <- gam(sv3 ~ e0+ed, data=dati)
#plot(m3)
#par(mfrow=c(1,1))
## fitted values
Vk.hat <- cbind(m1$fitted.values, 
                m2$fitted.values,
                m3$fitted.values)
## fitted values by linear model on first k svd-values 
ETA1s.hat.mod3 <- Ck%*%diag(sk)%*%t(Vk.hat)
ETA1s.hat.mod2 <- Ck[,1:2]%*%diag(sk[1:2]) %*% t(Vk.hat[,1:2])
ETA1s.hat.mod1 <- (Ck[,1]*sk[1]) %*% t(Vk.hat[,1])

matplot(xs, ETA1s.hat.mod3, col=colocou, t="l", lty=1)

ETA1array <- array(0, dim=c(length(x1:x2), p, 7),
                   dimnames = list(x1:x2, pop, c("Obs", "SVD1", "SVD2", "SVD3", "Mod1", "Mod2", "Mod3")))

ETA1array[,,1] <- ETA1s
ETA1array[,,2] <- ETA1s.hat1
ETA1array[,,3] <- ETA1s.hat2
ETA1array[,,4] <- ETA1s.hat3
ETA1array[,,5] <- ETA1s.hat.mod1
ETA1array[,,6] <- ETA1s.hat.mod2
ETA1array[,,7] <- ETA1s.hat.mod3

dim(ETA1array)
save(ETA1array,file = "ETA1fit.rda")
rm(list = ls())
bla <- load("ETA1fit.rda")


rany <- range(ETA1s)
par(mfrow=c(3,5))
i=1
par(mar=c(2,2,2,0))
for(i in 1:p){
  plot(xs, ETA1s[,i], col=1, t="l", lwd=4, ylim=rany, xlab="", ylab="rate-of-aging")
  #lines(xs, ETA1s.hat3[,i], col=2, lwd=4)
  lines(xs, ETA1s.hat.mod[,i], col=4, lwd=4)
  title(main=pop[i], col.main=colocou[i])
  legend("bottomleft", inset=0.1,
         legend=c("Actual Rate-of-aging", "SVD approx", "MODEL approx on SVD"),
         col=c(1,2,4), lwd=4)
  #locator(1)
}

dev.off()

plot(x2, Vk[,1])


gam1 <- gam(Vk[,1]~s(x1))

Vk1.e0.hat <- gam1$fitted.values

plot(x1, Vk[,1], xlab="e0", ylab="1st sing. val", t="n")
text(x1, Vk[,1], pop, col=colocou)
points(x1, Vk1.e0.hat, col=2, lwd=4)

summary(lm1)


