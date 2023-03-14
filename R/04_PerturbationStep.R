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
library(viridis)

## for a single country: Colombia (18)

## loading deaths
deaths <- read.csv("data_inter/deaths_sourced_infant_based.csv", header=TRUE)
## loading population
offset <- read.csv("data_inter/offsets.csv", header=TRUE)
## sex
sex <- "m"
sexLT <- "M"
deaths <- subset(deaths, Sex==sex)
offset <- subset(offset, Sex==sex)

## ages and dimensions
a <- 0:100
m <- length(a)
## pandemic years for the first forecasting part
t2 <- 2020:2021
n2 <- length(t2)
## countries
cou <- unique(deaths$Country)
nc <- length(cou)

j=34
cou.j <- cou[j]
## select the country from deaths
deaths.j <- subset(deaths, Country==cou.j)
## check available pre-pandemic years, country-dependent
t1.ava <- unique(deaths.j$Year)[unique(deaths.j$Year)<t2[1]]
## missing years
t1 <- min(t1.ava):max(t1.ava)
n1 <- length(t1)
whi.ava <- t1%in%t1.ava
## select t1 from deaths.j
Y.j <- subset(deaths.j, Year%in%t1.ava)
## all years
t <- c(t1,t2)
n <- n1+n2
ag.low <- unique(Y.j$Age)
ag.up1  <- ag.low-1
if(ag.low[2]==1){
  ag.up <- c(0, ag.up1[ag.up1>0], max(a))
}else{
  ag.up <- c(ag.up1[ag.up1>0], max(a))
}
cbind(ag.low, ag.up)
lg <- ag.up-ag.low+1
mg <- length(ag.low)
## only for plotting/aesthetic
ag.mid <- (ag.low+ag.up)/2 
ag.lab <- paste(ag.low, ag.up,sep="-")
Yg1 <- matrix(0, mg, n1)
Yg1[,whi.ava] <- matrix(Y.j$Deaths, mg, length(t1.ava))
## working on the offset
E.j <- subset(offset, Country==cou.j & Year%in%t1.ava)
E1 <- matrix(0, m, n1)
E1[,whi.ava] <- matrix(E.j$Population, m, length(t1.ava))
## matrix to map grouped to single age
G <- matrix(0, mg, m)
rownames(G) <- ag.mid
colnames(G) <- a
for(i in 1:mg){
  ag.low.i <- ag.low[i]
  ag.up.i <- ag.up[i]
  wc <- which(a>=ag.low.i & a<=ag.up.i)
  G[i,wc] <- 1
}
all(colSums(G)==1)
## grouping exposures only for plotting
Eg1 <- G%*%E1
## observed log-mortality, by age-groups, only for plotting and e0
ETAg1 <- log(Yg1/Eg1)


source("~/WORK/ConstrainedMortalityForecast/GroupedData/FUNCTIONS/GroupedSmoothConstrainedMortalityForecasting_Functions.R")

## weights for grouped-exposures
WEIg1 <- matrix(1, mg,n1)
WEIg1[,!whi.ava] <- 0
## there is age 0?
infantj <- ifelse(ag.up[1]==0, TRUE, FALSE)

BICfun <- function(par){
  FIT <- PSinfantGrouped(a.low = ag.low,
                         Yg=Yg1,
                         a=a,
                         E=E1,
                         WEIg=WEIg1,
                         lambdas=par,
                         verbose=TRUE,
                         kappa.shape=0,
                         infant=infantj)
  FIT$bic
}
## optimizing lambdas using greedy grid search
tim1 <- Sys.time()
OPT <- cleversearch(BICfun, lower=c(-1, 1), upper=c(3, 5),
                    ngrid=5, logscale=TRUE, startvalue=c(10^4,10^5),
                    verbose=TRUE)
tim2 <- Sys.time()
tim2-tim1


## estimating mortality with optimal lambdas
FITi <- PSinfantGrouped(a.low = ag.low,
                        Yg=Yg1,
                        a=a,
                        E=E1,
                        WEIg=WEIg1,
                        lambdas=OPT$par,
                        verbose=TRUE,
                        kappa.shape=0,
                        infant=infantj)

## estimated/smooth log-mortality on the pre-pandemic years
ETA1hat <- FITi$ETA

## FORECASTING STEP ## FORECASTING STEP ## FORECASTING STEP 
## where to apply the constraints
S <- matrix(1, m, n)
S[,1:n1] <- 0
## compute deltas : CIs of relative derivatives
deltasi <- deltasFUN(FITi)
## forecasting
CPSg.i <- CPSfunctionGrouped(a.low=ag.low, Yg1=Yg1, WEIg1=WEIg1,
                             a=a, E1=E1, obs.yrs=t1, 
                             for.hor=max(t),
                             ETA1hat=ETA1hat, deltas=deltasi, S=S, 
                             lambdas = OPT$par, verbose = TRUE,
                             infant = infantj,
                             kappa.shape=0)

ETAhat <- CPSg.i$ETA
str(CPSg.i)
## analytic confidence intervals
B <- kronecker(CPSg.i$Bt, CPSg.i$Ba)
Veta <- B %*% CPSg.i$Valpha %*% t(B)
SEeta <- matrix(sqrt(diag(Veta)), m, n)
psi2 <- CPSg.i$dev/(m * n1 - CPSg.i$ed)

ETAup <- ETAhat + 2*sqrt(psi2)*SEeta
ETAlow <- ETAhat - 2*sqrt(psi2)*SEeta


## plotting, starting from the beginning
## actual log-mortality
DFg <- expand.grid(list(ages=ag.low, years1=t1))
DFg$type <- "Actual grouped"
DFg$eta1 <- c(ETAg1)
DFg$ages.up <- ag.up+1
DFg$eta1.up <- c(ETAg1)
DFg$eta1.low <- NA
## fitted
DFhat <- expand.grid(list(ages=a, years1=t1))
DFhat$type <- "Fitted"
DFhat$eta1 <- c(ETA1hat)
DFhat$ages.up <- NA
DFhat$eta1.up <- NA
DFhat$eta1.low <- NA
## forecast
DFfor <- expand.grid(list(ages=a, years1=t))
DFfor$type <- "Fitted+Forecast"
DFfor$eta1 <- c(ETAhat)
DFfor$ages.up <- NA
DFfor$eta1.up <- c(ETAup)
DFfor$eta1.low <- c(ETAlow)

## all
DF <- rbind(DFg, DFhat, DFfor) 
DF$group <- factor(c(rep(ag.lab, n1),rep(rep(ag.lab, lg), n1),
                     rep(rep(ag.lab, lg), n)), levels = ag.lab)
DF$colo <- c(rep(ag.mid, n1), rep(a, n1), rep(a, n))

## over age
p <- ggplot(DF, aes(x=ages, y=eta1, color=type)) +
  geom_segment(data=filter(DF, type=="Actual grouped"),
               aes(x=ages, y=eta1, xend=ages.up, yend=eta1.up), size=1)+
  geom_line(data=filter(DF, type=="Fitted+Forecast"),
            aes(y=eta1),size=1)+
  geom_line(data=filter(DF, type=="Fitted+Forecast"),
            aes(y=eta1.up),size=1, linetype="dotted")+
  geom_line(data=filter(DF, type=="Fitted+Forecast"),
            aes(y=eta1.low),size=1, linetype="dotted")+
  facet_wrap(~years1, 2, 6, scales="free_y")+
  labs(x="age", y="log-mortality", title=cou.j)
print(p)


## over years
if(mg<=25){
  p <- ggplot(DF, aes(x=years1, y=eta1))+
    geom_line(data=filter(DF, type=="Fitted+Forecast"),
              aes(y=eta1, col=colo, group=ages),size=1)+
    geom_line(data=filter(DF, type=="Fitted+Forecast"),
              aes(y=eta1.up, col=colo, group=ages),size=1, linetype="dotted")+
    geom_line(data=filter(DF, type=="Fitted+Forecast"),
              aes(y=eta1.low, col=colo, group=ages),size=1, linetype="dotted")+
    geom_point(data=filter(DF, type=="Actual grouped"),
               aes(y=eta1, col=colo))+
    facet_wrap(~group, 5, 5, scales="free_y")+
    theme(legend.position = "none")+
    scale_color_viridis(discrete=FALSE, option="inferno")+
    labs(x="year", y="log-mortality", title=cou.j)+
    geom_rect(data=DF[1:mg,],
              aes(xmin = min(t2)-0.2, xmax = max(t2)+0.2,
                  ymin = -Inf, ymax = Inf),
              fill="cyan", alpha=0.05)
  p
}else{
  DFsel <- subset(DF, ages%in%seq(0,90,10))
  p <- ggplot(DFsel, aes(x=years1, y=eta1))+
    geom_line(data=filter(DFsel, type=="Fitted+Forecast"),
              aes(y=eta1, col=colo, group=ages),size=1)+
    geom_line(data=filter(DFsel, type=="Fitted+Forecast"),
             aes(y=eta1.up, col=colo, group=ages),size=1, linetype="dotted")+
    geom_line(data=filter(DFsel, type=="Fitted+Forecast"),
              aes(y=eta1.low, col=colo, group=ages),size=1, linetype="dotted")+
    geom_point(data=filter(DFsel, type=="Actual grouped"),
               aes(y=eta1, col=colo))+
    facet_wrap(~group, 2, 5, scales="free_y")+
    theme(legend.position = "none")+
    scale_color_viridis(discrete=FALSE, option="inferno")+
    labs(x="year", y="log-mortality", title=cou.j)+
    geom_rect(data=DFsel[1:mg,],
              aes(xmin = min(t2)-0.2, xmax = max(t2)+0.2,
                  ymin = -Inf, ymax = Inf),
              fill="cyan", alpha=0.05)
  p
}









## bootstrap deaths from fitted surface
Y1hat <- E1*exp(ETA1hat)
Yg1hat <- G%*%Y1hat
nb <- 100
## where to save all bootstrap+forecast log-mortality
ETAshat <- array(0, dim=c(m, n, nb))
b=1
for(b in 1:nb){
  Yg1.b <- matrix(rpois(mg*n1, c(Yg1)), mg, n1)
  FITi.b <- PSinfantGrouped(a.low = ag.low,
                            Yg=Yg1.b,
                            a=a,
                            E=E1,
                            WEIg=WEIg1,
                            lambdas=OPT$par,
                            verbose=FALSE,
                            kappa.shape=0,
                            infant=infantj)
  ## compute deltas : CIs of relative derivatives
  deltasi <- deltasFUN(FITi.b)
  ## forecasting
  CPSg.i.b <- CPSfunctionGrouped(a.low=ag.low, Yg1=Yg1.b, 
                                 WEIg1=WEIg1,
                                 a=a, E1=E1, obs.yrs=t1, 
                                 for.hor=max(t),
                                 ETA1hat=ETA1hat, deltas=deltasi, S=S, 
                                 lambdas = OPT$par, verbose = FALSE,
                                 infant = infantj,
                                 kappa.shape=0)
  ETAshat[,,b] <- CPSg.i.b$ETA
  cat(b, "\n")
}

## plotting, starting from the beginning
## actual log-mortality
DFg <- expand.grid(list(ages=ag.low, years1=t1))
DFg$type <- "Actual grouped"
DFg$eta1 <- c(ETAg1)
DFg$ages.up <- ag.up+1
DFg$eta1.up <- c(ETAg1)
## fitted
DFhat <- expand.grid(list(ages=a, years1=t1))
DFhat$type <- "Fitted"
DFhat$eta1 <- c(ETA1hat)
DFhat$ages.up <- NA
DFhat$eta1.up <- NA
## forecast
DFfor <- expand.grid(list(ages=a, years1=t))
DFfor$type <- "Fitted+Forecast"
DFfor$eta1 <- c(ETAhat)
DFfor$ages.up <- NA
DFfor$eta1.up <- NA
## all
DF <- rbind(DFg, DFhat, DFfor) 
DF$group <- factor(c(rep(ag.lab, n1),rep(rep(ag.lab, lg), n1),
                     rep(rep(ag.lab, lg), n)), levels = ag.lab)
DF$colo <- c(rep(ag.mid, n1), rep(a, n1), rep(a, n))
## bootstrap
DFboot <- expand.grid(list(ages=a, years1=t, boot=1:nb))
DFboot$eta <- c(ETAshat)

p <- ggplot(DFboot, aes(x=ages, y=eta))+
  geom_point(aes(y=eta), size=1, col="grey")+
  facet_wrap(~years1, 2, 6, scales="free_y")+
  labs(x="age", y="log-mortality", title=cou.j)
p

DFbootsel <- subset(DFboot, ages%in%seq(0,90,10))
p <- ggplot(DFbootsel, aes(x=years1, y=eta))+
  geom_point(aes(y=eta),size=1)+
  facet_wrap(~ages, 2, 5, scales="free_y")
p
  
  geom_line(data=filter(DF, type=="Fitted+Forecast"),
            aes(y=eta1, col=colo, group=ages),size=1)+
  geom_point(data=filter(DF, type=="Actual grouped"),
             aes(y=eta1, col=colo))+
  facet_wrap(~group, 5, 5, scales="free_y")+
  theme(legend.position = "none")+
  scale_color_viridis(discrete=FALSE, option="inferno")+
  labs(x="year", y="log-mortality", title=cou.j)+
  geom_rect(data=DF[1:mg,],
            aes(xmin = min(t2)-0.2, xmax = max(t2)+0.2,
                ymin = -Inf, ymax = Inf),
            fill="cyan", alpha=0.05)
p

## over age
p <- ggplot(DF, aes(x=ages, y=eta1, color=type)) +
  geom_segment(data=filter(DF, type=="Actual grouped"),
               aes(x=ages, y=eta1, xend=ages.up, yend=eta1.up), size=1)+
  geom_line(data=filter(DF, type=="Fitted+Forecast"),
            aes(y=eta1),size=1)+
  geom_line(data=filter(DF, type=="Fitted"),
            aes(y=eta1),size=1, linetype="dotted")+
  facet_wrap(~years1, 2, 6, scales="free_y")+
  labs(x="age", y="log-mortality", title=cou.j)
print(p)
## over years
if(mg<=25){
  p <- ggplot(DF, aes(x=years1, y=eta1))+
    geom_line(data=filter(DF, type=="Fitted+Forecast"),
              aes(y=eta1, col=colo, group=ages),size=1)+
    geom_point(data=filter(DF, type=="Actual grouped"),
               aes(y=eta1, col=colo))+
    facet_wrap(~group, 5, 5, scales="free_y")+
    theme(legend.position = "none")+
    scale_color_viridis(discrete=FALSE, option="inferno")+
    labs(x="year", y="log-mortality", title=cou.j)+
    geom_rect(data=DF[1:mg,],
              aes(xmin = min(t2)-0.2, xmax = max(t2)+0.2,
                  ymin = -Inf, ymax = Inf),
              fill="cyan", alpha=0.05)
  p
}else{
  DFsel <- subset(DF, ages%in%seq(0,90,10))
  p <- ggplot(DFsel, aes(x=years1, y=eta1))+
    geom_line(data=filter(DFsel, type=="Fitted+Forecast"),
              aes(y=eta1, col=colo, group=ages),size=1)+
    geom_point(data=filter(DFsel),
               aes(y=eta1, col=colo))+
    facet_wrap(~group, 2, 5, scales="free_y")+
    theme(legend.position = "none")+
    scale_color_viridis(discrete=FALSE, option="inferno")+
    labs(x="year", y="log-mortality", title=cou.j)+
    geom_rect(data=DFsel[1:mg,],
              aes(xmin = min(t2)-0.2, xmax = max(t2)+0.2,
                  ymin = -Inf, ymax = Inf),
              fill="cyan", alpha=0.05)
  p
}


## estimating perturbation 

## 2020 and 2021 are treated independently 

## y_2020 ~ Poi(C mu)
## log(mu) = eta^for_2020 + c + \delta
## smooth \delta and s.t. 1' \delta = 0

## let's start from 2020



## extract 2020 forecast log-mortality (~offset)
eta20 <- ETAhat[,ncol(ETAhat)-1]

## extract deaths and exposures for 2020
Y.j20 <- subset(deaths.j, Year==2020)
yg20 <- Y.j20$Deaths
e20 <- subset(offset, Country==cou.j & Year==2020)$Population
ag.low <- unique(Y.j20$Age)
ag.up1  <- ag.low-1
if(ag.low[2]==1){
  ag.up <- c(0, ag.up1[ag.up1>0], max(a))
}else{
  ag.up <- c(ag.up1[ag.up1>0], max(a))
}
cbind(ag.low, ag.up)
ag.mid <- (ag.low+ag.up)/2 
ag.lab <- paste(ag.low, ag.up,sep="-")

lg <- ag.up-ag.low+1
mg <- length(ag.low)
G20 <- matrix(0, mg, m)
rownames(G20) <- ag.mid
colnames(G20) <- a
for(i in 1:mg){
  ag.low.i <- ag.low[i]
  ag.up.i <- ag.up[i]
  wc <- which(a>=ag.low.i & a<=ag.up.i)
  G20[i,wc] <- 1
}
all(colSums(G20)==1)


## only for plotting
eg20 <- G20 %*% e20 
etag20 <- log(yg20/eg20)
## create C matrix with e20
C20 <- G20
C20[G20==1] <- c(e20)

## penalty stuff
D <- diff(diag(m), diff=2)
tDD <- t(D)%*%D
lambda <- 10^4
P0 <- lambda*tDD
P <- matrix(0,m+1,m+1)
P[-1,-1] <- P0

## design matrix
U <- cbind(1, diag(m))
kappa <- 0
## constraining \delta to sum up to 0
H <- matrix(c(0, rep(1,m)), 1, m+1)
## starting eta (log-mortality minus forecast log-mortality in 2020)
eta <- rep(0.01,m)
## PCLM regression with eta20 as offset
max.it <- 100
for(it in 1:max.it){
  gamma   <- exp(eta + eta20)
  mu      <- c(C20 %*% gamma)
  X       <- (C20 * ((1 / mu) %*% t(gamma)) ) %*% U
  w       <- as.vector(mu)
  r       <- yg20 - mu + C20 %*% (gamma * eta)
  tXWX    <- t(X) %*% (w * X) 
  tXWXpP  <- tXWX + P
  tXr     <- t(X) %*% r
  ## adding constraints
  LHS     <- rbind(cbind(tXWXpP, t(H)),
               cbind(H, 0))
  RHS     <- matrix(c(tXr, kappa), ncol=1)
  coeff   <- solve(LHS, RHS)
  ##
  betas   <- coeff[1:(m+1)]
  eta.old <- eta
  eta     <- U%*%betas
  dif.eta <- max(abs((eta - eta.old)/eta.old) )
  if(dif.eta < 1e-04 & it > 4) break
  cat(it, dif.eta, "\n")
}
c.hat <- coeff[1]
delta.hat <- coeff[1:m+1]
sum(delta.hat)
omega <- coeff[m+2]  
plot(a, delta.hat)
eta.hat <- eta20 + c.hat + delta.hat # U%*%betas+eta20
plot(a, eta.hat)

## observed
DFg <- data.frame(ages=ag.low, 
                  type="Actual grouped",
                  eta1=etag20,
                  ages.up=ag.up+1,
                  eta1.up=etag20)
## forecast
DFfor <- data.frame(ages=a, 
                    type="Forecast",
                    eta1=eta20,
                    ages.up=NA,
                    eta1.up=NA)
## fitted
DFhat <- data.frame(ages=a, 
                    type="Fitted",
                    eta1=eta.hat,
                    ages.up=NA,
                    eta1.up=NA)
## combine
DF <- rbind(DFg, DFfor, DFhat)

p <- ggplot(DF, aes(x=ages, y=eta1, color=type)) +
  geom_segment(data=filter(DF, type=="Actual grouped"),
               aes(x=ages, y=eta1, xend=ages.up, yend=eta1.up), size=1)+
  geom_line(data=filter(DF, type=="Fitted"),
            aes(y=eta1), size=1)+
  geom_line(data=filter(DF, type=="Forecast"),
            aes(y=eta1),size=1, linetype="dotted")+
  labs(x="age", y="log-mortality", title=cou.j)
p

DFdelta <- data.frame(ages=a, delta=delta.hat)
DFc <- data.frame(x=-2, c=c.hat)
p <- ggplot(DFdelta, aes(x=ages, y=delta))+
  geom_line(size=2, colour="darkgreen")+
  geom_hline(yintercept=0, linetype="dashed", color = "darkred")+
  labs(x="age", y=expression(paste("c, ", delta)), title=cou.j)+
  geom_point(aes(x=DFc$x,y=DFc$c),
             colour="darkblue", shape=3, size=3, stroke = 2)

p


## END




