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

## loading deaths
deaths <- read.csv("data_inter/annual_deaths_countries_selected_sources.csv", header=TRUE)
## loading population
offset <- read.csv("data_inter/offsets.csv", header=TRUE)
## sex
sex <- "f"
deaths <- subset(deaths, Sex==sex)
offset <- subset(offset, Sex==sex)

## ages and dimensions
a <- 0:100
m <- length(a)
## pre-pandemic years
t1 <- 2015:2019
n1 <- length(t1)
## pandemic years
t2 <- 2020:2021
n2 <- length(t2)
## all years
t <- c(t1,t2)
n <- n1+n2
## countries
cou <- unique(deaths$Country)
nc <- length(cou)

## questions:
##
## American Samoa: age 10 is missing in 2016 and year 2018 and off low value in the open-age groups
## Andorra: really weird observed age-pattern
## Kazakhstan sure upper age-group is 100
## Tunisia: how come we have zero deaths in the last age-group 85-100






j=3

#for(j in 1:2){
couj <- cou[j]
## working on the deaths and age-groups
Yj <- subset(deaths, Country==couj & Year%in%t1)
ag.low <- unique(Yj$Age)
ag.up1  <- ag.low-1
if(ag.low[2]==1){
  ag.up <- c(0, ag.up1[ag.up1>0], max(a))
}else{
  ag.up <- c(ag.up1[ag.up1>0], max(a))
}
lg <- ag.up-ag.low+1
mg <- length(ag.low)
## only for plotting/aesthetic
ag.mid <- (ag.low+ag.up)/2 
ag.lab <- paste(ag.low, ag.up,sep="-")
Yg1 <- matrix(Yj$Deaths, mg, n1)
## working on the offset
Ej <- subset(offset, Country==couj & Year%in%t1)
E1 <- matrix(Ej$Population, m, n1)

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
## observed log-mortality, by age-groups
ETAg1 <- log(Yg1/Eg1)

## plotting grouped actual log-mortality over ages
DFg <- expand.grid(list(ages=ag.low, years1=t1))
DFg$type <- "Actual grouped"
DFg$eta1 <- c(ETAg1)
DFg$ages.up <- ag.up+1
DFg$eta1.up <- c(ETAg1)
DFg$ag.lab <- factor(c(rep(ag.lab, n1)), levels = ag.lab)
## over age-groups
p <- ggplot(DFg, aes(x=ages, y=eta1, color=type)) +
  geom_segment(data=filter(DFg, type=="Actual grouped"),
               aes(x=ages, y=eta1, xend=ages.up, yend=eta1.up))+
  facet_wrap(~years1, 2, 7, scales="free_y")+
  labs(x="age", y="log-mortality", title=couj)
p
#locator(1)
}
## over years
p <- ggplot(DFg, aes(x=years1, y=eta1)) +
  geom_point(data=filter(DFg, type=="Actual grouped"),
             aes(y=eta1, col=ag.lab))+
  facet_wrap(~ag.lab, 5, 5, scales="free_y")+
  theme(legend.position = "none") +
  scale_color_viridis(discrete=TRUE, option="inferno")+
  labs(x="year", y="log-mortality")
p

source("~/WORK/ConstrainedMortalityForecast/GroupedData/FUNCTIONS/GroupedSmoothConstrainedMortalityForecasting_Functions.R")

## weights for grouped-exposures
WEIg1 <- matrix(1, mg,n1)
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
## all
DF <- rbind(DFg, DFhat) 
DF$group <- factor(c(rep(ag.lab, n1), rep(rep(ag.lab, lg), n1)),
                   levels = ag.lab)
DF$colo <- c(rep(ag.mid, n1), rep(a, n1))
## over age
p <- ggplot(DF, aes(x=ages, y=eta1, color=type)) +
  geom_segment(data=filter(DF, type=="Actual grouped"),
               aes(x=ages, y=eta1, xend=ages.up, yend=eta1.up), size=1)+
  geom_line(data=filter(DF, type=="Fitted"),
            aes(y=eta1),size=1)+
  facet_wrap(~years1, 2, 7)+
  labs(x="age", y="log-mortality")
p
## over years 
p <- ggplot(DF, aes(x=years1, y=eta1)) +
  geom_line(data=filter(DF, type=="Fitted"),
            aes(y=eta1, col=colo, group=ages),size=1)+
  geom_point(data=filter(DF, type=="Actual grouped"),
             aes(y=eta1, col=colo))+
  facet_wrap(~group, 5, 5, scales="free_y")+
  theme(legend.position = "none") +
  scale_color_viridis(discrete=FALSE, option="inferno")+
  labs(x="year", y="log-mortality")
p

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
DFfor$type <- "Forecast"
DFfor$eta1 <- c(ETAhat)
DFfor$ages.up <- NA
DFfor$eta1.up <- NA
## all
DF <- rbind(DFg, DFhat, DFfor) 
DF$group <- factor(c(rep(ag.lab, n1),rep(rep(ag.lab, lg), n1),
                     rep(rep(ag.lab, lg), n)), levels = ag.lab)
DF$colo <- c(rep(ag.mid, n1), rep(a, n1), rep(a, n))
## over age
p <- ggplot(DF, aes(x=ages, y=eta1, color=type)) +
  geom_segment(data=filter(DF, type=="Actual grouped"),
               aes(x=ages, y=eta1, xend=ages.up, yend=eta1.up), size=1)+
  geom_line(data=filter(DF, type=="Forecast"),
            aes(y=eta1),size=1)+
  geom_line(data=filter(DF, type=="Fitted"),
            aes(y=eta1),size=1, linetype="dotted")+
  facet_wrap(~years1, 5, 5, scales="free_y")+
  labs(x="age", y="log-mortality")
p
## over years 
p <- ggplot(DF, aes(x=years1, y=eta1))+
  geom_line(data=filter(DF, type=="Forecast"),
            aes(y=eta1, col=colo, group=ages),size=1)+
  geom_point(data=filter(DF, type=="Actual grouped"),
             aes(y=eta1, col=colo))+
  facet_wrap(~group, 5, 5, scales="free_y")+
  theme(legend.position = "none")+
  scale_color_viridis(discrete=FALSE, option="inferno")+
  labs(x="year", y="log-mortality")+
  geom_rect(data=DF[1:mg,], 
            aes(xmin = min(t2), xmax = max(t2), 
                ymin = -Inf, ymax = Inf),
            fill="cyan", alpha=0.05)
p


















## END




