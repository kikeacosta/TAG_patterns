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
deaths <- read.csv("data_inter/deaths_sourced_infant_based.csv", header=TRUE)
## loading population
offset <- read.csv("data_inter/offsets.csv", header=TRUE)
## sex
sex <- "f"
sexLT <- "F"
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

## questions:
## FEMALES

## Currently impossible to fit+forecast

## too few data
##  3) Andorra is missing age 30-35 in 2015
## 30) French Polynesia: 3 available years are too few (unless I modified the whole forecasting model)
## 54) Moldova: 4 available years are too few (unless I modified the whole forecasting model)
## 64) Peru: 3 available years are too few (unless I modified the whole forecasting model)
## 67) Quatar: 4 available years and old ones (2010 2011 2012 2014) are too few (unless I modified the whole forecasting model)

## 2019 not available (likely solvable if 2019 is interpolated as previous not-available years)
## 19) Cuba: missing 2019, what to do?
## 88) Uruguay: missing 2019, what to do?

## Issues with infant mortality

## Bermuda: no infant deaths for all pre-pandemic years but 2014 (~impossible to have a infant-soecialized coeff)
## Faroe Islands: no infant deaths in last years
## Liechtenstein: absence of infant deaths in the last period causes some issue (solvable?) 
## Montserrat: no infant deaths for all pre-pandemic years but 2015 (~impossible to have a infant-specialized coeff)
## Tuvalu: no infant deaths in the first available years (~impossible to have a infant-specialized coeff)


## Evident strange age-patterns (it could be solvable, but it denotes issues in actual data)

## Albania: strange age-pattern at highest ages and increasing infant mortality
## Aruba: strange increasing of infant mortality. Actually death at age 0 in 2012 and 2014
## Azerbaijan: odd rapid increase of oldest-age mortality
## Bosnia and Herzegovina: decreasing observed mortality at higher ages (I can force monotonicity over age, but weird data)
## Bulgaria: wiggling age-pattern in first years at old ages (force monotonicity over age?)
## Croatia: wiggling age-pattern in first years at old ages (force monotonicity over age?)
## Lesotho: special treatment is needed (would enforcement of monotonicity be enough?)
## Lithuania: wiggling age-pattern at old ages (force monotonicity over age?)
## Maldives: too strong leveling-off at odest ages and odd data in 2019 (missing 2018)
## Mauritius: very weird 2011 data
## Montenegro: decreasing observed mortality at higher ages (I can force monotonicity over age, but weird data)
## New Zeland: a bit too strong leveling-off at oldest ages in last observed years
## North Macedonia: wiggling age-pattern at old ages (force monotonicity over age?)
## Romania: odd rapid increase of oldest-age mortality
## Serbia: a bit too strong leveling-off at oldest ages in all observed years
## South Korea: wiggling/decreasing age-pattern at old ages (force monotonicity over age?)
## Taiwan: optimal smoothing parameter produces wiggling age-pattern, likely need to impose extra-smoothness
## Turks and Caicos Islands: 3 available years, the fit is possible because they cover 5 years and include 2019. However no infant deaths for all pre-pandemic years but 2015 (~impossible to have a infant-specialized coeff)
## Ukraine: optimal smoothing parameter produces wiggling age-pattern, likely need to impose extra-smoothness

## Issues only in computing e0

## due to large open-age group
## Belize: extremely large open-age interval, 65+ => issues in e0 computation
## Iran: extremely large open-age interval, 65+ => issues in e0 computation
## Thailand: extremely large open-age interval, 65+ => issues in e0 computation

## general issue (due to GC lifetable code?) also because the fit on log-rates seems OK
## Colombia: good fit on rates, but likely issues in computing e0
## Costa Rica: good fit on rates, but likely issues in computing e0
## Ecuador: good fit on rates, but likely issues in computing e0
## Honk-Kong: good fit on rates, but likely issues in computing e0
## Japan: good fit on rates, but likely issues in computing e0
## Malasya: good fit on rates, but likely issues in computing e0
## Oman: good fit on rates, but likely issues in computing e0
## Panama: good fit on rates, but likely issues in computing e0
## State of Palestine: good fit on rates, but likely issues in computing e0
## Suriname: good fit on rates, but likely issues in computing e0
## Uzbekistan: good fit on rates, but likely issues in computing e0


## General questions from GC:

## Iceland: data available only from 2013?
## Italy: only by age-group? and from 2011?













nocou <- c(3,30,54,64,67,19,88)
#unique(deaths$Country[which(deaths$age_spn==5)])
yescou <- 1:nc
yescou <- yescou[-nocou]

j=1
#pdf("ForecastFemales.pdf", width = 12, height = 10)
for(j in 1:nc){
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
}


# ## plotting grouped actual log-mortality over ages
# DFg <- expand.grid(list(ages=ag.low, years1=t1))
# DFg$type <- "Actual grouped"
# DFg$eta1 <- c(ETAg1)
# DFg$ages.up <- ag.up+1
# DFg$eta1.up <- c(ETAg1)
# DFg$ag.lab <- factor(c(rep(ag.lab, n1)), levels = ag.lab)
# ## over age-groups
# p <- ggplot(DFg, aes(x=ages, y=eta1, color=type)) +
#   geom_segment(data=filter(DFg, type=="Actual grouped"),
#                aes(x=ages, y=eta1, xend=ages.up, yend=eta1.up))+
#   facet_wrap(~years1, 2, 7, scales="free_y")+
#   labs(x="age", y="log-mortality", title=cou.j)
# p
# ## over years
# if(mg<=25){
#   p <- ggplot(DFg, aes(x=years1, y=eta1)) +
#     geom_point(data=filter(DFg, type=="Actual grouped"),
#                aes(y=eta1, col=ag.lab))+
#     facet_wrap(~ag.lab, 5, 5, scales="free_y")+
#     theme(legend.position = "none") +
#     scale_color_viridis(discrete=TRUE, option="inferno")+
#     labs(x="year", y="log-mortality", title=cou.j)
#   p
# }else{
#   DFgsel <- subset(DFg, ages%in%seq(0,90,10))
#   p <- ggplot(DFgsel, aes(x=years1, y=eta1)) +
#     geom_point(data=filter(DFgsel, type=="Actual grouped"),
#                aes(y=eta1, col=ag.lab))+
#     facet_wrap(~ag.lab, 2, 5, scales="free_y")+
#     theme(legend.position = "none") +
#     scale_color_viridis(discrete=TRUE, option="inferno")+
#     labs(x="year", y="log-mortality", title=cou.j)
#   p
# }

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

# ## plotting, starting from the beginning
# ## actual log-mortality
# DFg <- expand.grid(list(ages=ag.low, years1=t1))
# DFg$type <- "Actual grouped"
# DFg$eta1 <- c(ETAg1)
# DFg$ages.up <- ag.up+1
# DFg$eta1.up <- c(ETAg1)
# ## fitted
# DFhat <- expand.grid(list(ages=a, years1=t1))
# DFhat$type <- "Fitted"
# DFhat$eta1 <- c(ETA1hat)
# DFhat$ages.up <- NA
# DFhat$eta1.up <- NA
# ## all
# DF <- rbind(DFg, DFhat)
# DF$group <- factor(c(rep(ag.lab, n1), rep(rep(ag.lab, lg), n1)),
#                    levels = ag.lab)
# DF$colo <- c(rep(ag.mid, n1), rep(a, n1))
# ## over age
# p <- ggplot(DF, aes(x=ages, y=eta1, color=type)) +
#   geom_segment(data=filter(DF, type=="Actual grouped"),
#                aes(x=ages, y=eta1, xend=ages.up, yend=eta1.up))+
#   geom_line(data=filter(DF, type=="Fitted"),
#             aes(y=eta1),size=1)+
#   facet_wrap(~years1, 2, 5)+
#   labs(x="age", y="log-mortality", title=cou.j)
# p
# ## over years
# if(mg<=25){
#   p <- ggplot(DF, aes(x=years1, y=eta1)) +
#     geom_line(data=filter(DF, type=="Fitted"),
#               aes(y=eta1, col=colo, group=ages),size=1)+
#     geom_point(data=filter(DF, type=="Actual grouped"),
#                aes(y=eta1, col=colo))+
#     facet_wrap(~group, 5, 5, scales="free_y")+
#     theme(legend.position = "none") +
#     scale_color_viridis(discrete=FALSE, option="inferno")+
#     labs(x="year", y="log-mortality", title=cou.j)
#   p
# }else{
#   DFsel <- subset(DF, ages%in%seq(0,90,10))
#   p <- ggplot(DFsel, aes(x=years1, y=eta1)) +
#     geom_line(data=filter(DFsel, type=="Fitted"),
#               aes(y=eta1, col=colo, group=ages),size=1)+
#     geom_point(data=filter(DFsel, type=="Actual grouped"),
#                aes(y=eta1, col=colo))+
#     facet_wrap(~group, 2, 5, scales="free_y")+
#     theme(legend.position = "none") +
#     scale_color_viridis(discrete=FALSE, option="inferno")+
#     labs(x="year", y="log-mortality", title=cou.j)
#   p
# }
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
DFfor$type <- "Fitted+Forecast"
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
  geom_line(data=filter(DF, type=="Fitted+Forecast"),
            aes(y=eta1),size=1)+
  geom_line(data=filter(DF, type=="Fitted"),
            aes(y=eta1),size=1, linetype="dotted")+
  facet_wrap(~years1, 2, 6, scales="free_y")+
  labs(x="age", y="log-mortality", title=cou.j)
# print(p)
# ## over years
# if(mg<=25){
#   p <- ggplot(DF, aes(x=years1, y=eta1))+
#     geom_line(data=filter(DF, type=="Fitted+Forecast"),
#               aes(y=eta1, col=colo, group=ages),size=1)+
#     geom_point(data=filter(DF, type=="Actual grouped"),
#                aes(y=eta1, col=colo))+
#     facet_wrap(~group, 5, 5, scales="free_y")+
#     theme(legend.position = "none")+
#     scale_color_viridis(discrete=FALSE, option="inferno")+
#     labs(x="year", y="log-mortality", title=cou.j)+
#     geom_rect(data=DF[1:mg,],
#               aes(xmin = min(t2)-0.2, xmax = max(t2)+0.2,
#                   ymin = -Inf, ymax = Inf),
#               fill="cyan", alpha=0.05)
#   p
# }else{
#   DFsel <- subset(DF, ages%in%seq(0,90,10))
#   p <- ggplot(DFsel, aes(x=years1, y=eta1))+
#     geom_line(data=filter(DFsel, type=="Fitted+Forecast"),
#               aes(y=eta1, col=colo, group=ages),size=1)+
#     geom_point(data=filter(DFsel),
#                aes(y=eta1, col=colo))+
#     facet_wrap(~group, 2, 5, scales="free_y")+
#     theme(legend.position = "none")+
#     scale_color_viridis(discrete=FALSE, option="inferno")+
#     labs(x="year", y="log-mortality", title=cou.j)+
#     geom_rect(data=DFsel[1:mg,],
#               aes(xmin = min(t2)-0.2, xmax = max(t2)+0.2,
#                   ymin = -Inf, ymax = Inf),
#               fill="cyan", alpha=0.05)
#   p
# }
#   
# ## testing fit and forecast by life expectancy
# source("~/WORK/JUNK/LifeTableMatrixAlgebra/LifeTableFUN.R")
# e0hat <- numeric(n)
# e0act <- numeric(n1)
# i=4
# for(i in 1:n){
#   LT.i <- lifetable.mx(x=a, mx=exp(ETAhat)[,i], sex=sexLT)
#   e0hat[i] <- LT.i$ex[1]
#   if(i<=n1& whi.ava[i]){
#     agx.i <- G%*%LT.i$ax
#     LTact.i <- lifetable.mx(x=ag.low, mx=exp(ETAg1)[,i], 
#                             sex=sexLT)
#     e0act[i] <- LTact.i$ex[1]
#   }
# }
# 
# DFe0act <- data.frame(years1=t1, e0=e0act, type="Actual")
# DFe0act <- DFe0act[whi.ava,]
# DFe0hat <- data.frame(years1=t, e0=e0hat, type="Fitted+Forecast")
# DFe0 <- rbind(DFe0act, DFe0hat)
# 
# p <- ggplot(DFe0, aes(x=years1, y=e0, col=type))+
#   geom_line(data=filter(DFe0, type=="Fitted+Forecast"), 
#             aes(y=e0, col=type), size=2)+
#   geom_point(data=filter(DFe0, type=="Fitted+Forecast"), 
#             aes(y=e0, col=type), size=3.5)+
#   geom_point(data=filter(DFe0, type=="Actual"), 
#             aes(y=e0, col=type), size=3.5)+
#   theme(legend.position = "none")+
#   labs(x="year", y="life expectancy at birth", title=cou.j)+
#   geom_rect(data=DFe0[1,], 
#             aes(xmin = min(t2)-0.2, xmax = max(t2)+0.2, 
#                 ymin = -Inf, ymax = Inf),
#             fill="pink", color=NA, alpha=0.2)+
#   scale_x_continuous(breaks=seq(2010, 2021, 1), minor=FALSE)+
#   scale_y_continuous(breaks=seq(70, 90, 0.2), minor=FALSE)
# print(p)

}
#dev.off()












## END




