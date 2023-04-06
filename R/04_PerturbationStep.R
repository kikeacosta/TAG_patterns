library(MortalitySmooth)
library(magic)
library(colorspace)
library(ggplot2)
library(plotly)
library(viridis)
library(Matrix)
library(tidyverse)
source("R/04_PerturbationFunction.R")

deaths <- read.csv("data_inter/deaths_sourced_infant_based_99.csv", header=TRUE)
## loading population
offset <- read.csv("data_inter/offsets_99.csv", header=TRUE)

## sex
sex <- "f"
sexLT <- "F"
deaths <- subset(deaths, Sex==sex)
offset <- subset(offset, Sex==sex)

## ages and dimensions
a <- 0:99
m <- length(a)
## pandemic years for the first forecasting part
t2 <- 2020:2021
n2 <- length(t2)
## countries
cou <- unique(deaths$Country)
nc <- length(cou)

j=21
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


ETAhat <- CPSg.i$ETA
str(CPSg.i)
## analytic confidence intervals
B <- kronecker(CPSg.i$Bt, CPSg.i$Ba)
Veta.for <- B %*% CPSg.i$Valpha %*% t(B)
SEeta.for <- matrix(sqrt(diag(Veta.for)), m, n)

ETAup <- ETAhat + 2*SEeta.for
ETAlow <- ETAhat - 2*SEeta.for


## plotting, starting from the beginning
## actual log-mortality
DFg <- expand.grid(list(ages=ag.low, years1=t1))
DFg$type <- "Actual grouped"
DFg$eta1 <- c(ETAg1)
DFg$ages.up <- ag.up+1
DFg$eta1.up <- c(ETAg1)
DFg$eta1.low <- NA
## forecast
DFfor <- expand.grid(list(ages=a, years1=t))
DFfor$type <- "Fitted+Forecast"
DFfor$eta1 <- c(ETAhat)
DFfor$ages.up <- NA
DFfor$eta1.up <- c(ETAup)
DFfor$eta1.low <- c(ETAlow)

## all
DF <- rbind(DFg, DFfor) 
DF$group <- factor(c(rep(ag.lab, n1),
                     rep(rep(ag.lab, lg), n)), levels = ag.lab)
DF$colo <- c(rep(ag.mid, n1), rep(a, n))

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
  labs(x="age", y="log-mortality", title=paste(cou.j, sex))
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
    labs(x="year", y="log-mortality", title=paste(cou.j, sex))+
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
    labs(x="year", y="log-mortality", title=paste(cou.j, sex))+
    geom_rect(data=DFsel[1:mg,],
              aes(xmin = min(t2)-0.2, xmax = max(t2)+0.2,
                  ymin = -Inf, ymax = Inf),
              fill="cyan", alpha=0.05)
  p
}



## y_2020 ~ Poi(C mu)
## log(mu) = eta^for_2020 + c + \delta
## smooth \delta and s.t. 1' \delta = 0

## which pandemic years are available
tF.ava <- t.ava.dea[!t.ava.dea%in%t1.ava]
n.boot <- 500
j=1
list.out <- list()
for(j in 1:length(tF.ava)){
  
  ## for a given pandemic year
  
  ## extract deaths and exposures for tF.ava[j] 
  obs.y0 <- subset(DFTR, type=="Obs Deaths" & years==tF.ava[j])
  obs.y <- obs.y0$value
  obs.e0 <- subset(DFTR, type=="Obs Offset" & years==tF.ava[j])
  obs.e <- obs.e0$value
  
  ## extract tF.ava[j] forecast log-mortality (~offset in this analysis)
  eta.for0 <- subset(DFTR, type=="Baseline Logrates" & years==tF.ava[j])
  eta.for <- eta.for0$value 
  
  ## age-structure
  ag.low <- obs.y0$ages
  
  ## optimizing lambda for delta
  QIC_Pert <- function(par){
    FIT <- PertubFUN(ag.low=ag.low,
                     obs.y=obs.y,
                     obs.e=obs.e,
                     a=a,
                     eta.for=eta.for,
                     lambda=par)
    FIT$QIC
  }
  ## optimizing lambdas using greedy grid search
  OPT <- cleversearch(QIC_Pert, lower=-1, upper=8,
                      ngrid=19, logscale=TRUE, startvalue=10^4,
                      verbose=FALSE)
  lambda <- OPT$par
  
  ## estimating perturbation with optimal lambdas
  PertOPT <- PertubFUN(ag.low=ag.low,
                       obs.y=obs.y,
                       obs.e=obs.e,
                       a=a,
                       eta.for=eta.for,
                       lambda=lambda)
  nb <- ncol(PertOPT$Ba)
  ## fitted values
  c.hat <- PertOPT$coeff[1]
  expc.hat <- exp(c.hat)
  delta.hat <- PertOPT$Ba%*%PertOPT$coeff[1:nb+1]
  expdelta.hat <- exp(delta.hat)
  pert.hat <- c.hat + delta.hat
  exppert.hat <- exp(pert.hat)
  eta.hat <- eta.for + pert.hat
  
  ## given the optimized lambda bootstrapping
  cs.hat <- numeric(n.boot)
  DELTAs.hat <- PERTs.hat <- ETAs.hat <- matrix(0, m, n.boot)
  for(b in 1:n.boot){
    ## sample new observation
    obs.y.b <- rpois(n=length(obs.y), lambda = obs.y)
    ## simulate a new forecast log-mortality
    eta.for.b <- rnorm(m, mean=eta.for, sd=SEeta.for[, which(t==tF.ava[j])])
    ## permutation step with optimized lambda and bootstrapped data
    PertOPT.b <- PertubFUN(ag.low=ag.low,
                           obs.y=obs.y.b,
                           obs.e=obs.e,
                           a=a,
                           eta.for=eta.for.b,
                           lambda=lambda)
    nb <- ncol(PertOPT.b$Ba)
    ## save useful objects
    c.hat.b <- PertOPT.b$coeff[1]
    delta.hat.b <- PertOPT.b$Ba %*% PertOPT.b$coeff[1:nb+1]
    pert.hat.b <- c.hat.b + delta.hat.b
    eta.hat.b <- eta.for.b + pert.hat.b
    cs.hat[b] <- c.hat.b
    DELTAs.hat[,b] <- delta.hat.b
    PERTs.hat[,b] <- pert.hat.b
    ETAs.hat[,b] <- eta.hat.b
  }
  ## 95% CIs (now including both sources of uncertainty)
  c.up <- quantile(cs.hat, probs = 0.975)
  c.low <- quantile(cs.hat, probs = 0.025)
  delta.up <- apply(DELTAs.hat, 1, quantile, probs = 0.975)
  delta.low <- apply(DELTAs.hat, 1, quantile, probs = 0.025)
  expc.up <- quantile(exp(cs.hat), probs = 0.975)
  expc.low <- quantile(exp(cs.hat), probs = 0.025)
  expdelta.up <- apply(exp(DELTAs.hat), 1, quantile, probs = 0.975)
  expdelta.low <- apply(exp(DELTAs.hat), 1, quantile, probs = 0.025)
  pert.up <- apply(PERTs.hat, 1, quantile, probs = 0.975)
  pert.low <- apply(PERTs.hat, 1, quantile, probs = 0.025)
  exppert.up <- apply(exp(PERTs.hat), 1, quantile, probs = 0.975)
  exppert.low <- apply(exp(PERTs.hat), 1, quantile, probs = 0.025)
  eta.up <- apply(ETAs.hat, 1, quantile, probs = 0.975)
  eta.low <- apply(ETAs.hat, 1, quantile, probs = 0.025)
  ## for fitted values I take the median of the bootstrapped values?
  ## summation constrain for delta doesn't hold!
  #c.hat <- quantile(cs.hat, probs = 0.5)
  #delta.hat <- apply(DELTAs.hat, 1, quantile, probs = 0.5)
  #expc.hat <- quantile(exp(cs.hat), probs = 0.5)
  #expdelta.hat <- apply(exp(DELTAs.hat), 1, quantile, probs = 0.5)
  #pert.hat <- apply(PERTs.hat, 1, quantile, probs = 0.5)
  #exppert.hat <- apply(exp(PERTs.hat), 1, quantile, probs = 0.5)
  #eta.hat <- apply(ETAs.hat, 1, quantile, probs = 0.5)
  
  
  DFTRdelta <- expand.grid(ages=a, years=tF.ava[j],
                           type=c("Delta", "Exp Delta"))
  DFTRdelta$value <- c(delta.hat, expdelta.hat)
  DFTRdelta$up <- c(delta.up, expdelta.up)
  DFTRdelta$low <- c(delta.low, expdelta.low)
  
  DFTRc <- expand.grid(ages=NA, years=tF.ava[j],
                       type=c("c", "Exp c"))
  DFTRc$value <- c(c.hat, expc.hat)
  DFTRc$up <- c(c.up, expc.up)
  DFTRc$low <- c(c.low, expc.low)
  
  DFTRpert <- expand.grid(ages=a, years=tF.ava[j],
                           type=c("Perturbation", "Exp Perturbation"))
  DFTRpert$value <- c(pert.hat, exppert.hat)
  DFTRpert$up <- c(pert.up, exppert.up)
  DFTRpert$low <- c(pert.low, exppert.low)
  
  DFTReta <- expand.grid(ages=a, years=tF.ava[j],
                         type=c("Fitted Logrates"))
  DFTReta$value <- c(eta.hat)
  DFTReta$up <- c(eta.up)
  DFTReta$low <- c(eta.low)
  
  list.out[[j]] <- rbind(DFTRdelta, DFTRc, DFTRpert, DFTReta)
  
  ## observed for plotting 
  eta.obs <- subset(DFTR, type=="Obs Logrates" & years==tF.ava[j])$value
  
  DFg <- data.frame(ages=ag.low, 
                    type="Actual grouped",
                    eta1=eta.obs,
                    ages.up=ag.up+1,
                    eta1.up=eta.obs,
                    eta1.low=NA)
  ## forecast
  DFfor <- data.frame(ages=a, 
                      type="Forecast",
                      eta1=eta.for,
                      ages.up=NA,
                      eta1.up=NA,
                      eta1.low=NA)
  ## fitted
  DFhat <- data.frame(ages=a, 
                      type="Fitted",
                      eta1=eta.hat,
                      ages.up=NA,
                      eta1.up=eta.up,
                      eta1.low=eta.low)
  ## combine
  DF <- rbind(DFg, DFfor, DFhat)
  
  p <- ggplot(DF, aes(x=ages, y=eta1, color=type)) +
    geom_segment(data=filter(DF, type=="Actual grouped"),
                 aes(x=ages, y=eta1, xend=ages.up, yend=eta1.up), size=1)+
    geom_line(data=filter(DF, type=="Fitted"),
              aes(y=eta1), size=1)+
    geom_ribbon(data=filter(DF, type=="Fitted"),
                aes(ymin=eta1.low, ymax=eta1.up), alpha=.2)+
    geom_line(data=filter(DF, type=="Forecast"),
              aes(y=eta1),size=1.2)+
    labs(x="age", y="log-mortality", title=paste(cou.j, sex))
  p
  
  DFexpdelta <- data.frame(ages=a, expdelta=expdelta.hat,
                           expdelta.low=expdelta.low,
                           expdelta.up=expdelta.up)
  DFexpc <- data.frame(x=-2, expc=expc.hat, expc.low=expc.low, expc.up=expc.up)
  
  p <- ggplot(DFexpdelta, aes(x=ages))+
    geom_line(aes(y=expdelta), size=2, colour="darkgreen")+
    geom_ribbon(aes(ymin=expdelta.low, ymax=expdelta.up), alpha=0.5)+
    geom_pointrange(aes(x=DFexpc$x, y=DFexpc$expc, 
                        ymin = DFexpc$expc.low, 
                        ymax = DFexpc$expc.up),
                    colour="darkblue")+
    scale_y_log10(breaks=seq(0.1, 10, 0.1))+
    geom_hline(yintercept=1, linetype="dashed", color = "darkred")+
    labs(x="age", y=expression(paste("exp(c), exp(", delta, ")")), title=paste(cou.j, sex))
  
  p
  
}  


## old version in which forecast uncertainty is not accounted, 
## but everything is analytic

# tF.ava <- t.ava.dea[!t.ava.dea%in%t1.ava]
# j=1
# list.out <- list()
# for(j in 1:length(tF.ava)){
#   
#   ## extract tF.ava[j] forecast log-mortality (~offset in this analysis)
#   eta.for0 <- subset(DFTR, type=="Baseline Logrates" & years==tF.ava[j])
#   eta.for <- eta.for0$value
# 
#   ## extract deaths and exposures for tF.ava[j] 
#   obs.y0 <- subset(DFTR, type=="Obs Deaths" & years==tF.ava[j])
#   obs.y <- obs.y0$value
#   obs.e0 <- subset(DFTR, type=="Obs Offset" & years==tF.ava[j])
#   obs.e <- obs.e0$value
#   ## age-structure
#   ag.low <- obs.y0$ages
# 
#   QIC_Pert <- function(par){
#     FIT <- PertubFUN(ag.low=ag.low,
#                      obs.y=obs.y,
#                      obs.e=obs.e,
#                      a=a,
#                      eta.for=eta.for,
#                      lambda=par)
#     FIT$QIC
#   }
#   ## optimizing lambdas using greedy grid search
#   OPT <- cleversearch(QIC_Pert, lower=-1, upper=8,
#                       ngrid=19, logscale=TRUE, startvalue=10^4,
#                       verbose=TRUE)
# 
#   ## estimating perturbation with optimal lambdas
#   PertOPT <- PertubFUN(ag.low=ag.low,
#                        obs.y=obs.y,
#                        obs.e=obs.e,
#                        a=a,
#                        eta.for=eta.for,
#                        lambda=OPT$par)
#   nb <- ncol(PertOPT$Ba)
#   ## fitted values
#   c.hat <- PertOPT$coeff[1]
#   expc.hat <- exp(c.hat)
#   delta.hat <- PertOPT$Ba%*%PertOPT$coeff[1:nb+1]
#   expdelta.hat <- exp(delta.hat)
#   eta.hat <- eta.for + c.hat + delta.hat # U%*%betas+eta20
#   ## se for c
#   se.c <- sqrt(PertOPT$Vbetas[1,1])
#   ## se for delta
#   V.delta <- cbind(0, PertOPT$Ba) %*% PertOPT$Vbetas %*% t(cbind(0, PertOPT$Ba))
#   se.delta <- sqrt(diag(V.delta))
#   ## se for exp(c)
#   der.expc <- matrix(c(exp(c.hat), rep(0, nb)), 1, nb+1)
#   V.expc <- der.expc %*% PertOPT$Vbetas %*% t(der.expc)
#   se.expc <- sqrt(diag(V.expc))
#   ## se for exp(delta)
#   der.expdelta0 <- cbind(0, PertOPT$Ba)
#   der.expdelta <- diag(c(exp(delta.hat))) %*% der.expdelta0
#   V.expdelta <- der.expdelta %*% PertOPT$Vbetas %*% t(der.expdelta)
#   se.expdelta <- sqrt(diag(V.expdelta))
#   ## se for eta
#   Veta.pert <- PertOPT$U %*% PertOPT$Vbetas %*% t(PertOPT$U)
#   se.eta.pert <- sqrt(diag(Veta.pert))
#   
#   ## 95% CIs
#   c.up <- c.hat + 1.96*se.c
#   c.low <- c.hat - 1.96*se.c
#   delta.up <- delta.hat + 1.96*se.delta
#   delta.low <- delta.hat - 1.96*se.delta
#   expc.up <- expc.hat + 1.96*se.expc
#   expc.low <- expc.hat - 1.96*se.expc
#   expdelta.up <- expdelta.hat + 1.96*se.expdelta
#   expdelta.low <- expdelta.hat - 1.96*se.expdelta
#   eta.up <- eta.hat + 1.96*se.eta.pert
#   eta.low <- eta.hat - 1.96*se.eta.pert
# 
#   DFTRdelta <- expand.grid(ages=a, years=tF.ava[j],
#                            type=c("Delta", "Exp Delta"))
#   DFTRdelta$value <- c(delta.hat, expdelta.hat)
#   DFTRdelta$up <- c(delta.up, expdelta.up)
#   DFTRdelta$low <- c(delta.low, expdelta.low)
#   
#   DFTRc <- expand.grid(ages=NA, years=tF.ava[j],
#                        type=c("c", "Exp c"))
#   DFTRc$value <- c(c.hat, expc.hat)
#   DFTRc$up <- c(c.up, expc.up)
#   DFTRc$low <- c(c.low, expc.low)
#   
#   DFTReta <- expand.grid(ages=a, years=tF.ava[j],
#                            type=c("Fitted Logrates"))
#   DFTReta$value <- c(eta.hat)
#   DFTReta$up <- c(eta.up)
#   DFTReta$low <- c(eta.low)
#   
#   list.out[[j]] <- rbind(DFTRdelta, DFTRc, DFTReta)
#   
#   ## observed
#   
#   eta.obs <- subset(DFTR, type=="Obs Logrates" & years==tF.ava[j])$value
#   DFg <- data.frame(ages=ag.low, 
#                     type="Actual grouped",
#                     eta1=eta.obs,
#                     ages.up=ag.up+1,
#                     eta1.up=eta.obs,
#                     eta1.low=NA)
#   ## forecast
#   DFfor <- data.frame(ages=a, 
#                       type="Forecast",
#                       eta1=eta.for,
#                       ages.up=NA,
#                       eta1.up=NA,
#                       eta1.low=NA)
#   ## fitted
#   DFhat <- data.frame(ages=a, 
#                       type="Fitted",
#                       eta1=eta.hat,
#                       ages.up=NA,
#                       eta1.up=eta.up,
#                       eta1.low=eta.low)
#   ## combine
#   DF <- rbind(DFg, DFfor, DFhat)
#   
#   p <- ggplot(DF, aes(x=ages, y=eta1, color=type)) +
#     geom_segment(data=filter(DF, type=="Actual grouped"),
#                  aes(x=ages, y=eta1, xend=ages.up, yend=eta1.up), size=1)+
#     geom_line(data=filter(DF, type=="Fitted"),
#               aes(y=eta1), size=1)+
#     geom_ribbon(data=filter(DF, type=="Fitted"),
#                 aes(ymin=eta1.low, ymax=eta1.up), alpha=.2)+
#     geom_line(data=filter(DF, type=="Forecast"),
#               aes(y=eta1),size=1.2)+
#     labs(x="age", y="log-mortality", title=paste(cou.j, sex))
#   p
#   
#   DFexpdelta <- data.frame(ages=a, expdelta=expdelta.hat,
#                            expdelta.low=expdelta.low,
#                            expdelta.up=expdelta.up)
#   DFexpc <- data.frame(x=-2, expc=expc.hat, expc.low=expc.low, expc.up=expc.up)
#   
#   p <- ggplot(DFexpdelta, aes(x=ages))+
#     geom_line(aes(y=expdelta), size=2, colour="darkgreen")+
#     geom_ribbon(aes(ymin=expdelta.low, ymax=expdelta.up), alpha=0.5)+
#     geom_pointrange(aes(x=DFexpc$x, y=DFexpc$expc, 
#                         ymin = DFexpc$expc.low, 
#                         ymax = DFexpc$expc.up),
#                     colour="darkblue", size=0.8)+
#     scale_y_log10(breaks=seq(0.1, 10, 0.1))+
#     geom_hline(yintercept=1, linetype="dashed", color = "darkred")+
#     labs(x="age", y=expression(paste("exp(c), exp(", delta, ")")), title=paste(cou.j, sex))
#   
#   p
#   
# }
# 
# DFTRpert <- dplyr::bind_rows(list.out)
# 
# DFTR <- rbind(DFTR, DFTRpert)


## END
