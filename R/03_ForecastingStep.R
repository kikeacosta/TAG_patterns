## clearing workspace
# these lines only run on GC's
if (system("whoami", intern = TRUE)  == "gccamarda") {
  rm(list = ls())
  ## plotting in a difference device
  options(device = "X11")
  ## !!!! to be changed
  setwd("~/WORK/TAG_patterns/")
}


library(MortalitySmooth)
library(magic)
library(colorspace)
library(ggplot2)
library(plotly)
library(viridis)
library(Matrix)
source("R/GroupedSmoothConstrainedMortalityForecasting_Functions.R")
source("R/cleversearch.R")
BICfun <- function(par, ag.low, Yg1, a, E1, WEIg1, infantj) {
  FIT <- PSinfantGrouped(
    a.low = ag.low,
    Yg = Yg1,
    a = a,
    E = E1,
    WEIg = WEIg1,
    lambdas = par,
    verbose = TRUE,
    kappa.shape = 0,
    infant = infantj
  )
  FIT$bic
}

## loading deaths
deaths <-
  read.csv("data_inter/deaths_sourced_infant_based.csv", header = TRUE)
## loading population
offset <- read.csv("data_inter/offsets.csv", header = TRUE)
## sex


## ages and dimensions
a <- 0:100
m <- length(a)
## pandemic years for the first forecasting part
t2 <- 2020:2021
n2 <- length(t2)
## countries
cou <- unique(deaths$Country)
nc  <- length(cou)

## questions:
## FEMALES


nocou <-
  cou %in% c("Andorra", "Faroe Islands", "Montserrat", "Seychelles")

# nocou <- c(3,29,61,77)
# #unique(deaths$Country[which(deaths$age_spn==5)])
# yescou <- 1:nc
# yescou <- yescou[-nocou]

cou[nocou]
#unique(deaths$Country[which(deaths$age_spn==5)])
# yescou = cou; nocou = c(0)
yescou <- 1:nc
yescou <- yescou[!nocou]


get_forecast_rough <- function(deaths.j, offset) {
  cou.j    <- deaths.j$Country[1]
  sex.s    <- deaths.j$Sex[1]
  ## ages and dimensions
  a <- 0:100
  m <- length(a)
  ## pandemic years for the first forecasting part
  t2 <- 2020:2021
  n2 <- length(t2)
  
  ## check available pre-pandemic years, country-dependent
  t1.ava   <- unique(deaths.j$Year)[unique(deaths.j$Year) < t2[1]]
  
  ## missing years
  t1       <- min(t1.ava):max(t1.ava)
  n1       <- length(t1)
  whi.ava  <- t1 %in% t1.ava
  
  ## select t1 from deaths.j
  Y.j      <- subset(deaths.j, Year %in% t1.ava)
  
  ## all years
  t        <- c(t1, t2)
  n        <- n1 + n2
  ag.low   <- unique(Y.j$Age)
  ag.up1   <- ag.low - 1
  if (ag.low[2] == 1) {
    ag.up <- c(0, ag.up1[ag.up1 > 0], max(a))
  } else{
    ag.up <- c(ag.up1[ag.up1 > 0], max(a))
  }
  cbind(ag.low, ag.up)
  lg <- ag.up - ag.low + 1
  mg <- length(ag.low)
  ## only for plotting/aesthetic
  ag.mid <- (ag.low + ag.up) / 2
  ag.lab <- paste(ag.low, ag.up, sep = "-")
  Yg1 <- matrix(0, mg, n1)
  Yg1[, whi.ava] <- matrix(Y.j$Deaths, mg, length(t1.ava))
  
  ## working on the offset
  E.j <-
    subset(offset, Country == cou.j & Year %in% t1.ava & Sex == sex.s)
  E1 <- matrix(0, m, n1)
  E1[, whi.ava] <- matrix(E.j$Population, m, length(t1.ava))
  ## matrix to map grouped to single age
  G <- matrix(0, mg, m)
  rownames(G) <- ag.mid
  colnames(G) <- a
  for (i in 1:mg) {
    ag.low.i <- ag.low[i]
    ag.up.i <- ag.up[i]
    wc <- which(a >= ag.low.i & a <= ag.up.i)
    G[i, wc] <- 1
  }
  # all(colSums(G)==1)
  ## grouping exposures only for plotting
  Eg1 <- G %*% E1
  ## observed log-mortality, by age-groups, only for plotting and e0
  ETAg1 <- log(Yg1 / Eg1)
  
  # ## plotting grouped actual log-mortality over ages
  DFg <- expand.grid(list(ages = ag.low, years1 = t1))
  DFg$type    <- "Actual grouped"
  DFg$eta1    <- c(ETAg1)
  DFg$ages.up <- ag.up + 1
  DFg$eta1.up <- c(ETAg1)
  DFg$ag.lab  <- factor(c(rep(ag.lab, n1)), levels = ag.lab)
  
  ## weights for grouped-exposures
  WEIg1            <- matrix(1, mg, n1)
  WEIg1[, !whi.ava] <- 0
  ## is there an age 0?
  infantj <- ag.up[1] == 0
  
  
  ## optimizing lambdas using greedy grid search
  # tim1 <- Sys.time()
  OPT <- cleversearch2(
    BICfun,
    lower = c(-1, 1),
    upper = c(3, 5),
    ngrid = 5,
    logscale = TRUE,
    startvalue = c(10 ^ 4, 10 ^ 5),
    verbose = TRUE,
    ag.low = ag.low,
    Yg1 = Yg1,
    a = a,
    E1 = E1,
    WEIg1 = WEIg1,
    infantj = infantj
  )
  # tim2 <- Sys.time()
  # tim2-tim1
  
  
  ## estimating mortality with optimal lambdas
  FITi <- PSinfantGrouped(
    a.low = ag.low,
    Yg = Yg1,
    a = a,
    E = E1,
    WEIg = WEIg1,
    lambdas = OPT$par,
    verbose = TRUE,
    kappa.shape = 0,
    infant = infantj
  )
  
  ## estimated/smooth log-mortality on the pre-pandemic years
  ETA1hat <- FITi$ETA
  
  ## FORECASTING STEP ## FORECASTING STEP ## FORECASTING STEP
  ## where to apply the constraints
  S <- matrix(1, m, n)
  S[, 1:n1] <- 0
  ## compute deltas : CIs of relative derivatives
  deltasi <- deltasFUN(FITi)
  ## forecasting
  CPSg.i <- CPSfunctionGrouped(
    a.low = ag.low,
    Yg1 = Yg1,
    WEIg1 = WEIg1,
    a = a,
    E1 = E1,
    obs.yrs = t1,
    for.hor = max(t),
    ETA1hat = ETA1hat,
    deltas = deltasi,
    S = S,
    lambdas = OPT$par,
    verbose = TRUE,
    infant = infantj,
    kappa.shape = 0
  )
  
  ETAhat <- CPSg.i$ETA
  ## plotting, starting from the beginning
  ## actual log-mortality
  DFg <- expand.grid(list(ages = ag.low, years1 = t1))
  DFg$type <- "Actual grouped"
  DFg$eta1 <- c(ETAg1)
  DFg$ages.up <- ag.up + 1
  DFg$eta1.up <- c(ETAg1)
  ## fitted
  DFhat <- expand.grid(list(ages = a, years1 = t1))
  DFhat$type <- "Fitted"
  DFhat$eta1 <- c(ETA1hat)
  DFhat$ages.up <- NA
  DFhat$eta1.up <- NA
  ## forecast
  DFfor <- expand.grid(list(ages = a, years1 = t))
  DFfor$type <- "Fitted+Forecast"
  DFfor$eta1 <- c(ETAhat)
  DFfor$ages.up <- NA
  DFfor$eta1.up <- NA
  ## all
  DF <- rbind(DFg, DFhat, DFfor)
  DF$group <- factor(c(rep(ag.lab, n1), rep(rep(ag.lab, lg), n1),
                       rep(rep(ag.lab, lg), n)), levels = ag.lab)
  DF$colo <- c(rep(ag.mid, n1), rep(a, n1), rep(a, n))
  
  DF
}

deaths.j <- deaths |> subset(sex == "m" & Country == "Albania")
DFj <- get_forecast_rough(deaths.j, offset)

library(foreach)




# library(multidplyr)
# n_cores <- ceiling(parallel::detectCores() / 2)
# cluster <- new_cluster(n_cores)
# cluster_call(source("R/GroupedSmoothConstrainedMortalityForecasting_Functions.R"))
# cluster_library(cluster,c("tidyverse","Matrix","MortalitySmooth"))
# cluster_copy(cluster,c("get_forecast_rough",
#                        "offset",
#                        "cleversearch2",
#                        "BICfun","PSinfant","deltasFUN","CPSfunctionGrouped"))
# for (c in )
attempt <-
  deaths |>
  group_by(Country, Sex) %>%
  do(get_forecast_rough(., offset = offset))


## END
