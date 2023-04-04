# this is a selective redux of 04_PerturbationStep.R,
# aimed at wrapping everything in a function; or, if possible,
# using modular functions, although that's more ambitious.

library(MortalitySmooth)
library(magic)
library(colorspace)
library(ggplot2)
library(plotly)
library(viridis)
library(Matrix)
library(tidyverse)
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





deaths <- read.csv("data_inter/deaths_sourced_infant_based_99.csv", header=TRUE)
## loading population
offset <- read.csv("data_inter/offsets_99.csv", header=TRUE)

# start body here
fit_excess <- function(deaths.j, offset){
  
  
  cou.j    <- deaths.j$Country[1]
  sex.s    <- deaths.j$Sex[1]
  deaths.j <- deaths.j |> 
    filter(Year < 2022)
  cat(paste(cou.j, sex.s,"\n"))
  ## ages and dimensions
  a <- 0:99
  m <- length(a)
  ## pandemic years for the first forecasting part
  t2 <- 2020:2021
  n2 <- length(t2)
  
  ## check available pre-pandemic years, country-dependent
  t1.ava   <- deaths.j |> 
    filter(Year < t2[1]) |> 
    pull(Year) |> 
    unique()
  
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
  # cbind(ag.low, ag.up)
  lg     <- ag.up - ag.low + 1
  mg     <- length(ag.low)
  ## only for plotting/aesthetic
  ag.mid <- (ag.low + ag.up) / 2
  ag.lab <- paste(ag.low, ag.up, sep = "-")
  Yg1    <- matrix(0, mg, n1)
  Yg1[, whi.ava] <- matrix(Y.j$Deaths, mg, length(t1.ava))
  
  ## select offset
  E.j <-
    deaths.j |> 
    #filter(Year < t2[1]) |> 
    select(Country, Year, Sex) |> 
    distinct() |> 
    left_join(offset, by = c("Country","Year","Sex"),multiple = "all")
  
  E1.1 <- E.j |> 
    filter(Year < t2[1]) |> 
    pivot_wider(names_from = Year, values_from = "Population") |> 
    select(-Country, -Sex) |> 
    column_to_rownames("Age") |> 
    as.matrix()
  
  # robustness step to allow for gaps:
  E1 <- matrix(0,m,length(t1), dimnames=list(a,t1))
  E1[,colnames( E1.1 )] <-  E1.1 
  
  
  ## matrix to map grouped to single age
  
  A1 <- tibble(age = a, ones = rep(1, m))
  B  <- tibble(age = ag.low,value=1)
  G  <- left_join(A1, B, by = "age") |> 
    mutate(ag = age * value) |> 
    fill(ag,.direction = "down") |> 
    select(-value) |> 
    pivot_wider(names_from=age, values_from = ones, values_fill = 0) |> 
    column_to_rownames("ag") |> 
    as.matrix()
  # all(colSums(G)==1)
  ## grouping exposures only for plotting
  # Eg1 <- G %*% E1
  # ## observed log-mortality, by age-groups, only for plotting and e0
  # ETAg1 <- log(Yg1 / Eg1)
  # 
  # # ## plotting grouped actual log-mortality over ages
  # DFg <- expand.grid(list(ages = ag.low, years1 = t1))
  # DFg$type    <- "Actual grouped"
  # DFg$eta1    <- c(ETAg1)
  # DFg$ages.up <- ag.up + 1
  # DFg$eta1.up <- c(ETAg1)
  # DFg$ag.lab  <- factor(c(rep(ag.lab, n1)), levels = ag.lab)
  # 
  ## weights for grouped-exposures
  WEIg1            <- matrix(1, mg, n1)
  WEIg1[, !whi.ava] <- 0
  ## is there an age 0?
  infantj <- TRUE
  
  
  ## optimizing lambdas using greedy grid search
  # tim1 <- Sys.time()
  invisible(capture.output(
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
  ))
  # FITi is used to compute the derivative over age and time,
  # it's a smoothing of the prepandemic years
  invisible(capture.output(
    FITi <- PSinfantGrouped(a.low = ag.low,
                            Yg=Yg1,
                            a=a,
                            E=E1,
                            WEIg=WEIg1,
                            lambdas=OPT$par,
                            verbose=TRUE,
                            kappa.shape=10^4,
                            infant=infantj)
  ))
  ## estimated/smooth log-mortality on the pre-pandemic years
  ETA1hat <- FITi$ETA
  # over age and time of prepandemic years
  
  # TR: would be nice to add dimnames to this...
  
  ## FORECASTING STEP ## FORECASTING STEP ## FORECASTING STEP 
  ## where to apply the constraints
  S <- matrix(1, m, n)
  S[,1:n1] <- 0
  ## compute deltas : CIs of relative derivatives
  deltasi <- deltasFUN(FITi)
  ## forecasting
  invisible(capture.output(
    CPSg.i <- CPSfunctionGrouped(a.low=ag.low, Yg1=Yg1, WEIg1=WEIg1,
                                 a=a, E1=E1, obs.yrs=t1, 
                                 for.hor=max(t),
                                 ETA1hat=ETA1hat, deltas=deltasi, S=S, 
                                 lambdas = OPT$par, verbose = TRUE,
                                 infant = infantj,
                                 kappa.shape=10^4)
  ))
  
  # this covers prepandemic until 2020 or 2021
  ETAhat <- CPSg.i$ETA
  
  
  
  ## analytic confidence intervals
  B      <- kronecker(CPSg.i$Bt, CPSg.i$Ba)
  Veta   <- B %*% CPSg.i$Valpha %*% t(B)
  SEeta  <- matrix(sqrt(diag(Veta)), m, n)
  
  ETAup  <- ETAhat + 2*SEeta
  ETAlow <- ETAhat - 2*SEeta
  
  ## obs deaths
  t.ava.dea <- unique(deaths.j$Year)
  
  # obs deaths and log rates together
  DFTR.obsdeaths <-
    deaths.j |> 
    select(ages = Age, 
           years = Year, 
           `Obs Deaths` = Deaths, 
           Population) |> 
    mutate(`Obs Logrates` = log(`Obs Deaths`/Population),up = NA, low = NA) |> 
    select(-Population) |> 
    pivot_longer(c(`Obs Deaths` ,`Obs Logrates`),names_to ="type", values_to = "value") |> 
    arrange(type,years, ages)
  
  
  ## obs offset
  min.t.ava <- DFTR.obsdeaths |> pull(years) |> min()
  
  DFTR.offset <- offset |> 
    filter(Country == cou.j,
           between(Year, min.t.ava, 2021),
           Sex == sex.s) |> 
    select(-Country, -Sex) |> 
    rename(years = Year, ages = Age, value = Population) |> 
    mutate(type="Obs Offset",
           up = NA,
           low = NA)
  
  
  ## baseline logrates
  
  DFTR.baselogrates <- expand.grid(list(ages=a,
                                        years=t)) |> 
    mutate(type = "Baseline Logrates",
           value =  c(ETAhat),
           up =  c(ETAup),
           low = c(ETAlow))
  
  ## append all DFTR....
  DFTR <- bind_rows(DFTR.obsdeaths,
                    DFTR.offset,
                    DFTR.baselogrates)
  
  tF.ava <- t.ava.dea[!t.ava.dea%in%t1.ava]
  
  list.out <- list()
  for(j in 1:length(tF.ava)){
    
    ## extract tF.ava[j] forecast log-mortality (~offset in this analysis)
    eta.for0 <- subset(DFTR, type=="Baseline Logrates" & years==tF.ava[j])
    eta.for <- eta.for0$value
    
    ## extract deaths and exposures for tF.ava[j] 
    obs.y0 <- subset(DFTR, type=="Obs Deaths" & years==tF.ava[j])
    obs.y <- obs.y0$value
    obs.e0 <- subset(DFTR, type=="Obs Offset" & years==tF.ava[j])
    obs.e <- obs.e0$value
    ## age-structure
    ag.low <- obs.y0$ages
    
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
    invisible(capture.output(
      OPT <- cleversearch(QIC_Pert, lower=-1, upper=8,
                          ngrid=19, logscale=TRUE, startvalue=10^4,
                          verbose=TRUE)
    ))
    ## estimating perturbation with optimal lambdas
    invisible(capture.output(
      PertOPT <- PertubFUN(ag.low=ag.low,
                           obs.y=obs.y,
                           obs.e=obs.e,
                           a=a,
                           eta.for=eta.for,
                           lambda=OPT$par)
    ))
    nb            <- ncol(PertOPT$Ba)
    ## fitted values
    c.hat         <- PertOPT$coeff[1]
    expc.hat      <- exp(c.hat)
    delta.hat     <- PertOPT$Ba%*%PertOPT$coeff[1:nb+1]
    expdelta.hat  <- exp(delta.hat)
    eta.hat       <- eta.for + c.hat + delta.hat # U%*%betas+eta20
    ## se for c
    se.c          <- sqrt(PertOPT$Vbetas[1,1])
    ## se for delta
    V.delta       <- cbind(0, PertOPT$Ba) %*% PertOPT$Vbetas %*% t(cbind(0, PertOPT$Ba))
    se.delta      <- sqrt(diag(V.delta))
    ## se for exp(c)
    der.expc      <- matrix(c(exp(c.hat), rep(0, nb)), 1, nb+1)
    V.expc        <- der.expc %*% PertOPT$Vbetas %*% t(der.expc)
    se.expc       <- sqrt(diag(V.expc))
    ## se for exp(delta)
    der.expdelta0 <- cbind(0, PertOPT$Ba)
    der.expdelta  <- diag(c(exp(delta.hat))) %*% der.expdelta0
    V.expdelta    <- der.expdelta %*% PertOPT$Vbetas %*% t(der.expdelta)
    se.expdelta   <- sqrt(diag(V.expdelta))
    ## se for eta
    Veta          <- PertOPT$U %*% PertOPT$Vbetas %*% t(PertOPT$U)
    se.eta        <- sqrt(diag(Veta))
    
    ## 95% CIs
    c.up          <- c.hat + 1.96*se.c
    c.low         <- c.hat - 1.96*se.c
    delta.up      <- delta.hat + 1.96*se.delta
    delta.low     <- delta.hat - 1.96*se.delta
    expc.up       <- expc.hat + 1.96*se.expc
    expc.low      <- expc.hat - 1.96*se.expc
    expdelta.up   <- expdelta.hat + 1.96*se.expdelta
    expdelta.low  <- expdelta.hat - 1.96*se.expdelta
    eta.up        <- eta.hat + 1.96*se.eta
    eta.low       <- eta.hat - 1.96*se.eta
    
    DFTRdelta     <- expand.grid(ages=a, years=tF.ava[j],
                                 type=c("Delta", "Exp Delta"))
    DFTRdelta$value <- c(delta.hat, expdelta.hat)
    DFTRdelta$up  <- c(delta.up, expdelta.up)
    DFTRdelta$low <- c(delta.low, expdelta.low)
    
    DFTRc <- expand.grid(ages=NA, years=tF.ava[j],
                         type=c("c", "Exp c"))
    DFTRc$value   <- c(c.hat, expc.hat)
    DFTRc$up      <- c(c.up, expc.up)
    DFTRc$low     <- c(c.low, expc.low)
    
    DFTReta <- expand.grid(ages=a, years=tF.ava[j],
                           type=c("Fitted Logrates"))
    
    DFTReta$value <- c(eta.hat)
    DFTReta$up    <- c(eta.up)
    DFTReta$low   <- c(eta.low)
    
    list.out[[j]] <- rbind(DFTRdelta, DFTRc, DFTReta)
  }
  
  DFTRpert <- dplyr::bind_rows(list.out)
  
  DFTR <- rbind(DFTR, DFTRpert)
  DFTR
}


# source("R/temp_function.r")
# for spot testing
deaths.j <- 
deaths |> 
  filter(Country=="Andorra" & Sex == "m") 
  fit_excess(offset = offset)

# "Andorra" computationally singular
library(parallel)

fit_excess_wrap <- function(deaths.j, offset){
  out <- try(fit_excess(deaths.j = deaths.j,
                 offset = offset))
  if (class(out) == "try-error"){
    return(NULL)
  } else {
    return(out)
  }
}

big_L <-
deaths |> 
  filter(Sex != "t") |> 
  group_by(Country, Sex) |> 
  group_split()

N <- length(big_L)
out_L <- vector("list", N)
for (n in 1:N){
  ctr  <- big_L[[n]]$Country[1]
  sx   <- big_L[[n]]$Sex[1]
  outn <- try(fit_excess(big_L[[n]], offset = offset))
   if (class(outn)[1] == "try-error"){
     outn <- tibble(ages=integer(),
                    years=integer(),
                    up=double(),
                    low=double(),
                    type=character(),
                    value=double())
   }
   outn <-outn |> 
     mutate(Country = ctr,.before=1) |> 
     mutate(Sex = sx,.after=1)
     
   out_L[[n]] <- outn
}

big_test <- 
  out_L |> 
  bind_rows() |> 
  mutate(Sex = case_when(Sex == "m" ~ "male",
                         Sex == "f" ~ "female",
                         Sex == "t" ~ "total"))

big_test |> 
  write_csv("data_inter/eta_all.csv")


yr <- 2021
ctry <- big_test |> 
  filter(years == yr) |> 
  pull(Country) |> unique()

pdf(paste0("Figures/Fitted",yr,"all.pdf"))

for (cou.j in ctry){

  sexes <- big_test |> 
    filter(Country == cou.j,
           years == yr) |> 
    pull(Sex) |> 
    unique()
  for (sex in sexes){
    DF <-
      big_test |> 
      filter(Country == cou.j,
             Sex == sex,
             years == yr)
    
    DFobs <- DF |> 
      filter(type == "Obs Logrates") |> 
      mutate(ageup = ages + DemoTools::age2int(ages, OAvalue = 1)) |> 
      ungroup()
    
    p <-
      ggplot(DF, aes(x = ages, y = value, color = type)) +
      geom_segment(data = DFobs,
                   aes(x = ages, y = value, xend = ageup, yend = value), linewidth = 1) +
      geom_line(data = filter(DF, type == "Fitted Logrates"), linewidth = 1) +
      geom_ribbon(data = filter(DF, type == "Fitted Logrates"),
                  aes(ymin = low, ymax = up), alpha = .2) +
      geom_line(data = filter(DF, type == "Baseline Logrates"),
                aes(y = value), linewidth = 1.2) +
      labs(x = "age", y = "log-mortality", title = paste(cou.j,sex)) +
      theme_minimal()
    print(p)
  }
}

dev.off()

yr <- 2020
ctry <- big_test |> 
  filter(years == yr) |> 
  pull(Country) |> unique()

pdf(paste0("Figures/Deltas",yr,"all.pdf"))

for (cou.j in ctry){
  
  sexes <- big_test |> 
    filter(Country == cou.j,
           years == yr) |> 
    pull(Sex) |> 
    unique()
  for (sex in sexes){
    DFa <-
      big_test |> 
      filter(Country == cou.j,
             Sex == sex,
             years == yr,
             type %in% c("c","Delta"))
    
    .c  <-  DFa |> 
      filter(type == "c") |> 
      mutate(ages=0)
    
    DFd <- DFa |> 
      filter(type == "Delta") 
    
    p <-
      ggplot(DFa, aes(x = ages, y = value, color = type)) +
      geom_line(data=DFd,linewidth = 1) +
      geom_ribbon(data=DFd,mapping=aes(ymin = low, ymax = up), alpha = .2) +
      geom_point(data=.c, mapping=aes(x=ages,y=value),color="red",size=2) +
      geom_segment(data=.c, mapping=aes(x=ages,xend=ages,y=low,yend=up),color="red") +
      labs(x = "age", y = "log difference", title = paste(cou.j,sex)) +
      theme_minimal()
    print(p)
  }
}

dev.off()

