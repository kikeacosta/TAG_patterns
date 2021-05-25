
make_C <- function(Age, AgeInt, exposure){
  
 
  Agei  <- Age + 1
  m     <- sum(AgeInt)
  
  k     <- length(Age)
   
   if (length(exposure) == length(Age)){
     exposure <- rep(exposure, AgeInt)
   }
  if (length(exposure) > m){
    exposure[m] <- sum(exposure[m:length(exposure)])
    exposure <- exposure[1:m]
  }
   
   stopifnot(length(exposure) == m)
   
   
  C_out <- matrix(0,
                  nrow = k, 
                  ncol = m, 
                  dimnames = list(Age, 0:(sum(AgeInt) - 1)))
  
  Agei_top <- c(Agei[-1] - 1, m)  
  for (i in 1:k){
    ind <- Agei[i]:Agei_top[i]
  # cat(ind,"\n\n")  
   C_out[i,ind] <- exposure[ind]
  }
  
  C_out
}

# Note: repo name != package name!
# remotes::install_github("eshom/covid-age-data")
library(covidAgeData) # for COVerAGE-DB
library(tidyverse)
library(lubridate)
library(osfr)
# download the inputDB without loading it
download_covid("inputDB", dest = "Data", download_only = TRUE)


C19deaths <- read_subset_covid(zippath = "Data/inputDB.zip",
                  data = "inputDB",
                  Region = "All") %>% 
  dplyr::filter(Measure == "Deaths",
                !(Metric == "Fraction" & Age == "TOT")) %>% 
  pivot_wider(names_from = Sex, values_from = Value) %>% 
  mutate(b = case_when(is.na(b) & !is.na(f) & !is.na(m) ~ f + m,
                       TRUE ~ b)) %>% 
  select(Country, Date, Age, AgeInt, Deaths = b) %>% 
    mutate(Date = dmy(Date),
           DateDiff = abs(Date - ymd("2020-12-31"))) %>% 
    group_by(Country) %>% 
    dplyr::filter(DateDiff == min(DateDiff),
                  Date >= ymd("2020-08-01")) %>% 
  mutate(AgeInt = case_when(Country == "Slovenia" & Age == 0 & AgeInt == 45L ~ 5L,
                            TRUE ~AgeInt)) %>% 
  dplyr::filter(!(Country == "Slovenia" & Age == 0 & AgeInt == 35L)) %>% 
  mutate(Deaths = ifelse(is.na(Deaths),0,Deaths))
 
# to avoid some preprocessing just now, we can
# keep it to both-sex data


# No need to convert fractions since we'll rescale to Hopkins totals anyway
# unique(C19deaths$Metric)
# C19deaths %>% 
#   dplyr::filter(Metric == "Fraction")

# update coronavirus dataset in coronavirus package
OWD <- read_csv("https://covid.ourworldindata.org/data/owid-covid-data.csv") 
C19_total_deaths <-
  OWD %>% 
  dplyr::filter(date == as_date("2021-01-01"),
                !is.na(continent)) %>% 
  select(total_deaths, 
         Country = location) %>% 
  mutate(Country = ifelse(Country ==  "United States" ,"USA",Country))
countries_have <- unique(C19deaths$Country)

# Fine, we just do UK for now, haha
countries_have[!countries_have%in%C19_total_deaths$Country]

# Join
C19_use <-
  C19deaths %>% 
  filter(Age != "TOT", 
         Age != "UNK") %>% 
  left_join(C19_total_deaths, by = "Country") %>% 
  group_by(Country) %>% 
  mutate(total_deaths = case_when(is.na(total_deaths) ~ sum(Deaths),
                                  TRUE ~ total_deaths)) %>% 
  dplyr::filter(total_deaths > 200) %>% 
  group_by(Country) %>% 
  mutate(Deaths = Deaths / sum(Deaths) * total_deaths,
         Age = as.integer(Age))

# Finally need ideally single age pops for these countries,
# let's see what we got:

osf_code <-"unf6v"
osf_retrieve_file(osf_code) %>%
  osf_download(conflicts = "overwrite",
               path = "Data")
offsets <-  read_csv("Data/offsets.csv", skip = 1) %>% 
  dplyr::filter(Region == "All") %>% 
  pivot_wider(names_from = Sex, values_from = Population) %>% 
  mutate(b = case_when(is.na(b) & !is.na(f) & !is.na(m) ~ f + m,
                       TRUE ~ b)) %>% 
  select(-Region, Population = b, Age)

# TODO: get population counts for Moldova and UK!!!

C19_use <- 
  C19_use %>% 
  dplyr::filter(Country %in% offsets$Country)

countries_loop <- C19_use$Country %>% unique()


## clearing workspace
rm(list = ls())
## plotting in a difference device
options(device="X11")
## R-studio to get the same dir as for the .R file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

Y <- readRDS("C19_use.rds", refhook = NULL)
## delete last row for Argentina
Y <- Y[-115,]
E0 <- readRDS("offsets.rds", refhook = NULL)
## select E where we have data in D
E <- subset(E0, Country%in%unique(Y$Country))

names(table(Y$Country))
names(table(E$Country))
## divide data by regions
EUR <- c("Austria", "Belgium", "Czechia", "Denmark", 
         "Finland","France", "Germany", "Greece", "Hungary",
         "Ireland", "Israel", "Italy", "Netherlands", "Northern Ireland", "Norway",
         "Portugal","Romania", "Scotland", "Slovenia",
         "Spain", "Sweden","Switzerland", "Turkey","Ukraine")
AFR <- c("Chad", "Eswatini", "Kenya", "Malawi","Nigeria","Somalia", "South Africa")
ASI <- c("Afghanistan", "Bangladesh", "China", "India", "Indonesia", "Iraq", "Jordan",
         "Japan", "Nepal","Pakistan","Philippines", "South Korea")
SAM <- c("Argentina", "Brazil", "Chile", "Colombia", "Costa Rica", "Cuba", "Ecuador",
         "Mexico",
         "Panama", "Paraguay", "Peru", "Suriname", "Uruguay")
NAM <- c("Canada", "Jamaica", "USA")
OCE <- c("Australia")
length(EUR)+length(AFR)+length(ASI)+length(SAM)+length(NAM)+length(OCE)
length(unique(Y$Country))





## only for a given region
REG <- ASI
Y.REG <- subset(Y, Country%in%REG)
E.REG <- subset(E, Country%in%REG)
## number of age-group within each pop
n <- table(Y.REG$Country)
## number of pop
nc <- length(n)
## ages
x <- 0:104
m <- length(x)

## building country-specific composite matrix
## and place them in a list
Clist <- Glist <- LMX <- Egr <- list()
i=4
for(i in 1:nc){
  E.i <- subset(E.REG, Country==REG[i])
  Y.i <- subset(Y.REG, Country==REG[i])
  e.i <- E.i$Population
  age.low.i <- Y.i$Age
  age.up.i <- Y.i$Age+Y.i$AgeInt-1
  age.gr.i <- paste(age.low.i, age.up.i, sep="-")
  C.i <- matrix(0, n[i], m)
  rownames(C.i) <- age.gr.i
  colnames(C.i) <- x 
  G.i <- C.i
  for(j in 1:n[i]){
    age.low.ij <- age.low.i[j]
    age.up.ij <- age.up.i[j]
    whi <- which(x>=age.low.ij & x<=age.up.ij)
    C.i[j,whi] <- e.i[whi]
    G.i[j,whi] <- 1
  }
  Clist[[i]] <- C.i
  Glist[[i]] <- G.i
  egr.i <- G.i%*%e.i
  Egr[[i]] <- egr.i
  LMX[[i]] <- log(Y.i$Deaths/egr.i)
}
## composite matrix for the overall model
library(magic)
C <- do.call("adiag",Clist)

## response
y <- Y.REG$Deaths

## penalty stuff
D <- diff(diag(m), diff=2)
tDD <- t(D)%*%D
lambda <- 10^5
P0 <- lambda*tDD
P <- adiag(P0,0*diag(nc-1))

## model matrix
U0 <- kronecker(rep(1,nc), diag(m))
U1 <- kronecker(diag(nc), rep(1,m))
U1 <- U1[,-1]
U <- cbind(U0, U1)

## starting eta
library(MortalitySmooth)
eta.st.list <- list()
i=1
for(i in 1:nc){
  E.i <- subset(E.REG, Country==REG[i])
  Y.i <- subset(Y.REG, Country==REG[i])
  d.i <- Y.i$Deaths
  e.i <- E.i$Population
  e.i[e.i==0] <- 0.5
  len.i <- Y.i$AgeInt
  dst.i <- rep(d.i/len.i, len.i)
  fit0.i <- Mort1Dsmooth(x=x, y=dst.i,
                         offset=log(e.i),
                         method=3, lambda=10^5)
  eta.st.i <- fit0.i$logmortality
  eta.st.list[[i]] <- eta.st.i
}
eta.st <- unlist(eta.st.list)

## PCLM regression
eta <- eta.st
max.it <- 100
for(it in 1:max.it){
  gamma   <- exp(eta)
  mu      <- c(C %*% gamma)
  X       <- (C * ((1 / mu) %*% t(gamma)) ) %*% U
  w       <- as.vector(mu)
  r       <- y - mu + C %*% (gamma * eta)
  G       <- t(X) %*% (w * X) 
  GpP     <- G + P
  tXr     <- t(X) %*% r
  beta    <- solve(GpP, tXr) 
  eta.old <- eta
  eta     <- U%*%beta
  dif.eta <- max(abs((eta - eta.old)/eta.old) )
  if(dif.eta < 1e-04 & it > 4) break
  cat(it, dif.eta, "\n")
}
## fitted log-mortality
eta.hat <- list()
for(i in 1:nc){
  eta.hat[[i]] <- eta[1:m+(i-1)*m]
}

theta.hat <- beta[-c(1:m)]

## standard errors
H0 <- solve(GpP)
H1 <- H0 %*% G %*% H0
se.theta <- sqrt(diag(H1)[-c(1:m)])
se <- sqrt(diag(U %*% H1 %*% t(U)))

## upper and lower bound for eta.hat and theta.hat
eta.hatL <- eta.hatU <- list()
se.hat <- list()
for(i in 1:nc){
  se.hat[[i]] <- se[1:m+(i-1)*m]
  eta.hatL[[i]] <- eta.hat[[i]] - 2*se.hat[[i]]
  eta.hatU[[i]] <- eta.hat[[i]] + 2*se.hat[[i]]
}
theta.hatL <- theta.hat - 2*se.theta
theta.hatU <- theta.hat + 2*se.theta

library(colorspace)
colocou <- rainbow_hcl(nc)
## plotting outcomes
rany <- c(-15, -3)
ranx <- range(x)
plot(1,1, t="n", col=2, lwd=2, ylim=rany, xlim=ranx,
     axes=FALSE,
     xlab="age", ylab="mortality, log-scale",
     main="OCEANIA C19 mortality age-pattern")
yy <- 10^seq(-7, -2)
axis(2, at=log(yy), labels=yy)
axis(1);box()
xx <- c(x, rev(x))
for(i in 1:nc){
  Y.i <- subset(Y.REG, Country==REG[i])
  age.low.i <- Y.i$Age
  age.up.i <- Y.i$Age+Y.i$AgeInt
  lmx.i <- LMX[[i]]
  for(j in 1:n[i]){
    segments(x0=age.low.i[j], x1=age.up.i[j],
             y0=lmx.i[j], y1=lmx.i[j], col=colocou[i], lwd=3)
  }
  yy <- c(eta.hatL[[i]], rev(eta.hatU[[i]]))
  polygon(xx, yy, col=adjustcolor(colocou[i], 0.5), border=adjustcolor(colocou[i], 0.5))
  lines(x, eta.hat[[i]], col=colocou[i], lwd=4)
}
legend("topleft", inset=0.01,
       legend=REG, col=colocou, pch=15, lwd=2)

plot(c(0,theta.hat), 1:nc, pch=15, cex=1.2, col=colocou,
     ylim=c(1,nc+0.5), xlim=range(0, theta.hatL,theta.hatU), axes=FALSE,
     xlab="", ylab="")
axis(1);box()
abline(v=0, col=8, lwd=2, lty=3)
grid(ny=NA)
for(i in 2:nc){
  arrows(x0=theta.hatL[i-1], y0=i, x1 = theta.hatU[i-1], col=colocou[i], 
         lwd=2, code=3, length=0.1, angle=90)
}
for(i in 1:nc){
  text(c(0,theta.hat)[i], y=i+0.5, REG[i], col=colocou[i])
}







C_big <- list()
  for (i in countries_loop){
    
    C19d <- C19_use %>% dplyr::filter(Country == i, Age < 105)
    exposure <- offsets %>% 
      dplyr::filter(Country == i) %>% 
      dplyr::pull(Population)
    
    C_i <- try(make_C(Age = C19d$Age, AgeInt = C19d$AgeInt, exposure = exposure)) 
    C_big[[i]] <- C_i
  }
errors <- lapply(C_big,class) %>% unlist()
errors[errors == "try-error"]

# This is how we stack them, right?
C_big <- do.call("rbind",C_big)

# Now you take it from here!



