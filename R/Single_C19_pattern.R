
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
  dplyr::select(Country, Date, Age, AgeInt, Deaths = b,f,m) %>% 
    mutate(Date = dmy(Date),
           DateDiff = abs(Date - ymd("2020-12-31"))) %>% 
    group_by(Country) %>% 
    dplyr::filter(DateDiff == min(DateDiff)) %>% 
  mutate(AgeInt = case_when(Country == "Slovenia" & Age == 0 & AgeInt == 45L ~ 5L,
                            TRUE ~AgeInt)) %>% 
  dplyr::filter(!(Country == "Slovenia" & Age == 0 & AgeInt == 35L)) %>% 
  mutate(Deaths = ifelse(is.na(Deaths),0,Deaths)) %>% 
  group_by(Country) %>% 
  mutate(total_deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  dplyr::filter(total_deaths> 100) %>% 
  dplyr::select(-total_deaths)
 
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
  ungroup() %>% 
  dplyr::filter(total_deaths > 100) %>% 
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








