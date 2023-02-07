# Covariates

library(tidyverse)
library(readr)
library(HMDHFDplus)
library(MortalityLaws)

emp.der <- function(x, y){
  # simple function for "empirical" derivatives\
  # courtesy of GC
  m <- length(x)
  a <- c(diff(y), 1)
  b <- c(diff(x), 1)
  ab <- a/b
  wei <- c(rep(1, m-1), 0)
  
  a1 <- c(1, -1*diff(y))
  b1 <- c(1, -1*diff(x))
  ab1 <- a1/b1
  wei1 <- c(0, rep(1, m-1))
  
  y1emp <- (ab*wei + ab1*wei1)/(wei+wei1)
  return(y1emp)
}

WPP <- read_csv("Data/WPP2019_Life_Table_Medium.csv")
WPP$MidPeriod %>% unique()
WPP <- 
  WPP %>% 
  filter(MidPeriod == 2018)

gc()

WPP_calcs <-
  WPP %>% 
  arrange(LocID, SexID, AgeGrpStart) %>% 
  group_by(LocID, Location, Sex) %>% 
  mutate(ex_mid = approx(x = AgeGrpStart, y = ex, xout = (AgeGrpStart + ax), rule = 2)$y) %>% 
  summarize(e0 = ex[AgeGrpStart == 0],
            edagger = sum(ex_mid * dx, na.rm = TRUE) / sum(dx),
            .groups = "drop") %>% 
  mutate(Location = case_when(Location == "United States of America" ~ "USA",
                              Location == "Republic of Moldova" ~ "Moldova",
                              Location == "Republic of Korea" ~ "South Korea",
                              TRUE ~ Location),
         Sex = recode(Sex,
                      "Female" = "f",
                      "Male" = "m", 
                      "Total" = "b")) %>% 
  select(Country = Location, Sex, e0, edagger)

C19_use <- readRDS("Data/C19_use.rds")
c19countries <- C19_use$Country %>% unique()
wpp_locs <- WPP_calcs$Location %>% unique()

c19countries[!c19countries %in% wpp_locs]


# getHMDcountries()
get_these <- c("GBRTENW", "GBR_SCO", "GBR_NIR")
get_names <- c("England and Wales", "Scotland","Northern Ireland")
names(get_names) <- get_these
HMD_calcs <- tibble(Country = character(),
                    Sex = character(),
                    e0 = double(),
                    edagger = double(),
                    b = double())
for (xyz in get_these){
  for (s in c("f","m","b")){
    this_chunk<- readHMDweb(xyz,paste0(s,"ltper_1x1"),us,pw)
    HMD_calcs <- 
    this_chunk %>% 
      dplyr::filter(Year == 2018) %>% 
      mutate(ex_mid = approx(x = Age, y = ex, xout = (Age + ax), rule = 2)$y) %>% 
      summarize(e0 = ex[Age == 0],
                edagger = sum(dx*ex_mid) / sum(dx),
                b = get_b(chunk = .data)) %>% 
      mutate(Sex = s,
             Country = get_names[xyz]) %>% 
      bind_rows(HMD_calcs) 
  }
}

Covariates <- bind_rows(WPP_calcs, HMD_calcs)

saveRDS(Covariates, file = "Data/Covariates.rds")

# -------------------------------------------------------------
# get Gompertz b parameter for each country


WPP <- read_csv("Data/WPP2019_Life_Table_Medium.csv")
WPP$MidPeriod %>% unique()
WPP <- 
  WPP %>% 
  filter(MidPeriod == 2018)


  
get_rates_of_aging <- function(chunk){
  abc <- MortalityLaw(x = chunk$AgeGrpStart, 
               Dx = chunk$dx, 
               Ex = chunk$Lx,
               opt.method ="poissonL",
               law = "makeham", 
               fit.this.x = seq(30,85,by=5))$coefficients
  a5 <- seq(30,85,by=5)
  mx <- 
  chunk %>% 
    dplyr::filter(AgeGrpStart >= 30, AgeGrpStart <= 85) %>% 
    dplyr::pull(mx)
   ftd <- loess(log(mx) ~ a5)$fitted 
   
   rout <- emp.der(x = a5,y = ftd)
   
  tibble(A = abc[1], B = abc[2], C = abc[3], r30 = rout[1], r35 = rout[2], r40= rout[3],
         r45 = rout[4], r50= rout[5],
         r55 = rout[6], r60= rout[7],
         r65 = rout[8], r70= rout[9],
         r75 = rout[10], r80= rout[11],
         r85= rout[12])
}





adds_b <-
  WPP %>% 
  group_by(LocID, Location, Sex) %>% 
  do(get_rates_of_aging(chunk = .data)) %>% 
  ungroup() %>% 
  mutate(Location = case_when(Location == "United States of America" ~ "USA",
                              Location == "Republic of Moldova" ~ "Moldova",
                              Location == "Republic of Korea" ~ "South Korea",
                              TRUE ~ Location),
         Sex = recode(Sex,
                      "Female" = "f",
                      "Male" = "m", 
                      "Total" = "b")) %>% 
  rename(Country = Location)


Covariates <- readRDS("Data/Covariates.rds")
Covariates <- left_join(Covariates, adds_b) %>% 
  select(-LocID)

saveRDS(Covariates, file = "Data/WPP_Covariates.rds")





library(data.table, readxl, httr)

url1 <- 'https://population.un.org/wpp/Download/Files/4_Metadata/WPP2019_F02_METADATA.XLSX'

GET(url1, write_disk(myfile <- tempfile(fileext = ".xlsx")))



wpp_metadata_Overall_Mortality <- data.table(read_excel(path= myfile, sheet = "Overall_Mortality", col_names = TRUE, col_types = NULL, na = "", skip = 16))

wpp_metadata_Overall_Mortality <- unique(wpp_metadata_Overall_Mortality[, .(LocID, Location, Age_Specific_Type_of_Pattern, Age_Specific_Type_of_Model, Age_Specific_Type_of_Inputs, Derived_from_Child_Mortality, Use_Adult_Mortality, HIV_AIDS_Mortality_Impact)])
