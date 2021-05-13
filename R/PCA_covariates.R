# Covariates

library(tidyverse)
library(readr)

WPP <- read_csv("Data/WPP2019_Life_Table_Medium.csv")
WPP$MidPeriod %>% unique()
WPP <- 
  WPP %>% 
  filter(MidPeriod == 2018)

gc()
1:10 + lag(1:10)
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

library(HMDHFDplus)
getHMDcountries()
get_these <- c("GBRTENW", "GBR_SCO", "GBR_NIR")
get_names <- c("England and Wales", "Scotland","Northern Ireland")
names(get_names) <- get_these
HMD_calcs <- tibble(Country = character(),
                    Sex = character(),
                    e0 = double(),
                    edagger = double())
for (xyz in get_these){
  for (s in c("f","m","b")){
    this_chunk<- readHMDweb(xyz,paste0(s,"ltper_1x1"),us,pw)
    HMD_calcs <- 
    this_chunk %>% 
      dplyr::filter(Year == 2018) %>% 
      mutate(ex_mid = approx(x = Age, y = ex, xout = (Age + ax), rule = 2)$y) %>% 
      summarize(e0 = ex[Age == 0],
                edagger = sum(dx*ex_mid) / sum(dx)) %>% 
      mutate(Sex = s,
             Country = get_names[xyz]) %>% 
      bind_rows(HMD_calcs) 
  }
}

Covariates <- bind_rows(WPP_calcs, HMD_calcs)

saveRDS(Covariates, file = "Data/Covariates.rds")

