library(tidyverse)
library(readr)
library(countrycode)
WPP_hist <- read_csv("https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/CSV_FILES/WPP2019_PopulationBySingleAgeSex_1950-2019.csv")
WPPproj <- read_csv("https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/CSV_FILES/WPP2019_PopulationBySingleAgeSex_2020-2100.csv")
  # WPP <- read_csv("Data/WPP2019_PopulationBySingleAgeSex_1950-2019.csv") 

countries <- read_csv("Output/summary_selected_sources_by_country.csv") %>% 
  select(Country, Code) %>% 
  distinct()

colnames(WPP_hist)
colnames(WPPproj)

YrMin <- 2015
YrMax <- 2021

WPP <- bind_rows(WPP_hist,WPPproj)


locids <- WPP$LocID %>% unique()
sum(is.na( countrycode(locids, origin = "unpd", destination = "iso3c")))

offsets <-
  WPP %>% 
  dplyr::filter(between(Time, YrMin, YrMax)) %>% 
  pivot_longer(PopMale:PopTotal, names_to = "Sex", values_to = "Population") %>% 
  mutate(Sex = recode(Sex, 
                      "PopMale" = "m",
                      "PopFemale" = "f",
                      "PopTotal" = "t"),
         Region = "All") %>% 
  mutate(Population = Population * 1000,
         Code = suppressWarnings(countrycode(LocID, origin = "un", destination = "iso3c")),
         Code = ifelse(Location == "China, Taiwan Province of China","TWN",Code)) %>% 
  dplyr::filter(!is.na(Code)) %>% 
  select(Year = Time, Age = AgeGrpStart, Population, Code, Sex) %>% 
  arrange(Code, Year, Sex, Age) %>% 
  left_join(countries) %>% 
  dplyr::filter(!is.na(Country)) %>% 
  select(Country,Code,Year,Sex,Age,Population)

countries$Country
denom <-offsets$Country %>% unique()

countries$Country[!countries$Country%in% denom]

write_csv(offsets, file = "Output/offsets.csv")


# countries$Code[!countries$Code %in% offsets$Code]

# library(HMDHFDplus)
# get_these <- c("GBRTENW", "GBR_SCO", "GBR_NIR", "TWN")
# get_names <- c("England and Wales", "Scotland","Northern Ireland", "Taiwan")
# fir (i in 1:length(get_these)){
#   X <- readHMDweb(get_these[i],"Population",us,pw)
#   X$Country <-
# }


