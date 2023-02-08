# Note: repo name != package name!
# remotes::install_github("eshom/covid-age-data")
library(covidAgeData) # for COVerAGE-DB
library(tidyverse)
library(lubridate)
library(osfr)


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
  mutate(Age = as.integer(Age),
         Deaths = Deaths / sum(Deaths) * total_deaths,
         C19f = f / (f + m) * Deaths,
         C19m = m / (f + m) * Deaths) %>% 
  dplyr::filter(!is.na(f)) %>% 
  dplyr::select(-Deaths, -f, -m, -total_deaths, -DateDiff)
  

offsets <- readRDS("Data/offsets.rds")

ctries <- C19_use$Country %>% unique()
ctries[!ctries%in%offsets$Country]

all(C19_use$Country %in% offsets$Country)

cdb_use2 <- list()
for (i in ctries){
  off <- dplyr::filter(offsets, Country == i)
  cdb <- dplyr::filter(C19_use, Country == i)
  cdb <- mutate(cdb, ExpF = NA, ExpM = NA)
  if (nrow(off) > 0){
    maxn <- min(nrow(cdb),105)
    for (j in 1:maxn){
      a <- cdb$Age[j]
      int <-  cdb$AgeInt[j]
      ages <- a:(a+int-1)
      cdb$ExpF[j] <- off %>% dplyr::filter(Age %in% ages) %>% dplyr::pull(f) %>% sum()
      cdb$ExpM[j] <- off %>% dplyr::filter(Age %in% ages) %>% dplyr::pull(m) %>% sum()
    }
  }
  cdb_use2[[i]] <- cdb
}

C19_use_sex <- do.call("rbind", cdb_use2)

C19_use_sex %>% 
  select(Country, Age, AgeInt, C19F = C19f, C19M = C19m, ExpF, ExpM) %>% 
  mutate(C19F = ifelse(is.nan(C19F),0,C19F),
         C19M = ifelse(is.nan(C19M),0,C19M)) %>% 
  ungroup()

saveRDS(C19_use_sex,file = "Data/C19_use_sex.rds")
