rm(list=ls())
source("R/00_functions.R")
library(HMDHFDplus)

cds_hmd <- getHMDcountries()

# cds_hmd <- c("AUS", "AUT")
# Extract deaths and exposures in 2020 from the HMD

# HMD user and password
# ~~~~~~~~~~~~~~~~~~~~
hmd_us <- "kikepaila@gmail.com"
hmd_pw <- "secreto"

# identifying those with data for 2020
hmd <- tibble()
for(ct in cds_hmd){
  chunk_d <- 
    readHMDweb(ct, "Deaths_1x1", hmd_us, hmd_pw) %>%
    filter(Year >= 2015) %>%
    as_tibble() %>%
    mutate(Code = ct)

  hmd <- 
    hmd %>%
    bind_rows(chunk_d)
}

cods_exc <- c("GBRTENW", "GBRCENW", "GBR_SCO", "GBR_NIR",
              "DEUTE", "DEUTW", "FRACNP")

hmd2 <- 
  hmd %>%
  # only countries with data in 2020
  group_by(Code) %>% 
  filter(max(Year) >= 2020) %>% 
  ungroup() %>% 
  select(-OpenInterval) %>% 
  gather(Female, Male, Total, key = Sex, value = Deaths) %>% 
  filter(!Code %in% cods_exc) %>% 
  mutate(Sex = recode(Sex,
                      "Female" = "f",
                      "Male" = "m",
                      "Total" = "t"),
         Code = case_when(Code == "GBR_NP" ~ "GBR",
                          Code == "DEUTNP" ~ "DEU",
                          Code == "NZL_NP" ~ "NZL",
                          Code == "FRATNP" ~ "FRA",
                          TRUE ~ Code),
         Country = countrycode(Code, origin = "iso3c",
                               destination = "country.name"),
         Age = Age %>% as.character()) %>% 
  mutate(Source = "hmd")

unique(hmd2$Code)
unique(hmd2$Country)

# adding total age
hmd3 <- 
  hmd2 %>% 
  group_by(Year, Sex, Code, Country, Source) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Age = "TOT") %>% 
  bind_rows(hmd2) %>% 
  arrange(Country, Sex, Year, suppressWarnings(as.integer(Age))) %>% 
  select(Country, Code, Year, Sex, Age, Deaths, Source)

# re-scaling age and sex
hmd4 <- 
  hmd3 %>% 
  group_by(Country, Sex, Year) %>% 
  do(rescale_age(chunk = .data)) %>% 
  ungroup() %>%
  group_by(Country, Age, Year) %>%
  do(rescale_sex(chunk = .data)) %>% 
  ungroup() %>% 
  mutate(Age = Age %>% as.double()) %>% 
  arrange(Code, Year, Sex, Age)

write_csv(hmd4, "Output/hmd.csv")


