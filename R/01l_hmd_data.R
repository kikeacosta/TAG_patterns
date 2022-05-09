
library(HMDHFDplus)
library(tidyverse)
library(countrycode)

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

hmd2 <- 
  hmd %>%
  # only countries with data in 2020
  group_by(Code) %>% 
  filter(max(Year) >= 2020) %>% 
  ungroup() %>% 
  select(-OpenInterval) %>% 
  gather(Female, Male, Total, key = Sex, value = Deaths) %>% 
  mutate(Sex = recode(Sex,
                      "Female" = "f",
                      "Male" = "m",
                      "Total" = "t"),
         Country = countrycode(Code, origin = "iso3c",
                               destination = "country.name"),
         Age = Age %>% as.character()) %>% 
  mutate(Source = "hmd")

hmd_out <- 
  hmd2 %>% 
  group_by(Year, Sex, Code, Country, Source) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Age = "TOT") %>% 
  bind_rows(hmd2) %>% 
  arrange(Country, Sex, Year, suppressWarnings(as.integer(Age))) %>% 
  select(Country, Code, Year, Sex, Age, Deaths, Source)

write_csv(hmd_out, "Output/hmd.csv")


