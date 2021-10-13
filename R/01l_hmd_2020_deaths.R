
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
hmd_2020 <- tibble()
for(ct in cds_hmd){
  chunk_d <- readHMDweb(ct, "Deaths_1x1", hmd_us, hmd_pw) %>%
    filter(Year == 2020) %>%
    as_tibble() %>%
    mutate(Code = ct)

  hmd_2020 <- hmd_2020 %>%
    bind_rows(chunk_d)
}

cts_2020 <- unique(hmd_2020$Code)
hmd_long <- tibble()
for(ct in cts_2020){
  chunk_d <- readHMDweb(ct, "Deaths_1x1", hmd_us, hmd_pw) %>%
    filter(Year >= 2015) %>%
    as_tibble() %>%
    mutate(Code = ct)
  
  hmd_long <- hmd_long %>%
    bind_rows(chunk_d)
}

hmd_long2 <- 
  hmd_long %>% 
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

age_tot <- 
  hmd_long2 %>% 
  group_by(Year, Sex, Code, Country, Source) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Age = "TOT")

hmd_out <- 
  hmd_long2 %>% 
  bind_rows(age_tot) %>% 
  arrange(Country, Sex, Year, suppressWarnings(as.integer(Age)))

write_csv(hmd_out, "Output/hmd_2020.csv")



