rm(list=ls())
source("R/00_functions.R")
library(HMDHFDplus)

# exclude populations that are subsets of other populations
cts_exclude <- c("GBRCENW", "GBRTENW", "GBR_SCO", "GBR_NIR",
                 "DEUTE", "DEUTW", 
                 "FRACNP", 
                 "NZL_NM", "NZL_MA")

cds_hmd <- getHMDcountries() %>% pull(CNTRY)
cds_hmd2 <- cds_hmd[!cds_hmd %in% cts_exclude]


# Extract deaths and exposures between 2010 and 2021 from the HMD
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HMD user and password
# usethis::edit_r_environ()
# hmd_us="acosta@demogr.mpg.de"
# hmd_pw="Secreto_1"

# getting HMD username and password from the R environment
hmd_us <- Sys.getenv("hmd_us")
hmd_pw <- Sys.getenv("hmd_pw")

# identifying those with data for 2020
hmd <- tibble()
for(ct in cds_hmd2){
  cat(paste0(ct, "\n"))
  chunk_d <- 
    readHMDweb(ct, "Deaths_1x1", hmd_us, hmd_pw) %>%
    filter(Year >= 2010) %>%
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
  filter(max(Year) >= 2019) %>% 
  ungroup() %>% 
  select(-OpenInterval) %>% 
  gather(Female, Male, Total, key = Sex, value = Deaths) %>% 
  # filter(!Code %in% cods_exc) %>% 
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
         Country = recode(Country,
                          "United States" = "USA",
                          "Hong Kong SAR China" = "Hong Kong"),
         # Age = Age %>% as.character(),
         Age = ifelse(Age > 100, 100, Age)) %>% 
  group_by(Code, Country, Year, Sex, Age) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Source = "hmd") %>% 
  group_by(Country, Year, Sex) %>% 
  mutate(age_spn = ifelse(Age == max(Age), -1, lead(Age) - Age)) %>% 
  ungroup()

unique(hmd2$Code)
unique(hmd2$Country)

write_csv(hmd2, "data_inter/hmd.csv")




