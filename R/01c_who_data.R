library(readxl)
source("R/00_functions.R")
rm(list=ls())
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Updating WHO files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# downloading last version of WHO files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# these come from here: https://www.who.int/data/data-collection-tools/who-mortality-database
icd_base_url <-"https://cdn.who.int/media/docs/default-source/world-health-data-platform/mortality-raw-data/"
icd_files    <- c("mort_country_codes.zip","morticd10_part1.zip","morticd10_part2.zip",
                  "morticd10_part3.zip","morticd10_part4.zip","morticd10_part5.zip")
for (i in 1:length(icd_files)){
  url_i   <- paste0(icd_base_url,icd_files[i])
  local_i <- file.path("Data", "WHO", icd_files[i])
  download.file(url_i, destfile = local_i, overwrite = TRUE)
}

# a lookup table to match country names to codes
ctry_names <- read_csv(file.path("Data", "WHO", "mort_country_codes.zip")) %>% 
  rename(Country = country)

icd_all    <- list()
icd_files2 <- icd_files[-1]
#ICD download each of the 5 files
for (i in 1:length(icd_files2)){
  icd_i <- 
    read_csv(file.path("Data", "WHO", icd_files2[i]),
             col_types = cols(Admin1 = col_character(),SubDiv = col_character(),
                              List = col_character(), Cause = col_character(), 
                              Frmat = col_character(), IM_Frmat = col_character(),
                              .default = col_double())) %>% 
    left_join(ctry_names, by = "Country") %>% 
    dplyr::filter(Sex %in% c(1,2)) 
  
  icd_all[[i]] <- icd_i
}

# stick together
# ~~~~~~~~~~~~~~
icd_all <- 
  bind_rows(icd_all) %>% 
  select(name, everything())

# saving a consolidated file with all WHO data
write_rds(icd_all, "Data/WHO/who_raw.rds")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls()); gc()
source("R/00_functions.R")

db <- 
  read_rds("Data/WHO/who_raw.rds")

# selecting: 
# countries with data in 2020; data by age; all-cause mortality; data since 2015
db20 <- 
  db %>% 
  select(-Country) %>% 
  rename(Country = name) %>% 
  mutate(Year = Year %>% as.double()) %>% 
  group_by(Country) %>% 
  filter(max(Year) >= 2020 &  Frmat != "09") %>% 
  ungroup() %>% 
  filter(Cause %in% c("1000", "AAA"),
         Year >= 2010) %>% 
  select(-Admin1, -SubDiv, -Cause, -List) %>% 
  arrange(Country, Year)

unique(db20$Frmat)

db20_2 <- 
  db20 %>% 
  select(Country, Year, Sex, starts_with("Deaths")) %>% 
  gather(starts_with("Deaths"), key = Age, value = Deaths) %>% 
  filter(Age %in% paste0("Deaths", 1:25)) %>% 
  mutate(Age = recode(Age,
                      'Deaths1' = "TOT",
                      'Deaths2' = "0",
                      'Deaths3' = "1",
                      'Deaths4' = "1",
                      'Deaths5' = "1",
                      'Deaths6' = "1",
                      'Deaths7' = "5",
                      'Deaths8' = "10",
                      'Deaths9' = "15",
                      'Deaths10' = "20",
                      'Deaths11' = "25",
                      'Deaths12' = "30",
                      'Deaths13' = "35",
                      'Deaths14' = "40",
                      'Deaths15' = "45",
                      'Deaths16' = "50",
                      'Deaths17' = "55",
                      'Deaths18' = "60",
                      'Deaths19' = "65",
                      'Deaths20' = "70",
                      'Deaths21' = "75",
                      'Deaths22' = "80",
                      'Deaths23' = "85",
                      'Deaths24' = "90",
                      'Deaths25' = "95"),
         Deaths = Deaths %>% as.double(),
         Sex = recode(Sex,
                      "1" = "m",
                      "2" = "f")) %>%
  drop_na(Deaths) %>% 
  mutate(Code = countrycode(Country, origin = 'country.name', destination = 'iso3c'),
         Country = recode(Country,
                          "United Kingdom, England and Wales" = "England and Wales",
                          "United Kingdom, Northern Ireland" = "Northern Ireland",
                          "United Kingdom, Scotland" = "Scotland",
                          "United States of America" = "USA",
                          "Republic of Korea" = "South Korea",
                          "Czech Republic" = "Czechia"),
         Code = case_when(Country == "England and Wales" ~ "GBR-ENW",
                          Country == "Scotland" ~ "GBR-SCO",
                          Country == "Northern Ireland" ~ "GBR-NIR",
                          TRUE ~ Code)) %>% 
  group_by(Country, Code, Year, Sex, Age) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  arrange(Code, Year, Sex, Age) %>% 
  group_by(Country, Year, Sex) %>% 
  mutate(Source = "who_mort_db")

unique(db20_2$Code)
unique(db20_2$Age)
unique(db20_2$Sex)

# adding total sex
db20_3 <- 
  db20_2 %>% 
  group_by(Country, Code, Year, Age, Source) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  mutate(Sex = "t") %>% 
  bind_rows(db20_2) %>% 
  select(Country, Code, Year, Sex, Age, Deaths, Source) %>% 
  arrange(Country, Code, Year, Sex, suppressWarnings(as.integer(Age)))

db20_4 <- 
  db20_3 %>% 
  group_by(Country, Sex, Year) %>% 
  do(rescale_age(chunk = .data)) %>% 
  ungroup() %>%
  group_by(Country, Age, Year) %>%
  do(rescale_sex(chunk = .data)) %>% 
  ungroup() %>% 
  mutate(Age = Age %>% as.double()) %>% 
  arrange(Code, Year, Sex, Age) %>% 
  group_by(Country, Year, Sex) %>% 
  mutate(age_spn = ifelse(Age == max(Age), -1, lead(Age) - Age)) %>% 
  ungroup()

write_csv(db20_4, "data_inter/who.csv")
