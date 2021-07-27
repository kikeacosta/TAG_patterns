library(here)
source(here("R", "00_functions.R"))

# loading data
who <- read_csv("Output/who.csv")
unpd <- read_csv("Output/unpd.csv")
stmf <- read_csv("Output/stmf.csv") 
eurs <- read_csv("Output/eurs.csv") 
bra <- read_csv("Output/brazil.csv") 
mex <- read_csv("Output/mexico.csv") 
per <- read_csv("Output/peru.csv") 
zaf <- read_csv("Output/south_africa.csv")
col <- read_csv("Output/colombia.csv")
chl <- read_csv("Output/chile.csv")
ecu <- read_csv("Output/ecuador.csv",
                 col_types = "ccdccdc")

# adjusting sources for "country_public" data
latam <- 
  bind_rows(bra, 
            mex,
            per,
            ecu) %>% 
  mutate(Source = "country_public")


# excluding countries with missing ages from UNPD 
exc_unpd <- c("Costa Rica", "Oman", "Bermuda")
unpd2 <- 
  unpd %>% 
  filter(!Country %in% exc_unpd)


# Adjust for unknown ages and sex
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# WHO data already solved by Tim
# UNPD does not have missing data
# solving for stmf, eurostat, and brazil

db <- 
  bind_rows(stmf,
            eurs,
            latam,
            zaf,
            col,
            chl)

# imputing unknown ages and sexes
db2 <- 
  db %>% 
  group_by(Country, Sex, Year, Source) %>% 
  do(rescale_age(chunk = .data)) %>% 
  ungroup() %>% 
  group_by(Country, Age, Year, Source) %>% 
  do(rescale_sex(chunk = .data)) %>% 
  ungroup() %>% 
  mutate(Age = as.double(Age)) %>% 
  arrange(Source, Country, Year, Sex, Age)

# Adding data in 2020 for the UK
# grouping deaths for the UK in 2020 from the STMF
# 5-year age groups and 90+

uk2020 <- 
  db2 %>% 
  filter(Code %in% c("GBR_SCO", "GBR_NIR", "GBRTENW"),
         Year == 2020) %>% 
  mutate(Age = case_when(Age > 90 ~ 90, 
                         Age < 5 ~ 0,
                         TRUE ~ Age)) %>% 
  group_by(Year, Sex, Age) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Country = "United Kingdom", 
         Code = "GBR",
         Source = "stmf_adjusted")

# merge all data together
db3 <- 
  bind_rows(db2, 
            who,
            unpd2,
            uk2020) %>% 
  mutate(Deaths = ifelse(is.na(Deaths), 0, Deaths)) %>% 
  filter(!Code %in% c("GBR_SCO", "GBR_NIR", "GBRTENW"))


# create a table with data availability by population, including age groups and years

unique(db3$Age)
unique(db3$Sex)
unique(db3$Country) %>% sort()


# variable characteristics by country and source 
summ1 <- 
  db3 %>% 
  group_by(Country, Code, Year, Source, Sex) %>% 
  summarise(Age_groups = n(),
            Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  group_by(Country, Code, Year, Source, Age_groups) %>% 
  summarise(Sex_groups = n(),
            Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  group_by(Country, Code, Source, Age_groups, Sex_groups) %>% 
  summarise(Period = paste(min(Year), max(Year), sep = "-"),
         Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  group_by(Country) %>% 
  mutate(n_sources = n()) %>% 
  ungroup()

# variable configuration by source
summ2 <- 
  summ1 %>% 
  select(Country, Code, Source, n_sources) %>% 
  unique()


# selection of best source by country
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# based on data previous to 2020
# 

source_selec <- 
  db3 %>% 
  filter(Year <= 2019, Sex == "t") %>% 
  group_by(Country, Source, Year) %>% 
  summarise(Deaths = sum(Deaths),
            age_groups = n()) %>% 
  ungroup() %>% 
  group_by(Country) %>% 
  mutate(prop = age_groups / max(age_groups)) %>% 
  ungroup() %>% 
  filter(prop > 0.7) %>% 
  group_by(Country, Year) %>% 
  mutate(prop_d = Deaths / max(Deaths)) %>% 
  arrange(Country, Year, -Deaths) %>% 
  filter(prop_d > 0.98) %>% 
  ungroup() %>% 
  group_by(Country, Source) %>% 
  mutate(years = n()) %>% 
  ungroup() %>% 
  group_by(Country) %>% 
  filter(years == max(years)) %>% 
  ungroup() %>% 
  group_by(Country, Source) %>% 
  mutate(min_ages = min(age_groups)) %>% 
  group_by(Country) %>% 
  filter(min_ages == max(min_ages)) %>% 
  ungroup() %>% 
  group_by(Country, Source) %>% 
  mutate(all_deaths = sum(Deaths)) %>% 
  group_by(Country) %>% 
  filter(all_deaths == max(all_deaths)) %>% 
  ungroup() %>% 
  group_by(Country, Year) %>% 
  arrange(Country, Year, Source) %>% 
  mutate(order = 1:n()) %>% 
  ungroup() %>% 
  filter(order == 1) %>% 
  ungroup() %>% 
  select(Country, Source, Year, age_groups, min_ages)

best_source_by_country <- 
  source_selec %>% 
  select(Country, Source) %>% 
  unique()

# filtering deaths in best sources
db4 <- 
  db3 %>% 
  left_join(best_source_by_country %>% 
              mutate(best_source = 1)) %>% 
  filter(best_source == 1) %>% 
  select(-best_source)



# harmonization of ages before 2020
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# using the minimum common denominator

# selecting countries that need harmonization
cts_to_harmonize <- 
  source_selec %>% 
  select(Country, age_groups) %>% 
  unique() %>% 
  group_by(Country) %>% 
  summarise(age_configs = n()) %>% 
  filter(age_configs > 1) %>% 
  pull(Country)

# selecting the reference period 
age_groups_period <- 
  source_selec %>% 
  filter(Country %in% cts_to_harmonize) %>% 
  mutate(is_min_ages = ifelse(age_groups == min_ages, 1, 0)) %>% 
  filter(is_min_ages == 1) %>% 
  group_by(Country) %>% 
  mutate(is_reference = 1:n()) %>% 
  ungroup() %>% 
  filter(is_reference == 1) %>% 
  select(Country, Source, Year, is_reference) 
  
# pulling reference age groups
age_groups <- 
  db4 %>%
  left_join(age_groups_period) %>% 
  filter(Sex == "t",
         is_reference == 1) %>% 
  select(Country, Age)

harmonized_ages <- tibble()
for(ct in cts_to_harmonize){
  
  harmonized_ages <- 
    harmonized_ages %>% 
    bind_rows(assign_age_intervals(db4 , ct))
  
}

harmonized_ages2 <- 
  harmonized_ages %>% 
  select(-Age) %>% 
  rename(Age = Age_int) %>% 
  group_by(Country, Code, Source, Year, Sex, Age) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup()

db5 <- 
  db4 %>% 
  filter(!Country %in% cts_to_harmonize) %>% 
  bind_rows(harmonized_ages2)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# selecting best source for each country
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# criteria:
# 1. Three sex groups
# 2. More age groups
# 3. More observed periods
# 4. More deaths (only apply for Germany, so far)

# only Greece has two equal best sources: stmf and eurostat, lets choose stmf :)


best_source_year <- 
  db3 %>% 
  group_by(Country, Code, Year, Source, Sex) %>% 
  summarise(Age_groups = n(),
            Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  group_by(Country, Code, Year, Source, Age_groups) %>% 
  summarise(Sex_groups = n(),
            Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(priority = case_when(Source == "stmf" ~ 1,
                              Source == "eurs" ~ 2,
                              Source == "who" ~ 3,
                              TRUE ~ 4)) %>% 
  group_by(Country, Code, Source, Age_groups, Sex_groups) %>% 
  mutate(Period = paste(min(Year), max(Year), sep = "-"),
         Deaths = sum(Deaths), 
         Periods = max(Year) - min(Year) + 1) %>% 
  group_by(Country, Year) %>% 
  filter(Age_groups == max(Age_groups)) %>% 
  filter(Sex_groups == max(Sex_groups)) 



best_source_year <- 
  db3 %>% 
  group_by(Country, Code, Year, Source, Sex) %>% 
  summarise(Age_groups = n(),
            Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  group_by(Country, Code, Year, Source, Age_groups) %>% 
  summarise(Sex_groups = n(),
            Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(priority = case_when(Source == "stmf" ~ 1,
                              Source == "eurs" ~ 2,
                              Source == "who" ~ 3,
                              TRUE ~ 4)) %>% 
  group_by(Country, Code, Source, Age_groups, Sex_groups) %>% 
  mutate(Period = paste(min(Year), max(Year), sep = "-"),
         Deaths = sum(Deaths), 
         Periods = max(Year) - min(Year) + 1) %>% 
  group_by(Country, Year) %>% 
  filter(Age_groups == max(Age_groups)) %>% 
  filter(Sex_groups == max(Sex_groups)) %>% 
  filter(Periods == max(Periods)) %>% 
  filter(Deaths == max(Deaths)) %>% 
  filter(priority == min(priority)) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  group_by(Country, Code, Source, Age_groups) %>% 
  mutate(max_year = max(Year),
         min_year = min(Year)) %>% 
  ungroup() %>% 
  # filtering countries with at least data for 2019 and 2020  
  group_by(Country) %>% 
  mutate(max_year = max(max_year),
         min_year = min(min_year)) %>% 
  filter(max_year == 2020 & min_year <= 2019) %>% 
  mutate(Age_groups_2019 = Age_groups[Year == 2019],
         ratio_ages = Age_groups_2019 / Age_groups) %>% 
  ungroup() %>% 
  # excluding years with less than half of the age groups in 2019
  filter(ratio_ages < 2) %>% 
  select(Country, Code, Year, Source, Age_groups, Sex_groups) %>% 
  mutate(Best = 1)

best_source_summ <- 
  best_source_year %>% 
  group_by(Country, Code, Source, Age_groups, Sex_groups) %>% 
  mutate(Period = paste(min(Year), max(Year), sep = "-")) %>% 
  ungroup() %>% 
  select(Country, Code, Source, Age_groups, Sex_groups, Period) %>% 
  unique() %>% 
  group_by(Country) %>% 
  mutate(n_sources = n()) %>% 
  ungroup() %>% 
  arrange(-n_sources)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# filtering best sources by country
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

db_best <- 
  db3 %>% 
  left_join(best_source_year %>% 
              select(Code, Source, Year, Best)) %>% 
  filter(Best == 1) %>% 
  select(-Best)


# summary of selected sources by country 
available <- 
  db_best %>% 
  group_by(Country, Code, Year, Source, Sex) %>% 
  summarise(Age_groups = n(),
            Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  group_by(Country, Code, Year, Source, Age_groups) %>% 
  summarise(Sex_groups = n(),
            Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  group_by(Country, Code, Source, Age_groups, Sex_groups) %>% 
  summarise(Period = paste(min(Year), max(Year), sep = "-"),
            Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  group_by(Country) %>% 
  mutate(n_sources = n()) %>% 
  ungroup() 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# saving 3 output data objects
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1) summary of all available sources by country 
write_csv(summ1, "Output/summary_all_sources_by_country.csv")

# 2) summary of selected sources by country
write_csv(available, "Output/summary_selected_sources_by_country.csv")

# 3) output mortality data based on selected sources
write_csv(db_best, "Output/annual_deaths_countries_selected_sources.csv")

