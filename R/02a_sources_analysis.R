library(here)
source(here("R", "00_functions.R"))

# loading data
who <- read_csv("Output/who.csv")
stmf <- read_csv("Output/stmf.csv") 
eurs <- read_csv("Output/eurs.csv") 
brazil <- read_csv("Output/brazil.csv") 
mexico <- read_csv("Output/mexico.csv") 
peru <- read_csv("Output/peru.csv") 
safr <- read_csv("Output/south_africa.csv")


# solving issue with age in Eurostat data
unique(safr$Age)
unique(safr$Sex)
unique(safr$Year)

eurs2 <- 
  eurs %>% 
  mutate(Age = recode(Age,
                      "Y50-54" = "50",
                      "Y55-59" = "55"))

brazil2 <- 
  brazil %>% 
  mutate(Source = "brazil")

safr2 <- 
  safr %>% 
  tidyr::complete(Country, Code, Year, Sex, Age, fill = list(Deaths = 0)) %>% 
  group_by(Country, Code, Year, Sex) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Age = "TOT") %>% 
  bind_rows(safr %>% 
              mutate(Age = as.character(Age))) %>% 
  arrange(Year, Sex) %>% 
  mutate(Source = "south_africa")

unique(safr2$Age)

# Adjust for unknown ages and sex
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# WHO data already solved by Tim
# solving for stmf, eurostat, and brazil

db <- 
  bind_rows(stmf,
            eurs2,
            brazil2,
            peru,
            mexico,
            safr2)

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
         Source = "stmf_modified")

# merge all data together
db3 <- 
  bind_rows(db2, 
            who,
            uk2020) %>% 
  mutate(Deaths = ifelse(is.na(Deaths), 0, Deaths))


# create a table with data availability by population, including age groups and years

unique(db3$Age)
unique(db3$Sex)
unique(db3$Country) %>% sort()

eurs_test <- 
  db3 %>% 
  group_by() %>% 
  summarise(Deaths = sum(Deaths))


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
  group_by(Country, Code, Source, Age_groups, Sex_groups) %>% 
  mutate(Period = paste(min(Year), max(Year), sep = "-"),
         Deaths = sum(Deaths), 
         Periods = max(Year) - min(Year) + 1) %>% 
  group_by(Country, Year) %>% 
  filter(Sex_groups == max(Sex_groups),
         Age_groups == max(Age_groups)) %>% 
  group_by(Country, Year) %>% 
  filter(Periods == max(Periods)) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n == 1 | Source == "stmf") %>% 
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

