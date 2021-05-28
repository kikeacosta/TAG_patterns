source(here("R", "00_functions.R"))

# loading data
who <- read_csv("Output/who.csv")
stmf <- read_csv("Output/stmf.csv") 
eurs <- read_csv("Output/eurs.csv") 
brazil <- read_csv("Output/brazil.csv") 
mexico <- read_csv("Output/mexico.csv") 
peru <- read_csv("Output/peru.csv") 

# solving issue with age in Eurostat data
unique(eurs$Age)
unique(peru$Age)
unique(mexico$Age)

eurs2 <- 
  eurs %>% 
  mutate(Age = recode(Age,
                      "Y50-54" = "50",
                      "Y55-59" = "55"))

brazil2 <- 
  brazil %>% 
  mutate(Source = "brazil")

unique(eurs2$Age)

# Adjust for unknown ages and sex
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# WHO data already solved by Tim
# solving for stmf, eurostat, and brazil

db <- 
  bind_rows(stmf,
            eurs2,
            brazil2,
            peru,
            mexico)

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
         Source = "stmf modified")

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
  summarise(Ymin = min(Year),
         Ymax = max(Year),
         Deaths = sum(Deaths)) %>% 
  ungroup()
  
# variable configuration by source
summ2 <- 
  summ1 %>% 
  group_by(Country, Code, Source) %>% 
  summarise(n = n())

# best source by country
# ~~~~~~~~~~~~~~~~~~~~~~

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
  mutate(Period = paste(min(Year), max(Year), sep = "-")) %>% 
  ungroup() %>% 
  mutate(Best = 1)

best_source_summ <- 
  best_source_year %>% 
  select(Country, Code, Source, Period, Age_groups) %>% 
  unique() %>% 
  mutate(Best = 1)
  
# filtering best sources in each country
db_best <- 
  db3 %>% 
  left_join(best_source_year %>% 
              select(Code, Source, Year, Best)) %>% 
  filter(Best == 1)

unique(db_best$Country)

# data available
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
  summarise(Ymin = min(Year),
            Ymax = max(Year),
            Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  group_by(Country) %>% 
  mutate(n = n())

write_csv(db_best, "Output/annual_deaths_countries_best_source.csv")
write_csv(summ1, "Output/summary_sources_by_country.csv")
write_csv(best_source_summ, "Output/country_list_best_source.csv")
write_csv(available, "Output/annual_deaths_countries_best_source.csv")

