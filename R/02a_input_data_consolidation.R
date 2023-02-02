rm(list=ls())
source("R/00_functions.R")

# loading pre-processed data
hmd  <- read_csv("data_inter/hmd.csv")
who  <- read_csv("data_inter/who.csv")
unpd <- read_csv("data_inter/unpd.csv")
stmf <- read_csv("data_inter/stmf.csv") 
eurs <- read_csv("data_inter/eurs.csv") 
bra <- read_csv("data_inter/brazil.csv") %>% mutate(Source = "direct")
per <- read_csv("data_inter/peru.csv") %>% mutate(Source = "direct")
mex <- read_csv("data_inter/mexico.csv") %>% mutate(Source = "direct")
zaf <- read_csv("data_inter/south_africa.csv") %>% mutate(Source = "direct")
irn <- read_csv("data_inter/iran.csv") %>% mutate(Source = "direct")
mse <- read_csv("data_inter/msemburi_tag.csv") 


cts_exclude <- tibble(Country = 
                        c("Scotland", 
                          "England and Wales", 
                          "Northern Ireland",
                          "Kenya",
                          "Puerto Rico",
                          "Anguilla"))

# unique(out$Country)

unique(unpd$Country)
# Several direct sources have been superseded by UNPD and others

# putting all together
all_in <- 
  bind_rows(hmd,
            who,
            unpd,
            stmf,
            eurs,
            bra,
            per,
            mex,
            zaf,
            irn) %>% 
  replace_na(list(Deaths = 0)) %>% 
  group_by(Source, Country) %>% 
  filter(max(Year) >= 2020 & min(Year) <= 2017) %>% 
  ungroup() %>% 
  mutate(Age = ifelse(Age >= 100, 100, Age)) %>% 
  group_by(Country, Code, Source, Year, Sex, Age) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup()


# excluding sources with insufficient data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sources to exclude because of insufficient periods for baseline (minimum 3 periods)
insuf_pers <- 
  all_in %>% 
  select(Country, Source, Year, Sex) %>% 
  unique() %>% 
  filter(Year < 2020) %>% 
  group_by(Country, Source, Sex) %>% 
  summarise(yrs = n()) %>% 
  ungroup() %>% 
  filter(yrs < 3)  %>% 
  select(-yrs)

# sources to exclude because of insufficient age groups (minimum 4 intervals)
insuf_ages <- 
  all_in %>% 
  group_by(Country, Source, Year, Sex) %>% 
  summarise(ages = n()) %>% 
  ungroup() %>% 
  filter(ages < 4) %>% 
  select(-ages)

# excluding cases with insufficient periods and insufficient ages
all_in2 <- 
  all_in %>% 
  anti_join(insuf_pers) %>% 
  anti_join(insuf_ages) %>% 
  arrange(Source, Country, Year, Sex, Age)

# harmonizing ages before 2020
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
min_open_age <- 
  all_in2 %>% 
  filter(Year < 2020) %>% 
  select(Source, Country, Year, Sex, Age) %>% 
  unique() %>% 
  group_by(Source, Country, Year, Sex) %>% 
  filter(Age == max(Age)) %>% 
  ungroup() %>% 
  group_by(Source, Country, Sex) %>% 
  filter(Age == min(Age)) %>% 
  ungroup() %>% 
  rename(min_open = Age) %>% 
  select(-Year) %>% 
  unique()

# Largest age groups before 2020
ref_ages <- 
  all_in2 %>% 
  filter(Year < 2020) %>% 
  # no ages surpassing the minimum open age interval
  left_join(min_open_age) %>% 
  filter(Age <= min_open) %>% 
  select(Source, Country, Year, Age) %>% 
  unique() %>% 
  # calculate age span in each interval
  group_by(Source, Country, Year) %>% 
  mutate(span = ifelse(Age == max(Age),
                       Inf,
                       lead(Age) - Age)) %>%
  ungroup() %>% 
  # identify and choose the largest age span in all periods for each age
  group_by(Source, Country, Age) %>% 
  mutate(max_span = max(span)) %>% 
  ungroup() %>% 
  select(Source, Country, Age, max_span) %>% 
  unique() %>% 
  group_by(Source, Country) %>%
  mutate(test = ifelse(Age == 0, 0, lag(Age + max_span))) %>% 
  ungroup() %>% 
  filter(Age == test) %>% 
  select(-test, -max_span)

# harmonizing ages within countries before 2020 
harmonized_ages <- 
  all_in2 %>%
  filter(Year < 2020) %>% 
  group_by(Source, Country, Code, Year, Sex) %>% 
  do(harmon_age_intervals(chunk = .data)) %>% 
  ungroup()

# together harmonized ages before 2020, and original age groups since 2020
all_in3 <- 
  bind_rows(harmonized_ages, 
            all_in2 %>% 
              filter(Year >= 2020)) %>% 
  arrange(Country, Source, Year, Sex, Age) 
  

# Summarizing content
# ~~~~~~~~~~~~~~~~~~~
# summarizing content by source
sum_all <- 
  sum_source(all_in3) %>% 
  unique() %>% 
  arrange(Country)

write_excel(sum_all)

# selection based on source
sel_source <- 
  sum_all %>% 
  group_by(Country) %>% 
  filter((any(Source == "direct") & Source == "direct") | all(Source != "direct")) %>%
  filter((any(Source == "hmd") & Source == "hmd") | all(Source != "hmd")) %>%
  filter((any(Source == "stmf") & Source == "stmf") | all(Source != "stmf")) %>%
  filter((any(Source == "eurs") & Source == "eurs") | all(Source != "eurs")) %>%
  filter((any(Source == "unpd") & Source == "unpd") | all(Source != "unpd")) %>%
  filter((any(Source == "who_mort_db") & Source == "who_mort_db") | all(Source != "who_mort_db")) %>%
  unique() %>% 
  mutate(n = n()) %>% 
  ungroup()

# selection based on data
sel_criteria <- 
  sum_all %>% 
  group_by(Country) %>% 
  filter(ages >= max(ages)) %>% 
  filter(years == max(years)) %>% 
  filter(sexs == max(sexs)) %>% 
  filter(Deaths == max(Deaths)) %>% 
  unique() %>% 
  filter(1:n() == n()) %>% 
  mutate(n = n()) %>% 
  ungroup()

sel_mixed <- 
  sum_all %>% 
  group_by(Country) %>% 
  filter((any(Source == "direct") & Source == "direct") | all(Source != "direct")) %>%
  filter((any(Source == "hmd") & Source == "hmd") | all(Source != "hmd")) %>%
  filter(ages == max(ages)) %>% 
  filter(years == max(years)) %>% 
  filter(sexs == max(sexs)) %>% 
  filter(Deaths == max(Deaths)) %>% 
  unique() 


# comparison between both selection methods
comp <- 
  sel_source %>% 
  select(Country, ages_s = ages, years_s = years, Source_s = Source) %>%
  left_join(sel_criteria %>% 
              select(Country, ages_c = ages, years_c = years, Source_c = Source)) %>% 
  mutate(same = ifelse(Source_s == Source_c, 1, 0),
         ages_ratio = ages_c / ages_s,
         years_ratio = years_c / years_s)
  

# selecting best source for each population
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# prefer selection based on source, unless content-based selection has either
# 20% more ages groups or periods
best_source <- 
  comp %>% 
  mutate(Source = ifelse(ages_ratio > 1.2 | years_ratio > 1.2, Source_c, Source_s)) %>% 
  select(Country, Source)

out <- 
  all_in3 %>% 
  inner_join(best_source) %>% 
  arrange(Country, Year, Sex, Age) %>% 
  mutate(Country = ifelse(Country == "Faeroe Islands", "Faroe Islands", Country)) %>% 
  anti_join(cts_exclude)

unique(out$Country)

# 3) output mortality data based on selected sources
write_csv(out, "data_inter/annual_deaths_countries_selected_sources.csv")

# summary of selected sources by country 
available <- 
  out %>% 
  group_by(Country, Code, Year, Source, Sex) %>% 
  summarise(Age_groups = n(),
            Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  group_by(Country, Code, Year, Source, Age_groups) %>% 
  summarise(Sex_groups = n(),
            Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  group_by(Country, Code, Source) %>% 
  summarise(Period = paste(min(Year), max(Year), sep = "-"),
            Deaths = sum(Deaths),
            # putting together all the age and sex groups configurations between 2015 and 2021
            Age_groups = unique(Age_groups) %>% paste(collapse = ", "),
            Sex_groups = unique(Sex_groups) %>% paste(collapse = ", ")) %>% 
  ungroup() 

write_excel(available)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Harmonizing age groups in all periods
# =====================================
d1 <- 
  read_csv("data_inter/annual_deaths_countries_selected_sources.csv")

# harmonizing ages before 2020
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
min_open_age <- 
  d1 %>% 
  # filter(Year < 2020) %>% 
  select(Source, Country, Year, Sex, Age) %>% 
  unique() %>% 
  group_by(Source, Country, Year, Sex) %>% 
  filter(Age == max(Age)) %>% 
  ungroup() %>% 
  group_by(Source, Country, Sex) %>% 
  filter(Age == min(Age)) %>% 
  ungroup() %>% 
  rename(min_open = Age) %>% 
  select(-Year) %>% 
  unique()

# Largest age groups before 2020
ref_ages <- 
  d1 %>% 
  # filter(Year < 2020) %>% 
  # no ages surpassing the minimum open age interval
  left_join(min_open_age) %>% 
  filter(Age <= min_open) %>% 
  select(Source, Country, Year, Age) %>% 
  unique() %>% 
  # calculate age span in each interval
  group_by(Source, Country, Year) %>% 
  mutate(span = ifelse(Age == max(Age),
                       Inf,
                       lead(Age) - Age)) %>%
  ungroup() %>% 
  # identify and choose the largest age span in all periods for each age
  group_by(Source, Country, Age) %>% 
  mutate(max_span = max(span)) %>% 
  ungroup() %>% 
  select(Source, Country, Age, max_span) %>% 
  unique() %>% 
  group_by(Source, Country) %>%
  mutate(test = ifelse(Age == 0, 0, lag(Age + max_span))) %>% 
  ungroup() %>% 
  filter(Age == test) %>% 
  select(-test, -max_span)

# harmonizing ages within countries before 2020 
d2 <- 
  d1 %>%
  # filter(Year < 2020) %>% 
  group_by(Source, Country, Code, Year, Sex) %>% 
  do(harmon_age_intervals(chunk = .data)) %>% 
  ungroup()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Adding population exposures to each age interval
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

p1 <- 
  read_csv("data_inter/offsets.csv") %>% 
  mutate()

cts_yrs <- 
  d2 %>% 
  select(Country, Year, Sex) %>% 
  unique()

p2 <- 
  p1 %>% 
  inner_join(cts_yrs)


p3 <- 
  p2 %>% 
  group_by(Country, Year, Sex) %>% 
  do(assign_age_invals_pop(chunk = .data)) %>% 
  ungroup()

d3 <- 
  d2 %>% 
  left_join(p3)

no_pop <- 
  d3 %>% 
  filter(is.na(Population)) %>% 
  pull(Country) %>% 
  unique

write_csv(d3, "data_inter/annual_deaths_countries_selected_sources_harm.csv")

available2 <- 
  d3 %>% 
  group_by(Country, Code, Year, Source, Sex) %>% 
  summarise(Age_groups = n(),
            Deaths = sum(Deaths),
            Population = sum(Population)) %>% 
  ungroup() %>% 
  filter(Sex != "t") %>% 
  group_by(Country, Code, Year, Source, Age_groups) %>% 
  summarise(Sex_groups = n(),
            Deaths = sum(Deaths),
            Population = sum(Population)) %>% 
  ungroup() %>% 
  group_by(Country, Code, Source) %>% 
  summarise(Period = paste(min(Year), max(Year), sep = "-"),
            Deaths = sum(Deaths),
            Population = mean(Population),
            # putting together all the age and sex groups configurations between 2015 and 2021
            Age_groups = unique(Age_groups) %>% paste(collapse = ", "),
            Sex_groups = unique(Sex_groups) %>% paste(collapse = ", ")) %>% 
  ungroup() 

write_excel(available2)
