rm(list=ls())
source("R/00_functions.R")

# loading pre-processed data
hmd  <- read_csv("data_inter/hmd.csv")
who  <- read_csv("data_inter/who.csv")
unpd <- read_csv("data_inter/unpd.csv")
stmf <- read_csv("data_inter/stmf.csv") 
eurs_wk <- read_csv("data_inter/eurs_weekly.csv") 
eurs_an <- read_csv("data_inter/eurs_annual.csv") 

bra <- read_csv("data_inter/brazil.csv") %>% mutate(Source = "direct") %>% 
  group_by(Year, Sex) %>% 
  mutate(age_spn = ifelse(Age == max(Age), -1, lead(Age) - Age)) %>% 
  ungroup()

per <- read_csv("data_inter/peru.csv") %>% mutate(Source = "direct") %>% 
  group_by(Year, Sex) %>% 
  mutate(age_spn = ifelse(Age == max(Age), -1, lead(Age) - Age)) %>% 
  ungroup()

# mex <- read_csv("data_inter/mexico.csv") %>% mutate(Source = "direct")
zaf <- read_csv("data_inter/south_africa.csv") %>% mutate(Source = "direct") %>% 
  group_by(Year, Sex) %>% 
  mutate(age_spn = ifelse(Age == max(Age), -1, lead(Age) - Age)) %>% 
  ungroup()
# irn <- read_csv("data_inter/iran.csv") %>% mutate(Source = "direct")
# mse <- read_csv("data_inter/msemburi_tag.csv") 

# Exclude regions already included in countries, 
# as well as problematic data (e.g., Kenya)
cts_exclude <- tibble(Country = 
                        c("Scotland", 
                          "England and Wales", 
                          "Northern Ireland",
                          "Kenya",
                          "Puerto Rico",
                          "Anguilla"))

# unique(out$Country)

unique(unpd$Country)
unique(unpd$Source)
# Several direct sources have been superseded by UNPD and others

# putting all together
all_in <- 
  bind_rows(hmd,
            who,
            unpd,
            stmf,
            eurs_wk,
            eurs_an,
            bra,
            per,
            # mex,
            zaf,
            # irn
            ) %>% 
  replace_na(list(Deaths = 0)) %>% 
  group_by(Source, Country) %>% 
  filter(max(Year) >= 2020 & min(Year) <= 2017) %>% 
  # ungroup() %>% 
  # mutate(Age = ifelse(Age >= 100, 100, Age)) %>% 
  # group_by(Country, Code, Source, Year, Sex, Age, age_spn) %>% 
  # summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  unique()


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
  mutate(age_spn = ifelse(Age == max(Age), -1, lead(Age)-Age)) %>% 
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
  arrange(Country) %>% 
  mutate(last = str_sub(period, 6, 9) %>% as.integer(),
         yrs_bsn = 2020 - (str_sub(period, 1, 4) %>% as.integer()))

write_excel(sum_all)

unique(sum_all$Source)

# selection based on source
# ~~~~~~~~~~~~~~~~~~~~~~~~~

# source priority
s1 <- "direct"
s2 <- "unpd_crvs"
s3 <- "hmd"
s4 <- "eurs_annual"
s5 <- "stmf"
s6 <- "eurs_weekly"
s7 <- "unpd_hmd"
s8 <- "unpd_dy"
s9 <- "who_mort_db"

sel_source <- 
  sum_all %>% 
  group_by(Country) %>% 
  filter((any(Source == s1) & Source == s1) | all(Source != s1)) %>%
  filter((any(Source == s2) & Source == s2) | all(Source != s2)) %>%
  filter((any(Source == s3) & Source == s3) | all(Source != s3)) %>%
  filter((any(Source == s4) & Source == s4) | all(Source != s4)) %>%
  filter((any(Source == s5) & Source == s5) | all(Source != s5)) %>%
  filter((any(Source == s6) & Source == s6) | all(Source != s6)) %>%
  filter((any(Source == s7) & Source == s7) | all(Source != s7)) %>%
  filter((any(Source == s8) & Source == s8) | all(Source != s8)) %>%
  filter((any(Source == s9) & Source == s9) | all(Source != s9)) %>%
  unique() %>% 
  mutate(n = n()) %>% 
  ungroup()

# selection based on data
sel_criteria <- 
  sum_all %>% 
  group_by(Country) %>% 
  filter(last == max(last)) %>% 
  filter(ages >= max(ages)) %>% 
  filter(years == max(years)) %>% 
  filter(sexs == max(sexs)) %>% 
  filter(Deaths == max(Deaths)) %>% 
  unique() %>% 
  filter(1:n() == n()) %>% 
  mutate(n = n()) %>% 
  ungroup()

# comparison between both selection methods
comp <- 
  sel_source %>% 
  select(Country, ages_s = ages, years_s = years, Source_s = Source, last_s = last) %>%
  left_join(sel_criteria %>% 
              select(Country, ages_c = ages, years_c = years, Source_c = Source, last_c = last)) %>% 
  mutate(same = ifelse(Source_s == Source_c, 1, 0),
         ages_ratio = ages_c / ages_s,
         years_ratio = years_c / years_s)

# selection based on infant mortality
sel_crit_inf <- 
  sum_all %>% 
  group_by(Country) %>% 
  filter(infd == max(infd)) %>% 
  filter(last == max(last)) %>% 
  filter(ages >= max(ages)) %>%
  filter(yrs_bsn == max(yrs_bsn)) %>% 
  # filter(sexs == max(sexs)) %>% 
  filter(Deaths == max(Deaths)) %>%
  unique() %>% 
  # filter(1:n() == n()) %>% 
  mutate(n = n()) %>% 
  ungroup()

# selection based on 2021 data
sel_crit_pan <- 
  sum_all %>% 
  group_by(Country) %>% 
  filter(last == max(last)) %>% 
  filter(infd == max(infd)) %>% 
  filter(ages >= max(ages)) %>%
  filter(yrs_bsn == max(yrs_bsn)) %>% 
  # filter(sexs == max(sexs)) %>% 
  filter(Deaths == max(Deaths)) %>%
  unique() %>% 
  # filter(1:n() == n()) %>% 
  mutate(n = n()) %>% 
  ungroup()


# comparison between both selection methods
comp2 <- 
  sel_crit_inf %>% 
  select(Country, ages_s = ages, years_s = years, Source_s = Source, last_s = last, infd_s = infd) %>%
  left_join(sel_crit_pan %>% 
              select(Country, ages_c = ages, years_c = years, Source_c = Source, last_c = last, infd_c = infd)) %>% 
  mutate(same = ifelse(Source_s == Source_c, 1, 0),
         ages_ratio = ages_c / ages_s,
         years_ratio = years_c / years_s)


# sel_mixed <- 
#   sum_all %>% 
#   group_by(Country) %>% 
#   filter((any(Source == "direct") & Source == "direct") | all(Source != "direct")) %>%
#   filter((any(Source == "hmd") & Source == "hmd") | all(Source != "hmd")) %>%
#   filter(ages == max(ages)) %>% 
#   filter(years == max(years)) %>% 
#   filter(sexs == max(sexs)) %>% 
#   filter(Deaths == max(Deaths)) %>% 
#   unique() 




# selecting best source for each population
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# prefer selection based on source, unless content-based selection has either
# 20% more ages groups or periods
best_source <- 
  comp %>% 
  mutate(Source = ifelse(ages_ratio > 1.2 | years_ratio > 1.2, Source_c, Source_s)) %>% 
  select(Country, Source)

all_in4 <- 
  all_in3 %>% 
  inner_join(best_source) %>% 
  arrange(Country, Year, Sex, Age) %>% 
  mutate(Country = ifelse(Country == "Faeroe Islands", "Faroe Islands", Country)) %>% 
  anti_join(cts_exclude)

# unique(out$Country)
test <- 
  all_in4 %>% 
  select(Country, Age, age_spn) %>% 
  unique() %>% 
  group_by(Country, Age) %>% 
  summarise(n = n())


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Adding population exposures to each age interval
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ref_ages <- 
  all_in4 %>% 
  select(Country, Year, Age, age_spn) %>% 
  unique()

p1 <- 
  read_csv("data_inter/offsets.csv") %>% 
  mutate()

cts_yrs <- 
  all_in4 %>% 
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

all_in5 <- 
  all_in4 %>% 
  left_join(p3)

no_pop <- 
  all_in5 %>% 
  filter(is.na(Population)) %>% 
  pull(Country) %>% 
  unique

# 3) output mortality data based on selected sources
write_csv(all_in5, "data_inter/annual_deaths_countries_selected_sources.csv")

# summary of selected sources by country 
available <- 
  all_in5 %>% 
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

# # Harmonizing age groups in all periods
# # =====================================
# d1 <- 
#   read_csv("data_inter/annual_deaths_countries_selected_sources.csv")
# 
# span <- 
#   d1 %>% 
#   filter(Year < 2020) %>% 
#   select(Country, Age) %>%
#   unique() %>% 
#   group_by(Country) %>% 
#   mutate(age_int = lead(Age) - Age) %>% 
#   filter(Age == 0)
# 
# # harmonizing ages before 2020
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# min_open_age <- 
#   d1 %>% 
#   # filter(Year < 2020) %>% 
#   select(Source, Country, Year, Sex, Age) %>% 
#   unique() %>% 
#   group_by(Source, Country, Year, Sex) %>% 
#   filter(Age == max(Age)) %>% 
#   ungroup() %>% 
#   group_by(Source, Country, Sex) %>% 
#   filter(Age == min(Age)) %>% 
#   ungroup() %>% 
#   rename(min_open = Age) %>% 
#   select(-Year) %>% 
#   unique()
# 
# # Largest age groups before 2020
# ref_ages <- 
#   d1 %>% 
#   # filter(Year < 2020) %>% 
#   # no ages surpassing the minimum open age interval
#   left_join(min_open_age) %>% 
#   filter(Age <= min_open) %>% 
#   select(Source, Country, Year, Age) %>% 
#   unique() %>% 
#   # calculate age span in each interval
#   group_by(Source, Country, Year) %>% 
#   mutate(span = ifelse(Age == max(Age),
#                        Inf,
#                        lead(Age) - Age)) %>%
#   ungroup() %>% 
#   # identify and choose the largest age span in all periods for each age
#   group_by(Source, Country, Age) %>% 
#   mutate(max_span = max(span)) %>% 
#   ungroup() %>% 
#   select(Source, Country, Age, max_span) %>% 
#   unique() %>% 
#   group_by(Source, Country) %>%
#   mutate(test = ifelse(Age == 0, 0, lag(Age + max_span))) %>% 
#   ungroup() %>% 
#   filter(Age == test) %>% 
#   select(-test, -max_span)
# 
# # harmonizing ages within countries before 2020 
# d2 <- 
#   d1 %>%
#   # filter(Year < 2020) %>% 
#   group_by(Source, Country, Code, Year, Sex) %>% 
#   do(harmon_age_intervals(chunk = .data)) %>% 
#   ungroup()
# 
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Adding population exposures to each age interval
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# p1 <- 
#   read_csv("data_inter/offsets.csv") %>% 
#   mutate()
# 
# cts_yrs <- 
#   d2 %>% 
#   select(Country, Year, Sex) %>% 
#   unique()
# 
# p2 <- 
#   p1 %>% 
#   inner_join(cts_yrs)
# 
# 
# p3 <- 
#   p2 %>% 
#   group_by(Country, Year, Sex) %>% 
#   do(assign_age_invals_pop(chunk = .data)) %>% 
#   ungroup()
# 
# d3 <- 
#   d2 %>% 
#   left_join(p3)
# 
# no_pop <- 
#   d3 %>% 
#   filter(is.na(Population)) %>% 
#   pull(Country) %>% 
#   unique
# 
# write_csv(d3, "data_inter/annual_deaths_countries_selected_sources_harm.csv")
# 
# available2 <- 
#   d3 %>% 
#   group_by(Country, Code, Year, Source, Sex) %>% 
#   summarise(Age_groups = n(),
#             Deaths = sum(Deaths),
#             Population = sum(Population)) %>% 
#   ungroup() %>% 
#   filter(Sex != "t") %>% 
#   group_by(Country, Code, Year, Source, Age_groups) %>% 
#   summarise(Sex_groups = n(),
#             Deaths = sum(Deaths),
#             Population = sum(Population)) %>% 
#   ungroup() %>% 
#   group_by(Country, Code, Source) %>% 
#   summarise(Period = paste(min(Year), max(Year), sep = "-"),
#             Deaths = sum(Deaths),
#             Population = mean(Population),
#             # putting together all the age and sex groups configurations between 2015 and 2021
#             Age_groups = unique(Age_groups) %>% paste(collapse = ", "),
#             Sex_groups = unique(Sex_groups) %>% paste(collapse = ", ")) %>% 
#   ungroup() 
# 
# write_excel(available2)
