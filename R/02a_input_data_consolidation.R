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
  # mutate(Country = ifelse(Country == "China, Hong Kong SAR", "Hong Kong", Country)) %>% 
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
  select(Source, Country, Year, Sex, Age) %>% 
  unique() %>% 
  # calculate age span in each interval
  group_by(Source, Country, Year, Sex) %>% 
  mutate(span = ifelse(Age == max(Age),
                       Inf,
                       lead(Age) - Age)) %>%
  ungroup() %>% 
  # identify and choose the largest age span in all periods for each age
  group_by(Source, Country, Sex, Age) %>% 
  mutate(max_span = max(span)) %>% 
  ungroup() %>% 
  select(Source, Country, Sex, Age, max_span) %>% 
  unique() %>% 
  group_by(Source, Country, Sex) %>%
  mutate(test = ifelse(Age == 0, 0, lag(Age + max_span))) %>% 
  ungroup() %>% 
  filter(Age == test) %>% 
  select(-test, -max_span)

# harmonizing ages within countries before 2020 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # # a couple of tests
# chunk <-
#   all_in2 %>%
#   filter(Country == "Armenia",
#          Year == 2015,
#          Source == "unpd_dy")
# 
# chunk <- 
#   all_in2 %>%
#   filter(Year < 2020,
#          Country == "Germany",
#          Year == 2015,
#          Sex == "t",
#          Source == "eurs_weekly")

harmonized_ages <- 
  all_in2 %>%
  filter(Year < 2020) %>% 
  group_by(Source, Country, Code, Year, Sex) %>% 
  do(harmon_age_intervals(chunk = .data)) %>% 
  arrange(Age) %>% 
  mutate(age_spn = ifelse(Age == max(Age), -1, lead(Age)-Age)) %>% 
  ungroup() 

# together harmonized ages before 2020, and original age groups since 2020
all_in3 <- 
  bind_rows(harmonized_ages, 
            all_in2 %>% 
              filter(Year >= 2020)) %>% 
  arrange(Country, Source, Year, Sex, Age) %>% 
  mutate(Country = ifelse(Country == "Faeroe Islands", "Faroe Islands", Country)) %>% 
  anti_join(cts_exclude) %>% 
  select(Country, Code, Year, Sex, Age, age_spn, Source, Deaths)
  
  
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


# selection based on infant mortality
sel_crit_inf <- 
  sum_all %>% 
  group_by(Country) %>% 
  filter(infd == max(infd)) %>% 
  filter(ages >= max(ages)) %>%
  filter(last == max(last)) %>% 
  filter(yrs_bsn == max(yrs_bsn)) %>% 
  filter(Deaths == max(Deaths)) %>%
  unique() %>% 
  filter(1:n() == n()) %>%
  mutate(n = n()) %>% 
  ungroup()

# selection based on the last observed year
sel_crit_pan <- 
  sum_all %>% 
  group_by(Country) %>% 
  filter(last == max(last)) %>% 
  filter(infd == max(infd)) %>% 
  filter(ages >= max(ages)) %>%
  filter(yrs_bsn == max(yrs_bsn)) %>% 
  filter(Deaths == max(Deaths)) %>%
  unique() %>% 
  filter(1:n() == n()) %>%
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



# selecting best source for each population
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# prefer selection based on source, unless content-based selection has either
# 20% more ages groups or periods
best_inf <- 
  sel_crit_inf %>% 
  select(Country, Source)

best_pan <- 
  sel_crit_pan %>% 
  select(Country, Source)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# selecting best source for each country
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

inf <- 
  all_in3 %>% 
  inner_join(best_inf) %>% 
  arrange(Country, Year, Sex, Age)

pan <- 
  all_in3 %>% 
  inner_join(best_pan) %>% 
  arrange(Country, Year, Sex, Age)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Adding population exposures to each age interval
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p1 <- 
  read_csv("data_inter/offsets.csv") %>% 
  mutate()

# for data with infant priority
ref_ages <- 
  inf %>% 
  select(Country, Year, Sex, Age, age_spn) %>% 
  unique()

cts_yrs_i <- 
  inf %>% 
  select(Country, Year, Sex) %>% 
  unique()

pi <- 
  p1 %>% 
  inner_join(cts_yrs_i)

test <- 
  pi %>% 
  filter(Country == "Armenia",
         Sex == "t",
         Year == 2015)

pi2 <- 
  pi %>% 
  group_by(Country, Year, Sex) %>% 
  do(assign_age_invals_pop(chunk = .data)) %>% 
  ungroup()

inf2 <- 
  inf %>% 
  left_join(pi2)

# for data with pandemic period priority
ref_ages <- 
  pan %>% 
  select(Country, Year, Sex, Age, age_spn) %>% 
  unique()

cts_yrs_p <- 
  pan %>% 
  select(Country, Year, Sex) %>% 
  unique()

pp <- 
  p1 %>% 
  inner_join(cts_yrs_p)

pp2 <- 
  pp %>% 
  group_by(Country, Year, Sex) %>% 
  do(assign_age_invals_pop(chunk = .data)) %>% 
  ungroup()

pan2 <- 
  pan %>% 
  left_join(pp2)

no_pop <- 
  bind_rows(pan2, inf2) %>% 
  filter(is.na(Population)) %>% 
  pull(Country) %>% 
  unique

# 3) output mortality data based on selected sources
write_csv(inf2, "data_inter/deaths_sourced_infant_based.csv")
write_csv(pan2, "data_inter/deaths_sourced_period_based.csv")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# summary of selected sources by country 
available_inf <- 
  inf2 %>% 
  group_by(Country, Code, Year, Source, Sex) %>% 
  summarise(Age_groups = n(),
            infd = ifelse(any(Age == 0 & age_spn == 1), 1, 0),
            Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  group_by(Country, Code, Year, Source) %>% 
  summarise(Age_groups = min(Age_groups),
            infd = min(infd),
            Sex_groups = n(),
            Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  group_by(Country, Code, Source) %>% 
  summarise(Period = paste(min(Year), max(Year), sep = "-"),
            Deaths = sum(Deaths),
            infd = min(infd),
            # putting together all the age and sex groups configurations between 2015 and 2021
            Age_groups = unique(Age_groups) %>% paste(collapse = ", "),
            Sex_groups = unique(Sex_groups) %>% paste(collapse = ", ")) %>% 
  ungroup() %>% 
  mutate(last = str_sub(Period, 6, 9) %>% as.integer(),
         yrs_bsn = 2020 - (str_sub(Period, 1, 4) %>% as.integer()))

write_excel(available_inf)

# summary of selected sources by country 
available_pan <- 
  pan2 %>% 
  group_by(Country, Code, Year, Source, Sex) %>% 
  summarise(Age_groups = n(),
            infd = ifelse(any(Age == 0 & age_spn == 1), 1, 0),
            Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  group_by(Country, Code, Year, Source) %>% 
  summarise(Age_groups = min(Age_groups),
            infd = min(infd),
            Sex_groups = n(),
            Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  group_by(Country, Code, Source) %>% 
  summarise(Period = paste(min(Year), max(Year), sep = "-"),
            Deaths = sum(Deaths),
            infd = min(infd),
            # putting together all the age and sex groups configurations between 2015 and 2021
            Age_groups = unique(Age_groups) %>% paste(collapse = ", "),
            Sex_groups = unique(Sex_groups) %>% paste(collapse = ", ")) %>% 
  ungroup() %>% 
  mutate(last = str_sub(Period, 6, 9) %>% as.integer(),
         yrs_bsn = 2020 - (str_sub(Period, 1, 4) %>% as.integer()))

write_excel(available_pan)

