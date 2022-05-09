source("R/00_functions.R")

# loading pre-processed data
hmd  <- read_csv("Output/hmd.csv")
who  <- read_csv("Output/who.csv")
unpd <- read_csv("Output/unpd.csv")
stmf <- read_csv("Output/stmf.csv") 
eurs <- read_csv("Output/eurs.csv") 
zaf <- read_csv("Output/south_africa.csv") %>% mutate(Source = "direct")
irn <- read_csv("Output/iran.csv") %>% mutate(Source = "direct")
mse <- read_csv("Output/msemburi_tag.csv") 


# Several direct sources have been superseded by UNPD and others

# putting all together
all_in <- 
  bind_rows(hmd,
            who,
            unpd,
            stmf,
            eurs,
            zaf,
            irn)

# excluding sources with insufficient data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sources to exclude because of insufficient periods for baseline (minimum 3 periods)
insuf_pers <- 
  all_in %>% 
  select(Country, Source, Year) %>% 
  unique() %>% 
  filter(Year < 2020) %>% 
  group_by(Country, Source) %>% 
  summarise(yrs = n()) %>% 
  ungroup() %>% 
  filter(yrs < 3) %>% 
  mutate(exc_yrs = 1) %>% 
  select(-yrs)

# sources to exclude because of insufficient age groups (minimum 4 intervals)
insuf_ages <- 
  all_in %>% 
  group_by(Country, Source, Year, Sex) %>% 
  summarise(ages = n()) %>% 
  ungroup() %>% 
  filter(ages < 4) %>% 
  mutate(exc_ags = 1) %>% 
  select(-ages)


# imputing unknown ages and sexes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_in_adj <- 
  all_in %>% 
  left_join(insuf_pers) %>% 
  left_join(insuf_ages) %>% 
  filter(is.na(exc_yrs) & is.na(exc_ags)) %>% 
  select(-exc_yrs, -exc_ags) %>% 
  group_by(Country, Sex, Year, Source) %>% 
  do(rescale_age(chunk = .data)) %>% 
  ungroup() %>% 
  group_by(Country, Age, Year, Source) %>% 
  do(rescale_sex(chunk = .data)) %>% 
  ungroup() %>% 
  mutate(Age = as.double(Age)) %>% 
  arrange(Source, Country, Year, Sex, Age) %>% 
  replace_na(list(Deaths = 0))


# harmonizing ages before 2020
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# harmonizing to the minimum common amount of age intervals before 2020

age_groups <- 
  all_in_adj %>% 
  select(Source, Country, Year, Age) %>% 
  unique() %>% 
  group_by(Source, Country, Year) %>% 
  summarise(ages = n()) %>% 
  ungroup()

# identifying countries that need age harmonization
cts_to_harmonize <- 
  age_groups %>% 
  filter(Year < 2020) %>% 
  select(Source, Country, ages) %>% 
  unique() %>% 
  group_by(Source, Country) %>% 
  summarise(age_configs = n()) %>% 
  filter(age_configs > 1) %>% 
  select(Source, Country) %>% 
  mutate(to_harmon = 1) %>% 
  ungroup()

# selecting the reference period for age harmonization
ref_period <- 
  age_groups %>% 
  left_join(cts_to_harmonize) %>% 
  filter(Year < 2020 & to_harmon == 1) %>% 
  # selecting the periods with the minimum amount of age categories
  group_by(Source, Country) %>% 
  # mutate(min_age = min(ages))
  filter(ages == min(ages)) %>% 
  filter(Year == min(Year)) %>% 
  ungroup() %>% 
  mutate(is_reference = 1) %>% 
  select(Source, Country, Year, is_reference) 

# closing age groups
age_close <- 
  all_in_adj %>% 
  left_join(cts_to_harmonize) %>% 
  filter(Year < 2020 & to_harmon == 1) %>% 
  filter(Sex == "t") %>% 
  group_by(Country, Year, Source) %>% 
  filter(Age == max(Age)) %>% 
  group_by(Country, Source) %>% 
  filter(Age == min(Age)) %>% 
  select(Country, Source, Age_close = Age) %>% 
  unique()
  
# pulling reference age groups
ref_ages <- 
  all_in_adj %>%
  left_join(ref_period) %>% 
  filter(Sex == "t",
         is_reference == 1) %>% 
  select(Source, Country, Age) %>% 
  left_join(age_close) %>% 
  filter(Age <= Age_close) %>% 
  select(-Age_close)

harmonized_ages <- 
  all_in_adj %>%
  left_join(cts_to_harmonize) %>% 
  filter(Year < 2020 & to_harmon == 1) %>% 
  group_by(Source, Country, Code, Year, Sex) %>% 
  do(assign_age_intervals(chunk = .data)) %>% 
  ungroup()

all_in_adj2 <- 
  all_in_adj %>% 
  left_join(cts_to_harmonize) %>% 
  filter(is.na(to_harmon) | Year >= 2020) %>% 
  select(-to_harmon) %>% 
  bind_rows(harmonized_ages) %>% 
  arrange(Country, Source, Year, Sex, Age)

# Summarizing content
# ~~~~~~~~~~~~~~~~~~~

# summarizing content by source
sum_all <- 
  sum_source(all_in_adj2) %>% 
  unique() %>% 
  arrange(Country)

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
  mutate(n = n())

# selection based on data
sel_criteria <- 
  sum_all %>% 
  group_by(Country) %>% 
  filter(ages == max(ages)) %>% 
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
  select(Country, ages_s = ages, years_s = years, Source_s = Source) %>%
  left_join(sel_criteria %>% 
              select(Country, ages_c = ages, years_c = years, Source_c = Source)) %>% 
  mutate(ages_ratio = ages_c / ages_s,
         years_ratio = years_c / years_s)
  

# selecting best source for each population
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# prefer selection based on source, unless content-based selection has either
# 20% more ages groups or periods
best_source <- 
  comp %>% 
  mutate(Source = ifelse(ages_ratio > 1.2 | years_ratio > 1.2, Source_c, Source_s)) %>% 
  select(Country, Source) %>% 
  mutate(best = "y")

out <- 
  all_in_adj2 %>% 
  left_join(best_source) %>% 
  filter(best == "y") %>% 
  arrange(Country, Year, Sex, Age) %>% 
  select(-best)


# 3) output mortality data based on selected sources
write_csv(out, "Output/annual_deaths_countries_selected_sources.csv")

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
  group_by(Country, Code, Source, Age_groups, Sex_groups) %>% 
  summarise(Period = paste(min(Year), max(Year), sep = "-"),
            Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  group_by(Country) %>% 
  mutate(age_groupings = n()) %>% 
  ungroup() 

out <- read_csv("Output/annual_deaths_countries_selected_sources.csv")
unique(out$Country)
# 
# unique(mse$Country) %>% sort()
# length(unique(mse$Code))
# 
# 
# test <- 
#   out %>% 
#   select(Code) %>% 
#   unique() %>% 
#   mutate(our = 1) %>% 
#   full_join(mse %>% 
#               select(Code) %>% 
#               unique() %>% 
#               mutate(mse = 1))
# 
# test2 <- 
  


