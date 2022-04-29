library(here)
source(here("R", "00_functions.R"))

# loading pre-processed data
hmd <- read_csv("Output/hmd.csv")
who <- read_csv("Output/who.csv")
unpd <- read_csv("Output/unpd.csv")
stmf <- read_csv("Output/stmf.csv") 
eurs <- read_csv("Output/eurs.csv") 
zaf <- read_csv("Output/south_africa.csv") %>% mutate(Source = "direct")
irn <- read_csv("Output/iran.csv") %>% mutate(Source = "direct")

# These sources have been superceded by other sources (mostly by UNPD data)
# bra <- read_csv("Output/brazil.csv") 
# mex <- read_csv("Output/mexico.csv")
# per <- read_csv("Output/peru.csv")
# col <- read_csv("Output/colombia.csv")
# chl <- read_csv("Output/chile.csv")
# ecu <- read_csv("Output/ecuador.csv",
#                  col_types = "ccdccdc")
 
all_in <- 
  bind_rows(hmd,
            who,
            unpd,
            stmf,
            eurs,
            zaf,
            irn)

# imputing unknown ages and sexes
all_in_adj <- 
  all_in %>% 
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
# harmonization of ages before 2020
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# harmonizing to the minimum common amount of age intervals before 2020

age_groups <- 
  all_in_adj %>% 
  select(Source, Country, Year, Age) %>% 
  unique() %>% 
  group_by(Source, Country, Year) %>% 
  summarise(ages = n()) %>% 
  ungroup()

# identifying countries that need harmonization
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

# selecting the reference period 
ref_period <- 
  age_groups %>% 
  left_join(cts_to_harmonize) %>% 
  filter(Year < 2020 & to_harmon == 1) %>% 
  # selecting the periods with the minimum amount of age categories
  group_by(Source, Country) %>% 
  mutate(min_age = min(ages))
  
  filter(ages == min(ages)) %>% 
  filter(Year == min(Year)) %>% 
  ungroup() %>% 
  mutate(is_reference = 1) %>% 
  select(Source, Country, Year, is_reference) 

# pulling reference age groups
ref_ages <- 
  all_in_adj %>%
  left_join(ref_period) %>% 
  filter(Sex == "t",
         is_reference == 1) %>% 
  select(Source, Country, Age)

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

# unique(out$Sex)

# write_csv(out_all, "Output/tag_deaths_for_young.csv")
# write_csv(out_all, "C:/Users/kikep/OneDrive/Documents/gits/unicef_excess/Data/annual_deaths_countries_selected_sources_young.csv")

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




















# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # selecting best source for each country
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # based on data previous to 2020
# # criteria:
# # 1. Guarantee continuity of the same source
# # 2. Finer age categories
# # 3. More death count registration
# # 4. Longer period of observation
# 
# source_selec <- 
#   db3 %>% 
#   filter(Year <= 2019, Sex == "t") %>% 
#   group_by(Country, Source, Year) %>% 
#   # identifying amount of age groups
#   summarise(Deaths = sum(Deaths),
#             age_groups = n()) %>% 
#   ungroup() %>% 
#   group_by(Country) %>% 
#   # evaluating how different are sources in age categories 
#   mutate(prop = age_groups / max(age_groups)) %>% 
#   ungroup() %>% 
#   # exclude sources in which age groups reduced more than 30% 
#   filter(prop > 0.7) %>% 
#   group_by(Country, Year) %>% 
#   # comparing death registration across sources
#   mutate(prop_d = Deaths / max(Deaths)) %>% 
#   arrange(Country, Year, -Deaths) %>% 
#   # exclude those sources with a >= 2% reduction in deaths 
#   filter(prop_d > 0.98) %>% 
#   ungroup() %>% 
#   group_by(Country, Source) %>% 
#   # compare the period coverage across sources
#   mutate(years = n()) %>% 
#   ungroup() %>% 
#   group_by(Country) %>% 
#   # select the sources with longest period of observations
#   filter(years == max(years)) %>% 
#   ungroup() %>% 
#   group_by(Country, Source) %>% 
#   # identify the minimum common amount of age groups
#   mutate(min_ages = min(age_groups)) %>% 
#   group_by(Country) %>% 
#   filter(min_ages == max(min_ages)) %>% 
#   ungroup() %>% 
#   group_by(Country, Source) %>% 
#   # looking at total deaths during the observed periods
#   mutate(all_deaths = sum(Deaths)) %>% 
#   group_by(Country) %>% 
#   # selecting the source with highest count of deaths 
#   filter(all_deaths == max(all_deaths)) %>% 
#   ungroup() %>% 
#   # looking at countries with more than one source
#   group_by(Country, Year) %>% 
#   arrange(Country, Year, Source) %>% 
#   mutate(order = 1:n()) %>% 
#   ungroup() %>% 
#   # selecting the first listed source
#   filter(order == 1) %>% 
#   ungroup() %>% 
#   select(Country, Source, Year, age_groups, min_ages)
# 
# # list of best source by country
# best_source_by_country <- 
#   source_selec %>% 
#   select(Country, Source) %>% 
#   unique()
# 
# # filtering best sources by country
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# db4 <- 
#   db3 %>% 
#   left_join(source_selec %>% 
#               select(-age_groups, -min_ages) %>% 
#               mutate(Source = ifelse(Country == "United Kingdom", "eurs", Source),
#                      best_source = 1)) %>% 
#   filter(best_source == 1) %>% 
#   select(-best_source)


# harmonization of ages before 2020
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# harmonizing to the minimum common amount of age intervals before 2020

# selecting countries that need harmonization
cts_to_harmonize <- 
  source_selec %>% 
  select(Country, age_groups) %>% 
  unique() %>% 
  group_by(Country) %>% 
  summarise(age_configs = n()) %>% 
  filter(age_configs > 1,
         Country != "United Kingdom") %>% 
  pull(Country)

# selecting the reference period 
age_groups_period <- 
  source_selec %>% 
  filter(Country %in% cts_to_harmonize) %>% 
  # selecting the periods with the minimum amount of age categories
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

# ~~~~~~~~~~~~~~
# final database
# ~~~~~~~~~~~~~~

# data in 2020
db2020 <- 
  db3 %>% 
  filter(Year >= 2020) %>% 
  left_join(best_source_by_country %>% 
              mutate(Source = ifelse(Country == "United Kingdom", "stmf_adj", Source),
                     best_source = 1)) %>% 
  filter(best_source == 1) %>% 
  select(-best_source)

cts2020 <- db2020 %>% pull(Country) %>% unique() %>% sort()

# adding harmonized data with age groups previous to 2020 
# and all age groups in 2020
out <- 
  db4 %>% 
  # adding the harmonized ages
  filter(!Country %in% cts_to_harmonize) %>% 
  bind_rows(harmonized_ages2) %>%
  # selecting countries with data in 2020
  filter(Country %in% cts2020) %>% 
  # adding data for 2020
  bind_rows(db2020) %>%
  arrange(Country, Year, Sex, Age)

unique(out$Country) %>% sort()

# comparing countries in original and output databases
cts_diff <- 
  db3 %>% 
  select(Country) %>% 
  unique() %>% 
  left_join(out %>% 
                select(Country) %>% 
                unique() %>%
                mutate(is_end = 1))

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# saving 3 output data objects
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1) summary of all available sources by country 
write_csv(summ1, "Output/summary_all_sources_by_country.csv")

# 2) summary of selected sources by country
write_csv(available, "Output/summary_selected_sources_by_country.csv")

# 3) output mortality data based on selected sources
write_csv(out, "Output/annual_deaths_countries_selected_sources.csv")

unique(summ1$Country)
