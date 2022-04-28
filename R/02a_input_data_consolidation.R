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
irn <- read_csv("Output/iran.csv")
hmd <- read_csv("Output/hmd_2020.csv")

# adjusting sources for "country_public" data
latam <- 
  bind_rows(bra, 
            mex,
            per,
            ecu,
            col,
            chl,
            ecu) 

tot_latam <- 
  latam %>% 
  filter(Age == "TOT")

latam_gr5 <- 
  latam %>% 
  filter(Age != "TOT") %>% 
  mutate(Source = "country_public",
         Age = Age %>% as.double(),
         Age = case_when(Age == 0 ~ 0,
                         Age %in% 1:4 ~ 1,
                         Age >= 5 ~ Age - Age %% 5)) %>% 
  group_by(Country, Code, Year, Sex, Age, Source) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Age = Age %>% as.character()) %>% 
  bind_rows(tot_latam)
  
unique(latam_gr5$Age)


# excluding HMD countries from stmf and eurs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stmf2 <- 
  stmf %>% 
  filter(!Country %in% unique(hmd$Country))

eurs2 <- 
  eurs %>% 
  filter(!Country %in% unique(hmd$Country),
         !Country %in% unique(stmf2$Country))

unpd2 <- 
  unpd %>% 
  select(-AgeSpan) %>% 
  group_by(Country, Code, Year, Sex) %>% 
  mutate(max_age = max(Age)) %>% 
  group_by(Country, Code) %>% 
  mutate(min_max_age = min(max_age)) %>% 
  ungroup() %>% 
  mutate(Age = ifelse(Age >= min_max_age, min_max_age, Age)) %>% 
  group_by(Country, Code, Year, Sex, Age) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  # filter(!Country %in% unique(hmd$Country),
  #        !Country %in% unique(stmf2$Country),
  #        !Country %in% unique(eurs2$Country)) %>% 
  mutate(Source = "unpd") 

latam2 <- 
  latam_gr5 %>% 
  filter(!Country %in% unique(hmd$Country),
         !Country %in% unique(stmf2$Country),
         !Country %in% unique(eurs2$Country),
         !Country %in% unique(unpd2$Country))

who2 <- 
  who %>% 
  filter(!Country %in% unique(hmd$Country),
         !Country %in% unique(stmf2$Country),
         !Country %in% unique(eurs2$Country),
         !Country %in% unique(latam2$Country),
         !Country %in% unique(unpd2$Country)) %>% 
  select(-age_up)

# Adjust for unknown ages and sex
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# WHO data already solved by Tim
# UNPD does not have missing data
# solving for stmf, eurostat, and brazil

db <- 
  bind_rows(stmf,
            eurs,
            # latam2,
            hmd,
            zaf,
            irn)

unique(db$Sex)
unique(db$Age)

db_tots <- 
  db %>% 
  filter(Age == "TOT",
         Sex == "t")

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
  arrange(Source, Country, Year, Sex, Age) %>% 
  replace_na(list(Deaths = 0))

# Adding data in 2020 for the UK
# grouping deaths for the UK in 2020 from the STMF
# 5-year age groups and 90+

uk2020 <- 
  db2 %>% 
  filter(Code %in% c("GBR_SCO", "GBR_NIR", "GBRTENW"),
         Year == 2020) %>% 
  mutate(Age = case_when(Age > 90 ~ 90,
                         TRUE ~ Age)) %>% 
  group_by(Year, Sex, Age) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Country = "United Kingdom", 
         Code = "GBR",
         Source = "stmf_adj")

# merge all data together
db3 <- 
  bind_rows(db2, 
            who2,
            unpd2,
            # hmd,
            # ,
            # uk2020
            ) %>% 
  mutate(Deaths = ifelse(is.na(Deaths), 0, Deaths)) %>% 
  filter(!Code %in% c("GBR_SCO", "GBR_NIR", "GBRTENW")) %>% 
  drop_na(Sex)

# create a table with data availability by population, including age groups and years

unique(db3$Age)
unique(db3$Sex)
unique(db3$Country) %>% sort()


comp_all <- 
  db3 %>% 
  group_by(Country, Year, Sex, Source) %>% 
  mutate(ages = n()) %>% 
  ungroup() %>% 
  group_by(Country, Year, Age, Source) %>% 
  mutate(sexs = n()) %>% 
  ungroup() %>% 
  group_by(Country, Sex, Age, Source) %>% 
  mutate(years = n()) %>% 
  ungroup() %>% 
  group_by(Country, Source) %>% 
  filter(!(sexs == 3 & Sex == "t")) %>% 
  summarise(Deaths = sum(Deaths),
            ages = min(ages),
            sexs = min(sexs),
            years = min(years)) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(Country) %>% 
  # unique sources
  group_by(Country) %>% 
  # max age definition
  filter(ages == max(ages)) %>% 
  filter(years == max(years)) %>% 
  filter(sexs == max(sexs)) %>% 
  filter(Deaths == max(Deaths)) %>% 
  mutate(best = ifelse(n() == 1, Source, NA))

sel_all <- 
  comp_all %>% 
  select(Country, Source) %>% 
  mutate(keep = "y")

out_all <- 
  db3 %>% 
  left_join(sel_all) %>% 
  filter(keep == "y") %>% 
  arrange(Country, Sex, Age, Year) %>% 
  select(-keep)

unique(out_all$Country)

write_csv(out_all, "Output/tag_deaths_for_young.csv")
write_csv(out_all, "C:/Users/kikep/OneDrive/Documents/gits/unicef_excess/Data/annual_deaths_countries_selected_sources_young.csv")



# # variable characteristics by country and source 
# summ1 <- 
#   db3 %>% 
#   group_by(Country, Code, Year, Source, Sex) %>% 
#   summarise(Age_groups = n(),
#             Deaths = sum(Deaths)) %>% 
#   ungroup() %>% 
#   group_by(Country, Code, Year, Source, Age_groups) %>% 
#   summarise(Sex_groups = n(),
#             Deaths = sum(Deaths)) %>% 
#   ungroup() %>% 
#   group_by(Country, Code, Source, Age_groups, Sex_groups) %>% 
#   summarise(Period = paste(min(Year), max(Year), sep = "-"),
#          Deaths = sum(Deaths)) %>% 
#   ungroup() %>% 
#   group_by(Country) %>% 
#   mutate(n_sources = n()) %>% 
#   ungroup()
# 
# # variable configuration by source
# summ2 <- 
#   summ1 %>% 
#   select(Country, Code, Source, n_sources) %>% 
#   unique()
# 
# 
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
