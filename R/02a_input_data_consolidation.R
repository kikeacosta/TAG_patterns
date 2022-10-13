rm(list=ls())
source("R/00_functions.R")

# loading pre-processed data
hmd  <- read_csv("Output/hmd.csv")
who  <- read_csv("Output/who.csv")
unpd <- read_csv("Output/unpd.csv")
stmf <- read_csv("Output/stmf.csv") 
eurs <- read_csv("Output/eurs.csv") 
bra <- read_csv("Output/brazil.csv") %>% mutate(Source = "direct")
per <- read_csv("Output/peru.csv") %>% mutate(Source = "direct")
mex <- read_csv("Output/mexico.csv") %>% mutate(Source = "direct")
zaf <- read_csv("Output/south_africa.csv") %>% mutate(Source = "direct")
irn <- read_csv("Output/iran.csv") %>% mutate(Source = "direct")
mse <- read_csv("Output/msemburi_tag.csv") 

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
  filter(max(Year) >= 2020 & min(Year) <= 2017)

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

ref_ages <- 
  all_in2 %>% 
  filter(Year < 2020) %>% 
  # no ages overpassing the minimum open age interval
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
  # identify and choose the largest age span in al periods for each age
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
  arrange(Country, Year, Sex, Age)

unique(out$Country)

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
  group_by(Country, Code, Source) %>% 
  summarise(Period = paste(min(Year), max(Year), sep = "-"),
            Deaths = sum(Deaths),
            # putting together all the age and sex groups configurations between 2015 and 2021
            Age_groups = unique(Age_groups) %>% paste(collapse = ", "),
            Sex_groups = unique(Sex_groups) %>% paste(collapse = ", ")) %>% 
  ungroup() 

out <- read_csv("Output/annual_deaths_countries_selected_sources.csv")
unique(out$Country)

test <- 
out %>% 
  select(Country, Year) %>% 
  unique() %>% 
  mutate(id = 1) %>% 
  spread(Year, id)
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
  
# temp1 <- 
#   all_in2 %>% 
#   filter(Year < 2020)
# 
# srs <- unique(temp1$Source)
# for(sr in srs){
#   temp2 <- 
#     temp1 %>% 
#     filter(Source == sr)
#   cts <- unique(temp2$Country)
#   for(ct in cts){
#     temp3 <- 
#       temp2 %>% 
#       filter(Country == ct)
#     sxs <- unique(temp3$Sex)
#     for(sx in sxs){
#       temp4 <- 
#         temp3 %>% 
#         filter(Sex == sx)
#       yrs <- unique(temp4$Year)
#       for(yr in yrs){
#         temp5 <- 
#           temp4 %>% 
#           filter(Year == yr)
#         
#         test <- harmon_age_intervals(temp5)
#       }
#     }
#   }
# }



