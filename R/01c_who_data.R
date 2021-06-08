library(readxl)
library(tidyverse)
library(countrycode)

# assorted helper functions written on the fly
has_total <- function(age_vec){
  "total age" %in% age_vec
}
has_not_total <- function(age_vec){
  any(age_vec != "total age")
}

has_annual <- function(time_units){
  any(time_units == "Annual")
}
is_complete <- function(time, time_unit){
  if (any(time_unit == "Annual")){
    return(TRUE)
  }
  if (any(time_unit == "Month")){
    return(all(1:12 %in% time))
  }
  if (any(time_unit == "Week")){
    return(all(1:52 %in% time))
  }
  NA_integer_
}
rescale_age <- function(chunk){
  TOT <- chunk %>% dplyr::filter(age_cat_s == "total age") %>% dplyr::pull(deaths)
  chunk %>% 
    dplyr::filter(age_cat_s != "total age") %>% 
    mutate(deaths = (deaths / sum(deaths)) * TOT)
}
age2int2 <- function(Age){
  OA <- 105
  OAvalue <- 105 - max(Age)
  DemoTools::age2int(Age, OAG = TRUE, OAvalue = OAvalue)
}
# -------------------------------------
# read in May 10 data
WHOin <- read_excel("Data/WHO_Allcause_Mortality_Data_Call_01.06.2021.xlsx")


# filter down to incl age and 2020
WHO_age_2020 <-
  WHOin %>% 
  group_by(country, year, time_unit, time, sex) %>%
  mutate(ageq = has_not_total(age_cat_s)) %>% 
  ungroup() %>% 
  dplyr::filter(ageq) %>% 
  group_by(country, sex) %>%
  mutate(has2020 = any(year == 2020)) %>% 
  ungroup() %>% 
  dplyr::filter(has2020,
                year > 2015,
                year < 2021)


# WHO_age_2020 %>% 
#   group_by(country, year, sex, time_unit) %>% 
#   mutate(consistent_age = unique(format_age) %>% length()) %>% 
#   ungroup() 


# if has annual, prefer it.
# if not, ensure that the year is complete (all weeks or months)
# aggregate within year and age

# First aggregate based on complete weeks or months,
# ignoring annual
WHO_age_2020_annual <-
  WHO_age_2020 %>% 
  dplyr::filter(time_unit != "Annual") %>% 
  group_by(country, year, sex, time_unit) %>% 
  mutate(compl = is_complete(time, time_unit)) %>% 
  dplyr::filter(compl) %>% 
  group_by(country, year, sex, time_unit, age_cat_s) %>% 
  summarize(deaths = sum(deaths), .groups = "drop") 

# get annual data, ignoring georgia
WHO_age_2020_annual_pre <-
  WHO_age_2020 %>% 
  dplyr::filter(time_unit == "annual") %>% 
  dplyr::select(all_of(colnames(WHO_age_2020_annual))) %>% 
  dplyr::filter(!(country == "GEO" & year == 2020)) %>% 
  distinct()

# which are the combos we have annual data for already?
check1 <- 
  WHO_age_2020_annual_pre %>% 
  select(country, sex, year) %>% 
  distinct() %>% 
  mutate(check = paste(country, sex, year)) %>% 
  dplyr::pull(check) %>% 
  unique()

# remove self-aggregated combos, bind on annual data
WHO_age_2020_annual <-
  WHO_age_2020_annual %>% 
  mutate(check2 = paste(country, sex, year)) %>% 
  dplyr::filter(!check2 %in% check1) %>% 
  bind_rows(WHO_age_2020_annual_pre) %>% 
  dplyr::select(-check2)

# 01-06 version fixes this.
# Fix Armenia 0 + 1-4 oddity
# WHO_age_2020_annual <- 
# WHO_age_2020_annual %>% 
#   filter(!(age_cat_s == "0" & country == "ARM")) %>% 
#   mutate(age_cat_s = ifelse(country == "ARM" & age_cat_s == "1-4","0",age_cat_s))


# This produces a auxiliary selector dataset
# consisting in just subsets where  2020 deaths > 2019 deaths
# discarded in end.

# WHO_selector_excess <-
#   WHO_age_2020_annual %>% 
#   mutate(time_unit = "annual") %>% 
#   filter(age_cat_s != "total age",
#          sex != "Unknown") %>% 
#   group_by(country, sex, year) %>% 
#   summarize(tot = sum(deaths),
#             .groups = "drop") %>% 
#   pivot_wider(names_from = year, values_from = tot) %>% 
#   mutate(excess = `2020` > `2019`,
#          excess = ifelse(is.na(excess),TRUE, excess)) %>% 
#   filter(excess,
#          `2020` > 2000) %>% 
#   pivot_longer(`2016`:`2020`, names_to = "year", values_to = "tot") %>% 
#   filter(!is.na(tot)) %>% 
#   select(-excess, -tot) %>% 
#   mutate(year = as.double(year))

# leave off with rescale operations.

# Filter to useful categories,
# then rescale to total age
# then scale sex-specific data to sum to stated totals
WHO_selection <-
  WHO_age_2020_annual %>% 
  filter(age_cat_s != "Unknown",
         sex != "Unknown") %>% 
  select(-time_unit) %>% 
  group_by(country, sex, year) %>% 
  do(rescale_age(chunk = .data)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = sex, values_from = deaths) %>% 
  mutate(PM = Male / (Male + Female),
         PF = Female / (Male + Female),
         PM = ifelse(is.nan(PM),0,PM),
         PF = ifelse(is.nan(PF),0,PF),
         Total = ifelse(is.na(Total), Male + Female, Total),
         Male = PM * Total,
         Female = PF * Total) %>% 
  select(-PM, -PF) %>% 
  pivot_longer(Female:Total, 
               values_to = "deaths", 
               names_to = "sex",
               values_drop_na = TRUE) %>% 
  dplyr::filter(age_cat_s != "total age")

# need to summarize years < 2020 as mean to deliver just two
# time points per population / sex.
WHO_2020_save <- 
  WHO_selection %>% 
  dplyr::filter(year == 2020)

# merge, then recode values and colnames to standards 
db_who <-
  WHO_selection %>% 
  mutate(Age = recode(age_cat_s,
                        "0"     = 0L,           
                        "1-4"   = 1L,
                        "10-14" = 10L,
                        "15-19" = 15L,
                        "20-24" = 20L,
                        "25-29" = 25L,
                        "30-34" = 30L,
                        "35-39" = 35L,
                        "40-44" = 40L,
                        "45-49" = 45L,
                        "5-9"   = 5L,
                        "50-54" = 50L,
                        "55-59" = 55L,
                        "60-64" = 60L,
                        "65-69" = 65L,
                        "70-74" = 70L,
                        "75-79" = 75L,
                        "80-84" = 80L,
                        "85+"   = 85L,
                        "0-44"  = 0L,
                        "45-64" = 45L,
                        "65-74" = 65L,
                        "75-84" = 75L,
                        "< 5"   = 0L,     
                        "85-89" = 85L,
                        "90+"   = 90L,
                        "0-30"  = 0L ,
                        "85-89" = 85L,
                        "90-94" = 90L,
                        "95+"   = 95L,
                        "0-4"   = 0L ,
                        "100+"  = 100L ,
                        "95-100"= 95L,
                        "0-64"  = 0L ,
                        "65-79" = 65L,
                        "80+"   = 80L,
                        "100+"  = 100L,
                      "85-90" = 85L),
         Sex = recode(sex, "Male" = "m","Female" = "f", "Total" = "t") ) %>% 
  dplyr::select(Code = country, Year = year, Sex, Age, Deaths = deaths) %>% 
  dplyr::filter(!Code %in% c("AND")) %>% 
  mutate(Country =countrycode(sourcevar = Code, 
                              origin = "iso3c", 
                              destination = "country.name"),
         Country = ifelse(Country == "United States","USA",Country),
         Source = "who") %>% 
  arrange(Country, Year, Sex, Age) %>% 
  ungroup() %>% 
  select(Country, Code, Year, Sex, Age, Deaths, Source)

# save out
write_csv(db_who, file = "Output/who.csv")  


# who_countries <- WHO_compare$Country %>% unique()
# 
# who_countries[!who_countries%in%offsets$Country]
# offsets$Country %>% unique()
# "Armenia"   "Georgia"   "Mauritius"

# offsets <- readRDS("Data/Offsets.rds")  %>% 
#   dplyr::filter(Region == "All")
# 
# WPP <- read_csv("Data/WPP2019_PopulationBySingleAgeSex_1950-2019.csv") %>% 
#   dplyr::filter(MidPeriod == 2018.5,
#                 Location %in% c("Armenia",   "Georgia",   "Mauritius"))
# offsets <-
# WPP %>% 
#   pivot_longer(PopMale:PopTotal, names_to = "Sex", values_to = "Population") %>% 
#   mutate(Sex = recode(Sex, 
#                       "PopMale" = "m",
#                       "PopFemale" = "f",
#                       "PopTotal" = "b"),
#          Region = "All") %>% 
#   select(Age = AgeGrpStart, Population, Country = Location, Region,  Sex) %>% 
#   mutate(Population = Population * 1000) %>% 
#   bind_rows(offsets) %>% 
#   arrange(Country, Region, Sex, Age)
# 
# 
# saveRDS(offsets, "Data/Offsets.rds")
#   

exploratory <- FALSE
if (exploratory){
WHOin %>% 
  filter(country == "ZMB",
         sex == "Total") %>% 
  group_by(year) %>% 
  summarize(d = sum(deaths))
WHOin %>% 
  filter(country == "ZMB",
         sex == "Total",
         time_unit !="annual") %>% 
  group_by(year) %>% 
  summarize(d = sum(deaths))


WHOin %>% 
  filter(country == "GTM",
         age_cat_s == "total age",
         sex == "Total") %>% 
  View()
WHOin %>% 
  pull(time_unit) %>% unique()


WHOin %>% 
  filter(time_unit == "annual") %>% 
  pull(country)%>% unique()

WHOin %>% 
  filter(country == "ECU",
         sex == "Total") %>% 
  ggplot(aes(x = time, y = deaths, group = year)) + 
  geom_line(  ) 

WHOin %>% 
  filter(country == "MUS",
         age_cat_s != "total age",
         sex == "Total",
         time_unit == "annual") %>% 
  group_by(year) %>% 
  summarize(D = sum(deaths))


WHOin %>% 
  filter(country == "PER",
         sex == "Total",
         year == 2020) %>% 
  pull(age_cat_s)
}