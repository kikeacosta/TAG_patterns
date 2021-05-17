library(readxl)
library(tidyverse)
library(countrycode)
WHOin <- read_excel("Data/WHO_Allcause_Mortality_Data_Call_10.05.2021.xlsx")

str(WHOin)
iso3 <- WHOin$country %>% unique()

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
  geom_line(  ) %>% 
  View()
  
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


# %>% 
 # ggplot(aes(x=time, y=deaths,group = year))+
 #  geom_line()

# LBN is out
# MDV out, ZMB out

# programmatic filters
# throw out anyone without age (OMN)
# throw out anyone without 2020 (LBN, ECU, PAN)

# visual filters
# throw out anyone where 2020 doesn't look complete
# ZMB, MDV, 




WHOin %>% 
  filter(country == "LBN")


WHOin %>% 
  filter(year == 2020,
         sex == "Total",
         substr(age_cat_s, 1, 1) == "0",
         !country %in% c("MUS","ZMB","MDV")) %>% 
  pull(country) %>% unique()

# merge w WPP for 2020
# aggregate years prior to 2020

# get sex, country, deaths_prior, expos_prior, deaths, expos.

# Question, does each subset reliably have "total age" ?
WHOin %>% 
  group_by(country, year, time_unit, time, sex) %>% 
  summarize(tots = has_total(age_cat_s), .groups = "drop") %>% 
  pull(tots) %>% sum()

has_total <- function(age_vec){
  "total age" %in% age_vec
}
has_not_total <- function(age_vec){
  any(age_vec != "total age")
}

# WHOin %>% 
#   group_by(country, year, time_unit, time, sex) %>% 
#   summarize(ageq = has_not_total(age_cat_s), .groups = "drop") %>% 
#   pull(ageq) %>% sum()
#   #nrow()

# filter down to incl age
WHO_age <- WHOin %>% 
  group_by(country, year, time_unit, time, sex) %>%
  mutate(ageq = has_not_total(age_cat_s)) %>% 
  ungroup() %>% 
  dplyr::filter(ageq)

# only subsets with useful years
WHO_age_2020 <- WHO_age %>% 
  group_by(country, sex) %>%
  mutate(has2020 = any(year == 2020)) %>% 
  ungroup() %>% 
  dplyr::filter(has2020,
                year > 2015,
                year < 2021)

# get annual totals by sex

WHO_age_2020 %>% 
  group_by(country, year, sex, time_unit) %>% 
  mutate(consistent_age = unique(format_age) %>% length()) %>% 
  ungroup() 


# if it has annual, throw out the rest :-)
has_annual <- function(time_units){
  any(time_units == "annual")
}
WHO_age_2020_anuual1 <-
  WHO_age_2020 %>% 
  group_by(country, year, sex) %>% 
  mutate(annual = has_annual(time_unit)) %>% 
  ungroup() %>% 
  filter(!(annual & time_unit != "annual")) %>% 
  select(-ageq, -has2020, -annual)
  

is_complete <- function(time, time_unit){
  if (any(time_unit == "annual")){
    return(TRUE)
  }
  if (any(time_unit == "month")){
    return(all(1:12 %in% time))
  }
  if (any(time_unit == "week")){
    return(all(1:52 %in% time))
  }
  NA_integer_
}

WHO_age_2020_complete <-
  WHO_age_2020_anuual1 %>% 
  group_by(country, year, sex, time_unit) %>% 
  mutate(compl = is_complete(time, time_unit)) %>% 
  filter(compl) 

WHO_age_2020_annual <-
  WHO_age_2020_complete %>% 
  group_by(country, year, sex, time_unit, age_cat_s) %>% 
  summarize(deaths = sum(deaths), .groups = "drop") 

# before this, ensure

# 2020 at least as many deaths as 2019

WHO_selector <-
  WHO_age_2020_annual %>% 
  mutate(time_unit = "annual") %>% 
  filter(age_cat_s != "total age",
         sex != "Unknown") %>% 
  group_by(country, sex, year) %>% 
  summarize(tot = sum(deaths),
            .groups = "drop") %>% 
  pivot_wider(names_from = year, values_from = tot) %>% 
  mutate(excess = `2020` > `2019`,
         excess = ifelse(is.na(excess),TRUE, excess)) %>% 
  filter(excess,
         `2020` > 2000) %>% 
  pivot_longer(`2016`:`2020`, names_to = "year", values_to = "tot") %>% 
  filter(!is.na(tot)) %>% 
  select(-excess, -tot) %>% 
  mutate(year = as.double(year))

WHO_selection <-
  WHO_selector %>% 
  left_join(WHO_age_2020_annual) %>% 
  filter(age_cat_s != "Unknown") %>% 
  select(-time_unit)


# leave off with rescale operations.
rescale_age <- function(chunk){
  
}

WHO_selection %>% 
  pull(age_cat_s) %>% unique()


