rm(list=ls())
source("R/00_functions.R")

# t <- search_eurostat("mort", type = "dataset", fixed = TRUE)
# t <- search_eurostat("eaths", type = "dataset", fixed = TRUE)

d_y <- 	get_eurostat("demo_magec")

unique(d_y$age)  

db_eurs <-
  d_y %>% 
  mutate(year = year(time)) %>% 
  # separate(time, sep = "W", into = c("year","week"), convert = TRUE) %>% 
  filter(year %in% 2010:2022, 
                age != "UNK") %>% 
  group_by(geo) %>% 
  filter(max(year) >= 2020) %>% 
  ungroup() %>% 
  mutate(age = case_when(age == "Y_LT1" ~ "0",
                         age == "Y_OPEN" ~ "100",
                         age == "TOTAL" ~ "TOT",
                         TRUE ~ str_sub(age, 2, 3)),
         sex = tolower(sex),
         Country = suppressWarnings(countrycode(sourcevar = geo, 
                               origin = "iso2c", 
                               destination = "country.name")),
         Country = case_when(geo == "EL" ~ "Greece",
                             geo == "UK" ~ "United Kingdom",
                             TRUE ~ Country),
         Code = countryname(Country, destination = "iso3c"),
         Source = "eurs_annual") %>% 
  arrange(Country, year, sex, suppressWarnings(as.integer(age))) %>% 
  dplyr::select(Country, Code, Year = year, Sex = sex, Age = age, Deaths = values, Source) %>% 
  drop_na()
  
# re-scaling ages and sexes
db_eurs2 <- 
  db_eurs %>% 
  group_by(Country, Sex, Year) %>% 
  do(rescale_age(chunk = .data)) %>% 
  ungroup() 

# chunk <- 
#   db_eurs %>% 
#   filter(Country == "Azerbaijan",
#          Year == 2010,
#          Sex == "f")


to_scale_sex <- 
  db_eurs %>% 
  filter(Sex != "t") %>% 
  group_by(Country, Year, Age) %>% 
  summarise(dts_sum = sum(Deaths)) %>% 
  ungroup() %>% 
  left_join(db_eurs %>% 
              filter(Sex == "t") %>% 
              rename(dts_tot = Deaths)) %>% 
  filter(dts_sum != dts_tot) %>% 
  select(Country, Year, Sex) %>% unique()

if(dim(to_scale_sex)[1] == 0){
  db_eurs3 <- 
    db_eurs2 %>% 
    mutate(Age = Age %>% as.double()) %>% 
    group_by(Country, Year, Sex) %>% 
    arrange(Age) %>% 
    mutate(age_spn = ifelse(Age == max(Age), -1, lead(Age) - Age)) %>% 
    ungroup() 
}else{
  db_eurs3 <- 
    db_eurs2 %>% 
    inner_join(to_scale_sex) %>% 
    group_by(Country, Age, Year) %>%
    do(rescale_sex(chunk = .data)) %>% 
    ungroup() %>% 
    bind_rows(db_eurs2 %>% 
                anti_join(to_scale_sex)) %>% 
    mutate(Age = Age %>% as.double()) %>% 
    group_by(Country, Year, Sex) %>% 
    arrange(Age) %>% 
    mutate(age_spn = ifelse(Age == max(Age), -1, lead(Age) - Age)) %>% 
    ungroup()
  
  }

diff_age_grs <- 
  db_eurs3 %>% 
  group_by(Country, Code, Year, Sex) %>% 
  filter(n() < 101) %>% 
  ungroup() %>% 
  select(Country, Sex) %>% unique()

db_eurs4  <-  
  db_eurs3 %>% 
  anti_join(diff_age_grs)

write_csv(db_eurs4, "data_inter/eurs_annual.csv")
db_eurs4 <- read_csv("data_inter/eurs_annual.csv")





  

