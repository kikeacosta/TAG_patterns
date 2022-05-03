
source("R/00_functions.R")

TAG <- 
  read_excel("Data/msemburi_inputs_28.04.2022.xlsx") %>% 
  select(Country, Code = iso3, Year = year, Sex = sex, Age = age, Deaths = val) %>% 
  mutate(Source = "Msemburi(TAG)",
         Sex = ifelse(Sex == "Female","f","m"))  %>%
  pivot_wider(names_from = Sex, values_from = Deaths) %>% 
  mutate(t = f + m) %>% 
  pivot_longer(c(f,m,t), names_to = "Sex", values_to = "Deaths") %>% 
  mutate(Age = as.character(Age))
  

TAGt <- 
  TAG %>% 
  group_by(Country,  Code,   Year, Source, Sex) %>% 
  summarize(Deaths = sum(Deaths), .groups= "drop") %>% 
  mutate(Age = "TOT")

TAGout <- bind_rows(TAG, TAGt) %>% 
  arrange(Code,Year,Sex,Age) %>% 
  relocate(Source, .after = last_col())

# ------------------------------------ #
WHOold <- read_excel("Data/WHO_Allcause_Mortality_Data_Call_13.07.2021.xlsx")

# inherited custom functions 
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

WHO_age_pre <-
  WHOold %>% 
  group_by(country, year, time_unit, time, sex) %>%
  mutate(ageq = has_not_total(age_cat_s)) %>% 
  ungroup() %>% 
  dplyr::filter(ageq) %>% 
  group_by(country, sex) %>%
  ungroup() %>% 
  dplyr::filter(year > 2015,
                year < 2020)
# First aggregate based on complete weeks or months,
# ignoring annual
WHO_age_pre_annual <-
  WHO_age_pre %>% 
  dplyr::filter(time_unit != "Annual") %>% 
  group_by(country, year, sex, time_unit) %>% 
  mutate(compl = is_complete(time, time_unit)) %>% 
  dplyr::filter(compl) %>% 
  group_by(country, year, sex, time_unit, age_cat_s) %>% 
  summarize(deaths = sum(deaths), .groups = "drop") 

# get annual data, ignoring georgia
WHO_age_pre_annual2 <-
  WHO_age_pre %>% 
  dplyr::filter(time_unit == "Annual") %>% 
  dplyr::select(all_of(colnames(WHO_age_pre_annual))) %>% 
  dplyr::filter(country != "GEO") %>% 
  distinct()

# which are the combos we have annual data for already?
check1 <- 
  WHO_age_pre_annual2 %>% 
  select(country, sex, year) %>% 
  distinct() %>% 
  mutate(check = paste(country, sex, year)) %>% 
  dplyr::pull(check) %>% 
  unique()

# remove self-aggregated combos, bind on annual data
WHO_age_pre_annual3 <-
  WHO_age_pre_annual %>% 
  mutate(check2 = paste(country, sex, year)) %>% 
  dplyr::filter(!check2 %in% check1) %>% 
  bind_rows(WHO_age_pre_annual2) %>% 
  dplyr::select(-check2)

WHO_selection <-
  WHO_age_pre_annual3 %>% 
  filter(age_cat_s != "Unknown",
         sex != "Unknown") %>% 
  group_by(country, time_unit, sex, year) %>% 
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


WHO_selection_2 <- 
  WHO_selection %>% 
  group_by(country, sex, year) %>% 
  mutate(priority = case_when(time_unit == "Annual" ~ 1,
                              time_unit == "Month" ~ 2,
                              time_unit == "Week" ~ 3)) %>% 
  dplyr::filter(priority == min(priority)) %>% 
  ungroup() %>% 
  select(-priority) 

db_who <-
  WHO_selection_2 %>% 
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
                      "85-90" = 85L,
                      "15-30" = 15L,
                      "0-15" = 0L,
                      "95-99" = 95L),
         Sex = recode(sex, "Male" = "m","Female" = "f", "Total" = "t") ) %>% 
  dplyr::select(Code = country, Year = year, Sex, Age, Deaths = deaths) %>% 
  dplyr::filter(!Code %in% c("AND")) %>% 
  mutate(Country =countrycode(sourcevar = Code, 
                              origin = "iso3c", 
                              destination = "country.name"),
         Country = ifelse(Country == "United States","United States of America",Country),
         Source = "Msemburi(TAG)") %>% 
  arrange(Country, Year, Sex, Age) %>% 
  ungroup() %>% 
  select(Country, Code, Year, Sex, Age, Deaths, Source)

# merge w William's inputs for just 2020-2021
tag_countries <- TAGout$Country %>% unique()
who_countries <- db_who$Country %>% unique()

all(who_countries %in% tag_countries)
who_countries[!who_countries %in% tag_countries]
tag_countries[!tag_countries %in% who_countries]

who_totals <- 
  db_who %>% 
  group_by(Country, Code, Year, Sex,Source) %>% 
  summarize(Deaths = sum(Deaths),
            Age = "TOT", .groups = "drop")

colnames(db_who)


db_msemburi_tag_who <- 
  db_who %>% 
  dplyr::filter(Country %in% tag_countries) %>% 
  mutate(Source = "Msemburi(TAG)",
         Age = as.character(Age)) %>% 
  bind_rows(who_totals) %>% 
  bind_rows(TAGout)
# save out
write_csv(db_msemburi_tag_who, file = "Output/msemburi_tag.csv")  

db_msemburi_tag_who$Code %>% unique()
