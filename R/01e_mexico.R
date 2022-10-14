source("R/00_functions.R")
# Mexico mortality data
# ~~~~~~~~~~~~~~~~~~~~~
# data from INEGI
# https://www.inegi.org.mx/programas/mortalidad/#Microdatos

mex_files <- unzip("Data/mexico/muertes_2015-2020.zip", list = TRUE)

# all deaths from years 2016-2019
db_mx15_20 <- tibble()
for(i in 1:6){
  temp <- 
    read_csv(unz("Data/mexico/muertes_2015-2020.zip", mex_files[i,1]))
  
  temp2 <- 
    temp %>% 
    rename_all(tolower) %>% 
    select(Sex = sexo, 
           Age = edad, 
           Year = anio_ocur) %>% 
    mutate(Age = case_when(Age < 4000 ~ 0,
                           Age > 4000 & Age < 4130 ~ Age - 4000,
                           TRUE ~ NA_real_)) %>% 
    group_by(Year, Sex, Age) %>% 
    summarise(Deaths = n()) %>% 
    ungroup() %>% 
    mutate(file_orig = mex_files[i,1])
  
  db_mx15_20 <- 
    db_mx15_20 %>% 
    bind_rows(temp2)
}

# all Mexico deaths together
db_mx <- 
  db_mx15_20 %>% 
  select(Year, Sex, Age, Deaths) %>% 
  # drop_na(Date) %>%
  filter(Year >= 2015 & Year <= 2020) %>% 
  arrange(Year, Sex, Age) %>% 
  mutate(Age = case_when(Age > 100 & Age <= 130 ~ "100",
                               is.na(Age) | Age > 130 ~ "UNK", 
                               TRUE ~ as.character(Age)),
         Sex = case_when(Sex == 1 ~ "m",
                         Sex == 2 ~ "f",
                         TRUE ~ "UNK"),
         Country = "Mexico") %>% 
  group_by(Country, Year, Sex, Age) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  select(Country, Year, Sex, Age, Deaths)

# all sex and all age Peru and Mexico
db_mx_all_ages <- 
  db_mx %>% 
  group_by(Country, Year, Sex) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Age = "TOT") %>% 
  filter(Sex != "UNK")

db_mx_all_sex <- 
  db_mx %>% 
  group_by(Country, Year, Age) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Sex = "t") %>% 
  filter(Age != "UNK")

db_mx_all_sex_age <- 
  db_mx %>% 
  group_by(Country, Year) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Sex = "t",
         Age = "TOT")

db_mx2 <- 
  db_mx %>% 
  filter(Age != "UNK",
         Sex != "UNK") %>% 
  bind_rows(db_mx_all_ages,
            db_mx_all_sex,
            db_mx_all_sex_age) %>% 
  group_by(Sex, Year) %>% 
  do(rescale_age(chunk = .data)) %>% 
  ungroup() %>%
  group_by(Age, Year) %>%
  do(rescale_sex(chunk = .data)) %>%
  ungroup()

db_mx_out <- 
  db_mx2 %>% 
  arrange(Country, Year, Sex, suppressWarnings(as.numeric(Age))) %>% 
  mutate(Age = Age %>% as.double(),
         Code = "MEX", 
         Source = "country_public")

# saving annual deaths in Mexico
write_csv(db_mx_out, "Output/mexico.csv")



