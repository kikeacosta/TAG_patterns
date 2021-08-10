library(here)
source(here("R", "00_functions.R"))
# Mexico mortality data
# ~~~~~~~~~~~~~~~~~~~~~
# # files from OSF (Version 1) as of 15 March 2021 
# osf_retrieve_file("hbxkn") %>%
#   osf_download(conflicts = "overwrite",
#                path = "Data")
# 
mex_files <- unzip(here("Data", "Mexico", "mexico_deaths.zip"), list = TRUE)

# # all deaths from years 2020 and 2021
db_mx20 <- 
  read_csv(here("Data", "Mexico", "DDAAxsom2021SE14.csv"))

db_mx20_2 <- 
  db_mx20 %>% 
  select(Date = 5,
         Sexo = SEXO,
         Age = EDAD) %>% 
  mutate(Year = year(Date),
         # provisionally exchanging variable sex as its value was originally 
         # inverted
         Sex = case_when(Sexo == 1 ~ 2,
                         Sexo == 2 ~ 1,
                         TRUE ~ 3)) %>% 
  group_by(Year, Sex, Age) %>% 
  summarise(Deaths = n()) %>% 
  ungroup() 

# all deaths from years 2016-2019
db_mx12_19 <- tibble()
for(i in 1:8){
  temp <- 
    read_csv(unz(here("Data", "Mexico", "mexico_deaths.zip"), mex_files[i,1]))
  
  temp2 <- 
    temp %>% 
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
  
  db_mx12_19 <- 
    db_mx12_19 %>% 
    bind_rows(temp2)
}

# all Mexico deaths together
db_mx <- 
  db_mx12_19 %>% 
  select(Year, Sex, Age, Deaths) %>% 
  bind_rows(db_mx20_2) %>% 
  # drop_na(Date) %>%
  filter(Year >= 2016 & Year <= 2020) %>% 
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
  arrange(Country, Year, Sex, suppressWarnings(as.numeric(Age))) %>% 
  mutate(Code = "MEX", Source = "country_public")

# saving annual deaths in Mexico
write_csv(db_mx2, "Output/mexico.csv")



db_mx2 %>% 
  filter(Year == 2020,
         Age != "TOT") %>% 
  mutate(Age = as.integer(Age)) %>% 
  ggplot()+
  geom_line(aes(Age, Deaths, col = Sex, group = Sex))+
  scale_y_log10()

