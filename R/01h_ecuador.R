library(here)
source(here("R", "00_functions.R"))
library(lubridate)

# Ecuador mortality data
# ~~~~~~~~~~~~~~~~~~~~~
# # files from Ecuador (Version 1) as of 15 June 2021 
# 
ecu_files <- unzip(here("Data", "Ecuador", "BBD_EDG_2020_CSV_v1.zip"), list = TRUE)

db_2020 <- 
  read_csv2(unz(here("Data", "Ecuador", "BBD_EDG_2020_CSV_v1.zip"), ecu_files[1,1]))

# According to the Indice, file 1.2.6.csv:
# Número de defunciones en el año (t+1),por sexo, según grupos 
# de edad (quinquenales), según región y provincia de residencia 
# habitual de la persona fallecida.Periodo 2008 - 2020

# database in wide format, building variables names
names_grid <- 
  expand_grid(Year = 2008:2020, Sex = c("m", "f"), Age = c(0, 1, seq(5, 65, 5), "UNK")) %>% 
  mutate(names = paste(Sex, Age, Year, sep = "_"))

var_names <- c("Region", names_grid$names, "noth")

# loading data
db_ec <- 
  read_csv(here("Data", "Ecuador", "1.2.6.csv"),
           skip = 7,
           col_names = var_names,
           col_types = cols(.default = "c"))

# filtering only national level, naming variables, going to tidy format
db_ec2 <- 
  db_ec %>% 
  filter(Region == "Total Nacional:") %>% 
  select(-noth, -Region) %>% 
  gather(key = vars, value = Deaths) %>% 
  separate(vars, c("Sex", "Age", "Year"), sep = "_") %>% 
  mutate(Deaths = str_replace(Deaths, ",", "") %>% as.integer(),
         Year = Year %>% as.integer())


# for single-year ages without sex, according to the Indice, file 1.1.2D:
# "Número de defunciones en el año (t+1) por edades simples (años) a nivel 
# nacional. Periodo 1990 - 2020"

var_names_simple <- c("Year", "TOT", paste("a_", 1:120), "UNK", "noth", "noth2", "noth3") 
db_ec_simp <- 
  read_csv(here("Data", "Ecuador", "1.1.2D.csv"),
           skip = 4,
           col_names = var_names_simple,
           col_types = cols(.default = "c"))
           
db_ec_simp2 <- 
  db_ec_simp %>% 
  drop_na(TOT) %>% 
  select(-noth, -noth2, -noth3) %>% 
  gather(-Year, key = Age, value = Deaths) %>% 
  mutate(Deaths = str_replace(Deaths, ",", ""),
         Deaths = ifelse(Deaths == "-", "0", Deaths),
         Deaths = Deaths %>% as.integer(),
         Age = str_replace(Age, "a_", ""),
         Year = ifelse(str_detect(Year, "2020") , "2020", Year) %>% as.integer())

# ~~~~~~~~~~~~~~~~
# Consistency test
# ~~~~~~~~~~~~~~~~
# Age 0 is excluded from the simple age file
# testing that this is the case combining data from both datasets
# (it is the case)

db_ec_simp_t <- 
  db_ec_simp2 %>% 
  filter(Age == "TOT")

db_ec_simp_no_t <- 
  db_ec_simp2 %>% 
  filter(Age != "TOT") %>% 
  group_by(Year) %>% 
  summarise(Deaths_s = sum(Deaths))

db_ec_t <- 
  db_ec2 %>% 
  group_by(Year) %>% 
  summarise(Deaths_s2 = sum(Deaths))

db_ec_0 <- 
  db_ec2 %>% 
  filter(Age == "0") %>% 
  group_by(Year) %>% 
  summarise(Deaths_0 = sum(Deaths))

# testing that: Total deaths in the 5-year age group file = 
# total age in single-age file + ages 0-1 in the 5-year age group file  
# answer: yes
db_test <- 
  db_ec_simp_t %>% 
  left_join(db_ec_simp_no_t) %>% 
  left_join(db_ec_t) %>% 
  left_join(db_ec_0) %>% 
  drop_na() %>% 
  mutate(tot = Deaths_s + Deaths_0)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# merging data to have single year of age between 0 and 100+
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

db_0 <- 
  db_ec2 %>% 
  filter(Age == "0") %>% 
  group_by(Year) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Age = "0")

db_single <- 
  db_ec_simp2 %>% 
  mutate(Age = str_trim(Age)) %>% 
  filter(Age != "TOT",
         Age != "UNK") %>% 
  bind_rows(db_0) %>% 
  mutate(Age = Age %>% as.integer(),
         Age = ifelse(Age > 100, 100, Age)) %>% 
  group_by(Year, Age) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Age = Age %>% as.character())

db_unk <- 
  db_ec_simp2 %>% 
  filter(Age == "UNK")

unique(db_single$Age)

db_out <- 
  db_single %>% 
  bind_rows(db_unk) %>% 
  group_by(Year) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Age = "TOT") %>% 
  bind_rows(db_single) %>% 
  filter(Year >= 2015) %>% 
  mutate(Country = "Ecuador",
         Source = "country_public",
         Sex = "t",
         Code = "ECU") %>% 
  select(Country, Code, Year, Sex, Age, Deaths, Source) %>% 
  arrange(Year, suppressWarnings(as.numeric(Age)))
  
write_csv(db_out, "Output/ecuador.csv")


  

