source("R/00_functions.R")

d1 <- 
  read_csv("Data/ariel_age_sex_data_csv.csv")

unique(d1$age)
unique(d1$sex)




sum <-
  d1 %>%
  select(Country = country_name, Sex = sex, Year = year) %>%
  unique() %>%
  mutate(Sex = str_sub(Sex, 1, 1) %>% str_to_lower())

# ar2 <- 
#   d1 %>% 
#   select(Country = country_name, 
#          Sex = sex, 
#          Age = age,
#          Year = year, D_ariel = deaths)
  


our <- read_csv("Output/annual_deaths_countries_selected_sources.csv")

# # missing from Ariel
# # ~~~~~~~~~~~~~~~~~~
# our2 <- 
#   our %>% 
#   select(Country, Sex, Year) %>% 
#   anti_join(sum)
# 
# 
mis <-
  sum %>%
  anti_join(our %>%
                select(Country, Sex, Year))
#   
# 
# comp <- 
#   our %>% 
#   semi_join(sum) %>% 
#   left_join(ar2)


ar3 <- 
  d1 %>% 
  select(Country = country_name, 
         Sex = sex, 
         Age = age,
         Year = year, 
         Deaths = deaths) %>% 
  mutate(Sex = str_sub(Sex, 1, 1) %>% str_to_lower())

# %>% 
#   inner_join(mis)
  
unique(ar3$Age)



ar4 <- 
  ar3 %>% 
  mutate(age1 = Age) %>% 
  separate(age1, c("age3", "age4")) %>% 
  mutate(Country = ifelse(Country == "Bosnia", "Bosnia and Herzegovina", Country)) %>% 
  mutate(Deaths = Deaths %>% as.numeric(),
         age3 = case_when(Country == "Thailand" ~ str_sub(Age, 6,8),
                          Country == "Colombia" ~ str_sub(Age, 3,5),
                          Country == "Argentina" ~ str_sub(Age, 7,8),
                          TRUE ~ age3),
         age4 = case_when(age3 %in% c("Month", "Week", "Under", "Below", "Menos",
                                      "Menores", "MENORES", "Menor", "Up", "nor") ~ "0",
                          age3 %in% c("Unknown", "Not", "Edad", "Ignorado", 
                                      "unknown", "No", "ad ", "SIN") ~ "unk",
                          age3 %in% c("Above") ~ "85",
                          TRUE ~ age3),
         age4 = str_replace(age4, " ป", "") %>% str_trim,
         age4 = case_when(Age %in% c("不詳", "Неизвестно", "08.sin especificar") ~ "unk", 
                          Age %in% c("< Week", "<1", "1 >", "01.menor a 20 anios") ~ "0",
                          str_detect(age3, "อย") ~ "0", 
                          str_detect(age3, "กก") ~ "101",
                          TRUE ~ age4)) %>% 
  filter(Country != "Singapore")

unique(ar4$age4)


art <- 
  ar4 %>% 
  group_by(Country, Sex, Year) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Age = "TOT")

ar5 <- 
  ar4 %>% 
  filter(age4 != "unk") %>% 
  mutate(age4 = age4 %>% as.numeric(),
         age4 = ifelse(age4 > 100, 100, age4)) %>% 
  group_by(Country, Sex, age4, Year) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  rename(Age = age4) %>% 
  mutate(Age = Age %>% as.character())

ar6 <- 
  ar5 %>% 
  bind_rows(art) %>% 
  arrange(Country, Sex, Year, suppressWarnings(as.integer(Age))) %>% 
  group_by(Country, Sex, Year) %>% 
  do(rescale_age(chunk = .data)) %>% 
  ungroup() 

ar_ts <- 
  ar6 %>% 
  group_by(Country, Age, Year) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Sex = "t")

out <- 
  ar6 %>% 
  bind_rows(ar_ts) %>% 
  mutate(Age = Age %>% as.double()) %>% 
  arrange(Country, Sex, Year, Age) 


comp <- 
  our %>% 
  semi_join(out %>% select(Country, Sex, Age, Year) %>% unique) %>% 
  left_join(out %>% rename(D_ar = Deaths))



