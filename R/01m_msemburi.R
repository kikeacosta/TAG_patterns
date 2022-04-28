
source("R/00_functions.R")

TAG <- 
  read_excel("Data/msemburi_inputs_28.04.2022.xlsx") %>% 
  select(Country, Code = iso3, Year = year, Sex = sex, Age = age, Deaths = val) %>% 
  mutate(Source = "Msemburi(TAG)",
         Sex = ifelse(Sex == "Female","f","m"))  %>%
  pivot_wider(names_from = Sex, values_from = Deaths) %>% 
  mutate(t = f + m) %>% 
  pivot_longer(c(f,m,t), names_to = "Sex", values_to = "Deaths") 
  

TAGt <- 
  TAG %>% 
  group_by(Country, )


# chr (5): Country, Code, Sex, Age, Source
# dbl (2): Year, Deaths