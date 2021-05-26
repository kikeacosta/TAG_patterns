
library(readr)
library(tidyverse)
library(janitor)
BR2020 <-
  read_delim("https://opendatasus.saude.gov.br/dataset/57c4b23e-6ffa-43ef-a1e7-23103f73f376/resource/da17c5f6-aa89-4e7d-b7e2-ec4b77a5dc31/download/dobrano_.csv", delim = ";") %>% 
  group_by(IDADE, SEXO) %>% 
  summarize(Deaths = n(), .groups = "drop") %>% 
  mutate(age_type = substr(IDADE, 1,1),
         Age = substr(IDADE, 2,3),
         Age = case_when(age_type %in% c("0","9")~"UNK",
                         is.na(age_type)~"UNK",
                         age_type %in% c("1","2","3")~"00",
                         age_type == "5" ~ "100",
                         TRUE ~ Age),
         Sex = recode(SEXO,
                      "0" = "UNK",
                      "1" = "m",
                      "2" = "f"))%>% 
  group_by(Sex, Age) %>% 
  summarize(Deaths = sum(Deaths), .groups = "drop") %>% 
  pivot_wider(names_from = Sex, values_from = Deaths) %>%
  mutate(UNK = ifelse(is.na(UNK),0, UNK),
         b=f+m+UNK) %>% 
  select(-UNK) %>% 
  adorn_totals("row",name = "TOT") %>% 
  dplyr::filter(Age != "UNK") %>% 
  pivot_longer(f:b, names_to = "Sex", values_to = "Deaths") %>% 
  mutate(Age = ifelse(Age == "TOT", "TOT", Age %>% as.integer() %>% as.character())) %>% 
  arrange(Sex, Age) %>% 
  mutate(Country = "Brazil",
         Code = "BRA",
         Year = 2020,
         Source = "opendatasus.saude.gov.br") %>% 
  select(Country, Code, Year, Sex, Age, Deaths, Source)

# now do other years one at a time, then rbind and save.
# will need to do manual, though, since code might have inconsistencies?
# especially IDADE might be coded differently.





