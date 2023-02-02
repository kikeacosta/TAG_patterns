library(readr)
library(tidyverse)
library(janitor)
source("R/00_functions.R")

# files downloaded from Datasus, the source is the Ministry of Health
# data 2015-2020
# https://opendatasus.saude.gov.br/dataset/sim-1979-2019
# data 2021
# https://dados.gov.br/dataset/sistema-de-informacao-sobre-mortalidade


links <- paste0("Data/brazil/Mortalidade_Geral_", 2015:2022, ".csv")
i <- links[8]
out <- list()
for (i in links){
  cat(i)
  
  separ <- ifelse(str_detect(i, "2021|2022"), ",", ";")
  
  out[[i]] <- read_delim(i, 
                         delim = separ,
                         col_types = cols(.default = "c")) %>%  
    filter(TIPOBITO == "2") %>% 
    select(date_e = DTOBITO, 
           IDADE,
           SEXO) %>% 
    mutate(date_e = dmy(date_e),
           Year = year(date_e),
           Sex = recode(SEXO,
                        "0" = "UNK",
                        "1" = "m",
                        "2" = "f"),
           Age = case_when(IDADE <= 400 ~ "0",
                           IDADE > 400 & IDADE < 500 ~ str_sub(IDADE, 2, 3),
                           IDADE >= 500 & IDADE <= 600 ~ "100",
                           TRUE ~ "UNK")) %>% 
    group_by(Year, Sex, Age) %>% 
    summarize(Deaths = n(), .groups = "drop") %>% 
    spread(Sex, Deaths) %>% 
    replace_na(list(UNK = 0,
                    m = 0,
                    f = 0)) %>% 
    mutate(t = f + m + UNK) %>% 
    select(-UNK) %>% 
    gather(f, m, t, key = Sex, value = Deaths) %>% 
    group_by(Year, Sex, Age) %>% 
    summarize(Deaths = sum(Deaths), .groups = "drop")
}

dts <- 
  out %>% 
  bind_rows() %>% 
  ungroup() %>% 
  filter(Year <= 2021)

tot_age <- 
  dts %>% 
  group_by(Year, Sex) %>% 
  summarise(Deaths = sum(Deaths), .groups = "drop") %>% 
  mutate(Age = "TOT")

# re-scaling age and sex
dts2 <- 
  dts %>% 
  filter(Age != "UNK") %>% 
  bind_rows(tot_age) %>% 
  group_by(Sex, Year) %>% 
  do(rescale_age(chunk = .data)) %>% 
  ungroup() %>%
  group_by(Age, Year) %>%
  do(rescale_sex(chunk = .data)) %>%
  ungroup()

dts3 <- 
  dts2 %>% 
  mutate(Age = Age %>% as.double(),
         Country = "Brazil",
         Code = "BRA",
         Source = "brazil_sim") %>% 
  select(Country, Code, Year, Sex, Age, Deaths, Source) %>% 
  arrange(Year, Sex, Age)

# write_rds(dts3, "data_inter/brazil.rds")
# dts3 <- read_rds("Output/brazil.rds")
write_csv(dts3, "data_inter/brazil.csv")
