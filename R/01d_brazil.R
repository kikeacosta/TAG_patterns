library(readr)
library(tidyverse)
library(janitor)

links <- paste0("Data/brazil/Mortalidade_Geral_", 2015:2021, ".csv")
i <- links[1]
out <- list()
for (i in links){
  cat(i)
  out[[i]] <- read_delim(i, 
                         delim = ";",
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
  bind_rows()

tot_age <- 
  dts %>% 
  group_by(Year, Sex) %>% 
  summarise(Deaths = sum(Deaths), .groups = "drop") %>% 
  mutate(Age = "TOT")

# rescaling age and sex
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

write_rds(dts3, "Output/brazil.rds")
