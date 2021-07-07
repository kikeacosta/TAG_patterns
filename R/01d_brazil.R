
library(readr)
library(tidyverse)
library(janitor)

# doing this in one pass means that the input microdata doesn't eat up memory
# the only thing that gets saves is the tiny resulting table
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
         t=f+m+UNK) %>% 
  select(-UNK) %>% 
  adorn_totals("row",name = "TOT") %>% 
  dplyr::filter(Age != "UNK") %>% 
  pivot_longer(f:t, names_to = "Sex", values_to = "Deaths") %>% 
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


link2019 <- "https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SIM/Mortalidade_Geral_2019.csv"
link2018 <- "https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SIM/Mortalidade_Geral_2018.csv"
link2017 <- "https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SIM/Mortalidade_Geral_2017.csv"
link2016 <- "https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SIM/Mortalidade_Geral_2016.csv"
link2015 <- "https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SIM/Mortalidade_Geral_2015.csv"

years <- 2015:2019
out <- list()
for (i in 1:length(years)){
  cat(i,"\n")
  link <- paste0("https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SIM/Mortalidade_Geral_",years[i],".csv")
  out[[i]] <- read_delim(link, delim = ";") %>% 
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
           t=f+m+UNK) %>% 
    select(-UNK) %>% 
    adorn_totals("row",name = "TOT") %>% 
    dplyr::filter(Age != "UNK") %>% 
    pivot_longer(f:t, names_to = "Sex", values_to = "Deaths") %>% 
    mutate(Age = ifelse(Age == "TOT", "TOT", Age %>% as.integer() %>% as.character())) %>% 
    arrange(Sex, Age) %>% 
    mutate(Country = "Brazil",
           Code = "BRA",
           Year = years[i],
           Source = "country_public") %>% 
    select(Country, Code, Year, Sex, Age, Deaths, Source)
    
}

db_bra <-
  out %>% 
  bind_rows() %>% 
  bind_rows(BR2020) %>% 
  arrange(Year, Sex, as.integer(Age)) 
  
write_csv(db_bra, "Output/brazil.csv")

# db_bra %>% 
#   dplyr::filter(Age != "TOT", Sex != "t") %>% 
#   mutate(Age = as.integer(Age)) %>% 
#   ggplot(aes(x = Age, y = Deaths, color = Sex, group = interaction(Year, Sex))) + 
#   geom_line()



