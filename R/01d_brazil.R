
library(readr)
library(tidyverse)
library(janitor)

# doing this in one pass means that the input microdata doesn't eat up memory
# the only thing that gets saves is the tiny resulting table

# https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/SIM/sim_preliminar_2020.csv
# https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/SIM/DO21OPEN.csv

# links <- paste0(here("Data", "Brazil", "deaths"), "/Mortalidade_Geral_", 2015:2021, ".csv")
# 
# out <- list()
# for (i in links){
#   cat(i)
#   out[[i]] <- read_delim(i, 
#                          delim = ";",
#                          col_types = cols(.default = "c")) %>%  
#     filter(TIPOBITO == "2") %>% 
#     select(date_e = DTOBITO, mun_code = CODMUNOCOR) %>% 
#     mutate(date_e = dmy(date_e),
#            mth = month(date_e),
#            year = year(date_e),
#            date = make_date(d = 15, m = mth, y = year)) %>% 
#     group_by(date, mun_code) %>% 
#     summarise(dts = n()) %>% 
#     ungroup()
# }
# 
# dts <- 
#   out %>% 
#   bind_rows()






# 
# BR2020 <-
#   read_delim("https://opendatasus.saude.gov.br/dataset/57c4b23e-6ffa-43ef-a1e7-23103f73f376/resource/8d947ac1-addb-49f2-85ab-824a7408a432/download/dobrano_.csv", delim = ";") %>%
#   filter(TIPOBITO == 2) %>%
#   group_by(IDADE, SEXO) %>%
#   summarize(Deaths = n(), .groups = "drop") %>%
#   mutate(Age = case_when(IDADE < 400 ~ "0",
#                          IDADE > 400 & IDADE < 500 ~ str_sub(IDADE, 2, 3),
#                          IDADE >= 500 & IDADE <= 600 ~ "100",
#                          TRUE ~ "UNK"),
#          Sex = recode(SEXO,
#                       "0" = "UNK",
#                       "1" = "m",
#                       "2" = "f"))%>%
#   group_by(Sex, Age) %>%
#   summarize(Deaths = sum(Deaths), .groups = "drop") %>%
#   pivot_wider(names_from = Sex, values_from = Deaths) %>%
#   mutate(UNK = ifelse(is.na(UNK),0, UNK),
#          t=f+m+UNK) %>%
#   select(-UNK) %>%
#   adorn_totals("row",name = "TOT") %>%
#   dplyr::filter(Age != "UNK") %>%
#   pivot_longer(f:t, names_to = "Sex", values_to = "Deaths") %>%
#   mutate(Age = ifelse(Age == "TOT", "TOT", Age %>% as.integer() %>% as.character())) %>%
#   arrange(Sex, Age) %>%
#   mutate(Country = "Brazil",
#          Code = "BRA",
#          Year = 2020,
#          Source = "country_public") %>%
#   select(Country, Code, Year, Sex, Age, Deaths, Source)

# now do other years one at a time, then rbind and save.
# will need to do manual, though, since code might have inconsistencies?
# especially IDADE might be coded differently.
links <- tibble(year = 2021:2015, link = c("https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/SIM/DO21OPEN.csv",
                                          "https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/SIM/sim_preliminar_2020.csv",
                                          "https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SIM/Mortalidade_Geral_2019.csv",
                                          "https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SIM/Mortalidade_Geral_2018.csv",
                                          "https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SIM/Mortalidade_Geral_2017.csv",
                                          "https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SIM/Mortalidade_Geral_2016.csv",
                                          "https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SIM/Mortalidade_Geral_2015.csv"))
test <- read_delim("https://s3-sa-east-1.amazonaws.com/ckan.saude.gov.br/SIM/Mortalidade_Geral_2015.csv", delim = ";")
years <- 2015:2021
out <- list()
for (i in length(links$year)){
  y <- i
  cat(i,"\n")
  out[[i]] <- read_delim(links$link[i], delim = ";") %>% 
    group_by(IDADE, SEXO) %>% 
    summarize(Deaths = n(), .groups = "drop") %>% 
    mutate(Age = case_when(IDADE < 400 ~ "0",
                           IDADE > 400 & IDADE < 500 ~ str_sub(IDADE, 2, 3),
                           IDADE >= 500 & IDADE <= 600 ~ "100",
                           TRUE ~ "UNK"),
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
           Year = links$year[i],
           Source = "country_public") %>% 
    select(Country, Code, Year, Sex, Age, Deaths, Source)
    
}

db_bra <-
  out %>% 
  bind_rows() %>% 
  bind_rows(BR2020) %>% 
  arrange(Year, Sex, suppressWarnings(as.integer(Age)))
  
write_csv(db_bra, "Output/brazil.csv")

# db_bra %>% 
#   dplyr::filter(Age != "TOT", Sex != "t") %>% 
#   mutate(Age = as.integer(Age)) %>% 
#   ggplot(aes(x = Age, y = Deaths, color = Sex, group = interaction(Year, Sex))) + 
#   geom_line()



