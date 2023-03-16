tw_url <- "https://ws.moi.gov.tw/001/Upload/400/relfile/0/4405/48349492-6f8c-453b-a9d1-4a8f0593b979/year/y02-05-1.xls"

library(rio)
library(tidyverse)
IN <- import(tw_url,range = "C5:DA95")
head(IN)

dt <- 
  IN %>%  
  as_tibble() %>% 
  mutate(Year = rep(1992:2021,each=3), .before = 1) |> 
  filter(Year >= 2010) |> 
  rename(Sex = `...1`,
         tot = `Grand Total`) |> 
  mutate(Sex = case_when(Sex == "T."~3,
                         Sex == "M."~1,
                         Sex == "F."~2)) |> 
 mutate_if(is.character,as.double) |>
 pivot_longer(-c(Year,Sex,tot), names_to = "Age", values_to = "Deaths") |> 
 mutate(Age = gsub(Age, pattern = "Years...",replacement = ""),
        Age = case_when(Age == "Over" ~ 100,
                        TRUE ~ as.integer(Age)-3))  |> 
  select(-tot) 

dt2 <- 
  dt %>% 
  mutate(Sex = recode(Sex,
                      "2" = "f",
                      "1" = "m",
                      "3" = "t")) %>% 
  group_by(Year, Sex) %>% 
  mutate(age_spn = ifelse(Age == max(Age), -1, lead(Age) - Age)) %>% 
  ungroup() %>% 
  mutate(Code = "TWN",
         Country = "Taiwan",
         Source = "direct")

write_csv(dt2, "data_inter/taiwan.csv")


  

 
