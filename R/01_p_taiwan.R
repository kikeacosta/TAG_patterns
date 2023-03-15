tw_url <- "https://ws.moi.gov.tw/001/Upload/400/relfile/0/4405/48349492-6f8c-453b-a9d1-4a8f0593b979/year/y02-05-1.xls"

library(rio)
library(tidyverse)
IN <- import(tw_url,range = "C5:DA95")
head(IN)
IN |> 
  mutate(year = rep(1992:2021,each=3), .before = 1) |> 
  filter(year >= 2010) |> 
  rename(sex = `...1`,
         tot = `Grand Total`) |> 
  mutate(sex = case_when(sex == "T."~3,
                         sex == "M."~1,
                         sex == "F."~2)) |> 
 mutate_if(is.character,as.double) |> 
 pivot_longer(-c(year,sex,tot), names_to = "age", values_to = "deaths") |> 
 mutate(age = gsub(age, pattern = "Years...",replacement = ""),
        age = case_when(age == "Over" ~ 100L,
                        TRUE ~ as.integer(age)-3))  |> 
  select(-tot) 
  

 
