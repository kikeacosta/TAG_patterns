
out <- read_csv("Output/annual_deaths_countries_selected_sources.csv")


codes <- out$Code %>% unique()

code <- codes[36]
out %>% 
  filter(Code == code) %>% 
  ggplot(aes(x = Age, y = Deaths)) +
  geom_col() +
  facet_grid(vars(Year),vars(Sex)) +
  labs(title = code)


# AZE lower open age in 2021 only
# BLZ (Belize) low open age 65+
# BOL lower open age in 2020-2021
# BRA: lower open age in 2020 only
# CAN lowre open ages in 2020-2021
# CUB higher open age in 2021 only
# DOM higher open age in 2021 only
# EST tot crazy coincidence in 2021 ages 80 and 90
# PYF (French Polynesia) is small
# GEO (Georgia) possibly irregular in older ages
# DEU low open age in 2020-21
# GTM flat deaths dist, needs denoms to evaluate
# ISR low open age in 2021 
# KEN prepandemic data has very wide intervals and low open age
# LSO low open age 75+ has large fraction of deaths
# LIE (Liechtenstein) small pop
# MDA has higher open age only in 2020
# MCO small pop, high fraction in open age
# MNG (Mongolia) very low open age pre pandemic
# NZL lower open age in 2021 only
# MKD (North Macedonia): small population, and higher deaths in 2015-2018 ?
# OMN (Oman) lower open age pre-pandemic
# PRY higher open age only in 2021
# PER higher open age only in 2021
# QAT has low open age group in 2020
# SMR (San Marino) very small pop and low open ages in 2020-2021
# SYC (Seychelles) higher open age in 2020 only
# SGP (Singapore) has a lower open age in 2020
# KOR has lower open ages in 2020-2021
# SUR (Suriname) has lower open ages in 2020-21
# TUN has lower open age in 2020
# GBR has lower open ages in 2020-2021
# URY has higher open age in 2021 only
# USA has open age of 85+, but this might be an artifact of being derived from weekly data??

small_pops <-
  out %>% 
  filter(Sex == "t") %>% 
  group_by(Code, Year) %>% 
  summarize(D = sum(Deaths), .groups = "drop") %>% 
  group_by(Code) %>% 
  summarize(D = mean(D), .groups = "drop") %>% 
  mutate(size_cat = case_when(
    between(D, 0, 1000) ~ "very small",
    between(D, 1000, 5000) ~ "small",
    TRUE ~ "OK"
  ))
  
# TWO suspects: GBR_SCO and AUS, both from who_mort_database, do they need scaling?
small_pops %>% 
  filter(size_cat == 
           "very small") %>% 
  select(Code, D) %>% 
  left_join(out) %>% 
  select(Code, Source) %>% 
  distinct()

# library(LifeIneq)
# 
# ?ineq_var
# data(LT)
# V = ineq_var(age=LT$Age,dx=LT$dx,lx=LT$lx,ex=LT$ex,ax=LT$ax)
# V[1]
# 
# mu_a <- sum((LT$Age + LT$ax) * LT$Lx) / sum(LT$Lx)
# 
# e0 <- LT$ex[1]
# 
# (2 * mu_a) / LT$ex[1] - 1
# e0 * (2 * mu_a) - e0^2 # OK
# 
# .5 * (e0 + V[1]/e0) 
# 
# 
# (e0 + (V[1]*e0)) / 2
