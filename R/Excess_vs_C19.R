
library(covidAgeData) # for COVerAGE-DB
library(tidyverse)
library(lubridate)
library(osfr)
library(DemoTools)

Excess <- read_csv("https://raw.githubusercontent.com/christina-bohk-ewald/assess-total-infection-fatality-burden-of-COVID-19/main/input-data/cumulative_excess_age_2020_2021.csv")

# Pick date around end of 2020
Excess <-
  Excess %>% 
  group_by(Country, Sex) %>% 
  mutate(DateDiff = abs(Date - ymd("2020-12-31"))) %>% 
  group_by(Country, Sex) %>% 
  dplyr::filter(DateDiff == min(DateDiff),
                Date >= ymd("2020-08-01")) 


download_covid("Output_5", dest = "Data", download_only = TRUE)



OWD <- read_csv("https://covid.ourworldindata.org/data/owid-covid-data.csv") 


C19deaths <- read_subset_covid(zippath = "Data/Output_5.zip",
                               data = "Output_5",
                               Region = "All")

C19_total_deaths <-
  OWD %>% 
  dplyr::filter(date == as_date("2021-01-01"),
                !is.na(continent)) %>% 
  select(total_deaths, 
         Country = location,
         Continent = continent) %>% 
  mutate(Country = ifelse(Country ==  "United States" ,"USA",Country))

unique(C19deaths$Country) 
unique(C19_total_deaths$Country)

C19deaths_final <- 
  C19deaths %>% 
  dplyr::filter(!is.na(Deaths)) %>% 
  select(Country, Sex, Date, Age, AgeInt, C19_deaths = Deaths) %>% 
  mutate(Date = dmy(Date),
         DateDiff = abs(Date - ymd("2020-12-31"))) %>% 
  group_by(Country, Sex) %>% 
  dplyr::filter(DateDiff == min(DateDiff),
                Date >= ymd("2020-08-01")) %>% 
  ungroup() %>% 
  left_join(C19_total_deaths, by = "Country") %>% 
  pivot_wider(names_from = Sex, values_from = C19_deaths) %>% 
  group_by(Country) %>% 
  mutate(
    b = ifelse(is.na(b), f + m, b),
    total_deaths = ifelse(is.na(total_deaths),sum(b),total_deaths),
    b = b / sum(b, na.rm = TRUE) * total_deaths) %>% 
  ungroup() %>% 
  mutate(pm = m / (m + f),
         pm = ifelse(is.nan(pm), .5, pm),
         f = b * (1 - pm),
         m = b * pm
         ) %>% 
  pivot_longer(f:b, names_to = "Sex", values_to = "C19_deaths") %>% 
  select(Country, Continent, Sex, Age, C19_deaths) %>% 
  group_by(Country, Sex) %>% 
  mutate(C19_deaths = ifelse(Age == 90, 
                             sum(C19_deaths[Age >= 90]), 
                             C19_deaths)) %>% 
  ungroup() %>% 
  dplyr::filter(Age <= 90)


compare_countries <- intersect(unique(Excess$Country), unique(C19deaths$Country))

Compare <- inner_join(Excess, C19deaths_final, by = c("Country","Sex","Age")) %>% 
  arrange(Country, Sex, Age) %>% 
  select(Country, Sex, Age, CumEpi, CumExc, CumPos, C19_deaths, Exposure)


# CumPos vs C19 rates:

Compare %>% 
  mutate(Excess = CumPos / Exposure,
         C19 = C19_deaths / Exposure) %>% 
  select(Country, Sex, Age, Excess, C19) %>% 
  pivot_longer(Excess:C19, names_to = "Measure", values_to = "Rate") %>% 
  ggplot(aes(x = Age, 
             y = Rate, 
             group = interaction(Country, Measure), color = Measure)) + 
    geom_line(alpha = .5, size = 1) +
    facet_wrap(~Sex) + 
  scale_y_log10() +
  labs(title = "Positive Excess vs C19 rates",
       caption = "Excess by Enrique Acosta based on STMF
       Positive deviations only
       C19 data from COVerAGE-DB, dates close to Dec 31, 2020
       and scaled to Hopkins total") +
  theme(plot.title = element_text(size=35))

# CumExc vs C19 rates:
Compare %>% 
  mutate(Excess = CumExc / Exposure,
         C19 = C19_deaths / Exposure) %>% 
  select(Country, Sex, Age, Excess, C19) %>% 
  pivot_longer(Excess:C19, names_to = "Measure", values_to = "Rate") %>% 
  ggplot(aes(x = Age, 
             y = Rate, 
             group = interaction(Country, Measure), color = Measure)) + 
  geom_line(alpha = .5, size = 1) +
  facet_wrap(~Sex) + 
  scale_y_log10()  +
  labs(title = "Total Excess vs C19 rates",
       caption = "Excess by Enrique Acosta based on STMF
       includes negative and positive deviations
       C19 data from COVerAGE-DB, dates close to Dec 31, 2020
       and scaled to Hopkins total") +
  theme(plot.title = element_text(size=35))
  


# CumEpi vs C19 rates:
Compare %>% 
  mutate(Excess = CumEpi / Exposure,
         C19 = C19_deaths / Exposure) %>% 
  select(Country, Sex, Age, Excess, C19) %>% 
  mutate(Excess = ifelse(Excess == 0, NA_real_, Excess)) %>% 
  pivot_longer(Excess:C19, names_to = "Measure", values_to = "Rate") %>% 
  ggplot(aes(x = Age, 
             y = Rate, 
             group = interaction(Country, Measure), color = Measure)) + 
  geom_line(alpha = .5, size = 1) +
  facet_wrap(~Sex) + 
  scale_y_log10() +
  labs(title = "Epi Excess vs C19 rates",
       caption = "Excess by Enrique Acosta based on STMF
       includes only positive deviations beyond 95% upper,
       C19 data from COVerAGE-DB, dates close to Dec 31, 2020
       and scaled to Hopkins total")+
  theme(plot.title = element_text(size=35))


# CumExc vs C19 rates:
Compare %>% 
  dplyr::filter(Country == "Peru") %>% 
  mutate(Excess = CumExc / Exposure,
         C19 = C19_deaths / Exposure) %>% 
  select(Country, Sex, Age, Excess, C19) %>% 
  pivot_longer(Excess:C19, names_to = "Measure", values_to = "Rate") %>% 
  ggplot(aes(x = Age, 
             y = Rate, 
             group = interaction(Country, Measure), color = Measure)) + 
  geom_line(alpha = .5, size = 1) +
  facet_wrap(~Sex) + 
  scale_y_log10()  +
  labs(title = "Total Excess vs C19 rates",
       caption = "Excess by Enrique Acosta based on STMF
       includes negative and positive deviations
       C19 data from COVerAGE-DB, dates close to Dec 31, 2020
       and scaled to Hopkins total") +
  theme(plot.title = element_text(size=35))


Compare %>% 
  mutate(Excess = CumExc / Exposure,
         C19 = C19_deaths / Exposure) %>% 
  select(Country, Sex, Age, Excess, C19) %>% 
  pivot_longer(Excess:C19, names_to = "Measure", values_to = "Rate") %>% 
  ggplot(aes(x = Age, 
             y = Rate, 
             group = interaction(Sex, Measure), 
             color = Sex,
             linetype = Measure)) + 
  geom_line(alpha = .5, size = 1) +
  facet_wrap(~Country) + 
  scale_y_log10()  +
  labs(title = "Total Excess vs C19 rates",
       caption = "Excess by Enrique Acosta based on STMF
       includes negative and positive deviations
       C19 data from COVerAGE-DB, dates close to Dec 31, 2020
       and scaled to Hopkins total") +
  theme(plot.title = element_text(size=35))
# --------------------------------------- #
# Repeat with confidence intervals!       #
# --------------------------------------- #

ExcessCI <- readRDS("Data/excess_tim_enrique/baseline_mortality.rds")


C19_inputs <- readRDS("Data/C19_use.rds")


offsets <- readRDS("Data/offsets.rds")

countries <- C19_inputs$Country %>% unique()

library(ungroup)

mods <- list()

for (i in countries){
  c19i   <- filter(C19_inputs, Country == i, Age < 105) 
  y      <- c19i$Deaths
  x      <- c19i$Age
  nlast  <- rev(c19i$AgeInt)[1]
  off    <- filter(offsets, Country == i) %>% dplyr::pull(Population)
  off    <- off[1:105]
  ai     <- c19i$AgeInt
  y1     <- rep(y/ai,ai) 
  for (z in 0:3){
    Z <- 10^z
    mod <- pclm(x = x, 
                y = y * Z,
                nlast = nlast,
                offset = off, 
                control = list(deg = 3, lambda = 1e4))
    
    plot(0:104, y1 / off, log = 'y', main = Z)
    lines(mod$fitted / Z)
   # locator(1)
    
    if (all(mod$ci$lower > 0)){
      break
    }

     
  }
  mod[["c19_total"]] <- c19i$Deaths %>% sum()
  mod[["Z"]] <- Z
  mod[["Country"]] <- i
  mods[[i]] <- mod
  
}

# rescale loop:

C19_out <- lapply(mods, function(X){
  
  # back out counts and rescale to TOT
  TOT <- X$c19_total
  off <- X$input$offset
  d1  <- X$fitted * off
  MAR <- sum(d1)
  FAC <- TOT / MAR
  d1  <- d1 * FAC # seems OK
  
  # ------------------------------------ #
  # This doesn't work well everywhere..  #
  l   <- X$ci$lower * off * FAC          #
  u   <- X$ci$upper * off * FAC          #
  # ------------------------------------ #
  
  # Group to 5 year bins
  a   <- X$bin.definition$output$breaks["left",]
  d5  <- groupAges(d1, a, N = 5, OAnew = 90)
  l5  <- groupAges(l, a, N = 5, OAnew = 90)
  u5  <- groupAges(u, a, N = 5, OAnew = 90)
  p5  <- groupAges(off, a, N = 5, OAnew = 90)
  
  
  tibble(Country = X$Country[1],
         Age = DemoTools::names2age(d5),
         C19_deaths = d5,
         lower = l5,
         upper = u5,
         pop = p5)
}) %>% bind_rows()



C19_out %>% 
  group_by(Country) %>% 
  dplyr::filter(sum(C19_deaths) > 200) %>% 
  ggplot(aes(x = Age, 
             y = C19_deaths / pop, 
             ymin = lower/ pop, 
             ymax = upper/ pop, 
             group = Country)) +
  geom_line() + 
  geom_ribbon(alpha = .5) + 
  facet_wrap(~Country) + 
  scale_y_log10(limits = c(1e-6,1)) +
  xlim(30,90)


# C19_out %>% 
#   dplyr::filter(Age == 50) %>% 
#   arrange(-C19_deaths)

offsets %>% 
  dplyr::filter(Country == "India") %>% 
  ggplot(aes(x=Age,y=Population)) +
  geom_line()

Codes <- read_csv("https://raw.githubusercontent.com/kikeacosta/excess_mortality/master/Data/country_codes.csv")

ExcessCI_out <-
  ExcessCI %>% 
  dplyr::filter(Year == 2020) %>% 
  mutate(Excess = Deaths - Baseline,
         lower = Deaths - up,
         upper = Deaths - lp) %>% 
  group_by(PopCode,Age) %>% 
  summarize(point = sum(Excess[Excess > 0], na.rm = TRUE),
            lower = sum(lower[lower > 0], na.rm = TRUE),
            upper = sum(upper[upper > 0], na.rm = TRUE),
            .groups = "drop") %>% 
  left_join(Codes, by = "PopCode") %>% 
  pivot_longer(point:upper, names_to = "quantile", values_to = "Excess") %>% 
  select(-PopCode, -wpp_name)
 
ExcessCI_out %>% 
  pivot_wider(names_from = "quantile", values_from = "Excess") %>% 
  ggplot(aes(x = Age, y = point, ymin = lower, ymax = upper)) + 
  geom_line() + 
  geom_ribbon(alpha = .5) +
  facet_wrap(~Country)


excess_countries <- ExcessCI_out$Country %>% unique()

CompareCI <-
  C19_out %>% 
  rename(point = C19_deaths) %>% 
  pivot_longer(point:upper, names_to = "quantile", values_to = "C19 Deaths") %>% 
  inner_join(ExcessCI_out, by = c("Country","Age","quantile")) %>% 
  pivot_longer(`C19 Deaths`:Excess, names_to = "Metric", values_to = "Count") %>% 
  pivot_wider(names_from = "quantile", values_from = "Count")


CompareCI %>% 
  mutate(
    point = ifelse(point < 0,NA_real_, point),
    lower = ifelse(lower < 0, point / 2, lower)
         ) %>% 
  ggplot(aes(x = Age, 
             y = point / pop, 
             ymin = lower / pop, 
             ymax = upper / pop, 
             group = Metric, 
             color = Metric,
             fill = Metric)) + 
  geom_line() + 
  geom_ribbon(alpha = .2) +
  facet_wrap(~Country) +
  scale_y_log10(limits = c(1e-6,2)) +
  xlim(30,90)

ExcessCI %>% 
  ggplot(aes(x = Date, y = Deaths)) + 
  geom_line() +
  
