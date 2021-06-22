library(here)
source(here("R", "00_functions.R"))

col <- read_csv("Data/colombia_celade.csv")
chl <- read_csv("Data/chile_celade.csv")

deaths <- 
  bind_rows(std_db(col), 
            std_db(chl)) %>% 
  mutate(Sex = recode(Sex,
                      "Male" = "m",
                      "Female" = "f", 
                      "Both sexes" = "t"),
         Age = ifelse(Age > 100, 100, Age)) %>% 
  group_by(Country, Year, Sex, Age) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Age = as.character(Age))

deaths2 <- 
  deaths %>% 
  filter(Age != "-2") %>% 
  bind_rows(deaths %>% 
              group_by(Country, Year, Sex) %>% 
              summarise(Deaths = sum(Deaths)) %>% 
              ungroup() %>% 
              mutate(Age = "TOT")) %>% 
  arrange(Country, Year, Sex, suppressWarnings(as.numeric(Age)))

unique(deaths2$Age) %>% sort()

deaths3 <- 
  deaths2 %>% 
  filter(Year <= 2020) %>% 
  mutate(Code = countrycode(sourcevar = Country, 
                            origin = "country.name", 
                            destination = "iso3c"),
         Source = "country_public")

col <- 
  deaths3 %>% 
  filter(Country == "Colombia")

chl <- 
  deaths3 %>% 
  filter(Country == "Chile")

write_csv(col, "Output/colombia.csv")
write_csv(chl, "Output/chile.csv")
