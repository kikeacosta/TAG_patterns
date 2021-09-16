
library(tidyverse)
library(readr)
library(countrycode)
library(readxl)

WPP <- read_csv("Data/WPP2019_Life_Table_Medium.csv")

WPP
WPPout <- 
  WPP %>% 
  dplyr::filter(Sex != "Total",
                MidPeriod == 2018) %>% 
  select(LocID, Time, MidPeriod, Sex, Age = AgeGrpStart, mx) %>% 
  mutate(iso3 = countrycode::countrycode(LocID, "un", "iso3c")) %>% 
  dplyr::filter(!is.na(iso3)) 
WPPout %>% 
  write_csv("Output/WPP2019_baseline.csv")



WM_estimates<- read_excel("Data/EstimatesBySexAge.xlsx", sheet = 2, skip = 5)
clusters <- read_csv("Data/iso.clusters.csv")

WMout <- 
  WM_estimates %>% 
  dplyr::filter(measure == "deaths",
                `source year` == "Predicted 2020") %>% 
  group_by(iso3,
           sex) %>% 
  summarize(Yhat = sum(`mean`, na.rm = TRUE),
            .groups = "drop") %>% 
  ungroup() %>% 
  pivot_wider(names_from = sex, values_from = Yhat) %>% 
  filter(nchar(iso3) == 3) %>% 
  mutate(Yhat = Female + Male) %>% 
  left_join(clusters, by = "iso3") %>% 
  select(Country, iso3, WHO_region, Cluster, Yhat, Female, Male) 

WMout %>% 
  write.csv("Output/WM_constraints_and_clusters.csv")


read_csv("Output/offsets.csv")  %>% 
  mutate(AgeWM = Age - Age %% 5,
         AgeWM = ifelse(Age < 100, AgeWM, 100),
         AgeWM = ifelse(Age != 0 & AgeWM == 0, 1, AgeWM)) %>% 
  group_by(Country, Code, Year, Sex, AgeWM) %>% 
  summarize(Population = sum(Population)) %>% 
  rename("Age" = "AgeWM") %>% 
  write_csv("Output/offsets_WM.csv")


# ------------------------------------ #
dM <- read_delim("Data/deltasM.txt",delim = " ", col_names = c("Age",seq(1,8)), skip = 1) %>% 
  pivot_longer(-1, names_to = "Cluster", values_to = "delta") %>% 
  mutate(sex = "Male")
dF <- read_delim("Data/deltasF.txt",delim = " ", col_names = c("Age",seq(1,8)), skip = 1) %>% 
  pivot_longer(-1, names_to = "Cluster", values_to = "delta") %>% 
  mutate(sex = "Female")

Deltas <- 
  dM %>% 
  bind_rows(dF)
Deltas

WMout

WPPprime <- 
  WPPout %>% 
  filter(Age < 100) %>% 
  rename("sex" = "Sex") %>% 
  left_join(clusters, by = "iso3") %>% 
  mutate(Cluster = as.character(Cluster)) %>% 
  left_join(Deltas, by = c("Cluster", "sex","Age")) %>% 
  mutate(mx_prime = log(mx) + delta,
         mx_prime = exp(mx_prime)) %>% 
  select(-WHO_region, -LocID, -Country)

# for sex-specific constraints
# WPPprimeprime <-
#   WMout %>% 
#   select(-Yhat, -Cluster) %>% 
#   pivot_longer(Female:Male, names_to = "sex", values_to = "TOT") %>% 
#   right_join(WPPprime, by = c("iso3","sex"))

# for overall constraint
WPPprimeprime <-
  WMout %>% 
  select(iso3, TOT = Yhat) %>% 
  right_join(WPPprime, by = c("iso3"))

WPPpop <- read_csv("Data/WPP2019_PopulationBySingleAgeSex_1950-2019.csv")

WPP_pop_prep <-
  WPPpop %>% 
  filter(Time == max(Time)) %>% 
  select(LocID, Age = AgeGrp, Male = PopMale, Female = PopFemale) %>% 
  pivot_longer(Male:Female, names_to = "sex", values_to = "Population") %>% 
  mutate(iso3 =  countrycode::countrycode(LocID, "un", "iso3c")) %>% 
  filter(!is.na(iso3)) %>% 
  select(-LocID) %>% 
  mutate(AgeWM = Age - Age %% 5,
         AgeWM = ifelse(Age < 100, AgeWM, 100),
         AgeWM = ifelse(Age != 0 & AgeWM == 0, 1, AgeWM)) %>% 
  select(-Age) %>% 
  group_by(iso3, sex, Age = AgeWM) %>% 
  summarize(Population = sum(Population) * 1000, .groups = "drop") 

# our "final" estimate
Dxhat <-
  WPPprimeprime %>% 
  left_join(WPP_pop_prep, by = c("iso3", "sex","Age")) %>% 
  mutate(Dxhat = mx_prime * Population) %>% 
  group_by(iso3) %>% 
  mutate(Dxhat = Dxhat / sum(Dxhat) * TOT) %>% 
  ungroup() %>% 
  select(-TOT, -Time, -MidPeriod) %>% 
  mutate(age = ifelse(Age == 1, 0, Age),
         age = ifelse(age > 85, 85, age)) %>% 
  filter(!is.na(Dxhat)) %>% 
  group_by(iso3, sex, age) %>% 
  summarize(Dxhat_spinoff = sum(Dxhat), .groups = "drop")

obs <- 
  WM_estimates %>% 
  filter(measure == "deaths",
         `source year` == "Observed 2020") %>% 
  filter(!is.na(mean)) %>% 
  select(iso3, sex, age, Dx_observed = mean)

Dx_hat_out <-
  WM_estimates %>% 
  filter(measure == "deaths",
         `source year` == "Predicted 2020") %>% 
  select(Country, iso3, sex, age, Nx, Dxhat_WM = mean) %>% 
  left_join(Dxhat, by = c("iso3","sex","age")) %>% 
  left_join(obs) %>% 
  filter(!is.na(Dxhat_WM) & !is.na(Dxhat_spinoff))

write_csv(Dx_hat_out, file = "Output/age_sex_compare.csv")

mx_plot <-
  Dx_hat_out %>% 
  mutate(mx_wm = Dxhat_WM / Nx,
         mx_spinoff = Dxhat_spinoff / Nx,
         mx_observed = Dx_observed / Nx) %>% 
  select(Country, sex, age, mx_wm, mx_spinoff, mx_observed) %>% 
  pivot_longer(mx_wm:mx_observed, names_to = "variant", values_to = "mx") 

pdf("Figures/age_sex_compare.pdf")
for (i in unique(mx_plot$Country)){
  p <- 
    mx_plot %>% 
    filter(Country == i) %>% 
    ggplot(aes(x = age, y = mx, color = sex, linetype = variant)) +
    geom_line() +
    scale_y_log10() +
    labs(title = i)
  print(p)
}
dev.off()
pdf("Figures/sr_compare.pdf")
for (i in unique(mx_plot$Country)){
p <- 
  mx_plot %>% 
  filter(Country == i) %>% 
  pivot_wider(names_from = sex, values_from = mx) %>% 
  mutate(sr = Male / Female) %>% 
  ggplot(aes(x = age, y = sr, color = variant, linetype = variant)) +
  geom_line() +
  scale_y_log10() +
  labs(title = i)
print(p)
}
dev.off()
