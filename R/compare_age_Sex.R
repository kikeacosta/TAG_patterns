
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
  filter(measure == "deaths") %>% 
  select(-lwr, -uppr) %>% 
  pivot_wider(names_from = `source year`, values_from = `mean`) %>% 
  select(Country, iso3, sex, age, Nx, Dxhat_WM = `Predicted 2020`, 
         Dx_expected = `Expected 2020`,  Dx_observed = `Observed 2020`,
         GHE_2019 = `GHE 2019`) %>% 
  left_join(Dxhat, by = c("iso3","sex","age")) %>% 
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

cljoin <- 
clusters %>% 
  select(iso3, Cluster)
Dx_hat_out <- 
  Dx_hat_out %>% 
  left_join(cljoin, by = "iso3")


Dx_hat_out <- 
  Dx_hat_out %>% 
  arrange(Cluster, Country, age)

pdf("Figures/observed_vs_expected_compare.pdf")
for (ctry in unique(Dx_hat_out$Country)){
  chu <-
    Dx_hat_out %>% 
    dplyr::filter(Country == ctry)%>% 
    mutate(oe_wm = Dxhat_WM / Dx_expected,
           oe_spinoff = Dxhat_spinoff / Dx_expected,
           oe_obs = Dx_observed / Dx_expected) %>% 
    pivot_longer(oe_wm:oe_obs, names_to = "variant", values_to = "OE") 
  cl <- chu$Cluster[1]
  p <- 
    chu %>% 
    ggplot(aes(x = age, y = OE, color = variant, linetype = sex)) +
    geom_line() +
    scale_y_log10() +
    labs(title = paste(ctry, ", Cluster", cl))
 print(p)
}
dev.off()

# obs_clusters <- 
# mx_plot %>% 
#   filter(variant == "mx_observed",
#          !is.na(mx)) %>% 
#   select(Country) %>% 
#   distinct() %>% 
#   mutate(obs = TRUE) %>% 
#   right_join(clusters, by = "Country") %>% 
#   filter(obs) %>% 
#   group_by(Cluster) %>% 
#   tally() %>% 
#   rename("obs" = "n")
# 
# clusters %>% 
#   group_by(Cluster) %>% 
#   tally() %>% 
#   rename( "total" = "n") %>% 
#   left_join(obs_clusters) %>% 
#   mutate(prop = obs/ total)
# 
# 
# has_data <- 
# mx_plot %>% 
#   filter(variant == "mx_observed",
#          !is.na(mx)) %>% 
#   select(Country) %>% 
#   distinct() %>% 
#   mutate(obs = TRUE) %>% 
#   right_join(clusters, by = "Country") %>% 
#   filter(obs) %>% 
#   select(iso_a3 = iso3, obs)
# 
# dput(has_data)

deltasM <- read.table("Data/deltasMobs.txt", header =TRUE)%>% 
  pivot_longer(-1, names_to = "iso3", values_to = "delta") %>% 
  mutate(sex = "m")
deltasf <- read.table("Data/deltasFobs.txt", header =TRUE) %>% 
  pivot_longer(-1, names_to = "iso3", values_to = "delta")%>% 
  mutate(sex = "f")

deltas <-
  bind_rows(deltasM, deltasf) %>% 
  rename("age" = "Age")

deltas
WPP_pop_prep  %>% 
  left_join(deltas, by = c("iso3", "sex", "age")) %>% 
  filter(iso3 == "USA")
deltas %>% 
  filter(iso3 == "USA")


Dx_hat_out %>% 
  filter(!is.na(Dx_observed)) %>% 
  mutate(mx_obs = Dx_observed / Nx) %>% 
  filter(mx_obs == min(mx_obs))
  
  ggplot(aes(x = age, y = mx_obs, color = sex, group = interaction(Country, sex))) + 
  geom_line(alpha = .3) + 
  scale_y_log10()

  
  WPP2019 <- read_csv("Output/WPP2019_baseline.csv")
  WPP2019 %>% 
    ggplot(aes(x = Age, y = mx, color = Sex, group = interaction(iso3, Sex))) +
    geom_line(alpha = .3) + 
    scale_y_log10()
  
  
  Dx_hat_out %>% 
    mutate(mx_expected = GHE_2019 / Nx) %>% 
    ggplot(aes(x = age, y = mx_expected, color = sex, group = interaction(iso3, sex)))+
    geom_line(alpha = .3) + 
    scale_y_log10()

  
  pdf("Figures/spinoff_vs_wm.pdf")
  for (ctry in unique(Dx_hat_out$Country)){
    chu  <-
      Dx_hat_out %>% 
      dplyr::filter(Country == ctry)
    cl <- chu$Cluster[1]
    p <- 
      chu %>% 
      mutate(spinoff_vs_wm = Dxhat_spinoff / Dxhat_WM) %>% 
      ggplot(aes(x = age, y = spinoff_vs_wm, color = sex)) +
      geom_line() +
      scale_y_log10() +
      labs(title = paste(ctry, "Cluster =",cl))
    print(p)
  }
  dev.off()
  

# recalculate everything using GHE baseline.
  
# 1) plot baseline, spinoff, william separately by sex. 2-panel.
# 2) plot deltas by cluster
# 3) plot observed / expected by cluster
  
GHE2019_baseline <-
  WM_estimates %>% 
    filter(`source year` == "GHE 2019",
           measure == "deaths") %>% 
    arrange(iso3, sex, age) %>% 
  select(Country, iso3, sex, age, Nx, deaths = `mean`)
  
write_csv(GHE2019_baseline, file = "Output/GHE2019_baseline.csv")
WM2020_observed <-
  WM_estimates %>% 
    filter(`source year` == "Observed 2020",
           measure == "deaths") %>% 
    arrange(iso3, sex, age) %>% 
  select(Country, iso3, sex, age, Nx, deaths = `mean`)
write_csv(WM2020_observed, file = "Output/WM2020_observed.csv")
  
  out_iso3 <- deltasM$iso3 %>% unique()
  
  out_iso3[!out_iso3 %in% wm_iso3]
  wm_iso3[!wm_iso3 %in% out_iso3]

  
dat <- read_csv("Output/annual_deaths_countries_selected_sources.csv")  
allctry <- dat %>% pull(Country) %>% unique()

COL <-
dat %>% 
  filter(Country == "Colombia",
         Year == 2020) %>% 
  mutate(Sex = case_when(Sex == "f" ~ "Female",
                         Sex == "m" ~ "Male",
                         Sex == "t" ~ "Total")) %>% 
  rename("iso3" = "Code", "sex" = "Sex", "age" = "Age" ,"deaths" = "Deaths") %>% 
  select(Country, iso3, sex, age, deaths) %>% 
  mutate(age = ifelse(age >85,85,age)) %>% 
  group_by(Country, iso3, sex, age) %>% 
  summarize(deaths = sum(deaths), .groups = "drop") 

COL <- 
WM_estimates %>% 
  filter(`source year` == "Expected 2020",
         Country == "Colombia") %>% 
  select(Country, iso3, sex, age, Nx) %>% 
  left_join(COL) 

GHE2019_baseline %>% 
  select(Country, iso3, sex, age, Nx, deaths = `mean`)
WM2020_observed %>% 
  bind_rows(COL) %>% 
  write_csv("Output/WM2020_observed.csv")




DAT <- read_csv("Output/age_sex_compare.csv")
DAT$Country %>% unique()



# ----------------------------------------------------- #
# Here, working out new estimates and graphs for GC's new deltas that
# were based on WM's inputs
# ----------------------------------------------------- #
library(readr)
library(tidyverse)
deltas   <- read.table("Data/DELTAS_WM.txt", sep = " ", header = TRUE) %>% 
  select(-row)
clusters <- read_csv("Data/iso.clusters.csv") %>% 
  select(-Country)
baseline <- read_csv("Output/GHE2019_baseline.csv")

# merge clusters and deltas to baseline.

inputs <- 
baseline %>% 
  left_join(deltas, by = c("iso3","sex","age")) %>% 
  left_join(clusters,by = "iso3")

# visualize deltas by cluster and sex
p <-
inputs %>% 
  filter(!is.na(deltas)) %>% 
  ggplot(aes(x = age, y = deltas, color = as.factor(Cluster), group = iso3)) +
  geom_line() +
  facet_wrap(~Cluster + sex) +
  theme(legend.position = "none")

pdf("Figures/deltas_clusters.pdf")
print(p)
dev.off()
inputs %>% 
  filter(!is.na(deltas),
        Cluster == 2,
        sex == "Female",
        Country != "Albania") %>% 
  ggplot(aes(x = age, y = deltas, color = Country)) +
  geom_line()


inputs %>% 
  group_by(Cluster, sex, age) %>% 
  summarize(delta = median(deltas, na.rm = TRUE)) %>% 
  ggplot(aes(x = age, y= delta, group = Cluster))+
  geom_line() +
  facet_wrap(~sex)


inputs %>% 
  group_by(Cluster, sex, age) %>% 
  summarize(delta = median(deltas, na.rm = TRUE)) %>% 
  ggplot(aes(x = age, y= delta, group = Cluster))+
  geom_line() +
  facet_wrap(~sex)
# install.packages("quantreg")
library(quantreg)
library(splines)
dat <- 
inputs %>% 
  filter(Cluster == 1, sex == "Female") 
  rqfit <- rq(deltas ~ bs(age), data = dat,tau = .5)
  
predict(rqfit, newdata = data.frame(age = seq(0,85,by=5))) %>% 
  plot()

my_smooth_median <- function(chunk){
  newdata <- data.frame(age = seq(0,85,by=5))
  
  tibble(age = seq(0,85,by=5),
         delta_qm = rq(deltas ~ bs(age),data = chunk,tau = .5) %>% 
           predict(newdata = newdata)) 
}
qm_deltas <-
inputs %>% 
  dplyr::filter(!is.na(deltas),
                ! (iso3 == "ARM" & sex == "Male"),
                ! (iso3 == "USA" & sex == "Female")) %>% 
  group_by(Cluster, sex) %>% 
  do(my_smooth_median(chunk = .data)) %>% 
  ungroup()
qm_deltas %>% 
  ggplot(aes(x = age, y = delta_qm, color = sex)) +
  geom_line() +
  facet_wrap(~Cluster)

inputs %>% 
  filter(sex == "Female", Country == "Albania")

# constraints:
library(readxl)
WM_estimates<- read_excel("Data/EstimatesBySexAge.xlsx", sheet = 2, skip = 5)

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
  select(iso3, Cluster, Yhat) 


# constrained predictions:
mx_prime <-
inputs %>% 
  left_join(qm_deltas, by = c("Cluster","sex","age")) %>% 
  group_by(Cluster, sex) %>% 
  mutate(delta_mean = mean(deltas, na.rm = TRUE),
         delta_median = median(deltas, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(mx = deaths / Nx,
         mx_pred_qm = (log(mx) + delta_qm) %>% exp(),
         mx_pred_mean = (log(mx) + delta_mean) %>% exp(),
         mx_pred_median = (log(mx) + delta_median) %>% exp()) %>% 
  select(Cluster, Country, iso3, sex, age, Nx,mx:mx_pred_median) %>% 
  left_join(WMout, by = c("Cluster","iso3")) %>% 
  mutate(Dx_pred_qm = mx_pred_qm * Nx,
         Dx_pred_mean = mx_pred_mean * Nx,
         Dx_pred_median = mx_pred_median * Nx) %>% 
  group_by(iso3) %>% 
  mutate(Dx_pred_qm = Dx_pred_qm / sum(Dx_pred_qm) * Yhat,
         Dx_pred_mean = Dx_pred_mean / sum(Dx_pred_mean) * Yhat,
         Dx_pred_median = Dx_pred_median / sum(Dx_pred_median) * Yhat,
         mx_pred_qm = Dx_pred_qm / Nx,
         mx_pred_mean = Dx_pred_mean / Nx,
         mx_pred_median = Dx_pred_median / Nx) %>% 
  select(Cluster, iso3, sex, age, mx_pred_qm, mx_pred_mean, mx_pred_median) %>% 
  pivot_longer(mx_pred_qm:mx_pred_median, values_to = "mx", names_to = "variant")

Observed <- read_csv("Output/WM2020_observed.csv")

Joined <-
Observed %>% 
  mutate(variant = "observed",
         mx = deaths / Nx) %>% 
  select(-Nx, - deaths) %>% 
  filter(nchar(iso3) == 3) %>% 
  left_join(clusters, by = "iso3") %>% 
  select(-WHO_region) %>% 
  bind_rows(mx_prime) %>% 
  arrange(Cluster, iso3, sex, variant, age) %>% 
  select(-Country)

Joined <-
inputs %>% 
  select(Country, iso3) %>% 
  distinct() %>% 
  right_join(Joined) %>% 
  arrange(Cluster, Country, sex, age) %>% 
  filter(!is.na(Cluster))
  
clusters <-
  clusters %>% 
  select(-WHO_region)
Joined <- 
WM_estimates %>% 
  filter(`source year` == "Predicted 2020",
         measure == "deaths") %>% 
  mutate(mx = `mean` / Nx) %>% 
  select(Country, iso3, sex, age, mx) %>% 
  inner_join(clusters, by = c("iso3")) %>% 
  mutate(variant = "mx_WM") %>% 
  bind_rows(Joined) %>% 
  arrange(Cluster, iso3, sex, variant, age)

iso3 <- Joined %>% pull(iso3) %>% unique()

pdf("Figures/age_sex_compare_2.pdf")
for (i in iso3){
  chunk <- 
    Joined %>% 
    filter(iso3 == i)
  ctry <- chunk$Country %>% na.omit() %>% '['(1)
  cl <-  chunk$Cluster[1]
  p <-
    chunk %>% 
    ggplot(aes(x = age, y = mx, color = variant, linetype = variant))+
    geom_line() +
    scale_y_log10() +
    facet_grid(~sex) +
    labs(title = paste("Cluster",cl,ctry))
  print(p)
}
dev.off()

write_csv(Joined, file = "Output/age_sex_compare2.csv")
  
# 1) create SR plots 
sr_data <- 
Joined %>% 
  filter(variant %in% c("mx_pred_qm","mx_WM","observed")) %>% 
  distinct() %>% 
  pivot_wider(names_from = "sex", values_from = "mx") %>% 
  mutate(sr = Male / Female)


pdf("Figures/sr_compare_2.pdf")
for (i in iso3){
  chunk <- 
    sr_data %>% 
    filter(iso3 == i)
  ctry <- chunk$Country %>% na.omit() %>% '['(1)
  cl  <-  chunk$Cluster[1]
  p <-
    chunk %>% 
    ggplot(aes(x = age, y = sr, color = variant, linetype = variant))+
    geom_line() +
    scale_y_log10(limits = c(.1,10)) +
    labs(title = paste("Cluster",cl,ctry))
  print(p)
}
dev.off()

# 2) log(predicted / observed ) for us and william (in sample)

in_sample_data <- 
  Joined %>% 
  filter(variant %in% c("mx_pred_qm","mx_WM","observed")) %>% 
  distinct() %>% 
  pivot_wider(names_from = "variant", values_from = "mx") %>% 
  mutate(spinoff = mx_pred_qm / observed,
         wm = mx_WM / observed) %>% 
  filter(!is.na(observed)) %>% 
  pivot_longer(spinoff:wm, names_to = "variant", values_to = "ratio") %>% 
  arrange(Cluster, iso3, sex, variant,age)

iso3_2 <- in_sample_data %>% pull(iso3) %>% unique()
pdf("Figures/in_sample_compare.pdf")
for (i in iso3_2){
  chunk <- 
    in_sample_data %>% 
    filter(iso3 == i)
  ctry <- chunk$Country %>% na.omit() %>% '['(1)
  cl  <-  chunk$Cluster[1]
  p <-
    chunk %>% 
    ggplot(aes(x = age, y = ratio, color = variant, linetype = variant))+
    geom_line() +
    scale_y_log10(limits = c(.25,4)) +
    facet_grid(~sex) +
    labs(title = paste("Cluster",cl,ctry)) +
    geom_hline(yintercept=1,color = "red")
  print(p)
}
dev.off()

# 3) log(us / william)  (everything)

# map
