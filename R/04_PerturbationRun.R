library(MortalitySmooth)
library(magic)
library(colorspace)
library(ggplot2)
library(plotly)
library(viridis)
library(Matrix)
library(tidyverse)
source("R/04_PerturbationFunction.R")

deaths <- read.csv("data_inter/deaths_sourced_infant_based_99.csv", header=TRUE)
## loading population
offset <- read.csv("data_inter/offsets_99.csv", header=TRUE)



# source("R/temp_function.r")
# for spot testing
# deaths.j <- 
#   deaths |> 
#   filter(Country=="Andorra" & Sex == "m") 
# fit_excess(offset = offset)

# "Andorra" computationally singular

# this is a last-ditch effort since
# error-trapping in harder in parallel.
# approx 30-50 minutes execution time
do_big_fit <- FALSE
if (do_big_fit){
big_L <-
  deaths |> 
  filter(Sex != "t") |> 
  group_by(Country, Sex) |> 
  group_split()

N     <- length(big_L)
out_L <- vector("list", N)
for (n in 1:N){
  ctr  <- big_L[[n]]$Country[1]
  sx   <- big_L[[n]]$Sex[1]
  outn <- try(fit_excess(big_L[[n]], offset = offset))
  if (class(outn)[1] == "try-error"){
    outn <- tibble(ages=integer(),
                   years=integer(),
                   up=double(),
                   low=double(),
                   type=character(),
                   value=double())
  }
  outn <-outn |> 
    mutate(Country = ctr,.before=1) |> 
    mutate(Sex = sx,.after=1)
  
  out_L[[n]] <- outn
}

big_test <- 
  out_L |> 
  bind_rows() |> 
  mutate(Sex = case_when(Sex == "m" ~ "male",
                         Sex == "f" ~ "female",
                         Sex == "t" ~ "total"))

big_test |> 
  write_csv("data_inter/eta_all.csv")
}

big_test <- read_csv("data_inter/eta_all.csv")
# -------------------------------------------------- #
# Following lines are for pdf flip-books of results
yr <- 2021
ctry <- big_test |> 
  filter(years == yr) |> 
  pull(Country) |> unique()

pdf(paste0("Figures/Fitted",yr,"all.pdf"))

for (cou.j in ctry){
  
  sexes <- big_test |> 
    filter(Country == cou.j,
           years == yr) |> 
    pull(Sex) |> 
    unique()
  for (sex in sexes){
    DF <-
      big_test |> 
      filter(Country == cou.j,
             Sex == sex,
             years == yr)
    
    DFobs <- DF |> 
      filter(type == "Obs Logrates") |> 
      mutate(ageup = ages + DemoTools::age2int(ages, OAvalue = 1)) |> 
      ungroup()
    
    p <-
      ggplot(DF, aes(x = ages, y = value, color = type)) +
      geom_segment(data = DFobs,
                   aes(x = ages, y = value, xend = ageup, yend = value), linewidth = 1) +
      geom_line(data = filter(DF, type == "Fitted Logrates"), linewidth = 1) +
      geom_ribbon(data = filter(DF, type == "Fitted Logrates"),
                  aes(ymin = low, ymax = up), alpha = .2) +
      geom_line(data = filter(DF, type == "Baseline Logrates"),
                aes(y = value), linewidth = 1.2) +
      labs(x = "age", y = "log-mortality", title = paste(cou.j,sex)) +
      theme_minimal()
    print(p)
  }
}

dev.off()

yr <- 2021
ctry <- big_test |> 
  filter(years == yr) |> 
  pull(Country) |> unique()

pdf(paste0("Figures/Deltas",yr,"all.pdf"))

for (cou.j in ctry){
  
  sexes <- big_test |> 
    filter(Country == cou.j,
           years == yr) |> 
    pull(Sex) |> 
    unique()
  for (sex in sexes){
    DFa <-
      big_test |> 
      filter(Country == cou.j,
             Sex == sex,
             years == yr,
             type %in% c("Exp c","Exp Delta"))
    
    .c  <-  DFa |> 
      filter(type == "Exp c") |> 
      mutate(ages=0)
    
    DFd <- DFa |> 
      filter(type == "Exp Delta") 
    
    p <-
      ggplot(DFa, aes(x = ages, y = value, color = type)) +
      geom_hline(yintercept = 1,
                 linewidth = .5,
                 color = gray(.5)) +
      geom_line(data = DFd, linewidth = 1) +
      geom_ribbon(
        data = DFd,
        mapping = aes(ymin = low, ymax = up),
        alpha = .2,
        color = "transparent",
        fill = "red"
      ) +
      geom_point(
        data = .c,
        mapping = aes(x = ages, y = value),
        color = "red",
        size = 2
      ) +
      geom_segment(
        data = .c,
        mapping = aes(
          x = ages,
          xend = ages,
          y = low,
          yend = up
        ),
        color = "red"
      ) +
      labs(x = "age",
           y = "age perturbation",
           title = paste(cou.j, sex)) +
      theme_minimal()
    print(p)
  }
}

dev.off()

# deprecated
# library(parallel)
# 
# fit_excess_wrap <- function(deaths.j, offset){
#   out <- try(fit_excess(deaths.j = deaths.j,
#                         offset = offset))
#   if (class(out) == "try-error"){
#     return(NULL)
#   } else {
#     return(out)
#   }
# }



check_losses <- FALSE
if (check_losses){
big_test <- read_csv("data_inter/eta_all.csv")

chunks_out <-
  big_test |> 
  select(Country, Sex) |> 
  distinct()

chunks_in <-
  deaths |> 
  select(Country, Sex) |> 
  distinct() |> 
  filter(Sex!="t") |> 
  mutate(Sex = ifelse(Sex == "m","male","female"))

chunks_lost <- 
  anti_join(chunks_in,chunks_out,by = join_by(Country, Sex))

# Country                        Sex
# ------------------------------------
# Andorra                        male
# Liechtenstein                  male
# Montserrat                     male
# Saint Vincent and Grenadines   female
# Saint Vincent and Grenadines   male

}