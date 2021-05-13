
source("https://raw.githubusercontent.com/timriffe/covid_age/master/R/00_Functions.R")

library(googlesheets4)
library(tidyverse)
library(cartography)
library(rgdal)
library(tmap)
library(sf)
library(readr)
library(colorspace)
tidify_GC_matrix <- function(filename, valname){
  data.matrix(read.csv(paste0("Data/OutCluster/",filename,".txt"), 
                       header = TRUE, row.names = 1, sep = " ")) %>% 
    reshape2::melt(value.name = valname, varnames = c("Age - 29", "Country")) %>% 
    mutate(Age = `Age - 29` + 29) %>% 
    select(Country, Age, any_of(valname))
}
tidify_GC_center <- function(filename){
  read.csv(filename, header = TRUE, sep = " ") %>% 
    data.matrix() %>% 
    reshape2::melt(varnames = c("Cluster","x"), value.name = "center") %>% 
    mutate(Age = gsub(x, pattern = "p.", replacement = "") %>% as.integer() %>% '+'(29),
           Cluster = as.character(Cluster)) %>% 
    select(-x)
}
tidify_GC_memberships <- function(filename){
  read.csv(filename, header = TRUE, sep = " ") %>% 
    data.matrix() %>% 
    reshape2::melt(varnames = c("Country","Cluster"), value.name = "Membership") %>% 
    mutate(Country = as.character(Country),
           Cluster = as.character(Cluster),
           Cluster = gsub(Cluster, pattern = "\\.", replacement = " ")) %>% 
    as_tibble()
}



data(World)

# checking coordinate system
st_crs(World)

World$name <- as.character(World$name)

World <- 
  World %>% 
  mutate(name = case_when(name == "Czech Rep." ~ "Czechia",
                          name == "Swaziland" ~ "Eswatini",
                          name == "United States" ~ "USA",
                          TRUE ~ name))
# remove Antarctica
World <- World[World$name != "Antarctica",] %>% 
  st_transform("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

# Log mortality
LM <- tidify_GC_matrix("RateAgingClusteredData","log mortality") %>% 
  mutate(`C19 Mortality Rate` = exp(`log mortality`))

LM %>% 
  ggplot(aes(x= Age, y = `C19 Mortality Rate`, group = Country)) +
  geom_line(color = "#00000050") + 
  scale_y_log10() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="italic"),
        title = element_text(size = 16)) + 
  labs(title = "C19 mortality rates for 62 countries in COVerAGE-DB",
       subtitle = "Ages 30-70, year 2020 totals")

# Rate of aging values, following PCLM
RA <- tidify_GC_matrix("/RateAgingClusteredData","Rate of Aging")

RA %>% 
  ggplot(aes(x= Age, y = `Rate of Aging`, group = Country)) +
  geom_line(color = "#00000050") + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="italic"),
        title = element_text(size = 16)) + 
  labs(title = "C19 Rate of Aging for 62 countries in COVerAGE-DB",
       subtitle = "Ages 30-70, year 2020 totals")

tidify_GC_matrix("/RateAgingClusteredData","Rate of Aging")

all_files <- dir("Data/OutCluster")
fitted_files <- grep(all_files,pattern = "Fitted",value=TRUE)
centers_files <- grep(all_files,pattern = "Center",value=TRUE)



Centers <- tibble( Cluster = character(),     
                   center = double(), 
                   Age = double(), 
                   `Nr Clusters` = integer(), 
                   exponent = integer())

for (cc in 3:5){
  for (m in 1:3){
    fn <- paste0("Data/OutCluster/Centers",cc,"cluster",m,"m.txt")
    Memberships <- 
    tidify_GC_center(fn) %>% 
      mutate(`Nr Clusters` = cc,
             `exponent` = m) %>% 
      bind_rows(Memberships)
  }
}


Centers %>% 
    ggplot(aes(x = Age, y = center, color = Cluster, linetype = as.factor(exponent))) + 
    geom_line() + 
    facet_wrap(~`Nr Clusters`)
  

Memberships <- tibble(Country = character(),
                      Cluster = character(),
                      Membership = double())

for (cc in 3:5){
  for (m in 1:3){
    fn <- paste0("Data/OutCluster/Membership",cc,"cluster",m,"m.txt")
    Memberships <- 
      tidify_GC_memberships(fn) %>% 
      mutate(`Nr Clusters` = cc,
             `exponent` = m) %>% 
      bind_rows(Memberships)
  }
}

# Fuzzy Maps
# -----------------------------
tx=7

# Not working as expected just yet. grr  
World %>% 
  left_join(Memberships %>% 
              dplyr::filter(`Nr Clusters` == 3, 
                            exponent == 3,
                            Cluster == "Cluster 3") %>% 
              rename("name" = "Country") ) %>% 
  ggplot() + 
  geom_sf(aes(fill = Membership), col = "white", size = 0.2)+
  scale_x_continuous(expand=c(0.03,0.03)) +
  scale_y_continuous(expand=c(0.03,0.03)) +
  scale_fill_continuous_sequential("BrwnYl", na.value = gray(.8))
  guides(fill = guide_legend(title.position = "bottom",
                             keywidth = .5,
                             keyheight = .4))+
  theme(legend.text=element_text(size = tx + 5),
        legend.key.size = unit(1, "cm"),
        legend.title = element_blank(),
        # legend.position = c(0.1,.3),
        # legend.direction = "vertical",
        legend.position = "bottom",
        legend.direction = "horizontal",
        plot.margin=unit(c(0,0,0,0),"cm"),
        legend.spacing = unit(c(0,0,0,0),"cm"),
        legend.margin = margin(0, 0, 0, 0),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        #panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

