

library(tidyverse)
library(readr)
# library(googlesheets4)
library(cartography)
library(rgdal)
library(tmap)
library(sf)

map_these <- read_csv("Output/annual_deaths_countries_selected_sources.csv") %>% 
  distinct(Country, Source) %>% 
  mutate(Source= case_when(Source == "stmf modified" ~ "stmf",
                           Source %in% c("brazil","mexico","peru") ~ "nso",
                           TRUE ~Source)) %>% 
  select(name = Country, Source)

map_these$Source %>% unique()

source("https://raw.githubusercontent.com/timriffe/covid_age/master/R/00_Functions.R")


data(World)

# checking coordinate system
st_crs(World)

World$name <- as.character(World$name)

map_these$name
World <- 
  World %>% 
  mutate(name = case_when(name == "Czech Rep." ~ "Czechia",
                          name == "United States" ~ "USA",
                          name == "Korea" ~ "South Korea",
                          TRUE ~ name))
World$name
map_these$name[!map_these$name %in%World$name]
# lose Liechtenstein, Malta, Mauritius
# los E&W, N Ireland, Scotland

# remove Antarctica
World <- World[World$name != "Antarctica",]
World <- World %>% 
  left_join(map_these,by="name") 

# add Robinson projection
# world_rob<-
# world_rob %>% ggplot() + geom_sf()

tx        <- 7
# "C3" "C4"      "C5"      "Cadult3" "Cadult4" "Cadult5"
# "S3"      "S4"      "S5"      "Sadult3" "Sadult4" "Sadult5"
# "A3"      "A4"      "A5"

World %>% 
  st_transform("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") %>% 
  ggplot() + 
  geom_sf(aes(fill = as.factor(Source)), col = "white", size = 0.2) +
  scale_x_continuous(expand=c(0.03,0.03)) +
  scale_y_continuous(expand=c(0.03,0.03)) +
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

# ggsave("/home/tim/Desktop/stmf_map_pre.svg", map_out)




