

library(tidyverse)
library(readr)
# library(googlesheets4)
library(cartography)
library(rgdal)
library(tmap)
library(sf)

clusters <- read_csv("Data/iso.clusters.csv") %>% 
  select(iso_a3 = iso3, Cluster) %>% 
  left_join(

structure(list(iso_a3 = c("ALB", "ARM", "AUS", "AUT", "BEL", 
                          "BGR", "BRA", "CAN", "CHE", "CHL", "CRI", "CYP", "CZE", "DEU", 
                          "DNK", "ECU", "ESP", "EST", "FIN", "FRA", "GBR", "GEO", "GRC", 
                          "HRV", "HUN", "IRL", "IRQ", "ISL", "ISR", "ITA", "JPN", "KOR", 
                          "LTU", "LUX", "LVA", "MDA", "MEX", "MLT", "MNE", "MUS", "NLD", 
                          "NOR", "NZL", "PER", "POL", "PRT", "PRY", "ROU", "RUS", "SRB", 
                          "SVK", "SVN", "SWE", "TUN", "UKR", "URY", "USA", "ZAF"), obs = c(TRUE, 
                                                                                           TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
                                                                                           TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
                                                                                           TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
                                                                                           TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
                                                                                           TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
                                                                                           TRUE, TRUE)), row.names = c(NA, -58L), class = c("tbl_df", "tbl", 
                                                                                                                                            "data.frame"))) %>% 
  mutate(alpha = ifelse(is.na(obs),.4,1))


data(World)

# checking coordinate system
st_crs(World)


World$iso_a3 <- as.character(World$iso_a3)


# lose Liechtenstein, Malta, Mauritius
# los E&W, N Ireland, Scotland

# remove Antarctica
World <- World[World$name != "Antarctica",]
World <- World %>% 
  left_join(clusters,by="iso_a3") 

# add Robinson projection
# world_rob<-
# world_rob %>% ggplot() + geom_sf()

tx        <- 7
# "C3" "C4"      "C5"      "Cadult3" "Cadult4" "Cadult5"
# "S3"      "S4"      "S5"      "Sadult3" "Sadult4" "Sadult5"
# "A3"      "A4"      "A5"
map_out <- 
World %>% 
  st_transform("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") %>% 
  ggplot() + 
  geom_sf(aes(fill = as.factor(Cluster), alpha = `alpha`), col = "white", size = 0.2) +
  scale_x_continuous(expand=c(0.03,0.03)) +
  scale_y_continuous(expand=c(0.03,0.03)) +
  scale_alpha_identity() +
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

ggsave("Figures/wm_clusters.pdf", map_out)




