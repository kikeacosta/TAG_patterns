source("R/00_functions.R")


outs <- c("Armenia", "Cuba", "Dominican Republic", "Guatemala", "Mexico", "Mongolia")
outs <- c("Mongolia")

db_sca <- 
  read_delim("Output/allchat.txt", delim = '"') %>% 
  select(Country = 2, c = 3) %>% 
  mutate(c = exp(as.double(c)),
         out = ifelse(Country %in% outs, "out", "in"),
         c = ifelse(Country %in% outs, NA, c))

db_clu <- 
  read_delim("Output/Cluster.txt", delim = '"') %>% 
  select(Country = 2, Cluster = 3) %>% 
  mutate(out = ifelse(Country %in% outs, "out", "in"),
         Cluster = factor(str_trim(Cluster)))



# Map of data coverage
# ~~~~~~~~~~~~~~~~~~~~

library(cartography)
library(rgdal)
library(tmap)
library(sf)
data(World)
# checking coordinate system
st_crs(World)

World$name <- as.character(World$name)
World$name[World$name == "Swaziland"] <- "Eswatini"
World$name[World$name == "United States"] <- "USA" 
World$name[World$name == "Korea"] <- "South Korea"
World$name[World$name == "Dominican Rep."] <- "Dominican Republic"
World$name[World$name == "Czech Rep." ] <- "Czechia"
World$name[World$name == "Central African Rep." ] <- "Central African Republic"
World$name[World$name == "Eq. Guinea" ] <- "Equatorial Guinea"
World$name[World$name == "S. Sudan" ] <- "South Sudan"

# remove Antarctica
World <- World[!World$name == "Antarctica",]

# add Robinson projection
world_rob<-st_transform(World, "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
world_rob %>% ggplot() + geom_sf()

db_sca$Country[!db_sca$Country %in% world_rob$name]
# db_sca$Code[!db_sca$Code %in% world_rob$iso_a3]

map_joined <- 
  left_join(world_rob, db_sca, 
            by = c('name' = 'Country')) 

tx <- 5
library(ggnewscale)
# install.packages("ggnewscale")
# cols <- c("grey", "#52b69a", "#34a0a4", "#168aad", "#1a759f", "#1e6091", "#184e77")

map_joined %>% 
  ggplot() + 
  geom_sf(aes(fill = out), col = "black", size = 0.2) +
  scale_fill_manual(values = c("black", "transparent"), na.value = "grey80",
                    breaks = "out", labels = "")+
  new_scale("fill") +
  geom_sf(aes(fill = c), col = "black", size = 0.2) +
  # scale_x_continuous(expand=c(0.03,0.03)) +
  # scale_y_continuous(expand=c(0.03,0.03)) +
  # scale_fill_manual(values = cols) +
  # scale_fill_viridis_c()+
  scale_fill_gradient2(
    low = "#3a86ff",
    mid = "white",
    high = "#e63946",
    midpoint = 1,
    space = "Lab",
    na.value = "transparent",
    guide = "colourbar",
    aesthetics = "fill"
  )+
  # coord_sf(xlim = c(0, 4500))+
  theme(legend.text = element_text(size = tx),
        legend.title = element_text(size = tx),
        # legend.spacing = unit(c(0,0,0,0),"cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(0.5, "cm"),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank()
        # plot.background=element_blank()
        )

ggsave("Figures/scale_factor_map.png", dpi = 700, width = 6, height = 3)











map_joined <- 
  left_join(world_rob, db_clu, 
            by = c('name' = 'Country')) 

cols <- c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07")
cols <- c("#FC4E07", "#00AFBB", "#2E9FDF", "#E7B800")

tx <- 5
library(ggnewscale)
# install.packages("ggnewscale")
# cols <- c("grey", "#52b69a", "#34a0a4", "#168aad", "#1a759f", "#1e6091", "#184e77")

map_joined %>% 
  ggplot() + 
  # geom_sf(aes(fill = out), col = "black", size = 0.2) +
  # scale_fill_manual(values = c("black", "transparent"), na.value = "grey80",
  #                   breaks = "out", labels = "")+
  new_scale("fill") +
  geom_sf(aes(fill = Cluster), col = "black", size = 0.2) +
  scale_fill_manual(values = cols, breaks = c("1", "2", "3", "4"))+
  theme(legend.text = element_text(size = tx),
        legend.title = element_text(size = tx),
        # legend.spacing = unit(c(0,0,0,0),"cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(0.5, "cm"),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank()
        # plot.background = element_blank()
  )

ggsave("Figures/cluster_map.png", dpi = 700, width = 6, height = 3)

