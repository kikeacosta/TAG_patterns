

library(tidyverse)
library(readr)

A <- read_delim("Data/CLUSTERS2.txt",delim=" ")
B <- A[,-1] %>% as.matrix()
rownames(B) <-unlist(A[,1])
plot(cor(B))


D <- matrix(NA,43,43, dimnames = list(to = rownames(B), from = rownames(B)))

for (i in 1:43){
  for (j in 1:43){
    D[i,j] = sum(B[i,] == B[j, ])
  }
}
D = D / max(D)
image(D)

diag(D) <- NA

D[D < 8] <- NA
ctries <- colnames(D)
ctries[d$rowInd]
gp1 <- c("Brazil","Colombia","Belgium")
gp2 <- c("France","Czechia","Norway","Portugal","Italy","Canada","Netherlands","Israel","Scotland","Switzerland",
         "Spain","Slovenia","Australia","Germany")
gp3 <- c("Chad","Malawi","Pakistan","Nigeria","Ukraine","Romania","Iraq","Eswatini","India","Mexico")
gp4 <- c("Japan","Greece","Uruguay","Jordan","Cuba","Hungary","Argentina","USA","Philippines","Nepal","Peru","Kenya","Turkey","Chile","Panama","Paraguay")

length(gp1) + length(gp2) + length(gp3) + length(gp4)

MainClusters <- tibble(name = ctries) %>% 
  mutate(group = case_when(name %in% gp1 ~ "A",
                           name %in% gp2 ~ "B",
                           name %in% gp3 ~ "C",
                           name %in% gp4 ~ "D"))

source("https://raw.githubusercontent.com/timriffe/covid_age/master/R/00_Functions.R")

library(googlesheets4)
library(tidyverse)
library(cartography)
library(rgdal)
library(tmap)
library(sf)
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
World <- World[World$name != "Antarctica",]
World <- World %>% 
  left_join(MainClusters) %>% 
  left_join(A)

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
  geom_sf(aes(fill = as.factor(Cadult4)), col = "white", size = 0.2) +
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




