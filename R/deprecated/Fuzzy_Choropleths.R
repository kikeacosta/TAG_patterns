
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
LM <- tidify_GC_matrix("LogMortality","log mortality") %>% 
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
RA <- tidify_GC_matrix("RateAgingClusteredData","Rate of Aging")

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
    Centers <- 
    tidify_GC_center(fn) %>% 
      mutate(`Nr Clusters` = cc,
             `exponent` = m) %>% 
      bind_rows(Centers)
  }
}


Centers %>% 
  dplyr::filter(`Nr Clusters` == 3,
                exponent < 3) %>% 
    ggplot(aes(x = Age, y = center, color = Cluster, linetype = as.factor(exponent))) + 
    geom_line() 

# 1: red, 2: green, 3: blue




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


Clust3Maps <- list()

for (cl in c("Cluster 1", "Cluster 2", "Cluster 3")){
  for (expon in  1:3){
    Clust3Maps[[paste(expon, cl,sep = "_")]] <-
    World %>% 
      left_join(Memberships %>% 
                  dplyr::filter(`Nr Clusters` == 3, 
                                exponent == expon,
                                Cluster == cl) %>% 
                  rename("name" = "Country") ) %>% 
      ggplot() + 
      geom_sf(aes(fill = Membership), col = "white", size = 0.2)+
      scale_x_continuous(expand=c(0.03,0.03)) +
      scale_y_continuous(expand=c(0.03,0.03)) +
      scale_fill_continuous_sequential("BrwnYl", na.value = gray(.8), limits = c(0,1))
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
    
  }
}

Member3m2 <-
  Memberships %>% 
  dplyr::filter(`Nr Clusters` == 3, exponent == 2) %>% 
  rename("name" = "Country") %>% 
  left_join(World, by="name") %>% 
  drop_na(`Nr Clusters`) 

# MembershipWorld <-
#   Memberships %>% 
#   dplyr::filter(`Nr Clusters` == 3) %>% 
#   rename("name" = "Country") %>% 
#   left_join( World, by = "name")


Member3m2 %>% 
  ggplot(aes(x = Membership)) + 
  geom_density() +
  facet_grid(cols = vars(Cluster))

MainMember3m2 <- 
  Member3m2 %>% 
  group_by(name) %>% 
  mutate(maxMem = Membership == max(Membership)) %>% 
  ungroup() %>% 
  dplyr::filter(maxMem)

p1 <-
MainMember3m2 %>% 
  ggplot() +
  geom_sf(data = World,col = "white", size = 0.2, fill = gray(.8)) +
  geom_sf(data= MainMember3m2,aes(fill = factor(Cluster), geometry = geometry), col = "white", size = 0.2) + 
  scale_x_continuous(expand=c(0.03,0.03)) +
  scale_y_continuous(expand=c(0.03,0.03)) +
  guides(fill = guide_legend(title.position = "bottom",
                             keywidth = .5,
                             keyheight = .4,
                             title = "Primary cluster membership",
                             title.size = 14))+
  theme(legend.text=element_text(size = tx + 5),
        legend.key.size = unit(1, "cm"),
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

p2 <- 
Centers %>% 
  dplyr::filter(`Nr Clusters` == 3,
                exponent == 2) %>% 
  ggplot(aes(x = Age, y = center, color = Cluster)) + 
  geom_line(size = 2) +
  labs(title = "C19 Rate of aging, cluster centers",
       subtitle = "c-means, 3 clusters, membership degree = 2, ages 30-70") + 
  ylab("Rate of aging") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="italic"),
        title = element_text(size = 16))
# Cluster1
Cluster1 <- 
  Member3m2 %>% 
  ggplot() + 
  #facet_grid(rows = vars(Cluster)) +
  geom_sf(data = World,col = "white", size = 0.2, fill = gray(.8)) + 
  geom_sf(data= filter(Member3m2, Cluster == "Cluster 1"),aes(fill = Membership, geometry = geometry), col = "white", size = 0.2)+
  scale_x_continuous(expand=c(0.03,0.03)) +
  scale_y_continuous(expand=c(0.03,0.03)) +
  scale_fill_continuous_sequential("Reds 3", na.value = gray(.8), limits = c(0,1)) +
  guides(fill = guide_legend(title.position = "bottom",
                             keywidth = .5,
                             keyheight = .4))+
  theme(legend.text=element_text(size = tx + 5),
        legend.key.size = unit(1, "cm"),
        legend.title = element_blank(),
        # legend.position = c(0.1,.3),
        # legend.direction = "vertical",
        legend.position = "right",
        legend.direction = "vertical",
        plot.margin=unit(c(0,0,0,0),"cm"),
        legend.spacing = unit(c(0,0,0,0),"cm"),
        legend.margin = margin(0, 0, 0, 0),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

# Cluster 2
Cluster2 <- 
  Member3m2 %>% 
  ggplot() + 
  #facet_grid(rows = vars(Cluster)) +
  geom_sf(data = World,col = "white", size = 0.2, fill = gray(.8)) + 
  geom_sf(data= filter(Member3m2, Cluster == "Cluster 2"),aes(fill = Membership, geometry = geometry), col = "white", size = 0.2)+
  scale_x_continuous(expand=c(0.03,0.03)) +
  scale_y_continuous(expand=c(0.03,0.03)) +
  scale_fill_continuous_sequential("Greens 3", na.value = gray(.8), limits = c(0,1)) +
  guides(fill = guide_legend(title.position = "bottom",
                             keywidth = .5,
                             keyheight = .4))+
  theme(legend.text=element_text(size = tx + 5),
        legend.key.size = unit(1, "cm"),
        legend.title = element_blank(),
        # legend.position = c(0.1,.3),
        # legend.direction = "vertical",
        legend.position = "right",
        legend.direction = "vertical",
        plot.margin=unit(c(0,0,0,0),"cm"),
        legend.spacing = unit(c(0,0,0,0),"cm"),
        legend.margin = margin(0, 0, 0, 0),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

# Cluster 2
Cluster3 <- 
  Member3m2 %>% 
  ggplot() + 
  #facet_grid(rows = vars(Cluster)) +
  geom_sf(data = World,col = "white", size = 0.2, fill = gray(.8)) + 
  geom_sf(data= filter(Member3m2, Cluster == "Cluster 3"),aes(fill = Membership, geometry = geometry), col = "white", size = 0.2)+
  scale_x_continuous(expand=c(0.03,0.03)) +
  scale_y_continuous(expand=c(0.03,0.03)) +
  scale_fill_continuous_sequential("Blues 3", na.value = gray(.8), limits = c(0,1)) +
  guides(fill = guide_legend(title.position = "bottom",
                             keywidth = .5,
                             keyheight = .4))+
  theme(legend.text=element_text(size = tx + 5),
        legend.key.size = unit(1, "cm"),
        legend.title = element_blank(),
        # legend.position = c(0.1,.3),
        # legend.direction = "vertical",
        legend.position = "right",
        legend.direction = "vertical",
        plot.margin=unit(c(0,0,0,0),"cm"),
        legend.spacing = unit(c(0,0,0,0),"cm"),
        legend.margin = margin(0, 0, 0, 0),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())


(p1 / p2) | (Cluster1 / Cluster2 / Cluster3)

patch <- (p2 + plot_spacer()) 
  patch / p1

p1
p2


c3 <-
  Centers %>% 
  dplyr::filter(`Nr Clusters` == 3,
                exponent == 2)
RAp <-
RA %>% 
  ggplot(aes(x= Age, y = `Rate of Aging`, group = Country)) +
  geom_line(color = "#00000050") + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="italic"),
        title = element_text(size = 16)) + 
  labs(title = "C19 Rate of Aging for 62 countries in COVerAGE-DB",
       subtitle = "Ages 30-70, year 2020 totals") +
  ylim(-.02,.19)

p2 <- p2 + ylim(-.02,.19)


RAp | p2

p1 +
  theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="italic"),
          title = element_text(size = 16)) + 
  labs(title = "Primary cluster membership",
       subtitle = "c-means, 3 clusters, membership coef = 2\nAges 30-70, year 2020 totals")

Cluster1 / Cluster2 / Cluster3

library(cowplot)
plot_grid(Cluster1, Cluster2, Cluster3, ncol = 1, labels = c("Cluster 1","Cluster 2", "Cluster 3"))
Member3m2 %>% 
  ggplot() + 
  facet_grid(rows = vars(Cluster)) +
  geom_sf(data = World,col = "white", size = 0.2, fill = gray(.8)) + 
  geom_sf(data= Member3m2,aes(fill = Membership, geometry = geometry), col = "white", size = 0.2)+
  scale_x_continuous(expand=c(0.03,0.03)) +
  scale_y_continuous(expand=c(0.03,0.03)) +
  scale_fill_continuous_sequential("BrwnYl", na.value = gray(.8), limits = c(0,1)) +
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


  
  
  
  
library(patchwork)

cl3grid <- 
  (Clust3Maps[[1]] | Clust3Maps[[2]] | Clust3Maps[[3]])/
    (Clust3Maps[[4]] | Clust3Maps[[5]] | Clust3Maps[[6]]) / 
    (Clust3Maps[[7]] | Clust3Maps[[8]] | Clust3Maps[[9]])

cl3grid + plot_annotation("Cluster Membership by fuzziness exponent")

###########################3
# show some fits:

Fitted3c2m <- tidify_GC_matrix("FittedRateAging3cluster2m",valname = "fitted")

Compare3c2m <- 
RA %>% 
  left_join(Fitted3c2m, by = c("Country","Age")) %>% 
  pivot_longer(`Rate of Aging`:fitted, names_to = "Compare", values_to = "RA") 

Compare3c2m %>% 
  filter(Country == "USA") %>% 
    ggplot(aes(x = Age, y = RA, color = Compare)) + 
    geom_line() + 
  labs(title = "United States",
       subtitle ="c-means prediction based on weighting center patterns with membership coeficients")+
  ylab("Rate of Aging") + 
  theme(axis.text=element_text(size=12),
      axis.title=element_text(size=14,face="italic"),
      title = element_text(size = 16)) 

  


# from SVD:
# show singular vectors as lines.
# examples that 3rd component important
# scatterplot of covariates: could be done but doesn't work for these covariates
# show map of singular values to elicit suggestions of covariates.

SVmapdata <- 
read.table("Data/SingularValues.txt", 
           header = TRUE, 
           sep = " ") %>% data.matrix() %>% 
  reshape2::melt(varnames = c("name","component"), value.name = "singular value")

SVmapdata <- readr::read_delim("Data/ResponsesCovariates.txt",delim = " ") %>% 
  pivot_longer(sv1:roa.r85, names_to = "variable", values_to = "value") %>% 
  rename(name = Country)  
  
dev.new()
test <-
World %>% 
  left_join(SVmapdata, by = "name") %>% 
  dplyr::filter(variable == "sv1") %>% 
  ggplot() + 
  #geom_sf(data = World,col = "white", size = 0.2, fill = gray(.8)) + 
  geom_sf(mapping = aes(fill = value, geometry = geometry), col = "white", size = 0.2)+
  scale_x_continuous(expand=c(0.03,0.03)) +
  scale_y_continuous(expand=c(0.03,0.03)) +
  scale_fill_continuous_diverging("Purple-Green", na.value = gray(.8)) +
  guides(fill = guide_legend(title.position = "bottom",
                             keywidth = .5,
                             keyheight = .4))


SV1 <-
  World %>% 
  left_join(SVmapdata, by = "name") %>% 
  dplyr::filter(variable == "sv1") %>% 
  ggplot() + 
  geom_sf(data = World,col = "white", size = 0.2, fill = gray(.8)) + 
  geom_sf(mapping = aes(fill = -value, geometry = geometry), col = "white", size = 0.2) + 
  scale_fill_continuous_diverging("Purple-Green", na.value = gray(.8)) 

SV2 <-
  World %>% 
  left_join(SVmapdata, by = "name") %>% 
  dplyr::filter(variable == "sv2") %>% 
  ggplot() + 
  geom_sf(data = World,col = "white", size = 0.2, fill = gray(.8)) + 
  geom_sf(mapping = aes(fill = value, geometry = geometry), col = "white", size = 0.2) + 
  scale_fill_continuous_diverging("Purple-Green", na.value = gray(.8)) 


SV3 <-
World %>% 
  left_join(SVmapdata, by = "name") %>% 
  dplyr::filter(variable == "sv2") %>% 
  ggplot() + 
  geom_sf(data = World,col = "white", size = 0.2, fill = gray(.8)) + 
  geom_sf(mapping = aes(fill = value, geometry = geometry), col = "white", size = 0.2) + 
  scale_fill_continuous_diverging("Purple-Green", na.value = gray(.8), limits = c(-.25,.25)) 

plot_grid(SV1, SV2, SV3, ncol = 1, labels = c('1st','2nd','3rd'))


  
  load("Data/ETA1fit.rda")
  ETA1 <-
  ETA1array %>% 
    reshape2::melt(varnames = c("Age","Country","variable"), value.name = "value") %>% 
    mutate(Country = as.character(Country), variable = as.character(variable)) %>% 
    pivot_wider(names_from = "variable", values_from = "value")
  

  
  e0sv1 <- 
  SVmapdata %>% 
    pivot_wider(names_from = "variable", values_from = "value") %>% 
    ggplot(aes(x = e0, y = sv1)) +
    geom_text(mapping = aes(label = name)) + 
    geom_smooth()
  e0sv2 <- 
    SVmapdata %>% 
    pivot_wider(names_from = "variable", values_from = "value") %>% 
    ggplot(aes(x = e0, y = sv2)) +
    geom_text(mapping = aes(label = name)) + 
    geom_smooth()
  e0sv3 <- 
    SVmapdata %>% 
    pivot_wider(names_from = "variable", values_from = "value") %>% 
    ggplot(aes(x = e0, y = sv3)) +
    geom_text(mapping = aes(label = name)) + 
    geom_smooth()
  
  edsv1 <- 
    SVmapdata %>% 
    pivot_wider(names_from = "variable", values_from = "value") %>% 
    ggplot(aes(x = ed, y = sv1)) +
    geom_text(mapping = aes(label = name)) + 
    geom_smooth()
  edsv2 <- 
    SVmapdata %>% 
    pivot_wider(names_from = "variable", values_from = "value") %>% 
    ggplot(aes(x = ed, y = sv2)) +
    geom_text(mapping = aes(label = name)) + 
    geom_smooth()
  edsv3 <- 
    SVmapdata %>% 
    pivot_wider(names_from = "variable", values_from = "value") %>% 
    ggplot(aes(x = ed, y = sv3)) +
    geom_text(mapping = aes(label = name)) + 
    geom_smooth()
  
  
  bsv1 <- 
    SVmapdata %>% 
    pivot_wider(names_from = "variable", values_from = "value") %>% 
    ggplot(aes(x = b, y = sv1)) +
    geom_text(mapping = aes(label = name)) + 
    geom_smooth()
  bsv2 <- 
    SVmapdata %>% 
    pivot_wider(names_from = "variable", values_from = "value") %>% 
    ggplot(aes(x = b, y = sv2)) +
    geom_text(mapping = aes(label = name)) + 
    geom_smooth()
  bsv3 <- 
    SVmapdata %>% 
    pivot_wider(names_from = "variable", values_from = "value") %>% 
    ggplot(aes(x = b, y = sv3)) +
    geom_text(mapping = aes(label = name)) + 
    geom_smooth()
  
  
  (e0sv1 | e0sv2 | e0sv3) /
    (edsv1 | edsv2 | edsv3) /
    (bsv1 | bsv2 | bsv3)
  
 A <-  plot_grid(e0sv1,e0sv2,e0sv3, labels = c("1st","2nd","3rd"), nrow = 1)
 B <-  plot_grid(edsv1,edsv2,edsv3, labels = c("1st","2nd","3rd"), nrow = 1)
 C <-  plot_grid(bsv1,bsv2,bsv3, labels = c("1st","2nd","3rd"), nrow = 1)
 
 Fits <-
   ETA1 %>% 
   ggplot(aes(x = Age, y = Obs)) +
   geom_line() + 
   #geom_line(data = ETA1, mapping = aes(x = Age, y = Mod1), color = "orange")+
   #geom_line(data = ETA1, mapping = aes(x = Age, y = Mod2), color = "red")+
   geom_line(data = ETA1, mapping = aes(x = Age, y = Mod3), color = "red")+
   facet_wrap(~Country)
 
 A
 ggsave(filename = "Data/A.png")  
 B
 ggsave(filename = "Data/B.png")  
 C
 ggsave(filename = "Data/C.png")  
 
 Fits
 ggsave(filename = "Data/Fits.png") 
 