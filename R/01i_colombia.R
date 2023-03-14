source("R/00_functions.R")
# Colombia mortality data
# ~~~~~~~~~~~~~~~~~~~~~
# data from DANE
# https://www.dane.gov.co/index.php/estadisticas-por-tema/demografia-y-poblacion/nacimientos-y-defunciones

col_files <- unzip("Data/colombia/micro_deaths/", list = TRUE)

# all deaths from years 2016-2019
# db_mx15_20 <- tibble()
d1 <- 
  read_csv("Data/colombia/micro_deaths/nofetal2019.csv")
