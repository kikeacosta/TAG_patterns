rm (list = ls())
library(tidyverse)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Updating WHO files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# downloading last version of WHO files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# these come from here: https://www.who.int/data/data-collection-tools/who-mortality-database
icd_base_url <-"https://cdn.who.int/media/docs/default-source/world-health-data-platform/mortality-raw-data/"
icd_files    <- c("mort_country_codes.zip","morticd10_part1.zip","morticd10_part2.zip",
                  "morticd10_part3.zip","morticd10_part4.zip","morticd10_part5.zip")
for (i in 1:length(icd_files)){
  url_i   <- paste0(icd_base_url,icd_files[i])
  local_i <- file.path("Data", "WHO", icd_files[i])
  download.file(url_i, destfile = local_i, overwrite = TRUE)
}

# a lookup table to match country names to codes
ctry_names <- read_csv(file.path("Data", "WHO", "mort_country_codes.zip")) %>% 
  rename(Country = country)

icd_all    <- list()
icd_files2 <- icd_files[-1]
#ICD download each of the 5 files
for (i in 1:length(icd_files2)){
  icd_i <- 
    read_csv(file.path("Data", "WHO", icd_files2[i]),
             col_types = cols(Admin1 = col_character(),SubDiv = col_character(),
                              List = col_character(), Cause = col_character(), 
                              Frmat = col_character(), IM_Frmat = col_character(),
                              .default = col_double())) %>% 
    left_join(ctry_names, by = "Country") %>% 
    dplyr::filter(Sex %in% c(1,2)) 
  
  icd_all[[i]] <- icd_i
}

# stick together
# ~~~~~~~~~~~~~~
icd_all <- 
  bind_rows(icd_all) %>% 
  select(name, everything())

# saving a consolidated file with all WHO data
write_rds(icd_all, "data_inter/who_raw.rds")

