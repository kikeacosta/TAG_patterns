rm (list = ls())
source("R/00_functions.R")

## R wrapper to connect to DemoData (and optionally to DemoTools library)
## https://timriffe.github.io/DDSQLtools/
## https://timriffe.github.io/DDSQLtools/articles/Downloading-UNPD-data-into-R.html
## (optional) Tools for aggregate demographic analysis
## https://timriffe.github.io/DemoTools/

## Open API documentation for DemoData:
## https://popdiv.dfs.un.org//Demodata/swagger/ui/index#

## -----------------------------------------------------------------------
## install DemoTools and DDSQLtools (comments if already installed)
## devtools::install_github("timriffe/DemoTools", force=TRUE)
## -----------------------------------------------------------------------
## devtools::install_github("timriffe/DDSQLTools", force=TRUE)
## -----------------------------------------------------------------------

# List of packages for session
.packages = c("devtools", "data.table","tictoc","dplyr")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst],dependencies=TRUE)

# Load packages into session 
lapply(.packages, require, character.only=TRUE)


# devtools::install_github("timriffe/DDSQLTools", force=TRUE)
library(DDSQLtools)

options(unpd_server = "https://popdiv.dfs.un.org/DemoData/api/")
options(scipen=999999)

## get_indicators():     Get information about available indicators (IndicatorID)
## get_iitypes():     Get information about available indicators (IndicatorID) and indicatortypeids (IndicatorTypeId)
indicators <- data.table(get_iitypes())

indicators <- indicators[, .(IndicatorID=PK_IndicatorID, IndicatorName=Name, IndicatorShortName=ShortName, UnitShortLabel, VariableType, FormatString, ComponentID, ComponentName=IndicatorType.ComponentName, IndicatorTypeID, IndicatorTypeName=IndicatorType.Name, IsComplete, SortOrder)]
setorder(indicators, SortOrder)

unique(indicators$ComponentName)
indicators[ComponentName=="Population"]
indicators[ComponentName=="Fertility"]
indicators[ComponentName=="Mortality"]
indicators[ComponentName=="Life tables"]


## list of indicators related to Deaths by age and sex
indicators[IndicatorTypeName=="Deaths by age and sex"]

## Deaths by age and sex - abridged 
## Infant and child deaths by sex and age
indicators[IndicatorID %in% c(194, 314)]



## get_indicatortypes(): Get information about available indicators (IndicatorTypeID)

## get list of DataProcess
DataProcessType <- data.table(get_dataprocesstype())
DataProcess <- data.table(get_dataprocess())

DataProcess <- merge(DataProcess[, .(DataProcessID=PK_DataProcessID, DataProcessTypeID, Name, ShortName, SortOrder1)], DataProcessType[, .(PK_DataProcessTypeID, DataProcessTypeName=Name, DataProcessTypShortNamee=ShortName)], by.x="DataProcessTypeID", by.y="PK_DataProcessTypeID")

## example of selection for DataProcessTypeID
## 2=Census ; 11=Survey ; 12=Panel ; 8=PES
## 9=Register ; 7=Life Table (legacy UN DYB) ; 10=Sample Registration System (SRS)
## 6=Estimate
DataProcess[DataProcessTypeID %in% c(9), ]



## get list of locations from Server
Locations <- data.table(get_locations(addDefault = "false",
                                      includeDependencies = "false",
                                      includeFormerCountries = "false"))
Locations <- Locations[, .(LocID=PK_LocID, LocTypeID, LocName=Name)]

DataCatalog <- data.table(get_datacatalog(addDefault = "false"))
DataProcess <- data.table(get_dataprocess(addDefault = "false"))
DataProcessType <- data.table(get_dataprocesstype(addDefault = "false"))


DataCatalog <- subset(DataCatalog, select=c("DataCatalogID", "LocID", "LocTypeID", "LocName", "DataProcessTypeID", "DataProcessType", "DataProcessTypeShortName", "DataProcessID", "DataProcess", "DataProcessShortName", "Name", "ShortName", "OfficialName", "OfficialShortName", "ReferencePeriod", "ReferenceYearStart", "ReferenceYearEnd", "ReferenceYearMid", "FieldWorkStart", "FieldWorkEnd", "FieldWorkMiddle", "ParentDataCatalogID", "isSubnational"))
DataCatalog[, FieldWorkStart := as.Date(FieldWorkStart, format="%m/%d/%Y")]
DataCatalog[, FieldWorkEnd   := as.Date(FieldWorkEnd, format="%m/%d/%Y")]
DataCatalog <- DataCatalog[is.na(LocTypeID)==FALSE]
setorder(DataCatalog, LocName, ShortName, ReferenceYearStart)


## get list of all DataSources
all_ds <- data.table(get_datasources())
dyb_sd <- data.table(get_datasources(shortNames = "DYB"))
ipums_sd <- data.table(get_datasources(shortNames = "IPUMS"))

## sample of shortnames you can use to query only 1 data source:
# myDS <- "DYB"
# myDS <- "IPUMS"
## sample of shortnames you can use to query all data source:
myDS <- NULL


## example of test using Saudi Arabia, either using LocID or name
## myLocations <- c(682)
## myLocations <- "Saudi Arabia"

## example to query all locations in our DB
myLocations <- unique(Locations$LocID)


# Here replace with number of desired chunk of countries
n_chunks <- 20
chunk_groups <- rep(1:n_chunks, length.out = length(myLocations))
cnty_groups <- split(myLocations, chunk_groups)

# Loop through each location with `lapply`
# This is useful if you want to work with many locations because 
# the API can only handle a limited volume of queries at once.
myDT <- lapply(cnty_groups, function(x) {
  # Measure time of beginning
  tic()
  
  res <- get_recorddata(dataProcessTypeIds = c(9),    ## Register
                        startYear = 2010,
                        endYear = 2021,
                        indicatorIds = c(194),  ## Deaths by age and sex - abridged 
                        # isComplete = 1,       ## 0=Abridged or 1=Complete
                        locIds = x,             ## set of locations (M49/ISO3 numerical code or M49 names)
                        locAreaTypeIds = 2,     ## "Whole area"
                        subGroupIds = 2,        ## "Total or All groups"
                        dataSourceShortNames = myDS,  ## default = NULL for all sources of information
                        includeUncertainty = FALSE,
                        collapse_id_name = FALSE)
  
  # Print time it took to make the request
  cat("Country", x, ":")
  toc()
  
  # return the result
  return(res)
})

# Merge all separate country data frames into one data frame.
DT <- data.table(do.call(rbind, myDT))



## compute number of records per series/sex
DT[, count := .N, by=list(SeriesID, SexID)]

## filter out series that are too incomplete (i.e., a typical age distribution by 5-year age group should have > 10 records up to age 50)
## DT <- DT[count >= 10]

deduplicates <- function(myDT) {	
  
  ## sort records per location and year to order multiple observation by multi-criteria using sort orders
  setorder(myDT, LocID, TimeMid, DataCatalogShortName,
           StatisticalConceptSort, 
           DataStatusSort,
           DataProcessSort, DataProcessTypeSort, 
           DataSourceSort, -DataSourceYear, DataSourceShortName,
           -DataTypeSort,
           DataReliabilitySort,
           ModelPatternName, PeriodGroupName, PeriodStart, PeriodSpan,
           SexSort, AgeStart, AgeSpan)
  
  ## subset key attributes to rank most authoritative series
  mySeries <- unique(myDT[, .(SeriesID, LocID, DataCatalogShortName, TimeMid, DataSourceShortName, DataSourceYear, DataSourceSort, DataStatusName, DataStatusSort, DataProcessSort, DataProcessTypeSort, StatisticalConceptName, StatisticalConceptSort, DataTypeName, DataTypeSort, DataReliabilityName, DataReliabilitySort)])
  
  setorder(mySeries, LocID, DataCatalogShortName,
           StatisticalConceptSort, 
           DataStatusSort,
           DataProcessSort, DataProcessTypeSort, 
           DataSourceSort, -DataSourceYear, DataSourceShortName,
           -DataTypeSort,
           DataReliabilitySort)
  
  ## assign rank to each set of "dups"
  mySeries[, nrank := 1:.N, by=list(LocID, DataCatalogShortName, trunc(TimeMid))]
  mySeries <- mySeries[nrank==1]
  
  ## keep only the most authoritative version (top #1)
  myDT <- myDT[SeriesID %in% mySeries$SeriesID]
  return(myDT)
}

DT <- deduplicates(DT)


Series <- unique(DT$SeriesID)

## warning if you query both population and mortality, you'll get many series with population only -- no mortality data from those data sources

## create/keep a list of identifiers for the SeriesID we can use to merge back to upload output datasets back into SQL database
SeriesID_Characteristics <- unique(DT[, .(SeriesID, LocID, LocName, LocTypeName, LocAreaTypeName, SubGroupName,  SubGroupTypeName, SubGroupCombinationID, 
                                          DataCatalogID, DataCatalogName, DataCatalogShortName, FieldWorkStart, FieldWorkMiddle, DataProcess, DataSourceName, DataSourceAuthor, DataSourceYear, DataSourceShortName, 
                                          DataStatusName, StatisticalConceptName, DataTypeName, ModelPatternName, DataReliabilityName, PeriodTypeName, PeriodGroupName, TimeUnit, FootNoteID)])
setorder(SeriesID_Characteristics, LocName, FieldWorkMiddle, DataSourceAuthor, DataSourceYear, DataSourceShortName, DataStatusName, StatisticalConceptName)



# ==============================================================================
# saving raw unpd data 
# ~~~~~~~~~~~~~~~~~~~~
write_rds(DT, "Data/unpd_raw.rds")
# ==============================================================================
source("R/00_functions.R")
DT <- read_rds("Data/unpd_deaths_raw.rds")

un_data <- 
  DT %>% 
  as_tibble() %>% 
  select(LocID, 
         Country = LocName, 
         Date = TimeStart, 
         Sex = SexName, 
         Age = AgeStart, 
         AgeSpan, 
         AgeLabel,
         Deaths = DataValue,
         Source = DataSourceName) %>% 
  # filter(AgeLabel != "Total") %>%
  mutate(Date = dmy(Date),
         Year = year(Date),
         Code = countrycode::countrycode(LocID, "un", "iso3c"),
         Sex = case_when(Sex == "Male" ~ "m",
                         Sex == "Female" ~ "f",
                         Sex == "Both sexes" ~ "t")) %>% 
  select(-Date, -LocID) %>% 
  group_by(Country, Source) %>% 
  filter(max(Year) >= 2020,
         !is.na(Sex)) %>% 
  mutate(Country = recode(Country,
                          "Czech Republic" = "Czechia",
                          "Iran (Islamic Republic of)" = "Iran",
                          "Bolivia (Plurinational State of)" = "Bolivia",
                          "United States of America" = "USA",
                          "Republic of Korea" = "South Korea",
                          "China, Hong Kong SAR" = "Hong Kong",
                          "Russian Federation" = "Russia",
                          "Republic of Moldova" = "Moldova")) %>% 
  ungroup()


# Closing/open age interval
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# choosing only one open age interval
closing_age <- 
  un_data %>% 
  filter(AgeSpan == -1 & AgeLabel != "Total") %>% 
  group_by(Country, Code, Source, Year, Sex) %>% 
  filter(Age == max(Age)) %>% 
  ungroup()


# Child mortality harmonization
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# grouping ages 1 to 4
child_ages <-
  un_data %>% 
  filter(AgeLabel != "Total",
         Age %in% 0:4,
         AgeSpan %in% 0:5) %>% 
  mutate(Age = ifelse(Age %in% 1:4 & AgeSpan == 1, 1, Age),
         AgeLabel = ifelse(Age %in% 1:4 & AgeSpan == 1, "1-4", AgeLabel),
         AgeSpan = ifelse(Age %in% 1:4 & AgeSpan == 1, 4, AgeSpan)) %>% 
  group_by(Country, Code, Source, Year, Sex, Age, AgeLabel, AgeSpan) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() 

# countries with three data on infant-child mortality
child_ages2 <- 
  child_ages %>% 
  group_by(Country, Code, Source, Year, Sex) %>% 
  mutate(has_inf = ifelse(any(AgeLabel == "< 1"), 1, 0),
         has_1_4 = ifelse(any(AgeLabel == "1-4"), 1, 0),
         has_0_4 = ifelse(any(AgeLabel == "0-4"), 1, 0)) %>% 
  ungroup()  
  
# weird cases
# children deaths with no infant deaths
# having deaths 1-4 but without infant mortality (neither <1 or 0-4)
no_inf <- 
  child_ages2 %>% 
  filter(has_inf == 0  & has_1_4 == 1 & has_0_4 == 0)

# extracting only infant and 1-4 mortality
dts_inf <- 
  child_ages2 %>% 
  filter(has_inf == 1 & has_1_4 == 1 & ((Age == 0 & AgeSpan == 1) | (Age == 1))) %>% 
  select(Country, Code, Source, Year, Sex, Age, Deaths, AgeSpan) 

# extracting only 0-4 mortality where there is no data on infant deaths
dts_0_4 <- 
  child_ages2 %>% 
  filter(has_inf == 0 & has_1_4 == 0 & has_0_4 == 1) %>% 
  select(Country, Code, Source, Year, Sex, Age, Deaths, AgeSpan) 
  
dts_chd <- 
  bind_rows(dts_inf,
            dts_0_4) %>% 
  arrange(Country, Source, Year, Sex, Age)

# unknown ages
unks <- 
  un_data %>% 
  filter(Age == -2) %>% 
  select(-AgeLabel)


# putting together all ages
un_data2 <- 
  un_data %>% 
  filter(AgeSpan != -1,
         Age >= 5) %>% 
  # adding closing age groups
  left_join(closing_age %>% 
              select(Country, Source, Year, Sex, max_age = Age)) %>% 
  filter(Age < max_age,
         !is.na(max_age)) %>% 
  bind_rows(closing_age) %>% 
  # adding child mortality
  bind_rows(dts_chd) %>% 
  # adding unknown ages
  bind_rows(unks) %>% 
  arrange(Country, Code, Source, Year, Sex, Age) %>% 
  # closing all ages at 100
  mutate(Age = ifelse(Age > 100, 100, Age)) %>%
  group_by(Country, Code, Source, Year, Sex, Age) %>%
  summarise(Deaths = sum(Deaths)) %>%
  ungroup()
  
# Countries with incomplete data on age
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# incomplete ages
age_incomp <- 
  un_data2 %>% 
  filter(Age >= 0) %>% 
  group_by(Country, Source, Year, Sex) %>% 
  arrange(Country, Source, Year, Sex, Age) %>% 
  mutate(AgeSpan = lead(Age) - Age,
         open_int = max(Age)) %>% 
  mutate(sum_spans = sum(AgeSpan, na.rm = TRUE)) %>% 
  filter(sum_spans != open_int) %>% 
  select(Country, Source, Year, Sex) %>% 
  ungroup() %>% 
  unique()

# countries without data on age 0
no_zero <- 
  un_data2 %>% 
  filter(Age >= 0) %>% 
  group_by(Country, Source, Year, Sex) %>% 
  filter(min(Age) != 0) %>% 
  select(Country, Source, Year, Sex) %>% 
  ungroup() %>% 
  unique()

# countries without open age interval
no_open <- 
  un_data2 %>% 
  anti_join(closing_age %>% 
              select(Country, Source, Sex, Year) %>% 
              unique()) %>% 
  select(Country, Source, Year, Sex) %>% unique()

# total age 
# ~~~~~~~~~
all_ages <-
  un_data2 %>%
  anti_join(age_incomp) %>% 
  anti_join(no_zero) %>% 
  anti_join(no_open) %>% 
  group_by(Country, Source, Year, Code, Sex) %>%
  summarise(Deaths = sum(Deaths)) %>%
  ungroup() %>%
  mutate(Age = "TOT")

# adding all together
un_data3 <- 
  un_data2 %>% 
  anti_join(age_incomp) %>% 
  anti_join(no_zero) %>% 
  anti_join(no_open) %>% 
  mutate(Age = Age %>% as.character()) %>% 
  bind_rows(all_ages) %>%
  arrange(Country, Code, Source, Year, Sex, suppressWarnings(as.integer(Age))) %>% 
  # excluding unknown ages
  filter(Age != "-2")
  

# Total sex ====
# ~~~~~~~~~~~~~~
# identifying missing total sex values
only_t_sex <- 
  un_data3 %>% 
  select(Country, Source, Year, Sex, Age) %>% 
  unique() %>% 
  mutate(id = 1) %>% 
  spread(Sex, id) %>% 
  filter(!is.na(t) & is.na(m)) %>% 
  select(Country, Source, Year) %>% 
  unique()

# identifying missing total sex values
test_sex <- 
  un_data3 %>% 
  select(Country, Source, Year, Sex, Age) %>% 
  unique() %>% 
  mutate(id = 1) %>% 
  spread(Sex, id) %>% 
  filter(is.na(t)) %>% 
  select(Country, Source, Year, Age)

# total sex for those missing it
miss_tot_sex <- 
  un_data3 %>% 
  inner_join(test_sex, by = c("Country", "Source", "Year", "Age")) %>% 
  group_by(Country, Code, Source, Age, Year) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Sex = "t") %>%
  arrange(Country, Code, Source, Year, Sex, suppressWarnings(as.integer(Age))) 

# adding total sexes when missing and re-scaling age
un_data4 <- 
  un_data3 %>% 
  bind_rows(miss_tot_sex) %>% 
  arrange(Country, Code, Source, Year, Sex, suppressWarnings(as.integer(Age))) %>% 
  group_by(Country, Source, Sex, Year) %>% 
  do(rescale_age(chunk = .data)) %>% 
  ungroup()

# re-scaling sex
un_data5 <- 
  un_data4 %>% 
  anti_join(only_t_sex) %>% 
  group_by(Country, Source, Age, Year) %>%
  do(rescale_sex(chunk = .data)) %>% 
  ungroup() %>% 
  bind_rows(un_data4 %>% 
              inner_join(only_t_sex)) %>% 
  mutate(Age = Age %>% as.double()) %>% 
  arrange(Country, Code, Source, Year, Sex, Age) 

# countries with data in 2020 or 2021
pand <- 
  un_data5 %>% 
  select(Country, Source, Year, Sex) %>% 
  unique() %>% 
  group_by(Country, Source, Sex) %>% 
  filter(max(Year) >= 2020) %>% 
  ungroup() %>% 
  select(Country, Source, Sex) %>% unique()

# minimum three periods before 2020
three <- 
  un_data5 %>% 
  select(Country, Source, Year, Sex) %>% 
  unique() %>% 
  filter(Year < 2020) %>% 
  group_by(Country, Source, Sex) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  filter(n >= 3) %>% 
  select(-n)

unique(un_data5$Source)

un_data6 <- 
  un_data5 %>% 
  semi_join(pand, by = c("Country", "Source", "Sex")) %>% 
  semi_join(three, by = c("Country", "Source", "Sex")) %>% 
  group_by(Country, Source, Year, Sex) %>% 
  mutate(age_spn = ifelse(Age == max(Age), -1, lead(Age) - Age)) %>% 
  ungroup() %>% 
  mutate(Source = case_when(Source == "Demographic Yearbook" ~ "unpd_dy",
                            Source == "Human Mortality Database" ~ "unpd_hmd",
                            TRUE ~ "unpd_crvs"))
  
# # saving for TAG analysis
# readr::write_csv(un_data6, file = "data_inter/unpd.csv")

unique(un_data6$Source)


