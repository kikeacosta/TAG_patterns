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
write_rds(DT, "Data/unpd_deaths_raw.rds")
# ==============================================================================
rm (list = ls())
source("R/00_functions.R")
DT <- read_rds("Data/unpd_deaths_raw.rds")

dts <- 
  DT %>% 
  as_tibble() %>% 
  select(LocID, 
         Country = LocName, 
         Date = TimeStart, 
         Date2 = TimeEnd,
         unit = TimeUnit,
         dura = TimeDuration,
         Sex = SexName, 
         Age = AgeStart, 
         AgeSpan, 
         AgeLabel,
         Deaths = DataValue,
         Source = DataSourceName,
         src_yr = DataSourceYear) %>% 
  mutate(Date = dmy(Date),
         Date2 = dmy(Date2),
         Year = year(Date),
         Year2 = year(Date2),
         Code = countrycode::countrycode(LocID, "un", "iso3c"),
         Sex = case_when(Sex == "Male" ~ "m",
                         Sex == "Female" ~ "f",
                         Sex == "Both sexes" ~ "t",
                         Sex == "Unknown" ~ "unk"),
         Country = recode(Country,
                          "Czech Republic" = "Czechia",
                          "Iran (Islamic Republic of)" = "Iran",
                          "Bolivia (Plurinational State of)" = "Bolivia",
                          "United States of America" = "USA",
                          "China, Hong Kong SAR" = "Hong Kong",
                          "Republic of Korea" = "South Korea",
                          "Russian Federation" = "Russia",
                          "Republic of Moldova" = "Moldova")) %>% 
  drop_na(Sex) %>% 
  select(-Date, -LocID) %>% 
  # only annual data
  filter(unit == "year" & dura == 1) %>% 
  select(-Year2, -unit, -dura, -Date2) %>% 
  # only sources with data in 2020 or 2021
  group_by(Country, Source) %>% 
  filter(max(Year) >= 2020) %>% 
  ungroup() %>% 
  # adjusting wrong labels
  mutate(AgeLabel = case_when(
    Age == 100 & AgeSpan == 5 & AgeLabel == "99-104" ~ "100-104",
    Age == 105 & AgeSpan == 5 & AgeLabel == "104-109"  ~ "105-109",
    TRUE ~ AgeLabel)) %>% 
  # removing duplicate values
  unique() %>% 
  group_by(Country, Code, Source, Year, Sex, Age, AgeLabel, AgeSpan) %>% 
  # most recent update of the same source
  filter(src_yr == max(src_yr)) %>% 
  # highest value of deaths
  filter(Deaths == max(Deaths)) %>% 
  ungroup() %>% 
  # excluding HMD, as we are collecting them directly
  filter(Source != "Human Mortality Database") %>% 
  select(-src_yr)

# testing duplicates
dts %>% 
  group_by(Country, Code, Source, Year, Sex, Age, AgeLabel, AgeSpan) %>% 
  filter(n()>1) %>% 
  ungroup()

# test for data with wrong labels 
dts %>% 
  group_by(Country, Code, Source, Year, Sex, Age, AgeSpan) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n > 1)

unique(dts$AgeLabel) %>% sort

# sources and data series
# ~~~~~~~~~~~~~~~~~~~~~~~
srs <- 
  DT %>% 
  select(LocName, DataSourceAuthor, DataSourceName, DataSourceShortName) %>% 
  unique

# ~~~~~~~~~~~~~~~~~~~~~
# Data adjustments ====
# ~~~~~~~~~~~~~~~~~~~~~

# 1. adjust for open age interval (keeping the oldest)
# 2. adjust infant and child ages
# 3. removing incomplete ages
# 4. adjust for totals by age and re-scale
# 5. adjust for totals by sex and re-scale
# 6. harmonize ages (same intervals for all series, max closing at 100+)
# 7. remove series with insufficient periods for baseline (min 3 years 2015-2019)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. adjust for open age interval ==============================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# choosing only one open age interval (highest)
open_ages <- 
  dts %>% 
  filter(AgeSpan == -1 & AgeLabel != "Total") %>% 
  group_by(Country, Code, Source, Year, Sex) %>% 
  filter(Age == max(Age)) %>% 
  ungroup()

# no open age interval
no_open <- 
  dts %>% 
  anti_join(open_ages %>% 
              select(Country, Source, Sex, Year) %>% 
              unique()) %>% 
  select(Country, Source, Year, Sex) %>% 
  unique()

# unique open age
dts2 <- 
  dts %>% 
  # removing unknown age, total age, and open ages
  filter(AgeSpan >= 0) %>% 
  # adjusting open age interval
  left_join(open_ages %>% 
              select(Country, Source, Year, Sex, max_age = Age)) %>% 
  filter(Age < max_age | is.na(max_age)) %>% 
  # adding open age intervals 
  bind_rows(open_ages) %>% 
  arrange(Country, Source, Year, Sex, Age) %>% 
  select(-max_age) %>% 
  anti_join(no_open)

# test for duplicates 
dts2 %>% 
  group_by(Country, Code, Source, Year, Sex, Age, AgeSpan) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n > 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. adjust infant and child ages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# all child ages in infant and 1-4, or 0-4
child_ages <-
  dts2 %>% 
  filter(Age %in% 0:4,
         AgeSpan %in% 1:5) %>% 
  mutate(Age = ifelse(Age %in% 1:4 & AgeSpan == 1, 1, Age),
         AgeLabel = ifelse(Age %in% 1:4 & AgeSpan == 1, "1-4", AgeLabel),
         AgeSpan = ifelse(Age %in% 1:4 & AgeSpan == 1, 4, AgeSpan)) %>% 
  group_by(Country, Code, Source, Year, Sex, Age, AgeLabel, AgeSpan) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  group_by(Country, Code, Source, Year, Sex) %>% 
  mutate(has_inf = ifelse(any(AgeLabel == "< 1"), 1, 0),
         has_1_4 = ifelse(any(AgeLabel == "1-4"), 1, 0),
         has_0_4 = ifelse(any(AgeLabel == "0-4"), 1, 0)) %>% 
  ungroup()

child_ages %>% 
  group_by(Country, Code, Source, Year, Sex, Age, AgeSpan) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n > 1)


# combinations to exclude
# ~~~~~~~~~~~~~~~~~~~~~~~
# estimating age 1-4 for those with infant and 0-4 but no 1-4
needed_1_4 <- 
  child_ages %>% 
  filter(has_inf == 1 & has_1_4 == 0 & has_0_4 == 1) %>% 
  select(-starts_with("has_"), -AgeSpan) %>% 
  spread(AgeLabel, Deaths) %>% 
  mutate(`1-4` = `0-4` - `< 1`) %>% 
  gather(`1-4`, `0-4`, `< 1`, key = AgeLabel, value = Deaths) %>% 
  mutate(AgeSpan = case_when(AgeLabel == "< 1" ~ 1,
                             AgeLabel == "1-4" ~ 4,
                             AgeLabel == "0-4" ~ 5)) %>% 
  filter(AgeLabel == "1-4") %>% 
  mutate(Age = 1)

# to exclude with ages 1-4 but no infant nor 0-4
missing_inf <- 
  child_ages %>% 
  filter(has_inf == 0 & has_1_4 == 1 & has_0_4 == 0) %>% 
  select(Country, Source, Year, Sex) %>% 
  unique()

# to exclude with infant but no 1-4 nor 0-4
only_inf <- 
  child_ages %>% 
  filter(has_inf == 1 & has_1_4 == 0 & has_0_4 == 0) %>% 
  select(Country, Source, Year, Sex) %>% 
  unique()

# putting all together
dts_inf_chld <- 
  child_ages %>% 
  select(-starts_with("has")) %>% 
  bind_rows(needed_1_4) %>% 
  anti_join(missing_inf) %>% 
  anti_join(only_inf) %>% 
  arrange(Country, Source, Year, Sex, Age)

# unnecessary ages 0-4
unnecess_0_4 <- 
  dts_inf_chld %>% 
  group_by(Country, Code, Source, Year, Sex) %>% 
  mutate(has_inf = ifelse(any(AgeLabel == "< 1"), 1, 0),
         has_1_4 = ifelse(any(AgeLabel == "1-4"), 1, 0),
         has_0_4 = ifelse(any(AgeLabel == "0-4"), 1, 0)) %>% 
  ungroup() %>% 
  filter(has_inf == 1 & has_1_4 == 1 & has_0_4 == 1) %>% 
  filter(AgeLabel == "0-4") %>% 
  select(Country, Source, Year, Sex, Age, AgeLabel) %>% 
  unique()

dts_inf_chld2 <- 
  dts_inf_chld %>% 
  anti_join(unnecess_0_4) 

# test of completeness
dts_inf_chld2 %>% 
  group_by(Country, Code, Source, Year, Sex) %>% 
  mutate(span_inf = ifelse(any(AgeLabel == "< 1"), 1, 0),
         span_1_4 = ifelse(any(AgeLabel == "1-4"), 4, 0),
         span_0_4 = ifelse(any(AgeLabel == "0-4"), 5, 0),
         test = span_inf + span_1_4 + span_0_4) %>% 
  ungroup() %>% 
  filter(test != 5)

# no age zero
no_child <- 
  dts2 %>% 
  anti_join(dts_inf_chld2 %>% 
              select(Country, Source, Sex, Year) %>% 
              unique()) %>% 
  select(Country, Source, Year, Sex) %>% 
  unique()

# excluding countries without age 0
dts3 <- 
  dts2 %>% 
  filter(Age >= 5) %>% 
  bind_rows(dts_inf_chld2) %>% 
  arrange(Country, Source, Year, Sex, Age) %>% 
  anti_join(no_child)

# test for duplicates 
dts3 %>% 
  group_by(Country, Code, Source, Year, Sex, Age, AgeSpan) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n > 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. removing incomplete ages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test <- 
  dts3 %>% 
  filter(Year < 2020) %>% 
  group_by(Country, Code, Source, Sex, Age) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  select(Country, Code, Source, Sex, n) %>% 
  unique() %>% 
  group_by(Country, Code, Source, Sex) %>% 
  summarise(grs = n()) %>% ungroup() %>% 
  filter(grs>1)
  
incomp_ages <- 
  dts3 %>% 
  group_by(Country, Code, Source, Year, Sex) %>% 
  mutate(AgeSpan = ifelse(Age == max(Age), 0, AgeSpan)) %>% 
  filter(max(Age) != sum(AgeSpan)) %>% 
  select(Country, Source, Year, Sex) %>% 
  unique()
  
dts4 <- 
  dts3 %>% 
  anti_join(incomp_ages)

# test for duplicates 
dts4 %>% 
  group_by(Country, Code, Source, Year, Sex, Age, AgeSpan) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n > 1)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. adjust for totals by age and re-scale
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# included series
inc4 <- 
  dts4 %>% 
  select(Country, Code, Source, Year, Sex) %>% 
  unique()

# total ages
# ~~~~~~~~~~
# total ages already in data
age_tot_in <- 
  dts %>% 
  inner_join(inc4) %>%
  filter(AgeLabel == "Total") %>% 
  mutate(Age = "TOT") %>% 
  select(-AgeLabel)

# unknown ages (for combinations without totals)
# ~~~~~~~~~~~~
age_unk_in <- 
  dts %>% 
  inner_join(inc4) %>%
  filter(AgeLabel == "Unknown") %>% 
  filter(Deaths > 0) %>% 
  anti_join(age_tot_in %>% select(Country, Sex, Source, Year))

# no included in data
age_tot_no <- 
  dts4 %>% 
  anti_join(age_tot_in %>% select(Country, Source, Sex, Year) %>% unique()) %>% 
  bind_rows(age_unk_in) %>% 
  group_by(Country, Code, Source, Sex, Year) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  mutate(Age = "TOT",
         AgeSpan = -1) %>% 
  ungroup()

# putting total ages together
tots_all <- 
  bind_rows(age_tot_in, age_tot_no) %>% 
  arrange(Country, Code, Source, Sex, Year)

# putting all together
dts5 <- 
  dts4 %>% 
  mutate(Age = Age %>% as.character()) %>% 
  bind_rows(tots_all) %>% 
  arrange(Country, Code, Source, Year, Sex, suppressWarnings(as.integer(Age))) 

# re-scalling age
dts6 <- 
  dts5 %>% 
  group_by(Country, Source, Sex, Year) %>% 
  do(rescale_age(chunk = .data)) %>% 
  ungroup() %>% 
  mutate(Age = Age %>% as.double())

unique(dts6$Age)

# test for duplicates 
dts6 %>% 
  group_by(Country, Code, Source, Year, Sex, Age, AgeSpan) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n > 1)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. adjust for totals by sex and re-scale
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# total sex included series
tot_sex_in <- 
  dts6 %>% 
  filter(Sex == "t") 

# total sex for those missing
tot_sex_no <- 
  dts6 %>% 
  anti_join(tot_sex_in %>% 
              select(Country, Code, Source, Year, Age) %>% 
              unique()) %>% 
  # filter(Sex != "t") %>% 
  group_by(Country, Code, Source, Year, Age, AgeSpan, AgeLabel) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup() %>% 
  mutate(Sex = "t")

# adding total sex
dts7 <- 
  dts6 %>% 
  filter(Sex != "unk" & Sex != "t") %>% 
  bind_rows(tot_sex_in, tot_sex_no) %>% 
  arrange(Country, Code, Source, Year, Sex, suppressWarnings(as.integer(Age))) 

# only with total sex (to exclude)
only_t_sex <- 
  dts7 %>% 
  select(Country, Code, Source, Year, Sex) %>% 
  unique() %>% 
  group_by(Country, Code, Source, Year) %>% 
  filter(n() != 3) %>% 
  select(-Sex)

# imputation of unknown sex
# ~~~~~~~~~~~~~~~~~~~~~~~~~
dts8 <-
  dts7 %>% 
  anti_join(only_t_sex) %>% 
  group_by(Country, Source, Age, Year) %>%
  do(rescale_sex(chunk = .data)) %>%
  ungroup()

# adding those with only total sex
dts9 <- 
  dts8 %>% 
  bind_rows(dts7 %>% 
              inner_join(only_t_sex)) %>% 
  arrange(Country, Code, Source, Year, Sex, Age) %>% 
  replace_na(list(Deaths = 0))

unique(dts9$Sex)

# test for duplicates 
dts9 %>% 
  group_by(Country, Code, Source, Year, Sex, Age, AgeSpan) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n > 1)

write_rds(dts9, "data_inter/unpd_harmonized_test.rds")
dts9 <- read_rds("data_inter/unpd_harmonized_test.rds")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6. harmonize ages (same open age for all series, max closing at 100+)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
max_age <- 
  dts9 %>% 
  group_by(Country, Code, Source, Sex) %>% 
  summarise(max_age = max(Age)) %>% 
  mutate(max_age = min(max_age, 100)) %>% 
  ungroup() 

max_age <- 
  dts9 %>% 
  group_by(Country, Code, Source, Sex, Year) %>% 
  summarise(max_age = max(Age)) %>% 
  group_by(Country, Code, Source, Sex) %>% 
  summarise(max_age = min(max_age)) %>% 
  mutate(max_age = min(max_age, 100)) %>% 
  ungroup() 

dts10 <- 
  dts9 %>% 
  left_join(max_age) %>% 
  mutate(Age = ifelse(Age < max_age, Age, max_age),
         AgeLabel = ifelse(Age == max_age, paste0(max_age, "+"), AgeLabel),
         AgeSpan = ifelse(Age == max_age, -1, AgeSpan)) %>% 
  group_by(Country, Code, Source, Year, Sex, Age, AgeLabel, AgeSpan) %>% 
  summarise(Deaths = sum(Deaths)) %>% 
  ungroup()

# test for duplicates 
dts10 %>% 
  group_by(Country, Code, Source, Year, Sex, Age, AgeSpan) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n > 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7. looking at age groups differences
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

diff_age_grs <- 
  dts10 %>% 
  filter(Year < 2020) %>% 
  group_by(Country, Code, Source, Sex, Age) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  select(Country, Code, Source, Sex, n) %>% 
  unique() %>% 
  group_by(Country, Code, Source, Sex) %>% 
  summarise(grs = n()) %>% ungroup() %>% 
  filter(grs>1) %>% 
  select(-grs)
 
dts11 <- 
  dts10 %>% 
  anti_join(diff_age_grs)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 8. remove series with insufficient periods for baseline (min 3 years 2015-2019)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# series with sufficient years for estimating the baseline (at least 3 within 2015-2019)
# series with data during the pandemic 2020-2021 (at least 1)

unique(dts11$Source)

enough_pers <- 
  dts11 %>% 
  select(Country, Code, Source, Year, Sex) %>% 
  unique() %>% 
  mutate(pre = ifelse(Year %in% 2015:2019, 1, 0),
         pan = ifelse(Year %in% 2020:2021, 1, 0)) %>% 
  group_by(Country, Code, Source, Sex) %>% 
  filter(sum(pre) >= 3 & sum(pan) >= 1) %>% 
  ungroup() %>% 
  select(-pre, -pan)

dts12 <- 
  dts11 %>% 
  inner_join(enough_pers) %>% 
  group_by(Country, Source, Year, Sex) %>% 
  mutate(age_spn = ifelse(Age == max(Age), -1, lead(Age) - Age)) %>% 
  ungroup() %>% 
  mutate(Source = case_when(Source == "Demographic Yearbook" ~ "unpd_dy",
                            # Source == "WHO All-Cause Mortality Data Call" ~ "unpd_who",
                            str_detect(Source, "WHO") ~ "unpd_who",
                            Source == "Human Mortality Database" ~ "unpd_hmd",
                            TRUE ~ "unpd_crvs")) %>% 
  replace_na(list(Deaths = 0))

unique(dts12$Source)

# test for duplicates
dts12 %>%
  # select(Country, Code, Source2, Year, Sex, Age, age_up, Deaths) %>%
  group_by(Country, Code, Source, Year, Sex, Age, age_spn) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n > 1)

dts13 <-
  dts12 %>%
  group_by(Country, Code, Source, Year, Sex, Age, age_spn) %>%
  filter(Deaths == max(Deaths)) %>%
  select(-AgeLabel, -AgeSpan)

dts13 %>%
  group_by(Country, Code, Source, Year, Sex, Age, age_spn) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n > 1)

unique(dts13$Country)

incomp_ages <- 
  dts13 %>% 
  group_by(Country, Code, Source, Year, Sex) %>% 
  filter(max(Age) != sum(age_spn)+1) %>% 
  select(Country, Source, Year, Sex) %>% 
  unique()




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

readr::write_csv(dts13, file = "data_inter/unpd.csv")
readr::write_rds(dts13, file = "data_inter/unpd.rds")

# dts12 <- read_rds("data_inter/unpd.rds") %>%
#   select(-AgeSpan, -pre, -pan)
