### Prepping data files for NES-LTER data package -- this script manipulates the raw data to include correct column names and categories for downstream processing.
## Author: Diana Fontaine
# Input: 1) raw_npp_data.csv that has data from the analysis lab; 2) light_data.csv that includes % PAR data for each station per cruise
# Output: raw_data_table.csv that is used for downstream pp_rate calculations.


#load libaries
library(tidyverse)

#load data
light <- read_csv("light_data.csv")
npp_raw <- read_csv("raw_npp_data.csv")

#make npp_raw cast numeric
npp_raw$cast <- as.numeric(npp_raw$cast)

#make light cast numeric
light$cast <- as.numeric(light$cast)

#Merge to get light data in raw file
npp_raw <- npp_raw %>%
  left_join(light, by = c("cruise", "station", "cast", "niskin"))

#make surface irradiance in a % by multiplying by 100
npp_raw$percent_surface_irradiance <- npp_raw$percent_surface_irradiance*100

#separate alternate_id column
npp_raw <- npp_raw %>%
  separate(alternate_sample_category, into = c("depth_category", "filter_size"), sep = "_", remove = FALSE)

#all occurence of "GFF" under filter_size should be >0
npp_raw$filter_size[npp_raw$filter_size=='GFF'] <- ">0"
#adjust filter size for <5
npp_raw$filter_size[npp_raw$filter_size=='20'] <- ">0&<20"
#for 20um
npp_raw$filter_size[npp_raw$filter_size=='5'] <- ">0&<5"
#for 10um
npp_raw$filter_size[npp_raw$filter_size=='10'] <- ">0&<10"

#Change filter size for filter blanks
npp_raw$filter_size[npp_raw$filter_size=='blank'] <- ">0"

#Change depth category for filter blanks
npp_raw$depth_category[npp_raw$depth_category=='Filter'] <- "NA"

#Change filter size for Dark and NatAbun samples
npp_raw$filter_size[npp_raw$filter_size=='NatAbun'] <- ">0"

#Remove date column
npp_raw <- npp_raw %>%
  select(-date)

#Save as "raw_data_table.csv"
write_csv(npp_raw, "For_EDI/raw_data_table.csv")






