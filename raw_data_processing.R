### Prepping data files for NES-LTER data package -- this script manipulates the raw data to include correct column names and categories for downstream processing.
## Author: Diana Fontaine
# Input: 1) "raw_npp_data.csv" files dependent on version. These files have the raw data from the analysis lab; 2) light_data.csv that includes % PAR data for each station per cruise. Version 1: EN644 - EN649 and Version 2: EN655 - EN687
# Output: "raw_data_table.csv" with appropriate version appended to file name. These outputs are used for downstream npp_rate calculations.


#load libaries
library(tidyverse)

#load light data
light <- read_csv("light_data.csv")
npp_raw_vers1 <- read_csv("raw_npp_data_version1.csv")
npp_raw_vers2 <- read_csv("raw_npp_data_version2.csv")

#bind together two versions to be able to do the same thing
npp_raw <- rbind(npp_raw_vers1, npp_raw_vers2)


#make npp_raw cast numeric
npp_raw$cast <- as.numeric(npp_raw$cast)

#make light cast numeric
light$cast <- as.numeric(light$cast)
light$station <- as.character(light$station)

#Merge to get light data in raw file
npp_raw <- npp_raw %>%
  left_join(light, by = c("cruise", "station", "cast", "niskin"))

#make surface irradiance in a % by multiplying by 100
npp_raw$percent_surface_irradiance <- npp_raw$percent_surface_irradiance*100

#separate alternate_id column
npp_raw <- npp_raw %>%
  separate(alternate_sample_ID, into = c("depth_category", "filter_size"), sep = "_", remove = FALSE)

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

#save the two versions of the files
npp_table_version1 <- npp_raw %>%
  filter(cruise %in% c("EN644", "AR39B", "EN649"))

npp_table_version2 <- npp_raw %>%
  filter(!(cruise %in% c("EN644", "AR39B", "EN649")))

#Save as "raw_data_table.csv"-- make sure to save the correct versions
write_csv(npp_table_version1, "For_EDI/raw_data_table_version1.csv") 
write_csv(npp_table_version2, "For_EDI/raw_data_table_version2.csv") 






