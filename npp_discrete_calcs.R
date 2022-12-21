### Calculating discrete depth estimates of NPP
## Author: Diana Fontaine
# Input is raw_data_table.csv with version appended to file name. Version 1: EN644 - EN649 and Version 2: EN655 - EN687. These files were produced from the upstream "raw_data_processing.R" file
# Output is csv file with pp_rates and associated sample information: npp_discrete.csv with version appended to file name.

#load libraries
library(tidyverse)
library(geosphere)
#read in data
npp_discrete1 <- read_csv("For_EDI/raw_data_table_version1.csv")
npp_discrete2 <- read_csv("For_EDI/raw_data_table_version2.csv")

npp_discrete <- rbind(npp_discrete1,npp_discrete2)

#need to change L1.2 and L4.2  to be L1 and L4, respectively
npp_discrete$station[npp_discrete$station == "L1.2"] <- "L1"
npp_discrete$station[npp_discrete$station == "L4.2"] <- "L4"


#filter out "D4_10_exclude" row
npp_discrete <- npp_discrete %>%
  filter(alternate_sample_category != "D4_10_exclude")

#Make sample ID character
npp_discrete$alternate_sample_category <- as.character(npp_discrete$alternate_sample_category)

#Make cruise a factor
npp_discrete$cruise <- as.factor(npp_discrete$cruise)


#filter for the sample type of "filter blank" so that I can subtract the blank POC value from the measured values to get a "blank corrected POC" value
filt_blank <- npp_discrete %>%
  filter(alternate_sample_category == "Filter_blank")


#there was a contamination issue with the filter blanks from cruise EN644 so these filter blanks were not included in the average. Instead, I averaged the filter blank values from the other cruises and used this value for the EN644 cruise
#this does not include EN644 in the blank averages because of the contamination
filt_blank_noen644 <- filt_blank %>%
  filter(cruise != "EN644") %>%
  dplyr::group_by(cruise) %>%
  dplyr::summarise(mean_blank = mean(massC_mg, na.rm = TRUE)) 

#Add filter_blank values for EN644 cruise that had the contamination issue -- this is just an average of the other cruises without EN644
filt_blank <- filt_blank_noen644 %>%
  add_row(cruise = "EN644", mean_blank = mean(filt_blank_noen644$mean_blank))  

 
#add a column of blank corrected POC values
npp_discrete <- npp_discrete %>%
  filter(alternate_sample_category != "Filter_blank") %>%
  left_join(filt_blank, by = "cruise") %>%
  mutate(blank_cor_POC_mg = massC_mg - mean_blank) %>%
  mutate(POC_ug = blank_cor_POC_mg * 1000) #multiply the POC by 1000 to get into µg

#determine the amount of NaHCO3 per bottle
npp_discrete <- npp_discrete %>%
  mutate(mg_bottle = g_mL_NaHCO3*vol_added_mL*1000)

#get the POC in units of µg/L
npp_discrete <- npp_discrete %>%
  mutate(POC_ug_L = POC_ug/vol_filt_L) %>% 
  mutate(inc_time_day = inc_time_min/1440) #Convert incubation time from minutes to day

#need to convert d13C values to atom percent to ultimately calculate NPP
VPDB <- 0.0111803 #This is the absolute ratio for 13C according to the VPDB standard
npp_discrete <- npp_discrete %>%
  mutate(at_per = 100*VPDB*((d13C/1000+1))/(1+(VPDB*(d13C/1000+1))))

#need to get the natural abundance samples by themselves
nat <- npp_discrete[str_detect(npp_discrete$alternate_sample_category, "NatAbun"), ]

#need to add a row for L7bloom natabun sample since this wasn't collected (I decided to use the natural abundance sample from L6)
#nat_L6en668 <- nat %>%
 # filter(cruise == "EN668") %>%
  #filter(station == 6) %>%
  #select(at_per)
#nat_L6en668 <- as.numeric(nat_L6en668)

#select cast
#nat_L6en668_cast <- nat %>%
 # filter(cruise == "EN668") %>%
 # filter(station == 6) %>%
 # select(cast)
#nat_L6en668_cast <- as.numeric(nat_L6en668_cast)


#just need certain columns
nat <- nat %>%
 select(cruise, cast, station, at_per)

#nat <- nat %>%
  #add_row(cruise = "EN668", station = "7bloom", cast = nat_L6en668_cast, at_per =  nat_L6en668)



#join raw_pp df with nat df so that each sample has its respective natural abundance values (this value is used in the rate calculation)
pp <- left_join(npp_discrete, nat, by = c('cruise', 'cast', 'station')) %>%
  rename(at_per = at_per.x,
         at_nat = at_per.y)#rename columns so not to get confused between "x" and "y"

#for analysis, will need salinity data from the NES LTER data downloaded directly from the REST API
bottles_en644 <- read_csv('https://nes-lter-data.whoi.edu/api/ctd/en644/bottles.csv')
bottles_ar39 <- read_csv('https://nes-lter-data.whoi.edu/api/ctd/ar39b/bottles.csv')
bottles_en649 <- read_csv('https://nes-lter-data.whoi.edu/api/ctd/en649/bottles.csv')
bottles_en655 <- read_csv('https://nes-lter-data.whoi.edu/api/ctd/en655/bottles.csv')
bottles_en657 <- read_csv('https://nes-lter-data.whoi.edu/api/ctd/en657/bottles.csv')
bottles_en661 <- read.csv('https://nes-lter-data.whoi.edu/api/ctd/en661/bottles.csv')
bottles_en668 <- read.csv('https://nes-lter-data.whoi.edu/api/ctd/en668/bottles.csv')
bottles_ar61 <- read.csv('https://nes-lter-data.whoi.edu/api/ctd/ar61b/bottles.csv')
bottles_at46 <- read.csv('https://nes-lter-data.whoi.edu/api/ctd/at46/bottles.csv')
bottles_en687 <- read.csv('https://nes-lter-data.whoi.edu/api/ctd/en687/bottles.csv')



#just need certain columns from the bottles dataframe and some of these column names are hard to remember so change them to something better. do this for each cruise
bottles_en644 <- bottles_en644 %>%
  select(date, cruise, cast, niskin, latitude, longitude, prdm, sal00, t090c) %>%
  rename(depth = prdm,
         salinity = sal00,
         temp = t090c)

bottles_ar39 <- bottles_ar39 %>%
  select(date,cruise, cast, niskin, latitude, longitude, prdm, sal00, t090c) %>%
  rename(depth = prdm,
         salinity = sal00,
         temp = t090c)

bottles_en649 <- bottles_en649 %>%
  select(date,cruise, cast, niskin, latitude, longitude, prdm, sal00, t090c) %>%
  rename(depth = prdm,
         salinity = sal00,
         temp = t090c) 

bottles_en655 <- bottles_en655 %>%
  select(date,cruise, cast, niskin, latitude, longitude, prdm, sal00, t090c) %>%
  rename(depth = prdm,
         salinity = sal00,
         temp = t090c) 

bottles_en657 <- bottles_en657 %>%
  select(date,cruise, cast, niskin, latitude, longitude, prdm, sal00, t090c) %>%
  rename(depth = prdm,
         salinity = sal00,
         temp = t090c) 

bottles_en661 <- bottles_en661 %>%
  select(date, cruise, cast, niskin, latitude, longitude, prdm, sal00, t090c) %>%
  rename(depth = prdm,
         salinity = sal00,
         temp = t090c)

bottles_en668 <- bottles_en668 %>%
  select(date, cruise, cast, niskin, latitude, longitude, prdm, sal00, t090c) %>%
  rename(depth = prdm,
         salinity = sal00,
         temp = t090c)

bottles_ar61 <- bottles_ar61 %>%
  select(date, cruise, cast, niskin, latitude, longitude, prdm, sal00, t090c) %>%
  rename(depth = prdm,
         salinity = sal00,
         temp = t090c) 

bottles_at46 <- bottles_at46 %>%
  select(date, cruise, cast, niskin, latitude, longitude, prdm, sal00, t090c) %>%
  rename(depth = prdm,
         salinity = sal00,
         temp = t090c) 

bottles_en687 <- bottles_en687 %>%
  select(date, cruise, cast, niskin, latitude, longitude, prdm, sal00, t090c) %>%
  rename(depth = prdm,
         salinity = sal00,
         temp = t090c) 


#join bottle data for en644 and ar39 using rbind
bottles_join <- rbind(bottles_en644, bottles_ar39, bottles_en649, bottles_en655, bottles_en657, bottles_en661, bottles_en668,bottles_ar61, bottles_at46, bottles_en687)

#make cast for pp dataframe numerical for merging
pp$cast <- as.numeric(pp$cast)


#merge npp discrete and bottles_join df so that the salinity data are included in the pp data (need salinity for analysis)
#joining the pp df and the bottles_join df
pp_cruises <- pp %>%
  left_join(bottles_join, by = c('cruise', 'cast', 'niskin'))


#calculations to get to pp rates
npp_calcs <- pp_cruises %>%
  mutate(DIC = (((salinity*0.067)-0.05)*0.96)*12000) %>% #DIC units are mg/m3 or ug/L
  mutate(added_13C = (mg_bottle*1000*0.15294)/bottle_vol_L) %>% #Now add column for *******13C added***** in ugC/L by: (mg/sample column*1000 to get to ug*0.15294 which is the mass of 13C/mass of the entire compound (NaHCO3) and then divide all of this by the bottle volume to get added_13C in 
  mutate(at_13C = (added_13C*100)/(DIC + added_13C)) %>% #Add column for at%13C in the DIC pool by: (13Cadded column *100)/(DIC+13C added). the 100 is to get the value into a percent
  mutate(pp_rate = ((POC_ug_L*(at_per - at_nat))/(inc_time_day*(at_13C - at_nat)))*1.025) #now get the production rate in ugC/L/day

#add iode quality flags to negative values that were negative because of blank corrected POC values
npp_calcs$iode_quality_flag <-with(npp_calcs,
                                   ifelse(pp_rate < 0, "3",
                                          ifelse(pp_rate >0, "1",NA))) #Note there will be no quality flag for natural abundance samples

##D1 values from EN668 cast 20 (d6a) did not have 200um mesh (wasn't sampled properly). make the qual flag a 4 and make the values an NA
npp_en668_D1_forqual <- npp_calcs %>%
  filter(cruise == "EN668" & depth_category == "D1" & cast == "20")

npp_en668_D1_forqual$iode_quality_flag <- 4
npp_en668_D1_forqual$pp_rate <- NA


#rbind dataframe back together and filter out D4 10 exclude row and en668 d6a (cast 20) D1 samples since we binded them together (we don't want duplicate rows for the same samples)
npp_calcs_joined <- rbind(npp_calcs, npp_en668_D1_forqual)

npp_calcs_joined <- npp_calcs_joined %>%
  filter(!(cruise == "EN668" & cast == "20" & depth_category == "D1" & iode_quality_flag == 1))
  

### Using list made above to obtain nearest station info for discrete depth table
npp_dis <- npp_calcs_joined
# initialize nearest station and distance columns to NA
npp_dis$nearest_station <- NA_character_
npp_dis$distance <- NA_integer_
# read list csv into stations. using 668 station list because this has all the stations whereas the other master one doesn't
stations <- read_csv('https://nes-lter-data.whoi.edu/api/stations/en668.csv')
station_matrix <- matrix(data = c(stations$longitude, stations$latitude), nrow = 25, ncol = 2, byrow = FALSE,
                         dimnames = NULL)
# calculate distance per row of the data frame
for (df_row in 1:nrow(npp_dis)) {
  df_lon <- npp_dis$longitude[df_row]
  df_lat <- npp_dis$latitude[df_row]
  # add an if to skip the row if df lon and/or lat is NA
  if (!is.na(df_lon) & !is.na(df_lat)) {
    df_lon_lat <- c(df_lon,df_lat)
    km_from_df <- distHaversine(station_matrix, df_lon_lat, r=6378.137)
    # index the minimum distance
    index <- which.min(km_from_df)
    # use that index to pull the station name and its distance
    nearest_station_list <- stations[index,'name']
    # need to change this from a list to char
    nearest_station <- unname(unlist(nearest_station_list))
    distance <- km_from_df[index]
    # If distance less than 2 km, use base R to add to respective columns in the full data table within the for loop
    #if (distance < 3) {
      npp_dis$nearest_station[df_row] <- nearest_station
      npp_dis$distance[df_row] <- distance
  }
}
#Round distance column
npp_dis$distance <- round(npp_dis$distance, 2)

#Change date column and rename pp_rate columnn
npp_dis <- npp_dis %>%
  rename(date_time_utc = date) %>%
  rename(npp_rate = pp_rate)

#Round npp rate to 3 numbers after the decimal
npp_dis$npp_rate <- round(npp_dis$npp_rate, digits = 3)
#Round lat and long to 4 numbers after the decimal
npp_dis$latitude <- round(npp_dis$latitude, digits = 4)
npp_dis$longitude <- round(npp_dis$longitude, digits = 4)

#filter out Natabun samples
npp_dis <- npp_dis %>%
  filter(depth_category != "NatAbun")


#check to make sure station and nearest station match up.
npp_dis$stncheck <- ifelse(npp_dis$station==npp_dis$nearest_station,"TRUE","FALSE")
#Here, the false ones are for  L8 when ship drifted closer to u9a. The u9a stations will be manually changed in output table

#check station distance 
npp_dis$discheck <- ifelse(npp_dis$distance < 2.5,"TRUE","FALSE")
#range of rows classified as false is from 2.53 - 3.94. result of ship drifting


#select just what needs to be in the final datasheet
npp_calcs_final <- npp_dis %>%
  select(cruise, date_time_utc, cast, station, niskin,  latitude, longitude, depth, alternate_sample_category, depth_category, filter_size, replicate, percent_surface_irradiance, npp_rate, iode_quality_flag, nearest_station, distance)


#fix datetime
npp_calcs_final$date_time_utc <- str_replace(npp_calcs_final$date_time_utc, "T", " ")

#need to add new column called "incub_tank" for the experimental samples from en687
## for 26 degree experiment
en687exp_26 <- npp_calcs_final %>%
  filter(grepl("tank", alternate_sample_category)) %>%
  filter(grepl("tank3", alternate_sample_category))  %>%
  mutate(incub_type = "experiment_26deg")

#for 18 degree exp
en687exp_18 <- npp_calcs_final %>%
  filter(grepl("tank2", alternate_sample_category))  %>%
  mutate(incub_type = "experiment_18deg")

#for 14 degree exp
en687exp_14 <- npp_calcs_final %>%
  filter(grepl("tank1", alternate_sample_category)) %>%
  filter(nearest_station %in% c("L6", "L11")) %>%
  mutate(incub_type = "experiment_14deg")

#all the rest that need just "ambient_temp" in the incubation type column
no_experimental <- npp_calcs_final %>%
  filter(!grepl("tank", alternate_sample_category)) %>%
  mutate(incub_type = "ambient_temp")


#bind them back together
npp_calcs_final_withincub <- rbind(no_experimental, en687exp_14, en687exp_18, en687exp_26)


#selecting for version1
npp_calcs_final_vers1 <- npp_calcs_final_withincub %>%
  filter(cruise %in% c("EN644", "AR39B", "EN649"))

#selecting for version2
npp_calcs_final_vers2 <- npp_calcs_final_withincub%>%
  filter(cruise %in% c("EN655","EN657", "EN661", "EN668", "AR61B", "AT46", "EN687"))




#write_csv
write_csv(npp_calcs_final_vers1, "For_EDI/npp_discrete_version1.csv")
write_csv(npp_calcs_final_vers2, "For_EDI/npp_discrete_version2.csv")













