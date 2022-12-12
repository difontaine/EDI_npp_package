### Generate integrated npp rates
# Input: npp_discrete.csv that has the calculated pp_rates by depth and it an upstream product from the "npp_discrete_calcs.R" file. Note that there are two versions to this file. Version 1: EN644 - EN649 and Version 2: EN655 - EN687.
# Output: integrated_pp.csv (with versions) which is the size-fractionated integrated npp rates per station per cruise and has associated extinction coefficient data

#load libraries
library(tidyverse)
library(dplyr)


#read in data
pp_vers1 <- read_csv("For_EDI/npp_discrete_version1.csv")
pp_vers2 <- read_csv("For_EDI/npp_discrete_version2.csv")

pp_for_int <- rbind(pp_vers1, pp_vers2)

pp_for_int <- pp_for_int %>%
  filter(iode_quality_flag <4)

#Create a no reps df and averaged by replicate
pp_filt <- pp_for_int %>%
  filter(!grepl('NatAbun', alternate_sample_category)) %>% #filter out NatAbun rows 
  filter(!grepl('Dark', alternate_sample_category)) %>% #filter out Dark rows 
  select(date_time_utc, cruise, cast, nearest_station, latitude, longitude, niskin, alternate_sample_category, filter_size, replicate, depth, depth_category, npp_rate)

#npp_discrete data file doesn't have nearest station fo EN644 L8 so need to add that in somehow 
#can do it from light csv

light <- read_csv("light_data.csv")
#add "L" before each station value
light$cast <- as.numeric(light$cast)

#join light data and pp_filt together
pp_filt <- pp_filt %>%
  left_join(light, by = c("cruise", "cast", "niskin"))
  

#Now get averaging for those samples that have replicates
pp_filt_noreps <- pp_filt %>%
  group_by(cruise, cast, niskin, nearest_station,filter_size, depth, depth_category) %>%
  summarise(mean_pp = mean(npp_rate, na.rm = TRUE)) %>%
  pivot_wider(id_cols = cruise:depth_category,
              names_from = filter_size,
              values_from = mean_pp)


#Deselect the Dark sample so that can use just the case, stn, depth and different size classes for the production data
for_int_all <- pp_filt_noreps %>%
  #here you need to input the cruise name through a variable
  ungroup() %>%
  select(cruise, cast, nearest_station, depth, depth_category, `>0`,`>0&<5`,`>0&<20`) #Not including <10 here because these values only exist at the surface so they cannot be integrated according to this method

#Need to add a '0' for depth for each cruise for each station to be able to integrate later
for_int_all_subset <- for_int_all %>%
  filter(depth < 8.2)

#change depth to be 0 and depth category to be D0 and copy the pp values over so that pp is the same at 0 that it is at the first measured value 
for_int_all_subset$depth <- 0
for_int_all_subset$depth_category <- "D0"

#add the subsetted dataframe back into the original for_int_all
for_int_all <- rbind(for_int_all, for_int_all_subset)



########### Start of for loop process for calculating integrated production
cruisename <- for_int_all %>% #START with the cruise that has the longest # of stations
  select(cruise)%>%
  unique()
cruisecounts <- for_int_all %>% #START with the cruise that has the longest # of stations
  count(cruise)

stnlength <- for_int_all %>% #START with the cruise that has the longest # of stations
  select(cruise, nearest_station)%>%
  ungroup() %>%
  unique() %>%
  count(cruise)

ncruise <- length(cruisecounts$n)

#Defining dataframes for data to be placed
PPt_less5 <- data.frame(matrix(NA, ncol = ncruise, nrow = max(stnlength$n)))
colnames(PPt_less5) <- c(cruisename$cruise)
PPt_less_20 <- data.frame(matrix(NA, ncol = ncruise, nrow = max(stnlength$n)))
colnames(PPt_less_20) <- c(cruisename$cruise)
PPt_GFF <- data.frame(matrix(NA, ncol = ncruise, nrow = max(stnlength$n)))
colnames(PPt_GFF) <- c(cruisename$cruise)

#station per cruise dataframe to keep track of them
PPt_stations <- data.frame(matrix(NA, ncol = ncruise, nrow = max(stnlength$n)))
colnames(PPt_stations) <- c(cruisename$cruise)

PPt_stations <- PPt_stations %>%
  rename_all(paste0, "_stn")

#for loop goes here, it selects which cruise you're dealing with from 1 to 3
for (ii in 1:length(cruisename$cruise)){
  
  #here you define a variable, that changes with the for loop, that contains the name of the cruise to be worked on for this for-loop step
  for_int <- for_int_all %>%
    filter(cruise == cruisename$cruise[ii]) %>%
    arrange(nearest_station) #Arrange dataframe so that the stations are in order for the cruise.
  
  
  STATION <- for_int %>% #START with the cruise that has the longest # of stations
    select(nearest_station) %>%
    unique()
  
  
  PPt_stations[1:length(STATION$nearest_station),ii] = arrange(STATION, nearest_station) #Helps keep track of the order of the PP values
  
  n_station <- length(STATION$nearest_station) #Determine the length of STATION, in this case, it's 7 and assign that to "n_station"
  
  
  for (n1 in 1:n_station) { #this is the start of the loop--looping through the list of 1:n_station
    a1 = which(for_int$nearest_station == STATION$nearest_station[n1]); #searching the for_int$station column for instances where that column has values that match those in the STATION variable and assigning to a1, the index
    less_5_a1 = for_int$`>0&<5`[a1] #Selecting the row in the for_int$less_5 column (less_5 production) that matches to the a1 index which signifies each station and assigning this to "less_5_a1"
    less_20_a1 = for_int$`>0&<20`[a1]
    GFF_a1 = for_int$`>0`[a1]
    depth_a1 = for_int$depth[a1] #Selecting the row in the for_int$depth column (depth) that matches to the a1 index assigned above and assigning this to "depth_al"
    tmp=sort(depth_a1,index.return=TRUE); #this is sorting values (depths) in depth_a1 in ascending order
    depth_a2=tmp$x; a3=tmp$ix #and putting those sorted depth values in depth_a2 in ascending order and the respective indices in a3
    less_5_a2 = less_5_a1[a3]; #Rearrange the production data (less_5_a1) according to the new sorted depth indices and assign to less_5_a2
    less_20_a2 = less_20_a1[a3];
    GFF_a2 = GFF_a1[a3];
    
    PP_less5_int = (1:(length(a1)-1))*NA; #create a list from 1 to (length of a1) and fill that in with "NAs" to make it empty and ready
    PP_less_20_int = (1:(length(a1)-1))*NA;
    PP_GFF_int = (1:(length(a1)-1))*NA;
    

    for (n2 in 1:length(a1)-1) { #The start of another loop to calculate the integrated rates. This is looping through 1 to the (length of a1) -1 because there are 1 less depth intervals than the depths themselves
      PP_less5_int[n2] = ((less_5_a2[n2]+less_5_a2[n2+1])*(depth_a2[n2+1] - depth_a2[n2]))/2
      PP_less_20_int[n2] = ((less_20_a2[n2] + less_20_a2[n2+1]) * (depth_a2[n2+1] - depth_a2[n2]))/2
      PP_GFF_int[n2] = ((GFF_a2[n2] + GFF_a2[n2+1]) * (depth_a2[n2+1] - depth_a2[n2]))/2
      
    } #This is calculating the depth integrated production. So it's pulling out the first production value is less_5_a2 and adding that to the next production value. Then multiplying that whole thing by the difference between d2 and d1 if d2 and d1 are the depth intervals.
    PPt_less5[n1,ii] = sum(PP_less5_int) #This is taking the sum of PP_less5_int and storing it in a variable called PPt_less5. We need to take the sum because the PP_less5_int is the integrated value for each depth interval
    PPt_less_20[n1,ii] = sum(PP_less_20_int)
    PPt_GFF[n1,ii] = sum(PP_GFF_int)
    
    
  }
  
  if (ii == 1) {
    int_prod_less5 <- data.frame(matrix(cruisename$cruise[ii],ncol = 1, nrow = length(STATION$nearest_station)))
    colnames(int_prod_less5) <- "cruise"
    int_prod_less5$nearest_station  <- STATION$nearest_station
    int_prod_less5$pp_less5 <- PPt_less5[1:length(STATION$nearest_station),ii]
    
    int_prod_less_20 <- data.frame(matrix(cruisename$cruise[ii],ncol = 1, nrow = length(STATION$nearest_station)))
    colnames(int_prod_less_20) <- "cruise"
    int_prod_less_20$nearest_station  <- STATION$nearest_station
    int_prod_less_20$pp_less_20 <- PPt_less_20[1:length(STATION$nearest_station),ii]
    
    int_prod_GFF <- data.frame(matrix(cruisename$cruise[ii],ncol = 1, nrow = length(STATION$nearest_station)))
    colnames(int_prod_GFF) <- "cruise"
    int_prod_GFF$nearest_station  <- STATION$nearest_station
    int_prod_GFF$pp_GFF <- PPt_GFF[1:length(STATION$nearest_station),ii]
  }
  
  else {
    sacrificial <- data.frame(matrix(cruisename$cruise[ii],ncol = 1, nrow = length(STATION$nearest_station)))
    colnames(sacrificial) <- "cruise"
    sacrificial$nearest_station <- STATION$nearest_station
    sacrificial$pp_less5 <- PPt_less5[1:length(STATION$nearest_station),ii]
    
    int_prod_less5 <- rbind(int_prod_less5, sacrificial)
    
    sacrificial20 <- data.frame(matrix(cruisename$cruise[ii],ncol = 1, nrow = length(STATION$nearest_station)))
    colnames(sacrificial20) <- "cruise"
    sacrificial20$nearest_station  <- STATION$nearest_station
    sacrificial20$pp_less_20 <- PPt_less_20[1:length(STATION$nearest_station),ii]
    
    int_prod_less_20 <- rbind(int_prod_less_20, sacrificial20)
    
    sacrificial_GFF <- data.frame(matrix(cruisename$cruise[ii],ncol = 1, nrow = length(STATION$nearest_station)))
    colnames(sacrificial_GFF) <- "cruise"
    sacrificial_GFF$nearest_station <- STATION$nearest_station
    sacrificial_GFF$pp_GFF <- PPt_GFF[1:length(STATION$nearest_station),ii]
    
    int_prod_GFF <- rbind(int_prod_GFF, sacrificial_GFF)
  }
}


#Now cbind all size dataframes together
pp_int <- int_prod_GFF %>%
  left_join(int_prod_less_20, by = c('cruise', 'nearest_station')) %>%
  left_join(int_prod_less5, by = c('cruise', 'nearest_station')) %>%
  pivot_longer(cols = 'pp_GFF':'pp_less5',
               names_to = 'size',
               values_to = 'npp')

#Manipulate df to put into finalized form
pp_int <- pp_int %>%
  rename(filter_size = size,
         integrated_npp_mg_m2_day = npp)

pp_int$filter_size[pp_int$filter_size == "pp_less5"] <- ">0&<5"
pp_int$filter_size[pp_int$filter_size == "pp_GFF"] <- ">0"
pp_int$filter_size[pp_int$filter_size == "pp_less_20"] <- ">0&<20"

##Add cast information in
# Make stations dataframe to get cast info
stns <-  pp_filt_noreps %>%
  ungroup()%>%
  select(cruise, nearest_station, cast) %>%
  distinct()

pp_int <- pp_int %>%
  inner_join(stns, by = c("cruise","nearest_station"))

#Then need date/time and lat/long so use pp dataframe
lat_lon_date <- pp_for_int %>%
  select(cruise, cast, latitude, longitude, date_time_utc) %>%
  dplyr::group_by(cruise, cast) %>%
  dplyr::summarise(mean_date = mean(date_time_utc), #get mean info by cruise and cast
                  mean_long = mean(longitude),
                  mean_lat = mean(latitude))

#Join with integrated DF
pp_int <-  pp_int %>%
  left_join(lat_lon_date, by = c("cruise", "cast"))

#Rename for final datatable
pp_int <- pp_int %>%
  rename(date_time_utc = mean_date,
         latitude = mean_lat,
        longitude = mean_long)


#Add quality flag
pp_int$iode_quality_flag  <-with(pp_int,
                                 ifelse(integrated_npp_mg_m2_day< 0, "3",
                                                                      ifelse(integrated_npp_mg_m2_day >0, "1",NA))) 

#reorder and write to csv
pp_int <- pp_int %>%
  select(cruise, date_time_utc, latitude, longitude, nearest_station, cast, filter_size, integrated_npp_mg_m2_day, iode_quality_flag)

#add beam attenuation information from output files obtained using automated calculations in MATLBAB. See this Github page for more details: https://github.com/pmarrec/nes-lter-kd-calculation. The "extinction coefficient.csv" file is a file made by compiling output files from here: https://github.com/pmarrec/nes-lter-kd-calculation/tree/main/output_files
ext <- read_csv("extinction_coefficients.csv")
ext$station <- as.character(ext$station)

#join extinc coefficents with integrated table
pp_int <- pp_int %>%
  left_join(ext, by = c("cruise", "cast"))

#Round decimals
pp_int$integrated_npp_mg_m2_day <- round(pp_int$integrated_npp_mg_m2_day, digits = 3)
pp_int$latitude <- round(pp_int$latitude, digits = 4)
pp_int$longitude <- round(pp_int$longitude, digits = 4)


#change T and Z for time/date
pp_int$date_time_utc <- str_replace(pp_int$date_time_utc, "T", " ")

#reorder cruises chronological
pp_int$cruise <- fct_relevel(pp_int$cruise, c("EN644", "AR39B","EN649", "EN655", "EN657", "EN661", "EN668", "AR61B", "AT46", "EN687"))
pp_int <- pp_int %>%
  arrange(cruise)

#selecting for version1
pp_int_vers1 <- pp_int %>%
  filter(cruise %in% c("EN644", "AR39B", "EN649"))

pp_int_vers2 <- pp_int %>%
  filter(cruise %in% c("EN655", "EN657", "EN661", "EN668", "AR61B", "AT46", "EN687"))


#Trying to make the quality flag for these two size classes to be 9

pp_int_vers2_noL11_na <- pp_int_vers2 %>%
  filter(!(cruise == "AT46" & station == "L11" & filter_size == ">0&<20")) %>%
  filter(!(cruise == "AT46" & station == "L11" & filter_size == ">0&<5"))

pp_int_vers2_correctedL11 <- pp_int_vers2 %>%
  filter(cruise == "AT46" & station == "L11" & filter_size %in%c(">0&<20",">0&<5"))

pp_int_vers2_correctedL11$iode_quality_flag <- 9

pp_int_vers2_joined <- rbind(pp_int_vers2_noL11_na,pp_int_vers2_correctedL11)


#write_csv
write_csv(pp_int_vers1 , "For_EDI/npp_integrated_version1.csv")
write_csv(pp_int_vers2_joined, "For_EDI/npp_integrated_version2.csv")





