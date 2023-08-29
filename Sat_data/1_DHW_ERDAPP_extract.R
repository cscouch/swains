#DailyDHW  data extraction using CRW_dhw_v1_0 ERDAPP file and 10km buffers around each island for an 8 year period prior to the NCRMP mission. Runs region by region, then merge all the region files.

library(sf)
library(maps) 
library(mapdata) 
library(maptools) 
library(rerddap) 
library(rerddapXtracto)
library(dplyr)
library(lubridate)

rm(list=ls())

####Load the Island Boundary####
setwd("C:/github/swains/ArcGIS")


isl=st_read('SWA_10k_buff.shp') #can change this to be different buffer polygons 

names <- unique(isl$ISLAND) ## this shapefile has an attribute table that has the name of each island

years <- seq(from = 2008, to = 2010) # sequence of years i want to download data for

ERDDAP_Node="https://oceanwatch.pifsc.noaa.gov/erddap/"

out_df2 <- NULL # define dataframe that the for loops will populate




# nested for loops that cycle through each island and year extracts dhw data (CRW_dhw_v1_0), then merges it all together

for(i in c(1:length(names))){
  for(j in c(1:length(years))){
    isl_2 <- isl[ which(isl$ISLAND == names[i]),] # subset to island as specified by names[i]
    isl_2=st_coordinates(isl_2) # retrieves coordinates in matrix form for that island that we want data for
    
    poly <- as.data.frame(isl_2[,1:2]) # turn into a data frame
    names(poly) <- c("lon","lat") # add column names to coordinates
    
    # convert the longitudes to 0-360? to make it easier to work across the dateline
    I=which(poly$lon<0)
    poly$lon[I]=poly$lon[I]+360 
    
    ## extraction begins ##
    xcoord <- poly$lon # set boundaries
    ycoord <- poly$lat 
    tcoord <- c(paste(years[j],"-01-01", sep = ""), paste(years[j],"-12-31", sep = "")) # set the time period we pull data for
    
    # set data info for sst_cum
    dataInfo <- rerddap::info('CRW_dhw_v1_0', url=ERDDAP_Node) 
    dataInfo$variable$variable_name # choose from these variables
    parameter=dataInfo$variable$variable_name[1] # select sst as parameter of interest
    dhw <- rxtractogon(dataInfo, parameter=parameter, xcoord=xcoord, ycoord=ycoord, tcoord=tcoord)#extract dhw
    
    # calculate the spatial mean and max per island polygon:
    mean_dhw <- apply(dhw$degree_heating_week, 3, mean, na.rm = TRUE)# calculate spatial average -- the '3' refers to the date 
    max_dhw <- apply(dhw$degree_heating_week, 3, max, na.rm = TRUE)
    
    # return a dataframe with the spatial mean per timestep:
    this_df <- data.frame(date = dhw$time, isl = names[i], mean_dhw, max_dhw, stringsAsFactors = FALSE)  
    out_df2 <- rbind(out_df2, this_df)
    
    print(paste("DHW data added for ", names[i])) # print out status (it takes a while to download all this data!)
  }
}


#filter end of dataset to cruise date, if needed
out_df2$date <-as.Date(out_df2$date,"%Y/%m/%d")

#save data
setwd("C:/github/swains/Sat_data")
write.csv(out_df2, "DHW_daily_raw_extended_range.csv")





