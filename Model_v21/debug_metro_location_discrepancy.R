#find discrepancy between paris_subway_data.csv and metro_traffic_locations.csv.


metrotrafficloc_file_name <- paste0(dropbox_dir,"/../VelibData/ParisData/ParisTransportData/metro_traffic_locations.csv") 
metro_traffic_locations <- read.csv(metrotrafficloc_file_name)

paris_subway_data_file_name <- paste0(dropbox_dir,"/../VelibData/ParisData/ParisMetroLocationsData/paris_subway_data.csv") 
paris_subway_data <- read.csv(paris_subway_data_file_name)

#for all in metro_traffic_locations, find closest station in paris_subway_data
#and see the stats
metro_traffic_locations_disvec <- rep(NA,nrow(metro_traffic_locations))
for(i in c(1:nrow(metro_traffic_locations))) {
  dis_v <- latlondistance(metro_traffic_locations$lat[i],metro_traffic_locations$lon[i],
                          paris_subway_data$Lat,paris_subway_data$Lon)
  metro_traffic_locations_disvec[i] <- min(dis_v)
}

########################################################################
########################################################################
#check with new metro_stations_data_curated.csv
metro_stations_data_file_name <- paste0(dropbox_dir,"/../VelibData/ParisData/ParisMetroData/metro_stations_data_curated.csv") 
metro_stations_data <- read.csv(metro_stations_data_file_name)
paris_subway_data_file_name <- paste0(dropbox_dir,"/../VelibData/ParisData/Old_ParisMetroLocationsData/paris_subway_data.csv") 
paris_subway_data <- read.csv(paris_subway_data_file_name)

metro_stations_data_disvec <- rep(NA,nrow(metro_stations_data))
for(i in c(1:nrow(metro_stations_data))) {
  dis_v <- latlondistance(metro_stations_data$lat[i],metro_stations_data$lon[i],
                          paris_subway_data$Lat,paris_subway_data$Lon)
  metro_stations_data_disvec[i] <- min(dis_v)
}

summary(metro_stations_data_disvec)



