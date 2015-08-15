# #get local attributes to stationXstate
# get_local_attributes_st_state <- function(wdcMerged, points) {
#   #local density attributes of a station, and stockout indicator for n
#   #population density, 
#   #metro on, in 0 to 100mts and 100,300
#   #metro evening on, in 0 to 100mts and 100,300
#   #google places count, in 0 to 100mts and 100,300
#   #bus stop count, in 0 to 100mts and 100,300
#   station_data <- unique(wdcMerged[,c("station_id_index","tract","tw","lat","lon")])
#   colnames_density_dt <- c("census_density","metro_den_on_1","metro_den_on_2","metro_den_off_1","metro_den_off_2",
#                 "bus_den_1","bus_den_2","googleplaces_den_1","googleplaces_den_2")
#   station_data_density_dt <- as.data.frame(matrix(NA,nrow=nrow(station_data), ncol=length(colnames_density_dt)))
#   colnames(station_data_density_dt) <- colnames_density_dt
#   
#   for(i in 1:nrow(station_data)) {
#     lat1 = station_data$lat[i]
#     lon1 = station_data$lon[i]
#     dis_v <- latlondistance(lat1,lon1,points$lat,points$lon)  
#     points_sub_1 <- points[which(dis_v<0.1),]
#     points_sub_2 <- points[which(dis_v>=0.1 & dis_v<=0.3),]
#     #density 
#     census_density <- mean(points_sub_1$weight[which(points_sub_1$type==1)])
#     #metro on, in 0 to 100mts and 100,300
#     metro_den_1 <- sum(points_sub_1$weight[which(points_sub_1$type==2)])
#     metro_den_2 <- sum(points_sub_2$weight[which(points_sub_2$type==2)])  
#     metro_den_on_1 <- (station_data$tw[i]!=0)*metro_den_1
#     metro_den_on_2 <- (station_data$tw[i]!=0)*metro_den_2
#     metro_den_off_1 <- (station_data$tw[i]==0)*metro_den_1
#     metro_den_off_2 <- (station_data$tw[i]==0)*metro_den_2
# 
#     #bus on, in 0 to 100mts and 100,300
#     bus_den_1 <- sum(points_sub_1$weight[which(points_sub_1$type==3)])
#     bus_den_2 <- sum(points_sub_2$weight[which(points_sub_2$type==3)])  
# #     bus_den_on_1 <- (station_data$tw[i]!=0)*bus_den_1
# #     bus_den_on_2 <- (station_data$tw[i]!=0)*bus_den_2
# #     bus_den_off_1 <- (station_data$tw[i]==0)*bus_den_1
# #     bus_den_off_2 <- (station_data$tw[i]==0)*bus_den_2
#     
#     #google places count, in 0 to 100mts and 100,300
#     googleplaces_den_1 <- sum(points_sub_1$places_count)
#     googleplaces_den_2 <- sum(points_sub_2$places_count)  
#     station_data_density_dt[i,] <- c(census_density,metro_den_on_1,metro_den_on_2,metro_den_off_1,metro_den_off_2,
#                                     bus_den_1,bus_den_2,googleplaces_den_1,googleplaces_den_2)
#   }
#   station_data_density_dt <- cbind(station_data[,c("station_id_index","tw")],station_data_density_dt)
#   #expand station_data_density_dt to wdcMerged
# #  wdcMerged_density_dt <- wdcMerged[,c("station_id_index","tw")]
#   wdcMerged_density_dt <- wdcMerged
#   wdcMerged_density_dt <- merge(wdcMerged_density_dt, station_data_density_dt, by=c("station_id_index","tw"), all.x=T, sort=F)
#   print(identical(wdcMerged_density_dt$station_id_index, wdcMerged$station_id_index))
#   print(identical(wdcMerged_density_dt$tw, wdcMerged$tw))
#   #stockout state indicator (not necessary to exclude focal station stocked_out case as these entries are excluded anyways from computation of GMM objective)
#   wdcMerged_density_dt$sto_nearby <- as.numeric(sapply(wdcMerged$sto_state_local,FUN=splitchar_sum)>0)
#   return(wdcMerged_density_dt)
# }
# 

#get local attributes to stationXstate
get_local_attributes_st_state <- function(wdcMerged, points) {
  #local density attributes of a station, and stockout indicator for n
  #population density, 
  #metro on, in 0 to 100mts and 100,300
  #metro evening on, in 0 to 100mts and 100,300
  #google places count, in 0 to 100mts and 100,300
  #bus stop count, in 0 to 100mts and 100,300
  station_data <- unique(wdcMerged[,c("station_id_index","tract","tw","lat","lon")])
  colnames_density_dt <- c("census_density","metro_den_on_1","metro_den_on_2","metro_den_off_1","metro_den_off_2",
                           "bus_den_1","bus_den_2","googleplaces_den_1","googleplaces_den_2",
                           "metro_den_on_3","metro_den_on_4","metro_den_off_3","metro_den_off_4",
                           "log(metro_den_on_1+1)","log(metro_den_on_2+1)","log(metro_den_on_3+1)","log(metro_den_on_4+1)",
                           "log(metro_den_off_1+1)","log(metro_den_off_2+1)","log(metro_den_off_3+1)","log(metro_den_off_4+1)",
                           "bus_den_3","bus_den_4","log(bus_den_3+1)","log(bus_den_4+1)",
                           "googleplaces_den_3","googleplaces_den_4","log(googleplaces_den_3+1)","log(googleplaces_den_4+1)"
  )
  station_data_density_dt <- as.data.frame(matrix(NA,nrow=nrow(station_data), ncol=length(colnames_density_dt)))
  colnames(station_data_density_dt) <- colnames_density_dt
  
  for(i in 1:nrow(station_data)) {
    lat1 = station_data$lat[i]
    lon1 = station_data$lon[i]
    dis_v <- latlondistance(lat1,lon1,points$lat,points$lon)  
    points_sub_1 <- points[which(dis_v<0.1),]
    points_sub_2 <- points[which(dis_v>=0.1 & dis_v<0.3),]
    points_sub_3 <- points[which(dis_v>=0.3 & dis_v<0.6),]
    points_sub_4 <- points[which(dis_v>=0.6 & dis_v<1.0),]
    #density 
    census_density <- mean(points_sub_1$weight[which(points_sub_1$type==1)])
    #metro on, in 0 to 100mts and 100,300
    metro_den_1 <- sum(points_sub_1$weight[which(points_sub_1$type==2)])
    metro_den_2 <- sum(points_sub_2$weight[which(points_sub_2$type==2)])  
    metro_den_3 <- sum(points_sub_3$weight[which(points_sub_3$type==2)])  
    metro_den_4 <- sum(points_sub_4$weight[which(points_sub_4$type==2)])  
    metro_den_on_1 <- (station_data$tw[i]!=0)*metro_den_1
    metro_den_on_2 <- (station_data$tw[i]!=0)*metro_den_2
    metro_den_on_3 <- (station_data$tw[i]!=0)*metro_den_3
    metro_den_on_4 <- (station_data$tw[i]!=0)*metro_den_4
    metro_den_off_1 <- (station_data$tw[i]==0)*metro_den_1
    metro_den_off_2 <- (station_data$tw[i]==0)*metro_den_2
    metro_den_off_3 <- (station_data$tw[i]==0)*metro_den_3
    metro_den_off_4 <- (station_data$tw[i]==0)*metro_den_4
    
    #bus on, in 0 to 100mts and 100,300
    bus_den_1 <- sum(points_sub_1$weight[which(points_sub_1$type==3)])
    bus_den_2 <- sum(points_sub_2$weight[which(points_sub_2$type==3)])
    bus_den_3 <- sum(points_sub_3$weight[which(points_sub_3$type==3)])
    bus_den_4 <- sum(points_sub_4$weight[which(points_sub_4$type==3)])
    #     bus_den_on_1 <- (station_data$tw[i]!=0)*bus_den_1
    #     bus_den_on_2 <- (station_data$tw[i]!=0)*bus_den_2
    #     bus_den_off_1 <- (station_data$tw[i]==0)*bus_den_1
    #     bus_den_off_2 <- (station_data$tw[i]==0)*bus_den_2
    
    #google places count, in 0 to 100mts and 100,300
    googleplaces_den_1 <- sum(points_sub_1$places_count)
    googleplaces_den_2 <- sum(points_sub_2$places_count)
    googleplaces_den_3 <- sum(points_sub_3$places_count)
    googleplaces_den_4 <- sum(points_sub_4$places_count)
    station_data_density_dt[i,] <- c(census_density,metro_den_on_1,metro_den_on_2,metro_den_off_1,metro_den_off_2,
                                     bus_den_1,bus_den_2,googleplaces_den_1,googleplaces_den_2,
                                     metro_den_on_3,metro_den_on_4,metro_den_off_3,metro_den_off_4,
                                     log(metro_den_on_1+1),log(metro_den_on_2+1),log(metro_den_on_3+1),log(metro_den_on_4+1),
                                     log(metro_den_off_1+1),log(metro_den_off_2+1),log(metro_den_off_3+1),log(metro_den_off_4+1),
                                     bus_den_3,bus_den_4,log(bus_den_3+1),log(bus_den_4+1),
                                     googleplaces_den_3,googleplaces_den_4,log(googleplaces_den_3+1),log(googleplaces_den_4+1)
                                     
    )
  }
  station_data_density_dt <- cbind(station_data[,c("station_id_index","tw")],station_data_density_dt)
  #expand station_data_density_dt to wdcMerged
  #  wdcMerged_density_dt <- wdcMerged[,c("station_id_index","tw")]
  wdcMerged_density_dt <- wdcMerged
  wdcMerged_density_dt <- merge(wdcMerged_density_dt, station_data_density_dt, by=c("station_id_index","tw"), all.x=T, sort=F)
  print(identical(wdcMerged_density_dt$station_id_index, wdcMerged$station_id_index))
  print(identical(wdcMerged_density_dt$tw, wdcMerged$tw))
  #stockout state indicator (not necessary to exclude focal station stocked_out case as these entries are excluded anyways from computation of GMM objective)
  wdcMerged_density_dt$sto_nearby <- as.numeric(sapply(wdcMerged$sto_state_local,FUN=splitchar_sum)>0)
  return(wdcMerged_density_dt)
}

get_weather_data <- function(filelist_all_weather) {
  #read weather data
  weather_data <- c()
  weather_data_dir <- paste0(csv_dir,"/../ParisData/Weather/WeatherbaseCSV")
  for( i in 1:length(filelist_all_weather)) {
    filename <- paste0(weather_data_dir,"/",filelist_all_weather[i])
    weather_data_i <- read.csv(filename) 
    weather_data <- rbind(weather_data, weather_data_i)
  }
  
  weather_data$time <- paste(weather_data$Date,weather_data$Local.Time)
  weather_data$time <- as.POSIXct(strptime(weather_data$time, "%Y-%m-%d %I:%M %p"))
  weather_data <- subset(weather_data, !is.na(weather_data$time))
  #convert time in  weather_data to integer
  weather_data$time_halfhour_int <- round(as.numeric(weather_data$time)/1800)
  #remove duplicates
  weather_data <- weather_data[!duplicated(weather_data$time),]
  row.names(weather_data) <- NULL
  weather_data$Date <- NULL
  weather_data$Local.Time <- NULL
  weather_data$time_weather <- weather_data$time
  weather_data$time <- NULL
  
  
  weather_data$Conditions[which(weather_data$Conditions=="")] <- "Clear"
  #Group Conditions  
  weather_data$Conditions_org <- weather_data$Conditions
  Conditions_group_1 <- c("Thunderstorms and Rain","Light Thunderstorm", "Thunderstorm", "Light Thunderstorms and Rain",
                          "Light Rain Showers","Light Rain","Heavy Rain Showers","Rain", "Heavy Thunderstorms and Rain")
  Conditions_group_2 <- c("Patches of Fog","Light Fog","Mist", "Fog","Shallow Fog","Heavy Fog","Partial Fog")
  Conditions_group_3 <- c("Partly Cloudy", "Light Drizzle", "Mostly Cloudy", "Scattered Clouds", "Overcast")
  weather_data$Conditions <- as.character(weather_data$Conditions)
  weather_data$Conditions[which(weather_data$Conditions_org %in% Conditions_group_1)] <- "Thunderstorm and Rain"
  weather_data$Conditions[which(weather_data$Conditions_org %in% Conditions_group_2)] <- "Mist and Fog"
  weather_data$Conditions[which(weather_data$Conditions_org %in% Conditions_group_3)] <- "Cloudy"
  weather_data$Conditions <- as.factor(weather_data$Conditions)
  
  #Temperature <10 and >30
  weather_data$Temperature_group <- 0
  weather_data$Temperature_group[which(weather_data$Temperature<10)] <- 1
  weather_data$Temperature_group[which(weather_data$Temperature>=10 &
                                         weather_data$Temperature<=30)] <- 2
  weather_data$Temperature_group[which(weather_data$Temperature>30)] <- 3
  weather_data$Temperature_group  <- as.factor(weather_data$Temperature_group)
  
  weather_data$Humidity_high <- 0
  weather_data$Humidity_high[which(weather_data$Humidity>80)] <- 1
  weather_data$Humidity_high  <- as.factor(weather_data$Humidity_high)
  
  weather_data$Wind.Speed_high <- 0
  weather_data$Wind.Speed_high[which(weather_data$Wind.Speed>20)] <- 1
  weather_data$Wind.Speed_high  <- as.factor(weather_data$Wind.Speed_high)
  
  weather_data <- droplevels(weather_data)
  weather_data <- weather_data[,c("time_halfhour_int","Conditions","Temperature_group",
                                  "Humidity_high","Wind.Speed_high")]
  weather_data$Conditions_num <- as.numeric(weather_data$Conditions)
  weather_data$weather_state <- as.factor(apply(weather_data[,c("Conditions_num","Temperature_group",
                                                                "Humidity_high","Wind.Speed_high")],1, FUN=function(x) {paste(x, collapse="_")}))

  return(weather_data)
}

get_weather_data_junk <- function(filelist_all_weather) {
  #read weather data
  weather_data <- c()
  weather_data_dir <- paste0(csv_dir,"/../ParisData/Weather/WeatherbaseCSV")
  for( i in 1:length(filelist_all_weather)) {
    filename <- paste0(weather_data_dir,"/",filelist_all_weather[i])
    weather_data_i <- read.csv(filename) 
    weather_data <- rbind(weather_data, weather_data_i)
  }
  
  weather_data$time <- paste(weather_data$Date,weather_data$Local.Time)
  weather_data$time <- as.POSIXct(strptime(weather_data$time, "%Y-%m-%d %I:%M %p"))
  weather_data <- subset(weather_data, !is.na(weather_data$time))
  #convert time in  weather_data to integer
  weather_data$time_halfhour_int <- round(as.numeric(weather_data$time)/1800)
  #remove duplicates
  weather_data <- weather_data[!duplicated(weather_data$time),]
  row.names(weather_data) <- NULL
  weather_data$Date <- NULL
  weather_data$Local.Time <- NULL
  weather_data$time_weather <- weather_data$time
  weather_data$time <- NULL
  
  
  weather_data$Conditions[which(weather_data$Conditions=="")] <- "Clear"
  #Group Conditions  
  weather_data$Conditions_org <- weather_data$Conditions
  Conditions_group_1 <- c("Thunderstorms and Rain","Light Thunderstorm", "Thunderstorm", "Light Thunderstorms and Rain",
                          "Light Rain Showers","Light Rain","Heavy Rain Showers","Rain", "Heavy Thunderstorms and Rain")
  Conditions_group_2 <- c("Patches of Fog","Light Fog","Mist", "Fog","Shallow Fog","Heavy Fog","Partial Fog")
  Conditions_group_3 <- c("Partly Cloudy", "Light Drizzle", "Mostly Cloudy", "Scattered Clouds", "Overcast")
  weather_data$Conditions <- as.character(weather_data$Conditions)
  weather_data$Conditions[which(weather_data$Conditions_org %in% Conditions_group_1)] <- "Thunderstorm and Rain"
  weather_data$Conditions[which(weather_data$Conditions_org %in% Conditions_group_2)] <- "Mist and Fog"
  weather_data$Conditions[which(weather_data$Conditions_org %in% Conditions_group_3)] <- "Cloudy"
  weather_data$Conditions <- as.factor(weather_data$Conditions)
  
  #Temperature <10 and >30
  weather_data$Temperature_group <- 0
  weather_data$Temperature_group[which(weather_data$Temperature<10)] <- 1
  weather_data$Temperature_group[which(weather_data$Temperature>=10 &
                                         weather_data$Temperature<=30)] <- 2
  weather_data$Temperature_group[which(weather_data$Temperature>30)] <- 3
  weather_data$Temperature_group  <- as.factor(weather_data$Temperature_group)
  
  weather_data$Humidity_high <- 0
  weather_data$Humidity_high[which(weather_data$Humidity>80)] <- 1
  weather_data$Humidity_high  <- as.factor(weather_data$Humidity_high)
  
  weather_data$Wind.Speed_high <- 0
  weather_data$Wind.Speed_high[which(weather_data$Wind.Speed>20)] <- 1
  weather_data$Wind.Speed_high  <- as.factor(weather_data$Wind.Speed_high)
  
  weather_data <- droplevels(weather_data)
  weather_data <- weather_data[,c("time_halfhour_int","Conditions","Temperature_group",
                                  "Humidity_high","Wind.Speed_high")]
  weather_data$Conditions_num <- as.numeric(weather_data$Conditions)
  ############ Make Junk
  weather_data$Conditions <- as.factor(rep(0,nrow(weather_data)))
  weather_data$Temperature_group <- as.factor(rep(0,nrow(weather_data)))
  weather_data$Humidity_high <- as.factor(rep(0,nrow(weather_data)))
  weather_data$Wind.Speed_high <- as.factor(rep(0,nrow(weather_data)))
  weather_data$Conditions_num <- 0
  
  
  weather_data$weather_state <- as.factor(apply(weather_data[,c("Conditions_num","Temperature_group",
                                                                "Humidity_high","Wind.Speed_high")],1, FUN=function(x) {paste(x, collapse="_")}))
  
  return(weather_data)
}
