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
                "bus_den_1","bus_den_2","googleplaces_den_1","googleplaces_den_2")
  station_data_density_dt <- as.data.frame(matrix(NA,nrow=nrow(station_data), ncol=length(colnames_density_dt)))
  colnames(station_data_density_dt) <- colnames_density_dt
  
  for(i in 1:nrow(station_data)) {
    lat1 = station_data$lat[i]
    lon1 = station_data$lon[i]
    dis_v <- latlondistance(lat1,lon1,points$lat,points$lon)  
    points_sub_1 <- points[which(dis_v<0.1),]
    points_sub_2 <- points[which(dis_v>=0.1 & dis_v<=0.3),]
    #density 
    census_density <- mean(points_sub_1$weight[which(points_sub_1$type==1)])
    #metro on, in 0 to 100mts and 100,300
    metro_den_1 <- sum(points_sub_1$weight[which(points_sub_1$type==2)])
    metro_den_2 <- sum(points_sub_2$weight[which(points_sub_2$type==2)])  
    metro_den_on_1 <- (station_data$tw[i]!=0)*metro_den_1
    metro_den_on_2 <- (station_data$tw[i]!=0)*metro_den_2
    metro_den_off_1 <- (station_data$tw[i]==0)*metro_den_1
    metro_den_off_2 <- (station_data$tw[i]==0)*metro_den_2

    #bus on, in 0 to 100mts and 100,300
    bus_den_1 <- sum(points_sub_1$weight[which(points_sub_1$type==3)])
    bus_den_2 <- sum(points_sub_2$weight[which(points_sub_2$type==3)])  
#     bus_den_on_1 <- (station_data$tw[i]!=0)*bus_den_1
#     bus_den_on_2 <- (station_data$tw[i]!=0)*bus_den_2
#     bus_den_off_1 <- (station_data$tw[i]==0)*bus_den_1
#     bus_den_off_2 <- (station_data$tw[i]==0)*bus_den_2
    
    #google places count, in 0 to 100mts and 100,300
    googleplaces_den_1 <- sum(points_sub_1$places_count)
    googleplaces_den_2 <- sum(points_sub_2$places_count)  
    station_data_density_dt[i,] <- c(census_density,metro_den_on_1,metro_den_on_2,metro_den_off_1,metro_den_off_2,
                                    bus_den_1,bus_den_2,googleplaces_den_1,googleplaces_den_2)
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

