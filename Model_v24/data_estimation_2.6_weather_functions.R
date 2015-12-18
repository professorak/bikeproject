#get local attributes to stationXstate
get_local_attributes_st_state <- function(wdcMerged, points) {
  #local density attributes of a station, and stockout indicator for n
  #population density, 
  #metro on, _1 in in 0 to 100mts and _2 100,300
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
                           "googleplaces_den_3","googleplaces_den_4","log(googleplaces_den_3+1)","log(googleplaces_den_4+1)",
                           "metro_den_on_a","metro_den_on_b","metro_den_on_c",
                           "metro_den_off_a","metro_den_off_b","metro_den_off_c",
                           "googleplaces_den_a","googleplaces_den_b","googleplaces_den_c",
                           "googleplaces_food_den_1","googleplaces_food_den_2","googleplaces_food_den_3","googleplaces_food_den_4",
                           "googleplaces_food_den_a","googleplaces_food_den_b","googleplaces_food_den_c",
                           "googleplaces_grocery_den_1","googleplaces_grocery_den_2","googleplaces_grocery_den_3","googleplaces_grocery_den_4",
                           "googleplaces_grocery_den_a","googleplaces_grocery_den_b","googleplaces_grocery_den_c",
                           "googleplaces_government_den_1","googleplaces_government_den_2","googleplaces_government_den_3","googleplaces_government_den_4",
                           "googleplaces_government_den_a","googleplaces_government_den_b","googleplaces_government_den_c",
                           "touristlocs_den_1","touristlocs_den_2","touristlocs_den_3","touristlocs_den_4",
                           "touristlocs_den_a","touristlocs_den_b","touristlocs_den_c",
                           "googleplaces_lodging_den_1","googleplaces_lodging_den_2","googleplaces_lodging_den_3","googleplaces_lodging_den_4",
                           "googleplaces_museum_den_1","googleplaces_museum_den_2","googleplaces_museum_den_3","googleplaces_museum_den_4",
                           "googleplaces_movie_theater_den_1","googleplaces_movie_theater_den_2","googleplaces_movie_theater_den_3","googleplaces_movie_theater_den_4")
#  scale_cols <- c(0.00752302550265153, 4.83542598489239e-05, 0.000170904722330068, 2.50288841028638e-05, 7.95216548161101e-05, 0.00554762827586867, 0.0224567865081558, 0.137426485855502, 0.974440229038634, 0.000357492953470158, 0.000700384904760689, 0.000167675845000951, 0.000321261423444887, 4.10920914207535e-05, 0.000133062148061594, 0.000264545064954563, 0.00045220415925925, 2.09579084582854e-05, 6.28608402175063e-05, 0.000125158999915861, 0.000210143010126636, 0.0651212281617819, 0.137023149375346, 0.00346919852388493, 0.00412330376332579, 2.94024499593094, 6.20338857812594, 0.0067161687269596, 0.00733216596292547, 6.77715848359382e-06, 1.91722874965474e-05, 2.08411426680313e-05, 2.91859268962381e-06, 9.32564357736556e-06, 1.05593353165835e-05, 0.0109101284547938, 0.0286030639607752, 0.0460947986285681, 0.0352888128942108, 0.22726340549122, 0.655818203815696, 1.35006584425996, 0.00318496579384185, 0.00775049529457599, 0.0121574071509296, 0.00400790584661252, 0.0183371921054834, 0.0498963284549569, 0.102561024012223, 0.000451214063252186, 0.00115142969105167, 0.00148908987749543, 0.000938280137995079, 0.00352060727716425, 0.00925573164530426, 0.0171573370170951, 0.000145161405923873, 0.000334886635528987, 0.00035315295773006)
    
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
    
    points_sub_a <- points[which(dis_v<0.025),]
    points_sub_b <- points[which(dis_v>=0.025 & dis_v<0.05),]
    points_sub_c <- points[which(dis_v>=0.05 & dis_v<0.075),]

    #density 
    census_density <- mean(points_sub_1$weight[which(points_sub_1$type==1)])
    #metro on, in 0 to 100mts and 100,300
    metro_den_1 <- sum(points_sub_1$weight[which(points_sub_1$type==2)])
    metro_den_2 <- sum(points_sub_2$weight[which(points_sub_2$type==2)])  
    metro_den_3 <- sum(points_sub_3$weight[which(points_sub_3$type==2)])  
    metro_den_4 <- sum(points_sub_4$weight[which(points_sub_4$type==2)])  
    
    metro_den_a <- sum(points_sub_a$weight[which(points_sub_a$type==2)])  
    metro_den_b <- sum(points_sub_b$weight[which(points_sub_b$type==2)])  
    metro_den_c <- sum(points_sub_c$weight[which(points_sub_c$type==2)])  
  
    
    metro_den_on_1 <- (station_data$tw[i]!=0)*metro_den_1
    metro_den_on_2 <- (station_data$tw[i]!=0)*metro_den_2
    metro_den_on_3 <- (station_data$tw[i]!=0)*metro_den_3
    metro_den_on_4 <- (station_data$tw[i]!=0)*metro_den_4
    
    metro_den_on_a <- (station_data$tw[i]!=0)*metro_den_a
    metro_den_on_b <- (station_data$tw[i]!=0)*metro_den_b
    metro_den_on_c <- (station_data$tw[i]!=0)*metro_den_c


    metro_den_off_1 <- (station_data$tw[i]==0)*metro_den_1
    metro_den_off_2 <- (station_data$tw[i]==0)*metro_den_2
    metro_den_off_3 <- (station_data$tw[i]==0)*metro_den_3
    metro_den_off_4 <- (station_data$tw[i]==0)*metro_den_4
    
    metro_den_off_a <- (station_data$tw[i]==0)*metro_den_a
    metro_den_off_b <- (station_data$tw[i]==0)*metro_den_b
    metro_den_off_c <- (station_data$tw[i]==0)*metro_den_c
    
    #bus on, in 0 to 100mts and 100,300
    bus_den_1 <- sum(points_sub_1$weight[which(points_sub_1$type==3)])
    bus_den_2 <- sum(points_sub_2$weight[which(points_sub_2$type==3)])
    bus_den_3 <- sum(points_sub_3$weight[which(points_sub_3$type==3)])
    bus_den_4 <- sum(points_sub_4$weight[which(points_sub_4$type==3)])
    
    bus_den_a <- sum(points_sub_a$weight[which(points_sub_a$type==3)])  
    bus_den_b <- sum(points_sub_b$weight[which(points_sub_b$type==3)])  
    bus_den_c <- sum(points_sub_c$weight[which(points_sub_c$type==3)])  
    
    #touristlocs on, in 0 to 100mts and 100,300
    touristlocs_den_1 <- sum(points_sub_1$weight[which(points_sub_1$type==4)])
    touristlocs_den_2 <- sum(points_sub_2$weight[which(points_sub_2$type==4)])  
    touristlocs_den_3 <- sum(points_sub_3$weight[which(points_sub_3$type==4)])  
    touristlocs_den_4 <- sum(points_sub_4$weight[which(points_sub_4$type==4)])  
    
    touristlocs_den_a <- sum(points_sub_a$weight[which(points_sub_a$type==4)])  
    touristlocs_den_b <- sum(points_sub_b$weight[which(points_sub_b$type==4)])  
    touristlocs_den_c <- sum(points_sub_c$weight[which(points_sub_c$type==4)])  
    
    #google places count, in 0 to 100mts and 100,300
    googleplaces_den_1 <- sum(points_sub_1$places_count)
    googleplaces_den_2 <- sum(points_sub_2$places_count)
    googleplaces_den_3 <- sum(points_sub_3$places_count)
    googleplaces_den_4 <- sum(points_sub_4$places_count)
    
    googleplaces_den_a <- sum(points_sub_a$places_count)  
    googleplaces_den_b <- sum(points_sub_b$places_count)  
    googleplaces_den_c <- sum(points_sub_c$places_count)
    
    
    googleplaces_food_den_1 <- sum(points_sub_1$food)
    googleplaces_food_den_2 <- sum(points_sub_2$food)
    googleplaces_food_den_3 <- sum(points_sub_3$food)
    googleplaces_food_den_4 <- sum(points_sub_4$food)
    
    googleplaces_food_den_a <- sum(points_sub_a$food)  
    googleplaces_food_den_b <- sum(points_sub_b$food)  
    googleplaces_food_den_c <- sum(points_sub_c$food)
    
    googleplaces_grocery_den_1 <- sum(points_sub_1$grocery_or_supermarket)
    googleplaces_grocery_den_2 <- sum(points_sub_2$grocery_or_supermarket)
    googleplaces_grocery_den_3 <- sum(points_sub_3$grocery_or_supermarket)
    googleplaces_grocery_den_4 <- sum(points_sub_4$grocery_or_supermarket)
    
    googleplaces_grocery_den_a <- sum(points_sub_a$grocery_or_supermarket)  
    googleplaces_grocery_den_b <- sum(points_sub_b$grocery_or_supermarket)  
    googleplaces_grocery_den_c <- sum(points_sub_c$grocery_or_supermarket)
    
    
    googleplaces_government_den_1 <- sum(points_sub_1$local_government_office)
    googleplaces_government_den_2 <- sum(points_sub_2$local_government_office)
    googleplaces_government_den_3 <- sum(points_sub_3$local_government_office)
    googleplaces_government_den_4 <- sum(points_sub_4$local_government_office)
    
    googleplaces_government_den_a <- sum(points_sub_a$local_government_office)  
    googleplaces_government_den_b <- sum(points_sub_b$local_government_office)  
    googleplaces_government_den_c <- sum(points_sub_c$local_government_office)
    
    #lodging
    googleplaces_lodging_den_1 <- sum(points_sub_1$lodging)
    googleplaces_lodging_den_2 <- sum(points_sub_2$lodging)
    googleplaces_lodging_den_3 <- sum(points_sub_3$lodging)
    googleplaces_lodging_den_4 <- sum(points_sub_4$lodging)
    
    #museum
    googleplaces_museum_den_1 <- sum(points_sub_1$museum)
    googleplaces_museum_den_2 <- sum(points_sub_2$museum)
    googleplaces_museum_den_3 <- sum(points_sub_3$museum)
    googleplaces_museum_den_4 <- sum(points_sub_4$museum)
    
    #movie_theater  
    googleplaces_movie_theater_den_1 <- sum(points_sub_1$movie_theater)
    googleplaces_movie_theater_den_2 <- sum(points_sub_2$movie_theater)
    googleplaces_movie_theater_den_3 <- sum(points_sub_3$movie_theater)
    googleplaces_movie_theater_den_4 <- sum(points_sub_4$movie_theater)
    
    
    station_data_density_dt[i,] <- c(census_density,metro_den_on_1,metro_den_on_2,metro_den_off_1,metro_den_off_2,
                                     bus_den_1,bus_den_2,googleplaces_den_1,googleplaces_den_2,
                                     metro_den_on_3,metro_den_on_4,metro_den_off_3,metro_den_off_4,
                                     log(metro_den_on_1+1),log(metro_den_on_2+1),log(metro_den_on_3+1),log(metro_den_on_4+1),
                                     log(metro_den_off_1+1),log(metro_den_off_2+1),log(metro_den_off_3+1),log(metro_den_off_4+1),
                                     bus_den_3,bus_den_4,log(bus_den_3+1),log(bus_den_4+1),
                                     googleplaces_den_3,googleplaces_den_4,log(googleplaces_den_3+1),log(googleplaces_den_4+1),
                                     metro_den_on_a,metro_den_on_b,metro_den_on_c,
                                     metro_den_off_a,metro_den_off_b,metro_den_off_c,
                                     googleplaces_den_a,googleplaces_den_b,googleplaces_den_c,
                                     googleplaces_food_den_1, googleplaces_food_den_2, googleplaces_food_den_3, googleplaces_food_den_4,
                                     googleplaces_food_den_a, googleplaces_food_den_b, googleplaces_food_den_c,
                                     googleplaces_grocery_den_1, googleplaces_grocery_den_2, googleplaces_grocery_den_3, googleplaces_grocery_den_4,
                                     googleplaces_grocery_den_a, googleplaces_grocery_den_b, googleplaces_grocery_den_c,
                                     googleplaces_government_den_1, googleplaces_government_den_2, googleplaces_government_den_3, googleplaces_government_den_4,
                                     googleplaces_government_den_a, googleplaces_government_den_b, googleplaces_government_den_c,
                                     touristlocs_den_1,touristlocs_den_2,touristlocs_den_3,touristlocs_den_4,
                                     touristlocs_den_a,touristlocs_den_b,touristlocs_den_c,
                                     googleplaces_lodging_den_1,googleplaces_lodging_den_2,googleplaces_lodging_den_3,googleplaces_lodging_den_4,
                                     googleplaces_museum_den_1,googleplaces_museum_den_2,googleplaces_museum_den_3,googleplaces_museum_den_4,
                                     googleplaces_movie_theater_den_1,googleplaces_movie_theater_den_2,googleplaces_movie_theater_den_3,googleplaces_movie_theater_den_4
                                     )#/scale_cols
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

#   station_data_notw <- unique(wdcMerged[,c("station_id_index","tract","lat","lon")])
#   catchment_area_df <- get_catchment_area_data(station_data_notw)
#   catchment_area_df_full <- catchment_area_df[wdcMerged_density_dt$station_id_index,]
#   wdcMerged_density_dt <- cbind(wdcMerged_density_dt,catchment_area_df_full)  

  return(wdcMerged_density_dt)
}


get_local_attributes_st_state_alt_tw <- function(wdcMerged, points, offtw_list) {
  #local density attributes of a station, and stockout indicator for n
  #population density, 
  #metro on, _1 in in 0 to 100mts and _2 100,300
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
                           "googleplaces_den_3","googleplaces_den_4","log(googleplaces_den_3+1)","log(googleplaces_den_4+1)",
                           "metro_den_on_a","metro_den_on_b","metro_den_on_c",
                           "metro_den_off_a","metro_den_off_b","metro_den_off_c",
                           "googleplaces_den_a","googleplaces_den_b","googleplaces_den_c",
                           "googleplaces_food_den_1","googleplaces_food_den_2","googleplaces_food_den_3","googleplaces_food_den_4",
                           "googleplaces_food_den_a","googleplaces_food_den_b","googleplaces_food_den_c",
                           "googleplaces_grocery_den_1","googleplaces_grocery_den_2","googleplaces_grocery_den_3","googleplaces_grocery_den_4",
                           "googleplaces_grocery_den_a","googleplaces_grocery_den_b","googleplaces_grocery_den_c",
                           "googleplaces_government_den_1","googleplaces_government_den_2","googleplaces_government_den_3","googleplaces_government_den_4",
                           "googleplaces_government_den_a","googleplaces_government_den_b","googleplaces_government_den_c",
                           "touristlocs_den_1","touristlocs_den_2","touristlocs_den_3","touristlocs_den_4",
                           "touristlocs_den_a","touristlocs_den_b","touristlocs_den_c",
                           "googleplaces_lodging_den_1","googleplaces_lodging_den_2","googleplaces_lodging_den_3","googleplaces_lodging_den_4",
                           "googleplaces_museum_den_1","googleplaces_museum_den_2","googleplaces_museum_den_3","googleplaces_museum_den_4",
                           "googleplaces_movie_theater_den_1","googleplaces_movie_theater_den_2","googleplaces_movie_theater_den_3","googleplaces_movie_theater_den_4")
  #  scale_cols <- c(0.00752302550265153, 4.83542598489239e-05, 0.000170904722330068, 2.50288841028638e-05, 7.95216548161101e-05, 0.00554762827586867, 0.0224567865081558, 0.137426485855502, 0.974440229038634, 0.000357492953470158, 0.000700384904760689, 0.000167675845000951, 0.000321261423444887, 4.10920914207535e-05, 0.000133062148061594, 0.000264545064954563, 0.00045220415925925, 2.09579084582854e-05, 6.28608402175063e-05, 0.000125158999915861, 0.000210143010126636, 0.0651212281617819, 0.137023149375346, 0.00346919852388493, 0.00412330376332579, 2.94024499593094, 6.20338857812594, 0.0067161687269596, 0.00733216596292547, 6.77715848359382e-06, 1.91722874965474e-05, 2.08411426680313e-05, 2.91859268962381e-06, 9.32564357736556e-06, 1.05593353165835e-05, 0.0109101284547938, 0.0286030639607752, 0.0460947986285681, 0.0352888128942108, 0.22726340549122, 0.655818203815696, 1.35006584425996, 0.00318496579384185, 0.00775049529457599, 0.0121574071509296, 0.00400790584661252, 0.0183371921054834, 0.0498963284549569, 0.102561024012223, 0.000451214063252186, 0.00115142969105167, 0.00148908987749543, 0.000938280137995079, 0.00352060727716425, 0.00925573164530426, 0.0171573370170951, 0.000145161405923873, 0.000334886635528987, 0.00035315295773006)
  
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
    
    points_sub_a <- points[which(dis_v<0.025),]
    points_sub_b <- points[which(dis_v>=0.025 & dis_v<0.05),]
    points_sub_c <- points[which(dis_v>=0.05 & dis_v<0.075),]
    
    #density 
    census_density <- mean(points_sub_1$weight[which(points_sub_1$type==1)])
    #metro on, in 0 to 100mts and 100,300
    metro_den_1 <- sum(points_sub_1$weight[which(points_sub_1$type==2)])
    metro_den_2 <- sum(points_sub_2$weight[which(points_sub_2$type==2)])  
    metro_den_3 <- sum(points_sub_3$weight[which(points_sub_3$type==2)])  
    metro_den_4 <- sum(points_sub_4$weight[which(points_sub_4$type==2)])  
    
    metro_den_a <- sum(points_sub_a$weight[which(points_sub_a$type==2)])  
    metro_den_b <- sum(points_sub_b$weight[which(points_sub_b$type==2)])  
    metro_den_c <- sum(points_sub_c$weight[which(points_sub_c$type==2)])  
    
    
    metro_den_on_1 <- !(station_data$tw[i] %in% offtw_list)*metro_den_1
    metro_den_on_2 <- !(station_data$tw[i] %in% offtw_list)*metro_den_2
    metro_den_on_3 <- !(station_data$tw[i] %in% offtw_list)*metro_den_3
    metro_den_on_4 <- !(station_data$tw[i] %in% offtw_list)*metro_den_4
    
    metro_den_on_a <- !(station_data$tw[i] %in% offtw_list)*metro_den_a
    metro_den_on_b <- !(station_data$tw[i] %in% offtw_list)*metro_den_b
    metro_den_on_c <- !(station_data$tw[i] %in% offtw_list)*metro_den_c
    
    
    metro_den_off_1 <- (station_data$tw[i] %in% offtw_list)*metro_den_1
    metro_den_off_2 <- (station_data$tw[i] %in% offtw_list)*metro_den_2
    metro_den_off_3 <- (station_data$tw[i] %in% offtw_list)*metro_den_3
    metro_den_off_4 <- (station_data$tw[i] %in% offtw_list)*metro_den_4
    
    metro_den_off_a <- (station_data$tw[i] %in% offtw_list)*metro_den_a
    metro_den_off_b <- (station_data$tw[i] %in% offtw_list)*metro_den_b
    metro_den_off_c <- (station_data$tw[i] %in% offtw_list)*metro_den_c
    
    #bus on, in 0 to 100mts and 100,300
    bus_den_1 <- sum(points_sub_1$weight[which(points_sub_1$type==3)])
    bus_den_2 <- sum(points_sub_2$weight[which(points_sub_2$type==3)])
    bus_den_3 <- sum(points_sub_3$weight[which(points_sub_3$type==3)])
    bus_den_4 <- sum(points_sub_4$weight[which(points_sub_4$type==3)])
    
    bus_den_a <- sum(points_sub_a$weight[which(points_sub_a$type==3)])  
    bus_den_b <- sum(points_sub_b$weight[which(points_sub_b$type==3)])  
    bus_den_c <- sum(points_sub_c$weight[which(points_sub_c$type==3)])  
    
    #touristlocs on, in 0 to 100mts and 100,300
    touristlocs_den_1 <- sum(points_sub_1$weight[which(points_sub_1$type==4)])
    touristlocs_den_2 <- sum(points_sub_2$weight[which(points_sub_2$type==4)])  
    touristlocs_den_3 <- sum(points_sub_3$weight[which(points_sub_3$type==4)])  
    touristlocs_den_4 <- sum(points_sub_4$weight[which(points_sub_4$type==4)])  
    
    touristlocs_den_a <- sum(points_sub_a$weight[which(points_sub_a$type==4)])  
    touristlocs_den_b <- sum(points_sub_b$weight[which(points_sub_b$type==4)])  
    touristlocs_den_c <- sum(points_sub_c$weight[which(points_sub_c$type==4)])  
    
    #google places count, in 0 to 100mts and 100,300
    googleplaces_den_1 <- sum(points_sub_1$places_count)
    googleplaces_den_2 <- sum(points_sub_2$places_count)
    googleplaces_den_3 <- sum(points_sub_3$places_count)
    googleplaces_den_4 <- sum(points_sub_4$places_count)
    
    googleplaces_den_a <- sum(points_sub_a$places_count)  
    googleplaces_den_b <- sum(points_sub_b$places_count)  
    googleplaces_den_c <- sum(points_sub_c$places_count)
    
    
    googleplaces_food_den_1 <- sum(points_sub_1$food)
    googleplaces_food_den_2 <- sum(points_sub_2$food)
    googleplaces_food_den_3 <- sum(points_sub_3$food)
    googleplaces_food_den_4 <- sum(points_sub_4$food)
    
    googleplaces_food_den_a <- sum(points_sub_a$food)  
    googleplaces_food_den_b <- sum(points_sub_b$food)  
    googleplaces_food_den_c <- sum(points_sub_c$food)
    
    googleplaces_grocery_den_1 <- sum(points_sub_1$grocery_or_supermarket)
    googleplaces_grocery_den_2 <- sum(points_sub_2$grocery_or_supermarket)
    googleplaces_grocery_den_3 <- sum(points_sub_3$grocery_or_supermarket)
    googleplaces_grocery_den_4 <- sum(points_sub_4$grocery_or_supermarket)
    
    googleplaces_grocery_den_a <- sum(points_sub_a$grocery_or_supermarket)  
    googleplaces_grocery_den_b <- sum(points_sub_b$grocery_or_supermarket)  
    googleplaces_grocery_den_c <- sum(points_sub_c$grocery_or_supermarket)
    
    
    googleplaces_government_den_1 <- sum(points_sub_1$local_government_office)
    googleplaces_government_den_2 <- sum(points_sub_2$local_government_office)
    googleplaces_government_den_3 <- sum(points_sub_3$local_government_office)
    googleplaces_government_den_4 <- sum(points_sub_4$local_government_office)
    
    googleplaces_government_den_a <- sum(points_sub_a$local_government_office)  
    googleplaces_government_den_b <- sum(points_sub_b$local_government_office)  
    googleplaces_government_den_c <- sum(points_sub_c$local_government_office)
    
    #lodging
    googleplaces_lodging_den_1 <- sum(points_sub_1$lodging)
    googleplaces_lodging_den_2 <- sum(points_sub_2$lodging)
    googleplaces_lodging_den_3 <- sum(points_sub_3$lodging)
    googleplaces_lodging_den_4 <- sum(points_sub_4$lodging)
    
    #museum
    googleplaces_museum_den_1 <- sum(points_sub_1$museum)
    googleplaces_museum_den_2 <- sum(points_sub_2$museum)
    googleplaces_museum_den_3 <- sum(points_sub_3$museum)
    googleplaces_museum_den_4 <- sum(points_sub_4$museum)
    
    #movie_theater  
    googleplaces_movie_theater_den_1 <- sum(points_sub_1$movie_theater)
    googleplaces_movie_theater_den_2 <- sum(points_sub_2$movie_theater)
    googleplaces_movie_theater_den_3 <- sum(points_sub_3$movie_theater)
    googleplaces_movie_theater_den_4 <- sum(points_sub_4$movie_theater)
    
    
    station_data_density_dt[i,] <- c(census_density,metro_den_on_1,metro_den_on_2,metro_den_off_1,metro_den_off_2,
                                     bus_den_1,bus_den_2,googleplaces_den_1,googleplaces_den_2,
                                     metro_den_on_3,metro_den_on_4,metro_den_off_3,metro_den_off_4,
                                     log(metro_den_on_1+1),log(metro_den_on_2+1),log(metro_den_on_3+1),log(metro_den_on_4+1),
                                     log(metro_den_off_1+1),log(metro_den_off_2+1),log(metro_den_off_3+1),log(metro_den_off_4+1),
                                     bus_den_3,bus_den_4,log(bus_den_3+1),log(bus_den_4+1),
                                     googleplaces_den_3,googleplaces_den_4,log(googleplaces_den_3+1),log(googleplaces_den_4+1),
                                     metro_den_on_a,metro_den_on_b,metro_den_on_c,
                                     metro_den_off_a,metro_den_off_b,metro_den_off_c,
                                     googleplaces_den_a,googleplaces_den_b,googleplaces_den_c,
                                     googleplaces_food_den_1, googleplaces_food_den_2, googleplaces_food_den_3, googleplaces_food_den_4,
                                     googleplaces_food_den_a, googleplaces_food_den_b, googleplaces_food_den_c,
                                     googleplaces_grocery_den_1, googleplaces_grocery_den_2, googleplaces_grocery_den_3, googleplaces_grocery_den_4,
                                     googleplaces_grocery_den_a, googleplaces_grocery_den_b, googleplaces_grocery_den_c,
                                     googleplaces_government_den_1, googleplaces_government_den_2, googleplaces_government_den_3, googleplaces_government_den_4,
                                     googleplaces_government_den_a, googleplaces_government_den_b, googleplaces_government_den_c,
                                     touristlocs_den_1,touristlocs_den_2,touristlocs_den_3,touristlocs_den_4,
                                     touristlocs_den_a,touristlocs_den_b,touristlocs_den_c,
                                     googleplaces_lodging_den_1,googleplaces_lodging_den_2,googleplaces_lodging_den_3,googleplaces_lodging_den_4,
                                     googleplaces_museum_den_1,googleplaces_museum_den_2,googleplaces_museum_den_3,googleplaces_museum_den_4,
                                     googleplaces_movie_theater_den_1,googleplaces_movie_theater_den_2,googleplaces_movie_theater_den_3,googleplaces_movie_theater_den_4
    )#/scale_cols
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
  
  #   station_data_notw <- unique(wdcMerged[,c("station_id_index","tract","lat","lon")])
  #   catchment_area_df <- get_catchment_area_data(station_data_notw)
  #   catchment_area_df_full <- catchment_area_df[wdcMerged_density_dt$station_id_index,]
  #   wdcMerged_density_dt <- cbind(wdcMerged_density_dt,catchment_area_df_full)  
  
  return(wdcMerged_density_dt)
}

get_weather_data <- function(filelist_all_weather) {
  #read weather data
  weather_data <- read_weather_data(filelist_all_weather)
  
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
  weather_data <- weather_data[,c("day","time_halfhour_int","Conditions","Temperature_group",
                                  "Humidity_high","Wind.Speed_high")]
  weather_data$Conditions_num <- as.numeric(weather_data$Conditions)
  weather_data$weather_state <- as.factor(apply(weather_data[,c("Conditions_num","Temperature_group",
                                                                "Humidity_high","Wind.Speed_high")],1, FUN=function(x) {paste(x, collapse="_")}))

  return(weather_data)
}

read_weather_data <- function(filelist_all_weather) {
    weather_data <- c()
  weather_data_dir <- paste0(csv_dir,"/../ParisData/Weather/WeatherbaseCSV")
  for( i in 1:length(filelist_all_weather)) {
    filename <- paste0(weather_data_dir,"/",filelist_all_weather[i])
    weather_data_i <- read.csv(filename) 
    weather_data <- rbind(weather_data, weather_data_i)
  }
  
  weather_data$time <- paste(weather_data$Date,weather_data$Local.Time)
  weather_data$time <- as.POSIXct(strptime(weather_data$time, "%Y-%m-%d %I:%M %p", tz="Europe/Paris"))
  weather_data <- subset(weather_data, !is.na(weather_data$time))
  #convert time in  weather_data to day, tw, time_halfhour_int (closest half hour point)
  weather_data$tw <- as.POSIXlt(weather_data$time)$hour
  weather_data$seconds <- as.numeric(weather_data$time) %% 3600
  weather_data$day <- as.POSIXlt(weather_data$time)$yday
    
  weather_data$time_halfhour_int <- weather_data$tw + round(weather_data$seconds/1800)/2
  #remove duplicates
  weather_data <- weather_data[!duplicated(weather_data[,c("day","time_halfhour_int")]),]
  row.names(weather_data) <- NULL
  weather_data$Date <- NULL
  weather_data$Local.Time <- NULL
  weather_data$time_weather <- weather_data$time
  weather_data$time <- NULL
  
  #complete weather data table
  list_days <- unique(weather_data$day)  
  list_time_halfhour_int <- seq(0,23.5,by=0.5)
  weather_data_full <- merge(data.frame(day=list_days, dummy=1),
    data.frame(time_halfhour_int=list_time_halfhour_int, dummy=1), by="dummy")
  weather_data_full$dummy <- NULL
  
  weather_data_full <- merge(weather_data_full, weather_data, by=c("day","time_halfhour_int"), all.x=T)
  #for NA entries for Conditions, Temperature, Humidity, Wind.Speed fill up average of before and after entry
  for (i in which(is.na(weather_data_full$Conditions))) {
    j <- 1
    while(1) {
      if(i+j<=nrow(weather_data_full) & !is.na(weather_data_full$Conditions[i+j])) {
        weather_data_full$Conditions[i] <- weather_data_full$Conditions[i+j]
        break
      } else if(i-j>0 & !is.na(weather_data_full$Conditions[i+j])) {
        weather_data_full$Conditions[i] <- weather_data_full$Conditions[i+j]
        break
      }
      j <- j + 1
    }
  }
  for (i in which(is.na(weather_data_full$Temperature))) {
    j <- 1
    while(1) {
      if(i+j<=nrow(weather_data_full) & !is.na(weather_data_full$Temperature[i+j])) {
        weather_data_full$Temperature[i] <- weather_data_full$Temperature[i+j]
        break
      } else if(i-j>0 & !is.na(weather_data_full$Temperature[i+j])) {
        weather_data_full$Temperature[i] <- weather_data_full$Temperature[i+j]
        break
      }
      j <- j + 1
    }
  }
  for (i in which(is.na(weather_data_full$Humidity))) {
    j <- 1
    while(1) {
      if(i+j<=nrow(weather_data_full) & !is.na(weather_data_full$Humidity[i+j])) {
        weather_data_full$Humidity[i] <- weather_data_full$Humidity[i+j]
        break
      } else if(i-j>0 & !is.na(weather_data_full$Humidity[i+j])) {
        weather_data_full$Humidity[i] <- weather_data_full$Humidity[i+j]
        break
      }
      j <- j + 1
    }
  }
  for (i in which(is.na(weather_data_full$Wind.Speed))) {    
    j <- 1
    while(1) {
      if(i+j<=nrow(weather_data_full) & !is.na(weather_data_full$Wind.Speed[i+j])) {
        weather_data_full$Wind.Speed[i] <- weather_data_full$Wind.Speed[i+j]
        break
      } else if(i-j>0 & !is.na(weather_data_full$Wind.Speed[i+j])) {
        weather_data_full$Wind.Speed[i] <- weather_data_full$Wind.Speed[i+j]
        break
      }
      j <- j + 1
    }
  }
  weather_data <- weather_data_full
  weather_data_full <- NULL
  
  return(weather_data)
}

get_weather_data_detailed <- function(filelist_all_weather) {
  #read weather data
  weather_data <- read_weather_data(filelist_all_weather)
    
  weather_data$Conditions[which(weather_data$Conditions=="")] <- "Clear"
#   #Group Conditions  
#   weather_data$Conditions_org <- weather_data$Conditions
#   Conditions_group_1 <- c("Thunderstorms and Rain","Light Thunderstorm", "Thunderstorm", "Light Thunderstorms and Rain",
#                           "Light Rain Showers","Light Rain","Heavy Rain Showers","Rain", "Heavy Thunderstorms and Rain")
#   Conditions_group_2 <- c("Patches of Fog","Light Fog","Mist", "Fog","Shallow Fog","Heavy Fog","Partial Fog")
#   Conditions_group_3 <- c("Partly Cloudy", "Light Drizzle", "Mostly Cloudy", "Scattered Clouds", "Overcast")
#   weather_data$Conditions <- as.character(weather_data$Conditions)
#   weather_data$Conditions[which(weather_data$Conditions_org %in% Conditions_group_1)] <- "Thunderstorm and Rain"
#   weather_data$Conditions[which(weather_data$Conditions_org %in% Conditions_group_2)] <- "Mist and Fog"
#   weather_data$Conditions[which(weather_data$Conditions_org %in% Conditions_group_3)] <- "Cloudy"
#   weather_data$Conditions <- as.factor(weather_data$Conditions)
  
  #Temperature <10 and >30
  weather_data$Temperature_group <- 0
  weather_data$Temperature_group[which(weather_data$Temperature<10)] <- 1
  weather_data$Temperature_group[which(weather_data$Temperature>=10 &
                                         weather_data$Temperature<=30)] <- 2
  weather_data$Temperature_group[which(weather_data$Temperature>30)] <- 3
  weather_data$Temperature_group  <- as.factor(weather_data$Temperature_group)
  
  weather_data$Humidity_high <- 0
  weather_data$Humidity_high[which(weather_data$Humidity>40 & 
    weather_data$Humidity<=80)] <- 1
  weather_data$Humidity_high[which(weather_data$Humidity>80)] <- 2
  weather_data$Humidity_high  <- as.factor(weather_data$Humidity_high)
  
  weather_data$Wind.Speed_high <- 0
  weather_data$Wind.Speed_high[which(weather_data$Wind.Speed>20
    & weather_data$Wind.Speed<=30)] <- 1
  weather_data$Wind.Speed_high[which(weather_data$Wind.Speed>30)] <- 2
  weather_data$Wind.Speed_high  <- as.factor(weather_data$Wind.Speed_high)
  
  weather_data <- droplevels(weather_data)
  weather_data <- weather_data[,c("day","time_halfhour_int","Conditions","Temperature_group",
                                  "Humidity_high","Wind.Speed_high")]
  weather_data$Conditions_num <- as.numeric(weather_data$Conditions)
  weather_data$weather_state <- as.factor(apply(weather_data[,c("Conditions_num","Temperature_group",
                                                                "Humidity_high","Wind.Speed_high")],1, FUN=function(x) {paste(x, collapse="_")}))
  
  return(weather_data)
}

# get_weather_data_junk <- function(filelist_all_weather) {
#   #read weather data
#   weather_data <- c()
#   weather_data_dir <- paste0(csv_dir,"/../ParisData/Weather/WeatherbaseCSV")
#   for( i in 1:length(filelist_all_weather)) {
#     filename <- paste0(weather_data_dir,"/",filelist_all_weather[i])
#     weather_data_i <- read.csv(filename) 
#     weather_data <- rbind(weather_data, weather_data_i)
#   }
#   
#   weather_data$time <- paste(weather_data$Date,weather_data$Local.Time)
#   weather_data$time <- as.POSIXct(strptime(weather_data$time, "%Y-%m-%d %I:%M %p"))
#   weather_data <- subset(weather_data, !is.na(weather_data$time))
#   #convert time in  weather_data to integer
#   weather_data$time_halfhour_int <- round(as.numeric(weather_data$time)/1800)
#   #remove duplicates
#   weather_data <- weather_data[!duplicated(weather_data$time_halfhour_int),]
#   row.names(weather_data) <- NULL
#   weather_data$Date <- NULL
#   weather_data$Local.Time <- NULL
#   weather_data$time_weather <- weather_data$time
#   weather_data$time <- NULL
#   
#   
#   weather_data$Conditions[which(weather_data$Conditions=="")] <- "Clear"
#   #Group Conditions  
#   weather_data$Conditions_org <- weather_data$Conditions
#   Conditions_group_1 <- c("Thunderstorms and Rain","Light Thunderstorm", "Thunderstorm", "Light Thunderstorms and Rain",
#                           "Light Rain Showers","Light Rain","Heavy Rain Showers","Rain", "Heavy Thunderstorms and Rain")
#   Conditions_group_2 <- c("Patches of Fog","Light Fog","Mist", "Fog","Shallow Fog","Heavy Fog","Partial Fog")
#   Conditions_group_3 <- c("Partly Cloudy", "Light Drizzle", "Mostly Cloudy", "Scattered Clouds", "Overcast")
#   weather_data$Conditions <- as.character(weather_data$Conditions)
#   weather_data$Conditions[which(weather_data$Conditions_org %in% Conditions_group_1)] <- "Thunderstorm and Rain"
#   weather_data$Conditions[which(weather_data$Conditions_org %in% Conditions_group_2)] <- "Mist and Fog"
#   weather_data$Conditions[which(weather_data$Conditions_org %in% Conditions_group_3)] <- "Cloudy"
#   weather_data$Conditions <- as.factor(weather_data$Conditions)
#   
#   #Temperature <10 and >30
#   weather_data$Temperature_group <- 0
#   weather_data$Temperature_group[which(weather_data$Temperature<10)] <- 1
#   weather_data$Temperature_group[which(weather_data$Temperature>=10 &
#                                          weather_data$Temperature<=30)] <- 2
#   weather_data$Temperature_group[which(weather_data$Temperature>30)] <- 3
#   weather_data$Temperature_group  <- as.factor(weather_data$Temperature_group)
#   
#   weather_data$Humidity_high <- 0
#   weather_data$Humidity_high[which(weather_data$Humidity>80)] <- 1
#   weather_data$Humidity_high  <- as.factor(weather_data$Humidity_high)
#   
#   weather_data$Wind.Speed_high <- 0
#   weather_data$Wind.Speed_high[which(weather_data$Wind.Speed>20)] <- 1
#   weather_data$Wind.Speed_high  <- as.factor(weather_data$Wind.Speed_high)
#   
#   weather_data <- droplevels(weather_data)
#   weather_data <- weather_data[,c("time_halfhour_int","Conditions","Temperature_group",
#                                   "Humidity_high","Wind.Speed_high")]
#   weather_data$Conditions_num <- as.numeric(weather_data$Conditions)
#   ############ Make Junk
#   weather_data$Conditions <- as.factor(rep(0,nrow(weather_data)))
#   weather_data$Temperature_group <- as.factor(rep(0,nrow(weather_data)))
#   weather_data$Humidity_high <- as.factor(rep(0,nrow(weather_data)))
#   weather_data$Wind.Speed_high <- as.factor(rep(0,nrow(weather_data)))
#   weather_data$Conditions_num <- 0
#   
#   
#   weather_data$weather_state <- as.factor(apply(weather_data[,c("Conditions_num","Temperature_group",
#                                                                 "Humidity_high","Wind.Speed_high")],1, FUN=function(x) {paste(x, collapse="_")}))
#   
#   return(weather_data)
# }

get_weather_scale <- function(wdcMerged) {
  #generate half hourly, total number of trips, total observations,
  #total stocked out observations and total stocked in observations  
  
  agg_out_dem = binned_sum(wdcMerged$out_dem,wdcMerged$day_time_halfhour_int_fac)  
  agg_in_dem = binned_sum(wdcMerged$in_dem,wdcMerged$day_time_halfhour_int_fac)
  if(!identical(rownames(agg_out_dem),rownames(agg_in_dem))) stop("in aggdata, aggregated data not in same order")
  agg_stockedout = binned_sum(wdcMerged$stocked_out,wdcMerged$day_time_halfhour_int_fac)
  agg_stockedfull = binned_sum(wdcMerged$stocked_full,wdcMerged$day_time_halfhour_int_fac)
  
  wdcMerged_weather_df <- data.frame(day_time_halfhour_int_fac=rownames(agg_out_dem),
                                     out_dem_sum=agg_out_dem[,"sum"],in_dem_sum=agg_in_dem[,"sum"],
                                     total_obs=agg_out_dem[,"count"], stockedout_obs=agg_stockedout[,"sum"],
                                     stockedfull_obs=agg_stockedfull[,"sum"])
  
  wdcMerged_weather_df$non_stockedout_obs <- wdcMerged_weather_df$total_obs - 
    wdcMerged_weather_df$stockedout_obs
  wdcMerged_weather_df$non_stockedfull_obs <- wdcMerged_weather_df$total_obs - 
    wdcMerged_weather_df$stockedfull_obs
  
  wdcMerged_weather_df <- subset(wdcMerged_weather_df, non_stockedout_obs!=0 & 
                                   non_stockedfull_obs!=0)
  
  wdcMerged_weather_df$out_dem_rate <- wdcMerged_weather_df$out_dem_sum/wdcMerged_weather_df$non_stockedout_obs
  wdcMerged_weather_df$in_dem_rate <- wdcMerged_weather_df$in_dem_sum/wdcMerged_weather_df$non_stockedfull_obs
  #wdcMerged_weather_df$time_halfhour_int <- as.numeric(as.character(wdcMerged_weather_df$day_time_halfhour_int_fac))
  
  wdcMerged_weather_df_save <- wdcMerged_weather_df
  
  #merge with weather data to get weather variables.
  weather_data <- get_weather_data_detailed(filelist_all_weather)
  weather_data$day_time_halfhour_int_fac <- as.factor(paste(weather_data$day,
                                                            weather_data$time_halfhour_int))
  nrow_wdcMerged_weather_df <- nrow(wdcMerged_weather_df)
  wdcMerged_weather_df <- merge(wdcMerged_weather_df,weather_data, by="day_time_halfhour_int_fac", all.x=T)    
  if(nrow(wdcMerged_weather_df) < nrow_wdcMerged_weather_df) {stop("In get_weather_scale: ROWS missing in weather data")}
  
  #get tw at 1/2 hour level
  wdcMerged_weather_df$tw <- floor(wdcMerged_weather_df$time_halfhour_int)
  
  #run reg
  fit1 <- lm(log(out_dem_rate+1e-3) ~ factor(tw), data=wdcMerged_weather_df)
  summary(fit1)
  #displace 0 rates
  idx <- which(wdcMerged_weather_df$out_dem_rate==0)
  wdcMerged_weather_df$out_dem_rate[idx] <- min(wdcMerged_weather_df$out_dem_rate[-idx])
  idx <- which(wdcMerged_weather_df$in_dem_rate==0)
  wdcMerged_weather_df$in_dem_rate[idx] <- min(wdcMerged_weather_df$in_dem_rate[-idx])
  
  fit_out_dem <- lm(log(out_dem_rate) ~ factor(tw) + Conditions + Temperature_group +
                      Humidity_high + Wind.Speed_high, data=wdcMerged_weather_df,
                    weights=non_stockedout_obs)
  print("summary(fit_out_dem)")
  print(summary(fit_out_dem))
  
  fit_in_dem <- lm(log(in_dem_rate) ~ factor(tw) + Conditions + Temperature_group +
                  Humidity_high + Wind.Speed_high, data=wdcMerged_weather_df,
                  weights=non_stockedfull_obs)
  print("summary(fit_in_dem)")
  print(summary(fit_in_dem))
  
  #predict the effect of weather variables for each time_halfhour for out and indem
  wdcMerged_weather_predictdata_df <- wdcMerged_weather_df
  wdcMerged_weather_predictdata_df$tw <- 0 #setting tw=0 so that only weather effects are predicted
  
  outdem_weather_scale= predict(fit_out_dem, wdcMerged_weather_predictdata_df) #-coef(fit)["(Intercept)"]
  outdem_weather_scale <- outdem_weather_scale - mean(outdem_weather_scale, na.rm=T)
  outdem_weather_scale <- exp(outdem_weather_scale)
  outdem_weather_scale[which(is.na(outdem_weather_scale))] <- 1
  if(any(is.na(outdem_weather_scale))) stop("NA's in outdem_weather_scale")
  
  indem_weather_scale= predict(fit_in_dem, wdcMerged_weather_predictdata_df) #-coef(fit)["(Intercept)"]
  indem_weather_scale <- indem_weather_scale - mean(indem_weather_scale, na.rm=T)
  indem_weather_scale <- exp(indem_weather_scale)
  indem_weather_scale[which(is.na(indem_weather_scale))] <- 1
  if(any(is.na(indem_weather_scale))) stop("NA's in indem_weather_scale")
  
  weather_scale_df <- data.frame(day_time_halfhour_int_fac=wdcMerged_weather_predictdata_df$day_time_halfhour_int_fac,         
         outdem_weather_scale=outdem_weather_scale, indem_weather_scale=indem_weather_scale)
  
  return(weather_scale_df)
}

get_catchment_area_data <- function(station_data) {  
  dis_points_org <- dis_points
  dis_points <<- 0.010  
  file <- paste0("catchment_area_df_dis",paste0(tract,collapse="_"),"_dis_points_",dis_points)
  filepath <- paste0(csv_dir,"/Paris/",file,".RData")
  ######
  if(!file.exists(filepath)) {
    catchment_area_df <- get_catchment_area_data_in(station_data)    
    save(catchment_area_df,file=filepath)
  }
  ###
  load(filepath)  
  dis_points <<- dis_points_org
  return(catchment_area_df)
}

get_catchment_area_data_in <- function(station_data) {  
  
  points_dense <- generate_integration_points(station_data)
  points_dense <- subset(points_dense, type==1)
  #find closest station to each point
  closest_stid <- rep(0,nrow(points_dense))
  closest_stid_distance <- rep(0,nrow(points_dense))
  for(i in 1:nrow(points_dense)) {    
    lat1 = points_dense$lat[i]
    lon1 = points_dense$lon[i]
    dis_v <- latlondistance(lat1,lon1,station_data$lat,station_data$lon)  
    closest_stid[i] <- order(dis_v)[1]
    closest_stid_distance[i] <-  dis_v[order(dis_v)[1]]          
  }
  closest_stid_df <- data.frame(closest_stid, closest_stid_distance)
  catchment_area_df <- data.frame(catchment_area=rep(0,nrow(station_data)))
  #catchment area, for each station count number of points
  tab_catchment_area <- table(closest_stid_df$closest_stid)
  catchment_area_df$catchment_area[as.numeric(rownames(tab_catchment_area))] <- tab_catchment_area
  #catchment area within 0.076810 << 0.116700 << 0.165900;
  tab_catchment_area_step1 <- 
    table(closest_stid_df$closest_stid[which(closest_stid_df$closest_stid_distance<=0.076810)])
  catchment_area_df$catchment_area_step1[as.numeric(rownames(tab_catchment_area_step1))] <- 
    tab_catchment_area_step1
  tab_catchment_area_step2 <- 
    table(closest_stid_df$closest_stid[which(closest_stid_df$closest_stid_distance<=0.116700 &
                                               closest_stid_df$closest_stid_distance>0.076810)])
  tab_catchment_area_step3 <- 
    table(closest_stid_df$closest_stid[which(closest_stid_df$closest_stid_distance<=0.165900 &
                                               closest_stid_df$closest_stid_distance>0.116700)])
  catchment_area_df$catchment_area_step2 <- 0
  catchment_area_df$catchment_area_step2[as.numeric(rownames(tab_catchment_area_step2))] <- 
    tab_catchment_area_step2
  catchment_area_df$catchment_area_step3 <- 0
  catchment_area_df$catchment_area_step3[as.numeric(rownames(tab_catchment_area_step3))] <- 
    tab_catchment_area_step3
  tab_catchment_area_step4 <- 
    table(closest_stid_df$closest_stid[which(closest_stid_df$closest_stid_distance>0.165900)])
  catchment_area_df$catchment_area_step4 <- 0
  catchment_area_df$catchment_area_step4[as.numeric(rownames(tab_catchment_area_step4))] <- 
    tab_catchment_area_step4

  return(catchment_area_df)
}


