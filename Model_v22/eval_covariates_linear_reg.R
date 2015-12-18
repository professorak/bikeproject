eval_covariates_linear_reg <- function(wdcMerged) {
  
  stocked_list <- which(wdcMerged$stocked_out==FALSE)
  #compute attributes matrix X and instruments matrix Z
  #list of X attributes
  #1s, sevice level, diXw fixed effects, month fixed effects, weather fixed effects 
  #list of Z attributes
  #1s, diXw fixed effects, month fixed effects, weather fixed effects, two service level instruments, local density attributes of a station, 
  #stockout indicator for nearby station
  #of the Z's need to get --local density attributes of a station, and stockout indicator for nearby station, while generating data, 
  #one entry for each data row.
  
  #distance to nearest station
  #average distance to nearest 5 stations
  st_data <- wdcMerged[!duplicated(wdcMerged$station_id_index),
    c("station_id_index","lat","lon")]
  st_data$dis_nearest <- 0
  st_data$dis_nearest_5 <- 0
  for(i in c(1:nrow(st_data))) {
    dis_v <- latlondistance(st_data$lat[i], st_data$lon[i],
                            st_data$lat[-i], st_data$lon[-i])
    st_data$dis_nearest[i] <- min(dis_v)
    st_data$dis_nearest_5[i] <- mean(sort(dis_v)[c(1:5)])
  }
  st_data$dis_nearest_sq <- st_data$dis_nearest * st_data$dis_nearest
  st_data$dis_nearest_5_sq <- st_data$dis_nearest_5 * st_data$dis_nearest_5
  st_data <- st_data[,c("station_id_index","dis_nearest","dis_nearest_5",
                        "dis_nearest_sq","dis_nearest_5_sq")]
  wdcMerged <- merge(wdcMerged, st_data, by="station_id_index")  
  #   #metro_den_on_1_2 for 0 to 300mts metro impact
  #   wdcMerged$metro_den_on_1_2 <- wdcMerged$metro_den_on_1 + wdcMerged$metro_den_on_2
  #   wdcMerged$touristlocs_den_1_2 <- wdcMerged$touristlocs_den_1 + wdcMerged$touristlocs_den_2
  #   wdcMerged$googleplaces_food_den_1_2 <- wdcMerged$googleplaces_food_den_1 + wdcMerged$googleplaces_food_den_2
  #   wdcMerged$googleplaces_government_den_1_2 <- wdcMerged$googleplaces_government_den_1 + wdcMerged$googleplaces_government_den_2
  #   wdcMerged$googleplaces_lodging_den_1_2 <- wdcMerged$googleplaces_lodging_den_1 + wdcMerged$googleplaces_lodging_den_2
  #   wdcMerged$googleplaces_museum_den_1_2 <- wdcMerged$googleplaces_museum_den_1 + wdcMerged$googleplaces_museum_den_2
  #   wdcMerged$googleplaces_movie_theater_den_1_2 <- wdcMerged$googleplaces_movie_theater_den_1 + wdcMerged$googleplaces_movie_theater_den_2
  
  metro_vars <- c("metro_den_on_a","metro_den_on_b","metro_den_on_c","metro_den_on_1","metro_den_on_2",
                  "metro_den_off_a","metro_den_off_b","metro_den_off_c","metro_den_off_1","metro_den_off_2") #not including _3 and _4 as they are quite correlated with intercept.  
  metro_dummies_df <- wdcMerged[,metro_vars]
  metro_dummies_df <- metro_dummies_df>0
  colnames(metro_dummies_df) <- paste("dummy_",metro_vars, sep="")
  wdcMerged <- cbind(wdcMerged, metro_dummies_df)
  
  
  tract_list <- sort(unique(wdcMerged$tract)) 
  tw_list <- levels(wdcMerged$tw_fac)
  exclude_tw <- tw_list[1]
  tw_list <- tw_list[-1]
  
  Xbase <- matrix(1,length(stocked_list),1)
  colnames(Xbase) <- "Intercept"
  
  for(tract_i in tract_list) {
    for(tw_j in tw_list) {
      Xbase_col <- matrix(0,length(stocked_list),1)
      Xbase_col[which(wdcMerged$tract[stocked_list]==tract_i & 
                        wdcMerged$tw[stocked_list]==tw_j)] <- 1
      Xbase_col[which(wdcMerged$tract[stocked_list]==tract_i & 
                        wdcMerged$tw[stocked_list]==exclude_tw)] <- -1
      colnames(Xbase_col) <- paste0("tract_tw_fac",tract_i,"_",tw_j)
      Xbase_col <- Xbase_col
      Xbase <- cbind(Xbase,Xbase_col)
    }
  }
  
  #distance from city center
  city_center <- c(48.857453, 2.351759)
  dis_center <- latlondistance(wdcMerged$lat[stocked_list], wdcMerged$lon[stocked_list], city_center[1], city_center[2])
  Xbase <- cbind(Xbase, "dis_center"=dis_center)
  Zbase <- Xbase
  #add two service level instruments, local density attributes of a station, stockout indicator for nearby station
  #to get Z
  list_X_covariates <- c("dis_nearest","dis_nearest_sq","dis_nearest_5","dis_nearest_5_sq","metro_den_on_a","metro_den_on_b","metro_den_on_c","metro_den_on_1","metro_den_on_2","metro_den_on_3","metro_den_on_4",
                       "metro_den_off_a","metro_den_off_b","metro_den_off_c","metro_den_off_1","metro_den_off_2","metro_den_off_3","metro_den_off_4",
                       paste("dummy_",metro_vars, sep=""),
                       "touristlocs_den_1","touristlocs_den_2","touristlocs_den_3","touristlocs_den_4","touristlocs_den_a","touristlocs_den_b","touristlocs_den_c")
  list_Z_covariates <- c("dis_nearest","dis_nearest_sq","dis_nearest_5","dis_nearest_5_sq","metro_den_on_a","metro_den_on_b","metro_den_on_c","metro_den_on_1","metro_den_on_2","metro_den_on_3","metro_den_on_4",
                         "metro_den_off_a","metro_den_off_b","metro_den_off_c","metro_den_off_1","metro_den_off_2","metro_den_off_3","metro_den_off_4",
                         paste("dummy_",metro_vars, sep=""),
                         "touristlocs_den_1","touristlocs_den_2","touristlocs_den_3","touristlocs_den_4","touristlocs_den_a","touristlocs_den_b","touristlocs_den_c")
  
  if(length(unique(wdcMerged$tract)) > 1) {
    Xbase <- cbind(Xbase,
      as.matrix(wdcMerged[stocked_list,c("census_density",list_X_covariates)]))
    Zbase <- cbind(Zbase,
      as.matrix(wdcMerged[stocked_list,c("census_density",list_Z_covariates)]))
  } else {
    #not including census_density
    Xbase <- cbind(Xbase,
      as.matrix(wdcMerged[stocked_list,list_X_covariates]))
    Zbase <- cbind(Zbase,
      as.matrix(wdcMerged[stocked_list,list_Z_covariates]))
  }
  
  #drop from Xbase columns which have no non-zero entry
  #rowsum is quite sure way to test
  if(length(which(colSums(abs(Xbase))==0))) {
    Xbase <- Xbase[,-which(colSums(abs(Xbase))==0)]     
  }
    
  if(length(which(colSums(abs(Zbase))==0))) {
    Zbase <- Zbase[,-which(colSums(abs(Zbase))==0)]     
  }
  
  #X <- Xbase
  #add service level vector to convert Xbase to X
#   X <- cbind(Xbase, serv_lvl=wdcMerged$serv_lvl[stocked_list], 
#              serv_lvl_sq=wdcMerged$serv_lvl[stocked_list]*wdcMerged$serv_lvl[stocked_list])
  X <- Xbase
  Z <- Zbase
    
  return(list("X"=X,              
              "Z"=Z))
}
