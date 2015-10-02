eval_lambda_tw_groupin <- function(deltain, theta1, wdcMergedday, points, tw_groupin) {
  stocked_list <- which(wdcMergedday$stocked_out==FALSE)
  dem_T <- eval_lambda_new(deltain, theta1, wdcMergedday, points, tw_groupin)[stocked_list]  
  return(dem_T)
}


eval_lambda_full <- function(deltain_stin, theta1, wdcMerged, points) {
  
  #expand deltain to all observations, it is currently #stocked in observations
  deltain <- rep(-30, nrow(wdcMerged))
  deltain[which(wdcMerged$stocked_out==F)] <- deltain_stin
  
  tw_group_list <- unique(wdcMerged$tw_group)
  obj <- c()
  for(i in 1:length(tw_group_list)) {
    tw_groupin = tw_group_list[i]
    wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
    deltain_tw = deltain[which(wdcMerged$tw_group==tw_groupin)]    
    obj <- c(obj, 
             eval_lambda_tw_groupin(deltain_tw, theta1, wdcMergedday, 
                                    points, tw_groupin))
  }
  return(obj)
}




##debug_metro_effect.R functions

get_total_station_demand <- function() {
  source("data_estimation_2.6_weather_functions.R")
  #v3.0 two stockstate observations per stations, if results are consistent, can experiment with other fixed effects etc.
  point_range <- 3
  max_walking_dis <- 2.0
  dis_points <<- 0.050
  
  
  source("eval_func_4_GMM_new_temp_2.89.R")
  source("data_prep_functions_ffdf_limstostate.R")
  source("data_prep_inte_points_limstostate.R")
  source("eval_func_2_cpp_cntrt_map_2.86_2.R")
  source("eval_func_3_cpp_new_2.86_2.R")
  source("temp_read_google_places_data.R")
  
  #intergration points generation data
  min_lon <<- 2.293183
  max_lon <<- 2.374014
  min_lat <<- 48.83708
  max_lat <<- 48.884252
  
  #max_walking_dis = 0.1
  tract <- c(1:10)
  tract_in <- c(1:10)
  st_list <<- select_stations_in_tract(c(tract))
  
  #temporary - load st_6_7_8 file, should be able to generate this realtime, 
  save_dir <- paste0(csv_dir,"/ffdb/Paris/st_5_8")
  load(save_dir)
  st_list <<- intersect(st_5_8,st_list)
  
  filelist_all <- c("Paris/ind_paris_5.csv","Paris/ind_paris_6.csv",
                    "Paris/ind_paris_7.csv","Paris/ind_paris_8.csv")
  filelist_all_mons <- c(5:8)
  #   filelist_all <- c("Paris/ind_paris_7.csv")
  #   filelist_all_mons <- c(7)
  filelist_all_weather <- c("weatherdata_2013-5.csv","weatherdata_2013-6.csv",
                            "weatherdata_2013-7.csv","weatherdata_2013-8.csv")
  
  #filelist_servlvl <- c("Paris/ind_paris_4.csv","Paris/ind_paris_5.csv")
  wdcMerged <- c()
  for( i in 1:length(filelist_all_mons)) {
    filelist <- filelist_all[i]
    mon <- filelist_all_mons[i]
    save_file <- paste0(csv_dir,"/ffdb/wdcMerged_bootstrap_limstostate_5_8_den",
                        point_range,dis_points,"/mw6_tract_",paste0(tract,collapse="_"),"_mon_",mon)  
    print(save_file)
    #load files
    wdcMerged_bootstrap <- reg_run_bootstrap_simple(NULL,save_file)
    dim(wdcMerged_bootstrap)
    wdcMerged <- rbind(wdcMerged,wdcMerged_bootstrap)
    rm(wdcMerged_bootstrap)
  }  
  
  wdcMerged <- subset(wdcMerged, month>=5 & month<=8)
  wdcMerged <- droplevels(wdcMerged)
  
  wdcMerged_save <- wdcMerged
  wdcMerged <- subset(wdcMerged, stocked_out==F)
  wdcMerged_binned_sum <- binned_sum(wdcMerged$out_dem, bin=wdcMerged$station_id_index)
  wdcMerged_binned_sum <- as.data.frame(wdcMerged_binned_sum)
  wdcMerged_binned_sum$out_dem_sum <- wdcMerged_binned_sum$sum
  wdcMerged_binned_sum$obs_weight <- wdcMerged_binned_sum$count
  wdcMerged_binned_sum$out_dem_mean <- wdcMerged_binned_sum$out_dem_sum/wdcMerged_binned_sum$obs_weight
  wdcMerged_binned_sum$station_id_index <- as.numeric(row.names(wdcMerged_binned_sum))
  #get "station_id","lat","lon"
  station_data <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index")])
  station_data <- station_data[order(station_data$station_id_index),]
  wdcMerged_binned_sum <- wdcMerged_binned_sum[order(wdcMerged_binned_sum$station_id_index),]
  wdcMerged_binned_sum <- cbind(wdcMerged_binned_sum,station_data[,c("station_id","lat","lon")])
  
  
  #get points
  points <- generate_integration_points(wdcMerged_binned_sum)
  points <- subset(points, tract %in% c(1:10))
  points$weight <- points$density
  points$density <- NULL #changing density to weight to track where density is being used
  points$weight[which(points$type==2)] <- points$weight[which(points$type==2)]/max(points$weight[which(points$type==2)])
  points <- points[order(points$type, points$lat, points$lon),]  
  #google places data
  places_data <- read_googleplaces_data_dis_points()
  places_data <- subset(places_data, tract %in% tract_in)
  #places_cols_select <- c("lat","lon","cafe", "grocery_or_supermarket", "local_government_office")
  places_cols_select <- c("lat","lon","places_count","cafe", "grocery_or_supermarket","local_government_office","food")
  places_data <- places_data[,places_cols_select]
  places_data <- places_data[order(places_data$lat, places_data$lon),]
  if(!identical(round(points$lat[c(1:nrow(places_data))],4), round(places_data$lat,4))) stop("points lat dont match")
  if(!identical(round(points$lon[c(1:nrow(places_data))],4), round(places_data$lon,4))) stop("points lon dont match")
  places_data_temp <- places_data
  places_data_temp$lat <- NULL
  places_data_temp$lon <- NULL
  places_colnames <- colnames(places_data_temp)
  places_data_temp_full <- as.data.frame(matrix(0,nrow=nrow(points),ncol=ncol(places_data_temp)))
  colnames(places_data_temp_full) <- places_colnames
  places_data_temp_full[c(1:nrow(places_data_temp)),] <- places_data_temp
  points <- cbind(points, places_data_temp_full)
  
  
  #get local attributes  
  wdcCenTrParsed <- readtractboundaryfile()
  wdcMerged_binned_sum$tw <- 1  
  wdcMerged_binned_sum$tract  <- 
    assign_points_tract(wdcMerged_binned_sum$lat,wdcMerged_binned_sum$lon,wdcCenTrParsed)
  wdcMerged_binned_sum$sto_state_local <- 0
  wdcMerged_binned_sum <- get_local_attributes_st_state(wdcMerged_binned_sum, points)
  
  return(wdcMerged_binned_sum)    
}

generate_points_demand <- function(radius,wdcMerged_st_demand,points) {
  demand_vec <- rep(0,nrow(points))
  no_stations_vec <- rep(0,nrow(points))
  
  for(i in 1:nrow(points)) {
    lat1 = points$lat[i]
    lon1 = points$lon[i]
    dis_v <- latlondistance(lat1,lon1,wdcMerged_st_demand$lat,wdcMerged_st_demand$lon)  
    loc_st_i <- which(dis_v <= radius)
    demand_vec[i] <- sum(wdcMerged_st_demand$out_dem_mean[loc_st_i])
    no_stations_vec[i] <- length(loc_st_i)
  }
  return(cbind(demand_vec,no_stations_vec))
}



