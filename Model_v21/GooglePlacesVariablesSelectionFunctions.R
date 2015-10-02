source("eval_func_4_GMM_new_temp_2.89.R")
source("data_prep_functions_ffdf_limstostate.R")
source("data_prep_inte_points_limstostate.R")
source("eval_func_2_cpp_cntrt_map_2.86_2.R")
source("eval_func_3_cpp_new_2.86_2.R")

getTotalStationDemand <- function() {
  #v3.0 two stockstate observations per stations, if results are consistent, can experiment with other fixed effects etc.
  point_range <- 3
  no_bootstrap_runs <- 1
  reg_runs_save_file <- "reg_runs_discoef_2.89_8states_mw6_bootstrap.txt"
  
  
  
  #max_walking_dis = 0.1
  tract <- c(1:10)
  st_list = select_stations_in_tract(c(tract))
  
  #temporary - load st_6_7_8 file, should be able to generate this realtime, 
  save_dir <- paste0(csv_dir,"/ffdb/Paris/st_5_8")
  load(save_dir)
  st_list <- intersect(st_5_8,st_list)
  
  filelist_all <- c("Paris/ind_paris_6.csv",
                    "Paris/ind_paris_7.csv","Paris/ind_paris_8.csv")
  filelist_all_mons <- c(6:8)
  
  wdcMerged <- c()
  for( i in 1:length(filelist_all_mons)) {
    filelist <- filelist_all[i]
    mon <- filelist_all_mons[i]
    save_file <- paste0(csv_dir,"/ffdb/wdcMerged_bootstrap_limstostate_5_8_den",
                        point_range,"/mw6_tract_",paste0(tract,collapse="_"),"_mon_",mon)  
    print(save_file)
    #load files
    wdcMerged_bootstrap <- reg_run_bootstrap_simple(NULL,save_file)
    dim(wdcMerged_bootstrap)
    wdcMerged <- rbind(wdcMerged,wdcMerged_bootstrap)
    rm(wdcMerged_bootstrap)
  }    
  #use ffdfdply to find sum of out_dem variable by station_id_index
  station_demand <- binned_sum(wdcMerged$out_dem,wdcMerged$station_id_index, FUN=mean)  
  station_demand <- data.frame(station_id_index=as.numeric(rownames(station_demand)),
                  out_dem =station_demand[,"sum"]/station_demand[,"count"])
  
  station_data <- wdcMerged[ffwhich(wdcMerged, 
                                    !duplicated(wdcMerged$station_id_index)),]
  station_data <- station_data[,c("station_id_index","lat","lon")]
  station_demand <- merge(station_demand, station_data, by="station_id_index")
  return(station_demand)
}

getPointsData <- function(station_demand) {
  #for each points, find out stations that are in locality and store them as a string
  station_data <- unique(station_demand[,c("lat","lon","station_id_index")])
  station_data <- station_data[order(station_data$station_id_index),]  
#   station_data$local_points <- rep("",nrow(station_data))
  
  points <- generate_integration_points(station_data)  
#   points$local_stations <- rep(NA,nrow(points))
#   for(i in 1:nrow(points)) {
#     lat1 = points$lat[i]
#     lon1 = points$lon[i]
#     dis_v <- latlondistance(lat1,lon1,station_data$lat,station_data$lon)  
#     order_dis_v <- order(dis_v)
#     
#     points$local_stations[i] =paste(station_data$station_id_index[sort(order_dis_v[which(dis_v[order_dis_v[1:point_range]] <= max_walking_dis)])]
#                                     , collapse="_")
# #     for(stid in sort(order_dis_v[which(dis_v[order_dis_v[1:point_range]] <= max_walking_dis)])) {
# #       #add local points to stid
# #       if(station_data$local_points[stid]!="") {
# #         station_data$local_points[stid] <- paste(c(station_data$local_points[stid],i),collapse="_")
# #       } else {
# #         station_data$local_points[stid] <- i
# #       }    
# #     }
#   }
#   station_data <- station_data[,c("station_id_index","local_points")]
#   station_demand <- merge(station_demand,station_data, by=c("station_id_index"))
  return(list(station_demand=station_demand,
              points=points))
}


