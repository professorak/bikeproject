#demand 

source("constants_mw6.R")
source("data_estimation_2.6_weather_functions.R")
#v3.0 two stockstate observations per stations, if results are consistent, can experiment with other fixed effects etc.
point_range <- 6
max_walking_dis <- 2.0
max_points_distance <- 0.6 #600mts is the maximum a point can be far way from closest station.
dis_points <<- 0.050
no_bootstrap_runs <- 1
reg_runs_save_file <- "reg_runs_discoef_2.89_8states_mw6_bootstrap.txt"
num_top_states_select <- 8

source("eval_func_4_GMM_new_temp_2.89.R")
source("data_prep_functions_ffdf_limstostate.R")
source("data_prep_inte_points_limstostate.R")
source("eval_func_2_cpp_cntrt_map_2.86_2.R")
source("eval_func_3_cpp_new_2.86_2.R")
source("temp_read_google_places_data.R")

#intergration points generation data
min_lon <<- 2.249428
max_lon <<- 2.46408
min_lat <<- 48.81626
max_lat <<- 48.89964

#max_walking_dis = 0.1
tract <- c(1:20)

st_list <<- select_stations_in_tract(c(tract))

#temporary - load st_6_7_8 file, should be able to generate this realtime, 
save_dir <- paste0(csv_dir,"/ffdb/Paris/st_5_8")
load(save_dir)
st_list <<- intersect(st_5_8,st_list)

#stations which are in interior of paris
save_dir <- paste0(csv_dir,"/ffdb/Paris/paris_stations.RData")
load(save_dir)
st_list <<- intersect(paris_stations,st_list)


filelist_all <- c("Paris/ind_paris_5.csv","Paris/ind_paris_6.csv",
                  "Paris/ind_paris_7.csv","Paris/ind_paris_8.csv")
filelist_all_mons <- c(5:8)
#   filelist_all <- c("Paris/ind_paris_7.csv")
#   filelist_all_mons <- c(7)
filelist_all_weather <- c("weatherdata_2013-5.csv","weatherdata_2013-6.csv",
                          "weatherdata_2013-7.csv","weatherdata_2013-8.csv")

#filelist_servlvl <- c("Paris/ind_paris_4.c

#stid_list1 <- c(9004,9005,18041,9018,18042,18043,9027,18004,18114,18003,18001,18113,18111,18017,18016,18015,18002,18006,18005)
stid_list1 <- c(15118 ,15002, 9106, 9018, 8017, 4018, 10016, 
                10007, 8039, 8013, 8030, 8010)
st_list <<- stid_list1
station_data <- c()
for( i in 1:length(filelist_all_mons)) {
  filelist <- filelist_all[i]
  mon <- filelist_all_mons[i]
  save_file <- paste0(csv_dir,"/ffdb/wdcMerged_bootstrap_limstostate_5_8_den",
                      point_range,dis_points,"/mw6_tract_",paste0(tract,collapse="_"),"_mon_",mon)  
  print(save_file)
  
  
  #load files
  wdcRawOrg <- readcsvdemand(NULL,filelist) 
  #wdcRawOrg <- subset(wdcRawOrg, station_id==stid)
  station_data <- rbind(station_data, wdcRawOrg)
}

station_data <- as.data.frame(station_data)

save(station_data, file=paste0("station_data_stid_list2",".RData"))



