#source("constants_mw6.R")
#v3.0 two stockstate observations per stations, if results are consistent, can experiment with other fixed effects etc.
point_range <- 12
max_walking_dis <- 2.0
no_bootstrap_runs <- 1
reg_runs_save_file <- "reg_runs_discoef_2.89_8states_mw6_bootstrap.txt"


source("eval_func_4_GMM_new_temp_2.89.R")
source("data_prep_functions_ffdf_limstostate.R")
source("data_prep_inte_points_limstostate.R")
source("eval_func_2_cpp_cntrt_map_2.86_2.R")
source("eval_func_3_cpp_new_2.86_2.R")


load(file="wdcMerged_lin_weather_aggmonths_pr12.RData")
load(file="current_serv_lvl_lin_weather_aggmonths_pr12.RData")
load(file="user_serv_lvl_lin_weather_aggmonths_pr12.RData")
load(file="points_lin_weather_aggmonths_pr12.RData")
v0_vec <- c(0)

x0_start <<- NULL
prev_theta <<- NULL

#get local_stations
local_stations_vec <- wdcMerged[!duplicated(wdcMerged$station_id_index),c("station_id_index","local_stations")]
local_stations_vec <- local_stations_vec[order(local_stations_vec$station_id_index),]

tw_list <- unique(user_serv_lvl$tw)
instr_serv_lvl_1_df <- c()
no_st <- max(user_serv_lvl$station_id_index)  
for(tw_in in tw_list) {
  user_serv_lvl_subset <- user_serv_lvl[which(user_serv_lvl$tw==tw_in),]
  #if(!identical(order(user_serv_lvl_subset$station_id_index),c(1:nrow(user_serv_lvl_subset)))) stop("not in order")
  user_serv_lvl_subset <- user_serv_lvl_subset[order(user_serv_lvl_subset$station_id_index),]
  if(nrow(user_serv_lvl_subset)!=no_st) stop("incorrect number of rows in user_serv_lvl_subset")
  instr_serv_lvl_1_tw <- gen_serv_lvl_instr_neigh(user_serv_lvl_subset$serv_lvl,local_stations_vec$local_stations)
  instr_serv_lvl_1_tw_df <- data.frame(tw=user_serv_lvl_subset$tw, station_id_index=user_serv_lvl_subset$station_id_index,
                                       st_tw_index=user_serv_lvl_subset$st_tw_index,instr_serv_lvl=instr_serv_lvl_1_tw)
  instr_serv_lvl_1_df <- rbind(instr_serv_lvl_1_df, instr_serv_lvl_1_tw_df)
}
instr_serv_lvl_1_df <- instr_serv_lvl_1_df[order(instr_serv_lvl_1_df$tw,
                                                 instr_serv_lvl_1_df$station_id_index),]
if(!identical(user_serv_lvl$station_id_index, instr_serv_lvl_1_df$station_id_index)) stop("user_serv_lvl and instr_serv_lvl_1_df not in order")
if(!identical(user_serv_lvl$tw, instr_serv_lvl_1_df$tw)) stop("user_serv_lvl and instr_serv_lvl_1_df not in order")
user_serv_lvl$serv_lvl_neighbours <- instr_serv_lvl_1_df$instr_serv_lvl

#add serv_lvl instruments to wdcMerged
if(!identical(user_serv_lvl$st_tw_index, c(1:nrow(user_serv_lvl)))) stop("user_serv_lvl not in order")
wdcMerged$instr_serv_lvl <- user_serv_lvl$instr_serv_lvl[wdcMerged$st_tw_index]
wdcMerged$serv_lvl_neighbours <- user_serv_lvl$serv_lvl_neighbours[wdcMerged$st_tw_index]

rm(instr_serv_lvl_1_df)
rm(local_stations_vec)
rm(user_serv_lvl_subset)
rm(instr_serv_lvl_1_tw_df)

#save commands
# save(wdcMerged, file="wdcMerged_lin_weather.RData")
# save(current_serv_lvl, file="current_serv_lvl_lin_weather.RData")
# save(user_serv_lvl, file="user_serv_lvl_lin_weather.RData")
# save(points, file="points_lin_weather.RData")


