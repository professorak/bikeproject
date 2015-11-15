#This file generates data for observations weights=1 and assumed values of deltas and 
#then estimates them

#global variables
colnames_theta1 <<- c("dis coef","rand coef","density ridership","density metro","density intercept",
                      "density metro evening", "density google places count",
                      "density bus")
density_ridership_col <<- 3
density_metro_col <<- 4
density_intercept_col <<- 5
density_metro_evening_col <<- 6
density_google_places_count <<- 7
density_bus_col <<- 8

source("GetDataGooglePlaces_functions.R")
source("data_estimation_2.6_weather_20dis_aggmonth_hyperlocalstate_norealloc_stkoutthresh_five_averaged_pr3_saved.R")
#source("data_estimation_2.6_weather_small_saved.R")
source("eval_func_2_cpp_MPEC.R")

#make the weights to 1
#put a value of delta
#compute lambda


# #keeping "6_0" "6_3" and "7_0" "7_3"
# time_windows <- c("6_0","6_5","7_0","7_5")
# tws <- c("0","5")
# #time_windows <- c("6_3")
# wdcMerged <- subset(wdcMerged, tw_group_fac %in% time_windows)
# wdcMerged <- droplevels(wdcMerged)
# levels(wdcMerged$tw_group_fac)
# user_serv_lvl <- subset(user_serv_lvl, tw %in% tws)
# user_serv_lvl <- droplevels(user_serv_lvl)
# unique(user_serv_lvl$tw)
# 
# wdcMerged <- droplevels(wdcMerged)
# user_serv_lvl$st_tw_index <- as.numeric(user_serv_lvl$st_tw)
# wdcMerged$st_tw_index <- as.numeric(wdcMerged$st_tw)
rm(current_serv_lvl)

temp_cols <- t(apply(wdcMerged[,c("local_stations", "sto_state_local")],
      1, FUN=sort_sto_state_local))
wdcMerged$local_stations <- temp_cols[,1]
wdcMerged$sto_state_local <- temp_cols[,2]

# #making 0 demand observations insignificant
# wdcMerged$obs_weight[which(wdcMerged$stocked_out==F & wdcMerged$out_dem_sum<=0.1)] <- 0.1

# wdcMerged$out_dem_sum <- wdcMerged$out_dem_sum/wdcMerged$obs_weight
# wdcMerged$obs_weight <- 1

wdcMerged <<- removemodify_outlier_states(wdcMerged)

  


  
  
