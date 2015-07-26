#test time
#data
source("eval_obj_func_GMM_model4.R")
#source("GetDataGooglePlaces.R")
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


source("data_estimation_2.6_weather_saved.R")
#source("data_estimation_2.6_weather_small_saved.R")
source("eval_func_2_cpp_MPEC.R")

#make the weights to 1
#put a value of delta
#compute lambda


#keeping "6_0" "6_3" and "7_0" "7_3"
time_windows <- c("6_0","6_5","7_0","7_5")
tws <- c("0","5")
#time_windows <- c("6_3")
wdcMerged <- subset(wdcMerged, tw_group_fac %in% time_windows)
wdcMerged <- droplevels(wdcMerged)
levels(wdcMerged$tw_group_fac)
user_serv_lvl <- subset(user_serv_lvl, tw %in% tws)
user_serv_lvl <- droplevels(user_serv_lvl)
unique(user_serv_lvl$tw)
# 
# wdcMerged <- droplevels(wdcMerged)
# user_serv_lvl$st_tw_index <- as.numeric(user_serv_lvl$st_tw)
# wdcMerged$st_tw_index <- as.numeric(wdcMerged$st_tw)
rm(current_serv_lvl)

# #making 0 demand observations insignificant
# wdcMerged$obs_weight[which(wdcMerged$stocked_out==F & wdcMerged$out_dem_sum<=0.1)] <- 0.1

# wdcMerged$out_dem_sum <- wdcMerged$out_dem_sum/wdcMerged$obs_weight
# wdcMerged$obs_weight <- 1
#removing really low demand
wdcMerged$out_dem_sum[which(wdcMerged$out_dem_sum<=0.01 & 
                              wdcMerged$stocked_out==FALSE)] <- 0.01





############################################################
#run
theta <- c(-2.66047062428, 0.29046474623, 1681.87287960098, 3.43903446529, 
           836.49733133688, 0.31363526842, 0.01000001459)
length_theta <- length(theta)
set.seed(34675)
deltain_stin <- rnorm(length(which(!wdcMerged$stocked_out)),-3, 2)
params <- c(theta,deltain_stin)  

a1_eval_constraints <- eval_g(params,wdcMerged, points, length(theta))
a1_eval_jac_constraints <- eval_jac_g(params,wdcMerged, points, length(theta))

####################################################################
####################################################################
#check for a single tw group
deltain_stin <- params[-c(1:length_theta)]
theta1 <- c(theta[1],0,theta[-1])  

#expand deltain to all observations, it is currently #stocked in observations
deltain_all <- rep(-30, nrow(wdcMerged))
deltain_all[which(wdcMerged$stocked_out==F)] <- deltain_stin

tw_group_list <- unique(wdcMerged$tw_group)
grad_constraints <- c()
i <- 1
tw_groupin = tw_group_list[i]
wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
deltain = deltain_all[which(wdcMerged$tw_group==tw_groupin)]
# eval_grad_constraints_MPEC_tw_groupin(deltain_tw, theta1, 
#     wdcMergedday,points, tw_groupin)
system.time({
  grad_theta <- eval_grad_lambda_theta_new(theta1, deltain, wdcMergedday, points, tw_groupin)
  grad_theta <- grad_theta[,-2]  
})
system.time({
  grad_delta <- eval_lambda_delta_list_new(deltain, theta1, wdcMergedday, points, tw_groupin)[[2]]
})

