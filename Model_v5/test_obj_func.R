#This file generates data for observations weights=1 and assumed values of deltas and 
#then estimates them

source("data_estimation_discoef_2.89_limstostate_mw6_saved.R")
points_save <- points
wdcMerged_store <- wdcMerged
user_serv_lvl_store <- user_serv_lvl
#make the weights to 1
#put a value of delta
#compute lambda
wdcMerged <- wdcMerged_store
user_serv_lvl <- user_serv_lvl_store
points <- points_save
#keeping "6_0" "6_3" and "7_0" "7_3"
wdcMerged <- subset(wdcMerged, tw_group_fac %in% c("6_0","6_3","7_0","7_3"))
wdcMerged <- droplevels(wdcMerged)
levels(wdcMerged$tw_group_fac)
user_serv_lvl <- subset(user_serv_lvl, tw_group %in% c("6_0","6_3","7_0","7_3"))
user_serv_lvl <- droplevels(user_serv_lvl)
unique(user_serv_lvl$tw_group)
#making 0 demand observations insignificant
wdcMerged$obs_weight[which(wdcMerged$stocked_out==F & wdcMerged$out_dem_sum<=0.1)] <- 0.1

# wdcMerged$out_dem_sum <- wdcMerged$out_dem_sum/wdcMerged$obs_weight
# wdcMerged$obs_weight <- 1
#removing really low demand
wdcMerged$out_dem_sum[which(wdcMerged$out_dem_sum<=0.01 & 
                              wdcMerged$stocked_out==FALSE)] <- 0.01

theta <- c(-3.68,100,5688)
a <- eval_obj_GMM_list_extended_new (theta, wdcMerged, points)  

theta2 <- c(-3.68,1,0)
a2 <- eval_obj_GMM_list_extended_new (theta2, wdcMerged, points)  
