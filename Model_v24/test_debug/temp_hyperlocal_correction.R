dis_v <- latlondistance(wdcMerged$lat[which(wdcMerged$station_id_index==326)[1]],
                        wdcMerged$lon[which(wdcMerged$station_id_index==326)[1]],
                        points$lat,
                        points$lon)

which.min(dis_v)
dis_v[which.min(dis_v)]
points$type[which.min(dis_v)]



sourceCpp("eval_func_2_hyperlocal_nldis_v0weights_set2_debug.cpp")

ret <- eval_share_log_list_new(x0,theta1,wdcMergedday,points,tw_groupin)
dem_T <- ret$dem_T
dem_hat_T <- ret$dem_hat_T      

which(wdcMergedday$station_id_index==326)
dem_T[which(wdcMergedday$station_id_index==326)]

###
#use apply to get the ordered sto_state and stations

temp_cols <- t(apply(wdcMerged[,c("local_stations", "sto_state_local")],
                     1, FUN=sort_sto_state_local))
wdcMerged$local_stations <- temp_cols[,1]
wdcMerged$sto_state_local <- temp_cols[,2]

  



