target_use <- displace_use_vec[1]*use_actual
displace_serv_lvl <- displace_serv_lvl_vec[1]

optimal_distance_scaling(theta1, delta_list_cfl1, 
  wdcMerged_cfl1, points, displace_serv_lvl, serv_lvl_coef, target_use)

wdcMerged_save <- wdcMerged

delta_list <- delta_list_cfl1
wdcMerged <- wdcMerged_cfl1
demand_actual <- target_use
  
optimal_distance_scaling (theta1, delta_list, wdcMerged,points, 
                                     displace_serv_lvl, serv_lvl_coef, demand_actual)

a <- eval_lambda_full_unweighted(delta_list_cfl1, theta1,wdcMerged, points)

#because of obs_weights being 0, when wieghted average is taken, delta results in 0.

b <- wdcMerged$obs_weight[which(is.na(a))]
summary(b)
summary(wdcMerged$station_id_index[which(is.na(a))])
summary(wdcMerged$tw[which(is.na(a))])

c <- wdcMerged$obs_weight[-which(is.na(a))]
summary(c)



wdcMerged_0 <- wdcMerged
wdcMerged_0$obs_weight <- 0
d <- eval_lambda_full_unweighted(delta_list_cfl1, theta1,wdcMerged_0, points)









