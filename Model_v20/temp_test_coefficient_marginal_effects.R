source("counterfactual_functions.R")

#coutnerfactual level 1 with multiplicative change in distance, serv_lvl
theta1 <- c(-18.01914170296,0,1.12975255452827,1054.73245975706,0.0999999905134055)
delta_all <- compute_delta_list_cntrt_map_new(theta1,wdcMerged, points)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]

list_dis_factors <- c(-1,-2,-4,-8,-16,-32)/theta1[1]
list_dis_factors <- c(list_dis_factors,list_dis_factors*1.1)
ddc_vec <- c()
for(sc_factor in list_dis_factors) {
  ddc  <- demand_distance_change(sc_factor,theta1,delta_all, wdcMerged,points)  
  ddc_vec <- rbind(ddc_vec, c(sc_factor, ddc))
}

theta3 <- eval_error_xi_sl_MPEC(theta1,deltain_stin,wdcMerged,points)$theta3
serv_lvl_coef <- theta3[which(rownames(theta3)=="serv_lvl")]

list_sl_factors <- c(0.1,0.2,0.4,0.8,1.6,3.2)/serv_lvl_coef
list_sl_factors <- c(list_sl_factors,list_sl_factors*1.1)

sdc_vec <- c()
for(sl_factor in list_sl_factors) {
  sdc <- demand_serv_lvl_change(theta1,delta_all, wdcMerged,points,sl_factor,
                                serv_lvl_coef)  
  sdc_vec <- rbind(sdc_vec, c(sl_factor, sdc))
}

