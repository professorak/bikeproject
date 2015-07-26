source("counterfactual_functions.R")

#coutnerfactual level 1 with multiplicative change in distance, serv_lvl
theta1 <- c(-18.01914170296,0,1.12975255452827,1054.73245975706,0.0999999905134055,1054.732)
delta_all <- compute_delta_list_cntrt_map_new(theta1,wdcMerged, points)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]

theta3 <- eval_error_xi_sl_MPEC(theta1,deltain_stin,wdcMerged,points)$theta3
  
ddc  <- demand_distance_change(1.1,theta1,delta_all, wdcMerged,points)

sdc <- demand_serv_lvl_change(theta1,delta_all, wdcMerged,points,0.9,
                                   theta3[which(rownames(theta3)=="serv_lvl")])

#subs_perc <- average_demand_susbtitution(c(theta1[1],0),wdcMerged,points,delta_list)
#subs_perc <- mean(perc_frac_dem_lost)

print("ddc: ")
print(ddc)
print((ddc[1]-ddc[2])/ddc[2]*100)
print("sdc: ")
print(sdc)
print((sdc[1]-sdc[2])/sdc[2]*100)
# print("subs_perc: ")
# print(subs_perc)
