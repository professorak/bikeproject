source("counterfactual_functions.R")

#coutnerfactual level 1 with multiplicative change in distance, serv_lvl
theta1 <- c(-12.28769991999, 4.20013541466, 0.09999999003, 282.66078886268, 8.10824517514, 7505.59969617232)
delta_all <- compute_delta_list_cntrt_map_new(theta1,wdcMerged, points)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]

theta3 <- eval_error_xi_sl_MPEC(theta1,deltain_stin,wdcMerged,points)$theta3
  
ddc  <- demand_distance_change(1.1,theta1,delta_all, wdcMerged,points)

sdc <- demand_serv_lvl_change(theta1,delta_all, wdcMerged,points,0.9,
                                   theta3[which(rownames(theta3)=="serv_lvl")])

subs_perc <- average_demand_susbtitution_demModel(theta1,wdcMerged,points,delta_all)


print("ddc: ")
print(ddc)
print((ddc[1]-ddc[2])/ddc[2]*100)
print("sdc: ")
print(sdc)
print((sdc[1]-sdc[2])/sdc[2]*100)
# print("subs_perc: ")
# print(subs_perc)
