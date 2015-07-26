theta1 <- c(-4.637460685)
theta1 <- c(theta1[1],0)
theta2 <- eval_error_xi(theta1,wdcMerged,points)$theta2
delta_list <- compute_delta_list_cntrt_map_new(theta1,wdcMerged,points)


#here write the counterfactual fucntion for effect on net demand when one of the station is missing.
#what to do about stockouts, in first, do a simple thing, like you do for substituion and assume
#complete sotkcout in and average values of deltas





ret_demand_lost_st <- demand_lost_st(theta1,wdcMerged,points,delta_list)
tot_dem <- sum(ret_demand_lost_st[[2]])
#perc_frac_dem_lost <- (-ret_demand_lost_st[[1]]/tot_dem*100*max(wdcMerged$station_id_index))
perc_frac_dem_lost <- -(ret_demand_lost_st[[1]]/ret_demand_lost_st[[2]])*100

hist(perc_frac_dem_lost, ylab="Number of Stations",
     xlab="Lost demand as a percentage of before demand",
     main="")

write.csv(perc_frac_dem_lost, file="perc_frac_dem_lost.csv")

##########################################
##########################################
#isodemand curve

#increase service level by 0.1 and cap them at 1
# displace_serv_lvl <- 0.1
# ret <- optimal_distance_scaling (displace_serv_lvl, wdcMerged,points)
  

# displace_serv_lvl <- 0.1
# ret <- optimal_distance_scaling (displace_serv_lvl, wdcMerged,points)
#theta1 <- c(-4.27167090930182,0)
displace_serv_lvl_vec <- c(-0.3,-0.2,-0.1, 0, 0.1, 0.2, 0.3)
ret_vec <- c()
demand_actual <- get_base_demand(0, wdcMerged,points, theta2)

for(displace_serv_lvl in displace_serv_lvl_vec) {
   ret_vec <- c(ret_vec,optimal_distance_scaling (displace_serv_lvl, wdcMerged,points,theta2,demand_actual))
}

write.csv(cbind(displace_serv_lvl_vec,ret_vec), file="isodemand_linear.csv")


#write the code for isodemand for other levels 
displace_serv_lvl_vec <- c(-0.3,-0.2,-0.1, 0, 0.1, 0.2, 0.3)

ret_vec_mat <- c()
for(base_displace_serv_lvl_vec in displace_serv_lvl_vec) {
	ret_vec <- c()
	demand_actual <- get_base_demand(base_displace_serv_lvl_vec, wdcMerged,points, theta2)
	for(displace_serv_lvl in displace_serv_lvl_vec) {
		ret_vec <- c(ret_vec,optimal_distance_scaling (displace_serv_lvl, wdcMerged,points,theta2,demand_actual))
	}
	ret_vec_mat	<- rbind(ret_vec_mat, ret_vec)
}

write.csv(ret_vec_mat, file="isodemand_linear_mat.csv")

########################################################
########################################################
#coutnerfactual level 1 with multiplicative change in distance, serv_lvl

ddc  <- demand_distance_change_new(1.1,c(theta1[1],0),delta_list, wdcMerged,points)
sdc <- demand_serv_lvl_change_new2(c(theta1[1],0),wdcMerged,points,0.9,
                                   theta2[which(rownames(theta2)=="serv_lvl_covar")])
#subs_perc <- average_demand_susbtitution(c(theta1[1],0),wdcMerged,points,delta_list)
subs_perc <- mean(perc_frac_dem_lost)

print("ddc: ")
print(ddc)
print("sdc: ")
print(sdc)
print("subs_perc: ")
print(subs_perc)

#change in system-demand. above was system-use
sdc2 <- demand_serv_lvl_change_new_systemdemand(c(theta1[1],0),wdcMerged,points,0.9,
                                   theta2[which(rownames(theta2)=="serv_lvl_covar")])
print("sdc2: ")
print(sdc2)
######################################
#level 1 counterfactuals with additive changes in distance and service level

#additive changes in distance
displace_distance_vec <- c(-0.3,-0.2,-0.1, 0, 0.1, 0.2, 0.3)

dem_distance_vec <- c()
for(base_displace_distance in displace_distance_vec) {  
  dem_distance_vec <- c(dem_distance_vec,demand_distance_additive_change(base_displace_distance, theta1, wdcMerged, points, delta_list))
}
dem_distance_mat <- dem_distance_vec
no_st <- max(wdcMerged$station_id_index)
dem_distance_vec <- matrix(dem_distance_vec, length(dem_distance_vec)/no_st, no_st, byrow=T)
dem_distance_vec <- rowSums(dem_distance_vec)

  #to confirm the rate of fall is nearly same ;; it is!
for(i in 2:length(dem_distance_vec)) {
  print(dem_distance_vec[i-1]/dem_distance_vec[i])
}
  #mean perc change in demand when distance is change by amounts in displace_distance_vec
mean_perc <- c()
for(i in 2:length(dem_distance_vec)) {
  mean_perc <- c(mean_perc,1-dem_distance_vec[i]/dem_distance_vec[i-1])
}
print(mean(mean_perc)*100)
mean_perc_dis <- mean_perc

################################
#additive change in service level
displace_serv_lvl_vec <- c(-0.3,-0.2,-0.1, 0, 0.1, 0.2, 0.3)
dem_servlvl_vec <- c()
for(base_displace_serv_lvl_vec in displace_serv_lvl_vec) {  
  dem_servlvl_vec <- c(dem_servlvl_vec,get_base_demand(base_displace_serv_lvl_vec, wdcMerged,points, theta2))
}
#to confirm the rate of fall is nearly same ;; it is!
for(i in 2:length(dem_servlvl_vec)) {
  print(dem_servlvl_vec[i-1]/dem_servlvl_vec[i])
}
#mean perc change in demand when distance is change by amounts in displace_distance_vec
mean_perc <- c()
for(i in 2:length(dem_servlvl_vec)) {
  mean_perc <- c(mean_perc,1-dem_servlvl_vec[i-1]/dem_servlvl_vec[i])
}
print(mean(mean_perc)*100)
mean_perc_serv <- mean_perc
#there is really a huge change in demand when going from 0 to 0.1, demand goes 6 times unlike other cases where it goes less than 1.5 times

  
################################
#counterfactual level 2

dis_scaling <- c(0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.20)

serv_lvl_values <- c(0.499578416947, 0.667279607877, 0.750356654283, 0.798240923192, 0.832845203737, 0.857831671781, 0.889291282838, 
                     0.899833962270, 0.917188255333, 
                     0.929012726490, 0.936217702815, 0.944900500743, 0.951098883953, 0.957015177987, 0.960315335619, 0.964528097363, 0.967168731425, 0.976991182653)
serv_lvl_base <- serv_lvl_values[which(dis_scaling==1)]
ddc_vec <- c()

for(i in 1:length(dis_scaling)) {
  ddc_vec <- c(ddc_vec,demand_distance_change_new_linear_serv_lvl_adjusted_tract(1/dis_scaling[i], theta1, 
          wdcMerged,points, serv_lvl_values[i]/serv_lvl_base,theta2[which(rownames(theta2) == "serv_lvl_covar")]))  
}
ddc_mat <- matrix(ddc_vec, length(ddc_vec)/10,10, byrow=T)
ddc_mat_summed <- rowSums(ddc_mat)
dis_scaling[which.max(ddc_mat_summed)]
ddc_mat_summed_normalized <- ddc_mat_summed/ddc_mat_summed[which(dis_scaling==1)]
ddc_mat_summed_normalized

for(tract_in in tract) {
  net_dem <- ddc_mat[,tract_in]
  print(dis_scaling[which.max(net_dem)])
}
write.csv(cbind(dis_scaling, serv_lvl_values, ddc_mat_summed_normalized),
          file="counterfactual_l2.csv")

###############################################################
###############################################################
#compute the fraction of commuters coming from different distances.
#do in R implementation first
#take the average of exp(delta) values for different states instead of delta




#### finer points_new
points_new <- generate_integration_points(dis_points=0.025)
#for each points_new, find out stations that are in locality and store them as a string
points_new$local_stations <- rep(NA,nrow(points_new))
for(i in 1:nrow(points_new)) {
  lat1 = points_new$lat[i]
  lon1 = points_new$lon[i]
  dis_v <- latlondistance(lat1,lon1,station_data$lat,station_data$lon)  
  order_dis_v <- order(dis_v)
  
  points_new$local_stations[i] =paste(station_data$station_id_index[sort(order_dis_v[which(dis_v[order_dis_v[1:point_range]] <= max_walking_dis)])]
                                      , collapse="_")
}  
market_share = 0.10
tract_in <- tract
points_new <- subset(points_new, tract %in% tract_in)
for(tract_in in tract) {
  wdcMerged_tract <- wdcMerged[which(wdcMerged$tract == tract_in),]
  #net_demand = max(as.numeric(by(wdcMerged_tract$out_dem_mean, wdcMerged_tract$tw_group_fac, FUN=sum)))
  wdcMerged_tract <- subset(wdcMerged_tract,stocked_out==F)
  #at each stationXtw_group_fac level create a out_dem_mean
  wdcMerged_tract$st_tw_group_fac <- as.factor(paste(wdcMerged_tract$station_id_index, wdcMerged_tract$tw_group_fac))
  wdcMerged_tract$out_dem_mean_dem <- ave(wdcMerged_tract$out_dem_sum, wdcMerged_tract$st_tw_group_fac, FUN=sum)
  wdcMerged_tract$out_dem_mean_obs <- ave(wdcMerged_tract$obs_weight, wdcMerged_tract$st_tw_group_fac, FUN=sum)
  wdcMerged_tract$out_dem_mean <- wdcMerged_tract$out_dem_mean_dem/wdcMerged_tract$out_dem_mean_obs
  wdcMerged_tract$dup <- duplicated(wdcMerged_tract$st_tw_group_fac)
  wdcMerged_tract <- subset(wdcMerged_tract, dup==F)
  dem_sum <- as.numeric(by(wdcMerged_tract$out_dem_mean, wdcMerged_tract$tw_group_fac, FUN=sum))
  net_demand <- max(dem_sum)      
  points_new$density[which(points_new$tract %in% tract_in)] = net_demand/market_share/length(which(points_new$tract %in% tract_in))
}

#
ret <- distance_segment_demand(theta1, wdcMerged, points_new, delta_list)
print(ret)

################################################################

