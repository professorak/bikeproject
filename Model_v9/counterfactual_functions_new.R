

demand_distance_change_inner_new <- function(scale_distance_factor, theta1, delta_list, wdcMerged, points) {
  theta1[1] <- theta1[1]/scale_distance_factor
  tw_group_list <- c()
  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  
  lambda_t <- c()
  for(tw_groupin in tw_groupin_list) {
    wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
    deltain_tw <- delta_list[which(wdcMerged$tw_group==tw_groupin)]
    lambda_t <- c(lambda_t,
                              eval_lambda_delta_list_new(deltain_tw, theta1,wdcMergedday, points, tw_groupin)$objective)      
  }
  return(sum(lambda_t))
}

demand_distance_change_inner_new_serv_lvl_adjusted <- function(scale_distance_factor, theta1, delta_list, wdcMerged, points, 
                                                               serv_lvl_covar_new, serv_lvl_covar) {
  theta1[1] <- theta1[1]/scale_distance_factor
  tw_group_list <- c()
  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  
  lambda_t <- c()
  for(tw_groupin in tw_groupin_list) {
    wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
    deltain_tw <- delta_list[which(wdcMerged$tw_group==tw_groupin)]
    lambda_t <- c(lambda_t,
                  eval_lambda_delta_list_new(deltain_tw, theta1,wdcMergedday, points, tw_groupin)$objective)      
  }
  return(sum(lambda_t*serv_lvl_covar_new/serv_lvl_covar))
}

demand_distance_change_inner_new_nls <- function(scale_distance_factor, theta1, delta_list, wdcMerged, points) {
  theta1[1] <- theta1[1]/scale_distance_factor
  theta1[3] <- theta1[3]/(scale_distance_factor^2)
  
  tw_group_list <- c()
  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  
  lambda_t <- c()
  for(tw_groupin in tw_groupin_list) {
    wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
    deltain_tw <- delta_list[which(wdcMerged$tw_group==tw_groupin)]
    lambda_t <- c(lambda_t,
                  eval_lambda_delta_list_new(deltain_tw, theta1,wdcMergedday, points, tw_groupin)$objective)      
  }
  return(sum(lambda_t))
}

demand_distance_change_inner_new_nls_serv_lvl_adjusted <- function(scale_distance_factor, theta1, delta_list, wdcMerged, points, serv_lvl_covar_new, serv_lvl_covar) {
  theta1[1] <- theta1[1]/scale_distance_factor
  theta1[3] <- theta1[3]/(scale_distance_factor^2)
  
  tw_group_list <- c()
  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  
  lambda_t <- c()
  for(tw_groupin in tw_groupin_list) {
    wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
    deltain_tw <- delta_list[which(wdcMerged$tw_group==tw_groupin)]
    lambda_t <- c(lambda_t,
                  eval_lambda_delta_list_new(deltain_tw, theta1,wdcMergedday, points, tw_groupin)$objective)      
  }
  return(sum(lambda_t*serv_lvl_covar_new/serv_lvl_covar))
}

demand_distance_change_inner_new_nls_tract <- function(scale_distance_factor, theta1, delta_list, wdcMerged, points) {
  theta1[1] <- theta1[1]/scale_distance_factor
  theta1[3] <- theta1[3]/(scale_distance_factor^2)
  
  tw_group_list <- c()
  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  
  lambda_t <- c()
  for(tw_groupin in tw_groupin_list) {
    wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
    deltain_tw <- delta_list[which(wdcMerged$tw_group==tw_groupin)]
    lambda_t <- c(lambda_t,
                  eval_lambda_delta_list_new(deltain_tw, theta1,wdcMergedday, points, tw_groupin)$objective)      
  }
  return(c(by(lambda_t, wdcMerged$tract, FUN=sum)))
}

demand_distance_change_inner_new_nls_tract_serv_lvl_adjusted <- function(scale_distance_factor, theta1, delta_list, wdcMerged, points, 
                                                               serv_lvl_covar_new, serv_lvl_covar) {
  theta1[1] <- theta1[1]/scale_distance_factor
  theta1[3] <- theta1[3]/(scale_distance_factor^2)
  serv_lvl_covar_new <- serv_lvl_covar_new + 1e-4
  serv_lvl_covar <- serv_lvl_covar + 1e-4 
  
  tw_group_list <- c()
  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  
  lambda_t <- c()
  for(tw_groupin in tw_groupin_list) {
    wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
    deltain_tw <- delta_list[which(wdcMerged$tw_group==tw_groupin)]
    lambda_t <- c(lambda_t,
                  eval_lambda_delta_list_new(deltain_tw, theta1,wdcMergedday, points, tw_groupin)$objective)      
  }
  return(c(by(lambda_t*serv_lvl_covar_new/serv_lvl_covar, wdcMerged$tract, FUN=sum)))
}

demand_distance_change_inner_new_linear_tract <- function(scale_distance_factor, theta1, delta_list, wdcMerged, points) {
  theta1[1] <- theta1[1]/scale_distance_factor
  
  tw_group_list <- c()
  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  
  lambda_t <- c()
  for(tw_groupin in tw_groupin_list) {
    wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
    deltain_tw <- delta_list[which(wdcMerged$tw_group==tw_groupin)]
    lambda_t <- c(lambda_t,
                  eval_lambda_delta_list_new(deltain_tw, theta1,wdcMergedday, points, tw_groupin)$objective)      
  }
  return(c(by(lambda_t, wdcMerged$tract, FUN=sum)))
}

demand_distance_change_new <- function(scale_distance_factor, theta1, delta_list, wdcMerged, points) {
  demand_actual <- demand_distance_change_inner_new(1, theta1, delta_list, wdcMerged, points)
  demand_changed <- demand_distance_change_inner_new(scale_distance_factor, theta1, delta_list, wdcMerged, points)
  return(c(demand_changed,demand_actual))
}  

demand_distance_change_new_nls <- function(scale_distance_factor, theta1, delta_list, wdcMerged, points) {
  demand_actual <- demand_distance_change_inner_new_nls(1, theta1, delta_list, wdcMerged, points)
  demand_changed <- demand_distance_change_inner_new_nls(scale_distance_factor, theta1, delta_list, wdcMerged, points)
  return(c(demand_changed,demand_actual))
}  

demand_distance_change_new_nls_tract <- function(scale_distance_factor, theta1, delta_list, wdcMerged, points) {
  demand_changed <- demand_distance_change_inner_new_nls_tract(scale_distance_factor, theta1, delta_list, wdcMerged, points)
  return(c(demand_changed))
}

demand_distance_change_new_linear_tract <- function(scale_distance_factor, theta1, delta_list, wdcMerged, points) {
  demand_changed <- demand_distance_change_inner_new_linear_tract(scale_distance_factor, theta1, delta_list, wdcMerged, points)
  return(c(demand_changed))
}


demand_substitution <- function(focal_station_index, mode, theta1, wdcMerged, points,delta_list) {
  #find the number of people that substitute from a station to another
  #fix a station as focal station
  #select only points local to it. 
  #in first case select assume this station is stocked and find the demand at nearby stations
  #in next case take this station out from choice set and find demand and compare total demands
  #in each case
  
  no_st <- max(wdcMerged$station_id_index)
  tw_group_list <- unique(wdcMerged$tw_group)
  
  station_data <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index",
                                      "local_stations")])
  station_data <- station_data[order(station_data$station_id_index),]
  
  #   points_sub <- station_data$local_points[focal_station_index] 
  points_sub <- which(stid_in_localstations(focal_station_index, points$local_stations))
  points_sub <- points[points_sub,]
  
  
  
  beta1 <- theta1[1]
  sigma0 <- theta1[2]
  lat1 <- station_data[,"lat"]
  lon1 <- station_data[,"lon"]
  lambda_t <- rep(0,nrow(station_data)) 
  
  for(tw_groupin in tw_group_list) {
    deltain <- delta_list[which(wdcMerged$tw_group == tw_groupin)]
    deltain <- as.numeric(by(deltain, wdcMerged$station_id_index[which(wdcMerged$tw_group == tw_groupin)], FUN=mean))
    for(i in 1:nrow(points_sub)) {
      lat2 <- points_sub[i,"lat"]
      lon2 <- points_sub[i,"lon"]
      #dis_v <- latlondistance(lat1, lon1, lat2, lon2)
      st_point_list <- splitchar(points_sub$local_stations[i])
      list_obs <- st_point_list
      if(length(list_obs)==0) next
      
      station_data_i <- station_data[st_point_list,]
      dis_v <- latlondistance(points_sub[i,"lat"], points_sub[i,"lon"], 
                              station_data_i$lat, station_data_i$lon)
      sto_state = rep(0,nrow(station_data_i))
      if(mode==1) {
        sto_state[which(station_data_i$station_id_index==focal_station_index)]=1
      }
      util <- (beta1*dis_v) + deltain[st_point_list]
      exputil <- exp(util)
      util_st <- (!sto_state) * exputil
      den_util <- sum(util_st)
      lambda_st_t <- rep(0,length(util_st))
      for(k in 1:length(v0_vec)) {
        out <- exp(-v0_vec[k]*sigma0)
        prob_t <- util_st/(out+den_util)*(points_sub[i,"density"]/length(v0_vec))
        lambda_st_t <- lambda_st_t + prob_t
      }
      lambda_t[list_obs] <- lambda_t[list_obs] + lambda_st_t
    }
  }
  return(lambda_t)
}

demand_substitution_demModel <- function(focal_station_index, mode, theta1, wdcMerged, points,delta_list) {
  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  
  lambda_t <- c()
  #     no_st <- max(wdcMerged$station_id_index)
  #     tw_group_list <- unique(wdcMerged$tw_group)
  
  #     station_data <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index",
  #                                         "local_stations")])
  #     station_data <- station_data[order(station_data$station_id_index),]
  
  #   points_sub <- station_data$local_points[focal_station_index] 
  points_sub <- which(stid_in_localstations(focal_station_index, points$local_stations))
  points_sub <- points[points_sub,]    
  
  for(tw_groupin in tw_groupin_list) {
    wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
    deltain_tw <- delta_list[which(wdcMerged$tw_group==tw_groupin)]
    if(mode==1) {
      deltain_tw[which(wdcMergedday$station_id_index==focal_station_index)] <- -100
    }
    lambda_t <- c(lambda_t,
                  eval_lambda_delta_list_new(deltain_tw, theta1,wdcMergedday, points_sub, tw_groupin)$objective)      
  }
  return((lambda_t))
}

demand_computation <- function(theta1, wdcMerged, points,delta_list) {
  #find the number of people that substitute from a station to another
  #fix a station as focal station
  #select only points local to it. 
  #in first case select assume this station is stocked and find the demand at nearby stations
  #in next case take this station out from choice set and find demand and compare total demands
  #in each case
  
  no_st <- max(wdcMerged$station_id_index)
  tw_group_list <- unique(wdcMerged$tw_group)
  
  station_data <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index",
                                      "local_stations")])
  station_data <- station_data[order(station_data$station_id_index),]
  
  #   points_sub <- station_data$local_points[focal_station_index] 
  #points_sub <- which(stid_in_localstations(focal_station_index, points$local_stations))
  points_sub <- points
  
  beta1 <- theta1[1]
  sigma0 <- theta1[2]
  lat1 <- station_data[,"lat"]
  lon1 <- station_data[,"lon"]
  lambda_t <- rep(0,nrow(station_data)) 
  
  for(tw_groupin in tw_group_list) {
    deltain <- delta_list[which(wdcMerged$tw_group == tw_groupin)]
    deltain <- as.numeric(by(deltain, wdcMerged$station_id_index[which(wdcMerged$tw_group == tw_groupin)], FUN=mean))
    for(i in 1:nrow(points_sub)) {
      lat2 <- points_sub[i,"lat"]
      lon2 <- points_sub[i,"lon"]
      #dis_v <- latlondistance(lat1, lon1, lat2, lon2)
      st_point_list <- splitchar(points_sub$local_stations[i])
      list_obs <- st_point_list
      if(length(list_obs)==0) next
      
      station_data_i <- station_data[st_point_list,]
      dis_v <- latlondistance(points_sub[i,"lat"], points_sub[i,"lon"], 
                              station_data_i$lat, station_data_i$lon)
      sto_state = rep(0,nrow(station_data_i))
      util <- (beta1*dis_v) + deltain[st_point_list]
      exputil <- exp(util)
      util_st <- (!sto_state) * exputil
      den_util <- sum(util_st)
      lambda_st_t <- rep(0,length(util_st))
      for(k in 1:length(v0_vec)) {
        out <- exp(-v0_vec[k]*sigma0)
        prob_t <- util_st/(out+den_util)*(points_sub[i,"density"]/length(v0_vec))
        lambda_st_t <- lambda_st_t + prob_t
      }
      lambda_t[list_obs] <- lambda_t[list_obs] + lambda_st_t
    }
  }
  return(lambda_t)
}


demand_distance_additive_change <- function(distance_change, theta1, wdcMerged, points,delta_list) {
  #find the number of people that substitute from a station to another
  #fix a station as focal station
  #select only points local to it. 
  #in first case select assume this station is stocked and find the demand at nearby stations
  #in next case take this station out from choice set and find demand and compare total demands
  #in each case
  
  no_st <- max(wdcMerged$station_id_index)
  tw_group_list <- unique(wdcMerged$tw_group)
  
  station_data <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index",
                                      "local_stations")])
  station_data <- station_data[order(station_data$station_id_index),]
  
  #   points_sub <- station_data$local_points[focal_station_index] 
  #points_sub <- which(stid_in_localstations(focal_station_index, points$local_stations))
  points_sub <- points
  
  
  
  beta1 <- theta1[1]
  sigma0 <- theta1[2]
  lat1 <- station_data[,"lat"]
  lon1 <- station_data[,"lon"]
  lambda_t <- rep(0,nrow(station_data)) 
  
  for(tw_groupin in tw_group_list) {
    deltain <- average_delta(wdcMerged, tw_groupin, delta_list)
#	deltain <- delta_list[which(wdcMerged$tw_group == tw_groupin)]
#   deltain <- as.numeric(by(deltain, wdcMerged$station_id_index[which(wdcMerged$tw_group == tw_groupin)], FUN=mean))
    for(i in 1:nrow(points_sub)) {
      lat2 <- points_sub[i,"lat"]
      lon2 <- points_sub[i,"lon"]
      #dis_v <- latlondistance(lat1, lon1, lat2, lon2)
      st_point_list <- splitchar(points_sub$local_stations[i])
      list_obs <- st_point_list
      if(length(list_obs)==0) next
      
      station_data_i <- station_data[st_point_list,]
      dis_v <- latlondistance(points_sub[i,"lat"], points_sub[i,"lon"], 
                              station_data_i$lat, station_data_i$lon)
      dis_v <- dis_v + distance_change
      sto_state = rep(0,nrow(station_data_i))
      #       if(mode==1) {
      #         sto_state[which(station_data_i$station_id_index==focal_station_index)]=1
      #       }
      util <- (beta1*dis_v) + deltain[st_point_list]
      exputil <- exp(util)
      util_st <- (!sto_state) * exputil
      den_util <- sum(util_st)
      lambda_st_t <- rep(0,length(util_st))
      for(k in 1:length(v0_vec)) {
        out <- exp(-v0_vec[k]*sigma0)
        prob_t <- util_st/(out+den_util)*(points_sub[i,"density"]/length(v0_vec))
        lambda_st_t <- lambda_st_t + prob_t
      }
      lambda_t[list_obs] <- lambda_t[list_obs] + lambda_st_t
    }
  }
  return(lambda_t)
}


average_demand_susbtitution <- function(theta1,wdcMerged,points, delta_list) {
  perc_substituted = 0
  for(focal_station_index in c(1:max(wdcMerged$station_id_index))) {
    print(focal_station_index)
    lambda1 <- demand_substitution (focal_station_index, 0, theta1, wdcMerged, points,delta_list)     
    lambda2 <- demand_substitution (focal_station_index, 1, theta1, wdcMerged, points,delta_list)     
    focal_st_demand <- lambda1[focal_station_index]
    substituted_demand <- sum(lambda2[-focal_station_index]) -sum(lambda1[-focal_station_index])
    perc_substituted <- perc_substituted + (substituted_demand/focal_st_demand*100)
    print(substituted_demand/focal_st_demand*100)
  }
  perc_substituted <- perc_substituted/max(wdcMerged$station_id_index)
  return(perc_substituted)
}

average_demand_susbtitution_demModel <- function(theta1,wdcMerged,points, delta_list) {
  perc_substituted = 0
  for(focal_station_index in c(1:max(wdcMerged$station_id_index))) {
    print(focal_station_index)
    lambda1 <- demand_substitution_demModel (focal_station_index, 0, theta1, wdcMerged, points,delta_list)     
    lambda2 <- demand_substitution_demModel (focal_station_index, 1, theta1, wdcMerged, points,delta_list)     
    focal_st_demand <- sum(lambda1[which(wdcMerged$station_id_index==focal_station_index)])    
    substituted_demand <- sum(lambda2) - sum(lambda1) + focal_st_demand
    perc_substituted <- perc_substituted + (substituted_demand/focal_st_demand*100)
    print(substituted_demand/focal_st_demand*100)
  }
  perc_substituted <- perc_substituted/max(wdcMerged$station_id_index)
  return(perc_substituted)
}


demand_serv_lvl_change_new <- function(theta1, wdcMerged,points, serv_lvl_factor) {
  res <- eval_error_xi(theta1,wdcMerged,points)
  deltain <- (res$X1 %*% res$theta2 + res$xi)/(sqrt(wdcMerged$obs_weight))
  demand_actual <- demand_distance_change_inner_new(1, theta1, deltain, wdcMerged, points)
  X1 <- res$X1
  X1[,"serv_lvl_covar"] <- X1[,"serv_lvl_covar"]*serv_lvl_factor
  #cap service level by 1
  X1[which(X1[,"serv_lvl_covar"]>1),"serv_lvl_covar"] <- 1
  deltain_new <- (X1 %*% res$theta2 + res$xi)/(sqrt(wdcMerged$obs_weight))
  demand_new <- demand_distance_change_inner_new(1, theta1, deltain_new, wdcMerged, points)
  return(c(demand_new,demand_actual))
}

demand_serv_lvl_change_new_systemdemand <- function(theta1, wdcMerged,points, serv_lvl_factor,serv_lvl_coef) {
  if(!exists("prev_theta") || !identical(prev_theta,theta1)) {
    deltain <- compute_delta_list_cntrt_map_new(theta1,wdcMerged,points)
    prev_theta <<- theta1
    prev_deltain <<- deltain
  } else {
    deltain <- prev_deltain
  }
  demand_actual <- demand_distance_change_inner_new(1, theta1, deltain, wdcMerged, points)
  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  no_st <- max(wdcMerged$station_id_index)
  i=1
  serv_lvl_covar <- c()
  for(tw_groupin in tw_groupin_list) {
    wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
    tw_index <- which(tw_groupin_list==tw_groupin)
    deltain_tw <- deltain[which(wdcMerged$tw_group==tw_groupin)]
    low = (i-1)*no_st + 1
    high = i*no_st
    serv_lvl = user_serv_lvl$serv_lvl[c(low:high)]
    serv_lvl_covar <- c(serv_lvl_covar, serv_lvl[wdcMergedday$station_id_index])
    i= i+1
  }
  
  serv_lvl_covar_new <- serv_lvl_covar*serv_lvl_factor
  serv_lvl_covar_new[which(serv_lvl_covar_new>1)] <- 1
  serv_lvl_covar_new[which(serv_lvl_covar_new<0)] <- 0
  deltain_new <- (serv_lvl_covar_new-serv_lvl_covar)*serv_lvl_coef + deltain
	
  demand_new <- demand_distance_change_inner_new(1, theta1, deltain_new, wdcMerged, points)
  return(c(demand_new,demand_actual))
}


demand_serv_lvl_change_new2 <- function(theta1, wdcMerged,points, serv_lvl_factor,serv_lvl_coef) {
  if(!exists("prev_theta") || !identical(prev_theta,theta1)) {
    deltain <- compute_delta_list_cntrt_map_new(theta1,wdcMerged,points)
    prev_theta <<- theta1
    prev_deltain <<- deltain
  } else {
    deltain <- prev_deltain
  }
  demand_actual <- demand_distance_change_inner_new(1, theta1, deltain, wdcMerged, points)
  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  no_st <- max(wdcMerged$station_id_index)
  i=1
  serv_lvl_covar <- c()
  for(tw_groupin in tw_groupin_list) {
    wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
    tw_index <- which(tw_groupin_list==tw_groupin)
    deltain_tw <- deltain[which(wdcMerged$tw_group==tw_groupin)]
    low = (i-1)*no_st + 1
    high = i*no_st
    serv_lvl = user_serv_lvl$serv_lvl[c(low:high)]
    serv_lvl_covar <- c(serv_lvl_covar, serv_lvl[wdcMergedday$station_id_index])
    i= i+1
  }
  
  serv_lvl_covar_new <- serv_lvl_covar*serv_lvl_factor
  serv_lvl_covar_new[which(serv_lvl_covar_new>1)] <- 1
  serv_lvl_covar_new[which(serv_lvl_covar_new<0)] <- 0
  deltain_new <- (serv_lvl_covar_new-serv_lvl_covar)*serv_lvl_coef + deltain
	
  tiny	<- 1e-4
  demand_new <- demand_distance_change_inner_new_serv_lvl_adjusted(1, theta1, deltain_new, wdcMerged, points, 
                                                               serv_lvl_covar_new+tiny, serv_lvl_covar+tiny)
  return(c(demand_new,demand_actual))
}

demand_serv_lvl_change_new2_tract <- function(theta1, wdcMerged,points, serv_lvl_factor,serv_lvl_coef) {
  if(!exists("prev_theta") || !identical(prev_theta,theta1)) {
    deltain <- compute_delta_list_cntrt_map_new(theta1,wdcMerged,points)
    prev_theta <<- theta1
    prev_deltain <<- deltain
  } else {
    deltain <- prev_deltain
  }
  #demand_actual <- demand_distance_change_inner_new(1, theta1, deltain, wdcMerged, points)
  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  no_st <- max(wdcMerged$station_id_index)
  i=1
  serv_lvl_covar <- c()
  for(tw_groupin in tw_groupin_list) {
    wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
    tw_index <- which(tw_groupin_list==tw_groupin)
    deltain_tw <- deltain[which(wdcMerged$tw_group==tw_groupin)]
    low = (i-1)*no_st + 1
    high = i*no_st
    serv_lvl = user_serv_lvl$serv_lvl[c(low:high)]
    serv_lvl_covar <- c(serv_lvl_covar, serv_lvl[wdcMergedday$station_id_index])
    i= i+1
  }

  serv_lvl_covar_new <- serv_lvl_covar*serv_lvl_factor
  serv_lvl_covar_new[which(serv_lvl_covar_new>1)] <- 1
  serv_lvl_covar_new[which(serv_lvl_covar_new<0)] <- 0
  deltain_new <- (serv_lvl_covar_new-serv_lvl_covar)*serv_lvl_coef + deltain
  
  demand_new <- demand_distance_change_inner_new_nls_tract_serv_lvl_adjusted(1, theta1, deltain_new, wdcMerged, points,serv_lvl_covar_new, serv_lvl_covar)
  return(c(demand_new))
}

demand_dis_serv_lvl_change_new_tract <- function(scale_distance_factor, theta1, wdcMerged,points, serv_lvl_factor,serv_lvl_coef) {
  if(!exists("prev_theta") || !identical(prev_theta,theta1)) {
    deltain <- compute_delta_list_cntrt_map_new(theta1,wdcMerged,points)
    prev_theta <<- theta1
    prev_deltain <<- deltain
  } else {
    deltain <- prev_deltain
  }
  #demand_actual <- demand_distance_change_inner_new(1, theta1, deltain, wdcMerged, points)
  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  no_st <- max(wdcMerged$station_id_index)
  i=1
  serv_lvl_covar <- c()
  for(tw_groupin in tw_groupin_list) {
    wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
    tw_index <- which(tw_groupin_list==tw_groupin)
    deltain_tw <- deltain[which(wdcMerged$tw_group==tw_groupin)]
    low = (i-1)*no_st + 1
    high = i*no_st
    serv_lvl = user_serv_lvl$serv_lvl[c(low:high)]
    serv_lvl_covar <- c(serv_lvl_covar, serv_lvl[wdcMergedday$station_id_index])
    i= i+1
  }
  
  serv_lvl_covar_new <- serv_lvl_covar*serv_lvl_factor
  serv_lvl_covar_new[which(serv_lvl_covar_new>1)] <- 1
  serv_lvl_covar_new[which(serv_lvl_covar_new<0)] <- 0
  deltain_new <- (serv_lvl_covar_new-serv_lvl_covar)*serv_lvl_coef + deltain
  
  demand_new <- demand_distance_change_inner_new_nls_tract_serv_lvl_adjusted(scale_distance_factor,
        theta1, deltain_new, wdcMerged, points,serv_lvl_covar_new, serv_lvl_covar)
  return(c(demand_new))
}

demand_serv_lvl_change_new_stfixedeffects <- function(theta1, wdcMerged,points, serv_lvl_factor) {
  res <- eval_error_xi(theta1,wdcMerged,points)
  #construct delta from res, so that you can scale service level coefficient and recompute
  serv_lvl_covar_mat <- cbind(1,user_serv_lvl$serv_lvl)
  serv_lvl_error = res$theta2 - serv_lvl_covar_mat%*%res$serv_lvl_coef
  deltain <- (res$X1 %*% res$theta2 + res$xi)/(sqrt(wdcMerged$obs_weight))
  demand_actual <- demand_distance_change_inner_new(1, theta1, deltain, wdcMerged, points)
  
  serv_lvl_covar_mat[,2] <- serv_lvl_covar_mat[,2]*serv_lvl_factor
  serv_lvl_covar_mat[which(serv_lvl_covar_mat[,2]>1),2] <- 1
  theta2_cf <- serv_lvl_covar_mat%*%res$serv_lvl_coef + serv_lvl_error
  deltain_new <- (res$X1 %*% theta2_cf + res$xi)/(sqrt(wdcMerged$obs_weight))
  demand_new <- demand_distance_change_inner_new(1, theta1, deltain_new, wdcMerged, points)
  return(c(demand_new,demand_actual))
}












demand_distance_change_inner_ind_new <- function(scale_distance_factor, theta1, delta_list, wdcMerged, points) {
  theta1[1] <- theta1[1]/scale_distance_factor
  tw_group_list <- c()
  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  
  lambda_t <- c()
  for(tw_groupin in tw_groupin_list) {
    wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
    deltain_tw <- delta_list[which(wdcMerged$tw_group==tw_groupin)]
    lambda_t <- c(lambda_t,
                  eval_lambda_delta_list_new(deltain_tw, theta1,wdcMergedday, points, tw_groupin)$objective)      
  }

  return((lambda_t))
}

demand_distance_change_ind_new <- function(scale_distance_factor, theta1, delta_list, wdcMerged, points) {
  demand_actual <- demand_distance_change_inner_ind_new(1, theta1, delta_list, wdcMerged, points)
  demand_changed <- demand_distance_change_inner_ind_new(scale_distance_factor, theta1, delta_list, wdcMerged, points)
  return(list(demand_changed,demand_actual))
}  

average_demand_susbtitution_ind_new <- function(theta1,deltain, wdcMerged,points) {
  #combine deltain by station x tw and then apply the aggregate method.
  wdcMerged$st_tw <- as.factor(wdcMerged$station_id_index*100 + wdcMerged$tw_group)
  deltain_avg = as.numeric(by(deltain, wdcMerged$st_tw, FUN=mean))
  perc_substituted = 0
  substituted_demand_vec <- c()
  for(focal_station_index in c(1:max(wdcMerged$station_id_index))) {
    print(focal_station_index)
    lambda1 <- demand_substitution_new (focal_station_index, 0, theta1, deltain_avg, wdcMerged, points)     
    lambda2 <- demand_substitution_new (focal_station_index, 1, theta1, deltain_avg, wdcMerged, points)     
    focal_st_demand <- lambda1[focal_station_index]
    substituted_demand <- sum(lambda2[-focal_station_index]) -sum(lambda1[-focal_station_index])
    perc_substituted <- perc_substituted + (substituted_demand/focal_st_demand*100)
    print(substituted_demand/focal_st_demand*100)
    substituted_demand_vec <- c(substituted_demand_vec,substituted_demand/focal_st_demand*100)
  }
  return(substituted_demand_vec)
}


demand_serv_lvl_change_ind_new <- function(theta1, wdcMerged,points, serv_lvl_factor) {
  res <- eval_error_xi(theta1,wdcMerged,points)
  if(exists("weighted_GMM")) {
    deltain <- res$X1 %*% res$theta2 + res$xi
    deltain <- deltain/sqrt(wdcMerged$obs_weight)
    X1 <- res$X1
    X1 <- X1/sqrt(wdcMerged$obs_weight)
    res$xi <- res$xi/sqrt(wdcMerged$obs_weight)
  } else {
    deltain <- res$X1 %*% res$theta2 + res$xi
    X1 <- res$X1
    X1[,"serv_lvl_covar"] <- X1[,"serv_lvl_covar"]*serv_lvl_factor
    #cap service level by 1
    X1[which(X1[,"serv_lvl_covar"]>1),"serv_lvl_covar"] <- 1
  }
  X1[,"serv_lvl_covar"] <- X1[,"serv_lvl_covar"]*serv_lvl_factor
  #cap service level by 1
  X1[which(X1[,"serv_lvl_covar"]>1),"serv_lvl_covar"] <- 1
  deltain_new <- X1 %*% res$theta2 + res$xi    
  
  demand_actual <- demand_distance_change_inner_ind_new(1, theta1, deltain, wdcMerged, points)
  demand_new <- demand_distance_change_inner_ind_new(1, theta1, deltain_new, wdcMerged, points)
  return(list(demand_new,demand_actual))
}


demand_serv_lvl_change_ind_new_presentation <- function(theta1, theta2, wdcMerged,points, serv_lvl_factor) {
  res <- eval_error_xi(theta1,wdcMerged,points)
  if(!exists("prev_theta") || !identical(prev_theta,theta1)) {
    deltain <- compute_delta_list_cntrt_map_new(theta1,wdcMerged,points)
    prev_theta <<- theta1
    prev_deltain <<- deltain
  } else {
    deltain <- prev_deltain
  }  
  if(exists("weighted_GMM")) {
    #from res, only thing to use is X1, have to recompute xi.
    #scale deltain
    deltain <- deltain*sqrt(wdcMerged$obs_weight)
    res$xi <- deltain - res$X1 %*% theta2
    deltain <- deltain/sqrt(wdcMerged$obs_weight)
    X1 <- res$X1
    X1 <- X1/sqrt(wdcMerged$obs_weight)
    res$xi <- res$xi/sqrt(wdcMerged$obs_weight)
  } else {
    stop("correct this part")
    deltain <- res$X1 %*% theta2 + res$xi
    X1 <- res$X1
    X1[,"serv_lvl_covar"] <- X1[,"serv_lvl_covar"]*serv_lvl_factor
    #cap service level by 1
    X1[which(X1[,"serv_lvl_covar"]>1),"serv_lvl_covar"] <- 1
  }
  X1[,"serv_lvl_covar"] <- X1[,"serv_lvl_covar"]*serv_lvl_factor
  #cap service level by 1
  X1[which(X1[,"serv_lvl_covar"]>1),"serv_lvl_covar"] <- 1
  deltain_new <- X1 %*% theta2 + res$xi    
  
  demand_actual <- demand_distance_change_inner_ind_new(1, theta1, deltain, wdcMerged, points)
  demand_new <- demand_distance_change_inner_ind_new(1, theta1, deltain_new, wdcMerged, points)
  return(list(demand_new,demand_actual))
}


demand_distance_transform_new <- function(theta1, delta_list, wdcMerged, points, dis_change) {
  #find the number of people that substitute from a station to another
  #fix a station as focal station
  #select only points local to it. 
  #in first case select assume this station is stocked and find the demand at nearby stations
  #in next case take this station out from choice set and find demand and compare total demands
  #in each case
  
  no_st <- max(wdcMerged$station_id_index)
  tw_group_list <- unique(wdcMerged$tw_group)
  
  station_data <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index",
                                      "local_stations")])
  #     #compute local points to stations
  #for each station, find out points that are in locality and store them as a string
  station_data$local_points <- rep(NA,nrow(station_data))
  for(i in 1:nrow(station_data)) {
    lat1 = station_data$lat[i]
    lon1 = station_data$lon[i]
    dis_v <- latlondistance(lat1,lon1,points$lat,points$lon)
    station_data$local_points[i] = paste(which(dis_v <= max_walking_dis), collapse="_")
  }
  station_data <- station_data[,c("station_id_index","local_points")]
  wdcMerged <- merge(wdcMerged,station_data,by="station_id_index")
  
  station_data <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index",
                                      "local_stations","local_points")])
  station_data <- station_data[order(station_data$station_id_index),]
  
  points_sub <- points
  
  beta1 <- theta1[1]
  sigma0 <- theta1[2]
  lat1 <- station_data[,"lat"]
  lon1 <- station_data[,"lon"]
  lambda_t <- rep(0,nrow(station_data)) 
  
  for(tw_i in 1:length(tw_group_list)) {    
    tw_groupin <- tw_group_list[tw_i]
    deltain <- delta_list[((tw_i-1)*no_st+1):(tw_i*no_st)]
    
    for(i in 1:nrow(points_sub)) {
      lat2 <- points_sub[i,"lat"]
      lon2 <- points_sub[i,"lon"]
      #dis_v <- latlondistance(lat1, lon1, lat2, lon2)
      #stop("chage this to station ids in st_opint_list")
      #list_obs <- which(dis_v<=max_walking_dis )
      st_point_list <- as.numeric(strsplit(points_sub$local_stations[i], split="_")[[1]])
      list_obs = st_point_list
      if(length(list_obs)==0) next
      station_data_i <- station_data[st_point_list,]
      dis_v <- latlondistance(points_sub[i,"lat"], points_sub[i,"lon"], 
                              station_data_i$lat, station_data_i$lon)
		dis_v <- dis_v + dis_change
      sto_state = rep(0,nrow(station_data_i))      
      util <- (beta1*dis_v) + deltain[st_point_list]
      exputil <- exp(util)
      util_st <- (!sto_state) * exputil
      den_util <- sum(util_st)
      lambda_st_t <- rep(0,length(util_st))
      for(k in 1:length(v0_vec)) {
        out <- exp(-v0_vec[k]*sigma0)
        prob_t <- util_st/(out+den_util)*(points_sub[i,"density"]/length(v0_vec))
        lambda_st_t <- lambda_st_t + prob_t
      }
      lambda_t[list_obs] <- lambda_t[list_obs] + lambda_st_t
    }
  }
  return(lambda_t)
}

demand_distance_change_new_linear_serv_lvl_adjusted_tract <- function(scale_distance_factor, theta1, 
                                                          wdcMerged,points, serv_lvl_factor,serv_lvl_coef) {
  if(!exists("prev_theta") || !identical(prev_theta,theta1)) {
    deltain <- compute_delta_list_cntrt_map_new(theta1,wdcMerged,points)
    prev_theta <<- theta1
    prev_deltain <<- deltain
  } else {
    deltain <- prev_deltain
  }
  #demand_actual <- demand_distance_change_inner_new(1, theta1, deltain, wdcMerged, points)
  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  no_st <- max(wdcMerged$station_id_index)
  i=1
  serv_lvl_covar <- c()
  for(tw_groupin in tw_groupin_list) {
    wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
    tw_index <- which(tw_groupin_list==tw_groupin)
    deltain_tw <- deltain[which(wdcMerged$tw_group==tw_groupin)]
    low = (i-1)*no_st + 1
    high = i*no_st
    serv_lvl = user_serv_lvl$serv_lvl[c(low:high)]
    serv_lvl_covar <- c(serv_lvl_covar, serv_lvl[wdcMergedday$station_id_index])
    i= i+1
  }
  
  serv_lvl_covar_new <- serv_lvl_covar*serv_lvl_factor
  serv_lvl_covar_new[which(serv_lvl_covar_new>1)] <- 1
  serv_lvl_covar_new[which(serv_lvl_covar_new<0)] <- 0
  deltain_new <- (serv_lvl_covar_new-serv_lvl_covar)*serv_lvl_coef + deltain
  
  tiny <- 1e-4
  demand_new <- demand_distance_change_inner_new_linear_serv_lvl_adjusted_tract(scale_distance_factor, 
      theta1, deltain_new, wdcMerged, points, serv_lvl_covar_new+tiny, serv_lvl_covar+tiny)
  return(demand_new)
  
  
}


demand_distance_change_inner_new_linear_serv_lvl_adjusted_tract <- function(scale_distance_factor, theta1, deltain_new, wdcMerged, points, 
                                                                            serv_lvl_covar_new, serv_lvl_covar) {
  theta1[1] <- theta1[1]/scale_distance_factor
  tw_group_list <- c()
  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  
  lambda_t <- c()
  for(tw_groupin in tw_groupin_list) {
    wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
    deltain_tw <- deltain_new[which(wdcMerged$tw_group==tw_groupin)]
    lambda_t <- c(lambda_t,
                  eval_lambda_delta_list_new(deltain_tw, theta1,wdcMergedday, points, tw_groupin)$objective)      
  }
  return(c(by(lambda_t*serv_lvl_covar_new/serv_lvl_covar, wdcMerged$tract, FUN=sum)))
}

average_delta <- function(wdcMerged, tw_groupin, delta_list) {

  deltain <- exp(delta_list[which(wdcMerged$tw_group == tw_groupin)])
  deltain <- deltain*(wdcMerged$obs_weight[which(wdcMerged$tw_group == tw_groupin)])
  wt_sum <- ave(wdcMerged$obs_weight[which(wdcMerged$tw_group == tw_groupin)],
                wdcMerged$station_id_index[which(wdcMerged$tw_group == tw_groupin)], FUN=sum)
  deltain_sum <- ave(deltain,
                     wdcMerged$station_id_index[which(wdcMerged$tw_group == tw_groupin)], FUN=sum)
  deltain <- deltain_sum[!duplicated(wdcMerged$station_id_index[which(wdcMerged$tw_group == tw_groupin)])]/
    wt_sum[!duplicated(wdcMerged$station_id_index[which(wdcMerged$tw_group == tw_groupin)])]
  deltain <- log(deltain)
	return(deltain)
}

demand_lost_st <- function(theta1,wdcMerged,points, delta_list) {
  lost_dem_vec <- c()
  st_actual_dem <- c()
  for(focal_station_index in c(1:max(wdcMerged$station_id_index))) {
    print(focal_station_index)
    lambda1 <- demand_substitution_demModel (focal_station_index, 0, theta1, wdcMerged, points,delta_list)     
    lambda2 <- demand_substitution_demModel (focal_station_index, 1, theta1, wdcMerged, points,delta_list)     
    
    focal_st_demand <- sum(lambda1[which(wdcMerged$station_id_index==focal_station_index)])    
    #substituted_demand <- sum(lambda2) - sum(lambda1) + focal_st_demand
    
    lost_dem_vec <- c(lost_dem_vec,sum(lambda2) - sum(lambda1))
    st_actual_dem <- c(st_actual_dem, focal_st_demand)
    
  }
  return(list(lost_dem_vec, st_actual_dem))
}


get_base_demand <- function(displace_serv_lvl, wdcMerged,points, theta2) {
  serv_lvl_coef <- theta2[which(rownames(theta2)=="serv_lvl_covar")]
  
  if(!exists("prev_theta") || !identical(prev_theta,theta1)) {
    deltain <- compute_delta_list_cntrt_map_new(theta1,wdcMerged,points)
    prev_theta <<- theta1
    prev_deltain <<- deltain
  } else {
    deltain <- prev_deltain
  }
  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  no_st <- max(wdcMerged$station_id_index)
  i=1
  serv_lvl_covar <- c()
  for(tw_groupin in tw_groupin_list) {
    wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
    tw_index <- which(tw_groupin_list==tw_groupin)
    deltain_tw <- deltain[which(wdcMerged$tw_group==tw_groupin)]
    low = (i-1)*no_st + 1
    high = i*no_st
    serv_lvl = user_serv_lvl$serv_lvl[c(low:high)]
    serv_lvl_covar <- c(serv_lvl_covar, serv_lvl[wdcMergedday$station_id_index])
    i= i+1
  }
  serv_lvl_covar_new <- serv_lvl_covar + displace_serv_lvl
  serv_lvl_covar_new[which(serv_lvl_covar_new > 1)] <- 1
  serv_lvl_covar_new[which(serv_lvl_covar_new < 0)] <- 0
  deltain_new <- serv_lvl_covar_new*serv_lvl_coef + deltain - serv_lvl_covar*serv_lvl_coef
  tiny <- 1e-4
  demand_actual <- demand_distance_change_inner_new_serv_lvl_adjusted(1, theta1, deltain_new, wdcMerged, points, serv_lvl_covar_new+tiny, serv_lvl_covar+tiny)
  return(demand_actual)
  
}  


optimal_distance_scaling <- function(displace_serv_lvl, wdcMerged,points, theta2, demand_actual) {
  
  serv_lvl_coef <- theta2[which(rownames(theta2)=="serv_lvl_covar")]
  
  if(!exists("prev_theta") || !identical(prev_theta,theta1)) {
    deltain <- compute_delta_list_cntrt_map_new(theta1,wdcMerged,points)
    prev_theta <<- theta1
    prev_deltain <<- deltain
  } else {
    deltain <- prev_deltain
  }
  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  no_st <- max(wdcMerged$station_id_index)
  i=1
  serv_lvl_covar <- c()
  for(tw_groupin in tw_groupin_list) {
    wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
    tw_index <- which(tw_groupin_list==tw_groupin)
    deltain_tw <- deltain[which(wdcMerged$tw_group==tw_groupin)]
    low = (i-1)*no_st + 1
    high = i*no_st
    serv_lvl = user_serv_lvl$serv_lvl[c(low:high)]
    serv_lvl_covar <- c(serv_lvl_covar, serv_lvl[wdcMergedday$station_id_index])
    i= i+1
  }
  serv_lvl_covar_new <- serv_lvl_covar + displace_serv_lvl
  serv_lvl_covar_new[which(serv_lvl_covar_new > 1)] <- 1
  serv_lvl_covar_new[which(serv_lvl_covar_new < 0)] <- 0
  deltain_new <- serv_lvl_covar_new*serv_lvl_coef + deltain - serv_lvl_covar*serv_lvl_coef
  
  
  ftol <- 1
  tol <- 0.01
  
  a <- 0.05
  b <- 16
  N <- 1
  NMAX <- 100
  tiny <- 1e-4
  
  # func_a <- demand_distance_change_inner_new(1/a, theta1, deltain_new, wdcMerged, points)
  # func_b <- demand_distance_change_inner_new(1/b, theta1, deltain_new, wdcMerged, points)
  func_a <- demand_distance_change_inner_new_serv_lvl_adjusted(1/a, theta1, deltain_new, wdcMerged, points,
                                                               serv_lvl_covar_new+tiny, serv_lvl_covar+tiny)  
  func_b <- demand_distance_change_inner_new_serv_lvl_adjusted(1/b, theta1, deltain_new, wdcMerged, points, 
                                                               serv_lvl_covar_new+tiny, serv_lvl_covar+tiny)
  
  while(N < NMAX) {
    c <- (a+b)/2  
    print(c)
    func_c <- demand_distance_change_inner_new_serv_lvl_adjusted(1/c, theta1, deltain_new, wdcMerged, points, 
                                                                 serv_lvl_covar_new+tiny, serv_lvl_covar+tiny)
    if(!(func_c <= func_a & func_c >= func_b)) {
      stop("values are not in order")
    }	
    
    if(abs(func_c - demand_actual) < ftol | (b-a)/2 < tol) {
      print(c)
      break
    }
    N <- N + 1
    if(func_c > demand_actual) {
      a <- c
    } else {
      b <- c
    }
  }
  return(c)
}

distance_segment_demand <- function(theta1, wdcMerged, points, delta_list) {
  #compute the fraction of demand that goes to stations that is 10mts away, .. in buckets of 10mts
  #compute fraction of demand that goes to closest station
  #compute fraction of demand that goes to highest service level station
  
  dem_seg <- rep(0,max_walking_dis/0.01)
  outside_dem <- 0
  dem_closest_st <- 0
  dem_highserv_st <- 0
  density_mass <- 0
  no_st <- max(wdcMerged$station_id_index)
  tw_group_list <- unique(wdcMerged$tw_group)
  
  station_data <- unique(wdcMerged[,c("station_id","lat","lon","station_id_index",
                                      "local_stations")])
  station_data <- station_data[order(station_data$station_id_index),]
  
  #   points_sub <- station_data$local_points[focal_station_index] 
  #   points_sub <- which(stid_in_localstations(focal_station_index, points$local_stations))
  #   points_sub <- points[points_sub,]
  beta1 <- theta1[1]
  sigma0 <- theta1[2]
  lat1 <- station_data[,"lat"]
  lon1 <- station_data[,"lon"]
  lambda_t <- rep(0,nrow(station_data)) 
  
  for(tw_groupin in tw_group_list) {
    deltain <- exp(delta_list[which(wdcMerged$tw_group == tw_groupin)])
    deltain <- deltain*(wdcMerged$obs_weight[which(wdcMerged$tw_group == tw_groupin)])
    wt_sum <- ave(wdcMerged$obs_weight[which(wdcMerged$tw_group == tw_groupin)],
                  wdcMerged$station_id_index[which(wdcMerged$tw_group == tw_groupin)], FUN=sum)
    deltain_sum <- ave(deltain,
                       wdcMerged$station_id_index[which(wdcMerged$tw_group == tw_groupin)], FUN=sum)
    deltain <- deltain_sum[!duplicated(wdcMerged$station_id_index[which(wdcMerged$tw_group == tw_groupin)])]/
      wt_sum[!duplicated(wdcMerged$station_id_index[which(wdcMerged$tw_group == tw_groupin)])]
    deltain <- log(deltain)
    #     deltain <- delta_list[which(wdcMerged$tw_group == tw_groupin)]
    #     deltain_sum <- ave(deltain,
    #                        wdcMerged$station_id_index[which(wdcMerged$tw_group == tw_groupin)], FUN=max)
    #     deltain <- deltain_sum[!duplicated(wdcMerged$station_id_index[which(wdcMerged$tw_group == tw_groupin)])]
    
    for(i in 1:nrow(points)) {
      lat2 <- points[i,"lat"]
      lon2 <- points[i,"lon"]
      #dis_v <- latlondistance(lat1, lon1, lat2, lon2)
      st_point_list <- splitchar(points$local_stations[i])
      list_obs <- st_point_list
      if(length(list_obs)==0) next
      
      station_data_i <- station_data[st_point_list,]
      dis_v <- latlondistance(points[i,"lat"], points[i,"lon"], 
                              station_data_i$lat, station_data_i$lon)
      sto_state = rep(0,nrow(station_data_i))
      util <- (beta1*dis_v) + deltain[st_point_list]
      exputil <- exp(util)
      util_st <- (!sto_state) * exputil
      den_util <- sum(util_st)
      lambda_st_t <- rep(0,length(util_st))
      #       print("")
      #       print(st_point_list)
      #       print(dis_v)
      #       print(util_st/(1+den_util))
      
      for(k in 1:length(v0_vec)) {
        out <- exp(-v0_vec[k]*sigma0)
        prob_t <- util_st/(out+den_util)*(points[i,"density"]/length(v0_vec))
        dem_seg[ceiling(dis_v/0.01)] <- dem_seg[ceiling(dis_v/0.01)] + prob_t
        lambda_st_t <- lambda_st_t + prob_t
        outside_dem <- outside_dem + (points[i,"density"]/length(v0_vec)) - sum(prob_t)
        dem_closest_st <- dem_closest_st + prob_t[which.min(dis_v)]
        
        density_mass <- density_mass + points[i,"density"]/length(v0_vec)
        st_serv_lvl <- user_serv_lvl$serv_lvl[which(user_serv_lvl$tw_group==tw_groupin)][st_point_list]
        dem_highserv_st <- dem_highserv_st + prob_t[which.max(st_serv_lvl)]
      }
      lambda_t[list_obs] <- lambda_t[list_obs] + lambda_st_t
    }
  }
  return(list(lambda_t,dem_seg, outside_dem,dem_closest_st,dem_highserv_st,density_mass))
}

