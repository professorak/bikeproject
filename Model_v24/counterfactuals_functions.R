##data generations
generate_cf1data <- function(wdcMerged, delta_list, user_serv_lvl) {
  st_tw_data <- wdcMerged[!duplicated(wdcMerged$st_tw),]  
  st_tw_data$dummy <- 1
  table_sto_states <- data.frame(sto_state_no = c(0:(2^(hyperlocal_range-1)-1)))
  #hyperlocal_range-1 above so that focal stations is always stocked in in all states  
  table_sto_states$hyperlocal_sto_state_vec <- 
    sapply(table_sto_states$sto_state_no, FUN=convert_sto_state_no_to_str) 
  table_sto_states$dummy <- 1
  table_sto_states <- table_sto_states[,c("hyperlocal_sto_state_vec","dummy")]
  wdcMerged_cfl1 <- merge(st_tw_data, table_sto_states, by="dummy", all=T)
  wdcMerged_cfl1$sto_state_local <- NULL
  wdcMerged_cfl1$sto_state_local_fac <- NULL
  
  
  st_lat_lon_df <- wdcMerged_cfl1[!duplicated(wdcMerged_cfl1$station_id_index),
                                  c("station_id_index","lat","lon")]
  st_lat_lon_df <- st_lat_lon_df[order(st_lat_lon_df$station_id_index),]
  st_lat_vec <<- st_lat_lon_df$lat
  st_lon_vec <<- st_lat_lon_df$lon
  
  #generate sto_state_local based on hyperlocal_sto_state_vec and length of local_stations
  #   #rearrgange local stations in order of dis_v
  ret_list <- t(apply(wdcMerged_cfl1[,c("station_id_index","local_stations","hyperlocal_sto_state_vec")], 
                                          1, FUN= gen_local_sto_state_cf1))
  wdcMerged_cfl1$sto_state_local <- ret_list[,1]
  wdcMerged_cfl1$hyperlocal_stations <- ret_list[,2]
  wdcMerged_cfl1$stocked_out == F
  
  #take weighted average of delta_list by st_tw
  delta_list_df <- data.frame(st_tw=wdcMerged$st_tw[stocked_list],
                              st_tw_index=wdcMerged$st_tw_index[stocked_list],
                              obs_weight= wdcMerged$obs_weight[stocked_list],delta_list=delta_list)
  delta_list_df$tot_obs_weight <- 
    ave(delta_list_df$obs_weight, delta_list_df$st_tw, FUN=sum)
  delta_list_df$tot_obs_weight_delta <- 
    ave(delta_list_df$obs_weight*delta_list_df$delta_list, delta_list_df$st_tw, FUN=sum)
  delta_list_df$weightedavg_delta <- delta_list_df$tot_obs_weight_delta / delta_list_df$tot_obs_weight
  delta_list_df <- delta_list_df[!duplicated(delta_list_df$st_tw),]
  delta_list_df <- delta_list_df[,c("st_tw","st_tw_index","weightedavg_delta")]
  delta_list_df <- delta_list_df[order(delta_list_df$st_tw_index),]
  #use this to make delta_list_cfl1
  delta_list_cfl1 <- delta_list_df$weightedavg_delta[wdcMerged_cfl1$st_tw_index]
  
  wdcMerged_cfl1$obs_weight <- generate_weights(wdcMerged_cfl1, user_serv_lvl, 1)
    
  return(list(wdcMerged_cfl1, delta_list_cfl1))
}

gen_local_sto_state_cf1 <- function(in_vec) {
  stidindex <- as.numeric(in_vec[1]) 
  local_station_str <- in_vec[2] 
  hyperlocal_sto_state_str <- in_vec[3] 
  local_st_list <- as.numeric(strsplit(local_station_str, split="_")[[1]])
  dis_v_local_list <- latlondistance(st_lat_vec[stidindex], st_lon_vec[stidindex],
                                     st_lat_vec[local_st_list], st_lon_vec[local_st_list])
  local_sto_state_vec <- rep(0, length(local_st_list))
  local_sto_state_vec[(order(dis_v_local_list)[1:hyperlocal_range])] <- 
    as.numeric(strsplit(hyperlocal_sto_state_str, split="_")[[1]])    
  hyperlocal_stations_vec <- sort(local_st_list[(order(dis_v_local_list)[1:hyperlocal_range])])  
  return(c(paste(local_sto_state_vec, collapse="_"),
           paste(hyperlocal_stations_vec, collapse="_")))
}



convert_sto_state_no_to_str <- function(sto_state_no) {
  sto_state_vec <- rep(0, hyperlocal_range)
  if(sto_state_no!=0) {
    bin_vec <- digitsBase(sto_state_no, base=2)    
    sto_state_vec[(hyperlocal_range-length(bin_vec)+1):hyperlocal_range] <- bin_vec    
  }
  return(paste(sto_state_vec, collapse="_"))
}


generate_weights <- function(wdcMerged_cfl1, user_serv_lvl, scale_factor) {
  #scale serv_lvl in user_serv_lvl
  user_serv_lvl$serv_lvl <- user_serv_lvl$serv_lvl*scale_factor
  user_serv_lvl$serv_lvl[which(user_serv_lvl$serv_lvl>1)] <- 1
  user_serv_lvl$serv_lvl[which(user_serv_lvl$serv_lvl<0)] <- 0
    
  #for each tw run seperately
  #run an apply function which assgin the product of probabilities as weights.
  tw_group_list <- unique(wdcMerged_cfl1$tw_group)
  weights_vec <- rep(NA, nrow(wdcMerged_cfl1))
  for(tw_group_in in tw_group_list) {    
    wdcMerged_cfl1_tw <- wdcMerged_cfl1[which(wdcMerged_cfl1$tw_group==tw_group_in),]
    tw_in <- wdcMerged_cfl1_tw$tw[1]
    user_serv_lvl_tw <<- user_serv_lvl[which(user_serv_lvl$tw==tw_in),]
    user_serv_lvl_tw <<- user_serv_lvl_tw[order(user_serv_lvl_tw$station_id_index),]
    weights_vec[which(wdcMerged_cfl1$tw_group==tw_group_in)] <- 
      apply(wdcMerged_cfl1_tw[,c("local_stations","sto_state_local","hyperlocal_stations")],
      1,FUN=gen_state_prob_cf1)
  }
  if(length(which(is.na(weights_vec)))>0) {
    stop("NA present in weights_vec")
  }
  return(weights_vec)
}

generate_weights_additive_servlvl <- function(wdcMerged_cfl1, user_serv_lvl, displace_serv_lvl) {
  #scale serv_lvl in user_serv_lvl
  user_serv_lvl$serv_lvl <- user_serv_lvl$serv_lvl + displace_serv_lvl
  user_serv_lvl$serv_lvl[which(user_serv_lvl$serv_lvl>1)] <- 1
  user_serv_lvl$serv_lvl[which(user_serv_lvl$serv_lvl<0)] <- 0
  
  #for each tw run seperately
  #run an apply function which assgin the product of probabilities as weights.
  tw_group_list <- unique(wdcMerged_cfl1$tw_group)
  weights_vec <- rep(NA, nrow(wdcMerged_cfl1))
  for(tw_group_in in tw_group_list) {    
    wdcMerged_cfl1_tw <- wdcMerged_cfl1[which(wdcMerged_cfl1$tw_group==tw_group_in),]
    tw_in <- wdcMerged_cfl1_tw$tw[1]
    user_serv_lvl_tw <<- user_serv_lvl[which(user_serv_lvl$tw==tw_in),]
    user_serv_lvl_tw <<- user_serv_lvl_tw[order(user_serv_lvl_tw$station_id_index),]
    weights_vec[which(wdcMerged_cfl1$tw_group==tw_group_in)] <- 
      apply(wdcMerged_cfl1_tw[,c("local_stations","sto_state_local","hyperlocal_stations")],
            1,FUN=gen_state_prob_cf1)
  }
  if(length(which(is.na(weights_vec)))>0) {
    stop("NA present in weights_vec")
  }
  return(weights_vec)
}


gen_state_prob_cf1 <- function(in_vec) {
  local_station_str <- in_vec[1]
  sto_state_local <- in_vec[2] 
  hyperlocal_stations <- in_vec[3] 
  
  local_st_list <- as.numeric(strsplit(local_station_str, split="_")[[1]])
  sto_state_local_vec <- as.numeric(strsplit(sto_state_local, split="_")[[1]])
  hyperlocal_stations_vec <- as.numeric(strsplit(hyperlocal_stations, split="_")[[1]])
  
  hyperlocal_index <- 
    which(local_st_list %in% hyperlocal_stations_vec)
  if(!identical(local_st_list[hyperlocal_index],hyperlocal_stations_vec)) stop("hyperlocal_index incorrect")

  sto_state_hyperlocal_vec <- sto_state_local_vec[hyperlocal_index]
    
  return(prod(user_serv_lvl_tw$serv_lvl[hyperlocal_stations_vec] * (1-sto_state_hyperlocal_vec) + 
           (1-user_serv_lvl_tw$serv_lvl[hyperlocal_stations_vec]) * sto_state_hyperlocal_vec))
}




#counterfactual results computation

demand_distance_change <- function(scale_distance_factor_vec, theta1, delta_list, wdcMerged, points) {
  demand_vec <- c()
  for(scale_distance_factor in scale_distance_factor_vec) {
    demand_vec <- c(demand_vec,
      demand_distance_change_inner(scale_distance_factor, theta1, delta_list, wdcMerged, points))    
  }
  return(demand_vec)
}  

demand_distance_change_inner <- function(scale_distance_factor, theta1, delta_list, wdcMerged, points) {
  theta1[1] <- theta1[1]/scale_distance_factor
  theta1[3] <- theta1[3]/scale_distance_factor^2
  
  lambda_t <- eval_lambda_full_unweighted(delta_list, theta1, wdcMerged, points)

  return(sum(lambda_t*wdcMerged$obs_weight))
}

demand_serv_lvl_change <- function(theta1, delta_list, wdcMerged,points, serv_lvl_factor_vec, serv_lvl_coef) {
  demand_vec <- c()
  for(serv_lvl_factor in serv_lvl_factor_vec) {
    demand_vec <- c(demand_vec,
      demand_serv_lvl_change_inner(theta1, delta_list, wdcMerged,points, serv_lvl_factor, serv_lvl_coef))    
  }
  return(demand_vec)
}

demand_serv_lvl_change_inner <- function(theta1, delta_list, wdcMerged,points, serv_lvl_factor, serv_lvl_coef) {
  #construct new delta for serv_lvl_factor (assumed wdcMerged$serv_lvl is proper,
  #can also altenatively construct from user_serv_lvl)
  serv_lvl_covar_new <- serv_lvl_factor*wdcMerged$serv_lvl
  serv_lvl_covar_new[which(serv_lvl_covar_new>1)] <- 1
  serv_lvl_covar_new[which(serv_lvl_covar_new<0)] <- 0
  
  stocked_list <- which(wdcMerged$stocked_out==F)
  
  delta_list_new <- delta_list + serv_lvl_coef[1]*((serv_lvl_covar_new-wdcMerged$serv_lvl)[stocked_list]) +
    serv_lvl_coef[2]*((serv_lvl_covar_new^2-wdcMerged$serv_lvl^2)[stocked_list])
  
  #construct new weights
  wdcMerged$obs_weight <- generate_weights(wdcMerged, user_serv_lvl, serv_lvl_factor)
  
  lambda_t <- eval_lambda_full_unweighted(delta_list_new, theta1, wdcMerged, points)
  
  return(sum(lambda_t*wdcMerged$obs_weight))
  
}

demand_serv_lvl_change_shortterm <- function(theta1, delta_list, wdcMerged,points, serv_lvl_factor_vec, serv_lvl_coef) {
  demand_vec <- c()
  for(serv_lvl_factor in serv_lvl_factor_vec) {
    demand_vec <- c(demand_vec,
      demand_serv_lvl_change_inner_shortterm(theta1, delta_list, wdcMerged,points, serv_lvl_factor, serv_lvl_coef))    
  }
  return(demand_vec)
}

demand_serv_lvl_change_inner_shortterm <- function(theta1, delta_list, wdcMerged,points, serv_lvl_factor, serv_lvl_coef) {
  #not constructing new delta for short term effect, as service level effects
  #will remain the same, just weights of states will change
  
  #construct new weights
  wdcMerged$obs_weight <- generate_weights(wdcMerged, user_serv_lvl, serv_lvl_factor)
  
  lambda_t <- eval_lambda_full_unweighted(delta_list, theta1, wdcMerged, points)
  
  return(sum(lambda_t*wdcMerged$obs_weight))  
}

demand_distance_change_serv_lvl_change <- function(scale_distance_factor_vec, theta1, delta_list, wdcMerged,points, serv_lvl_factor_vec, serv_lvl_coef) {
  demand_vec <- c()
  for(i in c(1:length(serv_lvl_factor_vec))) {
    demand_vec <- c(demand_vec,
                    demand_distance_change_serv_lvl_change_inner(scale_distance_factor_vec[i],
                                                                 theta1, delta_list, wdcMerged,points, serv_lvl_factor_vec[i], serv_lvl_coef))    
  }
  return(demand_vec)
}

demand_distance_change_serv_lvl_change_inner <- function(scale_distance_factor,
                                                         theta1, delta_list, wdcMerged,points, serv_lvl_factor, serv_lvl_coef) {
  #distance factor change
  theta1[1] <- theta1[1]/scale_distance_factor
  theta1[3] <- theta1[3]/scale_distance_factor^2
  
  #construct new delta for serv_lvl_factor (assumed wdcMerged$serv_lvl is proper,
  #can also altenatively construct from user_serv_lvl)
  serv_lvl_covar_new <- serv_lvl_factor*wdcMerged$serv_lvl
  serv_lvl_covar_new[which(serv_lvl_covar_new>1)] <- 1
  serv_lvl_covar_new[which(serv_lvl_covar_new<0)] <- 0
  
  stocked_list <- which(wdcMerged$stocked_out==F)
  
  delta_list_new <- delta_list + serv_lvl_coef[1]*((serv_lvl_covar_new-wdcMerged$serv_lvl)[stocked_list]) +
    serv_lvl_coef[2]*((serv_lvl_covar_new^2-wdcMerged$serv_lvl^2)[stocked_list])
  
  #construct new weights
  wdcMerged$obs_weight <- generate_weights(wdcMerged, user_serv_lvl, serv_lvl_factor)
  
  lambda_t <- eval_lambda_full_unweighted(delta_list_new, theta1, wdcMerged, points)
  
  return(sum(lambda_t*wdcMerged$obs_weight))
  
}

# demand_serv_lvl_change_new2 <- function(theta1, wdcMerged,points, serv_lvl_factor,serv_lvl_coef) {
#   if(!exists("prev_theta") || !identical(prev_theta,theta1)) {
#     deltain <- compute_delta_list_cntrt_map_new(theta1,wdcMerged,points)
#     prev_theta <<- theta1
#     prev_deltain <<- deltain
#   } else {
#     deltain <- prev_deltain
#   }
#   demand_actual <- demand_distance_change_inner_new(1, theta1, deltain, wdcMerged, points)
#   
#   tw_groupin_list <- unique(wdcMerged$tw_group)    
#   no_st <- max(wdcMerged$station_id_index)
#   i=1
#   serv_lvl_covar <- c()
#   for(tw_groupin in tw_groupin_list) {
#     wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
#     tw_index <- which(tw_groupin_list==tw_groupin)
#     deltain_tw <- deltain[which(wdcMerged$tw_group==tw_groupin)]
#     low = (i-1)*no_st + 1
#     high = i*no_st
#     serv_lvl = user_serv_lvl$serv_lvl[c(low:high)]
#     serv_lvl_covar <- c(serv_lvl_covar, serv_lvl[wdcMergedday$station_id_index])
#     i= i+1
#   }
#   
#   serv_lvl_covar_new <- serv_lvl_covar*serv_lvl_factor
#   serv_lvl_covar_new[which(serv_lvl_covar_new>1)] <- 1
#   serv_lvl_covar_new[which(serv_lvl_covar_new<0)] <- 0
#   deltain_new <- (serv_lvl_covar_new-serv_lvl_covar)*serv_lvl_coef + deltain
#   
#   tiny  <- 1e-4
#   demand_new <- demand_distance_change_inner_new_serv_lvl_adjusted(1, theta1, deltain_new, wdcMerged, points, 
#                                                                    serv_lvl_covar_new+tiny, serv_lvl_covar+tiny)
#   return(c(demand_new,demand_actual))
# }
# 
# demand_distance_change_inner_new_serv_lvl_adjusted <- function(scale_distance_factor, theta1, delta_list, wdcMerged, points, 
#                                                                serv_lvl_covar_new, serv_lvl_covar) {
#   theta1[1] <- theta1[1]/scale_distance_factor
#   tw_group_list <- c()
#   
#   tw_groupin_list <- unique(wdcMerged$tw_group)    
#   
#   lambda_t <- c()
#   for(tw_groupin in tw_groupin_list) {
#     wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
#     deltain_tw <- delta_list[which(wdcMerged$tw_group==tw_groupin)]
#     lambda_t <- c(lambda_t,
#                   eval_lambda_delta_list_new(deltain_tw, theta1,wdcMergedday, points, tw_groupin)$objective)      
#   }
#   return(sum(lambda_t*serv_lvl_covar_new/serv_lvl_covar))
# }
# 

###do simulation for an hour length of time and see the fraction of time station is stocked out. <- not anymore, doing long term.
#we assume no queuing dynamics, so a demand lost due to stockout is lost forever.
#at time 0, a station starts with number of bikes given by uniform discrete distribution between 0 and station_size.

gen_simulated_serv_lvl_longterm <- function(mean_st_size, mean_out_dem, mean_in_dem,scaling_fac) { 
  mean_st_size <- mean_st_size * scaling_fac
  mean_out_dem <- mean_out_dem * scaling_fac
  mean_in_dem <- mean_in_dem * scaling_fac
  
  tmp <- .Random.seed
  #set.seed(100)
  time_limit <- 10000
  st_size <- round(mean_st_size)
  stockout_frac_list <- c()
  stockout_threshold <- 0 * scaling_fac
  for(starting_bikes in sample(0:st_size,10,replace=T)) {
    for(sims in 1:10) {
      num_bikes <- starting_bikes
      stockout_frac <- 0
      out_dem_epochs <- cumsum(rexp(mean_out_dem*time_limit*10, rate=mean_out_dem))
      if(length(which(out_dem_epochs>time_limit))==0) stop("short")
      out_dem_epochs <- out_dem_epochs[out_dem_epochs<time_limit]
      in_dem_epochs <- cumsum(rexp(mean_in_dem*time_limit*10, rate=mean_in_dem))
      if(length(which(in_dem_epochs>time_limit))==0) stop("short")
      in_dem_epochs <- in_dem_epochs[in_dem_epochs<time_limit]
      out_dem_epochs_df <- data.frame(time_epochs=out_dem_epochs, out_dem=rep(1,length(out_dem_epochs)))
      in_dem_epochs_df <- data.frame(time_epochs=in_dem_epochs, out_dem=rep(0,length(in_dem_epochs)))
      out_dem_epochs_df <- rbind(out_dem_epochs_df,in_dem_epochs_df)
      
      out_dem_epochs_df <- out_dem_epochs_df[order(out_dem_epochs_df$time_epochs),]
      out_dem_epochs_df <- rbind(out_dem_epochs_df,c(time_limit,0))
      
      if(num_bikes==stockout_threshold) {
        stockout_frac <- out_dem_epochs_df$time_epochs[1]
      }
      if(nrow(out_dem_epochs_df)>1) {
        for(i in 1:(nrow(out_dem_epochs_df)-1)) {
          if(out_dem_epochs_df$out_dem[i]==1) {
            num_bikes <- max(num_bikes-1,0)
          } else {
            num_bikes <- min(num_bikes+1,st_size)
          }
          if(num_bikes<=stockout_threshold) {
            stockout_frac <- stockout_frac + out_dem_epochs_df$time_epochs[i+1]-out_dem_epochs_df$time_epochs[i]
          }      
        }
      }
      stockout_frac_list <- c(stockout_frac_list,stockout_frac/time_limit) 
    }
  }
  print(1-stockout_frac_list)
  .Random.seed <- tmp
  return(1-mean(stockout_frac_list))
}

get_average_station_characteristics <- function(wdcMerged) {
  wdcMerged$station_size <- wdcMerged$bikes + wdcMerged$spaces_available
  mean_st_size <- mean(as.numeric(by(wdcMerged$station_size, 
    wdcMerged$station_id_index, FUN=get_median_non_NA)))
  
  #for each st_tw compute mean_dem|stocked_in and then weight by serv_lvl
  out_dem_df <- data.frame(st_tw=wdcMerged$st_tw, st_tw_index=wdcMerged$st_tw_index)
  out_dem_df$tot_dem_sum <- ave(wdcMerged$out_dem_sum,wdcMerged$st_tw, FUN=sum)
  out_dem_df$tot_obs_weight <- ave(wdcMerged$obs_weight,wdcMerged$st_tw, FUN=sum)
  out_dem_df$mean_outdem_stockedin <- out_dem_df$tot_dem_sum/out_dem_df$tot_obs_weight
  out_dem_df <- out_dem_df[!duplicated(out_dem_df$st_tw),]
#   #join with user_serv_lvl
#   out_dem_df <- out_dem_df[order(out_dem_df$st_tw_index),]
#   identical(out_dem_df$st_tw, user_serv_lvl$st_tw)
#   out_dem_df <- cbind(out_dem_df, serv_lvl=user_serv_lvl$serv_lvl)
#   out_dem_df$mean_outdem <- out_dem_df$mean_outdem_stockedin * out_dem_df$serv_lvl
  
  mean_out_dem <- mean(out_dem_df$mean_outdem_stockedin)*30
  mean_in_dem <- mean_out_dem
  if(is.na(mean_st_size+mean_out_dem+mean_in_dem)) stop("NA value for either of 
      mean_st_size,mean_out_dem,mean_in_dem")
  return(c(mean_st_size, mean_out_dem, mean_in_dem))
}


# #compute mean station size, mean out dem and mean in dem
run_counterfactual_l2 <- function(mean_st_size, mean_out_dem, mean_in_dem, theta1, delta_list, wdcMerged,points, serv_lvl_coef) {
  
  beta_list <- c(8,4,2,1/2,1/4,1/8) #high beta is high station density
  #beta scales distances linearly and service level aatributes (station size, demand) 
  #quadratically and inversely.
  #beta_list <- c(3.5,3,2.5,1.5)
  
  target_st_size <- c(ceiling(mean_st_size*beta_list^-2))
  beta_list <- (target_st_size/mean_st_size)^(-1/2)
  
  base_serv_lvl <- gen_simulated_serv_lvl_longterm(mean_st_size,mean_out_dem,mean_in_dem,1)  
  #generate serv_lvl_factor_vec from beta_list
  serv_lvl_factor_vec <- c()
  scale_distance_factor_vec <- c()
  for(beta in beta_list) {
    changed_serv_lvl <- gen_simulated_serv_lvl_longterm(mean_st_size,mean_out_dem,mean_in_dem,beta^-2)    
    serv_lvl_factor_vec <- c(serv_lvl_factor_vec, changed_serv_lvl/base_serv_lvl)  
    scale_distance_factor_vec <- c(scale_distance_factor_vec,beta)
  }
  serv_lvl_factor_vec <- c(1, serv_lvl_factor_vec)
  scale_distance_factor_vec <- c(1, scale_distance_factor_vec)
  beta_list <- c(1, beta_list)

  demand_vec <- demand_distance_change_serv_lvl_change(scale_distance_factor_vec, 
    theta1, delta_list, wdcMerged,points, serv_lvl_factor_vec, serv_lvl_coef)
    
  return(list("demand_vec"=demand_vec[order(beta_list)],
              "beta_list"=beta_list[order(beta_list)],
              "serv_lvl_factor_vec"=serv_lvl_factor_vec[order(beta_list)],
              "scale_distance_factor_vec"=scale_distance_factor_vec[order(beta_list)],
              "base_serv_lvl"=base_serv_lvl))
}



#compute average demand substitution
demand_substitution_demModel <- function(focal_station_index, mode, theta1, wdcMerged_sttw_ave, points,delta_list_sttw_ave) {
    
  points_sub <- which(stid_in_localstations(focal_station_index, points$local_stations))
  points_sub <- points[points_sub,]    
  
  if(mode==1) {
    delta_list_sttw_ave[which(wdcMerged_sttw_ave$station_id_index==focal_station_index)] <- -100
  }
  
  lambda_t <- eval_lambda_full_unweighted(delta_list_sttw_ave, theta1, wdcMerged_sttw_ave, points_sub)
  return((lambda_t))
}


average_demand_susbtitution_demModel <- function(theta1, wdcMerged_sttw_ave, points,delta_list_sttw_ave) {
  perc_substituted <- c()
  for(focal_station_index in c(1:max(wdcMerged_sttw_ave$station_id_index))) {
    print(focal_station_index)
    lambda1 <- demand_substitution_demModel (focal_station_index, 0, theta1, wdcMerged_sttw_ave, points,delta_list_sttw_ave)     
    lambda2 <- demand_substitution_demModel (focal_station_index, 1, theta1, wdcMerged_sttw_ave, points,delta_list_sttw_ave)     
    focal_st_demand <- sum(lambda1[which(wdcMerged_sttw_ave$station_id_index==focal_station_index)])    
    substituted_demand <- sum(lambda2) - sum(lambda1) + focal_st_demand
    perc_substituted <- c(perc_substituted, (substituted_demand/focal_st_demand*100))
    print(substituted_demand/focal_st_demand*100)
  }
  return(perc_substituted)
}

#pass wdcMerged_sttw_ave and delta_list_sttw_ave to this function
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
  
  beta1 <- theta1[1]
  sigma0 <- theta1[2]
  beta2 <- theta1[3]
  
  lat1 <- wdcMerged[,"lat"]
  lon1 <- wdcMerged[,"lon"]
  lambda_t <- rep(0,nrow(wdcMerged)) 
  
  tw_group_list <- unique(wdcMerged$tw_group)
  for(tw_groupin in tw_group_list) {
    delta_list_tw <- delta_list[which(wdcMerged$tw_group==tw_groupin)]
    wdcMerged_tw <- wdcMerged[which(wdcMerged$tw_group==tw_groupin),]
    
    points_density <- get_points_density(points, theta1, wdcMerged_tw$tw[1])
    
    for(i in 1:nrow(points)) {
      lat2 <- points[i,"lat"]
      lon2 <- points[i,"lon"]
      #dis_v <- latlondistance(lat1, lon1, lat2, lon2)
      st_point_list <- splitchar(points$local_stations[i])
      list_obs <- st_point_list
      if(length(list_obs)==0) next
      
      station_data_i <- wdcMerged_tw[st_point_list,]
      dis_v <- latlondistance(points[i,"lat"], points[i,"lon"], 
                              station_data_i$lat, station_data_i$lon)
      sto_state = rep(0,nrow(station_data_i))
      util <- (beta1*dis_v) + (beta2*dis_v*dis_v) + delta_list_tw[st_point_list]
      exputil <- exp(util)
      util_st <- (!sto_state) * exputil
      den_util <- sum(util_st)
      lambda_st_t <- rep(0,length(util_st))
      
      for(k in 1:length(v0_vec)) {
        out <- exp(-v0_vec[k]*sigma0)
        prob_t <- util_st/(out+den_util)*(points_density[i]/length(v0_vec))
        dem_seg[ceiling(dis_v/0.01)] <- dem_seg[ceiling(dis_v/0.01)] + prob_t
        lambda_st_t <- lambda_st_t + prob_t
        outside_dem <- outside_dem + (points_density[i]/length(v0_vec)) - sum(prob_t)
        dem_closest_st <- dem_closest_st + prob_t[which.min(dis_v)]
        
        density_mass <- density_mass + points_density[i]/length(v0_vec)
        st_serv_lvl <- user_serv_lvl$serv_lvl[which(user_serv_lvl$tw_group==tw_groupin)][st_point_list]
        dem_highserv_st <- dem_highserv_st + prob_t[which.max(st_serv_lvl)]
      }
      lambda_t[list_obs] <- lambda_t[list_obs] + lambda_st_t
    }
  }
  return(list(lambda_t,dem_seg, outside_dem,dem_closest_st,dem_highserv_st,density_mass))
}


#########################
#Isodemand functions

optimal_distance_scaling <- function(theta1, delta_list, wdcMerged,points, 
                                     displace_serv_lvl, serv_lvl_coef, demand_actual) {
  
  #construct new delta for displace_serv_lvl (assumed wdcMerged$serv_lvl is proper,
  #can also altenatively construct from user_serv_lvl)
  serv_lvl_covar_new <- displace_serv_lvl+wdcMerged$serv_lvl
  serv_lvl_covar_new[which(serv_lvl_covar_new>1)] <- 1
  serv_lvl_covar_new[which(serv_lvl_covar_new<0)] <- 0
  
  stocked_list <- which(wdcMerged$stocked_out==F)
  
  delta_list_new <- delta_list + serv_lvl_coef[1]*((serv_lvl_covar_new-wdcMerged$serv_lvl)[stocked_list]) +
    serv_lvl_coef[2]*((serv_lvl_covar_new^2-wdcMerged$serv_lvl^2)[stocked_list])
  
  #construct new weights
  wdcMerged$obs_weight <- generate_weights_additive_servlvl(wdcMerged, user_serv_lvl, displace_serv_lvl)
  wdcMerged$obs_weight[which(wdcMerged$obs_weight==0)] <- 
    min(wdcMerged$obs_weight[-which(wdcMerged$obs_weight==0)])*1e-5
  
  ftol <- 1
  tol <- 0.01
  
  a <- 0.05
  b <- 16
  N <- 1
  NMAX <- 10
  tiny <- 1e-4
  
  func_a <- demand_distance_change_inner(1/a, theta1, delta_list_new, wdcMerged, points)
  func_b <- demand_distance_change_inner(1/b, theta1, delta_list_new, wdcMerged, points)
  
  while(N < NMAX) {
    c <- (a+b)/2  
    print(c)
    func_c <- demand_distance_change_inner(1/c, theta1, delta_list_new, wdcMerged, points)
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





