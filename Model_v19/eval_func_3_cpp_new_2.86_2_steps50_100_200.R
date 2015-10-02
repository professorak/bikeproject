# sourceCpp("eval_func_3_new.cpp")
# 
# eval_lambda_theta <- function(theta1, deltain, wdcMergedday, points, tw_groupin) {
#   #for argument day compute the lambda values for all stations and times 
#   #in ascensing station_id n_time order
#   #this implies integrating over all points
#   #return will be a single vector of length no_stations X no_time in that day
#   return(eval_lambda_delta_list(deltain, theta1, wdcMergedday, points, tw_groupin)[[1]])    
# }  
# 
# eval_grad_lambda_theta_full <- function(theta1, deltain, wdcMergedday, points, tw_groupin) {
#   #for each point we have to generate the gradient wrt theta1 
#   #taking into account delta is function of theta
#   #for this will have to compute gradient of shares wrt delta and theta1
#   
#   #computing the gradients of station fixed effects wrt beta1
#   list_grad_share_delta <- eval_grad_share_delta (deltain, theta1, wdcMergedday, points, tw_groupin)
#   list_grad_share_theta <- eval_grad_share_theta (deltain, theta1, wdcMergedday, points, tw_groupin)
#   
#   grad_share_delta <- list_grad_share_delta[[1]]
#   grad_lambda_delta <- list_grad_share_delta[[2]]
#   grad_share_theta <- list_grad_share_theta[[1]]
#   grad_lambda_theta <- list_grad_share_theta[[2]]
#   
#   grad_delta_theta <- -solve(grad_share_delta) %*% as.matrix(grad_share_theta)
#   #the full gradient of lambda is given by gradient of lambda wrt theta 
#   #and that wrt delta is scaled by gradient of delta wrt theta 
#   grad_lambda_theta_full <- grad_lambda_theta + grad_lambda_delta %*% grad_delta_theta
#   
#   return (grad_lambda_theta_full)
# }  
# 
# 
eval_grad_share_delta_new <- function(deltain, theta1, wdcMergedday, points, tw_groupin) {
  #generates the gradient of shares sfT for dayin wrt \delta
  #calls eval_grad_lambda_delta and then aggregates.
  #output is of dimension no_stXno_st  
  grad_t <- eval_lambda_delta_list_new(deltain, theta1, wdcMergedday, points, tw_groupin)[[2]]
  
  grad_dem_T <- grad_t
  
  #returning both gradients of share and graident of lambda
  return(list(grad_dem_T,grad_t))
}

eval_grad_share_theta_new <- function(deltain, theta1, wdcMergedday, points, tw_groupin) {
  #generates the gradient of shares sfT for dayin wrt \delta
  #calls eval_grad_lambda_delta and then aggregates.
  #output is of dimension no_stXno_params 
  #note that this is partial gradient of share wrt theta and doesnt account delta is function of theta
  grad_t <- eval_grad_lambda_theta_new(theta1, deltain, wdcMergedday, points, tw_groupin)
  
  grad_dem_T <- grad_t
  return(list(grad_dem_T,grad_t))
}

eval_grad_lambda_theta_new <- function(theta1, deltain_tw, wdcMergedday, points, 
                                             tw_groupin) {
  #for argument day compute the lambda values for all stations and times 
  #in ascensing station_id n_time order
  #this implies integrating over all points
  #return will be a single vector of length no_stations X no_time in that day
  no_st <- max(wdcMergedday$station_id_index)
  tw_in <- wdcMergedday$tw[1]
  if(length(deltain_tw)!=nrow(wdcMergedday)) stop("error in eval_grad_lambda_theta")
  sto_state_local <- wdcMergedday$sto_state_local
  local_stations <- wdcMergedday$local_stations
  points_local_stations <- points$local_stations
  wdcMergedday  = wdcMergedday[,c("station_id",
                                  "stocked_out","station_id_index","lat","lon","obs_weight","out_dem_sum")]
  
  wdcMergedday = as.matrix(wdcMergedday)
  points_mat <- points
  points_mat$density <- get_points_density(points_mat, theta1, tw_in)
  points_mat = as.matrix(points_mat[,c("lat","lon","density")])

  grad_t <- eval_grad_lambda_theta_cpp_new(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
                                           as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations), nonden_ceoflength)

  #compute gradient wrt density_ridership and density_metro_col
  density_mat <- cbind(get_points_density_grad_ridership_col(points, tw_in)
                       ,
                       get_points_density_grad_metro_col(points, tw_in),  
                       get_points_density_grad_intercept_col(points, tw_in),  
                       get_points_density_grad_metro_evening_col(points, tw_in),
                       get_points_density_grad_places_count_col(points, tw_in),
                       get_points_density_grad_addden_cols(points, tw_in)
)
  points_mat <- cbind(points[,c("lat","lon")], density_mat)
  points_mat = as.matrix(points_mat)
  grad_t_densitycols <- eval_lambda_multiple_cpp_new(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
      as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations), nonden_ceoflength)  
  grad_t_all <- cbind(grad_t,grad_t_densitycols)  
  return(grad_t_all)  
}  

eval_grad_lambda_new <- function(theta1, deltain_tw, wdcMergedday, points, 
                                 tw_groupin) {
  #for argument day compute the lambda values for all stations and times 
  #in ascensing station_id n_time order
  #this implies integrating over all points
  #return will be a single vector of length no_stations X no_time in that day
  no_st <- max(wdcMergedday$station_id_index)
  tw_in <- wdcMergedday$tw[1]
  if(length(deltain_tw)!=nrow(wdcMergedday)) stop("error in eval_grad_lambda_theta")
  sto_state_local <- wdcMergedday$sto_state_local
  local_stations <- wdcMergedday$local_stations
  points_local_stations <- points$local_stations
  wdcMergedday  = wdcMergedday[,c("station_id",
                                  "stocked_out","station_id_index","lat","lon","obs_weight","out_dem_sum")]
  
  wdcMergedday = as.matrix(wdcMergedday)

  density_mat <- cbind(get_points_density_grad_ridership_col(points, tw_in)
                       ,
                       get_points_density_grad_metro_col(points, tw_in),  
                       get_points_density_grad_intercept_col(points, tw_in),  
                       get_points_density_grad_metro_evening_col(points, tw_in),
                       get_points_density_grad_places_count_col(points, tw_in),
                       get_points_density_grad_addden_cols(points, tw_in)
  )
  points_mat <- points
  points_mat$density <- get_points_density(points_mat, theta1, tw_in)
  points_mat <- cbind(points_mat[,c("lat","lon","density")], density_mat)
  points_mat = as.matrix(points_mat)
  grad_t_all <- eval_grad_lambda_cpp_new(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
      as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations),nonden_ceoflength)  
  return(grad_t_all)  
}  



# eval_grad_lambda_theta_new <- function(theta1, deltain_tw, wdcMergedday, points, tw_groupin) {
#   #for argument day compute the lambda values for all stations and times 
#   #in ascensing station_id n_time order
#   #this implies integrating over all points
#   #return will be a single vector of length no_stations X no_time in that day
#   no_st <- max(wdcMergedday$station_id_index)
#   tw_in <- wdcMergedday$tw[1]
#   if(length(deltain_tw)!=nrow(wdcMergedday)) stop("error in eval_grad_lambda_theta")
#   sto_state_local <- wdcMergedday$sto_state_local
#   local_stations <- wdcMergedday$local_stations
#   points_local_stations <- points$local_stations
#   wdcMergedday  = wdcMergedday[,c("station_id",
#                                   "stocked_out","station_id_index","lat","lon","obs_weight","out_dem_sum")]
#   
#   wdcMergedday = as.matrix(wdcMergedday)
#   points_mat <- points
#   points_mat$density <- get_points_density(points_mat, theta1, tw_in)
#   points_mat = as.matrix(points_mat[,c("lat","lon","density")])
#   
#   
#   grad_t <- eval_grad_lambda_theta_cpp_new(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
#                                     as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations))
#   
#   #compute gradient wrt density_ridership and density_metro_col
#   points_mat <- points
#   #points_mat$density <- ((points_mat$type==1)) *points_mat$weight
#   points_mat$density <- get_points_density_grad_ridership_col(points_mat, tw_in)
#   points_mat = as.matrix(points_mat[,c("lat","lon","density")])  
#   res <- eval_lambda_delta_list_cpp_new(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
#                                         as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations))  
#   grad_t_ridership <- res[,1]
#   points_mat <- points  
#   #points_mat$density <- ((points_mat$type==2)) *points_mat$weight
#   points_mat$density <- get_points_density_grad_metro_col(points_mat, tw_in)
#   points_mat = as.matrix(points_mat[,c("lat","lon","density")])  
#   res <- eval_lambda_delta_list_cpp_new(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
#                                         as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations))  
#   grad_t_metro <- res[,1]
#   points_mat <- points
#   #points_mat$density <- ((points_mat$type==1))
#   points_mat$density <- get_points_density_grad_intercept_col(points_mat, tw_in)
#   points_mat = as.matrix(points_mat[,c("lat","lon","density")])  
#   res <- eval_lambda_delta_list_cpp_new(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
#                                         as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations))  
#   grad_t_intercept <- res[,1]
#   
#   points_mat <- points  
#   #points_mat$density <- ((points_mat$type==2)) *points_mat$weight
#   points_mat$density <- get_points_density_grad_metro_evening_col(points_mat, tw_in)
#   points_mat = as.matrix(points_mat[,c("lat","lon","density")])  
#   res <- eval_lambda_delta_list_cpp_new(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
#                                         as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations))  
#   grad_t_metro_evening <- res[,1]
#   
#   points_mat <- points  
#   #points_mat$density <- ((points_mat$type==2)) *points_mat$weight
#   points_mat$density <- get_points_density_grad_places_count_col(points_mat, tw_in)
#   points_mat = as.matrix(points_mat[,c("lat","lon","density")])  
#   res <- eval_lambda_delta_list_cpp_new(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
#                                         as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations))  
#   grad_t_places_count <- res[,1]
#   
#   points_mat <- points  
#   #points_mat$density <- ((points_mat$type==2)) *points_mat$weight
#   points_mat$density <- get_points_density_grad_bus_col(points_mat, tw_in)
#   points_mat = as.matrix(points_mat[,c("lat","lon","density")])  
#   res <- eval_lambda_delta_list_cpp_new(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
#                                         as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations))  
#   grad_t_bus <- res[,1]
#   
#   grad_t_all <- cbind(grad_t,grad_t_ridership,grad_t_metro,grad_t_intercept,grad_t_metro_evening,
#                       grad_t_places_count,grad_t_bus)
#   
#   return(grad_t_all)  
# }

eval_grad_delta_theta_new <- function(theta1, deltain, wdcMergedday, points, tw_groupin) {
  #for each point we have to generate the gradient wrt theta1 
  #taking into account delta is function of theta
  #for this will have to compute gradient of shares wrt delta and theta1
  
  #computing the gradients of station fixed effects wrt beta1
  grad_share_delta <- eval_lambda_delta_list_new(deltain, theta1, wdcMergedday, points, tw_groupin)[[2]]
  grad_share_theta <- eval_grad_lambda_theta_new(theta1, deltain, wdcMergedday, points, tw_groupin)
  
  stocked_list <- which(wdcMergedday$stocked_out==FALSE)
  print("kappa grad_share_delta: ")
  print(kappa(grad_share_delta[stocked_list,stocked_list]))
  grad_delta_theta <- -solve(grad_share_delta[stocked_list,stocked_list]) %*% 
                            as.matrix(grad_share_theta[stocked_list,])
  grad_delta_theta_full <- matrix(0,nrow(wdcMergedday),ncol(grad_delta_theta))
  grad_delta_theta_full[stocked_list,] <- grad_delta_theta
  return(grad_delta_theta_full)
}


eval_grad_delta_theta_full <- function(theta1, deltain, wdcMerged, points) {
  #assumes wdcMerged is first sorted by tw_group
  grad_delta_theta <- c()
  tw_groupin_list <- unique(wdcMerged$tw_group)   
  for(tw_groupin in tw_groupin_list) {
    wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
    tw_index <- which(tw_groupin_list==tw_groupin)
    deltain_tw <- deltain[which(wdcMerged$tw_group==tw_groupin)]
    grad_delta_theta <- rbind(grad_delta_theta,
                              eval_grad_delta_theta_new(theta1, deltain_tw, wdcMergedday, points, tw_groupin))
  }  
  return(grad_delta_theta)
}

eval_hessian_delta_sq_tw <- function(deltain_tw, theta1, wdcMergedday, points, tw_groupin,lambda_multiplers_in) {
  no_st <- length(unique(wdcMergedday$station_id_index))
  no_obs <- nrow(wdcMergedday)
  tw_in <- wdcMergedday$tw[1]
  if(length(deltain_tw)!=no_obs) stop("error in eval_lambda_delta_list")  
  sto_state_local <- wdcMergedday$sto_state_local
  local_stations <- wdcMergedday$local_stations
  points_local_stations <- points$local_stations
  wdcMergedday  = wdcMergedday[,c("station_id",
                                  "stocked_out","station_id_index","lat","lon","obs_weight","out_dem_sum")]
  
  
  points_mat <- points
  points_mat$density <- get_points_density(points_mat, theta1, tw_in)
  points_mat = as.matrix(points_mat[,c("lat","lon","density")])
  wdcMergedday = as.matrix(wdcMergedday)
  
  hessian_lambda_delta_sq <- eval_hessian_lambda_delta_sq_cpp(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
                                                              as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations),
                                                              lambda_multiplers_in, nonden_ceoflength)
  
  return(hessian_lambda_delta_sq)  
  
}

eval_hessian_beta1_sq_tw <- function(deltain_tw, theta1, wdcMergedday, points, tw_groupin,lambda_multiplers_in) {
  no_st <- length(unique(wdcMergedday$station_id_index))
  no_obs <- nrow(wdcMergedday)
  tw_in <- wdcMergedday$tw[1]
  if(length(deltain_tw)!=no_obs) stop("error in eval_lambda_delta_list")  
  sto_state_local <- wdcMergedday$sto_state_local
  local_stations <- wdcMergedday$local_stations
  points_local_stations <- points$local_stations
  wdcMergedday  = wdcMergedday[,c("station_id",
                                  "stocked_out","station_id_index","lat","lon","obs_weight","out_dem_sum")]
  
  
  points_mat <- points
  points_mat$density <- get_points_density(points_mat, theta1, tw_in)
  points_mat = as.matrix(points_mat[,c("lat","lon","density")])
  wdcMergedday = as.matrix(wdcMergedday)
  
  hessian_lambda_beta1_sq <- eval_hessian_lambda_beta1_sq_cpp(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
                                                              as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations),
                                                              lambda_multiplers_in, nonden_ceoflength)
  
  return(hessian_lambda_beta1_sq)  
  
}

eval_hessian_theta1_sq_tw <- function(deltain_tw, theta1, wdcMergedday, points, tw_groupin,lambda_multiplers_in) {
  no_st <- length(unique(wdcMergedday$station_id_index))
  no_obs <- nrow(wdcMergedday)
  tw_in <- wdcMergedday$tw[1]
  if(length(deltain_tw)!=no_obs) stop("error in eval_lambda_delta_list")  
  sto_state_local <- wdcMergedday$sto_state_local
  local_stations <- wdcMergedday$local_stations
  points_local_stations <- points$local_stations
  wdcMergedday  = wdcMergedday[,c("station_id",
                                  "stocked_out","station_id_index","lat","lon","obs_weight","out_dem_sum")]
  
  density_mat <- cbind(get_points_density_grad_ridership_col(points, tw_in)
                       ,
                       get_points_density_grad_metro_col(points, tw_in),  
                       get_points_density_grad_intercept_col(points, tw_in),  
                       get_points_density_grad_metro_evening_col(points, tw_in),
                       get_points_density_grad_places_count_col(points, tw_in),
                       get_points_density_grad_addden_cols(points, tw_in)
  )
  points_mat <- points
  points_mat$density <- get_points_density(points_mat, theta1, tw_in)  
  points_mat <- cbind(points_mat[,c("lat","lon","density")], density_mat)
  points_mat = as.matrix(points_mat)
  
  wdcMergedday = as.matrix(wdcMergedday)
  
  hessian_lambda_theta1_sq <- eval_hessian_lambda_theta1_sq_cpp(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
                                                              as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations),
                                                              lambda_multiplers_in, nonden_ceoflength)
  print(hessian_lambda_theta1_sq[1,1])
  return(hessian_lambda_theta1_sq)  
  
}

eval_hessian_lambda_theta1_delta_tw <- function(deltain_tw, theta1, wdcMergedday, points, tw_groupin,lambda_multiplers_in) {
  no_st <- length(unique(wdcMergedday$station_id_index))
  no_obs <- nrow(wdcMergedday)
  tw_in <- wdcMergedday$tw[1]
  if(length(deltain_tw)!=no_obs) stop("error in eval_lambda_delta_list")  
  sto_state_local <- wdcMergedday$sto_state_local
  local_stations <- wdcMergedday$local_stations
  points_local_stations <- points$local_stations
  wdcMergedday  = wdcMergedday[,c("station_id",
                                  "stocked_out","station_id_index","lat","lon","obs_weight","out_dem_sum")]
  
  density_mat <- cbind(get_points_density_grad_ridership_col(points, tw_in)
                       ,
                       get_points_density_grad_metro_col(points, tw_in),  
                       get_points_density_grad_intercept_col(points, tw_in),  
                       get_points_density_grad_metro_evening_col(points, tw_in),
                       get_points_density_grad_places_count_col(points, tw_in),
                       get_points_density_grad_addden_cols(points, tw_in)
  )
  points_mat <- points
  points_mat$density <- get_points_density(points_mat, theta1, tw_in)  
  points_mat <- cbind(points_mat[,c("lat","lon","density")], density_mat)
  points_mat = as.matrix(points_mat)
  
  wdcMergedday = as.matrix(wdcMergedday)
  
  hessian_lambda_theta1_delta <- eval_hessian_lambda_theta1_delta_cpp(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
                                                                as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations),
                                                                lambda_multiplers_in, nonden_ceoflength)
  
  return(hessian_lambda_theta1_delta)  
  
}

eval_hessian_lambda_constraints_MPEC_tw_groupin <- function(deltain_tw, theta1, wdcMergedday, points, tw_groupin, lambda_multiplers_tw) {
  no_st <- length(unique(wdcMergedday$station_id_index))
  no_obs <- nrow(wdcMergedday)
  tw_in <- wdcMergedday$tw[1]
  stocked_list <- which(wdcMergedday$stocked_out==F)
  if(length(deltain_tw)!=no_obs) stop("error in eval_lambda_delta_list")  
  sto_state_local <- wdcMergedday$sto_state_local
  local_stations <- wdcMergedday$local_stations
  points_local_stations <- points$local_stations
  wdcMergedday  = wdcMergedday[,c("station_id",
                                  "stocked_out","station_id_index","lat","lon","obs_weight","out_dem_sum")]
  
  density_mat <- cbind(get_points_density_grad_ridership_col(points, tw_in)
                       ,
                       get_points_density_grad_metro_col(points, tw_in),  
                       get_points_density_grad_intercept_col(points, tw_in),  
                       get_points_density_grad_metro_evening_col(points, tw_in),
                       get_points_density_grad_places_count_col(points, tw_in),
                       get_points_density_grad_addden_cols(points, tw_in)
  )
  points_mat <- points
  points_mat$density <- get_points_density(points_mat, theta1, tw_in)  
  points_mat <- cbind(points_mat[,c("lat","lon","density")], density_mat)
  points_mat = as.matrix(points_mat)
  
  wdcMergedday = as.matrix(wdcMergedday)
  
  hessian_lambda_all <- eval_hessian_lambda_cpp(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
                                                as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations),
                                                lambda_multiplers_tw, nonden_ceoflength)
  theta_index <- c(1:length(theta1))
  delta_index <- length(theta1)+c(1:nrow(wdcMergedday))
  hessian_lambda_theta1_sq <- hessian_lambda_all[theta_index,theta_index]
  hessian_lambda_theta1_delta <- hessian_lambda_all[theta_index,delta_index]
  hessian_lambda_delta_sq <- hessian_lambda_all[delta_index,delta_index]
  
  #   hessian_lambda_constraints_tw_groupin <- matrix(0,length(theta1)+nrow(wdcMergedday),length(theta1)+nrow(wdcMergedday))
  #   theta_range <- c(1:length(theta1))
  #   delta_range <- c((1+length(theta1)):(length(theta1)+nrow(wdcMergedday)))
  #   #theta1_sq
  #   hessian_lambda_constraints_tw_groupin[theta_range,theta_range] <- hessian_lambda_theta1_sq
  #   
  #   #delta_sq
  #   hessian_lambda_constraints_tw_groupin[delta_range,delta_range] <- hessian_lambda_delta_sq
  #   
  #   #theta1_delta_sq
  #   hessian_lambda_constraints_tw_groupin[theta_range,delta_range] <- hessian_lambda_theta1_delta
  #   hessian_lambda_constraints_tw_groupin[delta_range,theta_range] <- t(hessian_lambda_theta1_delta)
  
  #   #remove 2 index of theta1
  #   hessian_lambda_constraints_tw_groupin <- hessian_lambda_constraints_tw_groupin[-2,-2]
  #   #keep only stocked_list of delta indexes
  #   stocked_list <- which(wdcMergedday$stocked_out==F)
  #   keep_index <- c(c(1:length(theta)),length(theta)+stocked_list)
  #   hessian_lambda_constraints_tw_groupin <- hessian_lambda_constraints_tw_groupin[keep_index,keep_index]
  #I cannot return hessian_lambda_constraints_tw_groupin seperately for each tw as a vector, since the 
  #input from other time windows will add to some of the same rows (corresponding to theta1)
  #so will return each component separately
  
  #remove 2 index of theta1
  hessian_lambda_theta1_sq <- hessian_lambda_theta1_sq[-2,-2]
  hessian_lambda_theta1_sq <- keep_lower_traingular_matrix(hessian_lambda_theta1_sq)
  hessian_lambda_theta1_delta <- hessian_lambda_theta1_delta[-2,]
  #keep only stocked_list of delta indexes    
  hessian_lambda_theta1_delta <- hessian_lambda_theta1_delta[,stocked_list]
  hessian_lambda_delta_sq <- hessian_lambda_delta_sq[stocked_list,stocked_list]
  hessian_lambda_delta_sq <- keep_lower_traingular_matrix(hessian_lambda_delta_sq)
  hessian_lambda_delta_theta1 <- t(hessian_lambda_theta1_delta)
  hessian_lambda_theta1_delta[,] <- 0 #since it is all in upper triangle portion
  
  #linearize and sparsify
  #not compressing hessian_lambda_theta1_sq and hessian_lambda_theta1_delta
  hessian_lambda_delta_vs_theta1_n_delta <- cbind(hessian_lambda_delta_theta1,hessian_lambda_delta_sq)
  hessian_lambda_delta_vs_theta1_n_delta <- linearize_sparsify(hessian_lambda_delta_vs_theta1_n_delta)
  
  return(list("hessian_lambda_theta1_sq"=hessian_lambda_theta1_sq,
              "hessian_lambda_theta1_delta"=hessian_lambda_theta1_delta,
              "hessian_lambda_delta_vs_theta1_n_delta_lin"=hessian_lambda_delta_vs_theta1_n_delta
  ))
}


eval_hessian_lambda_constraints_MPEC <- function(deltain_stin, theta1, eta, wdcMerged, points, lambda_multiplers_in) {  
  #expand deltain to all observations, it is currently #stocked in observations
  deltain <- rep(-30, nrow(wdcMerged))
  deltain[which(wdcMerged$stocked_out==F)] <- deltain_stin
  lambda_multiplers_full <- rep(0, nrow(wdcMerged))
  lambda_multiplers_full[which(wdcMerged$stocked_out==F)] <- lambda_multiplers_in  
  # theta <- theta1[-2]
  # stocked_list <- which(wdcMerged$stocked_out==F)
  
  tw_group_list <- unique(wdcMerged$tw_group)
  hessian_lambda_theta1_sq <- matrix(0,length(theta),length(theta))
  hessian_lambda_theta1_delta <- c()
  hessian_lambda_delta_vs_theta1_n_delta <- c()
  
  for(i in 1:length(tw_group_list)) {
    tw_groupin = tw_group_list[i]
    wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
    deltain_tw = deltain[which(wdcMerged$tw_group==tw_groupin)]
    lambda_multiplers_tw = lambda_multiplers_full[which(wdcMerged$tw_group==tw_groupin)]
    ret <- eval_hessian_lambda_constraints_MPEC_tw_groupin(deltain_tw, theta1, 
                                                           wdcMergedday, points, tw_groupin, lambda_multiplers_tw)
    hessian_lambda_theta1_sq <- hessian_lambda_theta1_sq + ret[[1]]
    hessian_lambda_theta1_delta <- cbind(hessian_lambda_theta1_delta,ret[[2]])
    hessian_lambda_delta_vs_theta1_n_delta <- c(hessian_lambda_delta_vs_theta1_n_delta,
                                                ret[[3]])
  }
  hessian_lambda_constraints_MPEC <- c(
    linearize_sparsify(cbind(hessian_lambda_theta1_sq,hessian_lambda_theta1_delta)),
    hessian_lambda_delta_vs_theta1_n_delta)
  
  return(hessian_lambda_constraints_MPEC)  
}

  
linearize_sparsify <- function(A) {
  B <- c(t(A))
  return(B[which(B!=0)])
}
