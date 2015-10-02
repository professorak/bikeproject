
test_eval_hessian_lambda_constraints_MPEC <- function(deltain_stin, theta1, eta, wdcMerged, points, lambda_multiplers_in) {  
  #expand deltain to all observations, it is currently #stocked in observations
  deltain <- rep(-30, nrow(wdcMerged))
  deltain[which(wdcMerged$stocked_out==F)] <- deltain_stin
  lambda_multiplers_full <- rep(0, nrow(wdcMerged))
  lambda_multiplers_full[which(wdcMerged$stocked_out==F)] <- lambda_multiplers_in
  theta <- theta1[-2]
  stocked_list <- which(wdcMerged$stocked_out==F)
  
  tw_group_list <- unique(wdcMerged$tw_group)
  hessian_lambda <- matrix(0, length(theta)+length(stocked_list),
                           length(theta)+length(stocked_list))  
  
  for(i in 1:length(tw_group_list)) {
    tw_groupin = tw_group_list[i]
    wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
    deltain_tw = deltain[which(wdcMerged$tw_group==tw_groupin)]
    lambda_multiplers_tw = lambda_multiplers_full[which(wdcMerged$tw_group==tw_groupin)]
    ret <- test_eval_hessian_lambda_constraints_MPEC_tw_groupin(deltain_tw, theta1, 
                                                                wdcMergedday, points, tw_groupin, lambda_multiplers_tw)
    
    delta_range <- which(wdcMerged$tw_group[stocked_list]==tw_groupin)
    row_range <- c(c(1:length(theta)), length(theta)+delta_range)
    hessian_lambda[row_range,row_range] <- hessian_lambda[row_range,row_range] + ret
    
  }
  hessian_lambda <- linearize_sparsify(hessian_lambda)
  
  return(hessian_lambda)  
}

test_eval_hessian_lambda_constraints_MPEC_tw_groupin <- function(deltain_tw, theta1, wdcMergedday, points, tw_groupin, lambda_multiplers_in) {
  no_st <- length(unique(wdcMergedday$station_id_index))
  no_obs <- nrow(wdcMergedday)
  stocked_list <- which(wdcMergedday$stocked_out==F)
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
                       get_points_density_grad_bus_col(points, tw_in)
  )
  points_mat <- points
  points_mat$density <- get_points_density(points_mat, theta1, tw_in)  
  points_mat <- cbind(points_mat[,c("lat","lon","density")], density_mat)
  points_mat = as.matrix(points_mat)
  
  wdcMergedday = as.matrix(wdcMergedday)
  
  hessian_lambda_theta1_delta <- eval_hessian_lambda_theta1_delta_cpp(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
                                                                      as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations),
                                                                      lambda_multiplers_in)
  hessian_lambda_delta_sq <- eval_hessian_lambda_delta_sq_cpp(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
                                                              as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations),
                                                              lambda_multiplers_in)
  hessian_lambda_theta1_sq <- eval_hessian_lambda_theta1_sq_cpp(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
                                                                as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations),
                                                                lambda_multiplers_in)
  hessian_lambda_constraints_tw_groupin <- matrix(0,length(theta1)+nrow(wdcMergedday),length(theta1)+nrow(wdcMergedday))
  theta_range <- c(1:length(theta1))
  delta_range <- c((1+length(theta1)):(length(theta1)+nrow(wdcMergedday)))
  #theta1_sq
  hessian_lambda_constraints_tw_groupin[theta_range,theta_range] <- hessian_lambda_theta1_sq
  
  #delta_sq
  hessian_lambda_constraints_tw_groupin[delta_range,delta_range] <- hessian_lambda_delta_sq
  
  #theta1_delta_sq
  hessian_lambda_constraints_tw_groupin[theta_range,delta_range] <- hessian_lambda_theta1_delta
  hessian_lambda_constraints_tw_groupin[delta_range,theta_range] <- t(hessian_lambda_theta1_delta)
  
  #remove 2 index of theta1
  hessian_lambda_constraints_tw_groupin <- hessian_lambda_constraints_tw_groupin[-2,-2]
  #keep only stocked_list of delta indexes
  
  keep_index <- c(c(1:length(theta)),length(theta)+stocked_list)
  hessian_lambda_constraints_tw_groupin <- hessian_lambda_constraints_tw_groupin[keep_index,keep_index]
  return(hessian_lambda_constraints_tw_groupin)
}

compute_lagrange  <- function(params, obj_factor, constraint_multipliers) {
  lagrange_function <- obj_factor * eval_obj_GMM_model5_obj(params,wdcMerged, points, 
                                                            length(theta), length(deltain_stin), length_eta) +
    (eval_g(params, wdcMerged, points, length_theta, length_delta, length_eta) %*% constraint_multipliers) 
  return(lagrange_function)
}

compute_hess_atindex <- function(params, index_1, index_2, obj_factor, constraint_multipliers) {
  
  a1 <- compute_lagrange(params, obj_factor, constraint_multipliers)
  diff1_vec <- rep(0, length(params))
  diff1_vec[index_1] <- diff
  diff2_vec <- rep(0, length(params))
  diff2_vec[index_2] <- diff
  #gradient wrt delta1_index at base value
  a1_2 <- compute_lagrange(params+diff1_vec, obj_factor, constraint_multipliers)
  a1_grad <- (a1_2-a1)/diff
  
  a2 <- compute_lagrange(params+diff2_vec, obj_factor, constraint_multipliers)
  #a2_hess <- eval_hessian_delta_list_new(deltain_tw_3, theta1, wdcMergedday, points, tw_groupin,lambda_multiplers_in)
  
  a2_2 <- compute_lagrange(params+diff1_vec+diff2_vec, obj_factor, constraint_multipliers)
  a2_grad <- (a2_2-a2)/diff
  
  ###
  a1_hess_num <- (a2_grad-a1_grad)/diff
  return(a1_hess_num)
  
}

