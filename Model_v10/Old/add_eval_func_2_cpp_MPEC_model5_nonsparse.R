
eval_hessian <- function(params, wdcMerged, points, length_theta, length_delta, length_eta, obj_factor, constraint_multipliers) {
  print("In eval_hessian: ")
  ptm <- proc.time()  
  theta <- params[c(1:length_theta)]
  deltain_stin <- params[c((length_theta+1):(length_theta+length_delta))]
  eta <- params[c((length_theta+length_delta+1):(length_theta+length_delta+length_eta))]
  theta1 <- c(theta[1],0,theta[-1])  
  #making multipliers non-zero so that the return values confirm to the sparse strcuture.
  obj_factor <- obj_factor+1e-16
  constraint_multipliers <- constraint_multipliers+1e-16
  lambda_multiplers_in <- constraint_multipliers[c(1:length_delta)]
  
  ret <- c(eval_hessian_lambda_constraints_MPEC(deltain_stin, theta1, 
                                                eta, wdcMerged, points, lambda_multiplers_in),
           eval_hessian_eta_sq(obj_factor, length_eta)
           )
  #make 0 the values less than 1e-16, they are only non-zero due to above shifting of multiplers and 
  #will otherwise cause issues with condition number of matrix.
  ret[which(abs(ret)<=1e-16)] <- 0
  print("eval time:")
  print(proc.time()-ptm)
  return(ret)
}

eval_hessian_eta_sq <- function(obj_factor, length_eta) {  
  if(is.null(weighing_GMM_mat)) {
    A_N <- diag(length_eta)    
  } else {
    A_N <- weighing_GMM_mat
  }
  return(c(t(2*A_N)))  
}





eval_hessian_structure <- function(params, wdcMerged, points, length_theta, length_delta, length_eta) {
  theta <- params[c(1:length_theta)]
  deltain_stin <- params[c((length_theta+1):(length_theta+length_delta))]
  eta <- params[c((length_theta+length_delta+1):(length_theta+length_delta+length_eta))]
  theta1 <- c(theta[1],0,theta[-1])  
  
  return(c(eval_hessian_structure_lambda_constraints_MPEC(deltain_stin, theta1, eta, wdcMerged, points),
           eval_hessian_structure_eta_sq(length_theta, length_delta,length_eta)
  ))
}

eval_hessian_structure_lambda_constraints_MPEC <- function(deltain_stin, theta1, eta, wdcMerged, points) {
  #expand deltain to all observations, it is currently #stocked in observations
  deltain <- rep(-30, nrow(wdcMerged))
  deltain[which(wdcMerged$stocked_out==F)] <- deltain_stin
  # theta <- theta1[-2]
  # stocked_list <- which(wdcMerged$stocked_out==F)
  
  tw_group_list <- unique(wdcMerged$tw_group)
  hessian_lambda_theta1_sq <- matrix(0,length(theta),length(theta))
  hessian_lambda_theta1_delta <- c()
  hessian_structure_lambda_delta_vs_theta1_n_delta <- c()
  span_start <- length(theta)
  
  for(i in 1:length(tw_group_list)) {
    tw_groupin = tw_group_list[i]
    wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
    deltain_tw = deltain[which(wdcMerged$tw_group==tw_groupin)]
    ret <- eval_hessian_structure_lambda_constraints_MPEC_tw_groupin(deltain_tw, theta1, 
              wdcMergedday, points, tw_groupin, span_start)
    hessian_lambda_theta1_sq <- hessian_lambda_theta1_sq + ret[[1]]
    hessian_lambda_theta1_delta <- cbind(hessian_lambda_theta1_delta,ret[[2]])
    hessian_structure_lambda_delta_vs_theta1_n_delta <- c(hessian_structure_lambda_delta_vs_theta1_n_delta,
                                                ret[[3]])
    span_start <- span_start + ret[[4]]
  }
  hessian_structure_lambda_constraints_MPEC <- c(
    make.sparse(cbind(hessian_lambda_theta1_sq,hessian_lambda_theta1_delta)),
    hessian_structure_lambda_delta_vs_theta1_n_delta)
  
  return(hessian_structure_lambda_constraints_MPEC)  
}

eval_hessian_structure_lambda_constraints_MPEC_tw_groupin <- function(deltain_tw, theta1, 
    wdcMergedday, points, tw_groupin, span_start) {
  no_st <- length(unique(wdcMergedday$station_id_index))
  no_obs <- nrow(wdcMergedday)
  tw_in <- wdcMergedday$tw[1]
  #lambda_multiplers_tw set to all 1's 
  lambda_multiplers_tw <- rep(1,nrow(wdcMergedday))
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
                       get_points_density_grad_bus_col(points, tw_in)
  )
  points_mat <- points
  points_mat$density <- get_points_density(points_mat, theta1, tw_in)  
  points_mat <- cbind(points_mat[,c("lat","lon","density")], density_mat)
  points_mat = as.matrix(points_mat)
  
  wdcMergedday = as.matrix(wdcMergedday)
  
  hessian_lambda_theta1_delta <- eval_hessian_lambda_theta1_delta_cpp(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
                                                                      as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations),
                                                                      lambda_multiplers_tw)
  hessian_lambda_delta_sq <- eval_hessian_lambda_delta_sq_cpp(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
                                                              as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations),
                                                              lambda_multiplers_tw)
  hessian_lambda_theta1_sq <- eval_hessian_lambda_theta1_sq_cpp(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
                                                                as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations),
                                                                lambda_multiplers_tw)
  
  #remove 2 index of theta1
  hessian_lambda_theta1_sq <- hessian_lambda_theta1_sq[-2,-2]
  hessian_lambda_theta1_delta <- hessian_lambda_theta1_delta[-2,]
  #keep only stocked_list of delta indexes    
  hessian_lambda_theta1_delta <- hessian_lambda_theta1_delta[,stocked_list]
  hessian_lambda_delta_sq <- hessian_lambda_delta_sq[stocked_list,stocked_list]
  hessian_lambda_delta_theta1 <- t(hessian_lambda_theta1_delta)
  
  #not compressing hessian_lambda_theta1_sq and hessian_lambda_theta1_delta
  hessian_structure_lambda_delta_vs_theta1_n_delta <- my_make_sparse(hessian_lambda_delta_theta1,
      hessian_lambda_delta_sq,span_start)
  
  return(list("hessian_lambda_theta1_sq"=hessian_lambda_theta1_sq,
              "hessian_lambda_theta1_delta"=hessian_lambda_theta1_delta,
              "hessian_lambda_delta_vs_theta1_n_delta_lin"=hessian_structure_lambda_delta_vs_theta1_n_delta,
              "span_shift"=ncol(hessian_lambda_delta_sq)
  ))
}


eval_hessian_structure_eta_sq <- function(length_theta, length_delta,length_eta) {
  hessian_eta_sq <- matrix(1,length_eta,length_eta)
  return(my_make_sparse_2(hessian_eta_sq,length_theta+length_delta))
}


