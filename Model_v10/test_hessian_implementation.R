source("constants_mw6.R")
sourceCpp("eval_func_2_new_deltaaveraged.cpp", verbose = T)

source("GetDataGooglePlaces_tiny.R")


weighing_GMM_mat <<- NULL

theta <- c(-4 ,1, 146.83, 0.052, 203.69, 0.42, 9.41)

set.seed(34675)
deltain_stin <- rnorm(length(which(!wdcMerged$stocked_out)),-3, 2)
eta <- rep(0,12)
length_eta <- length(eta)
length_theta <- length(theta)
length_delta <- length(deltain_stin)
params <- c(theta,deltain_stin, eta)

deltain <- rep(-30, nrow(wdcMerged))
deltain[which(wdcMerged$stocked_out==F)] <- deltain_stin
theta1 <- c(theta[1],0,theta[-1])
stocked_list <- which(wdcMerged$stocked_out==F)

lambda_multiplers_in <- rep(0,length(stocked_list))
lambda_multiplers_in[1] <- 1

# tw_group_list <- unique(wdcMerged$tw_group)
# i<- 1
# tw_groupin = tw_group_list[i]
# wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
# deltain_tw = deltain[which(wdcMerged$tw_group==tw_groupin)]
# lambda_multiplers_full <- rep(0, nrow(wdcMerged))
# lambda_multiplers_full[which(wdcMerged$stocked_out==F)] <- lambda_multiplers_in
# lambda_multiplers_tw = lambda_multiplers_full[which(wdcMerged$tw_group==tw_groupin)]


################################################################
################################################################

#check if the linearization of hessian is happening properly.
#Do this by first creating linear hessian in a full-proof way, ie. by creating
#entire matrix and then linearizing.

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

full_hess <- test_eval_hessian_lambda_constraints_MPEC(deltain_stin, theta1, 
                                                       eta, wdcMerged, points, lambda_multiplers_in)

a1_hess <- eval_hessian_lambda_constraints_MPEC(deltain_stin, theta1, 
                                     eta, wdcMerged, points, lambda_multiplers_in)


