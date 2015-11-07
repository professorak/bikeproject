#implementing a function which returns lambda_sum for each element of v0_vec seperately
#which is multiplied by corresponding v0_vec_weights
sourceCpp("eval_func_2_new_deltaaveraged_stepfunction_v0weights_set2.cpp") 

eval_lambdasum_v0_vec <- function(deltain_stin, theta1, wdcMerged, points) {
  
  #expand deltain to all observations, it is currently #stocked in observations
  deltain <- rep(-30, nrow(wdcMerged))
  deltain[which(wdcMerged$stocked_out==F)] <- deltain_stin
  
  tw_group_list <- unique(wdcMerged$tw_group)
  obj <- rep(0,length(v0_vec))
  for(i in 1:length(tw_group_list)) {
    tw_groupin = tw_group_list[i]
    wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
    deltain_tw = deltain[which(wdcMerged$tw_group==tw_groupin)]
    obj <- obj + 
             eval_lambdasum_v0_vec_tw_groupin(deltain_tw, theta1, wdcMergedday, 
                                    points, tw_groupin)
  }
  return(obj)
}



eval_lambdasum_v0_vec_tw_groupin <- function(deltain_tw, theta1, wdcMergedday, points, tw_groupin) {
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
  
  lambda_t <- eval_lambdasum_v0_vec_cpp_new(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec, v0_vec_weights,
                                  as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations), nonden_ceoflength)
  
  return(lambda_t)    
}



random_coef_acceptance_prob <- function(v_i,
                                        deltain_stin, theta1, wdcMerged, points) {
  #returns the probability that a random coefficient generated should be accepeted. 
  #its f_bar/s_bar
  #check ImportanceSampling.Lyx in Notes folder for more details
  
  #compute f_j, which is just predicted lambda for this v_i.
  v0_vec_org <- v0_vec
  v0_vec_weights_org <- v0_vec_weights
  v0_vec <<- v_i
  v0_vec_weights <<- rep(1,length(v0_vec))
  f_j <- eval_lambdasum_v0_vec(deltain_stin, theta1, wdcMerged, points)
  #share sum
  #fs_bar
  fs_bar <- f_j/get_total_density(theta1, wdcMerged, points)
  v0_vec <<- v0_vec_org
  v0_vec_weights <<- v0_vec_weights_org
  return(fs_bar)
}

generate_v0_importance_sampling <- function(no_vals = 40,
                                            deltain_stin, theta1, wdcMerged, points) {
  s_bar <- sum(wdcMerged$out_dem_mean)/get_total_density(theta1, wdcMerged, points)
  
  v0_list <- c()
  v0_weight_list <- c()
  while(length(v0_list)<no_vals) {
    v_i <- rnorm((no_vals-length(v0_list))*20)
    prob_i <- random_coef_acceptance_prob(v_i, deltain_stin, theta1, wdcMerged, points)
    unif_i <- runif(length(prob_i))
    if(length(which(unif_i<=prob_i))) {
      v0_list <- c(v0_list, v_i[which(unif_i<=prob_i)])
      v0_weight_list <- c(v0_weight_list, (s_bar/prob_i)[which(unif_i<=prob_i)])
    }
  }
  return(cbind(v0_list[c(1:no_vals)],v0_weight_list[c(1:no_vals)]))
}
