
eval_lambda_delta_list_new <- function(deltain_tw, theta1, wdcMergedday, points, tw_groupin) {
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
  
  res <- eval_lambda_delta_list_cpp_new(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec, v0_vec_weights,
      as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations), nonden_ceoflength)
  
  lambda_t <- res[,1]
  grad_t <- res[,c(2:ncol(res))]
  
  return(list("objective"=lambda_t,
              "gradient"=grad_t))  
  
}

eval_lambda_new <- function(deltain_tw, theta1, wdcMergedday, points, tw_groupin) {
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
  
  lambda_t <- eval_lambda_cpp_new(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec, v0_vec_weights,
      as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations), nonden_ceoflength)
  
  return(lambda_t)    
}

eval_lambda_full_weighted <- function(deltain_stin, theta1, wdcMerged, points) {
  return(eval_lambda_full_unweighted(deltain_stin, theta1, wdcMerged, points) * 
    wdcMerged$obs_weight)
}

eval_lambda_full_unweighted <- function(deltain_stin, theta1, wdcMerged, points) {
  
  #expand deltain to all observations, it is currently #stocked in observations
  deltain <- rep(-30, nrow(wdcMerged))
  deltain[which(wdcMerged$stocked_out==F)] <- deltain_stin
  
  tw_group_list <- unique(wdcMerged$tw_group)
  obj <- c()
  for(i in 1:length(tw_group_list)) {
    tw_groupin = tw_group_list[i]
    wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
    deltain_tw = deltain[which(wdcMerged$tw_group==tw_groupin)]    
    obj <- c(obj, 
             eval_lambda_tw_groupin(deltain_tw, theta1, wdcMergedday, 
                                    points, tw_groupin))
  }
  return(obj)
}

eval_lambda_tw_groupin <- function(deltain, theta1, wdcMergedday, points, tw_groupin) {
  stocked_list <- which(wdcMergedday$stocked_out==FALSE)
  dem_T <- eval_lambda_new(deltain, theta1, wdcMergedday, points, tw_groupin)[stocked_list]  
  return(dem_T)
}


