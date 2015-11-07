
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
