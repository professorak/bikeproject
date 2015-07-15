demand_distance_change_inner <- function(scale_distance_factor, theta1, delta_list, wdcMerged, points) {
  theta1[1] <- theta1[1]*scale_distance_factor
  tw_group_list <- c()
  if(length(delta_list)!=nrow(wdcMerged)) stop("demand_distance_change_inner error")
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

demand_distance_change <- function(scale_distance_factor, theta1, delta_list, wdcMerged, points) {
  demand_actual <- demand_distance_change_inner(1, theta1, delta_list, wdcMerged, points)
  demand_changed <- demand_distance_change_inner(scale_distance_factor, theta1, delta_list, wdcMerged, points)
  return(c(demand_changed,demand_actual))
}  

demand_distance_change_inner_serv_lvl_adjusted <- function(scale_distance_factor, theta1, delta_list, wdcMerged, points, 
                                                               serv_lvl_covar_new, serv_lvl_covar) {
  theta1[1] <- theta1[1]*scale_distance_factor
  tw_group_list <- c()
  if(length(delta_list)!=nrow(wdcMerged)) stop("demand_distance_change_inner_serv_lvl_adjusted error")
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

demand_serv_lvl_change <- function(theta1,deltain, wdcMerged,points, serv_lvl_factor,serv_lvl_coef) {
  demand_actual <- demand_distance_change_inner(1, theta1, deltain, wdcMerged, points)
  
  tw_list <- unique(wdcMerged$tw)    
  no_st <- max(wdcMerged$station_id_index)
  serv_lvl_covar <- c()
  for(tw_in in tw_list) {
    wdcMergedday = subset(wdcMerged, tw==tw_in)      
    deltain_tw <- deltain[which(wdcMerged$tw==tw_in)]
    user_serv_lvl_subset <- user_serv_lvl[which(user_serv_lvl$tw==tw_in),]
    if(length(which(user_serv_lvl_subset$station_id_index!=
                  c(1:nrow(user_serv_lvl_subset))))) stop("demand_serv_lvl_change error1")
    serv_lvl_covar <- c(serv_lvl_covar, 
                        user_serv_lvl_subset$serv_lvl[wdcMergedday$station_id_index])
  }
  
  serv_lvl_covar_new <- serv_lvl_covar*serv_lvl_factor
#   serv_lvl_covar_new[which(serv_lvl_covar_new>1)] <- 1
#   serv_lvl_covar_new[which(serv_lvl_covar_new<0)] <- 0
  deltain_new <- (serv_lvl_covar_new-serv_lvl_covar)*serv_lvl_coef + deltain
  
  tiny	<- 1e-4
  demand_new <- demand_distance_change_inner_serv_lvl_adjusted(1, theta1, deltain_new, wdcMerged, points, 
                                                                   serv_lvl_covar_new+tiny, serv_lvl_covar+tiny)
  return(c(demand_new,demand_actual))
}

