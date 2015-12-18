removemodify_outlier_states_keep_majoritystates <- function(wdcMerged) {
    
  #keeping obervations so that only enough states are kept to have 80% of the observations 
  print("Keeping states that make up 80% of observation weight in stXtw_group")
  
  perc_cumsum_threshold <- 0.8
  
  #for each station_id_index X tw_group.
  #order by obs_count
  
  wdcMerged <- wdcMerged[order(wdcMerged$tw_group,wdcMerged$station_id_index, wdcMerged$stocked_out,-wdcMerged$obs_weight),]
  wdcMerged$tot_obs_weight <- 
    ave(wdcMerged$obs_weight, wdcMerged$tw_group, wdcMerged$station_id_index, FUN=sum)
  wdcMerged$cumsum_obs_weight <- 
    ave(wdcMerged$obs_weight, wdcMerged$tw_group, wdcMerged$station_id_index, FUN=cumsum)
  wdcMerged$cumsum_percentile_obs_weight <- wdcMerged$cumsum_obs_weight/wdcMerged$tot_obs_weight
  wdcMerged$lag_cumsum_percentile_obs_weight <- 
    ave(wdcMerged$cumsum_percentile_obs_weight, wdcMerged$tw_group, wdcMerged$station_id_index, 
        FUN=function(x){return(c(0,x[-length(x)]))})
  
  wdcMerged <- subset(wdcMerged, lag_cumsum_percentile_obs_weight<perc_cumsum_threshold)
  
  
  
  #remove too low wieghts states
  print("Removing too low wieghts states:")
  print(paste0("Percentage of total observations removed: ",
               round(sum(wdcMerged$obs_weight[which(wdcMerged$obs_weight<=10 & wdcMerged$index>1)])/sum(wdcMerged$obs_weight)*100,3)," %"))
  wdcMerged <- subset(wdcMerged, obs_weight>10 | index==1)
  wdcMerged <- droplevels(wdcMerged)
  
  #truncating really low demand
  print("Removing too low demand states:")
  print(paste0("Percentage of 0 demand states: ",
               round(sum(wdcMerged$obs_weight[which(wdcMerged$out_dem_mean<=0 & 
                                                      wdcMerged$stocked_out==FALSE)])/sum(wdcMerged$obs_weight)*100,3)," %"))
  
  print(paste0("Percentage of too low demand states truncated at 1e-5: ",
               round(sum(wdcMerged$obs_weight[which(wdcMerged$out_dem_mean<=1e-5 & 
                                                      wdcMerged$stocked_out==FALSE)])/sum(wdcMerged$obs_weight)*100,3)," %"))
  wdcMerged$out_dem_mean[which(wdcMerged$out_dem_mean<=1e-5 & 
                                 wdcMerged$stocked_out==FALSE)] <- 1e-5
  
  #truncating really high demand
  print(paste0("Percentage of really high demand states truncated at 2 : ",
               round(sum(wdcMerged$obs_weight[which(wdcMerged$out_dem_mean>2 & 
                                                      wdcMerged$stocked_out==FALSE)])/sum(wdcMerged$obs_weight)*100,5)," %"))
  
  wdcMerged$out_dem_mean[which(wdcMerged$out_dem_mean>2 & 
                                 wdcMerged$stocked_out==FALSE)] <- 2
  
  wdcMerged$out_dem_sum <- wdcMerged$out_dem_mean * wdcMerged$obs_weight
  
  
  #check all tw's have all stations
  tw_group_list <- unique(wdcMerged$tw_group)
  no_st <- max(wdcMerged$station_id_index)
  for(i in 1:length(tw_group_list)) {
    tw_groupin = tw_group_list[i]
    idx <- which(wdcMerged$tw_group==tw_groupin)
    if(length(unique(wdcMerged$station_id_index[idx]))!=no_st) {
      print(paste0("in tw:", tw_groupin, "all stations are not present"))
    }
  }  
  return(wdcMerged)
}

removemodify_outlier_states_keep_majoritystates_no0dem <- function(wdcMerged) {
    
  #keeping obervations so that only enough states are kept to have 80% of the observations 
  print("Keeping states that make up 80% of observation weight in stXtw_group")
  
  perc_cumsum_threshold <- 0.8
  
  #for each station_id_index X tw_group.
  #order by obs_count
  
  wdcMerged <- wdcMerged[order(wdcMerged$tw_group,wdcMerged$station_id_index, wdcMerged$stocked_out,-wdcMerged$obs_weight),]
  wdcMerged$tot_obs_weight <- 
    ave(wdcMerged$obs_weight, wdcMerged$tw_group, wdcMerged$station_id_index, FUN=sum)
  wdcMerged$cumsum_obs_weight <- 
    ave(wdcMerged$obs_weight, wdcMerged$tw_group, wdcMerged$station_id_index, FUN=cumsum)
  wdcMerged$cumsum_percentile_obs_weight <- wdcMerged$cumsum_obs_weight/wdcMerged$tot_obs_weight
  wdcMerged$lag_cumsum_percentile_obs_weight <- 
    ave(wdcMerged$cumsum_percentile_obs_weight, wdcMerged$tw_group, wdcMerged$station_id_index, 
        FUN=function(x){return(c(0,x[-length(x)]))})
  
  wdcMerged <- subset(wdcMerged, lag_cumsum_percentile_obs_weight<perc_cumsum_threshold)
  
  
  
  #remove too low wieghts states
  print("Removing too low wieghts states:")
  print(paste0("Percentage of total observations removed: ",
               round(sum(wdcMerged$obs_weight[which(wdcMerged$obs_weight<=10 & wdcMerged$index>1)])/sum(wdcMerged$obs_weight)*100,3)," %"))
  wdcMerged <- subset(wdcMerged, obs_weight>10 | index==1)
  wdcMerged <- droplevels(wdcMerged)
  
  #truncating really low demand
  print("Removing too low demand states:")
  print(paste0("Percentage of 0 demand states: ",
               round(sum(wdcMerged$obs_weight[which(wdcMerged$out_dem_mean<=0 & 
                                                      wdcMerged$stocked_out==FALSE)])/sum(wdcMerged$obs_weight)*100,3)," %"))
  wdcMerged <- subset(wdcMerged, out_dem_mean>0 | index==1 | wdcMerged$stocked_out==T)
  
  print(paste0("Percentage of too low demand states truncated at ",
  round(min(wdcMerged$out_dem_mean[which(wdcMerged$out_dem_mean>0 & 
  wdcMerged$stocked_out==FALSE)]),5)
  , ": ",
  round(sum(wdcMerged$obs_weight[which(wdcMerged$out_dem_mean<=0 & 
  wdcMerged$stocked_out==FALSE)])/sum(wdcMerged$obs_weight)*100,3)," %"))
  
  wdcMerged$out_dem_mean[which(wdcMerged$out_dem_mean<=0 & 
                                 wdcMerged$stocked_out==FALSE)] <- 
    min(wdcMerged$out_dem_mean[which(wdcMerged$out_dem_mean>0 & 
                                       wdcMerged$stocked_out==FALSE)])
  
  #truncating really high demand
  print(paste0("Percentage of really high demand states truncated at 2 : ",
               round(sum(wdcMerged$obs_weight[which(wdcMerged$out_dem_mean>2 & 
                                                      wdcMerged$stocked_out==FALSE)])/sum(wdcMerged$obs_weight)*100,5)," %"))
  
  wdcMerged$out_dem_mean[which(wdcMerged$out_dem_mean>2 & 
                                 wdcMerged$stocked_out==FALSE)] <- 2
  
  wdcMerged$out_dem_sum <- wdcMerged$out_dem_mean * wdcMerged$obs_weight
  
  
  #check all tw's have all stations
  tw_group_list <- unique(wdcMerged$tw_group)
  no_st <- max(wdcMerged$station_id_index)
  for(i in 1:length(tw_group_list)) {
    tw_groupin = tw_group_list[i]
    idx <- which(wdcMerged$tw_group==tw_groupin)
    if(length(unique(wdcMerged$station_id_index[idx]))!=no_st) {
      print(paste0("in tw:", tw_groupin, "all stations are not present"))
    }
  }  
  return(wdcMerged)
}

removemodify_outlier_states <- function(wdcMerged) {
  #remove too low wieghts states
  print("Removing too low wieghts states:")
  print(paste0("Percentage of total observations removed: ",
               round(sum(wdcMerged$obs_weight[which(wdcMerged$obs_weight<=10 & wdcMerged$index>1)])/sum(wdcMerged$obs_weight)*100,3)," %"))
  wdcMerged <- subset(wdcMerged, obs_weight>10 | index==1)
  wdcMerged <- droplevels(wdcMerged)
  
  #truncating really low demand
  print("Removing too low demand states:")
  print(paste0("Percentage of 0 demand states: ",
               round(sum(wdcMerged$obs_weight[which(wdcMerged$out_dem_mean<=0 & 
                                                      wdcMerged$stocked_out==FALSE)])/sum(wdcMerged$obs_weight)*100,3)," %"))
  
  print(paste0("Percentage of too low demand states truncated at 1e-5: ",
               round(sum(wdcMerged$obs_weight[which(wdcMerged$out_dem_mean<=1e-5 & 
                                                      wdcMerged$stocked_out==FALSE)])/sum(wdcMerged$obs_weight)*100,3)," %"))
  wdcMerged$out_dem_mean[which(wdcMerged$out_dem_mean<=1e-5 & 
                                 wdcMerged$stocked_out==FALSE)] <- 1e-5
  
  #truncating really high demand
  print(paste0("Percentage of really high demand states truncated at 2 : ",
               round(sum(wdcMerged$obs_weight[which(wdcMerged$out_dem_mean>2 & 
                                                      wdcMerged$stocked_out==FALSE)])/sum(wdcMerged$obs_weight)*100,5)," %"))
  
  wdcMerged$out_dem_mean[which(wdcMerged$out_dem_mean>2 & 
                                 wdcMerged$stocked_out==FALSE)] <- 2
  
  wdcMerged$out_dem_sum <- wdcMerged$out_dem_mean * wdcMerged$obs_weight
  
  
  #check all tw's have all stations
  tw_group_list <- unique(wdcMerged$tw_group)
  no_st <- max(wdcMerged$station_id_index)
  for(i in 1:length(tw_group_list)) {
    tw_groupin = tw_group_list[i]
    idx <- which(wdcMerged$tw_group==tw_groupin)
    if(length(unique(wdcMerged$station_id_index[idx]))!=no_st) {
      print(paste0("in tw:", tw_groupin, "all stations are not present"))
    }
  }  
  return(wdcMerged)
}

get_servicelevel_instruments <- function() {
  #get local_stations
  local_stations_vec <- wdcMerged[!duplicated(wdcMerged$station_id_index),c("station_id_index","local_stations")]
  local_stations_vec <- local_stations_vec[order(local_stations_vec$station_id_index),]

  tw_list <- unique(user_serv_lvl$tw)
  instr_serv_lvl_1_df <- c()
  no_st <- max(user_serv_lvl$station_id_index)  
  for(tw_in in tw_list) {
    user_serv_lvl_subset <- user_serv_lvl[which(user_serv_lvl$tw==tw_in),]
    #if(!identical(order(user_serv_lvl_subset$station_id_index),c(1:nrow(user_serv_lvl_subset)))) stop("not in order")
    user_serv_lvl_subset <- user_serv_lvl_subset[order(user_serv_lvl_subset$station_id_index),]
    if(nrow(user_serv_lvl_subset)!=no_st) stop("incorrect number of rows in user_serv_lvl_subset")
    instr_serv_lvl_1_tw <- gen_serv_lvl_instr_neigh(user_serv_lvl_subset$serv_lvl,local_stations_vec$local_stations)
    instr_serv_lvl_1_tw_df <- data.frame(tw=user_serv_lvl_subset$tw, station_id_index=user_serv_lvl_subset$station_id_index,
                                         st_tw_index=user_serv_lvl_subset$st_tw_index,instr_serv_lvl=instr_serv_lvl_1_tw)
    instr_serv_lvl_1_df <- rbind(instr_serv_lvl_1_df, instr_serv_lvl_1_tw_df)
  }
  instr_serv_lvl_1_df <- instr_serv_lvl_1_df[order(instr_serv_lvl_1_df$tw,
                                                   instr_serv_lvl_1_df$station_id_index),]
  if(!identical(user_serv_lvl$station_id_index, instr_serv_lvl_1_df$station_id_index)) stop("user_serv_lvl and instr_serv_lvl_1_df not in order")
  if(!identical(user_serv_lvl$tw, instr_serv_lvl_1_df$tw)) stop("user_serv_lvl and instr_serv_lvl_1_df not in order")
  user_serv_lvl$serv_lvl_neighbours <<- instr_serv_lvl_1_df$instr_serv_lvl

  #add serv_lvl instruments to wdcMerged
  if(!identical(user_serv_lvl$st_tw_index, c(1:nrow(user_serv_lvl)))) stop("user_serv_lvl not in order")
  wdcMerged$instr_serv_lvl <<- user_serv_lvl$instr_serv_lvl[wdcMerged$st_tw_index]
  wdcMerged$serv_lvl_neighbours <<- user_serv_lvl$serv_lvl_neighbours[wdcMerged$st_tw_index]
  wdcMerged$instr_serv_lvl_sq <<- wdcMerged$instr_serv_lvl^2
  wdcMerged$serv_lvl_neighbours_sq <<- wdcMerged$serv_lvl_neighbours^2
  
  rm(instr_serv_lvl_1_df)
  rm(local_stations_vec)
  rm(user_serv_lvl_subset)
  rm(instr_serv_lvl_1_tw_df)
}

get_servicelevel_instruments_more <- function() {
  get_servicelevel_instruments()
  #add serv_lvl instruments to wdcMerged
  if(!identical(user_serv_lvl$st_tw_index, c(1:nrow(user_serv_lvl)))) stop("user_serv_lvl not in order")
  wdcMerged$instr_serv_lvl <<- user_serv_lvl$instr_serv_lvl[wdcMerged$st_tw_index]

  wdcMerged$in_dem_rate <<- user_serv_lvl$in_dem_rate[wdcMerged$st_tw_index]
  wdcMerged$out_dem_rate <<- user_serv_lvl$out_dem_rate[wdcMerged$st_tw_index]
  wdcMerged$in_dem_rate_lagtw <<- user_serv_lvl$in_dem_rate_lagtw[wdcMerged$st_tw_index]
  wdcMerged$out_dem_rate_lagtw <<- user_serv_lvl$out_dem_rate_lagtw[wdcMerged$st_tw_index]
  wdcMerged$diff_dem_rate_lagtw <<- user_serv_lvl$diff_dem_rate_lagtw[wdcMerged$st_tw_index]
  wdcMerged$out_dem_rate_lagtw <<- user_serv_lvl$out_dem_rate_lagtw[wdcMerged$st_tw_index]
  wdcMerged$diffdummy_dem_rate_lagtw <<- user_serv_lvl$diffdummy_dem_rate_lagtw[wdcMerged$st_tw_index]

  wdcMerged$in_dem_rate_sq <<- wdcMerged$in_dem_rate^2
  wdcMerged$in_dem_rate_lagtw_sq <<- wdcMerged$in_dem_rate_lagtw^2
  wdcMerged$diff_dem_rate_lagtw_sq <<- wdcMerged$diff_dem_rate_lagtw^2

}

