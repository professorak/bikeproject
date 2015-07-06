#check the theta gradient of random coef for a smaller instance
length_theta <- length(theta)

#theta1 <- c(theta[1],0,theta[-1])  
theta1 <- theta

#expand deltain to all observations, it is currently #stocked in observations
deltain <- rep(-30, nrow(wdcMerged))
deltain[which(wdcMerged$stocked_out==F)] <- deltain_stin

tw_group_list <- unique(wdcMerged$tw_group)
grad_constraints <- c()
i <- 1
tw_groupin = tw_group_list[i]
wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
deltain_tw = deltain[which(wdcMerged$tw_group==tw_groupin)]

grad_theta <- eval_grad_lambda_theta_new(theta1, deltain_tw, wdcMergedday, points, tw_groupin)

###temp
params_save <- params
theta_save <- params_save[c(1:length(theta))]
deltain_stin_save <- params_save[-c(1:length(theta))]

#compute density constraint
a1 <- eval_g(params,wdcMerged, points, length(theta))
a1[length_stklist+1]


#make rd to 0 and compute
params_2 <- params_save
params_2[2] <- 0
a2 <- eval_g(params_2,wdcMerged, points, length(theta))
a2[length_stklist+1]


#compute predicted lambda values to see if subsitution precentages computed in 
#counterfactuals are correct
theta1 <- c(-12.28769991999, 4.20013541466, 0.09999999003, 282.66078886268, 8.10824517514, 7505.59969617232)
delta_all <- compute_delta_list_cntrt_map_new(theta1,wdcMerged, points)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]

#test 1
tw_group_list <- unique(wdcMerged$tw_group)
lambda_t <- c()
for(i in 1:length(tw_group_list)) {
  tw_groupin = tw_group_list[i]
  wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
  deltain_tw = delta_all[which(wdcMerged$tw_group==tw_groupin)]    
  lambda_t <- c(lambda_t, 
                eval_lambda_delta_list_new(deltain_tw, theta1, wdcMergedday, 
                                            points, tw_groupin)[[1]])
}


test_dt <- wdcMerged[,c("station_id_index","tw_group","sto_state_local",
                        "weather_state","obs_weight")]
test_dt$lambda_t <- lambda_t
test_dt$mean_lambda_t <- test_dt$lambda_t/test_dt$obs_weight




calculate_subs_perc <- function(theta1, delta_all, wdcMerged, points) {
  #creating a data set where there is one entry for each station_tw_group and 
  #will manipulate stockout state and compute the substitution percewntage
  cf_data <- wdcMerged[,c("station_id_index","tw_group","sto_state_local","local_stations","obs_weight",
                          "st_tw","st_tw_index", "station_id",
                          "stocked_out","lat","lon","out_dem_sum")]
  #keep 1 entry for each st_twgroup.
  cf_data$st_twgroup <- paste0(cf_data$station_id_index,"_",cf_data$tw_group)
  delta_all_weighted <- delta_all*cf_data$obs_weight
  delta_sum <- by(delta_all_weighted, cf_data$st_twgroup, FUN=sum)
  weighted_sum <- by(cf_data$obs_weight, cf_data$st_twgroup, FUN=sum)
  delta_ave_df <- data.frame(delta_sum=as.numeric(delta_sum), st_twgroup=(names(delta_sum)), weighted_sum=as.numeric(weighted_sum))
  delta_ave_df$st_twgroup <- as.character(delta_ave_df$st_twgroup)
  delta_ave_df$delta_ave <- delta_ave_df$delta_sum/delta_ave_df$weighted_sum
  cf_data <- cf_data[which(!duplicated(cf_data$st_twgroup)),]
  
  delta_ave_df <- delta_ave_df[order(order(cf_data$st_twgroup)),]
  if(!identical(cf_data$st_twgroup, delta_ave_df$st_twgroup)) stop("st_twgroup in cf_data and delta_ave_df not equal")
  
  perc_substituted <- c()
  focal_station_index <- 1    
  for(focal_station_index in c(1:max(cf_data$station_id_index))) {
    print(focal_station_index)
    lambda1 <- calculate_subs_perc_focal_station(focal_station_index, 0, theta1, cf_data, points,delta_ave_df$delta_ave)
    lambda2 <- calculate_subs_perc_focal_station(focal_station_index, 1, theta1, cf_data, points,delta_ave_df$delta_ave)
    focal_st_demand <- sum(lambda1[which(cf_data$station_id_index==focal_station_index)])    
    substituted_demand <- sum(lambda2) - sum(lambda1) + focal_st_demand
    perc_substituted <- c(perc_substituted,(substituted_demand/focal_st_demand*100))
    print(substituted_demand/focal_st_demand*100)
  }
  return(perc_substituted)
    
}

calculate_subs_perc_focal_station <- function(focal_station_index, mode, theta1, cf_data, points,delta_ave) {
  so_vec_cf <- rep(0,max(cf_data$station_id_index))
  cf_data$stocked_out = FALSE
  if(mode==1) {
    so_vec_cf[focal_station_index] <- 1    
    cf_data$stocked_out[which(cf_data$station_id_index==focal_station_index)] = TRUE
  }
  cf_data$sto_state_local <- gen_sto_state_local_char(so_vec_cf, cf_data$local_stations)
  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  points_sub <- which(stid_in_localstations(focal_station_index, points$local_stations))
  points_sub <- points[points_sub,]   
  lambda_t <- c()
  for(tw_groupin in tw_groupin_list) {
    cf_dataday = subset(cf_data, tw_group==tw_groupin)    
    cf_dataday$obs_weight <- 1
    deltain_tw <- delta_ave[which(cf_data$tw_group==tw_groupin)]  
    lambda_t <- c(lambda_t,
                  eval_lambda_delta_list_new(deltain_tw, theta1,cf_dataday, points_sub, tw_groupin)$objective)      
  }
  return((lambda_t))
}


#test 3
theta1 <- c(-12.28769991999, 4.20013541466, 0.09999999003, 282.66078886268, 8.10824517514, 7505.59969617232)
delta_all <- compute_delta_list_cntrt_map_new(theta1,wdcMerged, points)

theta1_1 <- theta1
subs_perc_1 <- calculate_subs_perc(theta1_1,delta_all,wdcMerged,points)

theta1_2 <- theta1
theta1_2[2] <- 0 #rd coef = 0
subs_perc_2 <- calculate_subs_perc(theta1_2,delta_all,wdcMerged,points)

theta1_3 <- theta1
theta1_3[2] <- 0 #rd coef = 0
theta1_3[1] <- -4 #rd coef = 0
subs_perc_3 <- calculate_subs_perc(theta1_3,delta_all,wdcMerged,points)

theta1_4 <- theta1
theta1_4[2] <- 8 #rd coef = 0
theta1_4[1] <- -4 #rd coef = 0
subs_perc_4 <- calculate_subs_perc(theta1_4,delta_all,wdcMerged,points)

calculate_subs_perc_alldata <- function(theta1, delta_all, wdcMerged, points) {
  #creating a data set where there is one entry for each station_tw_group and 
  #will manipulate stockout state and compute the substitution percewntage
  cf_data <- wdcMerged[,c("station_id_index","tw_group","sto_state_local","local_stations","obs_weight",
                          "st_tw","st_tw_index", "station_id",
                          "stocked_out","lat","lon","out_dem_sum")]
##   keep 1 entry for each st_twgroup.
  cf_data$st_twgroup <- paste0(cf_data$station_id_index,"_",cf_data$tw_group)
#   delta_all_weighted <- delta_all*cf_data$obs_weight
#   delta_sum <- by(delta_all_weighted, cf_data$st_twgroup, FUN=sum)
#   weighted_sum <- by(cf_data$obs_weight, cf_data$st_twgroup, FUN=sum)
#   delta_ave_df <- data.frame(delta_sum=as.numeric(delta_sum), st_twgroup=(names(delta_sum)), weighted_sum=as.numeric(weighted_sum))
#   delta_ave_df$st_twgroup <- as.character(delta_ave_df$st_twgroup)
#   delta_ave_df$delta_ave <- delta_ave_df$delta_sum/delta_ave_df$weighted_sum
#   cf_data <- cf_data[which(!duplicated(cf_data$st_twgroup)),]
#   
#   delta_ave_df <- delta_ave_df[order(order(cf_data$st_twgroup)),]
#   if(!identical(cf_data$st_twgroup, delta_ave_df$st_twgroup)) stop("st_twgroup in cf_data and delta_ave_df not equal")
  
  perc_substituted <- c()
  focal_station_index <- 1    
  for(focal_station_index in c(1:max(cf_data$station_id_index))) {
    print(focal_station_index)
    lambda1 <- calculate_subs_perc_focal_station_alldata(focal_station_index, 0, theta1, cf_data, points,delta_all)    
    lambda1_dt <- getMeanDemandDf(lambda1, cf_data)
    lambda2 <- calculate_subs_perc_focal_station_alldata(focal_station_index, 1, theta1, cf_data, points,delta_all)
    lambda2_dt <- getMeanDemandDf(lambda2, cf_data)
    focal_st_demand <- sum(lambda1_dt$lambda_mean[which(lambda1_dt$station_id_index==focal_station_index)])    
    substituted_demand <- sum(lambda2_dt$lambda_mean) - sum(lambda1_dt$lambda_mean) + focal_st_demand
    perc_substituted <- c(perc_substituted,(substituted_demand/focal_st_demand*100))
    print(substituted_demand/focal_st_demand*100)
  }
  return(perc_substituted)
  
}

calculate_subs_perc_focal_station_alldata <- function(focal_station_index, mode, theta1, cf_data, points,delta_ave) {
  so_vec_cf <- rep(0,max(cf_data$station_id_index))
  cf_data$stocked_out = FALSE
  if(mode==1) {
    so_vec_cf[focal_station_index] <- 1    
    cf_data$stocked_out[which(cf_data$station_id_index==focal_station_index)] = TRUE
  }
  cf_data$sto_state_local <- gen_sto_state_local_char(so_vec_cf, cf_data$local_stations)
  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  points_sub <- which(stid_in_localstations(focal_station_index, points$local_stations))
  points_sub <- points[points_sub,]   
  lambda_t <- c()
  for(tw_groupin in tw_groupin_list) {
    cf_dataday = subset(cf_data, tw_group==tw_groupin)    
    deltain_tw <- delta_ave[which(cf_data$tw_group==tw_groupin)]  
    lambda_t <- c(lambda_t,
                  eval_lambda_delta_list_new(deltain_tw, theta1,cf_dataday, points_sub, tw_groupin)$objective)      
  }
  return((lambda_t))
}

getMeanDemandDf <- function(lambda1, cf_data) {
  #for each st_twgroup compute the total wighted lambdas
  #compute the total wieght of observations 
  #divide the two to compute average lambda and return in df, along with corresponding station_id_index of each
  cf_data_dt <- data.table(cf_data)  
  cf_data_dt$lambda <- lambda1
  lambda_dt <- cf_data_dt[,list(lambda_sum=sum(lambda),weight_sum=sum(obs_weight),
                                station_id_index=station_id_index[1]), by=c("st_twgroup")]
  lambda_dt$lambda_mean <- lambda_dt$lambda_sum/lambda_dt$weight_sum
  return(lambda_dt)
}


#test 3
theta1 <- c(-12.28769991999, 4.20013541466, 0.09999999003, 282.66078886268, 8.10824517514, 7505.59969617232)
delta_all <- compute_delta_list_cntrt_map_new(theta1,wdcMerged, points)

theta1_1 <- theta1
subs_perc_1 <- calculate_subs_perc_alldata(theta1_1,delta_all,wdcMerged,points)

theta1_2 <- theta1
theta1_2[2] <- 0 #rd coef = 0
subs_perc_2 <- calculate_subs_perc_alldata(theta1_2,delta_all,wdcMerged,points)

theta1_3 <- theta1
theta1_3[2] <- 0 #rd coef = 0
theta1_3[1] <- -4 #rd coef = 0
subs_perc_3 <- calculate_subs_perc_alldata(theta1_3,delta_all,wdcMerged,points)

theta1_4 <- theta1
theta1_4[2] <- 8 #rd coef = 0
theta1_4[1] <- -4 #rd coef = 0
subs_perc_4 <- calculate_subs_perc_alldata(theta1_4,delta_all,wdcMerged,points)
