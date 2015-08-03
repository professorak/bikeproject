sourceCpp("eval_func_2_new_deltaaveraged.cpp")

eval_hessian_delta_list_new <- function(deltain_tw, theta1, wdcMergedday, points, tw_groupin,lambda_multiplers_in) {
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
                                        lambda_multiplers_in)
  
  return(hessian_lambda_delta_sq)  
  
}

deltain <- rep(-30, nrow(wdcMerged))
deltain[which(wdcMerged$stocked_out==F)] <- deltain_stin
theta1 <- c(theta[1],0,theta[-1])
stocked_list <- which(wdcMerged$stocked_out==F)

tw_group_list <- unique(wdcMerged$tw_group)
i<- 1
tw_groupin = tw_group_list[i]
wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
deltain_tw = deltain[which(wdcMerged$tw_group==tw_groupin)]
ret <- eval_hessian_delta_list_new(deltain_tw, theta1, wdcMergedday, points, tw_groupin)
summary(c(ret))

################################################################
################################################################


compute_hess_atindex <- function(constraint_index, delta1_index, delta2_index, lambda_multiplers_in) {
  a1 <- eval_lambda_new(deltain_tw,theta1, wdcMergedday, points, tw_groupin)[constraint_index]
  #a1_hess <- eval_hessian_delta_list_new(deltain_tw, theta1, wdcMergedday, points, tw_groupin,lambda_multiplers_in)
  summary(c(a1_hess))
  
  deltain_tw_2 <- deltain_tw
  deltain_tw_2[delta1_index] <- deltain_tw_2[delta1_index] + diff
  
  #gradient wrt delta1_index at base value
  a1_2 <- eval_lambda_new(deltain_tw_2,theta1, wdcMergedday, points, tw_groupin)[constraint_index]
  a1_grad <- (a1_2-a1)/diff
  
  #gradient wrt delta1_index at delta2_index shifted value
  deltain_tw_3 <- deltain_tw
  deltain_tw_3[delta2_index] <- deltain_tw_3[delta2_index] + diff
  
  deltain_tw_4 <- deltain_tw_3
  deltain_tw_4[delta1_index] <- deltain_tw_4[delta1_index] + diff
  
  
  a2 <- eval_lambda_new(deltain_tw_3,theta1, wdcMergedday, points, tw_groupin)[constraint_index]
  #a2_hess <- eval_hessian_delta_list_new(deltain_tw_3, theta1, wdcMergedday, points, tw_groupin,lambda_multiplers_in)
  
  a2_2 <- eval_lambda_new(deltain_tw_4,theta1, wdcMergedday, points, tw_groupin)[constraint_index]
  a2_grad <- (a2_2-a2)/diff
  
  
  ###
  a1_hess_num <- (a2_grad-a1_grad)/diff
  return(a1_hess_num)
  
}

#test hessian 
constraint_index <- 2
delta1_index <- 13
delta2_index <- 13
diff <- 0.001
lambda_multiplers_in <- rep(0,nrow(wdcMergedday))
lambda_multiplers_in[constraint_index] <- 1
a1_hess <- eval_hessian_delta_list_new(deltain_tw, theta1, wdcMergedday, points, tw_groupin,lambda_multiplers_in)
compute_hess_atindex(constraint_index, delta1_index, delta2_index, lambda_multiplers_in)/wdcMergedday$obs_weight[delta1_index]/wdcMergedday$obs_weight[delta2_index]

a1_hess[delta1_index,delta2_index]/wdcMergedday$obs_weight[delta1_index]/wdcMergedday$obs_weight[delta2_index]
#a2_hess[delta1_index,delta2_index]

#compute entire numerical hessian

lambda_multiplers_in <- rep(0,nrow(wdcMergedday))
lambda_multiplers_in[constraint_index] <- 1
constraint_index <- 2
diff <- 0.001
a1_hess_num <- matrix(0,nrow(wdcMergedday),nrow(wdcMergedday))
for(i in 1:nrow(wdcMergedday)) {
  for(j in 1:i) {
    delta1_index <- i
    delta2_index <- j
    a1_hess_num[i,j] <- compute_hess_atindex(constraint_index, delta1_index, delta2_index, lambda_multiplers_in)
    a1_hess_num[j,i] <- a1_hess_num[i,j]
  }
}

idx <- which(a1_hess_num!=0)
a <- cbind(a1_hess_num[idx], a1_hess[idx])
summary(a1_hess_num[idx]-a1_hess[idx])
# #test lambda
# a1 <- eval_lambda_new(deltain_tw,theta1, wdcMergedday, points, tw_groupin)[constraint_index]
# deltain_tw_2 <- deltain_tw
# deltain_tw_2[9] <- deltain_tw_2[9] + 1/wdcMergedday$obs_weight[9]
# a2 <- eval_lambda_new(deltain_tw_2,theta1, wdcMergedday, points, tw_groupin)[constraint_index]
# deltain_tw_3 <- deltain_tw
# deltain_tw_3[10] <- deltain_tw_3[10] + 1/wdcMergedday$obs_weight[10]
# a3 <- eval_lambda_new(deltain_tw_3,theta1, wdcMergedday, points, tw_groupin)[constraint_index]

