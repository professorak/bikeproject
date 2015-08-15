eval_constraints_MPEC_tw_groupin <- function(deltain, theta1, wdcMergedday, points, tw_groupin) {
  stocked_list <- which(wdcMergedday$stocked_out==FALSE)
  dem_T <- eval_lambda_new(deltain, theta1, wdcMergedday, points, tw_groupin)[stocked_list]
  dem_hat_T <- wdcMergedday$out_dem_sum[stocked_list]
  obj <- dem_T-dem_hat_T
  return(obj)
}

eval_constraints_MPEC_tw_groupin_tiny <- function(deltain, theta1, wdcMergedday, points, tw_groupin) {
  stocked_list <- which(wdcMergedday$stocked_out==FALSE)
  dem_T <- eval_lambda_new(deltain, theta1, wdcMergedday, points, tw_groupin)[stocked_list]
  dem_hat_T <- wdcMergedday$out_dem_sum[stocked_list]
  obj <- (dem_T-dem_hat_T)/wdcMergedday$obs_weight[stocked_list]
  return(obj)
}


eval_constraints_MPEC <- function(params, wdcMerged, points, length_theta) {
  theta <- params[c(1:length_theta)]
  deltain_stin <- params[-c(1:length_theta)]
  theta1 <- c(theta[1],0,theta[-1])  
  
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
      eval_constraints_MPEC_tw_groupin(deltain_tw, theta1, wdcMergedday, 
                                       points, tw_groupin))
  }
  return(obj)
}

eval_constraints_MPEC_tiny <- function(params, wdcMerged, points, length_theta) {
  theta <- params[c(1:length_theta)]
  deltain_stin <- params[-c(1:length_theta)]
  theta1 <- c(theta[1],0,theta[-1])  
  
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
             eval_constraints_MPEC_tw_groupin_tiny(deltain_tw, theta1, wdcMergedday, 
                                              points, tw_groupin))
  }
  return(obj)
}

eval_g <- function(params, wdcMerged, points, length_theta) {
  if(print_iter_values) {print("In eval_g: ") }
  ptm <- proc.time()  
  ret <- c(eval_constraints_MPEC(params, wdcMerged, points, length_theta),
           get_total_density(params, wdcMerged, points))
  if(print_iter_values) {
    print("eval time:")
    print(proc.time()-ptm)
  }
  return(ret)  
}

eval_g_tiny <- function(params, wdcMerged, points, length_theta) {
  if(print_iter_values) {
    print("In eval_g: ")
  }
  ptm <- proc.time()  
  ret <- c(eval_constraints_MPEC_tiny(params, wdcMerged, points, length_theta),
           get_total_density(params, wdcMerged, points)/tot_density*target_density_constraint)
  if(print_iter_values) {
    print("eval time:")
    print(proc.time()-ptm)
  }
  return(ret)  
}

eval_grad_constraints_MPEC_tw_groupin <- function(deltain, theta1, wdcMergedday, points, tw_groupin) {
  grad_theta <- eval_grad_lambda_theta_new(theta1, deltain, wdcMergedday, points, tw_groupin)
  grad_theta <- grad_theta[,-2]
  grad_delta <- eval_lambda_delta_list_new(deltain, theta1, wdcMergedday, points, tw_groupin)[[2]]

  #reduce the rows and columns corresponding to only stocked in obervations of wdcMergedday
  stocked_list <- which(wdcMergedday$stocked_out==F)
  grad_theta <- grad_theta[stocked_list,]
  grad_delta <- grad_delta[stocked_list, stocked_list]
  grad_constraints_tw <- c(t(cbind(grad_theta,grad_delta)))  
  return(grad_constraints_tw[which(grad_constraints_tw!=0)])
}

eval_grad_constraints_MPEC_tw_groupin_tiny <- function(deltain, theta1, wdcMergedday, points, tw_groupin) {
  grad_theta <- eval_grad_lambda_theta_new(theta1, deltain, wdcMergedday, points, tw_groupin)
  grad_theta <- grad_theta[,-2]
  grad_delta <- eval_lambda_delta_list_new(deltain, theta1, wdcMergedday, points, tw_groupin)[[2]]
  
  #reduce the rows and columns corresponding to only stocked in obervations of wdcMergedday
  stocked_list <- which(wdcMergedday$stocked_out==F)
  grad_theta <- grad_theta[stocked_list,]/wdcMergedday$obs_weight[stocked_list]
  grad_delta <- grad_delta[stocked_list, stocked_list]/wdcMergedday$obs_weight[stocked_list]
  grad_constraints_tw <- c(t(cbind(grad_theta,grad_delta)))  
  return(grad_constraints_tw[which(grad_constraints_tw!=0)])
}

eval_grad_constraints_MPEC <- function(params, wdcMerged, points, length_theta) {
  theta <- params[c(1:length_theta)]
  deltain_stin <- params[-c(1:length_theta)]
  theta1 <- c(theta[1],0,theta[-1])  
  
  #expand deltain to all observations, it is currently #stocked in observations
  deltain <- rep(-30, nrow(wdcMerged))
  deltain[which(wdcMerged$stocked_out==F)] <- deltain_stin
  
  tw_group_list <- unique(wdcMerged$tw_group)
  grad_constraints <- c()
  for(i in 1:length(tw_group_list)) {
    tw_groupin = tw_group_list[i]
    wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
    deltain_tw = deltain[which(wdcMerged$tw_group==tw_groupin)]
    grad_constraints <- c(grad_constraints, 
                                    eval_grad_constraints_MPEC_tw_groupin(deltain_tw, theta1, 
                                          wdcMergedday,points, tw_groupin))
  }
  return(grad_constraints)  
}

eval_grad_constraints_MPEC_tiny <- function(params, wdcMerged, points, length_theta) {
  theta <- params[c(1:length_theta)]
  deltain_stin <- params[-c(1:length_theta)]
  theta1 <- c(theta[1],0,theta[-1])  
  
  #expand deltain to all observations, it is currently #stocked in observations
  deltain <- rep(-30, nrow(wdcMerged))
  deltain[which(wdcMerged$stocked_out==F)] <- deltain_stin
  
  tw_group_list <- unique(wdcMerged$tw_group)
  grad_constraints <- c()
  for(i in 1:length(tw_group_list)) {
    tw_groupin = tw_group_list[i]
    wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
    deltain_tw = deltain[which(wdcMerged$tw_group==tw_groupin)]
    grad_constraints <- c(grad_constraints, 
                          eval_grad_constraints_MPEC_tw_groupin_tiny(deltain_tw, theta1, 
                                                                wdcMergedday,points, tw_groupin))
  }
  return(grad_constraints)  
}

eval_jac_g <- function(params, wdcMerged, points, length_theta) {
  if(print_iter_values) {
    print("In eval_jac_g: ")
  }
  ptm <- proc.time()  
  ret <- c(eval_grad_constraints_MPEC(params, wdcMerged, points, length_theta),
           get_grad_total_density(params, wdcMerged, points))
  if(print_iter_values) {
    print("eval time:")
    print(proc.time()-ptm)
  }
  return(ret)
}

eval_jac_g_tiny <- function(params, wdcMerged, points, length_theta) {
  print("In eval_jac_g: ")
  ptm <- proc.time()  
  ret <- c(eval_grad_constraints_MPEC_tiny(params, wdcMerged, points, length_theta),
           get_grad_total_density(params, wdcMerged, points)/tot_density*target_density_constraint)
  print("eval time:")
  print(proc.time()-ptm)
  return(ret)
}

eval_grad_structure_constraints_MPEC_tw_groupin <- function(deltain, theta1, wdcMergedday, points, tw_groupin, span_start) {
  #create list with nrow(wdcMergedday) elements and
  #each element has vector from (1:lenth(deltain))+span_start.
  #not sure of a simpler way, creating in a loop
#   grad_structure_constraints_element <- list(c(c(1:length(theta1[-2])), span_start+c(1:length(deltain))))  
#   grad_structure_constraints <- rep(grad_structure_constraints_element, nrow(wdcMergedday))
#   return(grad_structure_constraints)
  grad_theta <- eval_grad_lambda_theta_new(theta1, deltain, wdcMergedday, points, tw_groupin)
  grad_theta <- grad_theta[,-2]
  grad_delta <- eval_lambda_delta_list_new(deltain, theta1, wdcMergedday, points, tw_groupin)[[2]]
  
  #reduce the rows and columns corresponding to only stocked in obervations of wdcMergedday
  stocked_list <- which(wdcMergedday$stocked_out==F)
  grad_theta <- grad_theta[stocked_list,]
  grad_delta <- grad_delta[stocked_list, stocked_list]
  
  grad_sparse <- my_make_sparse(grad_theta,grad_delta,span_start)
  return(list(grad_sparse,ncol(grad_delta)))
}

eval_grad_structure_constraints_MPEC <- function(params, wdcMerged, points, length_theta) {
  theta <- params[c(1:length_theta)]
  deltain_stin <- params[-c(1:length_theta)]
  theta1 <- c(theta[1],0,theta[-1])  
  
  #expand deltain to all observations, it is currently #stocked in observations
  deltain <- rep(-30, nrow(wdcMerged))
  deltain[which(wdcMerged$stocked_out==F)] <- deltain_stin
  
  tw_group_list <- unique(wdcMerged$tw_group)
  grad_structure_constraints <- list()
  span_start <- length(theta)
  for(i in 1:length(tw_group_list)) {
    tw_groupin = tw_group_list[i]
    wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
    deltain_tw = deltain[which(wdcMerged$tw_group==tw_groupin)]
    grad_sparse_list <- eval_grad_structure_constraints_MPEC_tw_groupin(deltain_tw, theta1, 
                                                    wdcMergedday,points, tw_groupin, span_start)
    grad_structure_constraints <- c(grad_structure_constraints, grad_sparse_list[[1]])
    span_start <- span_start + grad_sparse_list[[2]]
  }
  return(grad_structure_constraints)  
}  

eval_jac_g_structure <- function(params, wdcMerged, points, length_theta) {
  return(c(eval_grad_structure_constraints_MPEC(params, wdcMerged, points, length_theta),
           list(c(1:length_theta))
    ))
}
  
eval_obj_GMM_MPEC <- function(theta, deltain, wdcMerged, points) {  
  ptm <- proc.time()
  theta1 <- c(theta[1],0,theta[-1])  
  
  print("theta1:")
  print(paste(theta1))
  
  tryCatch({
    ret <- eval_error_xi_sl_MPEC (theta1, deltain, wdcMerged,points)
    gmm_objective <- ret$obj
    gmm_grad_delta <- ret$grad_obj_wrt_delta
    gmm_grad_theta <- rep(0,length(theta))
    gmm_gradient <- c(gmm_grad_theta,gmm_grad_delta)

  },error = function(e) {
    #error-handler-code
    print(paste("MY_ERROR:  ",e))
    print("theta:")
    print(paste(theta))
    ret = 1e12
    gradient = c(0,0)
  }      
  )
  print("obj:")
  print(gmm_objective)
  print("theta:")
  print(theta)
  print("gradient:")
  print(c(gmm_gradient[1:length(theta)],mean(gmm_gradient[-c(1:length(theta))])))
  print("eval_obj time:")
  print(proc.time()-ptm)
  return( list( "objective" = gmm_objective,
                "gradient"  = gmm_gradient) )
  
}


# eval_error_xi_MPEC <- function(theta1,deltain,wdcMerged,points) {
#   stocked_list <- which(wdcMerged$stocked_out==FALSE)
#   #deltain is the length of stocked_list
#   if(length(deltain)!=length(stocked_list)) stop("deltain is not of length stocked_list")
#   deltain_weighted <- deltain*(sqrt(wdcMerged$obs_weight)[stocked_list])
# 
#   if(length(levels(wdcMerged$tract_tw_fac))>1 & length(levels(wdcMerged$week_fac))>1) {
#     X <- model.matrix(~week_fac + tract_tw_fac , data=wdcMerged[stocked_list,])
#   } else if (length(levels(wdcMerged$tract_tw_fac))>1) {
#     X <- model.matrix(~tract_tw_fac , data=wdcMerged[stocked_list,])
#   } else if (length(levels(wdcMerged$week_fac))>1) {
#     X <- model.matrix(~week_fac , data=wdcMerged[stocked_list,])
#   } else {
#     X <- matrix(1,nrow=length(stocked_list), ncol=1)
#   }  
#   if(ncol(X)>0) {
#     X_weighted <- X*(sqrt(wdcMerged$obs_weight)[stocked_list])
#     theta2 = solve(t(X_weighted)%*%X_weighted) %*% (t(X_weighted)%*%deltain_weighted)
#     xi_norm <- deltain - X%*%theta2
#     
#   } else {
#     xi_norm <- deltain
#     theta2 <- c()
#   }
#   
#   #weigh now
#   xi_norm <- xi_norm * sqrt(wdcMerged$obs_weight[stocked_list])
#   xi_norm_full<- rep(0,nrow(wdcMerged))
#   xi_norm_full[stocked_list] <- xi_norm
# #   grad_xi_theta <- grad_xi_theta * sqrt(wdcMerged$obs_weight[stocked_list])
# #   grad_xi_theta_full <- matrix(0,nrow(wdcMerged),ncol(grad_xi_theta))
# #   grad_xi_theta_full[stocked_list,] <- grad_xi_theta
#   
#   #unwieghted xi
#   xi <-   xi_norm/ sqrt(wdcMerged$obs_weight[stocked_list])
#   
#   obj <- c(t(xi_norm ) %*% xi_norm )/length(xi_norm)
#   grad_obj_delta_1 <- xi*(wdcMerged$obs_weight[stocked_list])
#   grad_obj_delta_2 <- t(grad_obj_delta_1) %*% X %*% solve(t(X_weighted)%*%X_weighted) %*% t(X_weighted)
#   grad_obj_delta_2 <- grad_obj_delta_2 * sqrt(wdcMerged$obs_weight[stocked_list])
#   grad_obj_delta <- (c(grad_obj_delta_1)-c(grad_obj_delta_2))*2/length(xi_norm_full)
# #   grad_obj_delta_full <- rep(0,nrow(wdcMerged))
# #   grad_obj_delta_full[stocked_list] <- grad_obj_delta
#   
#   #
#   return(list("obj"=obj,
#               "grad_obj_wrt_delta"=grad_obj_delta))  
#   
# }

eval_error_xi_sl_MPEC <- function(theta1,deltain,wdcMerged,points) {
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  no_st <- max(wdcMerged$station_id_index)
  
#   i=1
#   serv_lvl_covar <- c()
#   instr_serv_lvl_covar <- c()
#   for(tw_groupin in tw_groupin_list) {
#     wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
#     tw_index <- which(tw_groupin_list==tw_groupin)
#     deltain_tw <- deltain[which(wdcMerged$tw_group==tw_groupin)]
#     low = (i-1)*no_st + 1
#     high = i*no_st
#     serv_lvl = user_serv_lvl$serv_lvl[c(low:high)]
#     serv_lvl_covar <- c(serv_lvl_covar, serv_lvl[wdcMergedday$station_id_index])
#     instr_serv_lvl = user_serv_lvl$instr_serv_lvl[c(low:high)]      
#     instr_serv_lvl_covar <- c(instr_serv_lvl_covar, instr_serv_lvl[wdcMergedday$station_id_index])
#     i= i+1
#   }
  if(!identical(order(user_serv_lvl$st_tw_index), c(1:nrow(user_serv_lvl)))) stop("user_serv_lvl not sorted")

  stocked_list <- which(wdcMerged$stocked_out==FALSE)
  
  
  #wdcMerged$serv_lvl_covar <- serv_lvl_covar  
  if(length(levels(wdcMerged$tract_tw_fac))>1 & length(levels(wdcMerged$week_fac))>1) {
    reg_formula <- "~week_fac + 0"
    if(length(levels(wdcMerged$Conditions))>1) paste(reg_formula, "+ Conditions")
    length_levels <- c(length(levels(wdcMerged$Conditions)), length(levels(wdcMerged$Temperature_group)),
                       length(levels(wdcMerged$Humidity_high)), length(levels(wdcMerged$Wind.Speed_high)))
    reg_terms <- c("Conditions", "Temperature_group", "Humidity_high", "Wind.Speed_high")[which(length_levels>1)]
    reg_formula <- paste(reg_formula, paste(reg_terms, collapse = " + "), sep = " + ")
    
    X <- model.matrix(as.formula(reg_formula), data=wdcMerged[stocked_list,])
    #to be able to estimate station fixed effects
    #else below has tract fixed effects which messes up.
    X <- X[,which(!((colnames(X) %like% "tract_tw_fac") & (colnames(X) %like% 
                                                             paste0("_",min(wdcMerged$tw),"$"))))]    
    #since there is no intercept, below in not automatically removed.
    X <- X[,which(!(colnames(X) %like% paste0("week_fac",min(wdcMerged$week))))]    
    if(length(which(colnames(X) %like% "Conditions")) != length(levels(wdcMerged$Conditions))-1) {
      stop("in eval_error_xi_sl_MPEC one conditions number not automatically removed")      
    }
    if(length(which(colnames(X) %like% "Temperature_group")) != length(levels(wdcMerged$Temperature_group))-1) {
      stop("in eval_error_xi_sl_MPEC one Temperature_group number not automatically removed")      
    }
    if(length(which(colnames(X) %like% "Humidity_high")) != length(levels(wdcMerged$Humidity_high))-1) {
      stop("in eval_error_xi_sl_MPEC one Humidity_high number not automatically removed")      
    }
    if(length(which(colnames(X) %like% "Wind.Speed_high")) != length(levels(wdcMerged$Wind.Speed_high))-1) {
      stop("in eval_error_xi_sl_MPEC one Wind.Speed_high number not automatically removed")      
    }
  } else if (length(levels(wdcMerged$tract_tw_fac))>1) {
    X <- model.matrix(~serv_lvl_covar + tract_tw_fac + 0, data=wdcMerged[stocked_list,])
    X <- X[,which(!((colnames(X) %like% "tract_tw_fac") & (colnames(X) %like% 
                                                             paste0("_",min(wdcMerged$tw),"$"))))]        
  } else if (length(levels(wdcMerged$week_fac))>1) {
    X <- model.matrix(~serv_lvl_covar + week_fac + 0, data=wdcMerged[stocked_list,])
    X <- X[,which(!(colnames(X) %like% paste0("week_fac",min(wdcMerged$week))))]    
  } else {
    X <- matrix(1,nrow=length(stocked_list), ncol=1)
  }  
  weight_sum <- ave(sqrt(wdcMerged$obs_weight)[stocked_list],
                    wdcMerged$st_tw_index[stocked_list], FUN=sum)
  
  X_weighted <- X*(sqrt(wdcMerged$obs_weight)[stocked_list])
  X_sum_short <- aggregate(X_weighted, list(Stid=wdcMerged$st_tw_index[stocked_list]), sum)
  X_sum_short <- X_sum_short[,-1]
  X_sum <- X_sum_short[wdcMerged$st_tw_index[stocked_list],]
  X_ave <- X_sum/weight_sum*(sqrt(wdcMerged$obs_weight)[stocked_list])
  X_ddot <- X_weighted - X_ave
  #remove one entry per stationXtw
  sttw_index <- ave(wdcMerged$station_id_index[stocked_list],
                    by=wdcMerged$st_tw_index[stocked_list], FUN=function(x)c(1:length(x)))
  X_ddot <- X_ddot[which(sttw_index!=1),]
  X_ddot <- as.matrix(X_ddot)
  #deltain is the length of stocked_list
  if(length(deltain)!=length(stocked_list)) stop("deltain is not of length stocked_list")
  deltain_weighted <- deltain*(sqrt(wdcMerged$obs_weight)[stocked_list])
  deltain_sum <- ave(deltain_weighted,wdcMerged$st_tw_index[stocked_list], FUN=sum)
  
  deltain_ave <- deltain_sum/weight_sum*(sqrt(wdcMerged$obs_weight)[stocked_list])
  deltain_ddot <- deltain_weighted - deltain_ave
  deltain_ddot <- deltain_ddot[which(sttw_index!=1)]
  
  #   ########
  #   #singlularity test
  #   a <- colSums(X_ddot)
  #   kappa(X_ddot)
  #   min(a)
  #   max(a)
  #   ########
  theta2 = solve(t(X_ddot)%*%X_ddot) %*% (t(X_ddot)%*%deltain_ddot)
  
  weight_sum <- aggregate(sqrt(wdcMerged$obs_weight)[stocked_list], list(St_tw=wdcMerged$st_tw_index[stocked_list]), sum)
  weight_sum <- weight_sum[order(weight_sum$St_tw),2]
  weight_sum_inv <- 1/weight_sum
  
  #X1 <- matrix(rep(1,length(stocked_list)))
  xi_bar <- deltain - X%*%theta2
  xi_bar_weighted <- xi_bar*(sqrt(wdcMerged$obs_weight)[stocked_list])
  gamma_f_w <- aggregate(xi_bar_weighted, list(St_tw=wdcMerged$st_tw_index[stocked_list]), sum)
  gamma_f_w <- gamma_f_w[order(gamma_f_w$St_tw),2]
  gamma_f_w <- gamma_f_w/weight_sum
  
#   stage_3_df <- data.frame(gamma_f_w=gamma_f_w)
#   stage_3_df$tw <- user_serv_lvl$tw
#   stage_3_df$serv_lvl <- user_serv_lvl$serv_lvl
#   stage_3_df$tract_tw <- user_serv_lvl$tract_tw
  stage_3_df <- user_serv_lvl
  #sum of weights at st_tw level and then take sqrt
  weights_st_tw_squared <- aggregate(wdcMerged$obs_weight[stocked_list], list(St_tw=wdcMerged$st_tw_index[stocked_list]), sum)
  weights_st_tw_squared <- weights_st_tw_squared[order(weights_st_tw_squared$St_tw),2]
  weights_st_tw <- sqrt(weights_st_tw_squared)

  if(!length(levels(stage_3_df$tract_tw))>1) stop("in eval_error_xi_sl_MPEC not enough tract_tw levels")
  X2 <- model.matrix(~serv_lvl + tract_tw, data=stage_3_df)
  X2 <- X2[,which(!((colnames(X2) %like% "tract_tw") & (colnames(X2) %like% 
                                                           paste0("_",min(stage_3_df$tw),"$"))))]
  Z2 <- model.matrix(~instr_serv_lvl + serv_lvl_neighbours + tract_tw, data=stage_3_df)
  #Z2 <- model.matrix(~instr_serv_lvl  + tract_tw, data=stage_3_df)
  Z2 <- Z2[,which(!((colnames(Z2) %like% "tract_tw") & (colnames(Z2) %like% 
                                                          paste0("_",min(stage_3_df$tw),"$"))))]
  Z2_weighted <- Z2 * weights_st_tw
  Z2_weighted_sq <- Z2_weighted * weights_st_tw
  X2_weighted <- X2 * weights_st_tw
  X2_full <- X2[wdcMerged$st_tw_index[stocked_list],]
  gamma_f_w_weighted <- gamma_f_w * weights_st_tw  

  #computing matrix P2 and its subparts before computing theta3
  P2_1 <- (t(X2_weighted) %*% Z2_weighted)
  P2_2 <- solve(t(Z2_weighted) %*% Z2_weighted)
  P2 <- solve(P2_1 %*% P2_2 %*% t(P2_1)) %*% P2_1 %*% P2_2 %*% t(Z2_weighted_sq)

  theta3 <- P2%*%gamma_f_w
  predicted_gamma_f_w <-  X2 %*% theta3
  predicted_gamma_f_w_full <- predicted_gamma_f_w[wdcMerged$st_tw_index[stocked_list]]

  eta = xi_bar - predicted_gamma_f_w_full
  eta_weighted <- eta*(sqrt(wdcMerged$obs_weight)[stocked_list])
  obj <- c(t(eta_weighted ) %*% eta_weighted )/length(eta_weighted)
#   xi_norm <- xi_bar - X1%*%(t(wdcMerged$obs_weight[stocked_list])%*%xi_bar)/(sum(wdcMerged$obs_weight[stocked_list]))  
#   #xi_norm is unweighted value of \eta corresponding to the document.
#   xi_weighted <- xi_norm*(sqrt(wdcMerged$obs_weight)[stocked_list])
#   
#   obj <- c(t(xi_weighted ) %*% xi_weighted )/length(xi_weighted)
  
  #gradient computation #refer to xi_gradient_formulation.lyx for analytical for below.
  
  wdcMerged_Sttw_mat <- model.matrix(~ factor(st_tw_index) + 0, data=wdcMerged[stocked_list,])
  row.names(wdcMerged_Sttw_mat) <- NULL 
  colnames(wdcMerged_Sttw_mat) <- NULL
#   W_a <- t(wdcMerged_Sttw_mat * sqrt(wdcMerged$obs_weight[stocked_list]))
#   W_a <- W_a/rowSums(W_a)

#   tx <- diag(nrow=max(wdcMerged$st_tw_index[stocked_list]))
#   tx <- tx*weight_sum_inv
#   tx <- tx[,wdcMerged$st_tw_index[stocked_list]]
#   tx <- t(t(tx) *  sqrt(wdcMerged$obs_weight[stocked_list]))

  grad_delta_hat_wrt_delta <- t(t(wdcMerged_Sttw_mat) * weight_sum_inv)
  grad_obj_wrt_delta_1 <- eta*(wdcMerged$obs_weight[stocked_list])
  grad_obj_wrt_delta_2 <- t(grad_obj_wrt_delta_1) %*% X %*% 
    solve(t(X_ddot)%*%X_ddot) %*% t(X_ddot)
  grad_obj_wrt_delta_2 <- grad_obj_wrt_delta_2 * (sqrt(wdcMerged$obs_weight[stocked_list])[which(sttw_index!=1)])
  grad_obj_wrt_delta_2_full <- rep(0,length(stocked_list))
  grad_obj_wrt_delta_2_full[which(sttw_index!=1)] <- grad_obj_wrt_delta_2
  rm(grad_obj_wrt_delta_2)
  grad_obj_wrt_delta_3 <- grad_obj_wrt_delta_2_full %*%  grad_delta_hat_wrt_delta
  grad_obj_wrt_delta_3_full <- grad_obj_wrt_delta_3[wdcMerged$st_tw_index[stocked_list]]
  grad_obj_wrt_delta_3_full <- grad_obj_wrt_delta_3_full* sqrt(wdcMerged$obs_weight[stocked_list])
  #rm(grad_obj_wrt_delta_3)
  grad_obj_wrt_delta <- (c(grad_obj_wrt_delta_1)-c(grad_obj_wrt_delta_2_full)+c(grad_obj_wrt_delta_3_full)
  )*2/length(eta)
  
# grad_obj_wrt_delta_4 <- t(eta*(wdcMerged$obs_weight[stocked_list])) %*% X2_full %*%
#   solve(t(Z2_weighted_sq) %*% X2) %*% t(Z2_weighted_sq)
# grad_obj_wrt_delta_4 <- grad_obj_wrt_delta_4 %*% W_a

Wa <- diag(nrow=max(wdcMerged$st_tw_index[stocked_list]))
Wa <- Wa*weight_sum_inv
grad_obj_wrt_delta_4 <- t(eta*(wdcMerged$obs_weight[stocked_list])) %*% X2_full %*% P2
grad_obj_wrt_delta_4 <- grad_obj_wrt_delta_4 %*% Wa
grad_obj_wrt_delta_4 <- grad_obj_wrt_delta_4[,wdcMerged$st_tw_index[stocked_list]]
grad_obj_wrt_delta_4 <- grad_obj_wrt_delta_4 *  sqrt(wdcMerged$obs_weight[stocked_list])

grad_obj_wrt_delta_5 <- grad_obj_wrt_delta_4 %*% X %*% 
  solve(t(X_ddot)%*%X_ddot) %*% t(X_ddot)
grad_obj_wrt_delta_5 <- grad_obj_wrt_delta_5 * (sqrt(wdcMerged$obs_weight[stocked_list])[which(sttw_index!=1)])
grad_obj_wrt_delta_5_full <- rep(0,length(stocked_list))
grad_obj_wrt_delta_5_full[which(sttw_index!=1)] <- grad_obj_wrt_delta_5

grad_obj_wrt_delta_6 <- grad_obj_wrt_delta_5_full %*%  grad_delta_hat_wrt_delta
grad_obj_wrt_delta_6_full <- grad_obj_wrt_delta_6[wdcMerged$st_tw_index[stocked_list]]
grad_obj_wrt_delta_6_full <- grad_obj_wrt_delta_6_full* sqrt(wdcMerged$obs_weight[stocked_list])  

grad_obj_wrt_delta <- (c(grad_obj_wrt_delta_1)-c(grad_obj_wrt_delta_2_full)+c(grad_obj_wrt_delta_3_full)
                       -c(grad_obj_wrt_delta_4)+c(grad_obj_wrt_delta_5_full)-c(grad_obj_wrt_delta_6_full))*2/length(eta)


  return(list("obj"=obj,
              "theta2"=theta2,
              "theta3"=theta3,
              "grad_obj_wrt_delta"=grad_obj_wrt_delta))  
  
}

eval_obj_GMM_MPEC_obj <- function(params, wdcMerged, points, length_theta) {    
  #obj_eval <<- 1
  theta <- params[c(1:length_theta)]
  deltain_stin <- params[-c(1:length_theta)]
  print("objective comput")
  return(eval_obj_GMM_MPEC(theta, deltain_stin, wdcMerged, points)$objective)
}  
eval_obj_GMM_MPEC_grad <- function(params, wdcMerged, points, length_theta) {    
  #obj_eval <<- 0
  theta <- params[c(1:length_theta)]
  deltain_stin <- params[-c(1:length_theta)]  
  print("gradint comput")
  return(eval_obj_GMM_MPEC(theta, deltain_stin, wdcMerged, points)$gradient)
} 

