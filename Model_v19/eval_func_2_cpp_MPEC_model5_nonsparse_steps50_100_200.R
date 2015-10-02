weighing_GMM_mat <<- NULL

eval_constraints_MPEC_tw_groupin <- function(deltain, theta1, wdcMergedday, points, tw_groupin) {
  stocked_list <- which(wdcMergedday$stocked_out==FALSE)
  dem_T <- eval_lambda_new(deltain, theta1, wdcMergedday, points, tw_groupin)[stocked_list]
  dem_hat_T <- wdcMergedday$out_dem_sum[stocked_list]/wdcMergedday$obs_weight[stocked_list]
  obj <- dem_T-dem_hat_T
  return(obj)
}


eval_constraints_MPEC <- function(deltain_stin, theta1, wdcMerged, points) {
  
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


eval_constraints_moment_conditions <- function(deltain, theta1, eta, wdcMerged, points) {  
  if(!identical(order(user_serv_lvl$st_tw_index), c(1:nrow(user_serv_lvl)))) stop("user_serv_lvl not sorted")  
  stocked_list <- which(wdcMerged$stocked_out==FALSE)
  #get covariates
  if(is.null(X) | is.null(Z)) {
    list_covariates <- eval_covariates_delta_reg(deltain,theta1,wdcMerged,points)
    X <<- list_covariates$X
    Z <<- list_covariates$Z
  }
  
  #theta_2 is given by, Inv(X'.W'.Z.A_N.Z'.W.X).(X'.W'.Z.A_N.Z'.W.delta)
  #W above can be replaced by W/|W|
  
  #Z_weighted equal to W'.Z
  Z_weighted <- (Z * (wdcMerged$obs_weight[stocked_list]))/sum(wdcMerged$obs_weight[stocked_list])
  if(is.null(weighing_GMM_mat)) {
    A_N <- diag(ncol(Z))    
  } else {
    A_N <- weighing_GMM_mat
  }
  if(print_iter_values) {
    print(paste0("kappa at theta_2 computation: ", kappa(t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted) %*% X)))  
  }  
  theta_2 <- solve(t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted) %*% X) %*%
    (t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted) %*% deltain)
  
  xi <- deltain - (X %*%  theta_2) 
  
  obj_moments <- t(Z_weighted) %*% xi
  
  return(obj_moments - eta)  
}

eval_grad_constraints_moment_conditions <- function(deltain, theta1, eta, wdcMerged, points) {
  if(!identical(order(user_serv_lvl$st_tw_index), c(1:nrow(user_serv_lvl)))) stop("user_serv_lvl not sorted")  
  stocked_list <- which(wdcMerged$stocked_out==FALSE)
  #get covariates
  if(is.null(X) | is.null(Z)) {
    list_covariates <- eval_covariates_delta_reg(deltain,theta1,wdcMerged,points)
    X <<- list_covariates$X
    Z <<- list_covariates$Z
  }
  #theta_2 is given by, Inv(X'.W'.Z.A_N.Z'.W.X).(X'.W'.Z.A_N.Z'.W.delta)
  #W above can be replaced by W/|W|
  
  #Z_weighted equal to W'.Z
  Z_weighted <- (Z * (wdcMerged$obs_weight[stocked_list]))/sum(wdcMerged$obs_weight[stocked_list])
  if(is.null(weighing_GMM_mat)) {
    A_N <- diag(ncol(Z))    
  } else {
    A_N <- weighing_GMM_mat
  }
  
  theta_2 <- solve(t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted) %*% X) %*%
    (t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted) %*% deltain)
  
  xi <- deltain - (X %*%  theta_2) 

  obj_moments <- t(Z_weighted) %*% xi
  ret <- obj_moments - eta
#   grad_theta <- matrix(0,length(obj_moments), length(theta1))
#   grad_theta <- grad_theta[,-2]
#   grad_xi_delta <- diag(1,nrow(X)) - X %*% solve(t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted) %*% X) %*%
#     (t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted))
#   grad_delta <- t(Z_weighted) %*% grad_xi_delta
  grad_delta <- t(Z_weighted) - t(Z_weighted) %*% X %*% solve(t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted) %*% X) %*%
    (t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted))

  #grad_delta <- grad_delta[, stocked_list]   #already done above.
  grad_eta <- diag(-1,length(eta))
  grad_all <- c(t(cbind(grad_delta,grad_eta)))
  return(grad_all[which(grad_all!=0)])
}
  



eval_grad_constraints_MPEC_tw_groupin <- function(deltain, theta1, wdcMergedday, points, tw_groupin) {
  grad_all <- eval_grad_lambda_new(theta1, deltain, wdcMergedday, points, tw_groupin)
  grad_all <- grad_all[,-2] #remove rand coef col
  
  #reduce the rows and columns corresponding to only stocked in obervations of wdcMergedday
  stocked_list <- which(wdcMergedday$stocked_out==F)
  length_theta <- length(theta1)-1
  keep_cols <- c(c(1:length_theta),length_theta+stocked_list)
  grad_all <- grad_all[stocked_list,keep_cols]
  return(linearize_sparsify(grad_all))
}




eval_grad_constraints_MPEC <- function(deltain_stin, theta1, eta, wdcMerged, points) {  
  #expand deltain to all observations, it is currently #stocked in observations
  deltain <- rep(-30, nrow(wdcMerged))
  deltain[which(wdcMerged$stocked_out==F)] <- deltain_stin
  # theta <- theta1[-2]
  # stocked_list <- which(wdcMerged$stocked_out==F)
  
  tw_group_list <- unique(wdcMerged$tw_group)
  grad_constraints <- c()
  for(i in 1:length(tw_group_list)) {
    tw_groupin = tw_group_list[i]
    wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
    deltain_tw = deltain[which(wdcMerged$tw_group==tw_groupin)]
    grad_constraints <- c(grad_constraints, 
                                    eval_grad_constraints_MPEC_tw_groupin(deltain_tw, theta1, 
                                          wdcMergedday,points, tw_groupin))
#     rows <- which(wdcMerged$tw_group[stocked_list]==tw_groupin)
#     cols <- c(c(1:length(theta)),length(theta)+rows)
#     grad_constraints[rows,cols] <- eval_grad_constraints_MPEC_tw_groupin(deltain_tw, theta1, 
#                                           wdcMergedday,points, tw_groupin)

  }
  return(grad_constraints)  
}



eval_g <- function(params, wdcMerged, points, length_theta, length_delta, length_eta) {
#  print("In eval_g: ")  
  ptm <- proc.time()  
  theta <- params[c(1:length_theta)]
  deltain_stin <- params[c((length_theta+1):(length_theta+length_delta))]
  eta <- params[c((length_theta+length_delta+1):(length_theta+length_delta+length_eta))]
  theta1 <- c(theta[1],0,theta[-1])  
  
  ret <- c(eval_constraints_MPEC(deltain_stin, theta1, wdcMerged, points),
           eval_constraints_moment_conditions(deltain_stin, theta1, eta, wdcMerged, points),
           get_total_density(params, wdcMerged, points))
  if(length(ret)!=(length_delta+length_eta+1)) stop ("eval_g error")
  if(print_iter_values) {
    print("In eval_g: summary of constraints. delta constraints, eta constraint, density constraint")
    print(as.numeric(summary(ret[1:length_delta],digits = 10)))
    print(as.numeric(summary(ret[length_delta+c(1:length_eta)],digits = 10)))
    print(as.numeric(ret[length_delta+1+length_eta]))
    print('')
    params_save <<- params
  }
#   print("eval time:")
#   print(proc.time()-ptm)
  return(ret)  
}

eval_jac_g <- function(params, wdcMerged, points, length_theta, length_delta, length_eta) {
#   print("In eval_jac_g: ")
  params <- params + rep(1e-32, length(params)) # to prevent from values in sparse jacobina left out due to 0 param values
  ptm <- proc.time()  
  theta <- params[c(1:length_theta)]
  deltain_stin <- params[c((length_theta+1):(length_theta+length_delta))]
  eta <- params[c((length_theta+length_delta+1):(length_theta+length_delta+length_eta))]
  theta1 <- c(theta[1],0,theta[-1])  
  
  ret <- c(eval_grad_constraints_MPEC(deltain_stin, theta1, eta, wdcMerged, points),
           eval_grad_constraints_moment_conditions(deltain_stin, theta1, eta, wdcMerged, points),
           get_grad_total_density(params, wdcMerged, points))

#   print("eval time:")
#   print(proc.time()-ptm)
  if(length(ret)!=length(unlist(eval_jac_g_structure_val))) stop("eval_jac_g length inconsistent with eval_jac_g_structure_val")
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

eval_grad_structure_constraints_MPEC <- function(deltain_stin, theta1, eta, wdcMerged, points) {
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

eval_grad_structure_moment_conditions <- function(deltain, theta1, eta, wdcMerged, points) {
  if(!identical(order(user_serv_lvl$st_tw_index), c(1:nrow(user_serv_lvl)))) stop("user_serv_lvl not sorted")  
  stocked_list <- which(wdcMerged$stocked_out==FALSE)
  #get covariates
  if(is.null(X) | is.null(Z)) {
    list_covariates <- eval_covariates_delta_reg(deltain,theta1,wdcMerged,points)
    X <<- list_covariates$X
    Z <<- list_covariates$Z
  }
  
  
  #theta_2 is given by, Inv(X'.W'.Z.A_N.Z'.W.X).(X'.W'.Z.A_N.Z'.W.delta)
  #W above can be replaced by W/|W|
  
  #Z_weighted equal to W'.Z
  Z_weighted <- (Z * (wdcMerged$obs_weight[stocked_list]))/sum(wdcMerged$obs_weight[stocked_list])
  if(is.null(weighing_GMM_mat)) {
    A_N <- diag(ncol(Z))    
  } else {
    A_N <- weighing_GMM_mat
  }
  
  theta_2 <- solve(t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted) %*% X) %*%
    (t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted) %*% deltain)
  
  xi <- deltain - (X %*%  theta_2) 
#   grad_xi_delta <- diag(1,nrow(X)) - X %*% solve(t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted) %*% X) %*%
#     (t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted))
  #seems like grad_xi_delta is not sparse.
 
  obj_moments <- t(Z_weighted) %*% xi
  ret <- obj_moments - eta
  grad_theta <- matrix(0,length(obj_moments), length(theta1))
  grad_theta <- grad_theta[,-2]
  #grad_delta <- t(Z_weighted) %*% grad_xi_delta
  grad_delta <- t(Z_weighted) - t(Z_weighted) %*% X %*% solve(t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted) %*% X) %*%
    (t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted))
  grad_eta <- diag(-1,length(eta))
  grad_all <- cbind(grad_theta, grad_delta,grad_eta)

  grad_sparse <- make.sparse(grad_all)
  return(grad_sparse)
}

eval_jac_g_structure <- function(params, wdcMerged, points, length_theta, length_delta, length_eta) {
  params <- params + rep(1e-32, length(params)) # to prevent from values in sparse jacobina left out due to 0 param values
  theta <- params[c(1:length_theta)]
  deltain_stin <- params[c((length_theta+1):(length_theta+length_delta))]
  eta <- params[c((length_theta+length_delta+1):(length_theta+length_delta+length_eta))]
  theta1 <- c(theta[1],0,theta[-1])  
  
  return(c(eval_grad_structure_constraints_MPEC(deltain_stin, theta1, eta, wdcMerged, points),
           eval_grad_structure_moment_conditions(deltain_stin, theta1, eta, wdcMerged, points),
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




# eval_error_xi_sl_MPEC <- function(theta1,deltain,wdcMerged,points) {
#   tw_groupin_list <- unique(wdcMerged$tw_group)    
#   no_st <- max(wdcMerged$station_id_index)
#   
# #   i=1
# #   serv_lvl_covar <- c()
# #   instr_serv_lvl_covar <- c()
# #   for(tw_groupin in tw_groupin_list) {
# #     wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
# #     tw_index <- which(tw_groupin_list==tw_groupin)
# #     deltain_tw <- deltain[which(wdcMerged$tw_group==tw_groupin)]
# #     low = (i-1)*no_st + 1
# #     high = i*no_st
# #     serv_lvl = user_serv_lvl$serv_lvl[c(low:high)]
# #     serv_lvl_covar <- c(serv_lvl_covar, serv_lvl[wdcMergedday$station_id_index])
# #     instr_serv_lvl = user_serv_lvl$instr_serv_lvl[c(low:high)]      
# #     instr_serv_lvl_covar <- c(instr_serv_lvl_covar, instr_serv_lvl[wdcMergedday$station_id_index])
# #     i= i+1
# #   }
#   if(!identical(order(user_serv_lvl$st_tw_index), c(1:nrow(user_serv_lvl)))) stop("user_serv_lvl not sorted")
# 
#   stocked_list <- which(wdcMerged$stocked_out==FALSE)
#   
#   
#   #wdcMerged$serv_lvl_covar <- serv_lvl_covar  
#   if(length(levels(wdcMerged$tract_tw_fac))>1 & length(levels(wdcMerged$week_fac))>1) {
#     reg_formula <- "~week_fac + 0"
#     if(length(levels(wdcMerged$Conditions))>1) paste(reg_formula, "+ Conditions")
#     length_levels <- c(length(levels(wdcMerged$Conditions)), length(levels(wdcMerged$Temperature_group)),
#                        length(levels(wdcMerged$Humidity_high)), length(levels(wdcMerged$Wind.Speed_high)))
#     reg_terms <- c("Conditions", "Temperature_group", "Humidity_high", "Wind.Speed_high")[which(length_levels>1)]
#     reg_formula <- paste(reg_formula, paste(reg_terms, collapse = " + "), sep = " + ")
#     
#     X <- model.matrix(as.formula(reg_formula), data=wdcMerged[stocked_list,])
#     #to be able to estimate station fixed effects
#     #else below has tract fixed effects which messes up.
#     X <- X[,which(!((colnames(X) %like% "tract_tw_fac") & (colnames(X) %like% 
#                                                              paste0("_",min(wdcMerged$tw),"$"))))]    
#     #since there is no intercept, below in not automatically removed.
#     X <- X[,which(!(colnames(X) %like% paste0("week_fac",min(wdcMerged$week))))]    
#     if(length(which(colnames(X) %like% "Conditions")) != length(levels(wdcMerged$Conditions))-1) {
#       stop("in eval_error_xi_sl_MPEC one conditions number not automatically removed")      
#     }
#     if(length(which(colnames(X) %like% "Temperature_group")) != length(levels(wdcMerged$Temperature_group))-1) {
#       stop("in eval_error_xi_sl_MPEC one Temperature_group number not automatically removed")      
#     }
#     if(length(which(colnames(X) %like% "Humidity_high")) != length(levels(wdcMerged$Humidity_high))-1) {
#       stop("in eval_error_xi_sl_MPEC one Humidity_high number not automatically removed")      
#     }
#     if(length(which(colnames(X) %like% "Wind.Speed_high")) != length(levels(wdcMerged$Wind.Speed_high))-1) {
#       stop("in eval_error_xi_sl_MPEC one Wind.Speed_high number not automatically removed")      
#     }
#   } else if (length(levels(wdcMerged$tract_tw_fac))>1) {
#     X <- model.matrix(~serv_lvl_covar + tract_tw_fac + 0, data=wdcMerged[stocked_list,])
#     X <- X[,which(!((colnames(X) %like% "tract_tw_fac") & (colnames(X) %like% 
#                                                              paste0("_",min(wdcMerged$tw),"$"))))]        
#   } else if (length(levels(wdcMerged$week_fac))>1) {
#     X <- model.matrix(~serv_lvl_covar + week_fac + 0, data=wdcMerged[stocked_list,])
#     X <- X[,which(!(colnames(X) %like% paste0("week_fac",min(wdcMerged$week))))]    
#   } else {
#     X <- matrix(1,nrow=length(stocked_list), ncol=1)
#   }  
#   weight_sum <- ave(sqrt(wdcMerged$obs_weight)[stocked_list],
#                     wdcMerged$st_tw_index[stocked_list], FUN=sum)
#   
#   X_weighted <- X*(sqrt(wdcMerged$obs_weight)[stocked_list])
#   X_sum_short <- aggregate(X_weighted, list(Stid=wdcMerged$st_tw_index[stocked_list]), sum)
#   X_sum_short <- X_sum_short[,-1]
#   X_sum <- X_sum_short[wdcMerged$st_tw_index[stocked_list],]
#   X_ave <- X_sum/weight_sum*(sqrt(wdcMerged$obs_weight)[stocked_list])
#   X_ddot <- X_weighted - X_ave
#   #remove one entry per stationXtw
#   sttw_index <- ave(wdcMerged$station_id_index[stocked_list],
#                     by=wdcMerged$st_tw_index[stocked_list], FUN=function(x)c(1:length(x)))
#   X_ddot <- X_ddot[which(sttw_index!=1),]
#   X_ddot <- as.matrix(X_ddot)
#   #deltain is the length of stocked_list
#   if(length(deltain)!=length(stocked_list)) stop("deltain is not of length stocked_list")
#   deltain_weighted <- deltain*(sqrt(wdcMerged$obs_weight)[stocked_list])
#   deltain_sum <- ave(deltain_weighted,wdcMerged$st_tw_index[stocked_list], FUN=sum)
#   
#   deltain_ave <- deltain_sum/weight_sum*(sqrt(wdcMerged$obs_weight)[stocked_list])
#   deltain_ddot <- deltain_weighted - deltain_ave
#   deltain_ddot <- deltain_ddot[which(sttw_index!=1)]
#   
#   #   ########
#   #   #singlularity test
#   #   a <- colSums(X_ddot)
#   #   kappa(X_ddot)
#   #   min(a)
#   #   max(a)
#   #   ########
#   theta2 = solve(t(X_ddot)%*%X_ddot) %*% (t(X_ddot)%*%deltain_ddot)
#   
#   weight_sum <- aggregate(sqrt(wdcMerged$obs_weight)[stocked_list], list(St_tw=wdcMerged$st_tw_index[stocked_list]), sum)
#   weight_sum <- weight_sum[order(weight_sum$St_tw),2]
#   weight_sum_inv <- 1/weight_sum
#   
#   #X1 <- matrix(rep(1,length(stocked_list)))
#   xi_bar <- deltain - X%*%theta2
#   xi_bar_weighted <- xi_bar*(sqrt(wdcMerged$obs_weight)[stocked_list])
#   gamma_f_w <- aggregate(xi_bar_weighted, list(St_tw=wdcMerged$st_tw_index[stocked_list]), sum)
#   gamma_f_w <- gamma_f_w[order(gamma_f_w$St_tw),2]
#   gamma_f_w <- gamma_f_w/weight_sum
#   
# #   stage_3_df <- data.frame(gamma_f_w=gamma_f_w)
# #   stage_3_df$tw <- user_serv_lvl$tw
# #   stage_3_df$serv_lvl <- user_serv_lvl$serv_lvl
# #   stage_3_df$tract_tw <- user_serv_lvl$tract_tw
#   stage_3_df <- user_serv_lvl
#   #sum of weights at st_tw level and then take sqrt
#   weights_st_tw_squared <- aggregate(wdcMerged$obs_weight[stocked_list], list(St_tw=wdcMerged$st_tw_index[stocked_list]), sum)
#   weights_st_tw_squared <- weights_st_tw_squared[order(weights_st_tw_squared$St_tw),2]
#   weights_st_tw <- sqrt(weights_st_tw_squared)
# 
#   if(!length(levels(stage_3_df$tract_tw))>1) stop("in eval_error_xi_sl_MPEC not enough tract_tw levels")
#   X2 <- model.matrix(~serv_lvl + tract_tw, data=stage_3_df)
#   X2 <- X2[,which(!((colnames(X2) %like% "tract_tw") & (colnames(X2) %like% 
#                                                            paste0("_",min(stage_3_df$tw),"$"))))]
#   Z2 <- model.matrix(~instr_serv_lvl + serv_lvl_neighbours + tract_tw, data=stage_3_df)
#   #Z2 <- model.matrix(~instr_serv_lvl  + tract_tw, data=stage_3_df)
#   Z2 <- Z2[,which(!((colnames(Z2) %like% "tract_tw") & (colnames(Z2) %like% 
#                                                           paste0("_",min(stage_3_df$tw),"$"))))]
#   Z2_weighted <- Z2 * weights_st_tw
#   Z2_weighted_sq <- Z2_weighted * weights_st_tw
#   X2_weighted <- X2 * weights_st_tw
#   X2_full <- X2[wdcMerged$st_tw_index[stocked_list],]
#   gamma_f_w_weighted <- gamma_f_w * weights_st_tw  
# 
#   #computing matrix P2 and its subparts before computing theta3
#   P2_1 <- (t(X2_weighted) %*% Z2_weighted)
#   P2_2 <- solve(t(Z2_weighted) %*% Z2_weighted)
#   P2 <- solve(P2_1 %*% P2_2 %*% t(P2_1)) %*% P2_1 %*% P2_2 %*% t(Z2_weighted_sq)
# 
#   theta3 <- P2%*%gamma_f_w
#   predicted_gamma_f_w <-  X2 %*% theta3
#   predicted_gamma_f_w_full <- predicted_gamma_f_w[wdcMerged$st_tw_index[stocked_list]]
# 
#   eta = xi_bar - predicted_gamma_f_w_full
#   eta_weighted <- eta*(sqrt(wdcMerged$obs_weight)[stocked_list])
#   obj <- c(t(eta_weighted ) %*% eta_weighted )/length(eta_weighted)
# #   xi_norm <- xi_bar - X1%*%(t(wdcMerged$obs_weight[stocked_list])%*%xi_bar)/(sum(wdcMerged$obs_weight[stocked_list]))  
# #   #xi_norm is unweighted value of \eta corresponding to the document.
# #   xi_weighted <- xi_norm*(sqrt(wdcMerged$obs_weight)[stocked_list])
# #   
# #   obj <- c(t(xi_weighted ) %*% xi_weighted )/length(xi_weighted)
#   
#   #gradient computation #refer to xi_gradient_formulation.lyx for analytical for below.
#   
#   wdcMerged_Sttw_mat <- model.matrix(~ factor(st_tw_index) + 0, data=wdcMerged[stocked_list,])
#   row.names(wdcMerged_Sttw_mat) <- NULL 
#   colnames(wdcMerged_Sttw_mat) <- NULL
# #   W_a <- t(wdcMerged_Sttw_mat * sqrt(wdcMerged$obs_weight[stocked_list]))
# #   W_a <- W_a/rowSums(W_a)
# 
# #   tx <- diag(nrow=max(wdcMerged$st_tw_index[stocked_list]))
# #   tx <- tx*weight_sum_inv
# #   tx <- tx[,wdcMerged$st_tw_index[stocked_list]]
# #   tx <- t(t(tx) *  sqrt(wdcMerged$obs_weight[stocked_list]))
# 
#   grad_delta_hat_wrt_delta <- t(t(wdcMerged_Sttw_mat) * weight_sum_inv)
#   grad_obj_wrt_delta_1 <- eta*(wdcMerged$obs_weight[stocked_list])
#   grad_obj_wrt_delta_2 <- t(grad_obj_wrt_delta_1) %*% X %*% 
#     solve(t(X_ddot)%*%X_ddot) %*% t(X_ddot)
#   grad_obj_wrt_delta_2 <- grad_obj_wrt_delta_2 * (sqrt(wdcMerged$obs_weight[stocked_list])[which(sttw_index!=1)])
#   grad_obj_wrt_delta_2_full <- rep(0,length(stocked_list))
#   grad_obj_wrt_delta_2_full[which(sttw_index!=1)] <- grad_obj_wrt_delta_2
#   rm(grad_obj_wrt_delta_2)
#   grad_obj_wrt_delta_3 <- grad_obj_wrt_delta_2_full %*%  grad_delta_hat_wrt_delta
#   grad_obj_wrt_delta_3_full <- grad_obj_wrt_delta_3[wdcMerged$st_tw_index[stocked_list]]
#   grad_obj_wrt_delta_3_full <- grad_obj_wrt_delta_3_full* sqrt(wdcMerged$obs_weight[stocked_list])
#   #rm(grad_obj_wrt_delta_3)
#   grad_obj_wrt_delta <- (c(grad_obj_wrt_delta_1)-c(grad_obj_wrt_delta_2_full)+c(grad_obj_wrt_delta_3_full)
#   )*2/length(eta)
#   
# # grad_obj_wrt_delta_4 <- t(eta*(wdcMerged$obs_weight[stocked_list])) %*% X2_full %*%
# #   solve(t(Z2_weighted_sq) %*% X2) %*% t(Z2_weighted_sq)
# # grad_obj_wrt_delta_4 <- grad_obj_wrt_delta_4 %*% W_a
# 
# Wa <- diag(nrow=max(wdcMerged$st_tw_index[stocked_list]))
# Wa <- Wa*weight_sum_inv
# grad_obj_wrt_delta_4 <- t(eta*(wdcMerged$obs_weight[stocked_list])) %*% X2_full %*% P2
# grad_obj_wrt_delta_4 <- grad_obj_wrt_delta_4 %*% Wa
# grad_obj_wrt_delta_4 <- grad_obj_wrt_delta_4[,wdcMerged$st_tw_index[stocked_list]]
# grad_obj_wrt_delta_4 <- grad_obj_wrt_delta_4 *  sqrt(wdcMerged$obs_weight[stocked_list])
# 
# grad_obj_wrt_delta_5 <- grad_obj_wrt_delta_4 %*% X %*% 
#   solve(t(X_ddot)%*%X_ddot) %*% t(X_ddot)
# grad_obj_wrt_delta_5 <- grad_obj_wrt_delta_5 * (sqrt(wdcMerged$obs_weight[stocked_list])[which(sttw_index!=1)])
# grad_obj_wrt_delta_5_full <- rep(0,length(stocked_list))
# grad_obj_wrt_delta_5_full[which(sttw_index!=1)] <- grad_obj_wrt_delta_5
# 
# grad_obj_wrt_delta_6 <- grad_obj_wrt_delta_5_full %*%  grad_delta_hat_wrt_delta
# grad_obj_wrt_delta_6_full <- grad_obj_wrt_delta_6[wdcMerged$st_tw_index[stocked_list]]
# grad_obj_wrt_delta_6_full <- grad_obj_wrt_delta_6_full* sqrt(wdcMerged$obs_weight[stocked_list])  
# 
# grad_obj_wrt_delta <- (c(grad_obj_wrt_delta_1)-c(grad_obj_wrt_delta_2_full)+c(grad_obj_wrt_delta_3_full)
#                        -c(grad_obj_wrt_delta_4)+c(grad_obj_wrt_delta_5_full)-c(grad_obj_wrt_delta_6_full))*2/length(eta)
# 
# 
#   return(list("obj"=obj,
#               "theta2"=theta2,
#               "theta3"=theta3,
#               "grad_obj_wrt_delta"=grad_obj_wrt_delta))  
#   
# }

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


eval_hessian <- function(params, wdcMerged, points, length_theta, length_delta, length_eta, obj_factor, hessian_lambda) {
#   print("In eval_hessian: ")
  params <- params + rep(1e-32, length(params)) # to prevent from values in sparse jacobina left out due to 0 param values
  ptm <- proc.time()  
  theta <- params[c(1:length_theta)]
  deltain_stin <- params[c((length_theta+1):(length_theta+length_delta))]
  eta <- params[c((length_theta+length_delta+1):(length_theta+length_delta+length_eta))]
  theta1 <- c(theta[1],0,theta[-1])  
  #making multipliers non-zero so that the return values confirm to the sparse strcuture.
  obj_factor <- obj_factor+1e-32
  constraint_multipliers <- hessian_lambda
  constraint_multipliers <- constraint_multipliers+1e-32
  lambda_multiplers_in <- constraint_multipliers[c(1:length_delta)]
  
  ret <- c(eval_hessian_lambda_constraints_MPEC(deltain_stin, theta1, 
                                                eta, wdcMerged, points, lambda_multiplers_in),
           eval_hessian_eta_sq(obj_factor, length_eta)
  )
  if(length(ret)!=length(unlist(eval_hess_struct))) stop("eval_hessian length inconsistent with eval_hess_struct")
  #make 0 the values less than 1e-16, they are only non-zero due to above shifting of multiplers and 
  #will otherwise cause issues with condition number of matrix.
  #ret[which(abs(ret)<=1e-32)] <- 0
#   print("eval time:")
#   print(proc.time()-ptm)
# params_save <<- params
# print(obj_factor)
# print(constraint_multipliers[1:10])
# print(summary(constraint_multipliers))
# print(summary(params))
# print(theta)
# print(ret[1:10])
  return(ret)
}

eval_hessian_eta_sq <- function(obj_factor, length_eta) {  
  if(is.null(weighing_GMM_mat)) {
    A_N <- diag(length_eta)    
  } else {
    A_N <- weighing_GMM_mat
  }
  A_N <- keep_lower_traingular_matrix(A_N)
  return(obj_factor*c(t(2*A_N)))  
}





eval_hessian_structure <- function(params, wdcMerged, points, length_theta, length_delta, length_eta) {
  params <- params + rep(1e-32, length(params)) # to prevent from values in sparse jacobina left out due to 0 param values
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
                       get_points_density_grad_addden_cols(points, tw_in)
  )
  points_mat <- points
  points_mat$density <- get_points_density(points_mat, theta1, tw_in)  
  points_mat <- cbind(points_mat[,c("lat","lon","density")], density_mat)
  points_mat = as.matrix(points_mat)
  
  wdcMergedday = as.matrix(wdcMergedday)
  
  hessian_lambda_theta1_delta <- eval_hessian_lambda_theta1_delta_cpp(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
                                                                      as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations),
                                                                      lambda_multiplers_tw, nonden_ceoflength)
  hessian_lambda_delta_sq <- eval_hessian_lambda_delta_sq_cpp(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
                                                              as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations),
                                                              lambda_multiplers_tw, nonden_ceoflength)
  hessian_lambda_theta1_sq <- eval_hessian_lambda_theta1_sq_cpp(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
                                                                as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations),
                                                                lambda_multiplers_tw, nonden_ceoflength)
  
  #remove 2 index of theta1
  hessian_lambda_theta1_sq <- hessian_lambda_theta1_sq[-2,-2]
  hessian_lambda_theta1_sq <- keep_lower_traingular_matrix(hessian_lambda_theta1_sq)
  hessian_lambda_theta1_delta <- hessian_lambda_theta1_delta[-2,]
  #keep only stocked_list of delta indexes    
  hessian_lambda_theta1_delta <- hessian_lambda_theta1_delta[,stocked_list]
  hessian_lambda_delta_sq <- hessian_lambda_delta_sq[stocked_list,stocked_list]
  hessian_lambda_delta_sq <- keep_lower_traingular_matrix(hessian_lambda_delta_sq)
  hessian_lambda_delta_theta1 <- t(hessian_lambda_theta1_delta)
  hessian_lambda_theta1_delta[,] <- 0 #since it is all in upper triangle portion
  
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

eval_error_xi_model5 <- function(deltain, theta1, wdcMerged, points) {  
  if(!identical(order(user_serv_lvl$st_tw_index), c(1:nrow(user_serv_lvl)))) stop("user_serv_lvl not sorted")  
  stocked_list <- which(wdcMerged$stocked_out==FALSE)
  #get covariates
  if(is.null(X) | is.null(Z)) {
    list_covariates <- eval_covariates_delta_reg(deltain,theta1,wdcMerged,points)
    X <<- list_covariates$X
    Z <<- list_covariates$Z
  }
    
  #theta_2 is given by, Inv(X'.W'.Z.A_N.Z'.W.X).(X'.W'.Z.A_N.Z'.W.delta)
  #W above can be replaced by W/|W|
  
  #Z_weighted equal to W'.Z
  Z_weighted <- (Z * (wdcMerged$obs_weight[stocked_list]))/sum(wdcMerged$obs_weight[stocked_list])
  if(is.null(weighing_GMM_mat)) {
    A_N <- diag(ncol(Z))    
  } else {
    A_N <- weighing_GMM_mat
  }
  
  theta_2 <- solve(t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted) %*% X) %*%
    (t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted) %*% deltain)
  
  xi <- deltain - (X %*%  theta_2) 
  
  return(list("theta2"=theta_2,
              "theta3"=NULL,              
              "Z"=Z,
              "weights"=wdcMerged$obs_weight[stocked_list],
              "xi"=xi))
}
