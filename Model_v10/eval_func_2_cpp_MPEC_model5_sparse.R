eval_g <- function(params, wdcMerged, points, length_theta, length_delta, length_eta) {
  print("In eval_g: ")
  ptm <- proc.time()  
  theta <- params[c(1:length_theta)]
  deltain_stin <- params[c((length_theta+1):(length_theta+length_delta))]
  eta <- params[c((length_theta+length_delta+1):(length_theta+length_delta+length_eta))]
  theta1 <- c(theta[1],0,theta[-1])  
  
  ret <- c(eval_constraints_MPEC(deltain_stin, theta1, wdcMerged, points),
           eval_constraints_moment_conditions(deltain_stin, theta1, eta, wdcMerged, points),
           get_total_density(params, wdcMerged, points))
  print("eval time:")
  print(proc.time()-ptm)
  return(ret)  
}

eval_constraints_MPEC <- function(deltain_stin, theta1, wdcMerged, points) {
  
  #expand deltain to all observations, it is currently #stocked in observations
  deltain_all <- rep(-30, nrow(wdcMerged))
  deltain_all[which(wdcMerged$stocked_out==F)] <- deltain_stin
  
  tw_group_list <- unique(wdcMerged$tw_group)
  obj <- c()
  for(i in 1:length(tw_group_list)) {
    tw_groupin = tw_group_list[i]
    wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
    deltain_tw = deltain_all[which(wdcMerged$tw_group==tw_groupin)]    
    obj <- c(obj, 
             eval_constraints_MPEC_tw_groupin(deltain_tw, theta1, wdcMergedday, 
                                              points, tw_groupin))
  }
  return(obj)
}

eval_constraints_MPEC_tw_groupin <- function(deltain, theta1, wdcMergedday, points, tw_groupin) {
  stocked_list <- which(wdcMergedday$stocked_out==FALSE)
  dem_T <- eval_lambda_new(deltain, theta1, wdcMergedday, points, tw_groupin)[stocked_list]
  dem_hat_T <- wdcMergedday$out_dem_sum[stocked_list]
  obj <- dem_T-dem_hat_T
  return(obj)
}

eval_constraints_moment_conditions <- function(deltain, theta1, eta, wdcMerged, points) {  
  if(!identical(order(user_serv_lvl$st_tw_index), c(1:nrow(user_serv_lvl)))) stop("user_serv_lvl not sorted")  
  stocked_list <- which(wdcMerged$stocked_out==FALSE)
  #get covariates
  list_covariates <- eval_covariates_delta_reg(deltain,theta1,wdcMerged,points)
  X <- list_covariates$X
  Z <- list_covariates$Z
  
  
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
  
  return(obj_moments - eta)  
}

eval_covariates_delta_reg <- function(deltain,theta1,wdcMerged,points) {
  stocked_list <- which(wdcMerged$stocked_out==FALSE)
  #compute attributes matrix X and instruments matrix Z
  #list of X attributes
  #1s, sevice level, diXw fixed effects, month fixed effects, weather fixed effects 
  #list of Z attributes
  #1s, diXw fixed effects, month fixed effects, weather fixed effects, two service level instruments, local density attributes of a station, 
  #stockout indicator for nearby station
  #of the Z's need to get --local density attributes of a station, and stockout indicator for nearby station, while generating data, 
  #one entry for each data row.
  
  #generating the common attributes of X and Z. Calling it Xbase
  #1s, diXw fixed effects, month fixed effects, weather fixed effects
  if(length(levels(wdcMerged$tract_tw_fac))>1 & length(levels(wdcMerged$week_fac))>1) {
    reg_formula <- "~week_fac + tract_tw_fac"
    if(length(levels(wdcMerged$Conditions))>1) paste(reg_formula, "+ Conditions")
    length_levels <- c(length(levels(wdcMerged$Conditions)), length(levels(wdcMerged$Temperature_group)),
                       length(levels(wdcMerged$Humidity_high)), length(levels(wdcMerged$Wind.Speed_high)))
    reg_terms <- c("Conditions", "Temperature_group", "Humidity_high", "Wind.Speed_high")[which(length_levels>1)]
    reg_formula <- paste(reg_formula, paste(reg_terms, collapse = " + "), sep = " + ")
    
    Xbase <- model.matrix(as.formula(reg_formula), data=wdcMerged[stocked_list,])    
    
    if(length(which(colnames(Xbase) %like% "Conditions")) != length(levels(wdcMerged$Conditions))-1) {
      stop("in eval_error_xi_sl_model4 one conditions number not automatically removed")      
    }
    if(length(which(colnames(Xbase) %like% "Temperature_group")) != length(levels(wdcMerged$Temperature_group))-1) {
      stop("in eval_error_xi_sl_model4 one Temperature_group number not automatically removed")      
    }
    if(length(which(colnames(Xbase) %like% "Humidity_high")) != length(levels(wdcMerged$Humidity_high))-1) {
      stop("in eval_error_xi_sl_model4 one Humidity_high number not automatically removed")      
    }
    if(length(which(colnames(Xbase) %like% "Wind.Speed_high")) != length(levels(wdcMerged$Wind.Speed_high))-1) {
      stop("in eval_error_xi_sl_model4 one Wind.Speed_high number not automatically removed")      
    }
    if(length(which(colnames(Xbase) %like% "tract_tw_fac")) != length(levels(wdcMerged$tract_tw_fac))-1) {
      stop("in eval_error_xi_sl_model4 one tract_tw_fac number not automatically removed")      
    }
  } else if (length(levels(wdcMerged$tract_tw_fac))>1) {
    stop("error in eval_error_xi")
  } else if (length(levels(wdcMerged$week_fac))>1) {
    stop("error in eval_error_xi")
  } else {
    stop("error in eval_error_xi")
  }  
  #drop from Xbase columns which have no non-zero entry
  #rowsum is quite sure way to test
  Xbase <- Xbase[,-which(colSums(Xbase)==0)] 
  
  #add service level vector to convert Xbase to X
  X <- cbind(Xbase, serv_lvl=wdcMerged$serv_lvl[stocked_list])
  
  #add two service level instruments, local density attributes of a station, stockout indicator for nearby station
  #to get Z
  Zbase <- as.matrix(wdcMerged[stocked_list,c("instr_serv_lvl","serv_lvl_neighbours","census_density","metro_den_on_1","metro_den_on_2","metro_den_off_1","metro_den_off_2"
                                              ,"bus_den_1","bus_den_2","googleplaces_den_1","googleplaces_den_2","sto_nearby")])
  Z <- cbind(Xbase, Zbase)
  
  Z <- colNormalize(Z) #so that column sums are equal to 1.
  
  return(list("X"=X,              
              "Z"=Z))
}

eval_jac_g <- function(params, wdcMerged, points, length_theta, length_delta, length_eta) {
  print("In eval_jac_g: ")
  ptm <- proc.time()  
  theta <- params[c(1:length_theta)]
  deltain_stin <- params[c((length_theta+1):(length_theta+length_delta))]
  eta <- params[c((length_theta+length_delta+1):(length_theta+length_delta+length_eta))]
  theta1 <- c(theta[1],0,theta[-1])  
  
  grad_total_density_full <- rep(0,length(params))
  grad_total_density_full[c(1:length_theta)] <- get_grad_total_density(params, wdcMerged, points)
  ret <- rbind(eval_grad_constraints_MPEC(deltain_stin, theta1, eta, wdcMerged, points),
               eval_grad_constraints_moment_conditions(deltain_stin, theta1, eta, wdcMerged, points),
               grad_total_density_full)
  
  print("eval time:")
  print(proc.time()-ptm)
  return(ret)
}
