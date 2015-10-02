eval_theta_variance_optimalGMM <- function(theta1,deltain,wdcMerged,points,active_coef=NULL) {  
  delta_all <- rep(-30, nrow(wdcMerged))
  delta_all[which(wdcMerged$stocked_out==F)] <- deltain 
  if(is.null(active_coef)) {
    active_coef <- c(1:(length(theta1)-1))
  }
  
  #gradient of moment conditions (G_hat)  
  list_covariates <- eval_covariates_delta_reg(deltain,theta1,wdcMerged,points)
  X <- list_covariates$X
  Z <- list_covariates$Z
  if(is.null(weighing_GMM_mat)) {
    A_N <- diag(ncol(Z))    
  } else {
    A_N <- weighing_GMM_mat
  }
  ret <- eval_error_xi_model5(deltain,theta1,wdcMerged,points)
  weights = ret$weights
  xi = ret$xi
  Z_scaled = Z * weights
  Z_weighted <- (Z * weights)/sum(weights)
  
  grad_delta_theta_full <- eval_grad_delta_theta_full(theta1, delta_all, wdcMerged, points)
  grad_delta_theta_full <- grad_delta_theta_full[which(wdcMerged$stocked_out==F),-2] #removing random coef.
  
  ############
  grad_total_density <- get_grad_total_density(params, wdcMerged, points)
  grad_total_density_multipliers <- -grad_total_density/grad_total_density[density_intercept_col-1]
  grad_delta_theta_full_constrained <- c()
  for(i in 1:ncol(grad_delta_theta_full)) {
    grad_delta_theta_full_constrained <- cbind(grad_delta_theta_full_constrained, 
        grad_delta_theta_full[,i]+grad_total_density_multipliers[i]*grad_delta_theta_full[,density_intercept_col-1])
  }
  grad_delta_theta_full <- grad_delta_theta_full_constrained
  
  ###########
  
  grad_delta_theta_full <- grad_delta_theta_full[,active_coef] #removing inactive den coef.
  grad_delta_theta_full <- grad_delta_theta_full[,-which(active_coef==(density_intercept_col-1))]
  
  #G_hat <- t(Z_scaled) %*% grad_delta_theta_full/sum(weights)  
  
  grad_delta <- t(Z_weighted) #this is G_hat when derivatives are taken wrt delta.
  #this takes theta2 as given, therefore gradient of xi wrt delta is identity.
  G_hat_theta1 <- grad_delta %*% grad_delta_theta_full
  G_hat_theta2 <- - t(Z_weighted) %*% X
  G_hat <- cbind(G_hat_theta1, G_hat_theta2)
  #check first order condition
  #   G <- t(Z_scaled) %*% xi/sum(weights)  
  #   grad_theta1 <- t(G_hat) %*% G #include wieghing matrix here.
  ############################  
  
  S0_inv <- eval_GMM_optimal_weighing_matrix(theta1,deltain,wdcMerged,points)  
  #theta1_variance <- solve(t(G_hat) %*% S0_inv %*% G_hat)/sum(weights) #returns singular value
  print(paste0("kappa t(G_hat) %*% S0_inv %*% G_hat: ",kappa(t(G_hat) %*% S0_inv %*% G_hat)))
  theta1_variance <- solve(t(G_hat_theta1) %*% S0_inv %*% G_hat_theta1)/sum(weights)
  return(theta1_variance)
}


eval_theta_variance <- function(theta1,deltain,wdcMerged,points,active_coef=NULL) {  
  thetafull_variance <- eval_thetafull_variance(theta1,deltain,wdcMerged,points,active_coef)
  theta1_variance <- thetafull_variance[c(1:length(active_coef_nointercept)),
                                        c(1:length(active_coef_nointercept))]
  
  return(theta1_variance)
}

eval_theta2_variance <- function(theta1,deltain,wdcMerged,points,active_coef=NULL) {  
  thetafull_variance <- eval_thetafull_variance(theta1,deltain,wdcMerged,points,active_coef)
  theta2_variance <- thetafull_variance[-c(1:length(active_coef_nointercept)),
                                        -c(1:length(active_coef_nointercept))]
  
  return(theta2_variance)
}

eval_thetafull_variance <- function(theta1,deltain,wdcMerged,points,active_coef=NULL) {  
  delta_all <- rep(-30, nrow(wdcMerged))
  delta_all[which(wdcMerged$stocked_out==F)] <- deltain 
  if(is.null(active_coef)) {
    active_coef <- c(1:(length(theta1)-1))
  }
  #gradient of moment conditions (G_hat)  
  list_covariates <- eval_covariates_delta_reg(deltain,theta1,wdcMerged,points)
  X <- list_covariates$X
  Z <- list_covariates$Z
  if(is.null(weighing_GMM_mat)) {
    A_N <- diag(ncol(Z))    
  } else {
    A_N <- weighing_GMM_mat
  }
  ret <- eval_error_xi_model5(deltain,theta1,wdcMerged,points)
  weights = ret$weights
  xi = ret$xi
  Z_scaled = Z * weights
  Z_weighted <- (Z * weights)/sum(weights)
  
  grad_delta_theta_full <- eval_grad_delta_theta_full(theta1, delta_all, wdcMerged, points)
  grad_delta_theta_full <- grad_delta_theta_full[which(wdcMerged$stocked_out==F),-2] #removing random coef.
  ############
  grad_total_density <- get_grad_total_density(params, wdcMerged, points)
  grad_total_density_multipliers <- -grad_total_density/grad_total_density[density_intercept_col-1]
  grad_delta_theta_full_constrained <- c()
  for(i in 1:ncol(grad_delta_theta_full)) {
    grad_delta_theta_full_constrained <- cbind(grad_delta_theta_full_constrained, 
                                               grad_delta_theta_full[,i]+grad_total_density_multipliers[i]*grad_delta_theta_full[,density_intercept_col-1])
  }
  grad_delta_theta_full <- grad_delta_theta_full_constrained
  
  ###########
  
  grad_delta_theta_full <- grad_delta_theta_full[,active_coef] #removing inactive den coef.
  grad_delta_theta_full <- grad_delta_theta_full[,-which(active_coef==(density_intercept_col-1))]
  
  #G_hat <- t(Z_scaled) %*% grad_delta_theta_full/sum(weights)  
  
  grad_delta <- t(Z_weighted) #this is G_hat when derivatives are taken wrt delta.
  #this takes theta2 as given, therefore gradient of xi wrt delta is identity.  #are taken wrt delta.
  G_hat_theta1 <- grad_delta %*% grad_delta_theta_full
  G_hat_theta2 <- - t(Z_weighted) %*% X
  G_hat <- cbind(G_hat_theta1, G_hat_theta2)
  
  scale_G_hat <<- sqrt(diag((t(G_hat) %*% G_hat)))
  G_hat <- scalecols(G_hat, scale_G_hat) #diving G_hat cols, equivalent to 
  #scaling up theta by  scale_G_hat_theta1.
  
  #check first order condition
  #   G <- t(Z_scaled) %*% xi/sum(weights)  
  #   grad_theta1 <- t(G_hat) %*% G #include wieghing matrix here.
  ############################  
  
  S0 <- eval_moments_variance(theta1,deltain,wdcMerged,points)  
  #   print(paste0("kappa t(G_hat) %*% A_N %*% G_hat: ",kappa(t(G_hat) %*% A_N %*% G_hat)))
  print(paste0("kappa t(G_hat) %*% A_N %*% G_hat: ",kappa(t(G_hat) %*% A_N %*% G_hat)))
  #below returns singular values
  #   var_p1 <- solve(t(G_hat) %*% A_N %*% G_hat)
  #   theta1_variance <- (var_p1 %*% t(G_hat) %*% A_N %*% S0 %*% A_N %*% G_hat %*% var_p1)/sum(weights)
  var_p1 <- solve(t(G_hat) %*% A_N %*% G_hat)
  theta1_variance <- (var_p1 %*% t(G_hat) %*% A_N %*% S0 %*% A_N %*% G_hat %*% var_p1)/sum(weights)
  
  theta1_variance <- scalecols(theta1_variance, scale_G_hat)
  theta1_variance <- theta1_variance/scale_G_hat #scale rows
  
  return(theta1_variance)
}

