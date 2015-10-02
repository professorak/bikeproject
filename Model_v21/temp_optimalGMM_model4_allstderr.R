eval_theta_variance_optimalGMM <- function(theta1,deltain,wdcMerged,points,active_coef=NULL) {  
  delta_all <- rep(-30, nrow(wdcMerged))
  delta_all[which(wdcMerged$stocked_out==F)] <- deltain 
  if(is.null(active_coef)) {
    active_coef <- c(1:(length(theta1)))
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
  ret <- eval_error_xi_model4(theta1,deltain,wdcMerged,points)
  weights = ret$weights
  xi = ret$xi
  Z_scaled = Z * weights
  Z_weighted <- (Z * weights)/sum(weights)
  
  grad_delta_theta_full <- eval_grad_delta_theta_full(theta1, delta_all, wdcMerged, points)
  grad_delta_theta_full <- grad_delta_theta_full[which(wdcMerged$stocked_out==F),-2] #removing random coef.
  
  ############
  grad_total_density <- get_grad_total_density(params, wdcMerged, points)
  grad_total_density_multipliers <- -grad_total_density/grad_total_density[density_intercept_col]
  grad_delta_theta_full_constrained <- c()
  for(i in 1:ncol(grad_delta_theta_full)) {
    grad_delta_theta_full_constrained <- cbind(grad_delta_theta_full_constrained, 
        grad_delta_theta_full[,i]+grad_total_density_multipliers[i]*grad_delta_theta_full[,density_intercept_col])
  }
  grad_delta_theta_full <- grad_delta_theta_full_constrained
  
  ###########
  
  grad_delta_theta_full <- grad_delta_theta_full[,active_coef] #removing inactive den coef.
  grad_delta_theta_full <- grad_delta_theta_full[,-which(active_coef==(density_intercept_col))]
  
  #G_hat <- t(Z_scaled) %*% grad_delta_theta_full/sum(weights)  
  
  grad_delta <- t(Z_weighted) - t(Z_weighted) %*% X %*% solve(t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted) %*% X) %*%
    (t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted))
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
  delta_all <- rep(-30, nrow(wdcMerged))
  delta_all[which(wdcMerged$stocked_out==F)] <- deltain 
  if(is.null(active_coef)) {
    active_coef <- c(1:(length(theta1)))
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
  ret <- eval_error_xi_model4(theta1,deltain,wdcMerged,points)
  weights = ret$weights
  xi = ret$xi
  Z_scaled = Z * weights
  Z_weighted <- (Z * weights)/sum(weights)
  
  grad_delta_theta_full <- eval_grad_delta_theta_full(theta1, delta_all, wdcMerged, points)
  grad_delta_theta_full <- grad_delta_theta_full[which(wdcMerged$stocked_out==F),] 
  ############
  grad_total_density <- get_grad_total_density(params, wdcMerged, points)
  grad_total_density_multipliers <- -grad_total_density/grad_total_density[density_intercept_col]
  grad_delta_theta_full_constrained <- c()
  for(i in 1:ncol(grad_delta_theta_full)) {
    grad_delta_theta_full_constrained <- cbind(grad_delta_theta_full_constrained, 
                                               grad_delta_theta_full[,i]+grad_total_density_multipliers[i]*grad_delta_theta_full[,density_intercept_col])
  }
  grad_delta_theta_full <- grad_delta_theta_full_constrained
  
  ###########
  
  grad_delta_theta_full <- grad_delta_theta_full[,active_coef] #removing inactive den coef.
  grad_delta_theta_full <- grad_delta_theta_full[,-which(active_coef==(density_intercept_col))]
  
  #G_hat <- t(Z_scaled) %*% grad_delta_theta_full/sum(weights)  
  
  grad_delta <- t(Z_weighted) - t(Z_weighted) %*% X %*% solve(t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted) %*% X) %*%
    (t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted))
  G_hat_theta1 <- grad_delta %*% grad_delta_theta_full
  G_hat_theta2 <- - t(Z_weighted) %*% X
  G_hat <- cbind(G_hat_theta1, G_hat_theta2)
  
  #check first order condition
  #   G <- t(Z_scaled) %*% xi/sum(weights)  
  #   grad_theta1 <- t(G_hat) %*% G #include wieghing matrix here.
  ############################  
  
  S0 <- eval_moments_variance(theta1,deltain,wdcMerged,points)  
  print(paste0("kappa t(G_hat) %*% A_N %*% G_hat: ",kappa(t(G_hat) %*% A_N %*% G_hat)))
  
  #below returns singular values
  #   var_p1 <- solve(t(G_hat) %*% A_N %*% G_hat)
  #   theta1_variance <- (var_p1 %*% t(G_hat) %*% A_N %*% S0 %*% A_N %*% G_hat %*% var_p1)/sum(weights)
  var_p1 <- solve(t(G_hat_theta1) %*% A_N %*% G_hat_theta1)
  theta1_variance <- (var_p1 %*% t(G_hat_theta1) %*% A_N %*% S0 %*% A_N %*% G_hat_theta1 %*% var_p1)/sum(weights)
  
  return(theta1_variance)
}

eval_theta2_variance <- function(theta1,deltain,wdcMerged,points,active_coef=NULL) {  
  delta_all <- rep(-30, nrow(wdcMerged))
  delta_all[which(wdcMerged$stocked_out==F)] <- deltain 
  if(is.null(active_coef)) {
    active_coef <- c(1:(length(theta1)))
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
  ret <- eval_error_xi_model4(theta1,deltain,wdcMerged,points)
  weights = ret$weights
  xi = ret$xi
  Z_scaled = Z * weights
  Z_weighted <- (Z * weights)/sum(weights)
  
  grad_delta_theta_full <- eval_grad_delta_theta_full(theta1, delta_all, wdcMerged, points)
  grad_delta_theta_full <- grad_delta_theta_full[which(wdcMerged$stocked_out==F),] 
  ############
  grad_total_density <- get_grad_total_density(params, wdcMerged, points)
  grad_total_density_multipliers <- -grad_total_density/grad_total_density[density_intercept_col]
  grad_delta_theta_full_constrained <- c()
  for(i in 1:ncol(grad_delta_theta_full)) {
    grad_delta_theta_full_constrained <- cbind(grad_delta_theta_full_constrained, 
                                               grad_delta_theta_full[,i]+grad_total_density_multipliers[i]*grad_delta_theta_full[,density_intercept_col])
  }
  grad_delta_theta_full <- grad_delta_theta_full_constrained
  
  ###########
  
  grad_delta_theta_full <- grad_delta_theta_full[,active_coef] #removing inactive den coef.
  grad_delta_theta_full <- grad_delta_theta_full[,-which(active_coef==(density_intercept_col))]
  
  #G_hat <- t(Z_scaled) %*% grad_delta_theta_full/sum(weights)  
  
  grad_delta <- t(Z_weighted) - t(Z_weighted) %*% X %*% solve(t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted) %*% X) %*%
    (t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted))
  G_hat_theta1 <- grad_delta %*% grad_delta_theta_full
  G_hat_theta2 <- - t(Z_weighted) %*% X
  G_hat <- cbind(G_hat_theta1, G_hat_theta2)
  
  #check first order condition
  #   G <- t(Z_scaled) %*% xi/sum(weights)  
  #   grad_theta1 <- t(G_hat) %*% G #include wieghing matrix here.
  ############################  
  
  S0 <- eval_moments_variance(theta1,deltain,wdcMerged,points)  
  print(paste0("kappa t(G_hat) %*% A_N %*% G_hat: ",kappa(t(G_hat) %*% A_N %*% G_hat)))
  
  var_p1 <- solve(t(G_hat_theta2) %*% A_N %*% G_hat_theta2)
  theta2_variance <- (var_p1 %*% t(G_hat_theta2) %*% A_N %*% S0 %*% A_N %*% G_hat_theta2 %*% var_p1)/sum(weights)
  
  return(theta2_variance)
}
