eval_theta_variance_temp <- function(theta1,deltain,wdcMerged,points,active_coef) {  
  delta_all <- rep(-30, nrow(wdcMerged))
  delta_all[which(wdcMerged$stocked_out==F)] <- deltain 
  
  #gradient of moment conditions (G_hat)  
  ret <- eval_error_xi_model4(theta1,deltain,wdcMerged,points)
  Z = ret$Z
  weights = ret$weights
  xi = ret$xi
  Z_scaled = Z * weights
  
  grad_delta_theta_full <- eval_grad_delta_theta_full(theta1, delta_all, wdcMerged, points)
  grad_delta_theta_full <- grad_delta_theta_full[which(wdcMerged$stocked_out==F),-2] #removing random coef.
  grad_delta_theta_full <- grad_delta_theta_full[,active_coef] #removing inactive den coef.
  
  G_hat <- t(Z_scaled) %*% grad_delta_theta_full/sum(weights)  
  
  #check first order condition
  #   G <- t(Z_scaled) %*% xi/sum(weights)  
  #   grad_theta1 <- t(G_hat) %*% G #include wieghing matrix here.
  ############################  
  
  S0_inv <- eval_GMM_optimal_weighing_matrix(theta1,deltain,wdcMerged,points)  
  theta1_variance <- solve(t(G_hat) %*% S0_inv %*% G_hat)/sum(weights)
  return(theta1_variance)
}


active_coef <- c(1,3,4,5,6,7)
var_covar_theta <- eval_theta_variance_temp(theta1,deltain_stin,wdcMerged,points, active_coef)
print(round(theta[active_coef]/sqrt(diag(var_covar_theta)),2))

st_err <- sqrt(diag(var_covar_theta))
correlation_theta_estimates <- diag(1/st_err) %*% var_covar_theta %*% diag(1/st_err)
print(round(correlation_theta_estimates,2))
