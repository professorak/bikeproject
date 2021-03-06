#optimal weighing matrix
eval_GMM_optimal_weighing_matrix <- function(theta1,deltain,wdcMerged,points) {
  ret <- eval_error_xi_model5(deltain,theta1,wdcMerged,points)
  Z = ret$Z
  weights = ret$weights
  xi = ret$xi
  Z_scaled = Z * weights
  #we need to calculate S_0^(-1)=1/N \sum(G_i*G_i') 
  #where G_i is column vector of moments for each observation
  #i.e. in my setting it is Z_scaled[i,]*xi[i]
  #\sum(G_i*G_i') can be also calculated as product of a 
  #matrix[G1,G2,..] i.e. where G1,G2 etc are stacked next to each on other in ith columns
  #..and its transpose.
  #we construct t(moment_vec) which is the matrix above.
  #Note the corresponding row weights are multipled to Z already in creating Z_scaled
  moment_vec <- c()
  for(i in 1:ncol(Z_scaled)) {
    moment_vec <- cbind(moment_vec,Z_scaled[,i]*xi)
  }  
  S0 = (t(moment_vec) %*% moment_vec) / sum(weights)
  S0_inv <- solve(S0)
  return(S0_inv)  
}

eval_moments_variance <- function(theta1,deltain,wdcMerged,points) {
  ret <- eval_error_xi_model5(deltain,theta1,wdcMerged,points)
  Z = ret$Z
  weights = ret$weights
  xi = ret$xi
  Z_scaled = Z * weights
  #we need to calculate S_0^(-1)=1/N \sum(G_i*G_i') 
  #where G_i is column vector of moments for each observation
  #i.e. in my setting it is Z_scaled[i,]*xi[i]
  #\sum(G_i*G_i') can be also calculated as product of a 
  #matrix[G1,G2,..] i.e. where G1,G2 etc are stacked next to each on other in ith columns
  #..and its transpose.
  #we construct t(moment_vec) which is the matrix above.
  #Note the corresponding row weights are multipled to Z already in creating Z_scaled
  moment_vec <- c()
  for(i in 1:ncol(Z_scaled)) {
    moment_vec <- cbind(moment_vec,Z_scaled[,i]*xi)
  }  
  S0 = (t(moment_vec) %*% moment_vec) / sum(weights)  
  return(S0)  
}

eval_theta_variance_optimalGMM <- function(theta1,deltain,wdcMerged,points,active_coef=NULL) {  
  delta_all <- rep(-30, nrow(wdcMerged))
  delta_all[which(wdcMerged$stocked_out==F)] <- deltain 
  if(is.null(active_coef)) {
    active_coef <- c(1:(length(theta1)-1))
  }
  #gradient of moment conditions (G_hat)  
  ret <- eval_error_xi_model5(deltain,theta1,wdcMerged,points)
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
  
  S0 <- eval_moments_variance(theta1,deltain,wdcMerged,points)  
  if(is.null(weighing_GMM_mat)) {
    A_N <- diag(length_eta)    
  } else {
    A_N <- weighing_GMM_mat
  }

  var_p1 <- solve(t(G_hat) %*% A_N %*% G_hat)
  theta1_variance <- (var_p1 %*% t(G_hat) %*% A_N %*% S0 %*% A_N %*% G_hat %*% var_p1)/sum(weights)

  return(theta1_variance)
}

eval_theta_variance <- function(theta1,deltain,wdcMerged,points,active_coef=NULL) {  
  delta_all <- rep(-30, nrow(wdcMerged))
  delta_all[which(wdcMerged$stocked_out==F)] <- deltain 
  if(is.null(active_coef)) {
    active_coef <- c(1:(length(theta1)-1))
  }
  
  #gradient of moment conditions (G_hat)  
  ret <- eval_error_xi_model5(deltain,theta1,wdcMerged,points)
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

