deltain <- deltain_stin
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

summary(grad_delta_theta_full[,1])
summary(grad_delta_theta_full[,2])



