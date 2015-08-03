# wdcMerged$out_dem_sum[which(wdcMerged$out_dem_sum<=0.01 & 
#                               wdcMerged$stocked_out==FALSE)] <- 0.01

wdcMerged_nstout <- wdcMerged[which(wdcMerged$stocked_out==F),]

length(which(wdcMerged_nstout$out_dem_sum<=0.01))
nrow(wdcMerged_nstout)
sum(wdcMerged_nstout$obs_weight[which(wdcMerged_nstout$out_dem_sum<=0.01)])
sum(wdcMerged_nstout$obs_weight)
sum(wdcMerged_nstout$obs_weight[which(wdcMerged_nstout$out_dem_sum<=0.01)])/sum(wdcMerged_nstout$obs_weight)

summary(wdcMerged_nstout$obs_weight[which(wdcMerged_nstout$out_dem_sum < 0.001)])
summary(wdcMerged_nstout$obs_weight[which(wdcMerged_nstout$out_dem_sum >= 0.001)])

####################
wdcMerged$out_dem_sum[which(wdcMerged$out_dem_sum<=0.1 & 
                              wdcMerged$stocked_out==FALSE)] <- 0.0001
params <- c(-3.458331564)
res <- check.derivatives(
  .x=params,
  func=eval_obj_GMM_list_extended_new_obj,
  func_grad=eval_obj_GMM_list_extended_new_grad,
  check_derivatives_print='all',  
  wdcMerged=wdcMerged, 
  points=points  
)
if(length(which(res$flag_derivative_warning==TRUE))>0) stop("eval_obj_GMM_list_extended_new test failed")  

wdcMerged$out_dem_sum[which(wdcMerged$out_dem_sum<=0.1 & 
                              wdcMerged$stocked_out==FALSE)] <- 0.1
params <- c(-3.458331564)
res <- check.derivatives(
  .x=params,
  func=eval_obj_GMM_list_extended_new_obj,
  func_grad=eval_obj_GMM_list_extended_new_grad,
  check_derivatives_print='all',  
  wdcMerged=wdcMerged, 
  points=points  
)
if(length(which(res$flag_derivative_warning==TRUE))>0) stop("eval_obj_GMM_list_extended_new test failed")  

##############################################
#testing effect on xi
theta_vec <- c(-3,-4)
xi_norm_vec <- c()
grad_xi_norm_vec <- c()
gmm_objective_vec <- c()
gmm_gradient_vec  <- c()

for(i in c(1:length(theta_vec))) {
  theta1 <- c(theta_vec[i],0)
  ret <- eval_error_xi(theta1,wdcMerged,points)
  xi_norm=ret$xi
  # theta2=ret$theta2
  # X1=ret$X1
  # Z=ret$Z
  # W=ret$W
  grad_xi_norm = ret$grad_delta_theta
  gmm_objective = c(t(xi_norm ) %*% xi_norm )/length(xi_norm)
  gmm_gradient <- 2*t(xi_norm ) %*% grad_xi_norm/length(xi_norm)  
  
  xi_norm_vec <- cbind(xi_norm_vec, xi_norm)
  grad_xi_norm_vec <- cbind(grad_xi_norm_vec, grad_xi_norm)
  gmm_objective_vec <- cbind(gmm_objective_vec,gmm_objective)
  gmm_gradient_vec <- cbind(gmm_gradient_vec,gmm_gradient)
}




