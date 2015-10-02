#deltain <- rep(-3, nrow(wdcMerged))
theta1 <- c(-4 ,0, 1, 1, 1)  
theta <- c(-4 , 1, 1, 1)  
# a <- eval_grad_constraints_MPEC(deltain,theta1, wdcMerged, points)
# b <- eval_grad_structure_constraints_MPEC(deltain,theta1, wdcMerged, points)

#
deltain_stin <- rnorm(length(which(wdcMerged$stocked_out==F)), -3,1)
params <- c(theta,deltain_stin)
length_theta=length(theta)
d_obj <- eval_obj_GMM_MPEC_obj(params, wdcMerged, points, length(theta))
d_grad <- eval_obj_GMM_MPEC_grad(params, wdcMerged, points, length(theta))
summary(d_grad)

######################################################
params <- c(theta,deltain_stin)
res <- check.derivatives(
  .x=params,
  func=eval_obj_GMM_MPEC_obj,
  func_grad=eval_obj_GMM_MPEC_grad,
  check_derivatives_print='all',  
  wdcMerged=wdcMerged, 
  points=points,
  length_theta=length(theta)
)
if(length(which(res$flag_derivative_warning==TRUE))>0) stop("eval_obj_GMM_list_extended_new test failed")

###################################################
#constraints
con_val <- eval_constraints_MPEC(params, wdcMerged, points, length(theta))
length(con_val)
con_val_grad <- eval_grad_constraints_MPEC(params, wdcMerged, points, length(theta))
length(con_val_grad)
length(which(con_val_grad!=0))
con_strcut <- eval_grad_structure_constraints_MPEC(params, wdcMerged, points, length(theta))
length(unlist(con_strcut))
length(con_strcut)
max(unlist(con_strcut))
#objective function
#eval_obj <- eval_obj_GMM_MPEC_obj(params, wdcMerged, points, length(theta))
eval_obj_grad <- eval_obj_GMM_MPEC_grad(params, wdcMerged, points, length(theta))
length(eval_obj_grad)

eval_grad_constraints_MPEC_test <- function(params, wdcMerged, points, length_theta) {
  con_val_grad <- eval_grad_constraints_MPEC(params, wdcMerged, points, length_theta)
  con_strcut <- eval_grad_structure_constraints_MPEC(params, wdcMerged, points, length_theta)
  p <- getfullfromsparsematrix(con_strcut, con_val_grad, NULL) 
  return(p)  
}

con_val_grad_test <- eval_grad_constraints_MPEC_test(params, wdcMerged, points, length(theta))

params <- c(theta,deltain_stin)
res <- check.derivatives(
  .x=params,
  func=eval_constraints_MPEC,
  func_grad=eval_grad_constraints_MPEC_test,
  check_derivatives_print='all',  
  wdcMerged=wdcMerged, 
  points=points,
  length_theta=length(theta)
)
if(length(which(res$flag_derivative_warning==TRUE))>0) stop("eval_obj_GMM_list_extended_new test failed")

################################

###


