
source("eval_func_2_cpp_MPEC.R")
length_stklist <- length(which(wdcMerged$stocked_out==F))
theta <- c(-4 , 1, 1, 1,1,1,1)  
#choose starting values for deltain_stin
delta_all <- compute_delta_list_cntrt_map_new(c(theta[1],0,theta[-1]),wdcMerged, points)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]

params <- c(theta,deltain_stin)
lb = c(-20,rep(0.01,6),rep(-30, length_stklist))
ub = c(0,rep(1000000,6),rep(0, length_stklist))
constraint_lb <- c(rep(0, length_stklist), 100000)
constraint_ub <- c(rep(0, length_stklist), 100000)

eval_jac_g_structure_val <- eval_jac_g_structure(params, wdcMerged, points, length(theta))
#W_optimal <<- eval_weighting_mat_new(c(-10),wdcMerged,points) 

opts <- list("tol"=1.0e-4,
             "print_info_string"='yes'
)

time <- proc.time()
res <- ipoptr( x0=params,
               lb=lb,
               ub=ub,
               eval_f=eval_obj_GMM_MPEC_obj, 
               eval_grad_f=eval_obj_GMM_MPEC_grad,
               eval_g=eval_g,
               eval_jac_g=eval_jac_g,
               eval_jac_g_structure=eval_jac_g_structure_val,
               constraint_lb=constraint_lb,
               constraint_ub=constraint_ub,
               wdcMerged=wdcMerged, 
               points=points, 
               length_theta=length(theta),
               opts=opts
) 
params <- res$solution 
#print(theta1)
time <- proc.time() - time
theta <- params[c(1:length(theta))]
deltain_stin <- params[-c(1:length(theta))]
theta1 <- c(theta[1],0,theta[-1])  
eval_xi <- eval_error_xi_sl_MPEC (theta1, deltain_stin, wdcMerged,points)
theta2 <- eval_xi$theta2  
theta2
theta3 <- eval_xi$theta3  
theta3


