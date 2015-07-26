#get data
source("eval_obj_func_GMM_model4.R")
source("GetDataGooglePlaces.R")

#run
source("temp_optimalGMM.R")

weighing_GMM_mat <<- NULL

source("eval_func_2_cpp_MPEC.R")
source("eval_obj_func_GMM_model4.R")
length_stklist <- length(which(wdcMerged$stocked_out==F))
theta <- c(-4 ,0.51, 146.83, 0.052, 203.69, 0.42, 9.41)
#choose starting values for deltain_stin
delta_all <- compute_delta_list_cntrt_map_new(c(theta[1],0,theta[-1]),wdcMerged, points)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
# #overriding above null setting as GMM values are quite small for identity weighing_GMM_mat setting
# weighing_GMM_mat <<- eval_GMM_optimal_weighing_matrix(c(theta[1],0,theta[-1]),deltain_stin,wdcMerged,
#                                                       points)
params <- c(theta,deltain_stin)
lb = c(-20,rep(-1,6),rep(-30, length_stklist))
ub = c(2,rep(1000000,6),rep(0, length_stklist))
constraint_lb <- c(rep(0, length_stklist), 100000)
constraint_ub <- c(rep(0, length_stklist), 100000)

eval_jac_g_structure_val <- eval_jac_g_structure(params, wdcMerged, points, length(theta))
#W_optimal <<- eval_weighting_mat_new(c(-10),wdcMerged,points) 

opts <- list("tol"=1.0e-6,
             "print_info_string"='yes',
             "file_print_level"=6,
             "output_file"="iterinfo_GMM_model4.out"
)

time <- proc.time()
res <- ipoptr( x0=params,
               lb=lb,
               ub=ub,
               eval_f=eval_obj_GMM_model4_obj, 
               eval_grad_f=eval_obj_GMM_model4_grad,
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
time <- proc.time() - time
print("time: ")
print(time)
params <- res$solution 
theta <- params[c(1:length(theta))]
deltain_stin <- params[-c(1:length(theta))]
theta1 <- c(theta[1],0,theta[-1])  

weighing_GMM_mat <<- eval_GMM_optimal_weighing_matrix(theta1,deltain_stin,wdcMerged,
                                                     points)

#print(theta1)
time <- proc.time()
res <- ipoptr( x0=params,
               lb=lb,
               ub=ub,
               eval_f=eval_obj_GMM_model4_obj, 
               eval_grad_f=eval_obj_GMM_model4_grad,
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
time <- proc.time() - time
print("time: ")
print(time)
params <- res$solution 


theta <- params[c(1:length(theta))]
deltain_stin <- params[-c(1:length(theta))]
theta1 <- c(theta[1],0,theta[-1])  
eval_xi <- eval_error_xi_model4 (theta1, deltain_stin, wdcMerged,points)
theta2 <- eval_xi$theta2  
theta2
theta3 <- eval_xi$theta3  
theta3


