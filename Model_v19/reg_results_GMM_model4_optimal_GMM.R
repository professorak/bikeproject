#get data
source("eval_obj_func_GMM_model4.R")
source("GetDataGooglePlaces.R")

#run
source("temp_optimalGMM.R")

weighing_GMM_mat <<- NULL

source("eval_func_2_cpp_MPEC.R")
source("eval_obj_func_GMM_model4.R")
length_stklist <- length(which(wdcMerged$stocked_out==F))
theta <- c(-4 ,0, 568.87, 4.13, 2532, 2.89, 9.08)
#choose starting values for deltain_stin
delta_all <- compute_delta_list_cntrt_map_new(c(theta[1],0,theta[-1]),wdcMerged, points)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
# #overriding above null setting as GMM values are quite small for identity weighing_GMM_mat setting
# weighing_GMM_mat <<- eval_GMM_optimal_weighing_matrix(c(theta[1],0,theta[-1]),deltain_stin,wdcMerged,
#                                                       points)

#########
Z <- eval_error_xi_model4(c(theta[1],0,theta[-1]),deltain_stin,wdcMerged,points)$Z
stocked_list <- which(wdcMerged$stocked_out==FALSE)
Z_sqweighted <- (Z * sqrt(wdcMerged$obs_weight[stocked_list]))/sum(wdcMerged$obs_weight[stocked_list])
print("kappa(t(Z_sqweighted) %*% Z_sqweighted): ")
print(kappa(t(Z_sqweighted) %*% Z_sqweighted))
weighing_GMM_mat <- solve(t(Z_sqweighted) %*% Z_sqweighted)
#########

params <- c(theta,deltain_stin)
lb = c(-100,0,rep(-1000,5),rep(-20, length_stklist))
ub = c(100,0,rep(1000000,5),rep(10, length_stklist))
constraint_lb <- c(rep(0, length_stklist), 100000)
constraint_ub <- c(rep(0, length_stklist), 100000)

time <- proc.time()
eval_jac_g_structure_val <- eval_jac_g_structure(params, wdcMerged, points, length(theta))
time <- proc.time() - time
print("time: ")
print(time)
#W_optimal <<- eval_weighting_mat_new(c(-10),wdcMerged,points) 

opts <- list("tol"=1.0e-6,
             "print_info_string"='yes'
#              ,"file_print_level"=6,
#              "output_file"="iterinfo_GMM_model4.out"
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
time_1 <- proc.time() - time
print("time: ")
print(time_1)
res_save <- res

params <- res$solution 
params_save <- params
res_save <- res
theta <- params[c(1:length(theta))]
theta_save <- theta
deltain_stin <- params[length(theta)+c(1:length_delta)]
theta1 <- c(theta[1],0,theta[-1])  
eval_xi <- eval_error_xi_model4 (theta1, deltain_stin, wdcMerged,points)
theta2_save <- eval_xi$theta2  
theta2_save

source("temp_optimalGMM_model4_allstderr.R")
active_coef <- c(1:7)[-c(2)]
active_coef_nointercept <- c(1:7)[-c(2,4)]
var_covar_theta_save <- eval_theta_variance(theta1,deltain_stin,wdcMerged,points,active_coef)
round(theta[active_coef_nointercept]/sqrt(diag(var_covar_theta_save)),2)


weighing_GMM_mat <<- eval_GMM_optimal_weighing_matrix(theta1,deltain_stin,wdcMerged,
                                                     points)
##########################
eval_jac_g_structure_val <- eval_jac_g_structure(params, wdcMerged, points, length(theta))
##########################
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

print(res_save$iterations + res$iterations)
print(time + time_1)


var_covar_theta <- eval_theta_variance(theta1,deltain_stin,wdcMerged,points,active_coef)
round(theta[active_coef_nointercept]/sqrt(diag(var_covar_theta)),2)

st_err <- sqrt(diag(var_covar_theta))
round(correlation_theta_estimates <- diag(1/st_err) %*% var_covar_theta %*% diag(1/st_err),3)

##########
print("Step 1 output: ")
#coef
theta_save
res_save$objective
theta2_save
round(theta_save[active_coef_nointercept]/sqrt(diag(var_covar_theta_save)),2)
res_save$iterations
time_1
print("Step 2 output: ")
theta
res$objective
theta2
round(theta[active_coef_nointercept]/sqrt(diag(var_covar_theta)),2)
res$iterations + res_save$iterations
time + time_1


