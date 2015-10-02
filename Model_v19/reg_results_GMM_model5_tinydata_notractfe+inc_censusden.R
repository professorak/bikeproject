#get data
source("eval_obj_GMM_model5.R")
source("GetDataGooglePlaces_tiny.R")
source("eval_func_2_cpp_MPEC_model5_nonsparse.R")

#run
print_iter_values <- 1
source("temp_optimalGMM_model5.R")


source("eval_covariates_delta_reg_notractfe+inc_censusden.R")
weighing_GMM_mat <<- NULL


length_stklist <- length(which(wdcMerged$stocked_out==F))
#theta <- c(-4 ,0, 146.83, 20, 203.69, 0.42, 9.41)
theta <- c(-4 ,1, 568.87, 4.13, 2532, 2.89, 9.08)
#choose starting values for deltain_stin
delta_all <- compute_delta_list_cntrt_map_new(c(theta[1],0,theta[-1]),wdcMerged, points)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
# #overriding above null setting as GMM values are quite small for identity weighing_GMM_mat setting
# weighing_GMM_mat <<- eval_GMM_optimal_weighing_matrix(c(theta[1],0,theta[-1]),deltain_stin,wdcMerged,
#                                                       points)
length_eta <- ncol(eval_covariates_delta_reg(deltain_stin,c(theta[1],0,theta[-1]),wdcMerged,points)$Z)

#########
Z <- eval_covariates_delta_reg(deltain_stin,c(theta[1],0,theta[-1]),wdcMerged,points)$Z
stocked_list <- which(wdcMerged$stocked_out==FALSE)
Z_sqweighted <- (Z * sqrt(wdcMerged$obs_weight[stocked_list]))/sum(wdcMerged$obs_weight[stocked_list])
print("kappa(t(Z_sqweighted) %*% Z_sqweighted): ")
print(kappa(t(Z_sqweighted) %*% Z_sqweighted))
weighing_GMM_mat <- solve(t(Z_sqweighted) %*% Z_sqweighted)
#########

eta <- rep(0,length_eta)

eta <- c(eval_constraints_moment_conditions(deltain_stin, c(theta[1],0,theta[-1]), eta, wdcMerged, points))
length_theta <- length(theta)
length_delta <- length(deltain_stin)
params <- c(theta,deltain_stin,eta)
get_total_density(params, wdcMerged, points)
tot_density <<- 100000
target_density_constraint <<- 100000
  
lb = c(-100,0,0,0,0,0,0,rep(-30, length_delta), rep(-1.0e19,length_eta))
ub = c(0,1000000,1000000,1000000,1000000,1000000,1000000,rep(0, length_delta), rep(1.0e19,length_eta))
constraint_lb <- c(rep(0, length_delta+length_eta), target_density_constraint) 
constraint_ub <- c(rep(0, length_delta+length_eta), target_density_constraint)

eval_jac_g_structure_val <<- eval_jac_g_structure(params, wdcMerged, points, length(theta), length(deltain_stin), length_eta)
eval_hess_struct <<- eval_hessian_structure(params,wdcMerged, points, length(theta), length(deltain_stin), length_eta)

#W_optimal <<- eval_weighting_mat_new(c(-10),wdcMerged,points) 

opts <- list("tol"=1.0e-6,
             "print_info_string"='yes'
#              ,"file_print_level"=6,
#              "output_file"="iterinfo_GMM_model4_tiny.out"
)

time <- proc.time()
res <- ipoptr( x0=params,
               lb=lb,
               ub=ub,
               eval_f=eval_obj_GMM_model5_obj, 
               eval_grad_f=eval_obj_GMM_model5_grad,
               eval_g=eval_g,
               eval_jac_g=eval_jac_g,
               eval_jac_g_structure=eval_jac_g_structure_val,
               eval_h=eval_hessian,
               eval_h_structure=eval_hess_struct,
               constraint_lb=constraint_lb,
               constraint_ub=constraint_ub,
               wdcMerged=wdcMerged, 
               points=points, 
               length_theta=length(theta),
               length_delta=length_delta,
               length_eta=length_eta,
               opts=opts
) 
time_1 <- proc.time() - time
print("time: ")
print(time_1)

params <- res$solution 
params_save <- params
res_save <- res
theta <- params[c(1:length(theta))]
deltain_stin <- params[length(theta)+c(1:length_delta)]
eta <- params[length(theta)+length_delta+c(1:length_eta)]
theta1 <- c(theta[1],0,theta[-1])  

# var_covar_theta <- eval_theta_variance(theta1,deltain_stin,wdcMerged,points)



###
#second step to do optimal GMM
weighing_GMM_mat <<- eval_GMM_optimal_weighing_matrix(theta1,deltain_stin,wdcMerged,
                                                      points)
##############
eval_jac_g_structure_val <<-  eval_jac_g_structure(params, wdcMerged, points, length(theta), length(deltain_stin), length_eta)
eval_hess_struct <<-  eval_hessian_structure(params,wdcMerged, points, length(theta), length(deltain_stin), length_eta)
##############


opts <- list("tol"=1.0e-6,
             "print_info_string"='yes'
             #              ,"file_print_level"=6,
             #              "output_file"="iterinfo_GMM_model4_tiny.out"
)
time <- proc.time()
res <- ipoptr( x0=params,
               lb=lb,
               ub=ub,
               eval_f=eval_obj_GMM_model5_obj, 
               eval_grad_f=eval_obj_GMM_model5_grad,
               eval_g=eval_g,
               eval_jac_g=eval_jac_g,
               eval_jac_g_structure=eval_jac_g_structure_val,
               eval_h=eval_hessian,
               eval_h_structure=eval_hess_struct,
               constraint_lb=constraint_lb,
               constraint_ub=constraint_ub,
               wdcMerged=wdcMerged, 
               points=points, 
               length_theta=length(theta),
               length_delta=length_delta,
               length_eta=length_eta,
               opts=opts
) 
time <- proc.time() - time
print("time: ")
print(time)

params <- res$solution 


theta <- params[c(1:length(theta))]
deltain_stin <- params[length(theta)+c(1:length_delta)]
theta1 <- c(theta[1],0,theta[-1])  
eval_xi <- eval_error_xi_model5 (deltain_stin, theta1, wdcMerged,points)
theta2 <- eval_xi$theta2  
theta2
theta3 <- eval_xi$theta3  
theta3

print(res_save$iterations + res$iterations)
print(time + time_1)



