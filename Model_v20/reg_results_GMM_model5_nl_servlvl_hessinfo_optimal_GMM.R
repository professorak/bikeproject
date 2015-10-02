
source("eval_obj_GMM_model5.R")
source("GetDataGooglePlaces_food.R")
source("eval_func_2_cpp_MPEC_model5_nonsparse.R")
print_iter_values <<- 1

#run
source("temp_optimalGMM_model5.R")
source("eval_covariates_delta_reg_nl_servlvl.R")
weighing_GMM_mat <<- NULL

length_stklist <- length(which(wdcMerged$stocked_out==F))
theta <- c(-4 ,1, 146.83, 20, 203.69, 0.42, 9.41)
#choose starting values for deltain_stin
#delta_all <- rnorm(nrow(wdcMerged),-10,1)
delta_all <- compute_delta_list_cntrt_map_new(c(theta[1],0,theta[-1]),wdcMerged, points)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
# #overriding above null setting as GMM values are quite small for identity weighing_GMM_mat setting
# weighing_GMM_mat <<- eval_GMM_optimal_weighing_matrix(c(theta[1],0,theta[-1]),deltain_stin,wdcMerged,
#                                                       points)
length_eta <- ncol(eval_covariates_delta_reg(deltain_stin,c(theta[1],0,theta[-1]),wdcMerged,points)$Z)

#########
list_covariates <- eval_covariates_delta_reg(deltain_stin,c(theta[1],0,theta[-1]),wdcMerged,points)
X <<- list_covariates$X
Z <<- list_covariates$Z

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

lb = c(-100,0,rep(0,5),rep(-20, length_delta), rep(-1.0e19,length_eta))
ub = c(100,1000000,1000000,1000000,1000000,1000000,1000000,rep(10, length_delta), rep(1.0e19,length_eta))
constraint_lb <- c(rep(0, length_delta+length_eta), target_density_constraint) 
constraint_ub <- c(rep(0, length_delta+length_eta), target_density_constraint)


eval_jac_g_structure_val <- eval_jac_g_structure(params, wdcMerged, points, length(theta), length(deltain_stin), length_eta)
eval_hess_struct <- eval_hessian_structure(params,wdcMerged, points, length(theta), length(deltain_stin), length_eta)
#W_optimal <<- eval_weighting_mat_new(c(-10),wdcMerged,points) 

opts <- list("tol"=1.0e-6,
             "print_info_string"='yes',
             "mumps_mem_percent"= 100
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
theta_save <- theta
deltain_stin <- params[length(theta)+c(1:length_delta)]
eta <- params[length(theta)+length_delta+c(1:length_eta)]
theta1 <- c(theta[1],0,theta[-1])  
eval_xi <- eval_error_xi_model5 (deltain_stin, theta1, wdcMerged,points)
theta2_save <- eval_xi$theta2  
theta2_save

# var_covar_theta <- eval_theta_variance(theta1,deltain_stin,wdcMerged,points)
source("temp_optimalGMM_model5_allstderr.R")
active_coef <- c(1:7)
active_coef_nointercept <- c(1:7)[-c(4)]
var_covar_theta_save <- eval_theta_variance(theta1,deltain_stin,wdcMerged,points,active_coef)
round(theta[active_coef_nointercept]/sqrt(diag(var_covar_theta_save)),2)


weighing_GMM_mat <<- eval_GMM_optimal_weighing_matrix(theta1,deltain_stin,wdcMerged,
                                                      points)

##############
eval_jac_g_structure_val <- eval_jac_g_structure(params, wdcMerged, points, length(theta), length(deltain_stin), length_eta)
eval_hess_struct <- eval_hessian_structure(params,wdcMerged, points, length(theta), length(deltain_stin), length_eta)
##############



#print(theta1)
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

var_covar_theta <- eval_theta_variance(theta1,deltain_stin,wdcMerged,points,active_coef)
round(theta[active_coef_nointercept]/sqrt(diag(var_covar_theta)),2)

st_err <- sqrt(diag(var_covar_theta))
round(correlation_theta_estimates <- diag(1/st_err) %*% var_covar_theta %*% diag(1/st_err),3)

##########print("Step 1 output: ")
#coef
print(theta_save)
print(res_save$objective)
print(theta2_save[which(rownames(theta2_save)=="serv_lvl")]/0.0006466922293)
print(theta2_save[which(rownames(theta2_save)=="serv_lvl_sq")]/0.0006466922293/0.8857722)
print(round(theta_save[active_coef_nointercept]/sqrt(diag(var_covar_theta_save)),2))
print(res_save$iterations)
print(time_1)
print("Step 2 output: ")
print(theta)
print(res$objective)
print(theta2[which(rownames(theta2)=="serv_lvl")]/0.0006466922293)
print(theta2[which(rownames(theta2)=="serv_lvl_sq")]/0.0006466922293/0.8857722)
print(round(theta[active_coef_nointercept]/sqrt(diag(var_covar_theta)),2))
print(res$iterations + res_save$iterations)
print(time + time_1)
print(round(correlation_theta_estimates <- diag(1/st_err) %*% var_covar_theta %*% diag(1/st_err),3))



