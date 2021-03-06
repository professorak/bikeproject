source("constants_mw6.R")
sink("rb17.1_1:10_full_stepf_pr6.txt",split=T)

#model with averaged station demand, no service level. Trying to get distance
#coefficient right.

get_points_density_grad_places_count_col <- NULL
source("eval_obj_GMM_model5.R")
source("GetDataGooglePlaces_aggmonth_averaged_pr6.R")

print_iter_values <<- 1

#run
source("temp_optimalGMM_model5.R")

sourceCpp("eval_func_2_new_deltaaveraged_stepfunction_steps_set7.cpp") 
#distance step at median quartile 0.300 kms
source("eval_func_2_cpp_steps50_100_200.R")
source("eval_func_2_cpp_MPEC_model5_nonsparse_steps50_100_200.R")
source("eval_func_3_cpp_new_2.86_2_steps50_100_200.R")
source("eval_covariates_delta_reg_densitycols_setb17.R") #all metro, googleplaces etc. density variables minus census density
#and catchement area instruments + service level single step
source("moredensitycols_functions.R")

colnames_theta1 <<- c("dis coef","rand coef","dis coef",
                      "density ridership","density metro","density intercept",
                      "density metro evening","local_government_office", "food",
                      "lodging","museum","movie_theater",
                      "density bus", "density tourist locations")
nonden_ceofrange <- c(1:3)
nonden_ceoflength <- 3
density_ridership_col <<- 4
density_metro_col <<- 5
density_intercept_col <<- 6
density_metro_evening_col <<- 7
#density_google_places_count <<- 10
#density_bus_col <<- 10
density_add_den_cols <<- c(13,14) #bus, tourist location


density_google_places_theta_cols <- c(8,9,10,11,12)
density_google_places_points_cols <- c(which(colnames(points)=="local_government_office"),
                                       which(colnames(points)=="food"),
                                       which(colnames(points)=="lodging"),
                                       which(colnames(points)=="museum"),
                                       which(colnames(points)=="movie_theater")
                                       )




weighing_GMM_mat <<- NULL

length_stklist <- length(which(wdcMerged$stocked_out==F))
theta <- c(-4 ,-4, 0, 146.83, 20, 0, rep(0.42, 5),0,1)

active_coef <- c(1:length(theta))[-c(density_add_den_cols[1]-1,density_metro_evening_col-1,density_ridership_col-1)]
active_coef_nointercept <- c(1:length(theta))[-c(density_intercept_col-1,density_add_den_cols[1]-1,density_metro_evening_col-1,
                                                 density_ridership_col-1)]



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

lb = c(-rep(300,2),
       0,rep(0,3),rep(0,5),rep(0,2),
       rep(-20, length_delta), rep(-1.0e19,length_eta))
ub = c(rep(300,2),
       0,rep(1000000,2),0,rep(1000000,5),0,1000000,
       rep(10, length_delta), rep(1.0e19,length_eta))
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
var_covar_theta_save <- eval_theta_variance(theta1,deltain_stin,wdcMerged,points,active_coef)
round(theta[active_coef_nointercept]/sqrt(diag(var_covar_theta_save)),2)
var_covar_theta_2_save <- eval_theta2_variance(theta1,deltain_stin,wdcMerged,points,active_coef)


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

var_covar_theta_2 <- eval_theta2_variance(theta1,deltain_stin,wdcMerged,points,active_coef)
serv_cols <- which(rownames(theta2) %like% "serv_lvl")
round(theta2_save[serv_cols]/sqrt(diag(var_covar_theta_2)[serv_cols]),2)
st_err_theta2 <- sqrt(diag(var_covar_theta_2))
round(correlation_theta_2_estimates <- diag(1/st_err_theta2[serv_cols]) %*% var_covar_theta_2[serv_cols,serv_cols] %*% 
        diag(1/st_err_theta2[serv_cols]),3)

serv_lvl_scale <- scalecols_X[rownames(theta2)[serv_cols]]
##########
theta1_save_t_stat <- rep("",length(theta_save))
theta1_save_t_stat[active_coef_nointercept] <- round(theta_save[active_coef_nointercept]/sqrt(diag(var_covar_theta_save)),2)
theta1_t_stat <- rep("",length(theta))
theta1_t_stat[active_coef_nointercept] <- round(theta[active_coef_nointercept]/sqrt(diag(var_covar_theta)),2)
##########
print("Step 1 output: ")
print("theta1 output: ")
print(paste0(theta_save[c(1:(nonden_ceoflength-1))],collapse = ", " ))
print(paste0(theta_save[-c(1:(nonden_ceoflength-1))],collapse = ", " ))
print("theta1 t-stat:")
print(paste0(theta1_save_t_stat[c(1:(nonden_ceoflength-1))],collapse = ", " ))
print(paste0(theta1_save_t_stat[-c(1:(nonden_ceoflength-1))],collapse = ", " ))
print("theta2 output: ")
print(theta2_save[which(rownames(theta2_save) %like% "serv_lvl")]/serv_lvl_scale)
print("theta2 t-stat:")
print(round(theta2_save[serv_cols]/sqrt(diag(var_covar_theta_2_save)[serv_cols]),2))
print(res_save$objective)
print(res_save$iterations)
print(time_1)

print("Step 2 output: ")
print("theta1 output: ")
print(paste0(theta[c(1:(nonden_ceoflength-1))],collapse = ", " ))
print(paste0(theta[-c(1:(nonden_ceoflength-1))],collapse = ", " ))
print("theta1 t-stat:")
print(paste0(theta1_t_stat[c(1:(nonden_ceoflength-1))],collapse = ", " ))
print(paste0(theta1_t_stat[-c(1:(nonden_ceoflength-1))],collapse = ", " ))
print("theta2 output: ")
print(theta2[which(rownames(theta2) %like% "serv_lvl")]/serv_lvl_scale)
print("theta2 t-stat:")
print(round(theta2[serv_cols]/sqrt(diag(var_covar_theta_2)[serv_cols]),2))
print(res$objective)
print(res$iterations + res_save$iterations)
print(time + time_1)
print(round(correlation_theta_estimates <- diag(1/st_err) %*% var_covar_theta %*% diag(1/st_err),3))
print(round(correlation_theta_2_estimates,3))

sink()
