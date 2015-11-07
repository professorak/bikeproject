source("constants_mw6.R")
sink("r70.1_noRD_bytw_metrodummies_logmetrocorrected_moments_1:20_full_nldis_pr6.txt",split=T)

#least squares implementation 

get_points_density_grad_places_count_col <- NULL
source("eval_obj_GMM_model5.R")
source("GetDataGooglePlaces_20dis_aggmonth_bytw_norealloc_stkoutthresh_five_averaged_pr6.R")
v0_vec <- c(0)
v0_vec_weights <- rep(1,length(v0_vec))
# load(file="v0_vec.RData")
# load(file="v0_vec_weights.RData")
points$lodging <- 0
points$museum <- 0
points$movie_theater <- 0
points$local_government_office <- 0
points$food <- 0
source("debug_functions.R")
points <- assign_log_metro_traffic_corrected(points)
#since we are adding dummy, we can bring lowest to 0 so that coefficient on this is positive.
points_metro_idx <- which(points$type==2)
points$weight[points_metro_idx] <- points$weight[points_metro_idx] - min(points$weight[points_metro_idx])

source("data_estimation_2.6_weather_functions.R")
#regenerate metro moment variables for wdcMerged
wdcMerged_newvars <- get_local_attributes_st_state(wdcMerged[,c("station_id_index","tract","tw","lat","lon","sto_state_local")], points)
list_metro_vars <- c("metro_den_on_1","metro_den_on_2","metro_den_off_1","metro_den_off_2",
                     "metro_den_on_3","metro_den_on_4","metro_den_off_3","metro_den_off_4",
                     "log(metro_den_on_1+1)","log(metro_den_on_2+1)","log(metro_den_on_3+1)","log(metro_den_on_4+1)",
                     "log(metro_den_off_1+1)","log(metro_den_off_2+1)","log(metro_den_off_3+1)","log(metro_den_off_4+1)",
                     "metro_den_on_a","metro_den_on_b","metro_den_on_c",
                     "metro_den_off_a","metro_den_off_b","metro_den_off_c")
wdcMerged[,list_metro_vars] <- wdcMerged_newvars[,list_metro_vars]
rm(wdcMerged_newvars)
print_iter_values <<- 1

#run
source("temp_optimalGMM_model5.R")

sourceCpp("eval_func_2_new_deltaaveraged_nldis_v0weights_set2.cpp")  
#distance step at median quartile 0.300 kms
source("eval_func_2_cpp_steps_v0weights.R")
source("eval_func_2_cpp_MPEC_model5_nonsparse_steps_v0weights.R")
source("eval_func_3_cpp_new_2.86_2_steps_v0weights.R")
source("eval_covariates_delta_reg_densitycols_setb40_wmetrooff_metrodummies.R") 
source("moredensitycols_functions_metrodummies.R")

colnames_theta1 <<- c("dis coef","rand coef","dis coef",
                      "density ridership","density metro","density intercept",
                      "density metro evening","local_government_office", "food",
                      "lodging","museum","movie_theater",
                      "density bus", "density tourist locations","density metro dummy","density metro evening dummy")
nonden_ceofrange <- c(1:3)
nonden_ceoflength <- 3
density_ridership_col <<- 4
density_metro_col <<- 5
density_intercept_col <<- 6
density_metro_evening_col <<- 7
#density_google_places_count <<- 10
#density_bus_col <<- 10
density_add_den_cols <<- c(13,14,15,16) #bus, tourist location, metro dummies


density_google_places_theta_cols <- c(8,9,10,11,12)
density_google_places_points_cols <- c(which(colnames(points)=="local_government_office"),
                                       which(colnames(points)=="food"),
                                       which(colnames(points)=="lodging"),
                                       which(colnames(points)=="museum"),
                                       which(colnames(points)=="movie_theater")
)




weighing_GMM_mat <<- NULL

length_stklist <- length(which(wdcMerged$stocked_out==F))
theta1 <- c(-2.309775e+00,0, -2.500129e+01,0,
            2.101448e-01,2.411672e-01,2.101448e-01,0.000000e+00,0.000000e+00,
            0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.587059e+00,0.1,0.1)

active_coef <- c(1:length(theta1))[-c(2,density_google_places_theta_cols,
                                      density_add_den_cols[1])]
active_coef_nointercept <- c(1:length(theta1))[-c(2,density_intercept_col,
                                                  density_google_places_theta_cols,density_add_den_cols[1])]



#choose starting values for deltain_stin
#delta_all <- rnorm(nrow(wdcMerged),-10,1)
delta_all <- compute_delta_list_cntrt_map_new(theta1,wdcMerged, points)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
# #overriding above null setting as GMM values are quite small for identity weighing_GMM_mat setting
# weighing_GMM_mat <<- eval_GMM_optimal_weighing_matrix(theta1,deltain_stin,wdcMerged,
#                                                       points)
length_eta <- ncol(eval_covariates_delta_reg(deltain_stin,theta1,wdcMerged,points)$Z)

#########
list_covariates <- eval_covariates_delta_reg(deltain_stin,theta1,wdcMerged,points)
X <<- list_covariates$X
Z <<- list_covariates$Z

stocked_list <- which(wdcMerged$stocked_out==FALSE)
Z_sqweighted <- (Z * sqrt(wdcMerged$obs_weight[stocked_list]))/sum(wdcMerged$obs_weight[stocked_list])
print("kappa(t(Z_sqweighted) %*% Z_sqweighted): ")
print(kappa(t(Z_sqweighted) %*% Z_sqweighted))
weighing_GMM_mat <- solve(t(Z_sqweighted) %*% Z_sqweighted)

#########

eta <- rep(0,length_eta)
eta <- c(eval_constraints_moment_conditions(deltain_stin, theta1, eta, wdcMerged, points))
length_theta1 <- length(theta1)
length_delta <- length(deltain_stin)
params <- c(theta1,deltain_stin,eta)
tot_density <<- 10000
target_density_constraint <<- 10000
scale_den <- target_density_constraint/get_total_density(theta1, wdcMerged, points)
theta1[-nonden_ceofrange] <- theta1[-nonden_ceofrange]*scale_den
get_total_density(theta1, wdcMerged, points)


lb = c(c(-100,0,-100),
       0,rep(0,3),rep(0,5),rep(0,4),
       rep(-20, length_delta), rep(-1.0e19,length_eta))
ub = c(c(50,0,50),
       1000000,rep(1000000,2),1000000,rep(0,5),0,1000000,1000000,1000000,
       rep(10, length_delta), rep(1.0e19,length_eta))
constraint_lb <- c(rep(0, length_delta+length_eta), target_density_constraint) 
constraint_ub <- c(rep(0, length_delta+length_eta), target_density_constraint)


eval_jac_g_structure_val <- eval_jac_g_structure(params, wdcMerged, points, length_theta1, length(deltain_stin), length_eta)
eval_hess_struct <- eval_hessian_structure(params,wdcMerged, points, length_theta1, length(deltain_stin), length_eta)
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
               length_theta=length_theta1,
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
theta <- params[c(1:length_theta1)]
theta_save <- theta
deltain_stin <- params[length_theta1+c(1:length_delta)]
eta <- params[length_theta1+length_delta+c(1:length_eta)]
theta1 <- theta_save  
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
eval_jac_g_structure_val <- eval_jac_g_structure(params, wdcMerged, points, length_theta1, length(deltain_stin), length_eta)
eval_hess_struct <- eval_hessian_structure(params,wdcMerged, points, length_theta1, length(deltain_stin), length_eta)
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
               length_theta=length_theta1,
               length_delta=length_delta,
               length_eta=length_eta,
               opts=opts
) 
time <- proc.time() - time
print("time: ")
print(time)

params <- res$solution 


theta <- params[c(1:length_theta1)]
deltain_stin <- params[length_theta1+c(1:length_delta)]
theta1 <- theta  
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

source("debug_functions.R")
demand_division_df <- c() 
for(i in c(4:length(theta1))) {
  theta_division <- rep(0,length(theta1))
  theta_division[c(1:3)] <-  theta1[c(1:3)]
  theta_division[i] <- theta1[i]
  demand_i <- eval_lambda_full(deltain_stin, theta_division,wdcMerged, points)
  demand_division_df <- cbind(demand_division_df, demand_i)
}
colnames(demand_division_df) <- colnames_theta1[4:length(colnames_theta1)]
summary(rowSums(demand_division_df)-wdcMerged$out_dem_mean)
print(t(t(colSums(demand_division_df)))/sum(demand_division_df)*100)

##########
theta1_save_t_stat <- rep("",length(theta_save))
theta1_save_t_stat[active_coef_nointercept] <- round(theta_save[active_coef_nointercept]/sqrt(diag(var_covar_theta_save)),2)
theta1_t_stat <- rep("",length_theta1)
theta1_t_stat[active_coef_nointercept] <- round(theta[active_coef_nointercept]/sqrt(diag(var_covar_theta)),2)
##########
print("Step 1 output: ")
print("theta1 output: ")
print(paste0(theta_save[c(1:(nonden_ceoflength))],collapse = ", " ))
print(paste0(theta_save[-c(1:(nonden_ceoflength))],collapse = ", " ))
print("theta1 t-stat:")
print(paste0(theta1_save_t_stat[c(1:(nonden_ceoflength))],collapse = ", " ))
print(paste0(theta1_save_t_stat[-c(1:(nonden_ceoflength))],collapse = ", " ))
print("theta2 output: ")
print(theta2_save[which(rownames(theta2_save) %like% "serv_lvl")]/serv_lvl_scale)
print("theta2 t-stat:")
print(round(theta2_save[serv_cols]/sqrt(diag(var_covar_theta_2_save)[serv_cols]),2))
print(res_save$objective)
print(res_save$iterations)
print(time_1)

print("Step 2 output: ")
print("theta1 output: ")
print(paste0(theta[c(1:(nonden_ceoflength))],collapse = ", " ))
print(paste0(theta[-c(1:(nonden_ceoflength))],collapse = ", " ))
print("theta1 t-stat:")
print(paste0(theta1_t_stat[c(1:(nonden_ceoflength))],collapse = ", " ))
print(paste0(theta1_t_stat[-c(1:(nonden_ceoflength))],collapse = ", " ))
print("theta2 output: ")
print(theta2[which(rownames(theta2) %like% "serv_lvl")]/serv_lvl_scale)
print("theta2 t-stat:")
print(round(theta2[serv_cols]/sqrt(diag(var_covar_theta_2)[serv_cols]),2))
print(res$objective)
print(res$iterations + res_save$iterations)
print(time + time_1)
print(round(correlation_theta_estimates <- diag(1/st_err) %*% var_covar_theta %*% diag(1/st_err),3))
print(round(correlation_theta_2_estimates,3))

thetaden_s1 <- theta_save[c(1:(nonden_ceoflength))]
thetanonden_s1 <- theta_save[-c(1:(nonden_ceoflength))]  
thetatstat_den_s1 <- theta1_save_t_stat[c(1:(nonden_ceoflength))]
thetatstat_nonden_s1 <- theta1_save_t_stat[-c(1:(nonden_ceoflength))]

theta2coef_s1 <- theta2_save[which(rownames(theta2_save) %like% "serv_lvl")]/serv_lvl_scale
theta2coeftstat_s1 <- round(theta2_save[serv_cols]/sqrt(diag(var_covar_theta_2_save)[serv_cols]),2)
obj_s1 <- res_save$objective
itr_s1 <- res_save$iterations
runtime_s1 <- time_1[1]

thetaden <- theta[c(1:(nonden_ceoflength))]
thetanonden <- theta[-c(1:(nonden_ceoflength))]  
thetatstat_den <- theta1_t_stat[c(1:(nonden_ceoflength))]
thetatstat_nonden <- theta1_t_stat[-c(1:(nonden_ceoflength))]

theta2coef <- theta2[which(rownames(theta2) %like% "serv_lvl")]/serv_lvl_scale
theta2coeftstat <- round(theta2[serv_cols]/sqrt(diag(var_covar_theta_2)[serv_cols]),2)
obj <- res$objective
itr <- res$iterations + res_save$iterations
runtime <- time[1] + time_1[1]

source("print_output_functions.R")
printoutput_theta_csvformattted (thetaden_s1, thetanonden_s1, thetatstat_den_s1, thetatstat_nonden_s1, theta2coef_s1, 
                                 theta2coeftstat_s1, obj_s1, itr_s1, runtime_s1, thetaden, thetanonden, thetatstat_den, thetatstat_nonden, theta2coef,
                                 theta2coeftstat, obj, itr, runtime)

sink()
