source("constants_mw6.R")

#least squares implementation #no instrument for service level 

get_points_density_grad_places_count_col <- NULL
source("eval_obj_GMM_model5.R")
#source("GetDataGooglePlaces_20dis_aggmonth_hyperlocalstate_norealloc_stkoutthresh_five_averaged_pr6.R")
source("GetDataGooglePlaces_20dis_aggmonth_bytw_norealloc_stkoutthresh_five_averaged_pr3.R")

wdcMerged <- subset(wdcMerged, index<=2)

v0_vec <- c(0)
v0_vec_weights <- rep(1,length(v0_vec))
# load(file="v0_vec.RData")
# load(file="v0_vec_weights.RData")
points$lodging <- 0
points$museum <- 0
points$movie_theater <- 0
points$local_government_office <- 0
points$food <- 0
#log of absolute metro traffic
points_metro_idx <- which(points$type==2)
points$weight[points_metro_idx] <- log(points$weight[points_metro_idx]*22468468)
source("data_estimation_2.6_weather_functions.R")

#regenerate metro moment variables for wdcMerged
# wdcMerged_newvars <- get_local_attributes_st_state(wdcMerged[,c("station_id_index","tract","tw","lat","lon","sto_state_local")], points)
# list_metro_vars <- c("metro_den_on_1","metro_den_on_2","metro_den_off_1","metro_den_off_2",
#                      "metro_den_on_3","metro_den_on_4","metro_den_off_3","metro_den_off_4",
#                      "log(metro_den_on_1+1)","log(metro_den_on_2+1)","log(metro_den_on_3+1)","log(metro_den_on_4+1)",
#                      "log(metro_den_off_1+1)","log(metro_den_off_2+1)","log(metro_den_off_3+1)","log(metro_den_off_4+1)",
#                      "metro_den_on_a","metro_den_on_b","metro_den_on_c",
#                      "metro_den_off_a","metro_den_off_b","metro_den_off_c")
# wdcMerged[,list_metro_vars] <- wdcMerged_newvars[,list_metro_vars]
# rm(wdcMerged_newvars)

#since we are adding dummy, we can bring lowest to 0 so that coefficient on this is positive.
points_metro_idx <- which(points$type==2)
points$weight[points_metro_idx] <- points$weight[points_metro_idx] - min(points$weight[points_metro_idx])

print_iter_values <<- 1

#run
source("temp_optimalGMM_model5.R")

sourceCpp("eval_func_2_hyperlocal_nldis_v0weights_set2.cpp")
#distance step at median quartile 0.300 kms
source("eval_func_2_cpp_steps_v0weights.R")
source("eval_func_2_cpp_MPEC_model5_nonsparse_steps_v0weights.R")
source("eval_func_3_cpp_new_2.86_2_steps_v0weights.R")
source("eval_covariates_delta_reg_densitycols_setb48_wservlvl_metrooff_metrodummies.R") 
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



tw_groupin <- wdcMerged$tw_group[1]

set.seed(1)
delta_all <- rnorm(nrow(wdcMerged),-2,1)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
deltain_tw <- deltain_stin[which(wdcMerged$tw_group==tw_groupin)]
wdcMergedday <- wdcMerged[which(wdcMerged$tw_group==tw_groupin),]


length_eta <- ncol(eval_covariates_delta_reg(deltain_stin,theta1,wdcMerged,points)$Z)

list_covariates <- eval_covariates_delta_reg(deltain_stin,theta1,wdcMerged,points)
X <<- list_covariates$X
Z <<- list_covariates$Z

stocked_list <- which(wdcMerged$stocked_out==FALSE)
Z_sqweighted <- (Z * sqrt(wdcMerged$obs_weight[stocked_list]))/sum(wdcMerged$obs_weight[stocked_list])
print("kappa(t(Z_sqweighted) %*% Z_sqweighted): ")
print(kappa(t(Z_sqweighted) %*% Z_sqweighted))
weighing_GMM_mat <- solve(t(Z_sqweighted) %*% Z_sqweighted)

eta <- rep(0,length_eta)
eta <- c(eval_constraints_moment_conditions(deltain_stin, theta1, eta, wdcMerged, points))
length_theta1 <- length(theta1)
length_delta <- length(deltain_stin)
params <- c(theta1,deltain_stin,eta)
tot_density <<- 10000
target_density_constraint <<- 10000

#sourceCpp("eval_func_2_new_deltaaveraged_nldis_v0weights_set2.cpp")

sourceCpp("eval_func_2_hyperlocal_nldis_v0weights_set2.cpp")
lambda_t_list <- eval_lambda_delta_list_new (deltain_tw, theta1, wdcMergedday, points, tw_groupin)
summary(lambda_t_list[[1]])
summary(c(lambda_t_list[[2]]))

eval_jac_g_structure_val <- eval_jac_g_structure(params, wdcMerged, points, length_theta1, length(deltain_stin), length_eta)
length(eval_jac_g_structure_val)
length(unlist(eval_jac_g_structure_val))

eval_hess_struct <- eval_hessian_structure(params,wdcMerged, points, length_theta1, length(deltain_stin), length_eta)
length(eval_hess_struct)
length(unlist(eval_hess_struct))

v_eval_g <- eval_g(params, wdcMerged, points, length_theta1, length(deltain_stin), length_eta)
summary(v_eval_g)
v_eval_jac_g <- eval_jac_g(params, wdcMerged, points, length_theta1, length(deltain_stin), length_eta)
v_eval_jac_g[1:10]
summary(v_eval_jac_g)
length(v_eval_jac_g)

v_eval_obj_GMM_model5_obj <- eval_obj_GMM_model5_obj(params, wdcMerged, points, length_theta1, length(deltain_stin), length_eta)
v_eval_obj_GMM_model5_obj

v_eval_obj_GMM_model5_grad <- eval_obj_GMM_model5_grad(params, wdcMerged, points, length_theta1, length(deltain_stin), length_eta)
summary(v_eval_obj_GMM_model5_grad)

set.seed(1)
obj_factor <- rnorm(1,0,1)  
hessian_lambda <- rnorm(length(v_eval_g),0,1)  
v_eval_hessian <- eval_hessian(params, wdcMerged, points, length_theta1, length(deltain_stin), length_eta,
                               obj_factor,hessian_lambda)
v_eval_hessian[1:10]
summary(v_eval_hessian)
length(v_eval_hessian)









# b <- eval_lambda_delta_list_new (deltain_tw, theta1, wdcMergedday, points, tw_groupin)
# 
# summary(a[[1]])
# summary(b[[1]])
# 
# summary(c(a[[2]]))
# summary(c(b[[2]]))
# 
# identical(c(a[[2]]), c(b[[2]]))


