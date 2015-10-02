source("eval_obj_GMM_model5.R")
source("GetDataGooglePlaces_tiny.R")

print_iter_values <<- 1

#run
source("temp_optimalGMM_model5.R")

sourceCpp("eval_func_2_new_deltaaveraged_stepfunction_steps50_100_200.cpp")
source("eval_func_2_cpp_steps50_100_200.R")
source("eval_func_2_cpp_MPEC_model5_nonsparse_steps50_100_200.R")
source("eval_func_3_cpp_new_2.86_2_steps50_100_200.R")
source("eval_covariates_delta_reg_step_servlvl.R")

weighing_GMM_mat <<- NULL
delta_all <- rnorm(nrow(wdcMerged),-10,1)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]


length_stklist <- length(which(wdcMerged$stocked_out==F))
theta <- c(-4 ,-4,-4,-4, 1, 146.83, 20, 203.69, 0.42, 1,2,3,9.41)
#choose starting values for deltain_stin
colnames_theta1 <<- c("dis coef","rand coef","dis coef","dis coef","dis coef",
                      "density ridership","density metro","density intercept",
                      "density metro evening", "density google places count",
                      "density bus")
nonden_ceofrange <- c(1:5)
nonden_ceoflength <- 5
density_ridership_col <<- 6
density_metro_col <<- 7
density_intercept_col <<- 8
density_metro_evening_col <<- 9
#density_google_places_count <<- 10
density_bus_col <<- 14





density_google_places_theta_cols <- c(10,11,12,13)
density_google_places_points_cols <- c(7,8,9,10)

get_points_density_grad_places_count_col <- NULL

get_points_density <- function(points_in, theta1, tw_in) {
  if(tw_in == 0) {
    density <- ((points_in$type==1)*theta1[density_ridership_col] + 
                  (points_in$type==2)*theta1[density_metro_evening_col]) *points_in$weight   + 
      (points_in$type==1)*theta1[density_intercept_col]    
  } else{
    density <- ((points_in$type==1)*theta1[density_ridership_col] + 
                  (points_in$type==2)*theta1[density_metro_col]) *points_in$weight   + 
      (points_in$type==1)*theta1[density_intercept_col]    
  }
  density <- density + 
    (as.matrix(points_in[,density_google_places_points_cols]) %*% 
    as.matrix(theta1[density_google_places_theta_cols])) +
    (points_in$type==3)*theta1[density_bus_col]
  
  return(density)
}

get_points_density_grad_places_count_col <- function(points_in, tw_in) {
  density <- points_in[,density_google_places_points_cols]
  return(density)
}

theta1 <- get_theta1fromtheta(theta)
tw_in = 1

a <- get_points_density(points_in, theta1, tw_in)

