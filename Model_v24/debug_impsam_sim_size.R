source("constants_mw6.R")
source("debug_functions.R")
#model with averaged station demand, no service level. Trying to get distance
#coefficient right.

get_points_density_grad_places_count_col <- NULL
source("eval_obj_GMM_model5.R")
source("GetDataGooglePlaces_aggmonth_averaged_pr6_tiny.R")

print_iter_values <<- 1

#run
source("temp_optimalGMM_model5.R")

sourceCpp("eval_func_2_new_deltaaveraged_stepfunction_v0weights_set2.cpp")
#distance step at median quartile 0.300 kms
source("eval_func_2_cpp_steps_v0weights.R")
source("eval_func_2_cpp_MPEC_model5_nonsparse_steps_v0weights.R")
source("eval_func_3_cpp_new_2.86_2_steps_v0weights.R")
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
theta1 <- c(-15.0363528438044, 4, -11.69613153711,
            c(0, 3.56400131121065, 0.0151809839374334, 0, 0.00822638613922138, 0.000751098339022488,
              9.47693298685854e-10, 0.00477686078154983, 8.46891338865097e-08, 0, 0.981490763648657)*10)

theta1_test <- theta1
theta1_test[c(1:3)] <- c(-12, 3, -18)

get_total_density(theta1, wdcMerged, points )
sum(wdcMerged$out_dem_mean)/get_total_density(theta1, wdcMerged, points )*100


#for lambda's
source("Importance_sampling_functions.R")

set.seed(2)
v0_vec <<- rnorm(100)
v0_vec_weights <<- rep(1,length(v0_vec))
delta_all <- compute_delta_list_cntrt_map_new(theta1,wdcMerged, points)

delta_all_test <- compute_delta_list_cntrt_map_new(theta1_test,wdcMerged, points)

lambda_df <- c()
size_vec <- c(1,2,4,8,20*2^c(0:7))
for(i in 1:length(size_vec)) {
  print(i)
  v0_vec_df <- generate_v0_importance_sampling(size_vec[i], delta_all, theta1, wdcMerged, points)
  v0_vec <<- v0_vec_df[,1]
  v0_vec_weights <<- v0_vec_df[,2]
  
  lambda_t <- eval_lambda_full(delta_all_test, theta1_test, wdcMerged, points)
  lambda_df <- cbind(lambda_df, lambda_t)  
  
}
colnames(lambda_df) <- size_vec

summary_df <- c()
for(i in 1:(ncol(lambda_df)-1)) {
  j <- i+1
  summary_df <- rbind(summary_df,
                      c(size_vec[i],as.numeric(summary(lambda_df[,i]-lambda_df[,j]-mean(lambda_df[,i]-lambda_df[,j])))))
}
colnames(summary_df) <- c("simsize",
                          "Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")

View(summary_df)

lambda_df_impsap <- lambda_df
summary_df_impsap <- summary_df


#############
#pseduo random
lambda_df <- c()

for(i in 1:length(size_vec)) {
  print(i)
  v0_vec <<- rnorm(size_vec[i])
  v0_vec_weights <<- rep(1, length(v0_vec))
  
  lambda_t <- eval_lambda_full(delta_all_test, theta1_test, wdcMerged, points)
  lambda_df <- cbind(lambda_df, lambda_t)  
  
}
colnames(lambda_df) <- size_vec

summary_df <- c()
for(i in 1:(ncol(lambda_df)-1)) {
  j <- i+1
  summary_df <- rbind(summary_df,
                      c(size_vec[i],as.numeric(summary(lambda_df[,i]-lambda_df[,j]-mean(lambda_df[,i]-lambda_df[,j])))))
}
colnames(summary_df) <- c("simsize",
                          "Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")

lambda_df_pseduor <- lambda_df
summary_df_pseduor <- summary_df



