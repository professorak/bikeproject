source("constants_mw6.R")
source("debug_functions.R")
#model with averaged station demand, no service level. Trying to get distance
#coefficient right.

get_points_density_grad_places_count_col <- NULL
source("eval_obj_GMM_model5.R")
source("GetDataGooglePlaces_aggmonth_averaged_pr6.R")
v0_vec <- generate_v0(100)

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
theta1 <- c(-15.7231085568937, 4, -6.36880404225048,
  c(0, 336.073466350849, 1.53276238137761, 0, 0.92488430899128, 0.0694077909057493, 1.08221956077344e-06, 
    0.462573747981363, 7.00834802982461e-06, 0, 105.623382345265)*1000/1e5)

active_coef <- c(1:length(theta1))[-c(density_add_den_cols[1],density_metro_evening_col,density_ridership_col)]
active_coef_nointercept <- c(1:length(theta1))[-c(density_intercept_col,density_add_den_cols[1],density_metro_evening_col,
                                                 density_ridership_col)]


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
length_theta <- length(theta1)
length_delta <- length(deltain_stin)
params <- c(theta1,deltain_stin,eta)


#eval obj
obj_1 <- eval_obj_GMM_model5_obj(params, wdcMerged, points, length_theta, length(deltain_stin), length_eta)

xi <- eval_error_xi_model5 (deltain_stin, 
                                  theta1, wdcMerged,points)$xi

#divide the demand coming at station from different density sources 
demand_division_df <- c() 
for(i in c(3:length(theta))) {
  theta_division <- rep(0,length(theta))
  theta_division[c(1:2)] <-  theta[c(1:2)]
  theta_division[i] <- theta[i]
  demand_i <- eval_lambda_full(deltain_stin, get_theta1fromtheta(theta_division),wdcMerged, points)
  demand_division_df <- cbind(demand_division_df, demand_i)
}
colnames(demand_division_df) <- colnames_theta1[4:length(colnames_theta1)]
summary(rowSums(demand_division_df)-wdcMerged$out_dem_mean)
View(t(t(colSums(demand_division_df)))/sum(demand_division_df)*100)

#finding closest type 1 points
points_intercept <- subset(points, type==1)
dis_st_closest_point <- rep(0,nrow(wdcMerged))
for(i in 1:nrow(wdcMerged)) {
  dis_v <- latlondistance(points_intercept$lat, points_intercept$lon, wdcMerged$lat[i],wdcMerged$lon[i])
  dis_st_closest_point[i] <- min(dis_v)
}


cor(dis_st_closest_point, xi-log(wdcMerged$out_dem_mean))

plot(sort(dis_st_closest_point),type="l")
par(new=T)
plot((deltain_stin-log(wdcMerged$out_dem_mean))[order(dis_st_closest_point)],type="l")

#########
#comparing with the case when theta discoef is -2,-2

# theta1_case2 <- c(-2,-2,
#   0, 340.198595762635, 1.55116482669844, 0, 0.49122522409031, 0.0606975905045715, 6.5052950946276e-08, 0.509453436659953, 4.4183408417253e-07, 0, 109.531656793844)
# theta1_case2 <- c(-2, -2,
#   0, 454.365056267057, 0.642354054985412, 0, 5.37565966472191, 0.264658314033594, 0.504334670936989, 2.21883319392775, 0.193297030056388, 0, 437.589412742883)
theta1_case2 <-  c(-15.7231085568937, 0, -6.36880404225048,
                    c(0, 336.073466350849, 1.53276238137761, 0, 0.92488430899128, 0.0694077909057493, 1.08221956077344e-06, 
                    0.462573747981363, 7.00834802982461e-06, 0, 105.623382345265)*1000/1e5)

get_total_density(theta1_case2,wdcMerged,points)
delta_all_case2 <- compute_delta_list_cntrt_map_new(theta1_case2,wdcMerged, points)
deltain_stin_case2 <- delta_all_case2[which(wdcMerged$stocked_out==F)]

eta_case2 <- rep(0,length_eta)
eta_case2 <- c(eval_constraints_moment_conditions(deltain_stin_case2, 
            theta1_case2, eta_case2, wdcMerged, points))
params_case2 <- c(theta1_case2,deltain_stin_case2,eta_case2)


obj_case2 <- eval_obj_GMM_model5_obj(params_case2, wdcMerged, points, length(theta1), length(deltain_stin), length_eta)

xi_case2 <- eval_error_xi_model5 (deltain_stin_case2, 
              theta1_case2, wdcMerged,points)$xi

#divide the demand coming at station from different density sources 
demand_division_case2_df <- c() 
for(i in c(3:length(theta1_case2))) {
  theta_division_case2 <- rep(0,length(theta1_case2))
  theta_division_case2[c(1:2)] <-  theta1_case2[c(1:2)]
  theta_division_case2[i] <- theta1_case2[i]
  demand_case2_i <- eval_lambda_full(deltain_stin_case2, get_theta1fromtheta(theta_division_case2),wdcMerged, points)
  demand_division_case2_df <- cbind(demand_division_case2_df, demand_case2_i)
}
colnames(demand_division_case2_df) <- colnames_theta1[4:length(colnames_theta1)]
summary(rowSums(demand_division_case2_df)-wdcMerged$out_dem_mean)

View(t(t(colSums(demand_division_case2_df)))/sum(demand_division_case2_df)*100)

print(theta1_case2,digits=2)
t(t(colSums(demand_division_case2_df)))
cbind(colnames(Z),eta_case2)
t(eta_case2) %*% weighing_GMM_mat %*% eta_case2
######################################################################
######################################################################
#comparing
eta_df <- cbind(eta, eta_case2, diag(weighing_GMM_mat))
rownames(eta_df) <- colnames(Z)
#are there outliers in xi values:
plot(xi, type="l")
plot(deltain_stin, type="l")
plot(xi_case2, type="l")
plot(deltain_stin_case2, type="l")
plot(xi_case2-xi, type="l")

outliers_1 <- wdcMerged[which(xi-xi_case2<=-0.3),
          c("lat","lon","station_id","station_id_index","tract","out_dem_mean","catchment_area_step2","metro_den_on_1","metro_den_on_2",
            "googleplaces_den_1","googleplaces_food_den_1","googleplaces_grocery_den_1")]

summary(wdcMerged$out_dem_sum)
summary(wdcMerged$metro_den_on_1)
summary(wdcMerged$googleplaces_den_1)
summary(wdcMerged$googleplaces_food_den_1)
summary(wdcMerged$googleplaces_grocery_den_1)

#find the quantile of density values of above outliers:
which(sort(wdcMerged$metro_den_on_1)>max(outliers_1$metro_den_on_1))[1]/nrow(wdcMerged)
which(sort(wdcMerged$googleplaces_food_den_1)>max(outliers_1$googleplaces_food_den_1))[1]/nrow(wdcMerged)
which(sort(wdcMerged$googleplaces_grocery_den_1)>max(outliers_1$googleplaces_grocery_den_1))[1]/nrow(wdcMerged)


#looking at values of xi which are high, these will generally hint what density variables are missing
top_xi_idx <- order(xi, decreasing = T)[1:10]
xi[top_xi_idx]
outliers_2 <- wdcMerged[top_xi_idx,
                        c("lat","lon","station_id","station_id_index","tract","out_dem_mean","catchment_area_step2","metro_den_on_1","metro_den_on_2",
                          "googleplaces_den_1","googleplaces_den_2","googleplaces_food_den_1","googleplaces_food_den_2","googleplaces_grocery_den_1",
                          "googleplaces_grocery_den_2")]
#find the quantile of density values of above outliers:
#max
which(sort(wdcMerged$out_dem_mean)>=max(outliers_2$out_dem_mean))[1]/nrow(wdcMerged)
which(sort(wdcMerged$metro_den_on_1)>=max(outliers_2$metro_den_on_1))[1]/nrow(wdcMerged)
which(sort(wdcMerged$googleplaces_food_den_1)>max(outliers_2$googleplaces_food_den_1))[1]/nrow(wdcMerged)
which(sort(wdcMerged$googleplaces_grocery_den_1)>max(outliers_2$googleplaces_grocery_den_1))[1]/nrow(wdcMerged)
#mean
which(sort(wdcMerged$out_dem_mean)>=mean(outliers_2$out_dem_mean))[1]/nrow(wdcMerged)
which(sort(wdcMerged$metro_den_on_1)>=mean(outliers_2$metro_den_on_1))[1]/nrow(wdcMerged)
which(sort(wdcMerged$googleplaces_food_den_1)>mean(outliers_2$googleplaces_food_den_1))[1]/nrow(wdcMerged)
which(sort(wdcMerged$googleplaces_grocery_den_1)>mean(outliers_2$googleplaces_grocery_den_1))[1]/nrow(wdcMerged)
#index i
for(i in 1:nrow(outliers_2)) {
  print(outliers_2$station_id_index[i])
  vals <- c()
  vals <- c(vals,which(sort(wdcMerged$out_dem_mean)>=outliers_2$out_dem_mean[i])[1]/nrow(wdcMerged))
  vals <- c(vals,which(sort(wdcMerged$metro_den_on_1)>=outliers_2$metro_den_on_1[i])[1]/nrow(wdcMerged))
  vals <- c(vals,which(sort(wdcMerged$googleplaces_food_den_1)>outliers_2$googleplaces_food_den_1[i])[1]/nrow(wdcMerged))
  vals <- c(vals,which(sort(wdcMerged$googleplaces_food_den_2)>outliers_2$googleplaces_food_den_2[i])[1]/nrow(wdcMerged))
  vals <- c(vals,which(sort(wdcMerged$googleplaces_grocery_den_1)>outliers_2$googleplaces_grocery_den_1[i])[1]/nrow(wdcMerged))
  vals <- c(vals,which(sort(wdcMerged$googleplaces_den_1)>outliers_2$googleplaces_den_1[i])[1]/nrow(wdcMerged))
  vals <- c(vals,which(sort(wdcMerged$googleplaces_den_2)>outliers_2$googleplaces_den_2[i])[1]/nrow(wdcMerged))
  print(cbind(c("out_dem_mean","metro_den_on_1","googleplaces_food_den_1","googleplaces_food_den_2",
                "googleplaces_grocery_den_1","googleplaces_den_1","googleplaces_den_2"),
              round(vals,3)))

}

####################################
#looking at values of xi which are high and low, these will generally hint what density variables are missing
top_xi_idx <- order(abs(xi), decreasing = T)[1:10]
outliers_2 <- wdcMerged[top_xi_idx,
                        c("lat","lon","station_id","station_id_index","tract","out_dem_mean","catchment_area_step2","metro_den_on_1","metro_den_on_2",
                          "googleplaces_den_1","googleplaces_den_2","googleplaces_food_den_1","googleplaces_food_den_2","googleplaces_grocery_den_1",
                          "googleplaces_grocery_den_2")]
outliers_2 <- cbind(outliers_2, xi[top_xi_idx],deltain_stin[top_xi_idx])
View(t(outliers_2))
points_type1 <- subset(points, type==1)
#index i
for(i in 1:nrow(outliers_2)) {
  print(outliers_2$station_id_index[i])
  vals <- c()
  vals <- c(vals,which(sort(wdcMerged$out_dem_mean)>=outliers_2$out_dem_mean[i])[1]/nrow(wdcMerged))
  vals <- c(vals,which(sort(wdcMerged$metro_den_on_1)>=outliers_2$metro_den_on_1[i])[1]/nrow(wdcMerged))
  vals <- c(vals,which(sort(wdcMerged$metro_den_on_2)>=outliers_2$metro_den_on_2[i])[1]/nrow(wdcMerged))
  vals <- c(vals,which(sort(wdcMerged$googleplaces_food_den_1)>=outliers_2$googleplaces_food_den_1[i])[1]/nrow(wdcMerged))
  vals <- c(vals,which(sort(wdcMerged$googleplaces_food_den_2)>=outliers_2$googleplaces_food_den_2[i])[1]/nrow(wdcMerged))
  vals <- c(vals,which(sort(wdcMerged$googleplaces_grocery_den_1)>=outliers_2$googleplaces_grocery_den_1[i])[1]/nrow(wdcMerged))
  vals <- c(vals,min(latlondistance(points_type1$lat,points_type1$lon, outliers_2$lat[i],outliers_2$lon[i])))
#   vals <- c(vals,which(sort(wdcMerged$googleplaces_den_1)>=outliers_2$googleplaces_den_1[i])[1]/nrow(wdcMerged))
#   vals <- c(vals,which(sort(wdcMerged$googleplaces_den_2)>=outliers_2$googleplaces_den_2[i])[1]/nrow(wdcMerged))
  print(cbind(c("out_dem_mean","metro_den_on_1","metro_den_on_2","googleplaces_food_den_1","googleplaces_food_den_2",
                "googleplaces_grocery_den_1","min_dis_point"
#                ,"googleplaces_den_1","googleplaces_den_2"
                ),
              round(vals,3)))
  
}

#looking at values of xi_case2 which are high, these will generally hint what density variables are missing
top_xi_idx <- order(xi_case2, decreasing = T)[1:10]
xi_case2[top_xi_idx]
outliers_2 <- wdcMerged[top_xi_idx,
                        c("lat","lon","station_id","station_id_index","tract","out_dem_mean","catchment_area_step2","metro_den_on_1","metro_den_on_2",
                          "googleplaces_den_1","googleplaces_den_2","googleplaces_food_den_1","googleplaces_food_den_2","googleplaces_grocery_den_1",
                          "googleplaces_grocery_den_2")]
outliers_2 <- cbind(outliers_2, xi_case2[top_xi_idx],deltain_stin_case2[top_xi_idx])
View(t(outliers_2))
#find the quantile of density values of above outliers:
#max
which(sort(wdcMerged$out_dem_mean)>=max(outliers_2$out_dem_mean))[1]/nrow(wdcMerged)
which(sort(wdcMerged$metro_den_on_1)>=max(outliers_2$metro_den_on_1))[1]/nrow(wdcMerged)
which(sort(wdcMerged$googleplaces_food_den_1)>max(outliers_2$googleplaces_food_den_1))[1]/nrow(wdcMerged)
which(sort(wdcMerged$googleplaces_grocery_den_1)>max(outliers_2$googleplaces_grocery_den_1))[1]/nrow(wdcMerged)
#mean
which(sort(wdcMerged$out_dem_mean)>=mean(outliers_2$out_dem_mean))[1]/nrow(wdcMerged)
which(sort(wdcMerged$metro_den_on_1)>=mean(outliers_2$metro_den_on_1))[1]/nrow(wdcMerged)
which(sort(wdcMerged$googleplaces_food_den_1)>mean(outliers_2$googleplaces_food_den_1))[1]/nrow(wdcMerged)
which(sort(wdcMerged$googleplaces_grocery_den_1)>mean(outliers_2$googleplaces_grocery_den_1))[1]/nrow(wdcMerged)
#index i
for(i in 1:nrow(outliers_2)) {
  print(outliers_2$station_id_index[i])
  vals <- c()
  vals <- c(vals,which(sort(wdcMerged$out_dem_mean)>=outliers_2$out_dem_mean[i])[1]/nrow(wdcMerged))
  vals <- c(vals,which(sort(wdcMerged$metro_den_on_1)>=outliers_2$metro_den_on_1[i])[1]/nrow(wdcMerged))
  vals <- c(vals,which(sort(wdcMerged$googleplaces_food_den_1)>outliers_2$googleplaces_food_den_1[i])[1]/nrow(wdcMerged))
  vals <- c(vals,which(sort(wdcMerged$googleplaces_food_den_2)>outliers_2$googleplaces_food_den_2[i])[1]/nrow(wdcMerged))
  vals <- c(vals,which(sort(wdcMerged$googleplaces_grocery_den_1)>outliers_2$googleplaces_grocery_den_1[i])[1]/nrow(wdcMerged))
  vals <- c(vals,which(sort(wdcMerged$googleplaces_den_1)>outliers_2$googleplaces_den_1[i])[1]/nrow(wdcMerged))
  vals <- c(vals,which(sort(wdcMerged$googleplaces_den_2)>outliers_2$googleplaces_den_2[i])[1]/nrow(wdcMerged))
  print(cbind(c("out_dem_mean","metro_den_on_1","googleplaces_food_den_1","googleplaces_food_den_2",
                "googleplaces_grocery_den_1","googleplaces_den_1","googleplaces_den_2"),
              round(vals,3)))
  
}


