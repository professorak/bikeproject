#This file generates data for observations weights=1 and assumed values of deltas and 
#then estimates them

source("data_estimation_2.6_weather_saved.R")
#source("data_estimation_2.6_weather_small_saved.R")
source("eval_func_2_cpp_MPEC.R")
source("temp_read_google_places_data.R")
points_save <- points
wdcMerged_store <- wdcMerged
user_serv_lvl_store <- user_serv_lvl

wdcMerged <- wdcMerged_store
user_serv_lvl <- user_serv_lvl_store
current_serv_lvl <- user_serv_lvl_store
points <- points_save

#merge points w google places data
#merge points with places data
points <- points[order(points$type, points$lat, points$lon),]
places_data <- read_googleplaces_data()
places_cols_select <- c("lat","lon","cafe", "grocery_or_supermarket", "local_government_office")
places_data <- places_data[,places_cols_select]
places_data <- places_data[order(places_data$lat, places_data$lon),]
if(!identical(round(points$lat[c(1:nrow(places_data))],4), round(places_data$lat,4))) stop("points lat dont match")
if(!identical(round(points$lon[c(1:nrow(places_data))],4), round(places_data$lon,4))) stop("points lon dont match")
places_data_temp <- places_data
places_data_temp$lat <- NULL
places_data_temp$lon <- NULL
places_colnames <- colnames(places_data_temp)
places_data_temp_full <- as.data.frame(matrix(0,nrow=nrow(points),ncol=ncol(places_data_temp)))
colnames(places_data_temp_full) <- places_colnames
places_data_temp_full[c(1:nrow(places_data_temp)),] <- places_data_temp
points <- cbind(points, places_data_temp_full)

#make the weights to 1
#put a value of delta
#compute lambda
colnames_theta1 <<- c("dis coef","rand coef","density ridership","density metro","density intercept",
                      "density metro evening", "density cafe", "density grocery_or_supermarket", "density local_government_office")
density_ridership_col <<- 3
density_metro_col <<- 4
density_intercept_col <<- 5
density_metro_evening_col <<- 6
density_cafe <<- 7
density_grocery <<- 8
density_gov_office <<- 9
density_bus_col <<- 10

# #keeping "6_0" "6_3" and "7_0" "7_3"
# time_windows <- c("6_0","6_5","7_0","7_5")
# tws <- c("0","5")
# #time_windows <- c("6_3")
# wdcMerged <- subset(wdcMerged, tw_group_fac %in% time_windows)
# wdcMerged <- droplevels(wdcMerged)
# levels(wdcMerged$tw_group_fac)
# user_serv_lvl <- subset(user_serv_lvl, tw %in% tws)
# user_serv_lvl <- droplevels(user_serv_lvl)
# unique(user_serv_lvl$tw)
# 
# wdcMerged <- droplevels(wdcMerged)
# user_serv_lvl$st_tw_index <- as.numeric(user_serv_lvl$st_tw)
# wdcMerged$st_tw_index <- as.numeric(wdcMerged$st_tw)
rm(current_serv_lvl)

# #making 0 demand observations insignificant
# wdcMerged$obs_weight[which(wdcMerged$stocked_out==F & wdcMerged$out_dem_sum<=0.1)] <- 0.1

# wdcMerged$out_dem_sum <- wdcMerged$out_dem_sum/wdcMerged$obs_weight
# wdcMerged$obs_weight <- 1
#removing really low demand
wdcMerged$out_dem_sum[which(wdcMerged$out_dem_sum<=0.01 & 
                              wdcMerged$stocked_out==FALSE)] <- 0.01





######
#test
#test gradient manually for fewer dimensions
theta <- c(-10.7374110745,0.8860478015,403.4177015742,2.8258972200, 200,10,10,10, 1)

set.seed(34675)
deltain_stin <- rnorm(length(which(!wdcMerged$stocked_out)),-3, 0.5)
params <- c(theta,deltain_stin)  

a1_eval_constraints <- eval_g(params,wdcMerged, points, length(theta))
get_total_density(params[1:length(theta)], wdcMerged, points)
a1_eval_jac_constraints <- eval_jac_g(params,wdcMerged, points, length(theta))
eval_jac_g_structure_val <- eval_jac_g_structure(params, wdcMerged, points, length(theta))
if(length(a1_eval_jac_constraints)!=length(unlist(eval_jac_g_structure_val))) stop("gradient lengths dont match a1_eval_jac_constraints")
a1_eval_jac_constraints <- getfullfromsparsematrix (eval_jac_g_structure_val, a1_eval_jac_constraints) 


idx <- 9
diff <- 0.0001

params_2 <- params
params_2[idx] <- params_2[idx] + diff  

a2_eval_constraints <- eval_g(params_2,wdcMerged, points, length(theta))
get_total_density(params_2[1:length(theta)], wdcMerged, points)
a2_eval_jac_constraints <- eval_jac_g(params_2,wdcMerged, points, length(theta))
a2_eval_jac_constraints <- getfullfromsparsematrix (eval_jac_g_structure_val, a2_eval_jac_constraints) 

(sum(a2_eval_constraints)-sum(a1_eval_constraints))/diff
sum(a2_eval_jac_constraints[,idx])
sum(a1_eval_jac_constraints[,idx])
identical(round((a2_eval_constraints-a1_eval_constraints)/diff,2),
          round(a2_eval_jac_constraints[,idx],2))

((a2_eval_constraints-a1_eval_constraints)/diff)[1:10]
(a2_eval_jac_constraints[1:10,idx])
a1_eval_jac_constraints[1:10,idx]
num_grad <- (a2_eval_constraints-a1_eval_constraints)/diff
identical(round((a2_eval_constraints-a1_eval_constraints)/diff,2),
          round(a2_eval_jac_constraints[,idx],2))
num_grad[which(round(num_grad,2)!=
                 round(a2_eval_jac_constraints[,idx],2))][1:10]
a2_eval_jac_constraints[,idx][which(round(num_grad,2)!=
                                      round(a2_eval_jac_constraints[,idx],2))][1:10]





