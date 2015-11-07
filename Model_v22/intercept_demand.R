source("constants_mw6.R")
source("debug_functions.R")
#model with averaged station demand, no service level. Trying to get distance
#coefficient right.

get_points_density_grad_places_count_col <- NULL
source("eval_obj_GMM_model5.R")
source("GetDataGooglePlaces_aggmonth_averaged_pr6.R")
#v0_vec <- c(0)
v0_vec <- generate_v0(20)

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
theta1 <- c(-15.3761554117194, 6, -7.75959563777998,
           0, 0, 1.817587, 0, 0, 0, 0, 0, 0, 0, 0)

active_coef <- c(1:length(theta1))[-c(density_add_den_cols[1],density_metro_evening_col,density_ridership_col)]
active_coef_nointercept <- c(1:length(theta1))[-c(density_intercept_col,density_add_den_cols[1],density_metro_evening_col,
                                                 density_ridership_col)]


#choose starting values for deltain_stin
delta_all <- rep(-6.424,nrow(wdcMerged))
#delta_all <- compute_delta_list_cntrt_map_new(c(theta[1],0,theta[-1]),wdcMerged, points)
deltain_stin <- delta_all[which(wdcMerged$stocked_out==F)]
# #overriding above null setting as GMM values are quite small for identity weighing_GMM_mat setting
# weighing_GMM_mat <<- eval_GMM_optimal_weighing_matrix(c(theta[1],0,theta[-1]),deltain_stin,wdcMerged,
#                                                       points)

get_total_density(theta1, wdcMerged, points)

dem <- eval_lambda_full(deltain_stin, theta1, wdcMerged, points)

#need to see how much demand contribution comes from each point for a given station.
#demand contribution for station 1:
#first get points which can contribute demand to this station, which will be points within max_walking_distance and 
points_type1 <- subset(points, type==1)

focal_station <- 1
demcontri_point1 <- rep(0, nrow(points_type1))

for(i in 1:nrow(points_type1)) {
  stids <- splitchar(points_type1$local_stations[i])  
  if(focal_station %in% stids) {    
    demand_i <- eval_lambda_full(deltain_stin, theta1, wdcMerged, points_type1[i,])
    demcontri_point1[i] <- demand_i[focal_station]
  }  
}
#plot these points and with bubbles represent their contribution
pointsid <- which(demcontri_point1>0)
summary(demcontri_point1[pointsid])
#nearby stations lat, lon to focal_station
local_stations_id <- splitchar(wdcMerged$local_stations[focal_station])
local_stations_id <- local_stations_id[-which(local_stations_id==focal_station)]

symbols(x=c(points_type1$lon[pointsid], wdcMerged$lon[c(focal_station,local_stations_id)]), 
        y=c(points_type1$lat[pointsid], wdcMerged$lat[c(focal_station,local_stations_id)]), 
        circles=c(sqrt(demcontri_point1[pointsid]),
                  rep(max(sqrt(demcontri_point1[pointsid])),1+length(local_stations_id))), 
        inches=0.05, fg=c(rep("gray30", length(pointsid)),"red",rep("green",length(local_stations_id))),
        xlab="lon",ylab="lat")

sum(demcontri_point1[order(demcontri_point1, decreasing = T)[1:400]])
sum(demcontri_point1)



