#load data from t_r82.7.1_noRD_hyperlocal_majoritystates_preweatherreg_p6_dis1to20_nldis_full_reg_results_GMM_moments_metrodummies.R
source("constants_mw6.R")

#least squares implementation 
#using incoming demand from previous tw as instrument
#using quadratic formulation for service level

get_points_density_grad_places_count_col <- NULL
source("eval_obj_GMM_model5.R")
source("GetDataGooglePlaces_20dis_aggmonth_hyperlocalstate_majoritystates_preweatherreg_moreinstr_withreallocthresh_stkoutthresh_five_averaged_pr6.R")
# {
#   source("GetDataGooglePlaces_2dis_aggmonth_hyperlocalstate_majoritystates_norealloc_stkoutthresh_five_averaged_pr6.R")
#   wdcMerged <- subset(wdcMerged, index<=2)
# }

v0_vec <- c(0)
v0_vec_weights <- rep(1,length(v0_vec))
# load(file="v0_vec.RData")
# load(file="v0_vec_weights.RData")

#log of absolute metro traffic
points_metro_idx <- which(points$type==2)
points$weight[points_metro_idx] <- log(points$weight[points_metro_idx]*22468468)


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
source("eval_covariates_delta_reg_densitycols_setb50.6_wservlvl_metrooff_metrodummies.R")  
# source("eval_covariates_delta_reg_densitycols_setb40.1_wmetrooff_metrodummies.R")

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
delta_all <- compute_delta_list_cntrt_map_new(theta1,wdcMerged, points, 1e-6, 1e-6)
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

#data loaded











source("counterfactuals_functions.R")

# #r82.7.1_noRD_hyperlocal_majoritystates_preweatherreg_metrodummies_moments_1/20_full_nldis_pr6.txt
# theta1 <- c(-2.20961045069201, 0, -19.6523343222199, 0.00523699070208834, 2.12663960984702, 
# 0.131087591513272, 1.53182163420034, 0, 0, 0, 0, 0, 
# 0, 86.6566512013304, 3.78614232898306, 3.12139107706004)
# source("eval_covariates_delta_reg_densitycols_setb50.6_wservlvl_metrooff_metrodummies.R")


#r82.7.1_noRD_googlevar1,4,5_hyperlocal_majoritystates_preweatherreg_metrodummies_moments_1/20_full_nldis_pr6.txt
theta1 <- c(-3.97824788960867, 0, -11.8291928329701, 0.00668806219342966, 
1.91813814237594, 0.0734336213016268, 1.19819644493666, 1.3554490413824, 
0, 0, 0.578749425785955, 1.01076691975751, 0, 58.37584271775, 3.75007002714952, 2.51769163717491)
source("eval_covariates_delta_reg_densitycols_setb50.6_googlevar1,4,5_wservlvl_metrooff_metrodummies.R")


delta_all <- compute_delta_list_cntrt_map_new(theta1,wdcMerged, points, 1e-14, 1e-14)
delta_list <- delta_all[which(wdcMerged$stocked_out==F)]

theta2 <- eval_error_xi_model5 (delta_list, theta1, wdcMerged,points)$theta2/scalecols_X
serv_lvl_coef <- c(theta2[which(row.names(theta2)=="serv_lvl")],
                   theta2[which(row.names(theta2)=="serv_lvl_sq")])

########################################################
########################################################
#generate wdcMerged_cfl1, delta_all_cfl1, and delta_list_cfl1 and generate
#wieghts based on service level. initially from user_serv_lvl
#delta_all_cfl1 are generated based on averages from delta_list
hyperlocal_range <<- 6

retlist <- generate_cf1data(wdcMerged, delta_list, user_serv_lvl)
wdcMerged_cfl1 <- retlist[[1]]
delta_list_cfl1 <- retlist[[2]]
  
#check on total obs_weight for st_tw (should be equal to tw service level)
check_serv_lvl <- data.frame(st_tw_serv_lvl=
  ave(wdcMerged_cfl1$obs_weight, wdcMerged_cfl1$st_tw, FUN=sum))
check_serv_lvl$st_tw <- wdcMerged_cfl1$st_tw
check_serv_lvl <- check_serv_lvl[!duplicated(check_serv_lvl$st_tw),]
if(!identical(check_serv_lvl$st_tw, user_serv_lvl$st_tw)) stop("ln 158: non identical")
if(!identical(round(check_serv_lvl$st_tw_serv_lvl,3), round(user_serv_lvl$serv_lvl,3))) stop("ln 159: non identical")
########################################################
########################################################
#coutnerfactual level 1 with multiplicative change in distance, serv_lvl

ddc  <- demand_distance_change(c(-5:5)/100+1,
                                   theta1,delta_list_cfl1, wdcMerged_cfl1,points)
print(ddc)

sdc <- demand_serv_lvl_change(theta1, delta_list_cfl1, wdcMerged_cfl1,points, 
          c(-5:5)/100+1, serv_lvl_coef)
print(sdc)
sdc_shortterm <- demand_serv_lvl_change_shortterm(theta1, delta_list_cfl1, 
          wdcMerged_cfl1,points, c(-5:5)/100+1, serv_lvl_coef)
print(sdc_shortterm)

# #subs_perc <- average_demand_susbtitution(c(theta1[1],0),wdcMerged,points,delta_list)
# subs_perc <- mean(perc_frac_dem_lost)
# 
# print("ddc: ")
# print(ddc)
# print("sdc: ")
# print(sdc)
# print("subs_perc: ")
# print(subs_perc)
# 
# #change in system-demand. above was system-use
# sdc2 <- demand_serv_lvl_change_new_systemdemand(c(theta1[1],0),wdcMerged,points,0.9,
#                                                 theta2[which(rownames(theta2)=="serv_lvl_covar")])
# print("sdc2: ")
# print(sdc2)



########################################################
########################################################
#perc demand substituted
wdcMerged_sttw_ave <- wdcMerged
wdcMerged_sttw_ave$out_dem_sum <- ave(wdcMerged_sttw_ave$out_dem_sum, 
  wdcMerged_sttw_ave$st_tw, FUN=sum)
wdcMerged_sttw_ave$obs_weight <- ave(wdcMerged_sttw_ave$obs_weight, 
                                      wdcMerged_sttw_ave$st_tw, FUN=sum)
wdcMerged_sttw_ave <- wdcMerged_sttw_ave[!duplicated(wdcMerged_sttw_ave$st_tw),]
delta_list_sttw_ave <- ave(delta_list*wdcMerged$obs_weight, wdcMerged$st_tw,
                           FUN=sum)
delta_list_sttw_ave <- delta_list_sttw_ave[!duplicated(wdcMerged$st_tw)]
delta_list_sttw_ave <- delta_list_sttw_ave/wdcMerged_sttw_ave$obs_weight

perc_demand_substituted <- average_demand_susbtitution_demModel(theta1, wdcMerged_sttw_ave, 
  points,delta_list_sttw_ave)

write.csv(perc_demand_substituted, "PaperPlots/perc_demand_substituted.csv")

perc_frac_dem_lost <- 100-perc_demand_substituted
perc_frac_dem_lost_hist <- hist(perc_frac_dem_lost)
pdf("PaperPlots/PercentageDemandLost.pdf",width=3.2,height=3.10,
    pointsize=9*0.7)
barplot(perc_frac_dem_lost_hist$density*100, xlab="Percentage Demand Lost", 
        ylab="Fraction of Stations (%)",font.main=1, col="grey", 
        main="", xpd=T, names.arg=perc_frac_dem_lost_hist$mids, 
        space=0.5,
        #             xaxt="n"
)
dev.off()

###
#Distance of users
results_distance_segment_demand <- distance_segment_demand(theta1, wdcMerged_sttw_ave, points, delta_list_sttw_ave)
write.csv(results_distance_segment_demand[[2]], "PaperPlots/results_distance_segment_demand.csv")
#keep only first 750mts  
demand_distance_segment <- results_distance_segment_demand[[2]][c(1:75)]
pdf("PaperPlots/Demand_distance_segment_10mts.pdf",width=3.65,height=3.25,
    pointsize=9*0.7)
x <- barplot(demand_distance_segment/sum(demand_distance_segment)*100, xlab="Distance to User", 
             ylab="Fraction of System-Use (%)",font.main=1, col="grey", 
             main="", xpd=T, names.arg=c(1:length(demand_distance_segment))*10, 
             space=0.5,
             #             xaxt="n"
)

dev.off()

########################################################
########################################################
#coutnerfactual level 2


average_station_characteristics_vec <- get_average_station_characteristics(wdcMerged)
mean_st_size <- average_station_characteristics_vec[1]
mean_out_dem <- average_station_characteristics_vec[2]
mean_in_dem <- average_station_characteristics_vec[3]
results_counterfactual_l2 <- run_counterfactual_l2(mean_st_size, mean_out_dem, mean_in_dem,
                                                   theta1, delta_list_cfl1, 
                                                   wdcMerged_cfl1,points, serv_lvl_coef)
write.csv(results_counterfactual_l2,"PaperPlots/results_counterfactual_l2_2.txt")


########################################################
#Isodemand curve

displace_use_vec <- c(0.8,0.9,1,1.1,1.2) 
displace_serv_lvl_vec <- c(-0.2, -0.15 ,-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2)

use_actual <- sum(eval_lambda_full_weighted(delta_list_cfl1, theta1, 
                                            wdcMerged_cfl1, points))

ret_vec_mat <- c()
for(displace_use in displace_use_vec) {
  target_use <- displace_use*use_actual
  ret_vec <- c()
  for(displace_serv_lvl in displace_serv_lvl_vec) {
    ret_vec <- c(ret_vec,optimal_distance_scaling(theta1, delta_list_cfl1, 
                                                  wdcMerged_cfl1, points, displace_serv_lvl, serv_lvl_coef, target_use))
  }  
  ret_vec_mat <- rbind(ret_vec_mat, ret_vec)
}
write.csv(ret_vec_mat, file="isodemand_mat.csv")


########################################################
########################################################
#coutnerfactual level 1 with additive change in distance, serv_lvl

#draw a line how use changes with distance
use_at_dis_x <- function(x, theta1, delta_val) {
  #x is in kms
  prob_x <- exp(theta1[1]*x+theta1[3]*x^2+delta_val)/(1+exp(theta1[1]*x+theta1[3]*x^2+delta_val))
  prob_0 <- exp(delta_val)/(1+exp(delta_val))
  return(prob_x/prob_0)
}

delta_val <- -10

scale_distance <- c(0:100)/100
use_vec <- c()
for(x in scale_distance) {
  use_vec <-  c(use_vec, use_at_dis_x(x,theta1,delta_val)) 
}
use_vec <- use_vec*100

pdf("PaperPlots/CounterfactualL1AdditiveDistanceEffect.pdf",width=3.25,height=3.25,
    pointsize=9*0.7)
plot(scale_distance,
     use_vec, col="white", bty="n", xlab="Distance to User", xaxt="n", ylab="Relative System-Use (%)")
lines(scale_distance,
      use_vec,col="black",lwd=1)
axis(1,at=scale_distance[seq(1,length(scale_distance),by=10)],labels=scale_distance[seq(1,length(scale_distance),by=10)])
dev.off()
idx <- seq(1,length(scale_distance),by=10)
perc_change <- c()
for(i in c(2:length(idx))) {
  perc_change <- c(perc_change, (use_vec[i]-use_vec[i-1])/(use_vec[i]+use_vec[i-1])*2*100)
}
perc_change


#draw a line how use changes with bike-availability
use_at_bikeavailability_x <- function(x, theta1, serv_lvl_coef, delta_val, dis) {
  #x is in kms
  util_x <- theta1[1]*dis+theta1[3]*dis^2+delta_val+serv_lvl_coef[1]*x+serv_lvl_coef[2]*x^2
  prob_x <- exp(util_x)/(1+exp(util_x))*x
  return(prob_x)
}

delta_val <- -10
dis <- 0.1
bike_availability_vec <- c(0:100)/100
use_vec <- c()
for(x in bike_availability_vec) {
  use_vec <-  c(use_vec, use_at_bikeavailability_x(x,theta1,serv_lvl_coef, delta_val, dis)) 
}
use_vec <- use_vec/max(use_vec)*100

pdf("PaperPlots/CounterfactualL1AdditiveBikeAvailabilityEffect.pdf",width=3.25,height=3.25,
    pointsize=9*0.7)
plot(bike_availability_vec,
     use_vec, col="white", bty="n", xlab="Bike-Availability", xaxt="n", ylab="Relative System-Use (%)")
lines(bike_availability_vec,
      use_vec,col="black",lwd=1)
axis(1,at=bike_availability_vec[seq(1,length(bike_availability_vec),by=10)],labels=bike_availability_vec[seq(1,length(bike_availability_vec),by=10)])
dev.off()
idx <- seq(1,length(bike_availability_vec),by=10)
perc_change <- c()
for(i in c(2:length(idx))) {
  perc_change <- c(perc_change, (use_vec[i]-use_vec[i-1])/(use_vec[i]+use_vec[i-1])*2*100)
}
perc_change








