#This file generates data for observations weights=1 and assumed values of deltas and 
#then estimates them

source("data_estimation_2.6_weather_saved.R")
source("eval_func_2_cpp_MPEC.R")
points_save <- points
wdcMerged_store <- wdcMerged
user_serv_lvl_store <- user_serv_lvl
#make the weights to 1
#put a value of delta
#compute lambda
colnames_theta1 <<- c("dis coef","rand coef","density ridership","density metro","density intercept")
density_ridership_col <<- 3
density_metro_col <<- 4
density_intercept_col <<- 5
density_metro_evening_col <<- 6

wdcMerged <- wdcMerged_store
user_serv_lvl <- user_serv_lvl_store
points <- points_save

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
v0_vec <- generate_v0(10)


# 
# # ptm <- proc.time()
# # theta1 <- c(-5,0)
# # 
# # 
# # deltain <- rnorm(nrow(wdcMerged),-10,5)
# # tw_groupin_list <- unique(wdcMerged$tw_group)
# # no_st <- max(wdcMerged$station_id_index)
# # grad_lambda_theta <- c()
# # lambda_t <- c()
# # for(i in 1:length(tw_groupin_list)) {
# #   tw_groupin <- tw_groupin_list[i]
# #   idx <- which(wdcMerged$tw_group==tw_groupin)
# #   delta_day <- deltain[idx]
# #   wdcMergedday <- wdcMerged[idx,]
# #   lambda_t <- c(lambda_t,
# #                 eval_lambda_delta_list_new(delta_day, theta1, wdcMergedday, points, tw_groupin)[[1]]
# #   )
# # }
# # print(proc.time()-ptm)
# # 
# # wdcMerged$out_dem_sum <- lambda_t
# ##market share rough
# wdcMerged_nostckout <- wdcMerged[which(wdcMerged$stocked_out==F),]
# st_demand <- by(wdcMerged_nostckout$out_dem_sum, wdcMerged_nostckout$station_id_index, FUN=sum)
# st_weight <- by(wdcMerged_nostckout$obs_weight, wdcMerged_nostckout$station_id_index, FUN=sum)
# st_demand <- st_demand/st_weight
# print(paste0("market share: ", round(sum(st_demand)/sum(points$density)*100),"%"))
# 
# 
# 
# ##estimate
# source("test_alternate_contraction_mapping.R")
# prev_theta <<- NULL
# source("reg_results_GMM_discoef_2.89.R")


##################################################
#run gradient test
run_gradient_test <- function() {
  x0_start <<- NULL
  prev_theta <<- NULL
  prev_deltain <<- NULL
  
  theta <- c(-10.7374110745,0.8860478015,403.4177015742,2.8258972200)
  deltain_stin <- rnorm(length(which(!wdcMerged$stocked_out)),-3, 0.5)
  params <- c(theta,deltain_stin)  
  
  res <- check.derivatives(
    .x=params,
    func=eval_obj_GMM_MPEC_obj,
    func_grad=eval_obj_GMM_MPEC_grad,
    check_derivatives_print='all',  
    wdcMerged=wdcMerged, 
    points=points  
  )
  if(length(which(res$flag_derivative_warning==TRUE))>0) stop("eval_obj_GMM_list_extended_new test failed")  
  
}

#test gradient manually for fewer dimensions
theta <- c(-10.7374110745,1,0.8860478015,403.4177015742,2.8258972200, 200)
set.seed(34675)
deltain_stin <- rnorm(length(which(!wdcMerged$stocked_out)),-3, 0.5)
params <- c(theta,deltain_stin)  

a1 <- eval_obj_GMM_MPEC_obj(params,wdcMerged, points, length(theta))
a1_g <- eval_obj_GMM_MPEC_grad(params,wdcMerged, points, length(theta))

idx <- 10
diff <- 0.0001

params_2 <- params
params_2[idx] <- params_2[idx] + diff
a2 <- eval_obj_GMM_MPEC_obj(params_2,wdcMerged, points, length(theta))
a2_g <- eval_obj_GMM_MPEC_grad(params_2,wdcMerged, points, length(theta))

(a2-a1)/diff
a2_g[idx]
a1_g[idx]

##~~~~~~
#gradient test for constrains to be able to check theta gradients as well
theta <- c(-10.7374110745,1,10.8860478015,403.4177015742,2.8258972200, 200)
set.seed(34675)
deltain_stin <- rnorm(length(which(!wdcMerged$stocked_out)),-3, 0.5)
params <- c(theta,deltain_stin)  

a1_eval_constraints <- eval_g(params,wdcMerged, points, length(theta))
a1_eval_jac_constraints <- eval_jac_g(params,wdcMerged, points, length(theta))
eval_jac_g_structure_val <- eval_jac_g_structure(params, wdcMerged, points, length(theta))
if(length(a1_eval_jac_constraints)!=length(unlist(eval_jac_g_structure_val))) stop("gradient lengths dont match a1_eval_jac_constraints")
a1_eval_jac_constraints <- getfullfromsparsematrix (eval_jac_g_structure_val, a1_eval_jac_constraints) 

idx <- 2
diff <- 0.0001

params_2 <- params
params_2[idx] <- params_2[idx] + diff  
  
a2_eval_constraints <- eval_g(params_2,wdcMerged, points, length(theta))
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

### test gradient for density

theta <- c(-10.7374110745,1,10.8860478015,403.4177015742,2.8258972200, 200)

a1_val <- get_total_density(theta, wdcMerged,points)
a1_grad <- get_grad_total_density(theta, wdcMerged,points)

idx <- 1
diff <- 0.0001
theta_2 <- theta
theta_2[idx] <- theta_2[idx] + diff  
a2_val <- get_total_density(theta_2, wdcMerged,points)
a2_grad <- get_grad_total_density(theta_2, wdcMerged,points)

(a2_val-a1_val)/diff
a2_grad[idx]

