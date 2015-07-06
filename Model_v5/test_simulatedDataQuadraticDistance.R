#testSimulatedData for quadratic distance

source("data_estimation_discoef_2.89_nldis_limstostate_mw6.R")

points_save <- points
wdcMerged_store <- wdcMerged
user_serv_lvl_store <- user_serv_lvl
#make the weights to 1
#put a value of delta
#compute lambda
wdcMerged <- wdcMerged_store
user_serv_lvl <- user_serv_lvl_store
points <- points_save
#keeping "6_0" "6_3" and "7_0" "7_3"
wdcMerged <- subset(wdcMerged, tw_group_fac %in% c("6_0","6_3","7_0","7_3"))
wdcMerged <- droplevels(wdcMerged)
levels(wdcMerged$tw_group_fac)
user_serv_lvl <- subset(user_serv_lvl, tw_group %in% c("6_0","6_3","7_0","7_3"))
user_serv_lvl <- droplevels(user_serv_lvl)
unique(user_serv_lvl$tw_group)

wdcMerged$out_dem_sum <- wdcMerged$out_dem_sum/wdcMerged$obs_weight
wdcMerged$obs_weight <- 1
#removing really low demand
wdcMerged$out_dem_sum[which(wdcMerged$out_dem_sum<=0.01 & 
                              wdcMerged$stocked_out==FALSE)] <- 0.01
ptm <- proc.time()
theta1 <- c(-5,-10,-10)


deltain <- rep(-3,nrow(wdcMerged))
tw_groupin_list <- unique(wdcMerged$tw_group)
no_st <- max(wdcMerged$station_id_index)
grad_lambda_theta <- c()
lambda_t <- c()
for(i in 1:length(tw_groupin_list)) {
  tw_groupin <- tw_groupin_list[i]
  idx <- which(wdcMerged$tw_group==tw_groupin)
  delta_day <- deltain[idx]
  wdcMergedday <- wdcMerged[idx,]
  lambda_t <- c(lambda_t,
                eval_lambda_delta_list_new(delta_day, theta1, wdcMergedday, points, tw_groupin)[[1]]
  )
}
print(proc.time()-ptm)

wdcMerged$out_dem_sum <- lambda_t
##market share rough
st_demand <- by(wdcMerged$out_dem_sum, wdcMerged$station_id_index, FUN=mean)
print(paste0("market share: ", round(sum(st_demand)/sum(points$density)*100),"%"))

##estimate
prev_theta <<- NULL
source("reg_results_GMM_discoef_2.86_2_nldis.R")



##################################################
#run gradient test
run_gradient_test <- function() {
  params <- c(-5.01,-9.90)
  res <- check.derivatives(
    .x=params,
    func=eval_obj_GMM_list_extended_new_obj,
    func_grad=eval_obj_GMM_list_extended_new_grad,
    check_derivatives_print='all',  
    wdcMerged=wdcMerged, 
    points=points  
  )
  if(length(which(res$flag_derivative_warning==TRUE))>0) stop("eval_obj_GMM_list_extended_new test failed")  
  
}



