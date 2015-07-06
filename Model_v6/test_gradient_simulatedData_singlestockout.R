#This file generates data for observations weights=1 and assumed values of deltas and 
#then estimates them

source("data_estimation_discoef_2.89_singlestostate_saved.R")
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

wdcMerged <- wdcMerged_store
user_serv_lvl <- user_serv_lvl_store
points <- points_save
#keeping "6_0" "6_3" and "7_0" "7_3"
#time_windows <- c("6_3","6_5","7_3","7_5")
time_windows <- c("6_0")
wdcMerged <- subset(wdcMerged, tw_group_fac %in% time_windows)
wdcMerged <- droplevels(wdcMerged)
levels(wdcMerged$tw_group_fac)
user_serv_lvl <- subset(user_serv_lvl, tw_group %in% time_windows)
user_serv_lvl <- droplevels(user_serv_lvl)
unique(user_serv_lvl$tw_group)
#making 0 demand observations insignificant
wdcMerged$obs_weight[which(wdcMerged$stocked_out==F & wdcMerged$out_dem_sum<=0.1)] <- 0.1

# wdcMerged$out_dem_sum <- wdcMerged$out_dem_sum/wdcMerged$obs_weight
# wdcMerged$obs_weight <- 1
#removing really low demand
wdcMerged$out_dem_sum[which(wdcMerged$out_dem_sum<=0.01 & 
                              wdcMerged$stocked_out==FALSE)] <- 0.01
# ptm <- proc.time()
# theta1 <- c(-5,0)
# 
# 
# deltain <- rnorm(nrow(wdcMerged),-10,5)
# tw_groupin_list <- unique(wdcMerged$tw_group)
# no_st <- max(wdcMerged$station_id_index)
# grad_lambda_theta <- c()
# lambda_t <- c()
# for(i in 1:length(tw_groupin_list)) {
#   tw_groupin <- tw_groupin_list[i]
#   idx <- which(wdcMerged$tw_group==tw_groupin)
#   delta_day <- deltain[idx]
#   wdcMergedday <- wdcMerged[idx,]
#   lambda_t <- c(lambda_t,
#                 eval_lambda_delta_list_new(delta_day, theta1, wdcMergedday, points, tw_groupin)[[1]]
#   )
# }
# print(proc.time()-ptm)
# 
# wdcMerged$out_dem_sum <- lambda_t
##market share rough
# wdcMerged_nostckout <- wdcMerged[which(wdcMerged$stocked_out==F),]
# st_demand <- by(wdcMerged_nostckout$out_dem_sum, wdcMerged_nostckout$station_id_index, FUN=sum)
# st_weight <- by(wdcMerged_nostckout$obs_weight, wdcMerged_nostckout$station_id_index, FUN=sum)
# st_demand <- st_demand/st_weight
# print(paste0("market share: ", round(sum(st_demand)/sum(points$density)*100),"%"))
# 


##estimate
#source("test_alternate_contraction_mapping.R")
prev_theta <<- NULL
source("reg_results_GMM_discoef_2.89.R")


##################################################
#run gradient test
run_gradient_test <- function() {
  x0_start <<- NULL
  prev_theta <<- NULL
  prev_deltain <<- NULL
  
  params <- c(-3,10,10,1)
  res <- check.derivatives(
    .x=params,
    func=eval_obj_GMM_list_extended_new_obj,
    func_grad=eval_obj_GMM_list_extended_new_grad,
    check_derivatives_print='all',  
    wdcMerged=wdcMerged, 
    points=points  
  )
  if(length(which(res$flag_derivative_warning==TRUE))>0) stop("eval_obj_GMM_list_extended_new test failed")  
  
  x0_start <<- NULL
  prev_theta <<- NULL
  prev_deltain <<- NULL
  
  params <- c(-3,10,10,1)
  res <- check.derivatives(
    .x=params,
    func=get_total_density,
    func_grad=get_grad_total_density,
    check_derivatives_print='all',  
    wdcMerged=wdcMerged, 
    points=points  
  )
  if(length(which(res$flag_derivative_warning==TRUE))>0) stop("eval_obj_GMM_list_extended_new test failed") 
}





##############################################
#contraction mapping algo here
run_contraction_mapping_test <- function() {
  theta1 <- c(-4.583810004,0)
  x0_start <- NULL
  no_st <- length(unique(wdcMerged$station_id_index))  
  no_obs <- nrow(wdcMerged)
  delta <- c()  
  
  ptm <- proc.time()
  lb_val <- -100
  ub_val <- 100
  
  tw_group_list <- unique(wdcMerged$tw_group)
  i <- 1
  
  tw_groupin = tw_group_list[i]
  wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
  
  if(!is.null(x0_start)) {
    x0 = x0_start[which(wdcMerged$tw_group==tw_groupin)]
  } else {
    wdcMergedday$out_dem_mean = wdcMergedday$out_dem_sum/wdcMergedday$obs_weight
    outside_share = sum(points$density) - mean(wdcMergedday$out_dem_mean)*max(wdcMergedday$station_id_index)
    x0 <- (wdcMergedday$out_dem_mean+ 1e-10)/(outside_share)
    
  }    
  stocked_list <- which(wdcMergedday$stocked_out==FALSE)
  
  version=2
  
  if(version==1) {
    #stocked_list <- stocked_list[1]
    ftol <- 1e-12
    htol <- 1e-12
    itr = 1
    repeat{
      #    for(j in 1:10) {
      ret <- eval_share_log_list_new(x0,theta1,wdcMergedday,points,tw_groupin)
      dem_T <- as.numeric(ret$dem_T + 1e-10)[stocked_list]
      dem_hat_T <- as.numeric(ret$dem_hat_T + 1e-10)[stocked_list]      
      #obj <- (norm(dem_T-dem_hat_T,"2"))/length(dem_hat_T)
      obj <- max(abs(dem_T-dem_hat_T))
      cat(obj," ")
      cat(mean(log(x0))," ")
      itr = itr + 1
      x0prev <- x0[stocked_list]
      x0[stocked_list] <- x0[stocked_list]*(dem_hat_T/dem_T)
      #if(norm(log(x0[stocked_list])-log(x0prev),"2")/length(x0[stocked_list]) < ftol | obj < htol) {
      if(max(abs(log(x0[stocked_list])-log(x0prev))) < ftol | obj < htol) {
        #      if(norm(log(x0)-log(x0prev),"2")/length(x0) < ftol ){
        break
      }
      #cat(norm(log(x0[stocked_list])-log(x0prev),"2")/length(x0[stocked_list])," ;; ")
      cat(max(abs(log(x0[stocked_list])-log(x0prev)))," ;; ")
    }
    if(obj>1e-5) stop("convergecne not reached")
    print(paste("itr: ", itr))
    #print( res )
    delta_day <- log(x0)       
  }
  
  if(version==2) {
    #stocked_list <- stocked_list[1]
    ftol <- 1e-9
    htol <- 1e-9
    itr = 1
    repeat{
      #    for(j in 1:10) {
      ret <- eval_share_log_list_new(x0,theta1,wdcMergedday,points,tw_groupin)
      dem_T <- as.numeric(ret$dem_T + 1e-10)[stocked_list]
      dem_hat_T <- as.numeric(ret$dem_hat_T + 1e-10)[stocked_list]      
      #obj <- (norm(dem_T-dem_hat_T,"2"))/length(dem_hat_T)
      obj <- max(abs(dem_T-dem_hat_T))
      cat(obj," ")
      
      cat(mean(log(x0))," ")
      itr = itr + 1
      x0prev <- x0[stocked_list]
      x0[stocked_list] <- x0[stocked_list]*(dem_hat_T/dem_T)
      t= abs(x0[stocked_list]-x0prev)
      norm=max(t)
      avgnorm = mean(t)
      
      if(norm <= ftol | avgnorm < ftol*1e-3) {
        break
      }
      #cat(norm(log(x0[stocked_list])-log(x0prev),"2")/length(x0[stocked_list])," ;; ")
      cat(max(abs(log(x0[stocked_list])-log(x0prev)))," ;; ")
    }
    if(obj>1e-5) stop("convergecne not reached")
    print(paste("itr: ", itr))
    #print( res )
    delta_day <- log(x0)       
  }  
}












