
# eval_lambda_delta_list <- function(deltain_tw, theta1, wdcMergedday, points, tw_groupin) {
#   no_st <- max(wdcMergedday$station_id_index)
#   if(length(deltain_tw)!=no_st) stop("error in eval_lambda_delta_list")  
#   sto_state_local <- wdcMergedday$sto_state_local
#   local_stations <- wdcMergedday$local_stations
#   points_local_stations <- points$local_stations
#   wdcMergedday  = wdcMergedday[,c("station_id",
#                                   "stocked_out","station_id_index","lat","lon","obs_weight","out_dem_sum")]
#   
#   wdcMergedday = as.matrix(wdcMergedday)
#   points_mat = as.matrix(points[,c("lat","lon","density")])
#   
#   
#   res <- eval_lambda_delta_list_cpp(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
#              as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations))
#   
#   lambda_t <- res[,1]
#   grad_t <- res[,c(2:ncol(res))]
#   return(list("objective"=lambda_t,
#               "gradient"=grad_t))  
# 
# }

eval_lambda_delta_list_new <- function(deltain_tw, theta1, wdcMergedday, points, tw_groupin) {
  no_st <- length(unique(wdcMergedday$station_id_index))
  no_obs <- nrow(wdcMergedday)
  tw_in <- wdcMergedday$tw[1]
  if(length(deltain_tw)!=no_obs) stop("error in eval_lambda_delta_list")  
  sto_state_local <- wdcMergedday$sto_state_local
  local_stations <- wdcMergedday$local_stations
  points_local_stations <- points$local_stations
  wdcMergedday  = wdcMergedday[,c("station_id",
                                  "stocked_out","station_id_index","lat","lon","obs_weight","out_dem_sum")]
  
  
  points_mat <- points
  points_mat$density <- get_points_density(points_mat, theta1, tw_in)
  points_mat = as.matrix(points_mat[,c("lat","lon","density")])
  wdcMergedday = as.matrix(wdcMergedday)
  
  res <- eval_lambda_delta_list_cpp_new(deltain_tw,theta1,wdcMergedday,points_mat,no_st,max_walking_dis,v0_vec,
                                    as.character(sto_state_local), as.character(local_stations), as.character(points_local_stations))
  
  lambda_t <- res[,1]
  grad_t <- res[,c(2:ncol(res))]
  
  return(list("objective"=lambda_t,
              "gradient"=grad_t))  
  
}

# eval_share_diff_list_single <- function(expdeltain_single, theta1, wdcMerged, points, tw_groupin) {
#   #calls eval_share_diff_list for same values of deltain for all stations
#   no_st <- max(wdcMerged$station_id_index)
#   deltain <- rep(deltain_single,no_st)
#   ret <- eval_share_diff_list(expdeltain, theta1, wdcMerged, points, tw_groupin)
#   
#   return( list( "objective" = ret$objective,
#                 "gradient"  = sum(ret$gradient)) )  
# }
# 
# eval_share_diff_list <- function(deltain, theta1, wdcMerged, points, tw_groupin) {
#   #generates the gradient wrt \delta
#   #of eval_share_diff objective function. 
#   #calls eval_grad_lambda_delta and then aggregates.
#   
#   lambda_list <- eval_lambda_delta_list(deltain, theta1, wdcMerged, points, tw_groupin)
#   lambda_t <- lambda_list[[1]]
#   grad_t <- lambda_list[[2]]
#   wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
#   dem_T <- by(lambda_t, wdcMergedday[,"station_id_index"], FUN=function(x) x[1])  
#   dem_hat_T <- by(wdcMergedday$out_dem_sum, wdcMergedday[,"station_id_index"], FUN=function(x) x[1])  
#   dem_T <- as.numeric(dem_T)
#   dem_hat_T <- as.numeric(dem_hat_T)
#   
#   diff_dem <-  dem_T - dem_hat_T  
#   diff_dem_norm <- norm(diff_dem, "2")
#   #normalize the objective.
#   no_st <- max(wdcMergedday$station_id_index)
#   avg_obs_weight <- mean(dem_hat_T)/10
#   diff_dem_norm <- diff_dem_norm/avg_obs_weight
#   
#   obj <- diff_dem_norm*diff_dem_norm
#   
#   grad_dem_T <- by(grad_t,wdcMergedday[,"station_id_index"], FUN=function(x) (x[1,]))  
#   #might replace this with 
#   #matrix(unlist(grad_dem_T),nrow=27,ncol=27)
#   #if the size of matrix grows and rbind is expensive
#   grad_dem_T <- do.call(rbind,grad_dem_T)
#   
#   #grad_dem_T <- as.numeric(grad_dem_T)
#   grad_share_diff <- 2*diff_dem*grad_dem_T/(avg_obs_weight*avg_obs_weight)
#   grad_share_diff <- colSums(grad_share_diff)
#   #   print("deltain: ")
#   #   print(deltain)
#   #   print("gradients: ")
#   #   print(grad_share_diff)
#   return( list( "objective" = obj,
#                 "gradient"  = grad_share_diff) )
#   
# }

# eval_share_diff_list_obj <- function(deltain, theta1, wdcMerged, points, tw_groupin) {
#   return(eval_share_diff_list(deltain, theta1, wdcMerged, points, tw_groupin)$objective)
# }
# 
# eval_share_diff_list_grad <- function(deltain, theta1, wdcMerged, points, tw_groupin) {
#   return(eval_share_diff_list(deltain, theta1, wdcMerged, points, tw_groupin)$gradient)
# }

  
# compute_delta_list <- function(theta1, wdcMerged, points) {
#   #just a wrapper around the optimization routine to compute optimal delta
#   #operated over eval_share_diff and eval_grad_share_diff 
#   print_lvl <- 5
#   no_st <- length(unique(wdcMerged$station_id_index))  
#   delta <- c()  
#   
#   ptm <- proc.time()
#   lb_val <- -200
#   ub_val <- 0
#   
#   compute_vals_seq <- function(val_seq,tw_groupin) {
#     val_list <- c()
#     for(val in val_seq) {      
#       val_list <- c(val_list,eval_share_logdiff_list_single(exp(val),theta1,wdcMerged,points,tw_groupin)$objective)
#     }
#     return(val_list)
#   }
#   
#   
#   tw_group_list <- unique(wdcMerged$tw_group)
#   for(tw_groupin in tw_group_list) {
# #     eval_h_structure <- list()
# #     for(i in 1:(no_st)) {
# #       eval_h_structure[[length(eval_h_structure)+1]] <- c(1:(no_st))
# #     }
#     
#     #compute at steps of 10 the function and seed neldermead with this value
#     val_seq <- seq(ub_val,lb_val,by=-20)
#     val_list <- compute_vals_seq(val_seq,tw_groupin)
#     val_seq <- seq(val_seq[max(which.min(val_list)-1,1)],
#                    val_seq[min(which.min(val_list)+1,length(val_list))],by=-5)
#     val_list <- compute_vals_seq(val_seq,tw_groupin)
#     val_seq <- seq(val_seq[max(which.min(val_list)-1,1)],
#                    val_seq[min(which.min(val_list)+1,length(val_list))],by=-1)
#     val_list <- compute_vals_seq(val_seq,tw_groupin)
#     
#     opts <- list("tol"=1.0e-8,
#                  #"print_options_documentation"='yes'
#                  #"print_info_string"='yes',
#                  "print_level"=print_lvl,
#                  "file_print_level"=12,
#                  "output_file"="debug_inloop.out")
#     
#     res <- ipoptr( x0=exp(rep(val_seq[which.min(val_list)],no_st)),
#                    lb=exp(rep(lb_val,no_st)),
#                    ub=exp(rep(ub_val,no_st)),
#                    eval_f=eval_share_logdiff_list_obj, 
#                    eval_grad_f=eval_share_logdiff_list_grad,
# #                    eval_h=eval_share_logdiff_hessian,
# #                    eval_h_structure=eval_h_structure,
#                    theta1=theta1, 
#                    wdcMerged=wdcMerged, 
#                    points=points, 
#                    tw_groupin=tw_groupin,
#                    opts=opts
#     ) 
#     if(res$status!=0 & res$status!=4) stop("convergence not reached")
#     if(res$objective>1e-4) {
#       print(paste("at theta1: ", theta1))
#       print(paste("at tw_groupin: ", tw_groupin))
#       print(paste("at tract: ", tract))
#       stop("objective value of balance equations non-zero")
#     }
#     #print( res )
#     delta_day <- log(res$solution)     
#     delta <- c(delta, delta_day)        
#   }
# 
#   print(proc.time()-ptm)
#   return(delta)
# }
