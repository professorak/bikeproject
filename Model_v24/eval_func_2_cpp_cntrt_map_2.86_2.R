
compute_delta_list_cntrt_map_new <- function(theta1, wdcMerged, points,
                                             ftol=1e-14, htol=1e-16) {
  #just a wrapper around the optimization routine to compute optimal delta
  #operated over eval_share_diff and eval_grad_share_diff 
  no_st <- length(unique(wdcMerged$station_id_index))  
  no_obs <- nrow(wdcMerged)
  delta <- c()  
  
  ptm <- proc.time()
  lb_val <- -100
  ub_val <- 100
  
  compute_vals_seq <- function(val_seq,wdcMergedday, tw_groupin) {
    val_list <- c()
    for(val in val_seq) {      
      val_list <- c(val_list,eval_share_logdiff_list_single_new(exp(val),theta1,wdcMergedday,points,tw_groupin)$objective)
    }
    return(val_list)
  }
  
  
  tw_group_list <- unique(wdcMerged$tw_group)
  for(i in 1:length(tw_group_list)) {
    tw_groupin = tw_group_list[i]
    wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
    
    if(!is.null(x0_start)) {
      x0 = x0_start[which(wdcMerged$tw_group==tw_groupin)]
    } else {
      wdcMergedday$out_dem_mean = wdcMergedday$out_dem_sum/wdcMergedday$obs_weight
      points_density <- get_points_density(points, theta1, wdcMergedday$tw[1])
      
      outside_share = sum(points_density) - mean(wdcMergedday$out_dem_mean)*max(wdcMergedday$station_id_index)
      x0 <- (wdcMergedday$out_dem_mean+ 1e-10)/(outside_share)
      
    }    
    stocked_list <- which(wdcMergedday$stocked_out==FALSE)
    
    itr = 1
    repeat{
      ret <- eval_share_log_list_new(x0,theta1,wdcMergedday,points,tw_groupin)
      dem_T <- ret$dem_T
      dem_hat_T <- ret$dem_hat_T      
      #obj <- (norm(dem_T-dem_hat_T,"2"))/length(dem_hat_T)
      obj <- max(abs(dem_T-dem_hat_T)/dem_T)
      cat(obj," ")
      
      cat(mean(log(x0))," ")
      itr = itr + 1
      x0prev <- x0[stocked_list]
      x0[stocked_list] <- x0[stocked_list]*(dem_hat_T/dem_T)
      t= abs(x0[stocked_list]-x0prev)
      norm=max(t)
      avgnorm = mean(t)
      
      #if(norm <= ftol | avgnorm < ftol*1e-3) {
      if(norm <= ftol) {
        break
      }
      #cat(norm(log(x0[stocked_list])-log(x0prev),"2")/length(x0[stocked_list])," ;; ")
      cat(max(abs(log(x0[stocked_list])-log(x0prev)))," ;; ")
    }
    if(ftol < 1e-10) {
      if(obj>1e-5) stop("convergecne not reached")  
    }    
    print(paste("itr: ", itr))
    #print( res )
    delta_day <- log(x0)     
    delta <- c(delta, delta_day)        
  }
  x0_start <<- exp(delta)
  print(proc.time()-ptm)
  return(delta)
}


eval_share_log_list <- function(expdeltain, theta1, wdcMergedday, points, tw_groupin) {
  #generates the gradient wrt \delta
  #of eval_share_diff objective function. 
  #calls eval_grad_lambda_delta and then aggregates.
  if(any(expdeltain<0)) {stop("expdeltain values <=0")}
  deltain = log(expdeltain)
  
  no_st <- max(wdcMergedday$station_id_index)
  lambda_list <- eval_lambda_delta_list(deltain, theta1, wdcMergedday, points, tw_groupin)
  lambda_t <- lambda_list[[1]]

  
  dem_T <- by(lambda_t, wdcMergedday[,"station_id_index"], FUN=sum)  
  dem_hat_T <- by(wdcMergedday$out_dem_sum, wdcMergedday[,"station_id_index"], FUN=sum)  
  
  #   print("deltain: ")
  #   print(deltain)
  #   print("gradients: ")
  #   print(grad_share_diff)
  return( list( "dem_T" = dem_T,
            "dem_hat_T"  = dem_hat_T) )
  
}

eval_share_log_list_new <- function(expdeltain, theta1, wdcMergedday, points, tw_groupin) {
  #generates the gradient wrt \delta
  #of eval_share_diff objective function. 
  #calls eval_grad_lambda_delta and then aggregates.
  if(any(expdeltain<0)) {stop("expdeltain values <=0")}
  deltain = log(expdeltain)
  
  stocked_list <- which(wdcMergedday$stocked_out==FALSE)
  
  dem_T <- eval_lambda_new(deltain, theta1, wdcMergedday, points, tw_groupin)[stocked_list]
  dem_hat_T <- wdcMergedday$out_dem_sum[stocked_list]/wdcMergedday$obs_weight[stocked_list]
  
  return( list( "dem_T" = dem_T,
                "dem_hat_T"  = dem_hat_T) )
}
