compute_delta_list_cntrt_map_new <- function(theta1, wdcMerged, points) {
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
      outside_share = sum(points$density) - mean(wdcMergedday$out_dem_mean)*max(wdcMergedday$station_id_index)
      x0 <- (wdcMergedday$out_dem_mean+ 1e-10)/(outside_share)
      
    }    
    stocked_list <- which(wdcMergedday$stocked_out==FALSE)
    ftol <- 1e-14
    htol <- 1e-16
    itr = 1
    repeat{
      ret <- eval_share_log_list_new(x0,theta1,wdcMergedday,points,tw_groupin)
      dem_T <- as.numeric(ret$dem_T)[stocked_list]
      dem_hat_T <- as.numeric(ret$dem_hat_T)[stocked_list]      
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
      
      #if(norm <= ftol | avgnorm < ftol*1e-3) {
      if(norm <= ftol) {
        break
      }
      #cat(norm(log(x0[stocked_list])-log(x0prev),"2")/length(x0[stocked_list])," ;; ")
      cat(max(abs(log(x0[stocked_list])-log(x0prev)))," ;; ")
    }
    if(obj>1e-5) stop("convergecne not reached")
    print(paste("itr: ", itr))
    #print( res )
    delta_day <- log(x0)     
    delta <- c(delta, delta_day)        
  }
  x0_start <<- exp(delta)
  print(proc.time()-ptm)
  return(delta)
}
