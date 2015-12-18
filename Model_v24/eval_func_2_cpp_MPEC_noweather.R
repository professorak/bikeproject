eval_error_xi_sl_MPEC <- function(theta1,deltain,wdcMerged,points) {
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  no_st <- max(wdcMerged$station_id_index)
  
  i=1
  serv_lvl_covar <- c()
  instr_serv_lvl_covar <- c()
  for(tw_groupin in tw_groupin_list) {
    wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
    tw_index <- which(tw_groupin_list==tw_groupin)
    deltain_tw <- deltain[which(wdcMerged$tw_group==tw_groupin)]
    low = (i-1)*no_st + 1
    high = i*no_st
    serv_lvl = user_serv_lvl$serv_lvl[c(low:high)]
    serv_lvl_covar <- c(serv_lvl_covar, serv_lvl[wdcMergedday$station_id_index])
    instr_serv_lvl = user_serv_lvl$instr_serv_lvl[c(low:high)]      
    instr_serv_lvl_covar <- c(instr_serv_lvl_covar, instr_serv_lvl[wdcMergedday$station_id_index])
    i= i+1
  }
  
  stocked_list <- which(wdcMerged$stocked_out==FALSE)
  #deltain is the length of stocked_list
  if(length(deltain)!=length(stocked_list)) stop("deltain is not of length stocked_list")
  deltain_weighted <- deltain*(sqrt(wdcMerged$obs_weight)[stocked_list])
  deltain_sum <- ave(deltain_weighted,wdcMerged$station_id_index[stocked_list], FUN=sum)
  weight_sum <- ave(sqrt(wdcMerged$obs_weight)[stocked_list],
                    wdcMerged$station_id_index[stocked_list], FUN=sum)
  deltain_ave <- deltain_sum/weight_sum*(sqrt(wdcMerged$obs_weight)[stocked_list])
  deltain_ddot <- deltain_weighted - deltain_ave
  
  wdcMerged$serv_lvl_covar <- serv_lvl_covar  
  if(length(levels(wdcMerged$tract_tw_fac))>1 & length(levels(wdcMerged$week_fac))>1) {
    X <- model.matrix(~serv_lvl_covar + week_fac + tract_tw_fac + 0, data=wdcMerged[stocked_list,])
    #to be able to estimate station fixed effects
    #else below has tract fixed effects which messes up.
    X <- X[,which(!((colnames(X) %like% "tract_tw_fac") & (colnames(X) %like% 
                                                             paste0("_",min(wdcMerged$tw),"$"))))]    
    #since there is no intercept, below in not automatically removed.
    X <- X[,which(!(colnames(X) %like% paste0("week_fac",min(wdcMerged$week))))]    
  } else if (length(levels(wdcMerged$tract_tw_fac))>1) {
    X <- model.matrix(~serv_lvl_covar + tract_tw_fac + 0, data=wdcMerged[stocked_list,])
    X <- X[,which(!((colnames(X) %like% "tract_tw_fac") & (colnames(X) %like% 
                                                             paste0("_",min(wdcMerged$tw),"$"))))]        
  } else if (length(levels(wdcMerged$week_fac))>1) {
    X <- model.matrix(~serv_lvl_covar + week_fac + 0, data=wdcMerged[stocked_list,])
    X <- X[,which(!(colnames(X) %like% paste0("week_fac",min(wdcMerged$week))))]    
  } else {
    X <- matrix(1,nrow=length(stocked_list), ncol=1)
  }  
  X_weighted <- X*(sqrt(wdcMerged$obs_weight)[stocked_list])
  X_sum_short <- aggregate(X_weighted, list(Stid=wdcMerged$station_id_index[stocked_list]), sum)
  X_sum_short <- X_sum_short[,-1]
  X_sum <- X_sum_short[wdcMerged$station_id_index[stocked_list],]
  X_ave <- X_sum/weight_sum*(sqrt(wdcMerged$obs_weight)[stocked_list])
  X_ddot <- X_weighted - X_ave
  #remove one entry per station
  st_index <- ave(wdcMerged$station_id_index[stocked_list],
                  by=wdcMerged$station_id_index[stocked_list], FUN=function(x)c(1:length(x)))
  X_ddot <- X_ddot[which(st_index!=1),]
  deltain_ddot <- deltain_ddot[which(st_index!=1)]
  X_ddot <- as.matrix(X_ddot)
  #   ########
  #   #singlularity test
  #   a <- colSums(X_ddot)
  #   each(min, max) (a)
  #   ########
  theta2 = solve(t(X_ddot)%*%X_ddot) %*% (t(X_ddot)%*%deltain_ddot)
  
  X1 <- matrix(rep(1,length(stocked_list)))
  xi_bar <- deltain - X%*%theta2
  xi_norm <- xi_bar - X1%*%(t(wdcMerged$obs_weight[stocked_list])%*%xi_bar)/(sum(wdcMerged$obs_weight[stocked_list]))  
  #xi_norm is unweighted value of \eta corresponding to the document.
  xi_weighted <- xi_norm*(sqrt(wdcMerged$obs_weight)[stocked_list])
  
  obj <- c(t(xi_weighted ) %*% xi_weighted )/length(xi_weighted)
  
  #gradient computation #refer to xi_gradient_formulation.lyx for analytical for below.
  weight_sum <- aggregate(sqrt(wdcMerged$obs_weight)[stocked_list], list(Stid=wdcMerged$station_id_index[stocked_list]), sum)
  weight_sum <- weight_sum[order(weight_sum$Stid),2]
  weight_sum_inv <- 1/weight_sum
  wdcMerged_Stid_mat <- model.matrix(~ factor(station_id_index) + 0, data=wdcMerged[stocked_list,])
  row.names(wdcMerged_Stid_mat) <- NULL 
  colnames(wdcMerged_Stid_mat) <- NULL
  grad_delta_hat_wrt_delta <- t(t(wdcMerged_Stid_mat) * weight_sum_inv)
  grad_obj_wrt_delta_1 <- xi_norm*(wdcMerged$obs_weight[stocked_list])
  grad_obj_wrt_delta_2 <- t(grad_obj_wrt_delta_1) %*% X %*% 
    solve(t(X_ddot)%*%X_ddot) %*% t(X_ddot)
  grad_obj_wrt_delta_2 <- grad_obj_wrt_delta_2 * (sqrt(wdcMerged$obs_weight[stocked_list])[which(st_index!=1)])
  grad_obj_wrt_delta_2_full <- rep(0,length(stocked_list))
  grad_obj_wrt_delta_2_full[which(st_index!=1)] <- grad_obj_wrt_delta_2
  rm(grad_obj_wrt_delta_2)
  grad_obj_wrt_delta_3 <- grad_obj_wrt_delta_2_full %*%  grad_delta_hat_wrt_delta
  grad_obj_wrt_delta_3_full <- grad_obj_wrt_delta_3[wdcMerged$station_id_index[stocked_list]]
  grad_obj_wrt_delta_3_full <- grad_obj_wrt_delta_3_full* sqrt(wdcMerged$obs_weight[stocked_list])
  rm(grad_obj_wrt_delta_3)
  grad_obj_wrt_delta <- (c(grad_obj_wrt_delta_1)-c(grad_obj_wrt_delta_2_full)+c(grad_obj_wrt_delta_3_full)
  )*2/length(xi_norm)
  
  
  return(list("obj"=obj,
              "theta2"=theta2,
              "grad_obj_wrt_delta"=grad_obj_wrt_delta))  
  
}
