weighted_GMM <<- 1
library("plm")
# eval_obj_GMM_list_extended_new <- function(theta, wdcMerged, points) {  
#   #theta=(theta1,\alpha)
#   ptm <- proc.time()
#   theta1 <- theta[1:2]  
#   print("theta:")
#   print(paste(theta))
#   
#   tryCatch({
#     tw_groupin_list <- unique(wdcMerged$tw_group)    
#     no_st <- max(wdcMerged$station_id_index)
#     
#     ret <- eval_error_xi(theta1,wdcMerged,points)
#     xi_norm=ret$xi
#     theta2=ret$theta2
#     X1=ret$X1
#     Z=ret$Z
#     W=ret$W
#     grad_delta_theta = ret$grad_delta_theta
#     
#     
#     moments = c(t(Z) %*% xi_norm)
#     
#     gmm_objective = moments %*% W %*% moments/length(xi_norm)
#     #need to compute gradients of these moment conditins as well, modify grad_delta_theta to grad_deltanorm_theta
#     #Since deltanorm is delta-mean(delta) (in a time window), grad_delta_theta_norm is also computed by substracting mean from gradients
#     #     grad_delta_ave <- cbind(ave(grad_delta_theta[,1],deltain_group,FUN=mean),
#     #                             ave(grad_delta_theta[,2],deltain_group,FUN=mean))
#     grad_xi_norm <- grad_delta_theta #- grad_delta_ave # as grad_delta_theta is same as grad_xi_theta because xi is linearly related to theta2
#     moments_grad <- t(Z) %*% grad_xi_norm
#     gmm_gradient <- 2*moments %*% W %*% moments_grad/length(xi_norm)
#     
#     #     gmm_objective = moment_2*moment_2
#     #     gmm_gradient = 2 * moment_2 * moments_2_grad
#   },error = function(e) {
#     #error-handler-code
#     print(paste("MY_ERROR:  ",e))
#     print("theta:")
#     print(paste(theta))
#     ret = 1e12
#     gradient = c(0,0)
#   }      
#   )
#   print("obj:")
#   print(gmm_objective)
#   print("theta:")
#   print(theta)
#   print("gradient:")
#   print(gmm_gradient)
#   print("eval_obj time:")
#   print(proc.time()-ptm)
#   return( list( "objective" = gmm_objective,
#                 "gradient"  = gmm_gradient ) )
#   
# }

eval_obj_GMM_list_extended_new <- function(theta1, wdcMerged, points) {  
  #theta=(theta1,\alpha)
  ptm <- proc.time()
  #theta1 <- c(theta[1],0)  

  print("theta1:")
  print(paste(theta1))
  
  tryCatch({
    tw_groupin_list <- unique(wdcMerged$tw_group)    
    no_st <- max(wdcMerged$station_id_index)
    
    ret <- eval_error_xi(theta1,wdcMerged,points)
    xi_norm=ret$xi
    theta2=ret$theta2
    X1=ret$X1
    Z=ret$Z
    W=ret$W
    grad_delta_theta = ret$grad_delta_theta
    
    gmm_objective = c(t(xi_norm ) %*% xi_norm )/length(xi_norm)
    grad_xi_norm <- grad_delta_theta #- grad_delta_ave # as grad_delta_theta is same as grad_xi_theta because xi is linearly related to theta2
    #moments_grad <- t(Z) %*% grad_xi_norm
    gmm_gradient <- 2*t(xi_norm ) %*% grad_xi_norm/length(xi_norm)
    
    #     gmm_objective = moment_2*moment_2
    #     gmm_gradient = 2 * moment_2 * moments_2_grad
  },error = function(e) {
    #error-handler-code
    print(paste("MY_ERROR:  ",e))
    print("theta:")
    print(paste(theta))
    ret = 1e12
    gradient = c(0,0)
  }      
  )
  print("obj:")
  print(gmm_objective)
  print("theta:")
  print(theta)
  print("gradient:")
  print(gmm_gradient)
  print("eval_obj time:")
  print(proc.time()-ptm)
  return( list( "objective" = gmm_objective,
                "gradient"  = gmm_gradient[-2] ) )
  
}

eval_obj_GMM_list_extended_new_obj <- function(theta, wdcMerged, points) {    
  #obj_eval <<- 1
  print("objective comput")
  return(eval_obj_GMM_list_extended_new(theta, wdcMerged, points)$objective)
}  
eval_obj_GMM_list_extended_new_grad <- function(theta, wdcMerged, points) {    
  #obj_eval <<- 0
  print("gradint comput")
  return(eval_obj_GMM_list_extended_new(theta, wdcMerged, points)$gradient)
} 


eval_weighting_mat_new <- function(theta, wdcMerged, points) {  
  
  #theta=(theta1,\alpha)
  ptm <- proc.time()
  theta1 <- c(theta[1],0)  
  #theta1 <- theta[1]  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  no_st <- max(wdcMerged$station_id_index)
  
  station_data <- unique(wdcMerged[,c("station_id_index","tract")])
  deltain_group = as.factor(paste(rep(tw_groupin_list,each=no_st),rep(station_data$tract,length(tw_groupin_list))))
  tw_group_vec <- as.factor(rep(tw_groupin_list,each=no_st))
  #deltain_ave=ave(deltain,by=deltain_tw_group,FUN=mean)
  #deltain_norm=deltain#-deltain_ave
  
  ret <- eval_error_xi(theta1,wdcMerged,points)
  xi_norm=ret$xi
  theta2=ret$theta2
  X1=ret$X1
  Z=ret$Z      
  W=ret$W
  
  moment_vec <- c()
  for(i in 1:ncol(Z)) {
    moment_vec <- cbind(moment_vec,Z[,i]*xi_norm)
  }
  
  Lambda_hat = (t(moment_vec) %*% moment_vec) / length(xi_norm)
  Lambda_hat_inv = solve(Lambda_hat)
  ###########################
  return(Lambda_hat_inv)
}


eval_standard_errors_new <- function(theta1, deltain, wdcMerged, points) {  
  #evaluates standard error for one step GMM.
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  no_st <- max(wdcMerged$station_id_index)
  
  station_data <- unique(wdcMerged[,c("station_id_index","tract")])
  deltain_group = as.factor(paste(rep(tw_groupin_list,each=no_st),rep(station_data$tract,length(tw_groupin_list))))
  tw_group_vec <- as.factor(rep(tw_groupin_list,each=no_st))
  
  ret <- eval_error_xi(theta1,wdcMerged,points)
  xi_norm=ret$xi
  theta2=ret$theta2
  X1=ret$X1
  Z=ret$Z      
  W=ret$W
  grad_delta_theta = ret$grad_delta_theta
  
  grad_xi_norm <- grad_delta_theta #- grad_delta_ave # as grad_delta_theta is same as grad_xi_theta because xi is linearly related to theta2
  moments_grad <- t(Z) %*% grad_xi_norm
  
  #st_error = solve(t(moments_grad) %*% W_optimal %*% moments_grad)/length(deltain)
  S_0 <- eval_moment_var_new(theta1,wdcMerged,points)      
  m1 = solve(t(moments_grad) %*% W %*% moments_grad)
  m2 = t(moments_grad) %*% W %*% S_0  %*% W %*% moments_grad
  st_error = (m1 %*% m2 %*% m1)/length(deltain)
  return(st_error)
}

eval_covariace_new <- function(theta1, wdcMerged, points) {  
  #evaluates based on M-estimator, using the homogeniety assumption.
  if(length(unique(wdcMerged$tw_group))>1) stop("standard error routine not implemented for multiple tw_group")     
  if(!exists("prev_theta") || !identical(prev_theta,theta1)) {
    deltain <- compute_delta_list_cntrt_map_new(theta1,wdcMerged,points)
    prev_theta <<- theta1
    prev_deltain <<- deltain
  } else {
    deltain <- prev_deltain
  }
  del_xi_theta1 = eval_grad_delta_theta_new(theta1, deltain, wdcMerged, points, 1)
  del_xi_theta1 <- del_xi_theta1[,1]
  ret <- eval_error_xi(theta1, wdcMerged, points)
  X1 <- ret$X1
  
  del_xi <- cbind(del_xi_theta1, -X1)
  
  #s_theta <- 2*deltain*del_xi  
  #B = t(s_theta) %*% s_theta / nrow(wdcMerged)
  B = 2*2* t(del_xi) %*% del_xi  / nrow(wdcMerged) * c(var(ret$xi))
  
  ##
  dis <- c( 0.0001, 0.000000000)
  deltain_dis <- compute_delta_list_cntrt_map_new(theta1 + dis, wdcMerged, points)
  del_xi_theta1_dis = eval_grad_delta_theta_new(theta1 + dis, deltain_dis, wdcMerged, points, 1)
  del_xi_theta1_dis <- del_xi_theta1_dis[,1]
  del2_xi_theta1 <- (del_xi_theta1_dis - del_xi_theta1)/sum(dis)
  
  A_int1 = 2*deltain*del2_xi_theta1 + 2*del_xi_theta1*del_xi_theta1
  A_int1 = sum(A_int1)
  
  Hess = 2* t(del_xi) %*% del_xi
  Hess[1,1] = Hess[1,1] + A_int1
  
  A =   Hess/nrow(wdcMerged)
  
  Var = solve(A)%*%B%*%solve(A)/(nrow(wdcMerged)-ncol(del_xi))
  return(Var)
}



eval_moment_var_new <- function(theta, wdcMerged, points) {  
  #theta=(theta1,\alpha)
  ptm <- proc.time()
  theta1 <- theta[1:2]  
  #theta1 <- theta[1]  
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  no_st <- max(wdcMerged$station_id_index)
  
  station_data <- unique(wdcMerged[,c("station_id_index","tract")])
  deltain_group = as.factor(paste(rep(tw_groupin_list,each=no_st),rep(station_data$tract,length(tw_groupin_list))))
  tw_group_vec <- as.factor(rep(tw_groupin_list,each=no_st))
  #deltain_ave=ave(deltain,by=deltain_tw_group,FUN=mean)
  #deltain_norm=deltain#-deltain_ave
  

  ret <- eval_error_xi(theta1,wdcMerged,points)
  xi_norm=ret$xi
  theta2=ret$theta2
  X1=ret$X1
  Z=ret$Z      
  W=ret$W
  
  moment_vec <- c()
  for(i in 1:ncol(Z)) {
    moment_vec <- cbind(moment_vec,Z[,i]*xi_norm)
  }
  
  Lambda_hat = (t(moment_vec) %*% moment_vec) / length(xi_norm)
  
  return(Lambda_hat)
}

eval_error_xi <- function(theta1,wdcMerged,points) {
  #returns xi, along with theta2, Z and X1
  if(!exists("prev_theta") || !identical(prev_theta,theta1)) {
    deltain <- compute_delta_list_cntrt_map_new(theta1,wdcMerged,points)
    prev_theta <<- theta1
    prev_deltain <<- deltain
  } else {
    deltain <- prev_deltain
  }
  tw_groupin_list <- unique(wdcMerged$tw_group)    
  no_st <- max(wdcMerged$station_id_index)
  
  i=1
  serv_lvl_covar <- c()
  instr_serv_lvl_covar <- c()
  grad_delta_theta <- c()
  for(tw_groupin in tw_groupin_list) {
    wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
    tw_index <- which(tw_groupin_list==tw_groupin)
    deltain_tw <- deltain[which(wdcMerged$tw_group==tw_groupin)]
    grad_delta_theta <- rbind(grad_delta_theta,
                              eval_grad_delta_theta_new(theta1, deltain_tw, wdcMergedday, points, tw_groupin))
    low = (i-1)*no_st + 1
    high = i*no_st
    serv_lvl = user_serv_lvl$serv_lvl[c(low:high)]
    serv_lvl_covar <- c(serv_lvl_covar, serv_lvl[wdcMergedday$station_id_index])
    instr_serv_lvl = user_serv_lvl$instr_serv_lvl[c(low:high)]      
    instr_serv_lvl_covar <- c(instr_serv_lvl_covar, instr_serv_lvl[wdcMergedday$station_id_index])
    i= i+1
  }
  #stop("will have to account for wieghts, create ddot such that fixed effet disapperas even when the observations are weighted")
  #   deltain = deltain * sqrt(wdcMerged$obs_weight)
  #   Z = Z * sqrt(wdcMerged$obs_weight)
  #   X1 = X1 * sqrt(wdcMerged$obs_weight)
  #   grad_delta_theta = grad_delta_theta * sqrt(wdcMerged$obs_weight)
  #   
  #   theta2 = solve(t(X1)%*%Z%*%W%*%t(Z)%*%X1)%*%(t(X1)%*%Z%*%W%*%t(Z)%*%deltain)      
  #   
  stocked_list <- which(wdcMerged$stocked_out==FALSE)
  deltain_weighted <- (deltain*sqrt(wdcMerged$obs_weight))[stocked_list]
  deltain_sum <- ave(deltain_weighted,wdcMerged$station_id_index[stocked_list], FUN=sum)
  weight_sum <- ave(sqrt(wdcMerged$obs_weight)[stocked_list],
                    wdcMerged$station_id_index[stocked_list], FUN=sum)
  deltain_ave <- deltain_sum/weight_sum*(sqrt(wdcMerged$obs_weight)[stocked_list])
  deltain_ddot <- deltain_weighted - deltain_ave
  
  wdcMerged$serv_lvl_covar <- serv_lvl_covar
  X <- model.matrix(~serv_lvl_covar + week_fac + tract_tw_fac + 0, data=wdcMerged[stocked_list,])
  Z <- solve(t(X)%*%X)
  X <- X[,which(!((colnames(X) %like% "tract_tw_fac") & (colnames(X) %like% 
                                                           paste0("_",min(wdcMerged$tw),"$"))))]
  
  X <- X[,which(!(colnames(X) %like% paste0("week_fac",min(wdcMerged$week))))]
  
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
  theta2 = solve(t(X_ddot)%*%X_ddot) %*% (t(X_ddot)%*%deltain_ddot)
  
  grad_delta_theta_weighted <- (grad_delta_theta*sqrt(wdcMerged$obs_weight))[stocked_list,]
  grad_delta_theta_sum1 <- ave(grad_delta_theta_weighted[,1],
                               wdcMerged$station_id_index[stocked_list], FUN=sum)
  grad_delta_theta_sum <- grad_delta_theta_sum1
  for(j in c(2:ncol(grad_delta_theta_weighted))) {
    grad_delta_theta_sum2 <- ave(grad_delta_theta_weighted[,j],
                                 wdcMerged$station_id_index[stocked_list], FUN=sum)
    grad_delta_theta_sum <- cbind(grad_delta_theta_sum,grad_delta_theta_sum2)
  }
  grad_delta_theta_ave <- grad_delta_theta_sum/weight_sum*(sqrt(wdcMerged$obs_weight)[stocked_list])
  grad_delta_theta_ddot <- grad_delta_theta_weighted - grad_delta_theta_ave
  grad_delta_theta_ddot <- grad_delta_theta_ddot[which(st_index!=1),]
  grad_xi_bar_theta <- grad_delta_theta[stocked_list,]  - 
    (X%*%solve(t(X_ddot)%*%X_ddot)%*%(t(X_ddot)%*%grad_delta_theta_ddot))
  
  X1 <- matrix(rep(1,length(stocked_list)))
  grad_xi_theta <- grad_xi_bar_theta - X1%*%(t(wdcMerged$obs_weight[stocked_list])%*%grad_xi_bar_theta)/(sum(wdcMerged$obs_weight[stocked_list]))
  xi_bar <- deltain[stocked_list] - X%*%theta2
  xi_norm <- xi_bar - X1%*%(t(wdcMerged$obs_weight[stocked_list])%*%xi_bar)/(sum(wdcMerged$obs_weight[stocked_list]))
  #weigh now
  xi_norm <- xi_norm * sqrt(wdcMerged$obs_weight[stocked_list])
  xi_norm_full<- rep(0,nrow(wdcMerged))
  xi_norm_full[stocked_list] <- xi_norm
  grad_xi_theta <- grad_xi_theta * sqrt(wdcMerged$obs_weight[stocked_list])
  grad_xi_theta_full <- matrix(0,nrow(wdcMerged),ncol(grad_xi_theta))
  grad_xi_theta_full[stocked_list,] <- grad_xi_theta
  #
  return(list("xi"=xi_norm_full,
              "theta2"=theta2,
              "grad_delta_theta"=grad_xi_theta_full))  
  
}


# eval_error_xi <- function(theta1,wdcMerged,points) {
#   #returns xi, along with theta2, Z and X1
#   if(!exists("prev_theta") || !identical(prev_theta,theta1)) {
#     deltain <- compute_delta_list_cntrt_map_new(theta1,wdcMerged,points)
#     prev_theta <<- theta1
#     prev_deltain <<- deltain
#   } else {
#     deltain <- prev_deltain
#   }
#   tw_groupin_list <- unique(wdcMerged$tw_group)    
#   no_st <- max(wdcMerged$station_id_index)
#   
#   i=1
#   serv_lvl_covar <- c()
#   instr_serv_lvl_covar <- c()
#   grad_delta_theta <- c()
#   for(tw_groupin in tw_groupin_list) {
#     wdcMergedday = subset(wdcMerged, tw_group==tw_groupin)      
#     tw_index <- which(tw_groupin_list==tw_groupin)
#     deltain_tw <- deltain[which(wdcMerged$tw_group==tw_groupin)]
#     grad_delta_theta <- rbind(grad_delta_theta,
#                               eval_grad_delta_theta_new(theta1, deltain_tw, wdcMergedday, points, tw_groupin))
#     low = (i-1)*no_st + 1
#     high = i*no_st
#     serv_lvl = user_serv_lvl$serv_lvl[c(low:high)]
#     serv_lvl_covar <- c(serv_lvl_covar, serv_lvl[wdcMergedday$station_id_index])
#     instr_serv_lvl = user_serv_lvl$instr_serv_lvl[c(low:high)]      
#     instr_serv_lvl_covar <- c(instr_serv_lvl_covar, instr_serv_lvl[wdcMergedday$station_id_index])
#     i= i+1
#   }
#   
#   #   station_data <- unique(wdcMerged[,c("station_id_index","tract")])    
#   #   deltain_group = as.factor(paste(wdcMerged$tw_group,wdcMerged$tract))
#   #   tw_group_vec <- as.factor(wdcMerged$tw_group)
#   #   tract_vec <-as.factor(wdcMerged$tract)
#   #   week_fac <- wdcMerged$week_fac
#   #   tw_fac <- wdcMerged$tw_fac
#   #   station_id_index_vec <- as.factor(wdcMerged$station_id_index)
#   
#   #   covariates_2 = covariates_st$density_stations[wdcMerged$station_id_index]    
#   #   Z = cbind(covariates_2)
#   #   
#   #   if(length(levels(tw_group_vec))<=1) {    
#   #     if(length(levels(wdcMerged$tract_fac))<=1) {
#   #       X1 = rep(1,nrow(Z))
#   #     } else {
#   #       X1 = model.matrix(~ wdcMerged$tract_fac + 0)
#   #     }
#   #   } else {    
#   #     if(length(levels(wdcMerged$tract_fac))<=1) {
#   #       X1 = model.matrix(~ tw_group_vec + 0)
#   #     } else {
#   #       X1 = model.matrix(~ tw_group_vec + wdcMerged$tract_fac + 0)
#   #     }
#   #   }   
#   #     
#   #   Z_proj = Z - X1%*%solve(t(X1)%*%X1)%*%(t(X1)%*%Z)
#   #   Z <- cbind(Z_proj,X1)
#   #X1 <- cbind(X1,serv_lvl_covar)
#   
#   # X1 <- matrix(rep(1, nrow(wdcMerged)))
#   # Z <- matrix(rep(1, nrow(wdcMerged)))
#   
#   #   if(length(levels(tw_group_vec))<=1) {    
#   #     if(length(levels(tract_vec))<=1) {    
#   #       X1 = matrix(rep(1, nrow(wdcMerged)))
#   #     } else {
#   #       X1 = model.matrix(~ tract_vec +0)
#   #     }
#   #   } else {    
#   #     if(length(levels(tract_vec))<=1) {    
#   #       X1 = model.matrix(~ tw_group_vec)
#   #     } else {
#   #       X1 = model.matrix(~ tw_group_vec + tract_vec +0)
#   #     }
#   #   } 
#   wdcMerged$deltain <- deltain
#   wdcMerged$serv_lvl_covar <- serv_lvl_covar
#   formulat <- "deltain ~ serv_lvl_covar"
#   #   if(length(levels(wdcMerged$week_fac))>1) formulat <- paste0(formulat, " + week_fac")
#   #   if(length(levels(wdcMerged$tw_fac))>1) formulat <- paste0(formulat, " + tw_fac")
#   
#   #   fit <- plm(as.formula(formulat), effect="individual",
#   #            model="within",index=c("station_id_index_fac","tw_group"),data=wdcMerged)
#   #   st_fixed_effects <- fixef(fit)[wdcMerged$station_id_index]
#   # 
#   # #   tract_fixed_effects <- ave(st_fixed_effects, wdcMerged$tract, FUN=mean)
#   # #   st_fixed_effects <- st_fixed_effects - tract_fixed_effects
#   # #   grad_delta_theta_1ave <- ave(grad_delta_theta[,1],wdcMerged$tract, FUN=mean)
#   # #   grad_delta_theta_2ave <- ave(grad_delta_theta[,2],wdcMerged$tract, FUN=mean)
#   # #   grad_delta_theta[,1] <- grad_delta_theta[,1] - grad_delta_theta_1ave
#   # #   grad_delta_theta[,2] <- grad_delta_theta[,2] - grad_delta_theta_2ave
#   #   
#   #   #xi_norm = fit$residuals + st_fixed_effects - mean(st_fixed_effects)
#   #   xi_norm = fit$residuals
#   
#   # #this linear method, has similar gradient buyt still not very precise
#   #   fit <- lm(deltain ~ serv_lvl_covar + station_id_index_fac + 0,data=wdcMerged)
#   #   coef_fit <- coef(fit)
#   #   station_effects <- as.numeric(coef_fit[which(names(coef_fit) %like% "station_id_index")])
#   # station_effects <- station_effects[wdcMerged$station_id_index]
#   #   xi_norm = fit$residuals #+ station_effects - mean(station_effects)
#   # 
#   #   theta2 = coef(summary(fit))[1,1]
#   
#   #the linear regression, but with corrected gradient
#   X1 <- matrix(rep(1, nrow(wdcMerged)))
#   X1 <- cbind(X1,serv_lvl_covar)
#   theta2 = solve(t(X1)%*%X1)%*%(t(X1)%*%deltain)
#   xi_norm = deltain - X1%*%theta2
#   grad_xi_theta = grad_delta_theta - X1%*%solve(t(X1)%*%X1)%*%(t(X1)%*%grad_delta_theta)
#   
#   #the fixed effect regression, 
#   
#   
#   
#   if(mean(wdcMerged$obs_weight)>1) stop("do plm with weights")
#   
#   return(list("xi"=xi_norm,
#               "theta2"=theta2,
#               "grad_delta_theta"=grad_delta_theta))  
#   
# }
