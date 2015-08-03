weighing_GMM_mat <<- NULL
eval_error_xi_model4 <- function(theta1,deltain,wdcMerged,points) {
  
  if(!identical(order(user_serv_lvl$st_tw_index), c(1:nrow(user_serv_lvl)))) stop("user_serv_lvl not sorted")  
  stocked_list <- which(wdcMerged$stocked_out==FALSE)
  
  #compute attributes matrix X and instruments matrix Z
  #list of X attributes
  #1s, sevice level, diXw fixed effects, month fixed effects, weather fixed effects 
  #list of Z attributes
  #1s, diXw fixed effects, month fixed effects, weather fixed effects, two service level instruments, local density attributes of a station, 
  #stockout indicator for nearby station
  #of the Z's need to get --local density attributes of a station, and stockout indicator for nearby station, while generating data, 
  #one entry for each data row.
  
  #generating the common attributes of X and Z. Calling it Xbase
  #1s, diXw fixed effects, month fixed effects, weather fixed effects
  if(length(levels(wdcMerged$tract_tw_fac))>1 & length(levels(wdcMerged$week_fac))>1) {
    reg_formula <- "~week_fac + tract_tw_fac"
    if(length(levels(wdcMerged$Conditions))>1) paste(reg_formula, "+ Conditions")
    length_levels <- c(length(levels(wdcMerged$Conditions)), length(levels(wdcMerged$Temperature_group)),
                       length(levels(wdcMerged$Humidity_high)), length(levels(wdcMerged$Wind.Speed_high)))
    reg_terms <- c("Conditions", "Temperature_group", "Humidity_high", "Wind.Speed_high")[which(length_levels>1)]
    reg_formula <- paste(reg_formula, paste(reg_terms, collapse = " + "), sep = " + ")
    
    Xbase <- model.matrix(as.formula(reg_formula), data=wdcMerged[stocked_list,])    

    if(length(which(colnames(Xbase) %like% "Conditions")) != length(levels(wdcMerged$Conditions))-1) {
      stop("in eval_error_xi_sl_model4 one conditions number not automatically removed")      
    }
    if(length(which(colnames(Xbase) %like% "Temperature_group")) != length(levels(wdcMerged$Temperature_group))-1) {
      stop("in eval_error_xi_sl_model4 one Temperature_group number not automatically removed")      
    }
    if(length(which(colnames(Xbase) %like% "Humidity_high")) != length(levels(wdcMerged$Humidity_high))-1) {
      stop("in eval_error_xi_sl_model4 one Humidity_high number not automatically removed")      
    }
    if(length(which(colnames(Xbase) %like% "Wind.Speed_high")) != length(levels(wdcMerged$Wind.Speed_high))-1) {
      stop("in eval_error_xi_sl_model4 one Wind.Speed_high number not automatically removed")      
    }
    if(length(which(colnames(Xbase) %like% "tract_tw_fac")) != length(levels(wdcMerged$tract_tw_fac))-1) {
      stop("in eval_error_xi_sl_model4 one tract_tw_fac number not automatically removed")      
    }
  } else if (length(levels(wdcMerged$tract_tw_fac))>1) {
    stop("error in eval_error_xi")
  } else if (length(levels(wdcMerged$week_fac))>1) {
    stop("error in eval_error_xi")
  } else {
    stop("error in eval_error_xi")
  }  
  #drop from Xbase columns which have no non-zero entry
  #rowsum is quite sure way to test
  if(length(which(colSums(Xbase)==0))) {
    Xbase <- Xbase[,-which(colSums(Xbase)==0)]     
  }
  
  #add service level vector to convert Xbase to X
  X <- cbind(Xbase, serv_lvl=wdcMerged$serv_lvl[stocked_list])

  #add two service level instruments, local density attributes of a station, stockout indicator for nearby station
  #to get Z
  Zbase <- as.matrix(wdcMerged[stocked_list,c("instr_serv_lvl","serv_lvl_neighbours","census_density","metro_den_on_1","metro_den_on_2","metro_den_off_1","metro_den_off_2"
                        ,"bus_den_1","bus_den_2","googleplaces_den_1","googleplaces_den_2","sto_nearby")])
  Z <- cbind(Xbase, Zbase)
  
  Z <- colNormalize(Z) #so that column sums are equal to 1.
  
  #theta_2 is given by, Inv(X'.W'.Z.A_N.Z'.W.X).(X'.W'.Z.A_N.Z'.W.delta)
  #W above can be replaced by W/|W|
  
  #Z_weighted equal to W'.Z
  Z_weighted <- (Z * (wdcMerged$obs_weight[stocked_list]))/sum(wdcMerged$obs_weight[stocked_list])
  if(is.null(weighing_GMM_mat)) {
    A_N <- diag(ncol(Z))    
  } else {
    A_N <- weighing_GMM_mat
  }

  theta_2 <- solve(t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted) %*% X) %*%
                (t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted) %*% deltain)
  
  xi <- deltain - (X %*%  theta_2) 
  
  obj_moments <- t(Z_weighted) %*% xi
  obj <- t(obj_moments) %*% A_N %*% obj_moments
  grad_obj_wrt_delta <- 2*c((Z_weighted) %*% A_N %*% obj_moments)
  
  return(list("obj"=obj,
              "theta2"=theta_2,
              "theta3"=NULL,
              "grad_obj_wrt_delta"=grad_obj_wrt_delta,
              "Z"=Z,
              "weights"=wdcMerged$obs_weight[stocked_list],
              "xi"=xi))
}

eval_obj_GMM_model4_obj <- function(params, wdcMerged, points, length_theta) {    
  #obj_eval <<- 1
  theta <- params[c(1:length_theta)]
  deltain_stin <- params[-c(1:length_theta)]
  print("objective comput")
  return(eval_obj_GMM_model4(theta, deltain_stin, wdcMerged, points)$objective)
}  
eval_obj_GMM_model4_grad <- function(params, wdcMerged, points, length_theta) {    
  #obj_eval <<- 0
  theta <- params[c(1:length_theta)]
  deltain_stin <- params[-c(1:length_theta)]  
  print("gradint comput")
  return(eval_obj_GMM_model4(theta, deltain_stin, wdcMerged, points)$gradient)
} 

eval_obj_GMM_model4 <- function(theta, deltain, wdcMerged, points) {  
  ptm <- proc.time()
  theta1 <- c(theta[1],0,theta[-1])  
  
  print("theta1:")
  print(paste(theta1))
  
  tryCatch({
    ret <- eval_error_xi_model4 (theta1, deltain, wdcMerged,points)
    gmm_objective <- ret$obj
    gmm_grad_delta <- ret$grad_obj_wrt_delta
    gmm_grad_theta <- rep(0,length(theta))
    gmm_gradient <- c(gmm_grad_theta,gmm_grad_delta)
    
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
  print(c(gmm_gradient[1:length(theta)],mean(gmm_gradient[-c(1:length(theta))])))
  print("eval_obj time:")
  print(proc.time()-ptm)
  return( list( "objective" = gmm_objective,
                "gradient"  = gmm_gradient) )
  
}



