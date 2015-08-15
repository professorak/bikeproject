

eval_covariates_delta_reg <- function(deltain,theta1,wdcMerged,points) {
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
    if(length(reg_terms)) {
      reg_formula <- paste(reg_formula, paste(reg_terms, collapse = " + "), sep = " + ")
    }
    
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
  list_covariates <- c("instr_serv_lvl","serv_lvl_neighbours","metro_den_on_1","metro_den_on_2","metro_den_off_1","metro_den_off_2"
                       ,"bus_den_1","bus_den_2","googleplaces_den_1","googleplaces_den_2","sto_nearby",
                       "metro_den_on_3","metro_den_on_4","metro_den_off_3","metro_den_off_4",
                       "log(metro_den_on_1+1)","log(metro_den_on_2+1)","log(metro_den_on_3+1)","log(metro_den_on_4+1)",
                       "log(metro_den_off_1+1)","log(metro_den_off_2+1)","log(metro_den_off_3+1)","log(metro_den_off_4+1)",
                       "bus_den_3","bus_den_4","log(bus_den_3+1)","log(bus_den_4+1)",
                       "googleplaces_den_3","googleplaces_den_4","log(googleplaces_den_3+1)","log(googleplaces_den_4+1)")
  if(length(unique(wdcMerged$tract)) > 1) {
    Zbase <- as.matrix(wdcMerged[stocked_list,c("census_density",list_covariates)])    
  } else {
    #not including census_density
    Zbase <- as.matrix(wdcMerged[stocked_list,list_covariates])
  }
  
  Z <- cbind(Xbase, Zbase)
  
  Z <- colNormalize(Z) #so that column sums are equal to 1.
  
  return(list("X"=X,              
              "Z"=Z))
}
