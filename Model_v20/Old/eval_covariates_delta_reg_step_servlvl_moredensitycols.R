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
  tract_list <- sort(unique(wdcMerged$tract)) 
  tw_list <- levels(wdcMerged$tw_fac)
  exclude_tw <- tw_list[1]
  tw_list <- tw_list[-1]
  
  Xbase <- matrix(1,length(stocked_list),1)
  colnames(Xbase) <- "Intercept"
  
  for(tract_i in tract_list) {
    for(tw_j in tw_list) {
      Xbase_col <- matrix(0,length(stocked_list),1)
      Xbase_col[which(wdcMerged$tract[stocked_list]==tract_i & 
                        wdcMerged$tw[stocked_list]==tw_j)] <- 1
      Xbase_col[which(wdcMerged$tract[stocked_list]==tract_i & 
                        wdcMerged$tw[stocked_list]==exclude_tw)] <- -1
      colnames(Xbase_col) <- paste0("tract_tw_fac",tract_i,"_",tw_j)
      Xbase_col <- Xbase_col
      Xbase <- cbind(Xbase,Xbase_col)
    }
  }
  
  #drop from Xbase columns which have no non-zero entry
  #rowsum is quite sure way to test
  if(length(which(colSums(abs(Xbase))==0))) {
    Xbase <- Xbase[,-which(colSums(Xbase)==0)]     
  }
  
  #add service level vector to convert Xbase to X
  serv_step1 <- c(0,0.5)
  serv_step2 <- c(0.5,0.8)
  serv_step3 <- c(0.8,1.1)  
  
  serv_lvl1 = trim_serv(wdcMerged$serv_lvl[stocked_list],serv_step1)
  serv_lvl2 = trim_serv(wdcMerged$serv_lvl[stocked_list],serv_step2)
  serv_lvl3 = trim_serv(wdcMerged$serv_lvl[stocked_list],serv_step3)
  
  X <- cbind(Xbase, serv_lvl1=serv_lvl1,
             serv_lvl2=serv_lvl2,
             serv_lvl3=serv_lvl3)
  
  #add two service level instruments, local density attributes of a station, stockout indicator for nearby station
  #to get Z
  list_covariates <- c("metro_den_on_1","metro_den_on_2","metro_den_on_3","metro_den_on_4",
                       "metro_den_off_1","metro_den_off_2","metro_den_off_3","metro_den_off_4",
                       "bus_den_1","bus_den_2","bus_den_3","bus_den_4",
                       "googleplaces_food_den_1","googleplaces_food_den_2","googleplaces_food_den_3","googleplaces_food_den_4",
                       "googleplaces_grocery_den_1","googleplaces_grocery_den_2","googleplaces_grocery_den_3","googleplaces_grocery_den_4",
                       "googleplaces_government_den_1","googleplaces_government_den_2","googleplaces_government_den_3","googleplaces_government_den_4"
                       )
  
  if(length(unique(wdcMerged$tract)) > 1) {
    Zbase <- as.matrix(wdcMerged[stocked_list,c("census_density",list_covariates)])    
    #Zbase <- as.matrix(wdcMerged[stocked_list,list_covariates])
  } else {
    #not including census_density
    Zbase <- as.matrix(wdcMerged[stocked_list,list_covariates])
  }
  Zserv_lvl <- cbind(trim_serv_frac(wdcMerged$instr_serv_lvl[stocked_list],serv_step1),
    trim_serv_frac(wdcMerged$instr_serv_lvl[stocked_list],serv_step2),
    trim_serv_frac(wdcMerged$instr_serv_lvl[stocked_list],serv_step3),
    trim_serv(wdcMerged$serv_lvl_neighbours[stocked_list],serv_step1),
    trim_serv(wdcMerged$serv_lvl_neighbours[stocked_list],serv_step2),
    trim_serv(wdcMerged$serv_lvl_neighbours[stocked_list],serv_step3)
    )
  colnames(Zserv_lvl) <- c("instr_serv_lvl1","instr_serv_lvl2","instr_serv_lvl3",
          "serv_lvl_neighbours1","serv_lvl_neighbours2","serv_lvl_neighbours3")
  Zbase <- cbind(Zbase, Zserv_lvl)
  
  if(length(which(colSums(Zbase)==0))) {
    Zbase <- Zbase[,-which(colSums(Zbase)==0)]     
  }
  
  
  Z <- cbind(Xbase, Zbase)
  
  #  Z <- colNormalize(Z) #so that column sums are equal to 1.
      
  stocked_list <- which(wdcMerged$stocked_out==FALSE)
  X_sqweighted <- (X * sqrt(wdcMerged$obs_weight[stocked_list]))/sum(wdcMerged$obs_weight[stocked_list])
  scalecols_X <<- sqrt(diag(t(X_sqweighted) %*% X_sqweighted))
  X <- scalecols(X,scalecols_X)
  
  Z_sqweighted <- (Z * sqrt(wdcMerged$obs_weight[stocked_list]))/sum(wdcMerged$obs_weight[stocked_list])
  scalecols_Z <<- sqrt(diag(t(Z_sqweighted) %*% Z_sqweighted))
  Z <- scalecols(Z,scalecols_Z)
  
  return(list("X"=X,              
              "Z"=Z))
}



