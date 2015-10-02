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
  
  Xbase <- matrix(1,length(stocked_list),1)/0.0001055609366/7.7874279
  colnames(Xbase) <- "Intercept"
  
  for(tract_i in tract_list) {
    for(tw_j in tw_list) {
      Xbase_col <- matrix(0,length(stocked_list),1)
      Xbase_col[which(wdcMerged$tract[stocked_list]==tract_i & 
                        wdcMerged$tw[stocked_list]==tw_j)] <- 1
      Xbase_col[which(wdcMerged$tract[stocked_list]==tract_i & 
                        wdcMerged$tw[stocked_list]==exclude_tw)] <- -1
      colnames(Xbase_col) <- paste0("tract_tw_fac",tract_i,"_",tw_j)
      Xbase_col <- Xbase_col/0.0001055609366/1.3628472
      Xbase <- cbind(Xbase,Xbase_col)
    }
  }
  
  #drop from Xbase columns which have no non-zero entry
  #rowsum is quite sure way to test
  if(length(which(colSums(abs(Xbase))==0))) {
    Xbase <- Xbase[,-which(colSums(Xbase)==0)]     
  }
  
  #add service level vector to convert Xbase to X
  X <- cbind(Xbase, serv_lvl=wdcMerged$serv_lvl[stocked_list]/0.0006466922293,
             serv_lvl_sq=(wdcMerged$serv_lvl[stocked_list]^2)/0.0006466922293/0.8857722)
  
  #add two service level instruments, local density attributes of a station, stockout indicator for nearby station
  #to get Z
  list_covariates <- c("instr_serv_lvl","serv_lvl_neighbours",
                       "metro_den_on_1","metro_den_on_2","metro_den_on_3","metro_den_on_4",
                       "metro_den_off_1","metro_den_off_2","metro_den_off_3","metro_den_off_4",
                       "bus_den_1","bus_den_2","bus_den_3","bus_den_4",
                       "googleplaces_food_den_1","googleplaces_food_den_2","googleplaces_food_den_3","googleplaces_food_den_4")
  
  if(length(unique(wdcMerged$tract)) > 1) {
    Zbase <- as.matrix(wdcMerged[stocked_list,c("census_density",list_covariates)])    
    #Zbase <- as.matrix(wdcMerged[stocked_list,list_covariates])
  } else {
    #not including census_density
    Zbase <- as.matrix(wdcMerged[stocked_list,list_covariates])
  }
  Zbase <- cbind(Zbase, "instr_serv_lvl_sq"=(wdcMerged$instr_serv_lvl[stocked_list]^2)/5.083953e-05)
  Zbase <- cbind(Zbase, "serv_lvl_neighbours_sq"=(wdcMerged$serv_lvl_neighbours[stocked_list]^2)/5.102270e-04) 
  
  if(length(which(colSums(Zbase)==0))) {
    Zbase <- Zbase[,-which(colSums(Zbase)==0)]     
  }
  
  Zbase[,which(colnames(Zbase)=="instr_serv_lvl")] <- Zbase[,which(colnames(Zbase)=="instr_serv_lvl")]/0.0001706792478
  Zbase[,which(colnames(Zbase)=="serv_lvl_neighbours")] <- Zbase[,which(colnames(Zbase)=="serv_lvl_neighbours")]/0.0006195223953            
  
  Z <- cbind(Xbase, Zbase)
  
  #  Z <- colNormalize(Z) #so that column sums are equal to 1.
  
  return(list("X"=X,              
              "Z"=Z))
}



