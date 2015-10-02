#tract fe, no service level, 

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
  
  #metro_den_on_1_2 for 0 to 300mts metro impact
  wdcMerged$metro_den_on_1_2 <- wdcMerged$metro_den_on_1 + wdcMerged$metro_den_on_2
  wdcMerged$touristlocs_den_1_2 <- wdcMerged$touristlocs_den_1 + wdcMerged$touristlocs_den_2
  wdcMerged$googleplaces_food_den_1_2 <- wdcMerged$googleplaces_food_den_1 + wdcMerged$googleplaces_food_den_2
  wdcMerged$googleplaces_government_den_1_2 <- wdcMerged$googleplaces_government_den_1 + wdcMerged$googleplaces_government_den_2
  wdcMerged$googleplaces_lodging_den_1_2 <- wdcMerged$googleplaces_lodging_den_1 + wdcMerged$googleplaces_lodging_den_2
  wdcMerged$googleplaces_museum_den_1_2 <- wdcMerged$googleplaces_museum_den_1 + wdcMerged$googleplaces_museum_den_2
  wdcMerged$googleplaces_movie_theater_den_1_2 <- wdcMerged$googleplaces_movie_theater_den_1 + wdcMerged$googleplaces_movie_theater_den_2

  
  #generating the common attributes of X and Z. Calling it Xbase  
  reg_formula <- "~tract_tw_fac + 0"

  Xbase <- model.matrix(as.formula(reg_formula), data=wdcMerged[stocked_list,])    
  
  
  #drop from Xbase columns which have no non-zero entry
  #rowsum is quite sure way to test
  if(length(which(colSums(abs(Xbase))==0))) {
    Xbase <- Xbase[,-which(colSums(Xbase)==0)]     
  }
  X <- Xbase
  
  #add two service level instruments, local density attributes of a station, stockout indicator for nearby station
  #to get Z
  list_covariates <- c("metro_den_on_1_2",
                       "googleplaces_food_den_1_2",
                       "googleplaces_government_den_1_2",
                       "catchment_area_step1","catchment_area_step2",
                       "touristlocs_den_1_2",
                       "googleplaces_lodging_den_1_2",
                       "googleplaces_museum_den_1_2",
                       "googleplaces_movie_theater_den_1_2")

  print(list_covariates)
  #not including census_density
  Zbase <- as.matrix(wdcMerged[stocked_list,list_covariates])
  
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



