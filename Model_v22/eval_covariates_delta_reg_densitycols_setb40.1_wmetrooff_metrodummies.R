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
  
#   #metro_den_on_1_2 for 0 to 300mts metro impact
#   wdcMerged$metro_den_on_1_2 <- wdcMerged$metro_den_on_1 + wdcMerged$metro_den_on_2
#   wdcMerged$touristlocs_den_1_2 <- wdcMerged$touristlocs_den_1 + wdcMerged$touristlocs_den_2
#   wdcMerged$googleplaces_food_den_1_2 <- wdcMerged$googleplaces_food_den_1 + wdcMerged$googleplaces_food_den_2
#   wdcMerged$googleplaces_government_den_1_2 <- wdcMerged$googleplaces_government_den_1 + wdcMerged$googleplaces_government_den_2
#   wdcMerged$googleplaces_lodging_den_1_2 <- wdcMerged$googleplaces_lodging_den_1 + wdcMerged$googleplaces_lodging_den_2
#   wdcMerged$googleplaces_museum_den_1_2 <- wdcMerged$googleplaces_museum_den_1 + wdcMerged$googleplaces_museum_den_2
#   wdcMerged$googleplaces_movie_theater_den_1_2 <- wdcMerged$googleplaces_movie_theater_den_1 + wdcMerged$googleplaces_movie_theater_den_2
  metro_vars <- c("metro_den_on_1","metro_den_on_2",
    "metro_den_off_1","metro_den_off_2") #not including _3 and _4 as they are quite correlated with intercept.  
  metro_dummies_df <- wdcMerged[,metro_vars]
  metro_dummies_df <- metro_dummies_df>0
  colnames(metro_dummies_df) <- paste("dummy_",metro_vars, sep="")
  wdcMerged <- cbind(wdcMerged, metro_dummies_df)


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
     
  #distance from city center
  city_center <- c(48.857453, 2.351759)
  dis_center <- latlondistance(wdcMerged$lat[stocked_list], wdcMerged$lon[stocked_list], city_center[1], city_center[2])
  Xbase <- cbind(Xbase, "dis_center"=dis_center)
 
  #drop from Xbase columns which have no non-zero entry
  #rowsum is quite sure way to test
  if(length(which(colSums(abs(Xbase))==0))) {
    Xbase <- Xbase[,-which(colSums(Xbase)==0)]     
  }
  X <- Xbase
  
  
  #add two service level instruments, local density attributes of a station, stockout indicator for nearby station
  #to get Z
  # list_covariates <- c("metro_den_on_1_2","metro_den_on_3","metro_den_on_4",
  #                      "googleplaces_food_den_1_2","googleplaces_food_den_3","googleplaces_food_den_4",                       
  #                      "googleplaces_government_den_1_2","googleplaces_government_den_3","googleplaces_government_den_4",
  #                      "catchment_area_step1","catchment_area_step2","catchment_area_step3","catchment_area_step4",
  #                      "touristlocs_den_1_2","touristlocs_den_3","touristlocs_den_4",
  #                      "googleplaces_lodging_den_1_2","googleplaces_lodging_den_3","googleplaces_lodging_den_4",
  #                      "googleplaces_museum_den_1_2","googleplaces_museum_den_3","googleplaces_museum_den_4",
  #                      "googleplaces_movie_theater_den_1_2","googleplaces_movie_theater_den_3","googleplaces_movie_theater_den_4")
  list_covariates <- c("metro_den_on_1","metro_den_on_2","metro_den_on_3","metro_den_on_4",
                     "metro_den_off_1","metro_den_off_2","metro_den_off_3","metro_den_off_4",
                     paste("dummy_",metro_vars, sep=""),
                     "touristlocs_den_1","touristlocs_den_2","touristlocs_den_3","touristlocs_den_4")

  if(length(unique(wdcMerged$tract)) > 1) {
    Zbase <- as.matrix(wdcMerged[stocked_list,c("census_density",list_covariates)])    
    #Zbase <- as.matrix(wdcMerged[stocked_list,list_covariates])
  } else {
    #not including census_density
    Zbase <- as.matrix(wdcMerged[stocked_list,list_covariates])
  }
    
  if(length(which(colSums(Zbase)==0))) {
    Zbase <- Zbase[,-which(colSums(Zbase)==0)]     
  }
    
  Z <- cbind(Xbase, Zbase)
  print(colnames(Z))
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



