#tract assignment functions

classify_to_tract <- function(data_point,wdcCenTrParsed, discrete_point_size,
                              intersection_gap, skip_upon_intersection) {
  #find directions to search in of size discrete_point_size
  min_lat_wdcCenTrParsed <- min(wdcCenTrParsed$lat)
  max_lat_wdcCenTrParsed <- max(wdcCenTrParsed$lat)
  min_lon_wdcCenTrParsed <- min(wdcCenTrParsed$lon)
  max_lon_wdcCenTrParsed <- max(wdcCenTrParsed$lon)

  if(data_point[1]>max_lat_wdcCenTrParsed | data_point[1]<min_lat_wdcCenTrParsed |
       data_point[2]>max_lon_wdcCenTrParsed | data_point[2]<min_lon_wdcCenTrParsed) {
    return(-1)
  }
  
  direction_list <- list(c(0.001,0),c(-0.001,0),c(0,0.001),c(0,-0.001))
  count_intersection_tracts_df <- c()
  for(direction in direction_list) {
    dis_v <- latlondistance(data_point[1], data_point[2], data_point[1]+direction[1], data_point[2]+direction[2])
    unit_step <- direction * discrete_point_size/dis_v
    
    ray <- draw_ray(data_point,unit_step,20,discrete_point_size) #20kms ray    
    # cut off the ray once it reaches rectangular boundary of wdcCenTrParsed
    ray <- ray[-which(ray$lat>max_lat_wdcCenTrParsed | ray$lat<min_lat_wdcCenTrParsed |
            ray$lon>max_lon_wdcCenTrParsed | ray$lon<min_lon_wdcCenTrParsed),]
        
    count_intersection_tracts_df <- rbind(count_intersection_tracts_df,
      count_intersection_tracts(ray, wdcCenTrParsed, discrete_point_size, 
                                intersection_gap, skip_upon_intersection))
  }
  
  count_intersection_tracts_odd_df <- count_intersection_tracts_df %% 2
  #taking prod columinwise
  count_intersection_tracts_odd <- apply(count_intersection_tracts_odd_df,2,FUN=prod)
  if(length(which(count_intersection_tracts_odd==1))) {
    return(which(count_intersection_tracts_odd==1))
  } else {
    return(-1)
  }
} 


draw_ray <- function(data_point,unit_step,length_ray,discrete_point_size) {
  #2. draw a ray of long enough length from this test_point towards this point.
  #number of points would be 4kms/50mts
  ray <- matrix(c(0:(length_ray/discrete_point_size)),,1) %*% matrix(unit_step,1,)
  ray[,1] <- ray[,1] + data_point[1]
  ray[,2] <- ray[,2] + data_point[2]
  ray <- as.data.frame(ray)
  colnames(ray) <- c("lat","lon")
  return(ray)
}

#count number of intersections with each of the tract
count_intersection_tracts <- function(ray, wdcCenTrParsed, discrete_point_size,
                                      intersection_gap, skip_upon_intersection) {
  #3. count how many times does this ray intersect. Drawling ray will be checking 
  #   on discrete point of this ray how close it is to the boundary points and define as intersection if it is 
  #   small enough. (the distance between boundary points is 10mts)

  
  count_intersections_vec <- c()
  for(tract_i in c(1:20)) {
    wdcCenTrParsed_tract <- subset(wdcCenTrParsed, tract==tract_i)
    min_dis_vec <- rep(0,nrow(ray))  
    for(j in c(1:nrow(ray))) {
      dis_v <- latlondistance(ray$lat[j],ray$lon[j],
                              wdcCenTrParsed_tract$lat, wdcCenTrParsed_tract$lon)        
      min_dis_vec[j] <- dis_v[which.min(dis_v)]
    }
    #if there is intersection (<intersection_gap ), then jump skip_upon_intersection and then start checking.
    count_intersections <- 0
    j <- 1
    while(j<=nrow(ray)) {
      if(min_dis_vec[j] <= intersection_gap) {
        count_intersections <- count_intersections + 1
        j <- j + (skip_upon_intersection/discrete_point_size)
      } else {
        j <- j + 1
      }
    }
    
    count_intersections_vec <- c(count_intersections_vec,count_intersections)  
  }
  return(count_intersections_vec)
}



