get_points_density <- function(points_in, theta1, tw_in) {
  if(tw_in == 0) {
    density <- ((points_in$type==1)*theta1[density_ridership_col] + 
                  (points_in$type==2)*theta1[density_metro_evening_col]) *points_in$weight   + 
      (points_in$type==1)*theta1[density_intercept_col]    
  } else{
    density <- ((points_in$type==1)*theta1[density_ridership_col] + 
                  (points_in$type==2)*theta1[density_metro_col]) *points_in$weight   + 
      (points_in$type==1)*theta1[density_intercept_col]    
  }
  density <- density + 
    (as.matrix(points_in[,density_google_places_points_cols]) %*% 
    as.matrix(theta1[density_google_places_theta_cols])) +
    ((points_in$type==3)*theta1[density_add_den_cols[1]] + #bus 
    (points_in$type==4)*theta1[density_add_den_cols[2]])*points_in$weight #tourist locations    
  
  return(density)
}

get_points_density_grad_places_count_col <- function(points_in, tw_in) {
  density <- points_in[,density_google_places_points_cols]
  return(density)
}


