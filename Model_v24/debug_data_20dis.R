source("constants_mw6.R")
source("debug_functions.R")
#model with averaged station demand, no service level. Trying to get distance
#coefficient right.

get_points_density_grad_places_count_col <- NULL
source("eval_obj_GMM_model5.R")
source("GetDataGooglePlaces_20dis_aggmonth_averaged_pr6.R")
v0_vec <- generate_v0(1000)
v0_vec_weights <- rep(1, length(v0_vec))

dim(wdcMerged)



wdcCenTrParsed <- readtractboundaryfile_presentation()

symbols(wdcCenTrParsed$lon, wdcCenTrParsed$lat, circles=rep(1,nrow(wdcCenTrParsed)),
        inches=0.08, fg="black", 
        cex.axis=0.7,bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
lon_vec <- c(wdcCenTrParsed$lon, wdcMerged$lon)
lat_vec <- c(wdcCenTrParsed$lat, wdcMerged$lat)
rad_vec <- c(rep(0.5,nrow(wdcCenTrParsed)), rep(1,nrow(wdcMerged)))

symbols(lon_vec, lat_vec, circles=rad_vec,
        inches=0.02, fg="black", 
        cex.axis=0.7,bty="n",xaxt="n",yaxt="n",xlab="",ylab="")


#testing if it is easy to implement ray casting algo
#take a point, draw a ray and see how many times it intersects boundary of each of the tract

#taking a point in 13 arr: 48.8355325  ,2.3525275
test_point_lat <- 48.8355325
test_point_lon <- 2.3525275

#steps
#1. find closest point on boundaries to a test_point
#2. draw a ray of long enough length from this test_point towards this point.
#3. count how many times does this ray intersect. Drawling ray will be checking 
#   on discrete point of this ray how close it is to the boundary points and define as intersection if it is 
#   small enough. (the distance between boundary points is 50mts)

#1
# dis_v <- latlondistance(test_point_lat, test_point_lon, wdcCenTrParsed$lat, wdcCenTrParsed$lon)
# dis_v[which.min(dis_v)]

dis_v <- latlondistance(test_point_lat, test_point_lon, test_point_lat+0.001, test_point_lon)
dis_v[which.min(dis_v)]

#draw a line which long enough. For now try 4kms
#using discrete points at roughly 50mts.
discrete_point_size <- 0.010

# unit_step <- c(wdcCenTrParsed$lat[which.min(dis_v)]-test_point_lat, wdcCenTrParsed$lon[which.min(dis_v)] - test_point_lon) * 
#   discrete_point_size/dis_v[which.min(dis_v)]
unit_step <-  c(0.001,0) * 
   discrete_point_size/dis_v[which.min(dis_v)]

#2. draw a ray of long enough length from this test_point towards this point.
#number of points would be 4kms/50mts
ray <- matrix(c(0:(10/discrete_point_size)),,1) %*% matrix(unit_step,1,)
ray[,1] <- ray[,1] + test_point_lat
ray[,2] <- ray[,2] + test_point_lon
ray <- as.data.frame(ray)
colnames(ray) <- c("lat","lon")

lon_vec <- c(wdcCenTrParsed$lon, ray[,2])
lat_vec <- c(wdcCenTrParsed$lat, ray[,1])
rad_vec <- c(rep(0.5,nrow(wdcCenTrParsed)), rep(1,nrow(ray)))

# symbols(lon_vec, lat_vec, circles=rad_vec,
#         inches=0.02, fg="black", 
#         cex.axis=0.7,bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
# 

#3. count how many times does this ray intersect. Drawling ray will be checking 
#   on discrete point of this ray how close it is to the boundary points and define as intersection if it is 
#   small enough. (the distance between boundary points is 50mts)
count_intersections_vec <- c()
for(tract_i in c(1:20)) {
  wdcCenTrParsed_tract <- subset(wdcCenTrParsed, tract==tract_i)
  min_dis_vec <- rep(0,nrow(ray))  
  for(j in c(1:nrow(ray))) {
    dis_v <- latlondistance(ray$lat[j],ray$lon[j],
      wdcCenTrParsed_tract$lat, wdcCenTrParsed_tract$lon)        
    min_dis_vec[j] <- dis_v[which.min(dis_v)]
  }
  #if there is intersection (<35mts), then jump 200 mts and then start checking.
  count_intersections <- 0
  j <- 1
  while(j<=nrow(ray)) {
    if(min_dis_vec[j] <= 0.035) {
      count_intersections <- count_intersections + 1
      j <- j + (0.200/discrete_point_size)
    } else {
      j <- j + 1
    }
  }
  
  count_intersections_vec <- c(count_intersections_vec,count_intersections)  
}


#find that there are two districts with 1 intersections. such situations will
#happen since there are cases of just going along boundary etc. best solution is to
#draw many lines (say 5 to some 5 centers of paris) and take intersection.

source("tract_assignment_functions.R")
# intersection_gap <- 0.020 # 20mts
# skip_upon_intersection <- 0.200 #200mts
# discrete_point_size <- 0.010 #distance between consecutive points of ray and also the boundary of tracts

discrete_point_size <- 0.025 #distance between consecutive points of ray and also the boundary of tracts
intersection_gap <- 0.025 # 20mts
skip_upon_intersection <- 0.200 #200mts

data_point <- c(48.835532 , 2.352527)
wdcCenTrParsed <- readtractboundaryfile_finegrain(discrete_point_size)
system.time({
  print(classify_to_tract(data_point,wdcCenTrParsed,discrete_point_size,
                          intersection_gap, skip_upon_intersection))  
})





