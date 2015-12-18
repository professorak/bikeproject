#finding the variables to inlcude from google places data
source("GooglePlacesVariablesSelectionFunctions.R")
source("temp_read_google_places_data.R")
getStationData <- function() {
  #Data generation:
  #constants
  tract_in <- c(1:10)
  number_local_points <- 14
  
  #1. generate station level data on demand
  station_demand <- getTotalStationDemand()
  ret <- getPointsData(station_demand)
  station_demand <- ret$station_demand
  points <- ret$points
  
  
  points <- subset(points, tract %in% tract_in)
  points$weight <- points$density
  points$density <- NULL #changing density to weight to track where density is being used
  points$weight[which(points$type==2)] <- points$weight[which(points$type==2)]/max(points$weight[which(points$type==2)])
  
  #merge points with places data
  places_data <- read_googleplaces_data()
  points <- points[which(points$type==1),]
  points <- points[order(points$lat, points$lon),]
  places_data <- places_data[order(places_data$lat, places_data$lon),]
  if(!identical(round(points$lat,4), round(places_data$lat,4))) stop("points lat dont match")
  if(!identical(round(points$lon,4), round(places_data$lon,4))) stop("points lon dont match")
  places_data_temp <- places_data
  places_data_temp$lat <- NULL
  places_data_temp$lon <- NULL
  places_colnames <- colnames(places_data_temp)
  points <- cbind(points, places_data_temp)
  
  #2. calculate nearby points characteristics from above google places data 
  #and previous data on metro stations ridership.
  #for each station consider the nearest 14 points. 14 because the number of grid points in 3 by 3 square is 14. (1+4+9)
  station_demand[,places_colnames] <- NA
  
  for(i in 1:nrow(station_demand)) {
    lat1 = station_demand$lat[i]
    lon1 = station_demand$lon[i]
    dis_v <- latlondistance(lat1,lon1,points$lat,points$lon)  
    order_dis_v <- order(dis_v)
    local_points <- order_dis_v[c(1:number_local_points)]
    station_demand[i,places_colnames] <- colSums(points[local_points,places_colnames])  
  }
  
  #distance to next station
  station_demand$dis_nearby_st <- NA
  #mean distance to nearest 5 station
  station_demand$dis_nearby_5_st <- NA
  
  for(i in 1:nrow(station_demand)) {
    lat1 = station_demand$lat[i]
    lon1 = station_demand$lon[i]
    dis_v <- latlondistance(lat1,lon1,station_demand$lat[-i],station_demand$lon[-i])  
    ordered_dis_v <- dis_v[order(dis_v)]
    station_demand$dis_nearby_st[i] <- ordered_dis_v[1]
    station_demand$dis_nearby_5_st[i] <- mean(ordered_dis_v[c(1:5)])
  }
  
  #assign tract
  wdcCenTrParsed <- readtractboundaryfile()
  tract  <- assign_points_tract(station_demand$lat,station_demand$lon,wdcCenTrParsed)
  X <- model.matrix(~factor(tract)+0)
  if(!identical(colnames(X), c("factor(tract)1", "factor(tract)2", "factor(tract)3", "factor(tract)4",
                               "factor(tract)5",  "factor(tract)6", "factor(tract)7", "factor(tract)8",
                               "factor(tract)9",  "factor(tract)10"))) stop("tract matrix colnames dont match")
  colnames(X) <- c("tract_1", "tract_2", "tract_3", "tract_4", "tract_5", "tract_6",
                   "tract_7", "tract_8","tract_9",  "tract_10")
  station_demand <- cbind(station_demand, X)
  
  
  return(list(station_demand=station_demand,
              places_colnames=places_colnames))
}

ret <- getStationData()
station_demand <- ret$station_demand
station_demand$dis_nearby_st_sq <- station_demand$dis_nearby_st^2
station_demand$dis_nearby_5_st_sq <- station_demand$dis_nearby_5_st^2
places_colnames <- ret$places_colnames
additional_controls <- c("dis_nearby_st","dis_nearby_st_sq","dis_nearby_5_st",
                         "dis_nearby_5_st_sq","tract_1",
                         "tract_2", "tract_3", "tract_4", "tract_5",
                         "tract_6", "tract_7", "tract_8", "tract_9", "tract_10")

#do lasso regression on station_demand with 

require("glmnet")
grid=exp(seq(10,-10,length=100))
x <- as.matrix(station_demand[,c(places_colnames,additional_controls)])

#If alpha=0 then a ridge regression model is fit, and if alpha=1 then a lasso model is fit
#lasso.mod=glmnet(x[train ,],y[train],alpha=1,lambda=grid)
cv.out=cv.glmnet(x,station_demand$out_dem,alpha=1)
plot(cv.out)
bestlam=cv.out$lambda.min

lambda_val <- 0.003 
out=glmnet(x,station_demand$out_dem,alpha=1,lambda=lambda_val)
lasso.coef=predict(out,type="coefficients",s=lambda_val)
lasso.coef[which(as.numeric((lasso.coef))!=0),]




