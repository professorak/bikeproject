latlondistance <- function(lat1, lon1, lat2, lon2) {
  #Equirectangular approximation
  lat1 <- lat1*3.14/180
  lon1 <- lon1*3.14/180
  lat2 <- lat2*3.14/180
  lon2 <- lon2*3.14/180
  R <- 6371
  x <- (lon2-lon1) * cos((lat1+lat2)/2)
  y <- (lat2-lat1)
  d <- sqrt(x*x + y*y) * R
  return(d)
}

# so_state <- function(so_vec) {
#   #convert binary so_vec of stockouts to decimal  
#   binary_vec <- c(0:(length(so_vec)-1))
#   binary_vec <- exp(binary_vec*log(2))
#   #return (as.character(sum(binary_vec*so_vec)))
#   return (sum(binary_vec*so_vec))
# }

gen_sto_state_local <- function(so_vec,local_stations) {
  #convert local_stations to a vector
  #extract corresponding elements from so_vec and convert them to bin and return
  #both arguments could be vectors for a single time observations.
  #have to apply entire vector so_vec for each local_station obs
  so_state_local_vec <- rep(NA,length(local_stations))
  for(i in 1:length(local_stations)) {
    indexvec <- binaryvec_to_indexes(local_stations[i])
    binary_vec <- c(0:(length(so_vec[indexvec])-1))
    binary_vec <- exp(binary_vec*log(2))
    so_state_local_vec[i] <- sum(binary_vec*so_vec[indexvec])
  }  
  return(so_state_local_vec)
}

gen_sto_state_local_char <- function(so_vec,local_stations) {
  
  so_state_local_vec <- rep(NA,length(local_stations))
  #for each of the local_stations entry, generate the projection from so_vec and create a char from it
  for(i in 1:length(local_stations)) {
    indexvec <- as.integer(strsplit(as.character(local_stations[i]),split="_")[[1]])
    so_state_local_vec[i] <- paste(as.numeric(so_vec[indexvec]),collapse="_")
  }
  return(so_state_local_vec)
}  

gen_bike_state_local_char <- function(bike_vec,local_stations) {
  
  bike_state_local_vec <- rep(NA,length(local_stations))
  #for each of the local_stations entry, generate the projection from so_vec and create a char from it
  for(i in 1:length(local_stations)) {
    indexvec <- as.integer(strsplit(as.character(local_stations[i]),split="_")[[1]])
    bike_state_local_vec[i] <- paste(as.numeric(bike_vec[indexvec]),collapse="_")
  }
  return(bike_state_local_vec)
}  

indexes_to_binaryvec <- function(indexes_vec) {
  binary_vec <- 2^indexes_vec
  return(sum(binary_vec))  
}

binaryvec_to_indexes <- function(bin_number) {
  return(which(rev(digitsBase(bin_number,base=2)==1))-1)
}


regout<-function(x){
  res<-c(paste(as.character(summary(x)$call),collapse=" "),
         x$coefficients[1],
         x$coefficients[2],
         length(x$model),
         summary(x)$coefficients[2,2],
         summary(x)$r.squared,
         summary(x)$adj.r.squared,
         summary(x)$fstatistic,
         pf(summary(x)$fstatistic[1],summary(x)$fstatistic[2],summary(x)$fstatistic[3],lower.tail=FALSE),
         summary(x$residuals))
  names(res)<-c("call","intercept","slope","n","slope.SE","r.squared","Adj. r.squared",
                "F-statistic","numdf","dendf","p.value",
                "Residuals-Min","1st Qu.","Median","Mean","3rd Qu.","Max.")
  return(res)
}  


splitchar <- function(str,splitchar="_"){
  return(as.numeric(strsplit(as.character(str)
                                  , split=splitchar)[[1]]))
}

splitchar_sum <- function(str,splitchar="_"){
  return(sum(as.numeric(strsplit(as.character(str)
                             , split=splitchar)[[1]])))
}

low_bikes <- function(str,splitchar="_"){
  a= as.numeric(strsplit(as.character(str), split=splitchar)[[1]])
  a = a[!is.na(a)]  
  return(length(which(a<=low_bikes_threshold))>=1)  
}


gen_serv_lvl_instr_local_char <- function(so_vec,local_stations) {  
  so_state_local_vec <- rep(NA,length(local_stations))
  #for each of the local_stations entry, generate the projection from so_vec and create a char from it
  for(i in 1:length(local_stations)) {
    indexvec <- as.integer(strsplit(as.character(local_stations[i]),split="_")[[1]])
    so_state_local_vec[i] <- mean(as.numeric(so_vec[indexvec]))
  }
  return(so_state_local_vec)
}  

gen_serv_lvl_instr_neigh <- function(so_vec,local_stations) {  
  so_state_local_vec <- rep(NA,length(local_stations))
  #for each of the local_stations entry, generate the projection from so_vec and create a char from it
  for(i in 1:length(local_stations)) {
    indexvec <- as.integer(strsplit(as.character(local_stations[i]),split="_")[[1]])
    indexvec <- indexvec[-which(indexvec==i)] #removing self statiod id to keep only neighbours
    so_state_local_vec[i] <- mean(as.numeric(so_vec[indexvec]))
  }
  return(so_state_local_vec)
} 

stid_in_localstations <- function(stid, localstations) {
  retvec <- c()
  for(localstations_str in localstations) {
    retvec <- c(retvec, stid %in% splitchar(localstations_str))  
  }
  return(retvec)
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)


getfullfromsparsematrix <- function (x, data = NULL, ncol = NULL) 
{
  #adopted from print.sparseness in ipoptr package
  stopifnot(is.list(x))
  if (is.null(ncol)) {
    ncol <- max(unlist(x))
  }
  p <- matrix(0, nrow = length(x), ncol)
  #names(p) <- 1:ncol
  cnt = 1
  for (row in 1:length(x)) {
    for (col in x[[row]]) {
      if (is.null(data)) {
        p[row, col] <- "x"
      }
      else {
        p[row, col] <- data[cnt]
      }
      cnt = cnt + 1    
    }
  }
  return(p)
}

make.sparse <- function (A) 
{
  S <- list()
  for (i in 1:nrow(A)) {
    indices <- c()
    for (j in 1:ncol(A)) {
      if (A[i, j]) {
        indices <- c(indices, j)
      }
    }
    S <- c(S, list(indices))
  }
  return(S)
}


my_make_sparse <- function (A,B,offset_B) 
{
  if(nrow(A)!=nrow(B)) stop("in my_make_sparse A and B have unequal rows")
  S <- list()
  for (i in 1:nrow(A)) {
    indices <- c()
    for (j in 1:ncol(A)) {
      if (A[i, j]) {
        indices <- c(indices, j)
      }
    }
    for (j in 1:ncol(B)) {
      if (B[i, j]) {
        indices <- c(indices, j+offset_B)
      }
    }
    S <- c(S, list(indices))
  }
  return(S)
}



##these functions are used to set constraints for total density
get_total_density <- function(x0_in, wdcMerged, points) {
  tw_list <- unique(wdcMerged$tw)
  tot_density <- 0
  for(tw_in in tw_list) {
    tot_density <- tot_density + sum(get_density_vec(x0_in, tw_in))
  }  
  return(tot_density/length(tw_list))
}
get_grad_total_density  <- function(x0_in, wdcMerged, points) {
  tw_list <- unique(wdcMerged$tw)
  grad <- c()
  for(tw_in in tw_list) {
    grad <- rbind(grad, c(0,colSums(get_grad_density_vec(x0_in, tw_in))))
  }
  return(colMeans(grad))
}
##

get_density_vec <- function(x0_in, tw_in) {
  x0 <- c(x0_in[1],0,x0_in[-1])
  return(get_points_density(points, x0, tw_in))
}


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
  density <- density + points_in$places_count*theta1[density_google_places_count] + 
    (points_in$type==3)*theta1[density_bus_col]
  
  return(density)
}

get_grad_density_vec <- function(x0_in, tw_in) {
  x0 <- c(x0_in[1],0,x0_in[-1])
  return(cbind(get_points_density_grad_ridership_col(points, tw_in),
               get_points_density_grad_metro_col(points, tw_in),
               get_points_density_grad_intercept_col(points, tw_in),
               get_points_density_grad_metro_evening_col(points, tw_in),
               get_points_density_grad_places_count_col(points, tw_in),
               get_points_density_grad_bus_col(points, tw_in))
         )
}

get_points_density_grad_ridership_col <- function(points_in, tw_in) {
  density <- ((points_in$type==1)) *points_in$weight
  return(density)
}

get_points_density_grad_metro_col <- function(points_in, tw_in) {
  if(tw_in==0) {
    density <- rep(0, nrow(points_in))
  } else {
    density <- ((points_in$type==2)) *points_in$weight      
  }
  return(density)
}

get_points_density_grad_intercept_col <- function(points_in, tw_in) {
  density <- ((points_in$type==1))
  return(density)
}

get_points_density_grad_metro_evening_col <- function(points_in, tw_in) {
  if(tw_in==0) {
    density <- ((points_in$type==2)) *points_in$weight      
  } else {
    density <- rep(0, nrow(points_in))
  }
  return(density)
}

get_points_density_grad_places_count_col <- function(points_in, tw_in) {
  density <- (points_in$places_count)
  return(density)
}

# get_points_density_grad_cafe_col <- function(points_in, tw_in) {
#   density <- (points_in$cafe)
#   return(density)
# }
# 
# get_points_density_grad_grocery_col <- function(points_in, tw_in) {
#   density <- ((points_in$grocery_or_supermarket))
#   return(density)
# }
# 
# get_points_density_grad_govoffice_col <- function(points_in, tw_in) {
#   density <- ((points_in$local_government_office))
#   return(density)
# }

get_points_density_grad_bus_col <- function(points_in, tw_in) {
  density <- ((points_in$type==3)) *points_in$weight    
  return(density)
}
