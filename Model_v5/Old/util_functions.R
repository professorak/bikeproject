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

get_density_vec <- function(x0_in) {
  x0 <- c(x0_in[1],0,x0_in[-1])
  points_density <- ((points$type==1)*x0[density_ridership_col] + 
                       (points$type==2)*x0[density_metro_col]) *points$weight + 
    (points$type==1)*x0[density_intercept_col]
  return(points_density)
}

get_grad_density_vec <- function(x0_in) {
  x0 <- c(x0_in[1],0,x0_in[-1])
  return(cbind((points$type==1)*points$weight, (points$type==2)*points$weight, (points$type==1)))
}

## these functions are used in setting the constraint of total density matching a prespecified amount to achieve identifiability
get_total_density <- function(x0_in, wdcMerged, points) {
  return(sum(get_density_vec(x0_in)))
}
get_grad_total_density  <- function(x0_in, wdcMerged, points) {
  grad <-c(0,colSums(get_grad_density_vec(x0_in))) 
  return(grad)
}
######

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

