#
dem_hat_T <- exp(-3)/(exp(-3)+1)
x0 <- exp(-2)
ftol <- 1e-12
htol <- 1e-12
itr = 1
repeat{
  dem_T <- x0/(x0+1)
  obj <- max(abs(dem_T-dem_hat_T))
  cat(obj," ")
  cat(mean(log(x0))," ")
  itr = itr + 1
  x0prev <- x0

  #x0 <- x0*(dem_hat_T/dem_T)
  x0 <- ((1+1/x0)*dem_T/dem_hat_T-1)^(-1)
  
  if(max(abs(log(x0)-log(x0prev))) < ftol | obj < htol) {
    #      if(norm(log(x0)-log(x0prev),"2")/length(x0) < ftol ){
    break
  }
  #cat(norm(log(x0[stocked_list])-log(x0prev),"2")/length(x0[stocked_list])," ;; ")
  cat(max(abs(log(x0)-log(x0prev)))," ;; ")
}
print(paste("itr: ", itr))




