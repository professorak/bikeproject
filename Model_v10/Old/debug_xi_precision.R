theta_2 <- solve(t(X) %*% X) %*% (t(X) %*% deltain)
theta_2 <- solve(t(X) %*% Z %*% A_N %*% t(Z) %*% X) %*%
  (t(X) %*% Z %*% A_N %*% t(Z) %*% deltain)

#check precision of solve.
a <- solve(t(X) %*% Z %*% A_N %*% t(Z) %*% X) %*%
  (t(X) %*% Z %*% A_N %*% t(Z) %*% X)
b <- solve(t(X) %*% X) %*% (t(X) %*% X)
c <- solve(t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted) %*% X) %*% 
  (t(X) %*% Z_weighted %*% A_N %*% t(Z_weighted) %*% X)
D_N <-  diag(ncol(X))
d <- solve(t(X) %*% X %*% D_N %*% t(X) %*% X) %*%
  (t(X) %*% X %*% D_N %*% t(X) %*% X)

Z_temp <- Z[,c(1:71)]
E_N <-  diag(ncol(Z_temp))
e <- solve(t(X) %*% Z_temp %*% E_N %*% t(Z_temp) %*% X) %*%
  (t(X) %*% Z_temp %*% E_N %*% t(Z_temp) %*% X)
summary(c(e))

#divide each column of Z by its colSum
Z_norm <- t(Z)
Z_norm <- Z_norm/rowSums(Z_norm)
Z_norm <- t(Z_norm)
F_N <-  diag(ncol(Z_norm))
f <- solve(t(X) %*% Z_norm %*% F_N %*% t(Z_norm) %*% X) %*%
  (t(X) %*% Z_norm %*% F_N %*% t(Z_norm) %*% X)
summary(c(f))
theta_2 <- solve(t(X) %*% Z_norm %*% A_N %*% t(Z_norm) %*% X) %*%
  (t(X) %*% Z_norm %*% A_N %*% t(Z_norm) %*% deltain)






