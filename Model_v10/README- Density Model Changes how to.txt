1. In test_gradients.R where have to specify the column name like density_metro_col <<- 4
2. eval_grad_lambda_theta_new: Here there are gradients coded wrt density coefficients. Have to be careful in specifying density, as it is gradient wrt density coefficients.
3. Changed implementation so that all relevant function are called from util_functions.R


