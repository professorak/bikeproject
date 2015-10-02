eval_jac_g_structure_val <- eval_jac_g_structure(params, wdcMerged, points, length(theta), length(deltain_stin), length_eta)
eval_hess_struct <- eval_hessian_structure(params,wdcMerged, points, length(theta), length(deltain_stin), length_eta)

v_eval_g <- eval_g(params, wdcMerged, points, length_theta, length(deltain_stin), length_eta)
summary(v_eval_g)

v_eval_jac_g <- eval_jac_g(params, wdcMerged, points, length_theta, length(deltain_stin), length_eta)
v_eval_jac_g[1:10]
summary(v_eval_jac_g)
length(v_eval_jac_g)

v_eval_obj_GMM_model5_obj <- eval_obj_GMM_model5_obj(params, wdcMerged, points, length_theta, length(deltain_stin), length_eta)
v_eval_obj_GMM_model5_obj

v_eval_obj_GMM_model5_grad <- eval_obj_GMM_model5_grad(params, wdcMerged, points, length_theta, length(deltain_stin), length_eta)
summary(v_eval_obj_GMM_model5_grad)

set.seed(1)
obj_factor <- rnorm(1,0,1)  
hessian_lambda <- rnorm(length(v_eval_g),0,1)  
v_eval_hessian <- eval_hessian(params, wdcMerged, points, length_theta, length(deltain_stin), length_eta,
                               obj_factor,hessian_lambda)
v_eval_hessian[1:10]
summary(v_eval_hessian)
length(v_eval_hessian)

