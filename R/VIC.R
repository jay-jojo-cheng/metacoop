# study-specific VIC
# step 1: compute VIC(lambda,gamma) for gamma sequence
# step 2: normalize curves along rows (along a fixed lambda)
# step 3: pick curve with greatest TV norm
# step 4: pick argmax of this curve
# returns:
# coefficients
# optimal gamma and lambda (and their idx)
# the particular curve for that gamma
# the normalized VIC curves
VIC <- function(y, x, weights, trt,
                larger_outcome_better, cutoff,
                aug_1, aug_0,
                coef_mat, gamma_scalingfactor = NULL){
  ############
  # # temp for testing
  # k=7
  # y = train_outcome_list[[k]]
  # x = train_scaled_data_list[[k]]
  # weights = train_weight_list[[k]]
  # trt = train_trt_list[[k]]
  # larger_outcome_better = FALSE
  # cutoff = 0
  # aug_1 = list_of_train_crossfit_augmentations_1[[k]]
  # aug_0 = list_of_train_crossfit_augmentations_0[[k]]
  # coef_mat = train_weighted_fit_glmnet$return_beta[((k-1)*42+1):(k*42),]
  # gamma_scalingfactor = NULL
  # # temp for testing
  ############
  n = length(y)
  
  # step 1: compute VIC(lambda,gamma) for gamma sequence
  orig_num_lambda = dim(coef_mat)[2]
  
  if(is.null(gamma_scalingfactor)){
    gamma_scalingfactor = 10^seq(from = 5, to = -5, length.out = 1000)
  }
  num_gamma = length(gamma_scalingfactor)
  
  beta_zero_norm = colSums(abs(coef_mat) > 0)
  first_idx_considered = which(beta_zero_norm > 0)[1] # hacky way to get rid
  # of zero coefficients. Not supposed to have them according to CIC/VIC theory
  # and the way lambda sequence is chosen.
  
  trunc_beta_zero_norm = beta_zero_norm[first_idx_considered:orig_num_lambda]
  trunc_coef_mat = coef_mat[,first_idx_considered:orig_num_lambda]
  trunc_num_lambda = dim(trunc_coef_mat)[2]
  
  V_hat_per_lambda = numeric(trunc_num_lambda)
  for (ell in 1:trunc_num_lambda) {
    V_hat_per_lambda[ell] = value_func_AIPW(y = y, x = x, weights = weights, trt = trt,
                                            coefficients = trunc_coef_mat[,ell],
                                            larger.outcome.better = larger_outcome_better,
                                            cutoff = cutoff,
                                            aug_1 = aug_1, aug_0 = aug_0)
  }
  
  # switch signs if smaller outcome is better
  if (larger_outcome_better == FALSE){
    V_hat_per_lambda = -V_hat_per_lambda
  }
  
  # V_hat_per_lambda is length num_lambda
  # gamma_scalingfactor is length num_gamma
  # beta_zero_norm is length num_lambda
  # note the second term, n*V is broadcasted into an array
  VIC_mat_unnormed = (- log(n) * gamma_scalingfactor %o% trunc_beta_zero_norm) + (n*matrix(replicate(num_gamma,V_hat_per_lambda),nrow=num_gamma, byrow = TRUE))
  
  # step 2: normalize curves along rows (along a fixed lambda)
  VIC_curves_normed = self_norm_rows(VIC_mat_unnormed)
  
  # step 3&4: pick curve with greatest TV norm, pick argmax of this curve
  VIC_res = lambdadiff(VIC_curves_normed)
  VIC_res$curves = VIC_curves_normed
  VIC_res$trunc_start = first_idx_considered
  VIC_res$max_untrunc_lambda = VIC_res$maxCIC_lambda + first_idx_considered - 1
  VIC_res$coef = coef_mat[,VIC_res$max_untrunc_lambda]
  
  return(VIC_res)
}


# which(VIC_curves_normed == max(VIC_curves_normed),arr.ind = TRUE)

# PLOTTING
# view plot
# image(1:ncol(VIC_mat_unnormed), 1:nrow(VIC_mat_unnormed), t(VIC_mat_unnormed), col = cm.colors(60), axes = FALSE)
# image(1:ncol(VIC_curves_normed), 1:nrow(VIC_curves_normed), t(VIC_curves_normed), col = cm.colors(60), axes = FALSE)
# plot(VIC_curves_normed[VIC_res$maxvariationrow_kappa,])


# aipw value function
value_func_AIPW <- function(y, x, propensity = NULL, weights = NULL, trt,
                            coefficients, larger.outcome.better = TRUE,
                            cutoff = 0, aug_1, aug_0){
  scorefunc = x %*% as.vector(coefficients)
  
  # print("larger outcome better:")
  # print(larger.outcome.better)
  if(larger.outcome.better){
    indicator_scorefunc_1 <- ifelse(scorefunc > cutoff, 1, 0)
    indicator_scorefunc_0 <- ifelse(scorefunc <= cutoff, 1 ,0)
  } else{
    indicator_scorefunc_1 <- ifelse(scorefunc < cutoff, 1, 0)
    indicator_scorefunc_0 <- ifelse(scorefunc >= cutoff, 1, 0)
  }
  
  
  if (is.null(weights)){
    aipw_weights = 1/(propensity*indicator_scorefunc_1 + (1-propensity)*indicator_scorefunc_0)
  } else{
    aipw_weights = weights
  }
  
  indicator_trt_1 = ifelse(trt == 1, 1, 0)
  indicator_trt_0 = ifelse(trt != 1, 1, 0)
  
  summands = indicator_scorefunc_1*(indicator_trt_1*(y - aug_1)*aipw_weights + aug_1) + indicator_scorefunc_0*(indicator_trt_0*(y - aug_0)*aipw_weights + aug_0)
  
  return(mean(summands, na.rm  = TRUE))
}
