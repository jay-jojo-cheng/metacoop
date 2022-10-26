# functions for estimating the value function over the test set
# adapted from p2p package
# ipw value function
estimate_test_val_func_ipw <- function(y_test, x_test, propensity_test = NULL,
                                       weight_test = NULL, trt_test,
                                       coefficients, larger.outcome.better = TRUE,
                                       cutoff = 0){
  if (is.null(weight_test)){
    ipw_weights = 1/(propensity_test*(trt_test == 1) + (1-propensity_test)*(trt_test != 1))
  } else{
    ipw_weights = weight_test
  }
  
  scorefunc = x_test %*% as.vector(coefficients)
  
  if(larger.outcome.better){
    agree1 <- (scorefunc > cutoff) & (trt_test == 1)
    agree0 <- (scorefunc <= cutoff) & (trt_test != 1)
  }
  else{
    agree1 <- (scorefunc < cutoff) & (trt_test == 1)
    agree0 <- (scorefunc >= cutoff) & (trt_test != 1)
  }
  
  return(weighted.mean(y_test[agree1 | agree0],
                       w = ipw_weights[agree1 | agree0],
                       na.rm  = TRUE))
}

# NOTE THAT THIS IS OBSOLETE, USE value_func_AIPW in VIC.R INSTEAD
# # aipw value function
# estimate_test_val_func_aipw <- function(y_test, x_test, propensity_test = NULL,
#                                         weight_test = NULL, trt_test,
#                                         coefficients, larger.outcome.better = TRUE,
#                                         cutoff = 0, aug_test_1, aug_test_0){
#   if (is.null(weight_test)){
#     aipw_weights = 1/(propensity_test*indicator_scorefunc_1 + (1-propensity_test)*indicator_scorefunc_0)
#   } else{
#     aipw_weights = weight_test
#   }
#   
#   scorefunc = x_test %*% as.vector(coefficients)
#   
#   if(larger.outcome.better){
#     indicator_scorefunc_1 <- ifelse(scorefunc > cutoff, 1, 0)
#     indicator_scorefunc_0 <- ifelse(scorefunc <= cutoff, 1 ,0)
#   }
#   else{
#     indicator_scorefunc_1 <- ifelse(scorefunc < cutoff, 1, 0)
#     indicator_scorefunc_0 <- ifelse(scorefunc >= cutoff, 1, 0)
#   }
#   
#   indicator_trt_1 = ifelse(trt_test == 1, 1, 0)
#   indicator_trt_0 = ifelse(trt_test == -1, 1, 0)
#   
#   summands = indicator_scorefunc_1*(indicator_trt_1*(y_test - aug_test_1)*aipw_weights + aug_test_1) + indicator_scorefunc_0*(indicator_trt_0*(y_test - aug_test_0)*aipw_weights + aug_test_0)
#   
#   return(mean(summands, na.rm  = TRUE))
# }