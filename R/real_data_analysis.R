# function that runs the analysis
# it takes in dataframe that has no missingness
# X_dat should be a "matrix" "array" with 1-hot encoding and an explicit intercept column
# If it is not already in this form (like factor encoding or a data.frame class),
# you can run X_dat = model.matrix(~., X_dat)

# k_dat, should be factor vector
# trt_dat should be numeric with -1, 1 coding

# numeric_var_names is a character vector of the names of the numeric columns
# we will use this for scaling them
# easiest way to obtain them is to call: colnames(X_dat)
# and then copy the string, delete the non-numeric vars

# wt_method can be either "energy weighting" or "propensity"
# aug_method can be either "random forest" or "glmnet" (LASSO)

# remember to check larger_outcome_better is the correct sign:
# If a smaller outcome is preferable, (e.g. a time to hospital discharge analysis)
# we would want FALSE
# If a larger outcome is preferable, (e.g. death time) we would want TRUE 
run_analysis <- function(X_dat, y_dat, k_dat, trt_dat, numeric_var_names,
                         aug_crossfit_num = 10, prop_crossfit_num = 10, wt_method,
                         aug_method, larger_outcome_better = FALSE,
                         solver = "weightit", subset_bal_var = NULL,
                         train_frac = 0.25){
  numcol = ncol(X_dat) # numcol is number of columns in 1 hot encoding, including intercept
  
  study_list = levels(k_dat)
  numstudies = length(study_list)
  
  data_list <- vector(mode="list", length=numstudies)
  trt_list <- vector(mode="list", length=numstudies)
  outcome_list <- vector(mode="list", length=numstudies)
  k_list <- vector(mode="list", length=numstudies)
  for (k in 1:numstudies){
    data_list[[k]] = X_dat[(k_dat == study_list[k]) & !is.na(k_dat),]
    trt_list[[k]] = trt_dat[(k_dat == study_list[k]) & !is.na(k_dat)]
    outcome_list[[k]] = y_dat[(k_dat == study_list[k]) & !is.na(k_dat)]
    k_list[[k]] = k_dat[(k_dat == study_list[k]) & !is.na(k_dat)]
  }
  
  train_data_list <- vector(mode="list", length=numstudies)
  train_trt_list <- vector(mode="list", length=numstudies)
  train_outcome_list <- vector(mode="list", length=numstudies)
  train_k_list <- vector(mode="list", length=numstudies)
  
  test_data_list <- vector(mode="list", length=numstudies)
  test_trt_list <- vector(mode="list", length=numstudies)
  test_outcome_list <- vector(mode="list", length=numstudies)
  test_k_list <- vector(mode="list", length=numstudies)
  
  # 3. train-test split
  print("splitting train-test sets")
  
  train_indices = vector(mode="list", length=numstudies)
  full_data_size = rep(NA, numstudies)
  train_size = rep(NA, numstudies)
  test_size = rep(NA, numstudies)
  for (k in 1:numstudies){
    full_data_size[k] = nrow(data_list[[k]])
    train_size[k] = ceiling(full_data_size[k]*train_frac)
    test_size[k] = full_data_size[k] - train_size[k]
    
    train = sample(c(rep(TRUE, train_size[k]), rep(FALSE, test_size[k])), size = full_data_size[k],
                   replace = FALSE)
    train_indices[[k]] = train
    
    train_data_list[[k]] = data_list[[k]][train,]
    train_trt_list[[k]] = trt_list[[k]][train]
    train_outcome_list[[k]] = outcome_list[[k]][train]
    train_k_list[[k]] = k_list[[k]][train]
    
    test_data_list[[k]] = data_list[[k]][!train,]
    test_trt_list[[k]] = trt_list[[k]][!train]
    test_outcome_list[[k]] = outcome_list[[k]][!train]
    test_k_list[[k]] = k_list[[k]][!train]
  }
  
  
  # NO SCALING OF THE BINARY OUTCOME HERE
  print("scaling X data")
  
  train_Xs = do.call(rbind,train_data_list)
  
  train_sx_nonweighted <- vector(mode="list", length=length(numeric_var_names))
  names(train_sx_nonweighted) <- numeric_var_names
  
  train_scaled_data_list = train_data_list
  test_scaled_data_list = test_data_list
  for (col in names(train_sx_nonweighted)){
    train_sx_nonweighted[[col]] <- sd(train_Xs[,col])
    for (k in 1:numstudies){
      train_scaled_data_list[[k]][,col] = train_scaled_data_list[[k]][,col] / (2*train_sx_nonweighted[[col]])
      
      test_scaled_data_list[[k]][,col] = test_scaled_data_list[[k]][,col] / (2*train_sx_nonweighted[[col]])
    }
  }
  
  # 3b. compute weights
  train_weight_list <- vector(mode="list", length=numstudies)
  test_weight_list <- vector(mode="list", length=numstudies)
  list_of_train_crossfit_propensities <- vector(mode="list", length=numstudies)
  list_of_test_crossfit_propensities <- vector(mode="list", length=numstudies)
  list_of_bal_tab <- vector(mode="list", length=numstudies)
  train_bal_dat_list <- vector(mode="list", length=numstudies)
  if (wt_method == "energy weighting"){
    print("fitting energy weights")
    for (k in 1:numstudies){
      if (is.null(subset_bal_var)){
        train_test_scaled_dat = data.frame(rbind(train_scaled_data_list[[k]], test_scaled_data_list[[k]]))
        
        train_bal_dat_list[[k]] = train_scaled_data_list[[k]]
        
        col_to_keep = rep(TRUE,numcol)
      } else {
        col_to_keep = colnames(train_scaled_data_list[[k]]) %in% subset_bal_var
        
        train_bal_dat_list[[k]] = train_scaled_data_list[[k]][,col_to_keep]
        
        train_test_scaled_dat = data.frame(rbind(
          train_bal_dat_list[[k]],
          test_scaled_data_list[[k]][,col_to_keep]))
      }
      
      train_test_trt = c(train_trt_list[[k]], test_trt_list[[k]])
      
      if (solver == "cccp"){
        # watch out, really slow! seems to take an hour for one of the studies
        # ebqs_cccp <- energy_balance(train_test_scaled_dat, train_test_trt, solver = "cccp")
        # train_weight_list[[k]] = ebqs_cccp$weights[1:train_size[k]]
        # test_weight_list[[k]] = ebqs_cccp$weights[(train_size[k]+1):(train_size[k]+test_size[k])]
      } else if (solver == "osqp"){
        ebqs_osqp <- energy_balance(train_test_scaled_dat, train_test_trt, solver = "osqp")
        train_weight_list[[k]] = ebqs_osqp$weights[1:train_size[k]]
        test_weight_list[[k]] = ebqs_osqp$weights[(train_size[k]+1):(train_size[k]+test_size[k])]
      } else if (solver == "weightit") {
        x_and_trt = data.frame(train_test_scaled_dat, trt = train_test_trt)
        n <- names(train_test_scaled_dat)
        f <- as.formula(paste("trt ~", paste(n[!n %in% "trt"], collapse = " + ")))
        weightit_weights = WeightIt::weightit(formula = f, data = x_and_trt,
                                              estimand = "ATE", method = "energy")
        train_weight_list[[k]] = weightit_weights$weights[1:train_size[k]]
        test_weight_list[[k]] = weightit_weights$weights[(train_size[k]+1):(train_size[k]+test_size[k])]
      }
      
    }
  } else if (wt_method == "propensity"){
    print("fitting propensities")
    for (k in 1:numstudies){
      kth_train_and_test_propens = propens_glmnet_crossfit_with_testset(x = train_scaled_data_list[[k]][,2:numcol],
                                                                        trt = train_trt_list[[k]],
                                                                        K = prop_crossfit_num,
                                                                        trtcoding = "minusplus",
                                                                        test_x = test_scaled_data_list[[k]][,2:numcol])
      list_of_train_crossfit_propensities[[k]] = kth_train_and_test_propens$propensvec
      list_of_test_crossfit_propensities[[k]] = kth_train_and_test_propens$testpropensvec
      
      train_weight_list[[k]] =  1/(train_trt_list[[k]]*list_of_train_crossfit_propensities[[k]] + (1-train_trt_list[[k]])/2)
      
      test_weight_list[[k]] =  1/(test_trt_list[[k]]*list_of_test_crossfit_propensities[[k]] + (1-test_trt_list[[k]])/2)
    }
  }
  
  # normalize the weights so that they sum up to the sample size (train and test separately)
  for (k in 1:numstudies){
    train_weight_list[[k]] = length(train_weight_list[[k]])*train_weight_list[[k]]/sum(train_weight_list[[k]])
    test_weight_list[[k]] = length(test_weight_list[[k]])*test_weight_list[[k]]/sum(test_weight_list[[k]])
    
  }
  
  
  # print the balance
  print("checking balance")
  for (k in 1:numstudies){
    train_scaled_data_list[[k]][,col_to_keep]
    
    # train balance    
    ind_optweight_bal = cobalt::bal.tab(train_bal_dat_list[[k]], treat = train_trt_list[[k]],
                                        weights = train_weight_list[[k]],
                                        disp=c("means"), un = TRUE, 
                                        stats = c("mean.diffs"))
    list_of_bal_tab[[k]] = ind_optweight_bal
    
    # double-checking the max std mean diff
    print(paste(k, "th study train size:", train_size[k]))
    print(paste(k, "th training set balance:",max(abs(list_of_bal_tab[[k]]$Balance$Diff.Adj))))
  }
  

  # 3d. crossfit augmentations using glmnet
  list_of_train_crossfit_augmentations <- vector(mode="list", length=numstudies)
  list_of_test_crossfit_augmentations_1 <- vector(mode="list", length=numstudies)
  list_of_test_crossfit_augmentations_0 <- vector(mode="list", length=numstudies)
  list_of_train_crossfit_augmentations_1 <- vector(mode="list", length=numstudies)
  list_of_train_crossfit_augmentations_0 <- vector(mode="list", length=numstudies)
  if (aug_method == "random forest"){
    print("fitting random forest augmentations")
    for (k in 1:numstudies){
      kth_train_and_test_aug = rf_aug_kfold_crossfit_with_testset(x = train_scaled_data_list[[k]][,2:numcol],
                                                                  y = train_outcome_list[[k]],
                                                                  trt = train_trt_list[[k]],
                                                                  prop_score_x = list_of_train_crossfit_propensities[[k]],
                                                                  use.crossfitting = TRUE,
                                                                  K = aug_crossfit_num,
                                                                  interactions = TRUE,
                                                                  test_x = test_scaled_data_list[[k]][,2:numcol],
                                                                  trtcoding = "minusplus")
      
      list_of_train_crossfit_augmentations[[k]] = kth_train_and_test_aug$predvec
      list_of_test_crossfit_augmentations_1[[k]] = kth_train_and_test_aug$testpredvec_1
      list_of_test_crossfit_augmentations_0[[k]] = kth_train_and_test_aug$testpredvec_0
      list_of_train_crossfit_augmentations_1[[k]] = kth_train_and_test_aug$train_predvec_1
      list_of_train_crossfit_augmentations_0[[k]] = kth_train_and_test_aug$train_predvec_0
    }
  } else if (aug_method == "glmnet"){
    print("fitting glmnet augmentations")
    for (k in 1:numstudies){
      kth_train_and_test_aug = glmnet_aug_kfold_crossfit_with_testset(x = train_scaled_data_list[[k]][,2:numcol],
                                                                      y = train_outcome_list[[k]],
                                                                      trt = train_trt_list[[k]],
                                                                      prop_score_x = list_of_train_crossfit_propensities[[k]],
                                                                      use.crossfitting = TRUE,
                                                                      K = aug_crossfit_num,
                                                                      predtype = "link",
                                                                      family = "gaussian",
                                                                      interactions = TRUE,
                                                                      cv.glmnet.args = NULL,
                                                                      test_x = test_scaled_data_list[[k]][,2:numcol])
      
      list_of_train_crossfit_augmentations[[k]] = kth_train_and_test_aug$predvec
      list_of_test_crossfit_augmentations_1[[k]] = kth_train_and_test_aug$testpredvec_1
      list_of_test_crossfit_augmentations_0[[k]] = kth_train_and_test_aug$testpredvec_0
      list_of_train_crossfit_augmentations_1[[k]] = kth_train_and_test_aug$train_predvec_1
      list_of_train_crossfit_augmentations_0[[k]] = kth_train_and_test_aug$train_predvec_0
    }
  }

  
  # 3e. Misc. precomputing / data prep
  train_Yaug_list_weighted <- vector(mode="list", length=numstudies)
  for (k in 1:numstudies){
    train_Yaug_list_weighted[[k]] =  train_outcome_list[[k]] - list_of_train_crossfit_augmentations[[k]]
  }
  
  # for weighted learning
  train_x_times_trt <- vector(mode="list", length=numstudies)
  for (k in 1:numstudies){
    train_x_times_trt[[k]] =  as.vector(train_trt_list[[k]]) * train_scaled_data_list[[k]]
  }
  
  numcovariates = dim(data_list[[1]])[2] - 1
  
  # 3f. Fit models to get solution paths
  
  # fit weighted learning -> regress y_aug on x_times_trt, weighted by wi
  print("fitting metaanalysis model")
  train_weighted_fit_glmnet <- fit_metacoop_cpp(train_x_times_trt,
                                                train_Yaug_list_weighted,
                                                numcovariates + 1, 100,
                                                train_size, train_weight_list,
                                                c(0,1))
  
  # 3g. Use various CIC to get particular coefficients along paths
  # Use CIC to pick lambda for weighted_fit and Alearning_fit
  # note that this version picks them jointly
  # (a single lambda over all studies rather than K lambdas)
  train_metacoop_w_CIC_glmnet <- all_CICs_pick_lambda_beta(metacoop_fit = train_weighted_fit_glmnet,
                                                           data_list = train_scaled_data_list,
                                                           numobsperstudy = train_size,
                                                           numstudies = numstudies,
                                                           outcome_list = train_outcome_list,
                                                           list_of_crossfit_augmentations = list_of_train_crossfit_augmentations,
                                                           list_of_weights = train_weight_list,
                                                           treatment_list = train_trt_list)
  
  # Modified CIC picker
  # SELF-NORMALIZED CIC
  snCIC_W_curves_glmnet = self_normalized_CIC(train_metacoop_w_CIC_glmnet$concordance_per_lambda_study,
                                              train_metacoop_w_CIC_glmnet$numobsperstudy,
                                              train_metacoop_w_CIC_glmnet$betasize_per_lambda_study)
  snW_lambdadiff_glmnet = lambdadiff(snCIC_W_curves_glmnet)
  snW_coords_glmnet = c(snW_lambdadiff_glmnet$maxvariationrow_kappa,
                        snW_lambdadiff_glmnet$maxCIC_lambda)
  metacoop_snW_beta_glmnet = train_weighted_fit_glmnet$return_beta[,snW_coords_glmnet[2]]
  
  # USE VIC TO SELECT
  VIC_results <- vector(mode="list", length=numstudies)
  metacoop_beta_VIC_list <- vector(mode="list", length=numstudies)
  for (k in 1:numstudies){
    VIC_results[[k]] = VIC(y=train_outcome_list[[k]], x=train_scaled_data_list[[k]],
                           weights=train_weight_list[[k]], trt=train_trt_list[[k]],
                           larger_outcome_better=larger_outcome_better, cutoff=0,
                           aug_1=list_of_train_crossfit_augmentations_1[[k]],
                           aug_0=list_of_train_crossfit_augmentations_0[[k]],
                           coef_mat=train_weighted_fit_glmnet$return_beta[((k-1)*numcol+1):(k*numcol),],
                           gamma_scalingfactor = NULL)
    
    metacoop_beta_VIC_list[[k]] = VIC_results[[k]]$coef
  }
  
  metacoop_beta_VIC = unlist(metacoop_beta_VIC_list)
  
  # FIT POOLED MODEL
  print("fitting pooled model")
  pooled_x_times_trt = do.call(rbind,train_x_times_trt)
  pooled_Yaug_list_weighted_glmnet = unlist(train_Yaug_list_weighted)
  pooled_weights = unlist(train_weight_list)
  
  stacked_weighted_fit <- glmnet::cv.glmnet(pooled_x_times_trt, pooled_Yaug_list_weighted_glmnet,
                                    weights = pooled_weights, intercept = FALSE,
                                    standardize = FALSE,
                                    penalty.factor = c(0,rep(1, (numcol-1))))
  stacked_W_coef = coef(stacked_weighted_fit, s="lambda.min")
  
  # FIT SEPARATE LEARNING MODEL (weighted learning on each dataset)
  print("fitting separate models")
  separate_weighted_fits <- vector(mode="list", length=numstudies)
  for (k in 1:numstudies){
    # note that intercept is already included
    entire_list <- glmnet::cv.glmnet(train_x_times_trt[[k]],
                             train_Yaug_list_weighted[[k]],
                             weights = train_weight_list[[k]], intercept = FALSE,
                             standardize = FALSE,
                             penalty.factor = c(0,rep(1, (numcol-1))))
    
    separate_weighted_fits[[k]] = coef(entire_list,
                                       s="lambda.min")[2:(numcovariates+2)]
  }
  separate_W_coef = unlist(separate_weighted_fits)
  
  # FIT POOLED MODEL WITH INDICATORS
  print("fitting pooled model with indicators")
  train_x_times_trt_with_indicators <- vector(mode="list", length=numstudies)
  for (k in 1:numstudies){
    k_list_times_trt = as.vector(train_trt_list[[k]]) * model.matrix(~.-1, data.frame(train_k_list[[k]]))
    train_x_times_trt_with_indicators[[k]] = cbind(train_x_times_trt[[k]], k_list_times_trt)
  }
  
  pooled_x_times_trt_with_indicators = do.call(rbind,train_x_times_trt_with_indicators)
  stacked_weighted_with_indicators_fit <- glmnet::cv.glmnet(pooled_x_times_trt_with_indicators,
                                                    pooled_Yaug_list_weighted_glmnet,
                                                    weights = pooled_weights, intercept = FALSE,
                                                    standardize = FALSE,
                                                    penalty.factor = c(0, rep(1, ncol(pooled_x_times_trt_with_indicators) - 1)))
  stacked_W_indicators_coef = coef(stacked_weighted_with_indicators_fit, s="lambda.min")
  
  # 4. Estimate value function on test data
  print("computing value function on k test sets")
  stacked_test_valfunc_ipw = rep(NA,numstudies)
  separate_test_valfunc_ipw = rep(NA,numstudies)
  CIC_test_valfunc_ipw = rep(NA,numstudies)
  VIC_test_valfunc_ipw = rep(NA,numstudies)
  stacked_indicators_test_valfunc_ipw = rep(NA,numstudies)
  
  stacked_test_valfunc_aipw = rep(NA,numstudies)
  separate_test_valfunc_aipw = rep(NA,numstudies)
  CIC_test_valfunc_aipw = rep(NA,numstudies)
  VIC_test_valfunc_aipw = rep(NA,numstudies)
  stacked_indicators_test_valfunc_aipw = rep(NA,numstudies)
  # getting the estimate_test_val_func
  for (k in 1:numstudies){
    startidx = ((k-1)*(numcovariates+1)+1)
    endidx = startidx + numcovariates
    
    # for the indicator learning,
    # bind the hospital indicator to the test set
    indicator_testset = cbind(test_scaled_data_list[[k]], model.matrix(~.-1, data.frame(test_k_list[[k]])))
    
    stacked_test_valfunc_ipw[k] = estimate_test_val_func_ipw(y_test = test_outcome_list[[k]],
                                                             x_test = test_scaled_data_list[[k]],
                                                             weight_test = test_weight_list[[k]],
                                                             trt_test = test_trt_list[[k]],
                                                             coefficients = stacked_W_coef[2:(numcovariates + 2)],
                                                             larger.outcome.better = larger_outcome_better,
                                                             cutoff = 0)
    separate_test_valfunc_ipw[k] = estimate_test_val_func_ipw(y_test = test_outcome_list[[k]],
                                                              x_test = test_scaled_data_list[[k]],
                                                              weight_test = test_weight_list[[k]],
                                                              trt_test = test_trt_list[[k]],
                                                              coefficients = separate_W_coef[startidx:endidx],
                                                              larger.outcome.better = larger_outcome_better,
                                                              cutoff = 0)
    CIC_test_valfunc_ipw[k] = estimate_test_val_func_ipw(y_test = test_outcome_list[[k]],
                                                         x_test = test_scaled_data_list[[k]],
                                                         weight_test = test_weight_list[[k]],
                                                         trt_test = test_trt_list[[k]],
                                                         coefficients = metacoop_snW_beta_glmnet[startidx:endidx],
                                                         larger.outcome.better = larger_outcome_better,
                                                         cutoff = 0)
    VIC_test_valfunc_ipw[k] = estimate_test_val_func_ipw(y_test = test_outcome_list[[k]],
                                                         x_test = test_scaled_data_list[[k]],
                                                         weight_test = test_weight_list[[k]],
                                                         trt_test = test_trt_list[[k]],
                                                         coefficients = metacoop_beta_VIC[startidx:endidx],
                                                         larger.outcome.better = larger_outcome_better,
                                                         cutoff = 0)
    stacked_indicators_test_valfunc_ipw[k] = estimate_test_val_func_ipw(y_test = test_outcome_list[[k]],
                                                                        x_test = indicator_testset,
                                                                        weight_test = test_weight_list[[k]],
                                                                        trt_test = test_trt_list[[k]],
                                                                        coefficients = stacked_W_indicators_coef[2:length(stacked_W_indicators_coef)],
                                                                        larger.outcome.better = larger_outcome_better,
                                                                        cutoff = 0)
    
    stacked_test_valfunc_aipw[k] = value_func_AIPW(y=test_outcome_list[[k]],
                                                   x=test_scaled_data_list[[k]],
                                                   weight = test_weight_list[[k]],
                                                   trt=test_trt_list[[k]],
                                                   coefficients=stacked_W_coef[2:(numcovariates + 2)],
                                                   larger.outcome.better = larger_outcome_better,
                                                   cutoff = 0,
                                                   aug_1=list_of_test_crossfit_augmentations_1[[k]],
                                                   aug_0=list_of_test_crossfit_augmentations_0[[k]])
    
    separate_test_valfunc_aipw[k] = value_func_AIPW(y=test_outcome_list[[k]],
                                                    x=test_scaled_data_list[[k]],
                                                    weight = test_weight_list[[k]],
                                                    trt=test_trt_list[[k]],
                                                    coefficients=separate_W_coef[startidx:endidx],
                                                    larger.outcome.better = larger_outcome_better,
                                                    cutoff = 0,
                                                    aug_1=list_of_test_crossfit_augmentations_1[[k]],
                                                    aug_0=list_of_test_crossfit_augmentations_0[[k]])
    
    CIC_test_valfunc_aipw[k] = value_func_AIPW(y=test_outcome_list[[k]],
                                               x=test_scaled_data_list[[k]],
                                               weight = test_weight_list[[k]],
                                               trt=test_trt_list[[k]],
                                               coefficients=metacoop_snW_beta_glmnet[startidx:endidx],
                                               larger.outcome.better = larger_outcome_better,
                                               cutoff = 0,
                                               aug_1=list_of_test_crossfit_augmentations_1[[k]],
                                               aug_0=list_of_test_crossfit_augmentations_0[[k]])
    
    VIC_test_valfunc_aipw[k] = value_func_AIPW(y=test_outcome_list[[k]],
                                               x=test_scaled_data_list[[k]],
                                               weight = test_weight_list[[k]],
                                               trt=test_trt_list[[k]],
                                               coefficients=metacoop_beta_VIC[startidx:endidx],
                                               larger.outcome.better = larger_outcome_better,
                                               cutoff = 0,
                                               aug_1=list_of_test_crossfit_augmentations_1[[k]],
                                               aug_0=list_of_test_crossfit_augmentations_0[[k]])
    
    stacked_indicators_test_valfunc_aipw[k] = value_func_AIPW(y=test_outcome_list[[k]],
                                                              x=indicator_testset,
                                                              weight= test_weight_list[[k]],
                                                              trt=test_trt_list[[k]],
                                                              coefficients=stacked_W_indicators_coef[2:length(stacked_W_indicators_coef)],
                                                              larger.outcome.better = larger_outcome_better,
                                                              cutoff = 0,
                                                              aug_1=list_of_test_crossfit_augmentations_1[[k]],
                                                              aug_0=list_of_test_crossfit_augmentations_0[[k]])
    
  }
  
  return(list(pooled_aipw = stacked_test_valfunc_aipw,
              separate_aipw = separate_test_valfunc_aipw,
              pooled_indicator_aipw = stacked_indicators_test_valfunc_aipw,
              CIC_meta_aipw = CIC_test_valfunc_aipw,
              VIC_meta_aipw = VIC_test_valfunc_aipw,
              pooled_coef = stacked_W_coef,
              separate_coef = separate_W_coef,
              pooled_indicator_coef = stacked_W_indicators_coef,
              CIC_meta_coef = metacoop_snW_beta_glmnet,
              VIC_meta_coef = metacoop_beta_VIC,
              k_list = k_list))
}
