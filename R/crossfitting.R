# propensity crossfitting function adapted from personalized package
# this version uses LOOCV and creates balanced folds in terms of 
# treatment and control
kfold_crossfit_propensity_glmnet_loocv <- function(x, trt,
                                                   use.crossfitting = TRUE,
                                                   K = 5,
                                                   trtcoding = "minusplus")
{
  if (trtcoding == "minusplus"){
    trt[trt == -1] = 0
  }
  
  n_1 = NROW(trt[trt == 1])
  n_0 = NROW(trt[trt == 0])
  
  rowswith1 = which(trt == 1)
  rowswith0 = which(trt == 0)
  
  total_n <- NROW(x)
  
  tm <- "deviance"
  
  reorderedrowswith1 = sample(rowswith1)
  reorderedrowswith0 = sample(rowswith0)
  
  foldid = rep(NA, total_n)
  foldid[reorderedrowswith1] = sample(rep(seq(K), length = n_1))
  foldid[reorderedrowswith0] = sample(rep(seq(K), length = n_0))
  
  propensvec <- numeric(total_n)
  
  for (i in seq(K))
  {
    ithfold <- foldid == i
    
    glmnetfolds = sum(!ithfold)
    
    cv.glmnet.args <- list(type.measure = tm, nfolds = glmnetfolds, grouped = FALSE)
    
    glmfit_propens <- do.call(glmnet::cv.glmnet, c(list(y = trt[!ithfold], x = x[!ithfold,,drop=FALSE],
                                                family = "binomial"), cv.glmnet.args))
    
    ## get propensity scores for the held out fold
    propensvec[ithfold] <- unname(drop(stats::predict(glmfit_propens, newx = x[ithfold,,drop=FALSE],
                                               s = "lambda.1se", type = "response")))
  }
  
  ## propensity scores will never be outside of 0 or 1 and
  ## shouldn't have missing values, but this code is a safety
  ## check just in case
  propensvec[is.na(propensvec)] <- mean(propensvec[!is.na(propensvec)])
  propensvec[propensvec <= 0] <- 1e-5
  propensvec[propensvec >= 1] <- 1 - 1e-5
  
  propensvec
}

# similar but k folds. needs to be merged with above
# propensity crossfitting function adapted from personalized package
glmnet_propensity_kfold_crossfit <- function(x, trt, use.crossfitting = TRUE,
                                             K = 5, cv.glmnet.args = NULL)
{
  n <- NROW(x)
  tm <- "auc"
  
  if (is.null(cv.glmnet.args))
  {
    cv.glmnet.args <- list(type.measure = tm, nfolds = 5)
  }
  
  cv.glmnet.args[c("x", "y", "family", "parallel")] <- NULL
  cv.glmnet.args$parallel <- FALSE
  
  if (!("type.measure" %in% names(cv.glmnet.args) ))
  {
    cv.glmnet.args$type.measure <- tm
  }
  
  propensvec <- numeric(n)
  
  foldid = sample(rep(seq(K), length = n))
  
  for (i in seq(K))
  {
    which <- foldid == i
    
    glmfit_propens <- do.call(glmnet::cv.glmnet, c(list(y = trt[!which], x = x[!which,,drop=FALSE],
                                                family = "binomial"), cv.glmnet.args))
    
    ## get propensity scores for the held out fold
    propensvec[which] <- unname(drop(
      stats::predict(glmfit_propens, newx = x[which,,drop=FALSE],
                                             s = "lambda.1se", type = "response")))
  }
  
  ## propensity scores will never be outside of 0 or 1 and
  ## shouldn't have missing values, but this code is a safety
  ## check just in case
  propensvec[is.na(propensvec)] <- mean(propensvec[!is.na(propensvec)])
  propensvec[propensvec <= 0] <- 1e-5
  propensvec[propensvec >= 1] <- 1 - 1e-5
  
  propensvec
}

# new rewritten augmentation crossfitting function by JC 2021 April
crossfit_augmentation_rf <- function(x, y, trt, prop_score_x,
                                     use.crossfitting = TRUE, K=5,
                                     interactions = TRUE,
                                     trtcoding = "minusplus"){
  
  if (trtcoding == "minusplus"){
    trt[trt == -1] = 0
  }
  
  dat_all = data.frame(outcome = y,trt = trt,x)
  dat_1 = data.frame(trt = 1,x)
  dat_0 = data.frame(trt = 0,x)
  
  
  n <- NROW(dat_all)
  # check that all rows are either 0 or 1
  # trt_is_one = (trt == 1)
  # trt_is_zero = (trt == 0)
  # print(n == sum(trt_is_one) + sum(trt_is_zero))
  
  predvec <- numeric(n)
  
  if (use.crossfitting)
  {
    foldid = sample(rep(seq(K), length = n))
    
    for (i in seq(K))
    {
      ithfold <- foldid == i
      
      v.grow <- randomForestSRC::rfsrc(outcome ~ . , data = dat_all[-ithfold, ])
      
      
      if (interactions)
      {
        ## get predictions for trt = 1 and then trt = 0
        v.pred1 = randomForestSRC::predict.rfsrc(v.grow, newdata = dat_1[ithfold,])
        v.pred0 = randomForestSRC::predict.rfsrc(v.grow, newdata = dat_0[ithfold,])
        
        predvec[ithfold] <- (1-prop_score_x[ithfold]) * (v.pred1$predicted) + (prop_score_x[ithfold]) * (v.pred0$predicted)
      } else
      {
        # get predictions for all in ith fold
        v.pred = randomForestSRC::predict.rfsrc(v.grow, newdata = dat_all[ithfold,])
        
        predvec[ithfold] <- v.pred$predicted
      }
      
    }
  } else
  {
    
  }
  
  predvec
}

# rf augmentation for the whole data Feb 2
rf_aug <- function(x, y, trt, prop_score_x = NULL, use.crossfitting = TRUE,
                   K = 5, interactions = TRUE, trtcoding = "minusplus"){
  if (trtcoding == "minusplus"){
    trt[trt == -1] = 0
  }
  
  dat_all = data.frame(outcome = y,trt = trt,x)
  dat_1 = data.frame(trt = 1,x)
  dat_0 = data.frame(trt = 0,x)
  
  n <- NROW(dat_all)
  # check that all rows are either 0 or 1
  # trt_is_one = (trt == 1)
  # trt_is_zero = (trt == 0)
  # print(n == sum(trt_is_one) + sum(trt_is_zero))
  
  predvec <- numeric(n)
  train_predvec_1 <- numeric(n)
  train_predvec_0 <- numeric(n)
  
  if (use.crossfitting)
  {
    foldid = sample(rep(seq(K), length = n))
    
    for (i in seq(K))
    {
      ithfold <- foldid == i
      
      v.grow <- randomForestSRC::rfsrc(outcome ~ . , data = dat_all[-ithfold, ])
      
      
      if (interactions)
      {
        ## get predictions for trt = 1 and then trt = 0
        v.pred1 = randomForestSRC::predict.rfsrc(v.grow, newdata = dat_1[ithfold,])
        v.pred0 = randomForestSRC::predict.rfsrc(v.grow, newdata = dat_0[ithfold,])
        
        if (is.null(prop_score_x)){
          predvec[ithfold] <- (0.5) * (v.pred1$predicted) + (0.5) * (v.pred0$predicted)
        } else{
          predvec[ithfold] <- (1-prop_score_x[ithfold]) * (v.pred1$predicted) + (prop_score_x[ithfold]) * (v.pred0$predicted)
        }
        
        train_predvec_1[ithfold] = v.pred1$predicted
        train_predvec_0[ithfold] = v.pred0$predicted
        
      } else
      {
        stop("No.interactions not implemented yet for this function.
         Not meant for production.")
        
        # # get predictions for all in ith fold
        # v.pred = randomForestSRC::predict.rfsrc(v.grow, newdata = dat_all[ithfold,])
        # 
        # predvec[ithfold] <- v.pred$predicted
      }
      
    }
  } else
  {
    stop("No.crossfitting not implemented yet for this function.
         Not meant for production.")
  }
  
  return(list(predvec=predvec,
              train_predvec_1 = train_predvec_1,
              train_predvec_0 = train_predvec_0))
}

# rf augmentation with testset Oct 21
rf_aug_kfold_crossfit_with_testset <- function(x, y, trt,
                                               prop_score_x = NULL,
                                               use.crossfitting = TRUE,
                                               K = 5,
                                               interactions = TRUE,
                                               test_x,
                                               trtcoding = "minusplus"){
  if (trtcoding == "minusplus"){
    trt[trt == -1] = 0
  }
  
  dat_all = data.frame(outcome = y,trt = trt,x)
  dat_1 = data.frame(trt = 1,x)
  dat_0 = data.frame(trt = 0,x)
  
  # test set, create two datasets with counterfactual treatments
  df_test_1 <- data.frame(trt = 1, test_x)
  df_test_0 <- data.frame(trt = -1, test_x)
  
  n <- NROW(dat_all)
  # check that all rows are either 0 or 1
  # trt_is_one = (trt == 1)
  # trt_is_zero = (trt == 0)
  # print(n == sum(trt_is_one) + sum(trt_is_zero))
  
  predvec <- numeric(n)
  train_predvec_1 <- numeric(n)
  train_predvec_0 <- numeric(n)
  testaugmat_1 = matrix(data = NA, nrow = nrow(test_x), ncol = K)
  testaugmat_0 = matrix(data = NA, nrow = nrow(test_x), ncol = K)
  
  if (use.crossfitting)
  {
    foldid = sample(rep(seq(K), length = n))
    
    for (i in seq(K))
    {
      ithfold <- foldid == i
      
      v.grow <- randomForestSRC::rfsrc(outcome ~ . , data = dat_all[-ithfold, ])
      
      
      if (interactions)
      {
        ## get predictions for trt = 1 and then trt = 0
        v.pred1 = randomForestSRC::predict.rfsrc(v.grow, newdata = dat_1[ithfold,])
        v.pred0 = randomForestSRC::predict.rfsrc(v.grow, newdata = dat_0[ithfold,])
        
        if (is.null(prop_score_x)){
          predvec[ithfold] <- (0.5) * (v.pred1$predicted) + (0.5) * (v.pred0$predicted)
        } else{
          predvec[ithfold] <- (1-prop_score_x[ithfold]) * (v.pred1$predicted) + (prop_score_x[ithfold]) * (v.pred0$predicted)
        }
        
        train_predvec_1[ithfold] = v.pred1$predicted
        train_predvec_0[ithfold] = v.pred0$predicted
        
      } else
      {
        stop("No.interactions not implemented yet for this function.
         Not meant for production.")
        
        # # get predictions for all in ith fold
        # v.pred = randomForestSRC::predict.rfsrc(v.grow, newdata = dat_all[ithfold,])
        # 
        # predvec[ithfold] <- v.pred$predicted
      }
      
      
      ## get test set predictions for trt = 1 and -1
      test_pred1 <- randomForestSRC::predict.rfsrc(v.grow, newdata = df_test_1)
      test_pred0 <- randomForestSRC::predict.rfsrc(v.grow, newdata = df_test_0)
      
      testaugmat_1[,i] <- test_pred1$predicted
      testaugmat_0[,i] <- test_pred0$predicted
      
    }
  } else
  {
    stop("No.crossfitting not implemented yet for this function.
         Not meant for production.")
  }
  
  testpredvec_1 = rowMeans(testaugmat_1)
  testpredvec_0 = rowMeans(testaugmat_0)
  
  return(list(predvec=predvec,
              testpredvec_1 = testpredvec_1,
              testpredvec_0 = testpredvec_0,
              train_predvec_1 = train_predvec_1,
              train_predvec_0 = train_predvec_0))
}


# augmentation crossfitting function adapted from personalized package
glmnet_aug_kfold_crossfit <- function(x, y, trt, prop_score_x, wts = NULL,
                                      use.crossfitting = TRUE,
                                      K = 5,
                                      predtype = c("link", "response"),
                                      family = c("gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian"),
                                      interactions = TRUE, cv.glmnet.args = NULL)
{
  
  predtype <- match.arg(predtype)
  family   <- match.arg(family)
  
  
  if (family == "binomial")
  {
    tm <- "auc"
  } else
  {
    tm <- "mse"
  }
  
  if (is.null(cv.glmnet.args))
  {
    
    cv.glmnet.args <- list(type.measure = tm, nfolds = 5)
  }
  
  cv.glmnet.args[c("x", "y", "family", "weights", "parallel")] <- NULL
  cv.glmnet.args$parallel <- FALSE
  
  
  if (!("type.measure" %in% names(cv.glmnet.args) ))
  {
    cv.glmnet.args$type.measure <- tm
  }
  
  if (is.null(wts))
  {
    wts <- rep(1, NROW(x))
  }
  
  
  if (interactions)
  {
    ## full model for nonzeroness
    df_all <- data.frame(x, trt = trt)
    df_1   <- data.frame(x, trt = 1)
    df_0   <- data.frame(x, trt = -1)
    
    mm_all <- stats::model.matrix(~x*trt-1, data = df_all)
    mm_1   <- stats::model.matrix(~x*trt-1, data = df_1)
    mm_0   <- stats::model.matrix(~x*trt-1, data = df_0)
  } else
  {
    mm_all <- x
  }
  
  n <- NROW(mm_all)
  
  predvec <- numeric(n)
  
  if (use.crossfitting)
  {
    foldid = sample(rep(seq(K), length = n))
    
    for (i in seq(K))
    {
      which <- foldid == i
      
      glmfit_zero_main <- do.call(glmnet::cv.glmnet, c(list(y = y[!which], x = mm_all[!which,,drop=FALSE],
                                                    weights = wts[!which], family = family), cv.glmnet.args))
      
      if (interactions)
      {
        ## get predictions for trt = 1 & -1
        pred1_zerr <- unname(drop(stats::predict(glmfit_zero_main, newx = mm_1[which,,drop=FALSE], s = "lambda.min", type = predtype)))
        pred0_zerr <- unname(drop(stats::predict(glmfit_zero_main, newx = mm_0[which,,drop=FALSE], s = "lambda.min", type = predtype)))
        
        predvec[which] <- (1-prop_score_x[which]) * (pred1_zerr) + (prop_score_x[which]) * (pred0_zerr)
      } else
      {
        ## get predictions for trt = 1 & -1
        pred_zerr <- unname(drop(stats::predict(glmfit_zero_main, newx = mm_all[which,,drop=FALSE], s = "lambda.min", type = predtype)))
        
        predvec[which] <- pred_zerr
      }
      
    }
  } else
  {
    glmfit_zero_main <- do.call(glmnet::cv.glmnet, c(list(y = y, x = mm_all,
                                                  weights = wts, family = family), cv.glmnet.args))
    
    if (interactions)
    {
      ## get predictions for trt = 1 & -1
      pred1_zerr <- unname(drop(stats::predict(glmfit_zero_main, newx = mm_1, s = "lambda.min", type = predtype)))
      pred0_zerr <- unname(drop(stats::predict(glmfit_zero_main, newx = mm_0, s = "lambda.min", type = predtype)))
      
      predvec <- (1-prop_score_x) * pred1_zerr + prop_score_x * pred0_zerr
    } else
    {
      ## get predictions for trt = 1 & -1
      pred_zerr <- unname(drop(stats::predict(glmfit_zero_main, newx = mm_all, s = "lambda.min", type = predtype)))
      
      predvec <- pred_zerr
    }
  }
  
  predvec
}

# version of propensity function where we also return the test set propensities
propens_glmnet_crossfit_with_testset <- function(x, trt,
                                                 use.crossfitting = TRUE,
                                                 K = 5,
                                                 trtcoding = "minusplus",
                                                 test_x,
                                                 glmnet_cvfolds = 5)
{
  if (trtcoding == "minusplus"){
    trt[trt == -1] = 0
  }
  
  n_1 = NROW(trt[trt == 1])
  n_0 = NROW(trt[trt == 0])
  
  rowswith1 = which(trt == 1)
  rowswith0 = which(trt == 0)
  
  total_n <- NROW(x)
  
  tm <- "deviance"
  
  reorderedrowswith1 = sample(rowswith1)
  reorderedrowswith0 = sample(rowswith0)
  
  foldid = rep(NA, total_n)
  foldid[reorderedrowswith1] = sample(rep(seq(K), length = n_1))
  foldid[reorderedrowswith0] = sample(rep(seq(K), length = n_0))
  
  propensvec <- numeric(total_n)
  testpropensmat = matrix(data = NA, nrow = nrow(test_x), ncol = K)
  
  for (i in seq(K))
  {
    ithfold <- foldid == i
    
    if (glmnet_cvfolds == "loocv"){
      glmnetfolds = sum(!ithfold)
    } else{
      glmnetfolds = glmnet_cvfolds
    }
    
    cv.glmnet.args <- list(type.measure = tm, nfolds = glmnetfolds,
                           grouped = FALSE,
                           standardize = FALSE)
    
    glmfit_propens <- do.call(glmnet::cv.glmnet, c(list(y = trt[!ithfold], x = x[!ithfold,,drop=FALSE],
                                                        family = "binomial"), cv.glmnet.args))
    
    ## get propensity scores for the held out fold
    propensvec[ithfold] <- unname(drop(stats::predict(glmfit_propens, newx = x[ithfold,,drop=FALSE],
                                                      s = "lambda.min", type = "response")))
    
    ## get propensity scores for the test set
    testpropensmat[,i] = unname(drop(stats::predict(glmfit_propens, newx = test_x,
                                                    s = "lambda.min", type = "response")))
  }
  
  ## propensity scores will never be outside of 0 or 1 and
  ## shouldn't have missing values, but this code is a safety
  ## check just in case
  propensvec[is.na(propensvec)] <- mean(propensvec[!is.na(propensvec)])
  propensvec[propensvec <= 0] <- 1e-5
  propensvec[propensvec >= 1] <- 1 - 1e-5
  
  testpropensvec = rowMeans(testpropensmat)
  
  
  return(list(propensvec=propensvec, testpropensvec = testpropensvec))
}

# version of augmentation crossfitting function where we also return the test set augmentations
glmnet_aug_kfold_crossfit_with_testset <- function(x, y, trt,
                                                   prop_score_x = NULL,
                                                   wts = NULL,
                                                   use.crossfitting = TRUE,
                                                   K = 5,
                                                   predtype = c("link", "response"),
                                                   family = c("gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian"),
                                                   interactions = TRUE, cv.glmnet.args = NULL,
                                                   test_x)
{
  predtype <- match.arg(predtype)
  family   <- match.arg(family)
  
  if (family == "binomial")
  {
    tm <- "auc"
  } else
  {
    tm <- "mse"
  }
  
  if (is.null(cv.glmnet.args))
  {
    cv.glmnet.args <- list(type.measure = tm, nfolds = 5)
  }
  
  cv.glmnet.args[c("x", "y", "family", "weights", "parallel")] <- NULL
  cv.glmnet.args$parallel <- FALSE
  cv.glmnet.args$standardize = FALSE
  cv.glmnet.args$intercept = TRUE
  
  if (!("type.measure" %in% names(cv.glmnet.args) ))
  {
    cv.glmnet.args$type.measure <- tm
  }
  
  if (is.null(wts))
  {
    wts <- rep(1, NROW(x))
  }
  
  if (interactions)
  {
    ## full model for nonzeroness
    # first one makes a treatment column with actual treatments
    # second and third make a treatment column with counterfactual treatments
    df_all <- data.frame(x, trt = trt)
    df_1   <- data.frame(x, trt = 1)
    df_0   <- data.frame(x, trt = -1)
    
    # -1 is for no intercept
    mm_all <- stats::model.matrix(~x*trt-1, data = df_all)
    mm_1   <- stats::model.matrix(~x*trt-1, data = df_1)
    mm_0   <- stats::model.matrix(~x*trt-1, data = df_0)
    
    # test set, create two datasets with counterfactual treatments
    df_test_1 <- data.frame(test_x, trt = 1)
    df_test_0 <- data.frame(test_x, trt = -1)
    
    mm_test_1 <- stats::model.matrix(~test_x*trt-1, data = df_test_1)
    mm_test_0 <- stats::model.matrix(~test_x*trt-1, data = df_test_0)
  } else
  {
    stop("Error: no-interactions not implemented for this function.
         Not meant for production, just auxiliary function for paper.")
  }
  
  n <- NROW(mm_all)
  
  predvec <- numeric(n)
  testaugmat_1 = matrix(data = NA, nrow = nrow(test_x), ncol = K)
  testaugmat_0 = matrix(data = NA, nrow = nrow(test_x), ncol = K)
  
  if (use.crossfitting)
  {
    foldid = sample(rep(seq(K), length = n))
    
    for (i in seq(K))
    {
      which <- foldid == i
      
      glmfit_zero_main <- do.call(glmnet::cv.glmnet, c(list(y = y[!which], x = mm_all[!which,,drop=FALSE],
                                                            weights = wts[!which], family = family), cv.glmnet.args))
      
      if (interactions)
      {
        ## get predictions for trt = 1 & -1
        pred1_zerr <- unname(drop(stats::predict(glmfit_zero_main,
                                                 newx = mm_1[which,,drop=FALSE],
                                                 s = "lambda.min",
                                                 type = predtype)))
        pred0_zerr <- unname(drop(stats::predict(glmfit_zero_main,
                                                 newx = mm_0[which,,drop=FALSE],
                                                 s = "lambda.min",
                                                 type = predtype)))
        
        if (is.null(prop_score_x)){
          predvec[which] <- 0.5 * (pred1_zerr) + 0.5 * (pred0_zerr)
        } else{
          predvec[which] <- (1-prop_score_x[which]) * (pred1_zerr) + (prop_score_x[which]) * (pred0_zerr)
        }
        
        
        ## get test set predictions for trt = 1 and -1
        test_pred1 <- unname(drop(stats::predict(glmfit_zero_main,
                                                 newx = mm_test_1,
                                                 s = "lambda.min",
                                                 type = predtype)))
        test_pred0 <- unname(drop(stats::predict(glmfit_zero_main,
                                                 newx = mm_test_0,
                                                 s = "lambda.min",
                                                 type = predtype)))
        
        testaugmat_1[,i] <- test_pred1
        testaugmat_0[,i] <- test_pred0
        
      } else
      {
        stop("Error: no-interactions not implemented for this function.
         Not meant for production, just auxiliary function for paper.")
      }
      
    }
  }
  else{
    stop("Error: non-crossfitting not implemented for this function.
         Not meant for production, just auxiliary function for paper.")
  }
  
  testpredvec_1 = rowMeans(testaugmat_1)
  testpredvec_0 = rowMeans(testaugmat_0)
  
  return(list(predvec=predvec,
              testpredvec_1 = testpredvec_1,
              testpredvec_0 = testpredvec_0))
}

# version of RF propensity where we also return the test set propensities
propens_rf_crossfit_with_testset <- function(x, trt,
                                                 use.crossfitting = TRUE,
                                                 K = 5,
                                                 trtcoding = "minusplus",
                                                 test_x)
{
  if (trtcoding == "minusplus"){
    trt[trt == -1] = 0
  }
  
  n_1 = NROW(trt[trt == 1])
  n_0 = NROW(trt[trt == 0])
  
  rowswith1 = which(trt == 1)
  rowswith0 = which(trt == 0)
  
  total_n <- NROW(x)
  
  tm <- "deviance"
  
  reorderedrowswith1 = sample(rowswith1)
  reorderedrowswith0 = sample(rowswith0)
  
  foldid = rep(NA, total_n)
  foldid[reorderedrowswith1] = sample(rep(seq(K), length = n_1))
  foldid[reorderedrowswith0] = sample(rep(seq(K), length = n_0))
  
  propensvec <- numeric(total_n)
  testpropensmat = matrix(data = NA, nrow = nrow(test_x), ncol = K)
  
  for (i in seq(K))
  {
    ithfold <- foldid == i
    
    ##################################### TODO ###############################
    ##################################### TODO ###############################
    ## FIX THE call to the rfsrc
    ##################################### TODO ###############################
    ##################################### TODO ###############################
    
    v.grow <- randomForestSRC::rfsrc(outcome ~ . , data = dat_all[-ithfold, ])
    
    cv.glmnet.args <- list(type.measure = tm, nfolds = glmnetfolds, grouped = FALSE)
    
    glmfit_propens <- do.call(glmnet::cv.glmnet, c(list(y = trt[!ithfold], x = x[!ithfold,,drop=FALSE],
                                                        family = "binomial"), cv.glmnet.args))
    
    ## get propensity scores for the held out fold
    propensvec[ithfold] <- unname(drop(stats::predict(glmfit_propens, newx = x[ithfold,,drop=FALSE],
                                                      s = "lambda.1se", type = "response")))
    
    ## get propensity scores for the test set
    testpropensmat[,i] = unname(drop(stats::predict(glmfit_propens, newx = test_x,
                                                    s = "lambda.1se", type = "response")))
  }
  
  ## propensity scores will never be outside of 0 or 1 and
  ## shouldn't have missing values, but this code is a safety
  ## check just in case
  propensvec[is.na(propensvec)] <- mean(propensvec[!is.na(propensvec)])
  propensvec[propensvec <= 0] <- 1e-5
  propensvec[propensvec >= 1] <- 1 - 1e-5
  
  testpropensvec = rowMeans(testpropensmat)
  
  
  return(list(propensvec=propensvec, testpropensvec = testpropensvec))
}