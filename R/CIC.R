all_CICs_pick_lambda_beta <- function(metacoop_fit, data_list, numobsperstudy,
                                      numstudies, outcome_list,
                                      list_of_crossfit_augmentations,
                                      list_of_crossfit_propensities = NULL,
                                      list_of_weights = NULL,
                                      treatment_list){
  lambda_length = dim(metacoop_fit[[1]])[2]
  numcovariates = dim(data_list[[1]])[2] - 1
  
  concordance_per_lambda_study = matrix(data=NA,nrow=lambda_length,ncol=numstudies)
  # sign_concordance_per_lambda_study = matrix(data=NA,nrow=lambda_length,ncol=numstudies)
  # concordance_kendall_tau = matrix(data=NA,nrow=lambda_length,ncol=numstudies)
  betasize_per_lambda_study = matrix(data=NA,nrow=lambda_length,ncol=numstudies)
  
  # outer ell loop for lambda
  for (ell in 1:lambda_length){
    for (k in 1:numstudies){
      startidx = ((k-1)*(numcovariates+1)+1)
      endidx = startidx + numcovariates
      
      score_func = data_list[[k]] %*% metacoop_fit$return_beta[startidx:endidx,ell]
      
      # for CIC 1-4
      if (is.null(list_of_crossfit_propensities)){
        concordance_per_lambda_study[ell,k] = concordance_calc_ties(itr = score_func,
                                                                    y = outcome_list[[k]],
                                                                    mu_hat_x = list_of_crossfit_augmentations[[k]],
                                                                    weights = list_of_weights[[k]],
                                                                    trt = treatment_list[[k]],
                                                                    double.robust = FALSE,
                                                                    trt_coding = "minusplus",
                                                                    larger.outcome.better = FALSE)
      }else{
        concordance_per_lambda_study[ell,k] = concordance_calc_ties(itr = score_func,
                                                                    y = outcome_list[[k]],
                                                                    mu_hat_x = list_of_crossfit_augmentations[[k]],
                                                                    propens = list_of_crossfit_propensities[[k]],
                                                                    trt = treatment_list[[k]],
                                                                    double.robust = FALSE,
                                                                    trt_coding = "minusplus",
                                                                    larger.outcome.better = FALSE)
      }
      
      # for CIC 5-6
      # sign_concordance_per_lambda_study[ell,k] = concordance_calc_ties(sign(score_func),
      #                                                                  outcome_list[[k]],
      #                                                                  list_of_crossfit_augmentations[[k]],
      #                                                                  list_of_crossfit_propensities[[k]],
      #                                                                  treatment_list[[k]],
      #                                                                  double.robust=FALSE,
      #                                                                  trt_coding = "minusplus")
      
      betasize_per_lambda_study[ell,k] = sum(metacoop_fit$return_beta[startidx:endidx,ell] != 0)
      
      # # for CIC_kendall
      # concordance_kendall_tau[ell,k] = stats::cor(x=score_func,
      #                                      y= (outcome_list[[k]] - list_of_crossfit_augmentations[[k]]) * ((treatment_list[[k]]+1/2)-list_of_crossfit_propensities[[k]]) / (list_of_crossfit_propensities[[k]] * (1 - list_of_crossfit_propensities[[k]])),
      #                                      method = "kendall")
    }
  }
  
  # Compute CIC1
  # First version of CIC: sum_{k=1}^{K}.......
  # set this concordances to be the average of concordances
  avgconcordance_perlambda = rowSums(concordance_per_lambda_study) / numstudies
  
  # compute 0-norm of beta
  zeronorm_beta_perlambda = colSums(metacoop_fit$return_beta != 0)
  nprime = sqrt(sum(numobsperstudy^2))
  
  CIC1 = nprime*avgconcordance_perlambda - log(nprime)*zeronorm_beta_perlambda
  CIC1_best_lambda_idx = which.max(CIC1)
  CIC1_beta = metacoop_fit$return_beta[,CIC1_best_lambda_idx]
  
  # Compute CIC2
  # Second version of CIC: sum_{k=1}^{K}[n_{k}C_{k}(\beta) - log[n_{k}]||\beta||_{0}
  # compute 0-norm of beta_ks
  CIC2 = concordance_per_lambda_study %*% as.vector(numobsperstudy) - betasize_per_lambda_study %*% as.vector(log(numobsperstudy)) 
  CIC2_best_lambda_idx = which.max(CIC2)
  CIC2_beta = metacoop_fit$return_beta[,CIC2_best_lambda_idx]
  
  # Compute CIC3
  # Third version of CIC: sum_{k=1}^{K}[C_{k}(\beta)/K - log[n_{k}]||\beta||_{0}/(Kn_{k})
  # (scaled by K)
  CIC3 = rowSums(concordance_per_lambda_study) - betasize_per_lambda_study %*% as.vector(log(numobsperstudy)/numobsperstudy)
  CIC3_best_lambda_idx = which.max(CIC3)
  CIC3_beta = metacoop_fit$return_beta[,CIC3_best_lambda_idx]
  
  # Compute CIC4
  # Fourth version of CIC: sum_{k=1}^{K}[C_{k}(\beta) - 2||\beta||_{0}/(n_{k})]
  # (fixed kappa of 2)
  CIC4 = rowSums(concordance_per_lambda_study) - betasize_per_lambda_study %*% as.vector(2/numobsperstudy)
  CIC4_best_lambda_idx = which.max(CIC4)
  CIC4_beta = metacoop_fit$return_beta[,CIC4_best_lambda_idx]
  
  # # Compute CIC5
  # # Fifth version of CIC: same as CIC1 but with signconcordance
  # avg_signconcordance_perlambda = rowSums(sign_concordance_per_lambda_study) / numstudies
  # 
  # # compute 0-norm of beta
  # zeronorm_beta_perlambda = colSums(metacoop_fit$return_beta != 0)
  # nprime = sqrt(sum(numobsperstudy^2))
  # 
  # CIC5 = nprime*avg_signconcordance_perlambda - log(nprime)*zeronorm_beta_perlambda
  # CIC5_best_lambda_idx = which.max(CIC5)
  # CIC5_beta = metacoop_fit$return_beta[,CIC5_best_lambda_idx]
  # 
  # # Compute CIC6
  # # Sixth version of CIC: same as CIC4 but with signconcordance
  # CIC6 = rowSums(sign_concordance_per_lambda_study) - betasize_per_lambda_study %*% as.vector(2/numobsperstudy)
  # CIC6_best_lambda_idx = which.max(CIC6)
  # CIC6_beta = metacoop_fit$return_beta[,CIC6_best_lambda_idx]
  # 
  # # Compute CIC_kendall
  # CIC_kendall = concordance_kendall_tau %*% as.vector(numobsperstudy) - betasize_per_lambda_study %*% as.vector(log(numobsperstudy))
  # CIC_kendall_best_lambda_idx = which.max(CIC_kendall)
  # CIC_kendall_beta = metacoop_fit$return_beta[,CIC_kendall_best_lambda_idx]
  
  idx_list = c(CIC1_best_lambda_idx,
               CIC2_best_lambda_idx,
               CIC3_best_lambda_idx,
               CIC4_best_lambda_idx #,
               # CIC5_best_lambda_idx,
               # CIC6_best_lambda_idx,
               # CIC_kendall_best_lambda_idx
               )
  
  parameter_matrix = rbind(CIC1_beta,CIC2_beta,CIC3_beta,CIC4_beta
                           # ,CIC5_beta,CIC6_beta,CIC_kendall_beta
                           )
  
  
  return(list(parameters = parameter_matrix, lambda_idx = idx_list,
              CIC1 = CIC1,
              CIC2 = CIC2,
              CIC3 = CIC3,
              CIC4 = CIC4,
              # CIC5 = CIC5,
              # CIC6 = CIC6,
              # CIC_kendall = CIC_kendall,
              concordance_per_lambda_study = concordance_per_lambda_study,
              # concordance_kendall_tau = concordance_kendall_tau,
              numobsperstudy = numobsperstudy,
              betasize_per_lambda_study = betasize_per_lambda_study))
}

# return matrix of self-normalized CIC curves (first in x then in y)
self_normalized_CIC <- function(concordance_per_lambda_study, numobsperstudy,
                                betasize_per_lambda_study, scalingfactor = NULL){
  
  if(is.null(scalingfactor)){
    # scalingfactor = 1/seq(from = 1, to = 10, length.out = 1000)
    scalingfactor = 10^seq(from = 5, to = -5, length.out = 1000)
  }
  
  num_curves = length(scalingfactor)
  num_models = dim(concordance_per_lambda_study)[1]
  
  CIC_curves = matrix(data = NA, nrow = num_curves, ncol = num_models)
  for (i in 1:num_curves){
    CIC_curves[i,] = concordance_per_lambda_study %*%
      as.vector(numobsperstudy) -
      scalingfactor[i] *
      betasize_per_lambda_study %*%
      as.vector(log(numobsperstudy))
  }
  
  CIC_curves = self_norm_rows(CIC_curves)
  
  return(CIC_curves)
}

lambdadiff <- function(CIC_curves){
  num_curves = dim(CIC_curves)[1]
  num_models = dim(CIC_curves)[2]
  
  lambda_differences = matrix(data = NA, nrow = num_curves, ncol = (num_models - 1))
  for (i in 1:num_curves){
    lambda_differences[i,] = diff(CIC_curves[i,])
  }
  
  TVbyrow = rowSums(abs(lambda_differences))
  
  maxvariationrow_kappa = which.max(TVbyrow)
  
  maxCIC_lambda = which.max(CIC_curves[maxvariationrow_kappa, ])
  
  return(list(lambda_differences=lambda_differences,
              TVbyrow=TVbyrow,
              maxvariationrow_kappa=maxvariationrow_kappa,
              maxCIC_lambda=maxCIC_lambda))
}



# basically we're using total variation (normalized) now to compare the curves

# new method:
# take the lambda differences (x direction) and sum up the absolute values.
# take the curve with largest total variation distance

# Take a rolling average of total variation distances e.g. take the row abs lambda diff
# sum and sum up the ones on either side.