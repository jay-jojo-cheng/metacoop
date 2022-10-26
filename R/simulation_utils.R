gen_mvn_train_data <- function(numstudies, numobsperstudy, numcovariates){
  train_data_list <- vector(mode="list", length=numstudies)
  for (k in 1:numstudies){
    train_data_list[[k]] = cbind(1,matrix(rnorm((numobsperstudy[k])*numcovariates, mean=0, sd=1),
                                          (numobsperstudy[k]), numcovariates))
  }
  
  return(train_data_list)
}


gen_mvn_test_data <- function(numstudies, testsetsize, numcovariates){
  test_data_list <- vector(mode="list", length=numstudies)
  for (k in 1:numstudies){
    test_data_list[[k]] = cbind(1,matrix(rnorm(testsetsize*numcovariates, mean=0, sd=1),
                                         testsetsize, numcovariates))
  }
  
  return(test_data_list)
}

gen_propensity_scores <- function(propensity_sparsity, numcovariates, train_data_list){
  propensity_coefs = sample(c(-0.75,-0.5,-0.25,0.25,0.5,0.75), propensity_sparsity, replace = TRUE)
  propensity_coefs = c(0, propensity_coefs, rep(0, numcovariates - propensity_sparsity))
  
  train_propensity_list <- vector(mode="list", length=numstudies)
  for (k in 1:numstudies){
    propensity_coef_perturbed = propensity_coefs + c(0, rnorm(propensity_sparsity,mean=0,sd=0.1), rep(0, numcovariates - propensity_sparsity))
    #plogis is the logistic function
    train_propensity_list[[k]] = plogis(train_data_list[[k]] %*% propensity_coef_perturbed)
  }
  
  return(train_propensity_list)
}

gen_treatments <- function(numstudies, train_propensity_list){
  # generate treatments from Bernoulli trial of propensities
  # -1 is control, 1 is treatment
  train_treatment_list <- vector(mode="list", length=numstudies)
  for (k in 1:numstudies){
    train_treatment_list[[k]] = ifelse(
      train_propensity_list[[k]] > runif(length(train_propensity_list[[k]])),
      1, -1)
  }
  
  return(train_treatment_list)
}

gen_lin_main_effects <- function(numstudies, maineffect_sparsity, numcovariates){
  # generate base linear main effects
  # and study-specific linear main effects (perturbed from main effect coefs)
  maineffect_coef_list <- vector(mode="list", length=numstudies)
  
  # the signs will be fixed across studies
  maineffect_sign = sample(c(-1,1), maineffect_sparsity + 1, replace=TRUE)
  maineffect_sign = c(maineffect_sign, rep(0, numcovariates - maineffect_sparsity))
  
  # base + perturb, base is shared across studies, perturb is independent
  maineffect_coef = rlnorm(n=(maineffect_sparsity + 1), meanlog = 0)
  maineffect_coef = c(maineffect_coef, rep(0, numcovariates - maineffect_sparsity))
  
  for (k in 1:numstudies){
    maineffect_perturb = rlnorm(n=(maineffect_sparsity + 1), meanlog = 0)
    maineffect_perturb = c(maineffect_perturb, rep(0, numcovariates - maineffect_sparsity))
    maineffect_coef_list[[k]] = (maineffect_coef + maineffect_perturb)*maineffect_sign
  }
  
  return(maineffect_coef_list)
}

gen_lin_int_effects <- function(numstudies, inteffect_sparsity, numcovariates,
                                INTBASE, INTPERTURB){
  # generate treatment interaction effect per study (perturbed from interaction
  # effect coefs)
  inteffect_coef_list <- vector(mode="list", length=numstudies)
  
  # the signs will be fixed across studies
  inteffect_sign = sample(c(-1,1), inteffect_sparsity + 1, replace=TRUE)
  inteffect_sign = c(inteffect_sign, rep(0, numcovariates - inteffect_sparsity))
  
  # interaction is base + perturb, base is shared across studies
  # perturb is independent
  inteffect_coef = rlnorm(n=(inteffect_sparsity + 1), meanlog = INTBASE)
  inteffect_coef = c(inteffect_coef, rep(0, numcovariates - inteffect_sparsity))
  
  for (k in 1:numstudies){
    # generate treatment interaction effect coefficients
    # the signs are the same, but the sizes are different
    inteffect_perturb = rlnorm(n=(inteffect_sparsity + 1), meanlog = INTPERTURB)
    inteffect_perturb = c(inteffect_perturb, rep(0, numcovariates - inteffect_sparsity))
    
    inteffect_coef_list[[k]] = (inteffect_coef+inteffect_perturb)*inteffect_sign
  }
  
  return(inteffect_coef_list)
}

gen_nonlin_main_effects <- function(nonlinearmaineffect_sparsity,
                                    maineffect_sparsity, numstudies,
                                    train_data_list, test_data_list,
                                    nonlinear_ME_on){
  # generating nonlinear main effects
  nonlinearmaineffectsigns = sample(c(-1,1), nonlinearmaineffect_sparsity, replace = TRUE)
  nonlinearmaineffect = rlnorm(n=nonlinearmaineffect_sparsity, meanlog = 0) * nonlinearmaineffectsigns
  
  nonlinear_ME_index1 = sample(1:maineffect_sparsity, nonlinearmaineffect_sparsity, replace=TRUE)
  nonlinear_ME_index2 = sample(1:maineffect_sparsity, nonlinearmaineffect_sparsity, replace=TRUE)
  
  train_nonlinear_ME_list <- vector(mode="list", length=numstudies)
  test_nonlinear_ME_list <- vector(mode="list", length=numstudies)
  
  for (k in 1:numstudies){
    train_nonlinear_ME_list[[k]] = rep(0, dim(train_data_list[[k]])[1])
    test_nonlinear_ME_list[[k]] = rep(0, dim(test_data_list[[k]])[1])
    
    if (nonlinear_ME_on){
      for (j in 1:nonlinearmaineffect_sparsity){
        train_nonlinear_ME_list[[k]] = train_nonlinear_ME_list[[k]] +
          nonlinearmaineffect[j]*train_data_list[[k]][,nonlinear_ME_index1[j]]*train_data_list[[k]][,nonlinear_ME_index2[j]]
        
        test_nonlinear_ME_list[[k]] = test_nonlinear_ME_list[[k]] +
          nonlinearmaineffect[j]*test_data_list[[k]][,nonlinear_ME_index1[j]]*test_data_list[[k]][,nonlinear_ME_index2[j]]
      }
    }
    
  }
  
  return(list(train_nonlinear_ME_list = train_nonlinear_ME_list,
              test_nonlinear_ME_list = test_nonlinear_ME_list))
  
}

gen_nonlin_int_effects <- function(nonlinearinteffect_sparsity, inteffect_sparsity,
                                   numstudies, train_data_list,
                                   test_data_list, nonlinear_TE_on){
  nonlinearinteffectsigns = sample(c(-1,1), nonlinearinteffect_sparsity, replace = TRUE)
  nonlinearinteffect = rlnorm(n=nonlinearinteffect_sparsity, meanlog = 0) * nonlinearinteffectsigns
  
  nonlinear_TE_index1 = sample(1:inteffect_sparsity, nonlinearinteffect_sparsity, replace=TRUE)
  nonlinear_TE_index2 = sample(1:inteffect_sparsity, nonlinearinteffect_sparsity, replace=TRUE)
  
  # note that these are the actual additive terms to the effects and not coefficients
  
  train_nonlinear_TE_list <- vector(mode="list", length=numstudies)
  test_nonlinear_TE_list <- vector(mode="list", length=numstudies)
  
  for (k in 1:numstudies){
    train_nonlinear_TE_list[[k]] = rep(0, dim(train_data_list[[k]])[1])
    test_nonlinear_TE_list[[k]] = rep(0, dim(test_data_list[[k]])[1])
    
    if (nonlinear_TE_on){
      for (j in 1:nonlinearinteffect_sparsity){
        train_nonlinear_TE_list[[k]] = train_nonlinear_TE_list[[k]] +
          nonlinearinteffect[j]*train_data_list[[k]][,nonlinear_TE_index1[j]]*train_data_list[[k]][,nonlinear_TE_index2[j]]
        
        test_nonlinear_TE_list[[k]] = test_nonlinear_TE_list[[k]] +
          nonlinearinteffect[j]*test_data_list[[k]][,nonlinear_TE_index1[j]]*test_data_list[[k]][,nonlinear_TE_index2[j]]
      }
    }
    
  }
  
  return(list(train_nonlinear_TE_list = train_nonlinear_TE_list,
              test_nonlinear_TE_list = test_nonlinear_TE_list,
              nonlinearinteffect = nonlinearinteffect,
              nonlinear_TE_index1 = nonlinear_TE_index1,
              nonlinear_TE_index2 = nonlinear_TE_index2))
}

gen_train_outcome_list <- function(numstudies, train_data_list, maineffect_coef_list,
                                   train_treatment_list, inteffect_coef_list,
                                   train_nonlinear_ME_list, train_nonlinear_TE_list){
  # generate training outcomes by adding maineffect, nonlinear maineffect,
  # interaction effect, and error term
  train_outcome_list <- vector(mode="list", length=numstudies)
  for (k in 1:numstudies){
    train_outcome_list[[k]] =  train_data_list[[k]] %*% maineffect_coef_list[[k]] +
      train_treatment_list[[k]] * train_data_list[[k]] %*% inteffect_coef_list[[k]] +
      train_nonlinear_ME_list[[k]] + # turned off automatically when nonlinear_ME_on is FALSE
      train_treatment_list[[k]] * train_nonlinear_TE_list[[k]] + # turned off automatically when nonlinear_TE_on is FALSE
      rnorm(dim(train_data_list[[k]])[1], mean=0, sd=1)
  }
  
  return(train_outcome_list)
}

gen_test_maineffect <- function(numstudies, test_data_list,
                                maineffect_coef_list, test_nonlinear_ME_list){
  test_maineffect <- vector(mode = "list", length=numstudies)
  for (k in 1:numstudies){
    test_maineffect[[k]] = test_data_list[[k]] %*% maineffect_coef_list[[k]] +
      test_nonlinear_ME_list[[k]]
  }
  
  return(test_maineffect)
}

gen_test_inteffect <- function(numstudies, test_data_list, inteffect_coef_list,
                               test_nonliear_TE_list){
  test_interactioneffect <- vector(mode = "list", length=numstudies)
  for (k in 1:numstudies){
    test_interactioneffect[[k]] = test_data_list[[k]] %*% inteffect_coef_list[[k]] +
      test_nonlinear_TE_list[[k]]
  }
  
  return(test_interactioneffect)
}

gen_test_pot_outcomes <- function(numstudies, test_data_list,
                                  maineffect_coef_list, test_nonliear_ME_list,
                                  inteffect_coef_list, test_nonliear_TE_list){
  test_maineffect <- gen_test_maineffect(numstudies=numstudies,
                                         test_data_list=test_data_list,
                                         maineffect_coef_list=maineffect_coef_list,
                                         test_nonlinear_ME_list=test_nonlinear_ME_list)
  
  test_interactioneffect <- gen_test_inteffect(numstudies=numstudies,
                                               test_data_list=test_data_list,
                                               inteffect_coef_list=inteffect_coef_list,
                                               test_nonliear_TE_list=test_nonliear_TE_list)
  
  # generate test potential outcomes by adding maineffect, nonlinear maineffect,
  # interaction effect, and error term
  test_pot_outcome_pos_list <- vector(mode="list", length=numstudies)
  test_pot_outcome_neg_list <- vector(mode="list", length=numstudies)
  test_max_pot_outcome_list <- vector(mode="list", length=numstudies)
  for (k in 1:numstudies){
    test_pot_outcome_pos_list[[k]] =  test_maineffect[[k]] + test_interactioneffect[[k]] +
      rnorm(dim(test_data_list[[k]])[1], mean=0, sd=1)
    test_pot_outcome_neg_list[[k]] =  test_maineffect[[k]] - test_interactioneffect[[k]] +
      rnorm(dim(test_data_list[[k]])[1], mean=0, sd=1)
    
    # this is for computing oracle value function
    test_max_pot_outcome_list[[k]] = pmax(test_pot_outcome_pos_list[[k]],test_pot_outcome_neg_list[[k]])
  }
  
  return(list(test_pot_outcome_pos_list=test_pot_outcome_pos_list,
              test_pot_outcome_neg_list=test_pot_outcome_neg_list,
              test_max_pot_outcome_list=test_max_pot_outcome_list,
              test_interactioneffect=test_interactioneffect))
}