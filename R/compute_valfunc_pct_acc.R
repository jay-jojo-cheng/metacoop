compute_valfunc_pct_acc <- function(numstudies,numcovariates,test_data_list,
                                    simrep,beta,pos_pot_outc,neg_pot_outc,
                                    test_interactioneffect,
                                    oracle_value_function){
  
  valfunc = 0
  sumcorrect = 0
  for (k in 1:numstudies){
    startidx = ((k-1)*(numcovariates+1)+1)
    endidx = startidx + numcovariates
    
    scorefunc = test_data_list[[k]] %*% as.vector(beta[simrep,startidx:endidx])
    valfunc = valfunc + sum(pos_pot_outc[[k]][scorefunc >= 0]) + sum(neg_pot_outc[[k]][scorefunc < 0])
    
    sumcorrect = sumcorrect + sum(sign(scorefunc) == sign(test_interactioneffect[[k]]))
  }
  valfunc_pct = valfunc / abs(oracle_value_function)
  acc = sumcorrect / (length(unlist(pos_pot_outc)))
  
  return(list(valfunc_pct=valfunc_pct, acc=acc))
}

