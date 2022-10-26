concordance_calc_ties <- function(itr, y, mu_hat_x, propens = NULL,
                                  weights = NULL, trt,
                                  double.robust=FALSE,
                                  trt_coding = "minusplus",
                                  larger.outcome.better = TRUE) {
  # Original function by Guanhua Chen
  # Functionality for accommodating ties in 'itr',
  # weights instead of propensities, and different preferred
  # outcome direction by J. Jojo Cheng
  
  # switch signs if smaller outcome is better
  if (larger.outcome.better == FALSE){
    itr = -itr
    y = -y
    mu_hat_x = -mu_hat_x
  }
  
  # Fix treatment notation
  if(trt_coding == "minusplus"){
    # it takes -1,1 treatment coding and changes it to 0,1 for this func.
    trt = (trt+1)/2
  }
  else if(trt_coding == "zeroone"){
    # we already have 0,1 trt coding - do nothing
  }
  
  n       <- length(y)
  y.resid <- y - mu_hat_x
  
  # use weights if supplied directly, otherwise use the propensity score
  # for weighting
  if (is.null(weights)){
    weight <- (trt - propens) / (propens * (1 - propens))
  } else{
    weight <- weights
  }
  
  values <- y.resid * weight
  
  # we sort by ITR small to big
  order_idx = order(itr)
  itr_ordered = itr[order_idx]
  values_ordered <- values[order_idx]
  
  # each value is larger than minrank - 1 values (pos contribution)
  # and smaller than n - maxrank values (neg contribution)
  # note that rank returns 1 for the smallest value, etc.
  max_itr_rank = rank(itr_ordered, ties.method = "max")
  min_itr_rank = rank(itr_ordered, ties.method = "min")
  
  if(double.robust == FALSE) {
    # we add the number of pos terms and subtract the number of neg terms
    numberofterms = (min_itr_rank - 1) - (n - max_itr_rank)
    conc.contribution = numberofterms * values_ordered
  }
  else if(double.robust == TRUE) {
    dr.term              <- trt / propens
    dr.term_ordered <- dr.term[order_idx] 
    
    # the ith element is the sum of i smallest dr values
    cumsum.dr.term = cumsum(dr.term_ordered)
    # the jth element is the sum of j largest dr values
    rev_cumsum.dr.term = cumsum(rev(dr.term_ordered))
    
    # we explicitly create an element for zero contribution
    # so vec[0] gives 0 and vec[i+1] gives the sum of i smallest (or largest)
    cumsum.dr.term_zeroidx = c(0, cumsum.dr.term)
    rev_cumsum.dr.term_zeroidx = c(0, rev_cumsum.dr.term)
    
    # pos contribution is sum of dr terms from smaller itrs,
    # which is cumsum.dr.term_zeroidx[(n - max_itr_rank + 1)]
    # neg contribution is sum of dr terms from larger itrs,
    # which is rev_cumsum.dr.term_zeroidx[(min_itr_rank - 1 + 1)]
    dr_scaling = cumsum.dr.term_zeroidx[(n - max_itr_rank + 1)] -
      rev_cumsum.dr.term_zeroidx[(min_itr_rank - 1 + 1)]
    conc.contribution = dr_scaling * values_ordered
  }
  
  concordance <- (1/(n * (n-1))) * sum(conc.contribution)
  return(concordance)
}