energy_balance <- function(x,  ## matrix of covariates for all RCT and non-RCT patients
                           trt, ## vector of treatment indicators, 1 for each patient with s=1,
                           ## in the order that the 1's appear in the vector 's'
                           standardize_x = FALSE, ## standardize the covariates?
                           verbose = FALSE,
                           sigma = 0.5,
                           solver = c("osqp", "cccp"),
                           type = c("energy",
                                    "gaussian_kernel")) ## what kernel/distance to use?
{
  if (standardize_x)
  {
    x <- scale(x)
  }
  
  p <- NCOL(x)
  
  
  type    <- match.arg(type)
  solver  <- match.arg(solver)
  
  dist_x  <- as.matrix(dist(x))
  
  median_dist <- median(dist_x)
  
  
  ## we will optimize weighted energy distance
  ## by minimizing f(w) = w'Qw + a'w subject to constraints on w.
  ## constraints are that w >= 0 and Aw = b, where b is a vector
  
  ## the following code constructs the Q matrix and 'a' vector
  ## corresponding to the weighted energy distance.
  
  ## if the treatment is not randomized, we will
  ## balance the treatment arm and the control arm separately
  ## to the target population (instead of the combined trt+control sample).
  ## this makes Q and 'a' slightly different.
  
  
  N <- nrow(x)
  
  if (grepl("gaussian_kernel", type))
  {
    Q_x   <- -exp(-sigma * dist_x ^ 2 / median_dist ^ 2)
  } else
  {
    Q_x   <- dist_x
  }
  
  
  Q_mat <- -trt * t( trt * t(Q_x)) / sum(trt == 1) ^ 2 -
    (1 - trt) * t( (1 - trt) * t(Q_x)) / sum(trt != 1) ^ 2
  
  aa <- 2 * as.vector(rowSums(trt * Q_x)) / (sum(trt == 1) * N) +
    2 * as.vector(rowSums((1-trt) * Q_x)) / (sum(trt != 1) * N)
  
  
  #evals <- eigen(Q_mat)$values
  
  ## only add nugget if we plan to use the solve.QP() function
  # if (min(evals) < 0 & (grepl("gaussian_kernel", type) | grepl("gaussian_anova_kernel", type)))
  # {
  #   nugget <- -min(evals) + 1e-10
  # } else
  # {
  #   nugget <- 0
  # }
  
  nugget <- 0
  
  ## add the smallest value to diagonal such that
  ## the Q matrix is positive definite
  Q_mat <- Q_mat + (nugget) * diag(nrow(Q_mat))
  
  
  ## set up constraints on weights
  AA1           <- matrix(0, ncol = N, nrow = 1)
  AA0           <- matrix(0, ncol = N, nrow = 1)
  
  AA1[,trt == 1] <- 1
  AA0[,trt != 1] <- 1
  
  
  A           <- rbind(AA1, AA0)
  rownames(A) <- paste0("eq", 1:nrow(A))
  sum.constr   <- c(sum(trt == 1), sum(trt != 1))
  
  
  ## A*w = sum.constr forces the weights
  ## to sum to the sample size of each treatment group.
  ## ie w_i : trt_i = 1 will sum to the number of trt == 1
  
  
  if (solver == "osqp")
  {
    ## minimizing f(w) = 0.5*w'Qw - a'w
    
    ## we will set up the quadratic program to be solved by
    ## the solve.QP function accessed from quadprog. this solver
    ## is very fast but cannot handle negative definite Q_mat
    
    
    ## set up the linear (in)equality constraints. The first 1 or 2
    ## constraints are the equality constraints. the remaining are the
    ## inequality constraints (ie making sure all weights are positive)
    
    #Amat <- t( rbind(A, diag(nrow(Q_mat)), -diag(nrow(Q_mat)) ) )
    
    Amat <- rbind(diag(N), AA1, AA0)
    
    lvec <- c(rep(0, N), sum.constr)
    uvec <- c(rep(Inf, N), sum.constr)
    
    ## Amat * w (>)= bvec
    #bvec <- c(sum.constr, rep(1e-14, nrow(Q_mat)), rep(-10*(nrow(Q_mat) ^ (1/3)), nrow(Q_mat)))
    
    #res <- solve.QP(Dmat = Q_mat, dvec = -aa/2, Amat = Amat,
    #                bvec = bvec, meq = length(sum.constr))
    
    
    lagrangian <- function(x, lambda)
    {
      print("aa")
      print(str(aa))
      print("lam")
      print(str(lambda))
      print("A")
      print(str(Amat))
      
      term_1 <- 0.5 * t(x) %*% Q_mat %*% x + (t(aa)/2 - t(lambda[c(1:N)]) %*% Amat[c(1:N),]) %*% x
      term_2 <- (t(lambda[c(1:N)]) %*% lvec[c(1:N)])  # - (t(lambda[-c(1:n)]) %*% uvec[-c(1:n)])
      term_3 <- t(lambda[-c(1:N)]) %*% Amat[-c(1:N),,drop=FALSE] %*% x
      
      term_1 + term_2 + term_3
    }
    
    opt.out <- osqp::solve_osqp(Q_mat,
                                q = aa/2,
                                A = Amat, l = lvec, u = uvec,
                                pars = osqp::osqpSettings(max_iter = 2e5,
                                                          eps_abs = 1e-8,
                                                          eps_rel = 1e-8,
                                                          verbose = FALSE))
    
    # if (!identical(opt.out$info$status, "maximum iterations reached") & !(any(opt.out$x > 1e5)))
    # {
    #   break
    # }
    
    energy_wts <- unname(opt.out$x)
    energy_wts[energy_wts < 0] <- 0 #due to numerical imprecision
    
  } else
  {
    ## minimizing f(w) = w'Qw + a'w
    
    ## we will set up the quadratic program to be solved by
    ## the cccp function accessed from optiSolve. this solver
    ## can handle when Q_mat is negative definite very well!
    
    rownames(Q_mat) <- paste(1:NROW(Q_mat))
    
    qf <- optiSolve::quadfun(Q = Q_mat, a = aa, id = rownames(Q_mat)) #quadratic obj.
    
    lb <- optiSolve::lbcon(0.0, id = rownames(Q_mat)) #lower bound
    
    ub <- optiSolve::ubcon(10*(nrow(Q_mat) ^ (1/3)), id = rownames(Q_mat)) #upper bound. not really needed
    
    ## sup up linear constraints
    lc <- optiSolve::lincon(A   = A,
                 dir = rep("==", length(sum.constr)),
                 val = sum.constr,
                 id  = rownames(Q_mat)) #linear constraint
    
    ## set up the QP
    lcqp <- optiSolve::cop( f = qf, lb = lb, ub = ub, lc = lc )
    
    opt.out <- optiSolve::solvecop(lcqp, solver = "cccp", quiet = !verbose)
    
    energy_wts <- unname(opt.out$x)
    
    lagrangian <- NULL
    
  }
  
  energy_wts[energy_wts < 0] <- 0
  
  list(weights = energy_wts, opt = opt.out, lagrangian = lagrangian)
}