// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "metacoop.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::Map;

// constructor
metacoop::metacoop(Rcpp::List Xmat_list,
                   Rcpp::List Y_list,
                   const Eigen::Ref<const Eigen::VectorXi> & NOBSPERSTUDY_,
                   Rcpp::List weights_list,
                   const Eigen::Ref<const Eigen::VectorXd> & LAMBDA_,
                   const int & totalnvars,
                   const int & nlambda,
                   const int & nvars,
                   const int & nstudies,
                   const double & tol,
                   const double & lambda_min_ratio) 
  : // the form below is how to construct a Map Matrix variable
    // you need a pointer to the region of memory of the coefficients
    // and the desired shape (rows and cols)
    // X(X_.data(),X_.rows(),X_.cols()),
    
    // the form below is how to construct a Map Vector variable
    // you need a pointer to the region of memory of the coefficients
    // and the desired size, note the matrix style constructor would
    // also work since Eigen vectors are Eigen matrices.
    // Y(Y_.data(),Y_.size()),
    nobs_per_study(NOBSPERSTUDY_.data(),NOBSPERSTUDY_.size()),
    // weights(WEIGHTS_.data(),WEIGHTS_.size()),
    // will make a way for user to input lambda sequence in future
    lambda_OLD(LAMBDA_.data(),LAMBDA_.size()),
    nvars(nvars),
    totalnvars(totalnvars),
    nstudies(nstudies),
    nlambda(nlambda),
    tol(tol),
    lambda_min_ratio(lambda_min_ratio),
    
    // strongrule datamembers
    ever_active_set(VectorXd::Zero(totalnvars)),
    strong_set(VectorXd::Zero(totalnvars)),
    eligible_predictor_set(VectorXd::Zero(totalnvars)),
    kkt_violators(VectorXd::Zero(totalnvars)),
    all_set(VectorXd::Ones(totalnvars)),
    eligible_predictor_mask(MatrixXd::Identity(totalnvars,totalnvars)),
    
    // datamembers for algorithm
    // beta and beta_old are nvars*nstudies long vectors
    beta(VectorXd::Zero(totalnvars)),
    beta_old(VectorXd::Zero(totalnvars)),
    
    // gammajs is nvars long vector
    gammajs(VectorXd::Zero(nvars)),
    Ujs(MatrixXd::Zero(nstudies,nvars)),
    
    // return data
    return_beta(MatrixXd::Zero(totalnvars, nlambda))
    {
  //Rcpp::Rcout << "runtime_constructor_debug" << std::endl;

  // create list of X matrices - checked
  for (int k=0; k < nstudies; k++) {
    Rcpp::NumericMatrix tmpcv = Xmat_list[k];
    double *pcv = &tmpcv(0,0);

    Eigen::Map<Eigen::MatrixXd> tmpmapd(pcv, tmpcv.nrow(), tmpcv.ncol());
    vec_of_Xmats.push_back(tmpmapd);
  }
  // Rcpp::Rcout << "got here a" << std::endl;

  // create list of weights - checked
  for (int k=0; k < nstudies; k++) {
    Rcpp::NumericVector tmpcv = weights_list[k];
    Eigen::Map<Eigen::VectorXd> tmpmapd(Rcpp::as<Map<VectorXd> >(tmpcv));
    vec_of_weights.push_back(tmpmapd);
  }

  // Rcpp::Rcout << "got here b" << std::endl;

  // create list of wX matrices - checked
  for (int k=0; k < nstudies; k++) {
    MatrixXd tmpmat = vec_of_weights[k].asDiagonal() * vec_of_Xmats[k];
    vec_of_wXs.push_back(tmpmat);
  }

  // create list of Ys - checked
  for (int k=0; k < nstudies; k++) {
    Rcpp::NumericVector tmpcv = Y_list[k];
    Eigen::Map<Eigen::VectorXd> tmpmapd(Rcpp::as<Map<VectorXd> >(tmpcv));
    vec_of_Ys.push_back(tmpmapd);
  }

  // initialize residuals OLD
  // residuals = Y-(X * beta);
  // initialize residuals xlist version
  vec_of_residuals = vec_of_Ys;
  // Rcpp::Rcout << "initial residuals[0] = " << vec_of_residuals[0] << std::endl;

  // initialize Ujs - checked
  for (int k=0; k < nstudies; k++){
    // factor of 2?
    Ujs.row(k) = -1.0 * vec_of_Ys[k].transpose() * vec_of_weights[k].asDiagonal() * vec_of_Xmats[k] / nobs_per_study(k);
  }
  // Rcpp::Rcout << "initial Ujs " << Ujs << std::endl;

  // initialize ever_active_set to start with main effects
  for (int k = 0; k < nstudies; ++k){
    ever_active_set(k*nvars) = 1;
  }
  // Rcpp::Rcout << "starting ever_active_set make sure main eff is 1" << ever_active_set << std::endl;
  
  
  // initialize lambdamax
  lambdamax = compute_lambdamax();

  // initialize lambda sequence
  double lambdamin = lambda_min_ratio*lambdamax;
  double logmax = std::log10(lambdamax);
  double logmin = std::log10(lambdamin);
  double increment = (logmax - logmin) / (nlambda - 1);
  lambda = VectorXd::Zero(nlambda);
  for (int ell = 0; ell < nlambda; ell++){
    lambda(ell) = std::pow(10, logmax - ell*increment);
  }
  
  // Rcpp::Rcout << "lambda is" << lambda << std::endl;
  
  // Rcpp::Rcout << "LOOKINGFORTHISPRINT" << std::endl;
}

// Algorithm 1 with x and weight lists and strong rule
void metacoop::fit_path_xlist() {
  // Rcpp::Rcout << "debug 0" << std::endl;
  // update gammaj just once in our method
  gammajs = update_gammajs_xlist();

  // loop over lambdas
  for (int i = 0; i < nlambda; ++i) {
    // Rcpp::Rcout << "debug 1" << std::endl;
    // set current lambda
    lambdacur = lambda[i];
    
    // updated 2020/11/17
    update_strong_set(i);

    // Rcpp::Rcout << "debug 2" << std::endl;
    // strong outer loop iteration counter
    int itercount_strong_outer = 0;
    do {
      // Rcpp::Rcout << "debug 3" << std::endl;
      // add KKT violators to ever_active_set
      for(int j=0; j < totalnvars; j++){
        // Rcpp::Rcout << "debug 4" << std::endl;
        if (kkt_violators(j) && (ever_active_set(j) < 1) ) {
          // Rcpp::Rcout << "debug 5" << std::endl;
          ever_active_set(j) = 1;
        }
      }
      
      // strong step 1. set E = everactiveset(lambda)
      eligible_predictor_set = ever_active_set;

      // Rcpp::Rcout << "debug 6" << std::endl;

      // strong step 3. while there exists Strong KKT violators
      // add them into E and go back to step 1 using current soln as warm start
      int itercount_strong_inner = 0;
      do {
        // Rcpp::Rcout << "debug 7" << std::endl;

        // add strong violators to eligible set from prev
        // add strong violators
        for(int j=0; j < totalnvars; j++){
          // Rcpp::Rcout << "debug 8" << std::endl;
          
          if (kkt_violators(j) && (eligible_predictor_set(j) < 1)) {
            // Rcpp::Rcout << "debug 9" << std::endl;
            
            eligible_predictor_set(j) = 1;
          }
        }

        // strong step 2. solve problem at lambda using only E
        // Rcpp::Rcout << "debug 10" << std::endl;

        // don't need masking for this version
        // update_eligible_predictor_mask(eligible_predictor_set);

        // Rcpp::Rcout << "debug 11" << std::endl;

        // Rcpp::Rcout << eligible_predictor_set << std::endl;
        fit_beta2(eligible_predictor_set);

        // Rcpp::Rcout << "debug 12" << std::endl;

        itercount_strong_inner++;
        //TODO fix the checking set
        // while there is a kkt violator in strong set \ eligible predictor set
        // add into eligible_predictor_set
      }while( exists_kkt_violators( set_minus(strong_set, eligible_predictor_set) , eligible_predictor_set) );
      
      itercount_strong_outer++;
      //TODO fix the checking set
      // while there exists kkt violator for all predictor \ ()eligiblepredictor union strong set)
      // add into ever_active_set
    }while( exists_kkt_violators( set_minus(all_set, eligible_predictor_set.array().max(strong_set.array())) , ever_active_set) );

    // save the current betas to the corresponding column of return_beta
    return_beta.col(i) = beta;
  }
}

// I think phi_m_cpp should be a global function, I think this is faster?
// maybe need to try it as a member function too and benchmark it.
Eigen::VectorXd phi_m_cpp(const Eigen::Ref<const Eigen::VectorXd> & v, const int & m) {
  // const Map<MatrixXd> v(as<Map<MatrixXd> >(VV));
  // int m = Rcpp::as<int >(MM);
  
  int vlen = v.size();
  VectorXd retvec(vlen);
  retvec.setZero();

  if (v(m) > 0.0)
  {
    // calculate v plus
    for (int k = 0; k < vlen; ++k)
    {
      retvec(k) = std::max(0.0, v(k));
    }
  } else if (v(m) < 0.0)
  {
    // calculate v minus
    for (int k = 0; k < vlen; ++k)
    {
      retvec(k) = std::max(0.0, -1.0 * v(k));
    }
  }
  return(retvec);
}

bool metacoop::converged(const Eigen::VectorXd & cur, const Eigen::VectorXd & prev, const double & tolerance) {
  for (int i = 0; i < cur.rows(); i++)
  {
    if ( (std::abs(cur(i)) > 1e-13 && std::abs(prev(i)) <= 1e-13) ||
         (std::abs(cur(i)) <= 1e-13 && std::abs(prev(i)) > 1e-13) ) {
      return 0;
    }
    if (std::abs(cur(i)) > 1e-13 && std::abs(prev(i)) > 1e-13 &&
        std::pow( (cur(i) - prev(i)) / prev(i), 2) > tolerance) {
      return 0;
    }
  }
  return 1;
}

// Eigen::VectorXd metacoop::takeBreturnB(const Eigen::Ref<const Eigen::VectorXd> & B) {
//   return B;
// }

// Eigen::VectorXd metacoop::xbeta_cpp(const Eigen::Ref<const Eigen::MatrixXd> & X, const Eigen::Ref<const Eigen::VectorXd> & B) {
//   // const Map<MatrixXd> X(as<Map<MatrixXd> >(XX));
//   // const Map<VectorXd> B(as<Map<VectorXd> >(BB));
// 
//   VectorXd XB = X * B;
// 
//   return XB;
// }

// DONE
Eigen::MatrixXd metacoop::update_gammajs_xlist() {
  // create matrix for comparison
  MatrixXd Hjmatrix(MatrixXd::Zero(nvars, nstudies));
  for (int k = 0; k < nstudies; ++k) {
    // compute component-wise square
    MatrixXd retvec_temp1 = vec_of_Xmats[k].array().square();
    
    // matrix multiplication, get a row vector with
    // elements that are the blocks in S.18b
    Hjmatrix.col(k) = vec_of_weights[k].transpose() * retvec_temp1 / nobs_per_study(k);
  }
  // then take row-wise (max over each study) max
  MatrixXd retvec = Hjmatrix.rowwise().maxCoeff();
  
  // we need to transpose it to get the right form
  return(retvec);
}

// DONE needs a test
void metacoop::update_residuals_xlist(const int & j) {
  for(int k=0; k < nstudies; k++){
    // residual update is Xj column * (betaoldj-betanewj) scalar
    double beta_dif = beta_old((k)*nvars + j) - beta((k)*nvars + j);
    // Rcpp::Rcout << "beta_old-beta_new in update_resids is " << beta_dif << std::endl;
    vec_of_residuals[k] += vec_of_Xmats[k].col(j) * beta_dif;
    // Rcpp::Rcout << "WE START PRINTING IN THE RESIDUAL UPDATE FN" << std::endl;
    // 
    // Rcpp::Rcout << "beta_old" << beta_old((k)*nvars + j) << std::endl;
    // Rcpp::Rcout << "beta_new" << beta((k)*nvars + j) << std::endl;
    // Rcpp::Rcout << "vec_of_Xmats[k].col(j)" << vec_of_Xmats[k].col(j) << std::endl;
    // Rcpp::Rcout << "residuals[0] = " << vec_of_residuals[k] << std::endl;
    // 
    // Rcpp::Rcout << "WE ARE DONE PRINTING IN THE RESIDUAL UPDATE FN" << std::endl;
  }
}

// DONE needs a test
void metacoop::update_Ujs_cpp_xlist(const int & j) {
  for(int k=0; k < nstudies; k++){
    
    // factor of 2?
    Ujs(k,j) = -1.0 * (vec_of_residuals[k].array() * vec_of_wXs[k].col(j).array()).sum() / nobs_per_study(k);
    // update all js at once
    // Ujs.row(k) = -2.0 * (vec_of_wXs[k] * vec_of_residuals[k].asDiagonal()).colwise().sum() / nobs_per_study(k);
    
    // old ones that don't work
    // Ujs.row(k) = -2.0 * (vec_of_weights[k].array() * vec_of_Xmats[k].col(j).array() * vec_of_residuals[k].array()).sum() / nobs_per_study(k);
    // Ujs(k,j) -= 2.0 * vec_of_wXtXs[k](j,j) * beta_dif / nobs_per_study(k);
  }
}

double metacoop::update_Bjs_cpp(const int & j, const int & m) {
  // we assume gj is 1 according to p18
  
  // pick the jth gamma scalar from gammaj
  double curgammaj = gammajs[j];
  
  // create the betaj subvector
  // the following takes the jth slot of beta, and then takes mutliples of nvars,
  // shifted by j to get the jth variable of each of the nstudies.
  Eigen::Map<Eigen::VectorXd,0,Eigen::InnerStride<> > betaj_sub(beta.segment(j,totalnvars - j).data(),nstudies,Eigen::InnerStride<>(nvars));
  
  // pick the Uj subvector
  VectorXd Uj_sub = Ujs.col(j);
  
  double denom = phi_m_cpp(curgammaj * betaj_sub - Uj_sub, m).norm();
  
  double rightterm = std::max(0.0, 1 - lambdacur * 1 / denom);
  
  double leftterm = betaj_sub(m) - Uj_sub(m)/curgammaj;
  
  // Rcpp::Rcout << "beta-U/gamma is " << leftterm << std::endl;
  // Rcpp::Rcout << "max(0,1-blob) is" << rightterm << std::endl;
  
  return leftterm * rightterm;
}

double metacoop::compute_lambdamax() {
  double max_val = 0;
  
  VectorXd XgroupTy = VectorXd::Zero(nstudies);
  for (int j=1; j < nvars; j++){
    for(int k=0; k < nstudies; k++){
      // todo: this computation is repeated, we can just save a running copy
      XgroupTy(k) = vec_of_Xmats[k].col(j).adjoint()*vec_of_Ys[k];
    }
    
    for(int k=0; k < nstudies; k++){
      // old without scaling / I_k 
      // max_val = std::max(max_val, phi_m_cpp(XgroupTy,k).norm());
      // new with scaling / I_k
      max_val = std::max(max_val, phi_m_cpp(XgroupTy,k).norm()/nobs_per_study(k));
    }
  }
  
  return(max_val);
}

void metacoop::update_strong_set(int lambda_index) {
  // Rcpp::Rcout << "update_strong_debug_2" << std::endl;
  double rhs;
  if (lambda_index == 0) {
    rhs = 2*lambda(lambda_index) - lambdamax;
  } else {
    rhs = 2*lambda(lambda_index) - lambda(lambda_index-1);
  }
  
  double lhs;
  
  // Rcpp::Rcout << "update_strong_debug_3" << std::endl;
  for(int j=1; j < nvars; j++){
    VectorXd XjgroupTresid = VectorXd::Zero(nstudies);
    for(int k=0; k < nstudies; k++){
      XjgroupTresid(k) = vec_of_Xmats[k].col(j).adjoint()*vec_of_residuals[k];
    }
    
    // comparison of each strongrule for jgroup
    for(int k=0; k < nstudies; k++){
      lhs = phi_m_cpp(XjgroupTresid, k).norm();
      if (lhs < rhs) {
        strong_set(k*nvars + j)=0;
      } else {
        strong_set(k*nvars + j)=1;
      }
    }
  }
}

// returns complement e.g. [0,1,0,0] -> [1,0,1,1]
Eigen::VectorXd metacoop::complement_set(Eigen::VectorXd set_to_check) {
  Eigen::VectorXd comp_set = VectorXd::Ones(totalnvars);
  
  for (int j = 0; j < totalnvars; j++){
    if (set_to_check(j)){
      comp_set(j) = 0;
    }
  }
  
  return comp_set;
}

Eigen::VectorXd metacoop::set_minus(Eigen::VectorXd set_to_subtract_from, Eigen::VectorXd set_to_subtract) {
  Eigen::VectorXd result_set = (set_to_subtract_from - set_to_subtract).array().max(0);
  
  return result_set;
}

// updated 11/18 with set_to_add, may need to move xjgroupTresid into a data member
bool metacoop::exists_kkt_violators(Eigen::VectorXd set_to_check, Eigen::VectorXd set_to_add) {
  bool exists = false;
  kkt_violators = VectorXd::Zero(totalnvars);
  
  // add the violator set
  double lhs;
  
  // Rcpp::Rcout << "exists_kkt_violators_debug_3" << std::endl;
  // note that we never need to check the main effects
  for(int j=1; j < nvars; j++){
    VectorXd XjgroupTresid = VectorXd::Zero(nstudies);
    for(int k=0; k < nstudies; k++){
      // todo: this computation is repeated, we can just save a running copy
      XjgroupTresid(k) = vec_of_Xmats[k].col(j).adjoint()*vec_of_residuals[k];
    }
    
    // check kkt
    for(int k=0; k < nstudies; k++){
      if(set_to_check(k*nvars + j)){
        lhs = phi_m_cpp(XjgroupTresid, k).norm();
        if (lhs > lambdacur) {
          kkt_violators(k*nvars + j) = 1;
          exists = true;
        }
      }
      
    }
  }
  
  // add violators into set
  set_to_add = set_to_add.array().max(kkt_violators.array());
  
  return exists;
}


bool metacoop::approx_zero(double value) {
  if (std::abs(value) < tol) {
    return true;
  }
  else {
    return false;
    }
}

void metacoop::fit_beta2(Eigen::VectorXd eligible_predictors) {
  // iteration counter for debugging
  int itercount_beta = 0;
  
  // while beta not converged, keep iterating
  // limiting iters to debug:
  // while( itercount < 3) {
  // original:
  while( !converged(beta, beta_old, tol) || itercount_beta < 10) {
    // Can turn off this print line when not needed
    // Rcpp::Rcout << "itercount = " << itercount << std::endl;
    // Rcpp::Rcout << "beta = " << beta << std::endl;
    // Rcpp::Rcout << "Uj = " << Ujs << std::endl;
    
    // update each beta_jm (all elements of beta)
    for (int j = 0; j < nvars; ++j) {
      if( eligible_predictors(j)) {
        update_Ujs_cpp_xlist(j);
        beta_old = beta;
        // replace with nstudies
        for (int k = 0; k < nstudies; ++k) {
          beta((k)*nvars + j) = update_Bjs_cpp(j,k);
        }
        update_residuals_xlist(j);
      }
    }
    // Can turn off this print line when not needed
    // Rcpp::Rcout << "itercount_beta = " << itercount_beta << std::endl;
    itercount_beta++;
  }
}