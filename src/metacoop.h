#ifndef GUARD_metacoop
#define GUARD_metacoop

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>


// [[Rcpp::depends(RcppEigen)]]

Eigen::VectorXd phi_m_cpp(const Eigen::Ref<const Eigen::VectorXd> & v, const int & m);

class metacoop {
public:
  // constructor
  metacoop(Rcpp::List Xmat_list,
           Rcpp::List Y_list,
           const Eigen::Ref<const Eigen::VectorXi> & NOBSPERSTUDY_,
           Rcpp::List weights_list,
           const Eigen::Ref<const Eigen::VectorXd> & LAMBDA_,
           const int & totalnvars,
           const int & nlambda,
           const int & nvars,
           const int & nstudies,
           const double & tol,
           const double & lambda_min_ratio);
  
  // fit methods
  // void fit_path();
  // void fit_path_no_SR();
  void fit_path_xlist();
  // void fit_path_xlist_no_SR();
  
  // subroutines
  // Eigen::VectorXd takeBreturnB(const Eigen::Ref<const Eigen::VectorXd> & B);
  // Eigen::VectorXd xbeta_cpp(const Eigen::Ref<const Eigen::MatrixXd> & X, const Eigen::Ref<const Eigen::VectorXd> & B);
  Eigen::MatrixXd update_gammajs_cpp();
  Eigen::MatrixXd update_Ujs_cpp();
  double update_Bjs_cpp(const int & j, const int & m);
  bool converged(const Eigen::VectorXd & cur, const Eigen::VectorXd & prev, const double & tolerance);
  
  Eigen::MatrixXd update_gammajs_xlist();
  void update_residuals_xlist(const int & j);
  void update_Ujs_cpp_xlist(const int & j);
  
  // strongrule methods
  double compute_lambdamax();
  void update_strong_set(int lambda_index);
  // void update_eligible_predictor_mask(Eigen::VectorXd eligible_predictors);
  bool exists_kkt_violators(Eigen::VectorXd set_to_check, Eigen::VectorXd set_to_add);
  // bool fails_kkt(int totalnvar_idx);
  // bool fails_kkt_xlist(int totalnvar_idx);
  bool approx_zero(double value);
  // void fit_beta(Eigen::VectorXd eligible_predictor_set);
  void fit_beta2(Eigen::VectorXd eligible_predictor_set);
  
  // strongrule set utilities
  Eigen::VectorXd set_minus(Eigen::VectorXd set_to_subtract_from, Eigen::VectorXd set_to_subtract);
  Eigen::VectorXd complement_set(Eigen::VectorXd set_to_check);
  
  // changing Alg1 data
  Eigen::VectorXd beta;
  Eigen::VectorXd beta_old;
  // Eigen::VectorXd residuals;
  Eigen::VectorXd gammajs;
  Eigen::VectorXd xbeta;
  Eigen::MatrixXd Ujs;
  double lambdacur;
  double lambdaprev;
  
  // changing strongrule data
  // these have dimension determined by totalnvars (not nvars)
  Eigen::VectorXd ever_active_set;
  Eigen::VectorXd strong_set;
  Eigen::VectorXd eligible_predictor_set;
  Eigen::MatrixXd eligible_predictor_mask;
  Eigen::VectorXd kkt_violators;
  
  // return data
  Eigen::MatrixXd return_beta;
  
  // temp for debugging - move back to protected
  const int nvars;
  const int nstudies;
  const int totalnvars;
  const int nlambda;
  
  // lambda sequence data
  double lambdamax;
  Eigen::VectorXd lambda;
  
  // fixed xlist data
  std::vector<Eigen::MatrixXd> vec_of_Xmats;
  std::vector<Eigen::VectorXd> vec_of_weights;
  std::vector<Eigen::MatrixXd> vec_of_wXs;
  std::vector<Eigen::VectorXd> vec_of_Ys;
  
  // changing xlist data
  std::vector<Eigen::VectorXd> vec_of_residuals;
  
protected:
  // fixed data
  // const Eigen::Map<const Eigen::MatrixXd> X;
  const Eigen::Map<const Eigen::VectorXi> nobs_per_study;
  const Eigen::Map<const Eigen::VectorXd> lambda_OLD;
  const double tol;
  const double lambda_min_ratio;
  
  // fixed strongrule data
  Eigen::VectorXd all_set;
};


#endif