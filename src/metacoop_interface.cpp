#include "metacoop.h"

using Rcpp::as;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::Map;

// [[Rcpp::export]]
Rcpp::List fit_metacoop_cpp(Rcpp::List Xmat_list,
                            Rcpp::List Y_list,
                            int nvars,
                            int nlambda,
                            const Rcpp::IntegerVector &NOBS_PER_STUDY_,
                            Rcpp::List weights_list,
                            const Rcpp::NumericVector &LAMBDA_) // we will take in lambda seq but not do anything with it
{
  // Rcpp::Rcout << "doweenter" << std::endl;
  // type conversion
  // const Map<MatrixXd> X(as<Map<MatrixXd> >(XBLOCK_));
  const Map<VectorXi> nobs_per_study(as<Map<VectorXi> >(NOBS_PER_STUDY_));
  // const Map<VectorXd> weights(as<Map<VectorXd> >(WEIGHTS_));
  Map<VectorXd> lambda(as<Map<VectorXd> >(LAMBDA_));
  
  // helpful const integers
  int nstudies = Xmat_list.size();
  int totalnvars = nstudies*nvars;
  
  // hardcoded for now, add functionality in R to change these
  // int nlambda = 100;
  double lambda_min_ratio = 1e-4;
  
  // Rcpp::Rcout << "nstudies,nvars,totalnvars,nlambda" << std::endl;
  // Rcpp::Rcout << nstudies << std::endl;
  // Rcpp::Rcout << nvars << std::endl;
  // Rcpp::Rcout << totalnvars << std::endl;
  // Rcpp::Rcout << nlambda << std::endl;
  // Rcpp::Rcout << lambda_min_ratio << std::endl;
  
  // rewrite so this is a default not hardcoded
  double tol = 1e-5;
  
  // sort lambda in descending order (maybe just do this in R)
  // std::vector<double> vec(lambda.data(), lambda.data()+lambda.size());
  // std::sort(vec.begin(), vec.end());
  // Map<VectorXd> lambda_sorted(*vec, vec.size());
  
  // vector<MatrixXd> matrices_vector;
  
  // Rcpp::Rcout << "runtime_debug_a" << std::endl;
  
  // getting rid of weights
  metacoop metacoop_obj(Xmat_list, Y_list, nobs_per_study, weights_list, lambda, // fixed data
                        totalnvars, nlambda, nvars, nstudies, // const integers
                          tol, lambda_min_ratio // const double
                          );
  
  // Rcpp::Rcout << "runtime_debug_b" << std::endl;
  
  metacoop_obj.fit_path_xlist();
  
  // Rcpp::Rcout << "runtime_debug_c" << std::endl;
  
  // return Rcpp::List::create( Rcpp::Named("xblockfeed")=X);
    
  return Rcpp::List::create( Rcpp::Named("return_beta")=metacoop_obj.return_beta);
  
  // NOTE:
  // THERE IS AN ISSUE IN C++ WHERE IF WE WANT A ZERO ARGUMENT CONSTRUCTOR
  // WE USE metacoop mymcobj;
  // NOT metacoop mymcobj();
  // metacoop mymcobj;
}