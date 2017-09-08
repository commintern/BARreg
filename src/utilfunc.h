// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

arma::mat diagmult_c(arma::vec diagv, arma::mat X);
arma::vec rss_c(arma::mat A, arma::vec B, arma::mat beta);
arma::vec rssf_c(arma::mat A, arma::vec B, arma::mat beta, double tss);
List standard_c(arma::mat A, arma::vec B, int n);