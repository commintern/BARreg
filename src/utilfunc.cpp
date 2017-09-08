// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

arma::mat diagmult_c(arma::vec diagv, arma::mat X){
  int n = X.n_cols;
  arma::mat res = mat(size(X));
  int i,j;
  for(i=0;i<n;i++){
    for(j=0;j<i;j++){
      res(i,j)=diagv(i) * diagv(j) * X(i,j);
      res(j,i)=res(i,j);
    }
  }
  for(i=0;i<n;i++){
    res(i,i)=diagv(i) * diagv(i) * X(i,i);
  }
  return res;
}

// [[Rcpp::export]]
arma::vec rss_c(arma::mat A, arma::vec B, arma::mat beta){
  int p =beta.n_elem;
  //double p1 = as_scalar(B.t() * inv_sympd(A) * B);
  double p1=0;
  arma::vec p2 = -2 * beta.t() *  B;
  arma::vec p3 =vec(beta.n_cols);
  for(int i = 0; i<beta.n_cols; i++){
    p3(i) = as_scalar(beta.col(i).t() * A * beta.col(i));
  }
  arma::vec res = p1+p2+p3;
  return(res);
}

// [[Rcpp::export]]
arma::vec rssf_c(arma::mat A, arma::vec B, arma::mat beta, double tss){
  int p =beta.n_elem;
  arma::vec p2 = 2 * beta.t() *  B;
  arma::vec p3 =vec(beta.n_cols);
  for(int i = 0; i<beta.n_cols; i++){
    p3(i) = as_scalar(beta.col(i).t() * A * beta.col(i));
  }
  arma::vec res = tss+p2+p3;
  return(res);
}

// [[Rcpp::export]]
List standard_c(arma::mat A, arma::vec B, int n){
  List out(3);
  arma::vec stdft = n/sqrt(A.diag());
  arma::mat W = diagmat(stdft);
  arma::mat As = W * A * W;
  arma::mat Bs = W * B;
  out[0] = As;
  out[1] = Bs;
  out[2] = stdft;
  return(out);
}