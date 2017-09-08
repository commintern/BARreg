// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;
// // [[Rcpp::export]]
// List cv_test1_c(arma::mat X, arma::vec Y, arma::vec lambdav, arma::vec betao, double abstol, int nfold){
//   List out(3);
//   int p = X.n_cols;
//   int n = X.n_rows;
//   int nq = n / nfold;
//   
//   arma::vec bi = vec(size(betao));
//   int l,i;
//   arma::mat res = mat(nfold,lambdav.n_elem);
//   arma::mat trainx;
//   arma::mat testx;
//   arma::vec trainy;
//   arma::vec testy;
//   
//   for(i=0; i < nfold-1; i++){
//     trainx = X;
//     trainx.shed_rows(i*nq,nq*(i+1)-1);
//     testx = X.rows(i*nq,nq*(i+1)-1);
//     trainy = Y;
//     trainy.shed_rows(i*nq,nq*(i+1)-1);
//     testy = Y.rows(i*nq,nq*(i+1)-1);
//     for(l = 0; l<lambdav.n_elem;l++){
//       bi = barrep_c(trainx,trainy,lambdav(l),betao,abstol);
//       res(i,l) = pow(norm(testy - testx * bi,2),2);
//     }
//   }
//   trainx = X;
//   trainx.shed_rows((nfold-1)*nq,n-1);
//   testx = X.rows((nfold-1)*nq,n-1);
//   trainy = Y;
//   trainy.shed_rows((nfold-1)*nq,n-1);
//   testy = Y.rows((nfold-1)*nq,n-1);
//   for(l = 0; l<lambdav.n_elem;l++){
//     bi = barrep_c(trainx,trainy,lambdav(l),betao,abstol);
//     res(nfold-1,l) = pow(norm(testy - testx * bi,2),2);
//   }
//   
//   arma::vec errv = trans(mean(res,0));
//   double lopt = as_scalar(lambdav(errv.index_min())); 
//   
//   out[0] = barrep_c(X,Y,lopt,betao,abstol);
//   out[1] = lopt;
//   out[2] = errv;
//   
//   return(out);
// }
// 
// // [[Rcpp::export]]
// List cv_test2_c(arma::mat X, arma::vec Y, arma::vec lambdav, arma::vec betao, double abstol, int nfold){
//   List out(4);
//   int p = X.n_cols;
//   int n = X.n_rows;
//   int nq = n / nfold;
//   
//   arma::vec bi = vec(size(betao));
//   int l,i;
//   arma::mat res = mat(nfold,lambdav.n_elem);
//   arma::mat trainx;
//   arma::mat testx;
//   arma::vec trainy;
//   arma::vec testy;
//   
//   for(i=0; i < nfold-1; i++){
//     trainx = X;
//     trainx.shed_rows(i*nq,nq*(i+1)-1);
//     testx = X.rows(i*nq,nq*(i+1)-1);
//     trainy = Y;
//     trainy.shed_rows(i*nq,nq*(i+1)-1);
//     testy = Y.rows(i*nq,nq*(i+1)-1);
//     
//     for(l = 0; l<lambdav.n_elem;l++){
//       bi = itercd_c(trainx,trainy,p,n,lambdav(l),betao,abstol);
//       res(i,l) = pow(norm(testy - testx * bi,2),2);
//     }
//   }
//   trainx = X;
//   trainx.shed_rows((nfold-1)*nq,n-1);
//   testx = X.rows((nfold-1)*nq,n-1);
//   trainy = Y;
//   trainy.shed_rows((nfold-1)*nq,n-1);
//   testy = Y.rows((nfold-1)*nq,n-1);
//   for(l = 0; l<lambdav.n_elem;l++){
//     bi = itercd_c(trainx,trainy,p,n,lambdav(l),betao,abstol);
//     res(nfold-1,l) = pow(norm(testy - testx * bi,2),2);
//   }
//   
//   arma::vec errv = trans(mean(res,0));
//   double lopt = as_scalar(lambdav(errv.index_min())); 
//   
//   out[0] = itercd_c(X,Y,p,n,lopt,betao,abstol);
//   out[1] = lopt;
//   out[2] = errv;
//   
//   return(out);
// }
// 
// // [[Rcpp::export]]
// List cv_test3_c(arma::mat X, arma::vec Y, arma::vec lambdav, arma::vec xiv, arma::vec betao, double abstol, int nfold){
//   List out(6);
//   int p = X.n_cols;
//   int n = X.n_rows;
//   int nq = n / nfold;
//   
//   arma::vec bi = vec(size(betao));
//   arma::vec bin0 = vec(size(betao));
//   int l,i,j;
//   arma::mat res = mat(nfold,lambdav.n_elem*xiv.n_elem);
//   arma::mat res1 = mat(nfold,lambdav.n_elem*xiv.n_elem);
//   arma::mat bmat = mat(p,lambdav.n_elem);
//   arma::mat bmat1 = mat(p,lambdav.n_elem);
//   arma::mat trainx;
//   arma::mat testx;
//   arma::vec trainy;
//   arma::vec testy;
//   for(i=0; i < nfold-1; i++){
//     trainx = X;
//     trainx.shed_rows(i*nq,nq*(i+1)-1);
//     testx = X.rows(i*nq,nq*(i+1)-1);
//     trainy = Y;
//     trainy.shed_rows(i*nq,nq*(i+1)-1);
//     testy = Y.rows(i*nq,nq*(i+1)-1);
//     for( j =0; j < xiv.n_elem; j++){
//       bin0 = inv_sympd(X.t() * X + xiv(j)*eye(p,p)) * X.t() * Y;
//       for(l = 0; l<lambdav.n_elem;l++){
//         
//         bi = itercd_c(trainx,trainy,p,n,lambdav(l),betao,abstol);
//         res(i,l+j*lambdav.n_elem) = pow(norm(testy - testx * bi,2),2);
//         bi = barrep_c(trainx,trainy,lambdav(l),betao,abstol);
//         res1(i,l+j*lambdav.n_elem) = pow(norm(testy - testx * bi,2),2);
//       }
//     }
//   }
//   trainx = X;
//   trainx.shed_rows((nfold-1)*nq,n-1);
//   testx = X.rows((nfold-1)*nq,n-1);
//   trainy = Y;
//   trainy.shed_rows((nfold-1)*nq,n-1);
//   testy = Y.rows((nfold-1)*nq,n-1);
//   for( j =0; j < xiv.n_elem; j++){
//     bin0 = inv_sympd(X.t() * X + xiv(j)*eye(p,p)) * X.t() * Y;
//     for(l = 0; l<lambdav.n_elem;l++){
//       bi = itercd_c(trainx,trainy,p,n,lambdav(l),betao,abstol);
//       res(nfold-1,l+j*lambdav.n_elem) = pow(norm(testy - testx * bi,2),2);
//       bi = barrep_c(trainx,trainy,lambdav(l),betao,abstol);
//       res1(nfold-1,l+j*lambdav.n_elem) = pow(norm(testy - testx * bi,2),2);
//     }
//   }
//   
//   arma::vec errv = trans(mean(res1,0));
//   double lopt = as_scalar(lambdav(errv.index_min() % lambdav.n_elem)); 
//   double xiopt = as_scalar(xiv(errv.index_min() / lambdav.n_elem));
//   out[0] = barrep_c(X,Y,lopt,inv_sympd(X.t() * X + xiopt*eye(p,p)) * X.t() * Y,abstol);
//   out[1] = lopt; 
//   out[2] = xiopt;
//   
//   errv = trans(mean(res,0));
//   lopt = as_scalar(lambdav(errv.index_min() % lambdav.n_elem)); 
//   xiopt = as_scalar(xiv(errv.index_min() / lambdav.n_elem));
//   out[3] = itercd_c(X,Y,p,n,lopt,inv_sympd(X.t() * X + xiopt*eye(p,p)) * X.t() * Y,abstol);
//   out[4] = lopt;
//   out[5] = xiopt;
//   
//   
//   
//   return(out);
// }

// [[Rcpp::export]]
arma::vec test(arma::vec a, double b){
  arma::vec temp = vec(1);
  temp(0) = b;
  a.insert_rows(a.n_elem,temp);
  return a;
}
