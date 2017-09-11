
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stdbool.h>
#include "utilfunc.h"
using namespace Rcpp;
using namespace arma;
using namespace std;



// All kinds of update functions

// [[Rcpp::export]]
arma::vec barrdglupd_c(arma::mat A, arma::vec B, double lambda, arma::vec beta0){
  arma::mat res = diagmult_c(beta0,inv(diagmult_c(beta0,A) + lambda * eye(size(A))))  * B;
  return(res);
}

// [[Rcpp::export]]
arma::mat barrdghupd_c(arma::mat X, arma::vec B, double lambda, arma::vec weight){
  int p = X.n_cols;
  int r = X.n_rows;
  arma::mat xga = mat(r,p);
  xga = X * diagmat(weight);
  arma::mat M= xga * xga.t();
  M.diag() = M.diag() + lambda;
  arma::mat M1= xga.t() * inv_sympd(M) * xga;
  arma::mat M2 = mat(p,p);
  double dd;
  for(int i=0;i<p;i++){
    for(int j=0;j<=i;j++){
      if(j==i){
        dd=1;
      } else {
        dd=0;
      }
      M2(i,j) = weight(i)*weight(j)*(dd-M1(i,j));
      M2(j,i) = M2(i,j);
    }
  }
  arma::vec res = M2 * B /lambda;
  return(res);
}


arma::vec barrdgH_c(arma::mat X, arma::mat A, arma::vec B, double lambda, arma::vec betao, double abstol, int* totiter){
  arma::vec beta1 = vec(betao.n_elem);
  arma::vec beta0 = betao;
  int p = X.n_cols;
  unsigned int r = X.n_rows;
  arma::uvec actset =conv_to< uvec >::from(linspace(0, p-1,p));
  int iter = 0;
  do{
    beta1.fill(0);
    //if(actset.n_elem <= r){
    if(actset.n_elem >= r){
      beta1.rows(actset) = barrdghupd_c(X.cols(actset), B.rows(actset), lambda, beta0.rows(actset));
    } else {
      beta1.rows(actset) = barrdglupd_c(A.submat(actset,actset),B.rows(actset),lambda,beta0.rows(actset));
    }
    actset = find(abs(beta1)>1e-20);
    //std::cout.setf(std::ios::fixed);
    //std::cout.precision(10);
    //beta1.t().raw_print(cout);
    iter++;
    beta1( find(abs(beta1) < 1e-20) ).zeros();
    //for(int i=0; i < betao.n_elem; i++){
    //  if(fabs(beta1(i))<1e-20){
    //    beta1(i) = 0;
    //  }
    //}
    if(norm(beta1-beta0)<abstol){
        break;
      } else {
        beta0 = beta1;
        
    }
  } while((iter<3000));
  
  //cout << "The iteration time is " << iter << endl;
  if(iter>=3000) cout << iter << ": " << "Ridge update fails to converge" << endl;
  *totiter = iter;
  return(beta1);
}



arma::vec barrdgL_c(arma::mat A, arma::vec B, double lambda, arma::vec betao, double abstol, int* totiter){
  int p = A.n_cols;
  //int n = A.n_rows;
  arma::vec beta1 = vec(betao.n_elem);
  arma::vec beta0 = betao;
  arma::uvec actset =conv_to< uvec >::from(linspace(0, p-1,p));
  int iter = 0;
  do{
      beta1.fill(0);
      beta1.rows(actset) = barrdglupd_c(A.submat(actset,actset),B.rows(actset),lambda,beta0.rows(actset));
      actset = find(abs(beta1)>1e-20);
      //std::cout.setf(std::ios::fixed);
      //std::cout.precision(10);
      //beta1.t().raw_print(cout);
      iter++;
      beta1( find(abs(beta1) < 1e-20) ).zeros();
      //for(int i=0; i < betao.n_elem; i++){
      //  if(fabs(beta1(i))<1e-20){
      //    beta1(i) = 0;
      //  }
      //}
      if(norm(beta1-beta0)<abstol){
        break;
      } else {
        beta0 = beta1;
      }
    } while((iter<3000));
  
  
  //cout << "The iteration time is " << iter << endl;
  if(iter>=3000) cout << iter << ": " << "Ridge update fails to converge" << endl;
  *totiter = iter;
  return(beta1);
}


// Add different lambdav for differen xi
// [[Rcpp::export]]
List barridge_c(arma::mat A, arma::vec B, int n, int p, arma::vec xiv, arma::vec lambdav, double abstol,int stand){
  int iter = 0;
  unsigned int i,j;
  List out(3);
  arma::vec beta0;
  arma::mat res = mat(p,xiv.n_elem*lambdav.n_elem);
  arma::rowvec totaliter = rowvec(xiv.n_elem*lambdav.n_elem);
  int current_index=0;
  arma::mat xil = mat(2,xiv.n_elem*lambdav.n_elem);
  arma::vec std;
  // Standardize
  if(stand==1){
    List temp(3);
    temp = standard_c(A,B,1);
    A = Rcpp::as<arma::mat>(temp[0]);
    B = Rcpp::as<arma::mat>(temp[1]);
    std = Rcpp::as<arma::vec>(temp[2]);
  }
  
  
  int r=p;
  if(n<p){
    r = rank(A);
  }
  if(r >= p){
    for(i=0; i < xiv.n_elem; i++){
      beta0 = inv_sympd(A+xiv(i)*eye(size(A))) * B;
      for(j=0;j < lambdav.n_elem; j++){
        current_index =i*lambdav.n_elem+j;
        res.col(current_index) = barrdgL_c(A,B,lambdav(j),beta0,abstol,&iter);
        totaliter(current_index) = iter;
        xil(0,current_index) = xiv(i);
        xil(1,current_index) = lambdav(j);
        if(iter>=3000) cout << i*lambdav.n_elem+j << ": " << "Ridge update fails to converge" << endl;
      }
    }
  } else {
    arma::mat V;
    arma::vec d;
    eigs_sym(d,V,sp_mat(A),r);
    arma::mat Xt = diagmat(sqrt(d)) * V.t();
    for(i=0; i < xiv.n_elem; i++){
      beta0 = inv_sympd(A+xiv(i)*eye(size(A))) * B;
      for(j=0;j < lambdav.n_elem; j++){
        current_index =i*lambdav.n_elem+j;
        res.col(current_index) = barrdgH_c(Xt, A, B, lambdav(j), beta0, abstol, &iter);
        totaliter(current_index) = iter;
        xil(0,current_index) = xiv(i);
        xil(1,current_index) = lambdav(j);
        if(iter>=3000) cout << i*lambdav.n_elem+j << ": " << "Matrix update fails to converge" << endl;
      }
    }
    
  }
  
  if(stand==1){
    res = diagmat(std) * res; 
  }
  out[0]=res;
  out[1]=xil;
  out[2]=totaliter;
  return(out);
}

// iterative Coordinate descent
arma::vec barcdint_c(arma::mat A, arma::vec B, double lambda, arma::vec beta0, double abstol, int* totiter){
  int p = A.n_cols;
  arma::vec b1= vec(p);
  arma::vec b0=beta0;
  arma::uvec actset =conv_to< uvec >::from(linspace(0, p-1,p));
  //double maxchange = 0.0;
  //arma::vec z = vec(p);
  //arma::vec r = Y - X * b0;
  int iter = 0;
  int idx;
  arma::uvec temp(1);
  //arma::vec temp;
  do
  {
    //cout << iter << ": " <<endl;
    //maxchange =0;
    b1 = b0;
    for(unsigned int j=0; j < actset.n_elem; j++){
      idx = actset(j);
      temp(0) =idx;
      b1(idx) = (B(idx)-as_scalar(dot(A.submat(actset,temp) , b1(actset)))  +  A(idx,idx)* b0(idx))*pow(b0(idx),2)/(A(idx,idx)*pow(b0(idx),2)+lambda);
    }
    actset = find(abs(b1)>1e-20);
    if(norm(b1-b0)<abstol){
      break;
    } else {
      b0.fill(0);
      b0(actset) = b1(actset);
      iter++;
    }
    
  } while((iter<3000));
  if(iter==3000) cout << "Iterative Coordinate fails to converge" << endl;
  *totiter = iter;
  return(b1);
}


// add time count, add active set
// [[Rcpp::export]]
List barcd_c(arma::mat A, arma::vec B, int n, int p, arma::vec xiv, arma::vec lambdav, double abstol){
  int iter = 0;
  unsigned int i,j;
  List out(3);
  arma::vec beta0;
  arma::mat res = mat(p,xiv.n_elem*lambdav.n_elem);
  arma::rowvec totaliter = rowvec(xiv.n_elem*lambdav.n_elem);
  int current_index=0;
  arma::mat xil = mat(2,xiv.n_elem*lambdav.n_elem);
  for(i=0; i < xiv.n_elem; i++){
      beta0 = inv_sympd(A+xiv(i)*eye(size(A))) * B;
      for(j=0;j < lambdav.n_elem; j++){
        current_index =i*lambdav.n_elem+j;
        res.col(current_index) = barcdint_c(A, B, lambdav(j), beta0, abstol, &iter);
        totaliter(current_index) = iter;
        xil(0,current_index) = xiv(i);
        xil(1,current_index) = lambdav(j);
        if(iter>=3000) cout << i*lambdav.n_elem+j << ": " << "CD update fails to converge" << endl;
      }
    }
  
  out[0]=res;
  out[1]=xil;
  out[2]=totaliter;
  return(out);
}




// [[Rcpp::export]]
List barcv_c(arma::mat A, arma::vec B, List trainL, List testL, int n, int p, arma::vec xiv, arma::vec lambdav, int nfold, std::string method, double abstol, int stand){
  //int iter = 0;
  int i;
  arma::mat beta = mat(p,xiv.n_elem*lambdav.n_elem);
  arma::uvec status = uvec(xiv.n_elem*lambdav.n_elem);
  arma::mat msemat = mat(xiv.n_elem*lambdav.n_elem,nfold);
  arma::mat betamat = mat(p,xiv.n_elem*lambdav.n_elem*nfold);
  //List out(5);
  List reseach(3);
  List orires(3);
  List oneres(3);
  List train(2);
  List test(2);
  if(!method.compare("cd")){
    orires = barcd_c(A, B, n, p, xiv, lambdav, abstol);
    for(i=0; i<nfold; i++){
      train = Rcpp::as<Rcpp::List>(trainL[i]);
      test = Rcpp::as<Rcpp::List>(testL[i]);
      oneres = barcd_c(train[0], train[1], n-n/nfold, p, xiv, lambdav, abstol);
      beta = Rcpp::as<arma::mat>(oneres[0]);
      betamat.cols(xiv.n_elem*lambdav.n_elem*i,xiv.n_elem*lambdav.n_elem*(i+1)-1) = beta;
      msemat.col(i) = rss_c(test[0],test[1],beta);
    }
  } else { 
    orires = barridge_c(A, B, n, p, xiv, lambdav, abstol,stand);
    for(i=0; i<nfold; i++){
      train = Rcpp::as<Rcpp::List>(trainL[i]);
      test = Rcpp::as<Rcpp::List>(testL[i]);
      oneres = barridge_c(train[0], train[1], n-n/nfold, p, xiv, lambdav, abstol,stand);
      beta = Rcpp::as<arma::mat>(oneres[0]);
      betamat.cols(xiv.n_elem*lambdav.n_elem*i,xiv.n_elem*lambdav.n_elem*(i+1)-1) = beta;
      msemat.col(i) = rss_c(test[0],test[1],beta);
    }
  }
  status = find(all(abs(Rcpp::as<arma::mat>(orires[0]))<1e-10));
  //status.t().print();
  
  
  arma::vec amse = mean(msemat,1);
  int idx;
  if(status.n_elem>0) {
    amse(status).fill(1e100);
    idx = amse.index_min();
  } else {
    idx = amse.index_min();
  }
   
  arma::vec betaopt = Rcpp::as<arma::mat>(orires[0]).col(idx);
  arma::vec paropt = Rcpp::as<arma::mat>(orires[1]).col(idx);
  return List::create(Named("betaopt") = betaopt,
                      Named("paropt") = paropt,
                      Named("barregpath") = orires,
                      Named("amse") = amse,
                      Named("msemat") = msemat,
                      Named("optidx") = idx,
                      Named("betamat") = betamat);

}

// arma::mat barlseq_c(arma::mat X, arma::vec Y, arma::vec lambdav, arma::vec betao, double abstol, int if_hd){
//   arma::mat res= arma::mat(X.n_cols,lambdav.n_elem);
//   for (std::size_t i=0; i < lambdav.n_elem; i++){
//     res.col(i) = barrep_c(X, Y, lambdav(i), betao, abstol);
//   }
//   return(res);
// }



