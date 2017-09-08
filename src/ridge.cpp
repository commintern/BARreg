// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma; 

// [[Rcpp::export]]
arma::vec ridreg_c(arma::mat A, arma::vec B, arma::vec binital, double lambda, arma::vec weight, double abstol){
  int p = B.n_rows;
  //int n = Y.n_elem;
  //arma::vec b1= vec(p);
  arma::vec b0= binital;
  //double sse0, sse1;
  //arma::mat A = X.t() * X;
  //arma::mat B = X.t() * Y;
  //b0 = inv(X.t() * X) * X.t() * Y;
  //b0.fill(0);
  double b1; 
  //sse0 = pow(norm(Y- X*b0,2),2);
  
  double maxchange = 0.0;
  arma::vec z = vec(p);
  //arma::vec r = Y - X * b0;
  int iter = 0;
  do
  {
    maxchange =0;
    for(int j=0; j < p; j++){
      z(j) = (B(j)-as_scalar(dot(A.row(j), b0)))/ 100 + b0(j);
      //cout << z(j) << " j is " <<j << endl;
      //z(j) = as_scalar(B(j)-A.row(j) * b1) / n + b0(j);
      b1 = z(j)/(1.0+lambda/weight(j));
      if(fabs(b1-b0(j))>maxchange){ 
        maxchange = fabs(b1-b0(j));
      }
      b0(j) = b1;
    }
    //sse1 = pow(norm(Y- X*b1,2),2);
    if(maxchange < abstol){
      break;
    } else {
      iter++;
    }
    
  } while((iter<1000));
  if(iter==1000) cout << "Fail to converge" << endl;
  cout << "The iteration time is " << iter << endl;
  return(b0);
}



// Coordinate descent
// [[Rcpp::export]]
arma::vec ridreg0_c(arma::mat X, arma::vec Y, arma::vec binital, double lambda, arma::vec weight, double abstol){
  int p = X.n_cols;
  int n = X.n_rows;
  //cout << p << n << endl;
  arma::vec b1= vec(p);
  arma::vec b0= binital;
  //double sse0, sse1;
  //arma::mat A = X.t() * X;
  //arma::mat B = X.t() * Y;
  //b0 = inv(X.t() * X) * X.t() * Y;
  //b0.fill(0);
  
  //double b1; 
  //sse0 = pow(norm(Y- X*b0,2),2);
  
  double maxchange = 0.0;
  arma::vec z = vec(p);
  arma::vec r = Y - X * b0;
  int iter = 0;
  do
  {
    maxchange =0;
    for(int j=0; j < p; j++){
      //z(j) = (B(j)-as_scalar(dot(A.row(j), b0)))/ 100 + b0(j);
      z(j) = as_scalar(dot(X.col(j) , r))  +  n* b0(j); 
      //cout << z(j) << " j is " <<j << endl;
      //z(j) = as_scalar(B(j)-A.row(j) * b1) / n + b0(j);
      b1(j) = z(j)/(n+lambda/pow(weight(j),2));
      r = r - (b1(j)-b0(j)) * X.col(j);
      if(fabs(b1(j)-b0(j))>maxchange){ 
        maxchange = fabs(b1(j)-b0(j));
      }
      //b0(j) = b1;
    }
    //sse1 = pow(norm(Y- X*b1,2),2);
    if(maxchange < abstol){
     
      //cout << b0(1) << endl;
      break;
    } else {
      b0 = b1;
      //cout << b0(1) << endl;
      iter++;
    }
    
  } while((iter<100000));
  if(iter==100000) cout << "Fail to converge" << endl;
  cout << "The iteration time is " << iter << endl;
  return(b1);
}


// iterative Coordinate descent
// [[Rcpp::export]]
arma::vec ridregi_c(arma::mat X, arma::vec Y, double lambda, arma::vec weight, double abstol){
  int p = X.n_cols;
  int n = X.n_rows;
  arma::vec b1= vec(p);
  arma::vec b0= weight;
  
  
  double maxchange = 0.0;
  arma::vec z = vec(p);
  arma::vec r = Y - X * b0;
  int iter = 0;
  do
  {
    maxchange =0;
    for(int j=0; j < p; j++){
      z(j) = as_scalar(dot(X.col(j) , r))  +  n* b0(j); 
      //b1(j) = z(j)*pow(b0(j),2)/(n*pow(b0(j),2)+lambda);
      b1(j) = z(j)/(n+lambda/pow(b0(j),2));
      r = r - (b1(j)-b0(j)) * X.col(j);
      if(fabs(b1(j)-b0(j))>maxchange){ 
        maxchange = fabs(b1(j)-b0(j));
      }
    }
    if(maxchange < abstol){
      break;
    } else {
      b0 = b1;
      iter++;
    }
    
  } while((iter<1000));
  if(iter==1000) cout << "Fail to converge" << endl;
  cout << "The iteration time is " << iter << endl;
  return(b1);
}

// [[Rcpp::export]]
arma::vec ridreg1_c(arma::mat X, arma::vec Y, double lambda, arma::vec weight, double abstol){

 
 int n = Y.n_elem;
  arma::vec res = inv(X.t() * X +  lambda * diagmat(weight) * n) * X.t() * Y;
  
  return(res);
}

// [[Rcpp::export]]
arma::vec ridreg2_c(arma::mat A, arma::vec B, double lambda, arma::vec weight, double abstol){
  
  arma::mat w = diagmat(weight);
  int n = B.n_elem;
  arma::mat xg = w * A * w;
  arma::mat V;
  arma::vec d;
  eig_sym(d,V,xg);
  arma::mat gv = w * V;
  arma::mat cen = diagmat(1/(pow(d,2)+lambda));
  arma::vec res = gv * cen * gv.t() * B;
  
  return(res);
}

// [[Rcpp::export]]
arma::vec ridreg3_c(arma::mat A, arma::vec B, double lambda, arma::vec weight, double abstol){
  
  
  int p = A.n_cols;
  arma::mat gat = mat(p,p);
  
  for(int i=0;i<p;i++){
    for(int j=0;j<=i;j++){
      
      gat(i,j)=weight(i) * weight(j);
      gat(j,i)=gat(i,j);
    }
  }
  //arma::mat gat = weight * weight.t();
  arma::mat M = gat % A;
  M.diag() = M.diag() + lambda;
  arma::vec res = gat % inv_sympd(M) * B;
  
  return(res);
}
// [[Rcpp::export]]
arma::vec ridreg4_c(arma::mat A, arma::vec B, double lambda, arma::vec weight, double abstol){
  
  //arma::mat ga = diagmat(pow(weight,-2));
  arma::vec res = solve(A +  lambda * diagmat(pow(weight,-2)) , B);
  
  return(res);
}

// [[Rcpp::export]]
arma::mat ridregM1_c(arma::mat A, arma::vec B, double lambda, arma::vec weight, double abstol){
  
  
  int p = A.n_cols;
  arma::mat gat = mat(p,p);
  
  for(int i=0;i<p;i++){
    for(int j=0;j<=i;j++){
      
      gat(i,j)=weight(i) * weight(j) * A(i,j);
      gat(j,i)=gat(i,j);
    }
  }
  //arma::mat gat = weight * weight.t();
  
  gat.diag() = gat.diag() + lambda;
  arma::mat M = inv_sympd(gat); 
  //cout << M(1,2) <<endl;
  for(int i=0;i<p;i++){
    for(int j=0;j<=i;j++){
      
      M(i,j)=weight(i) * weight(j) * M(i,j);
      M(j,i)=M(i,j);
    }
  }
  //cout << M.row(1) <<endl;
  arma::vec res = M * B;
  
  return(res);
}

// [[Rcpp::export]]
arma::vec ridregM2_c(arma::mat A, arma::vec B, double lambda, arma::vec weight, int r,double abstol){
  
  
  int p = A.n_cols;
  arma::mat gat = mat(p,p);
  
  for(int i=0;i<p;i++){
    for(int j=0;j<=i;j++){
      
      gat(i,j)=A(i,j)*weight(i) * weight(j);
      gat(j,i)=gat(i,j);
    }
  }
  arma::vec delta;
  arma::mat V;
  arma::eigs_sym( delta, V, sp_mat(gat), r);
  arma::mat gaV = diagmat(weight) * V;
  //arma::mat gat = weight * weight.t();
  arma::mat res = V * diagmat(1/(delta+lambda)) * V.t() * B;
  
  return(res);
}

// [[Rcpp::export]]
arma::vec ridregM3_c(arma::mat X, arma::vec B, double lambda, arma::vec weight, int r,double abstol){
  
  
  int p = X.n_cols;
  arma::mat xga = mat(r,p);
  //for(int i=0;i<r;i++){
  //  for(int j=0;j<p;j++){
  //    xga(i,j) = X(i,j) * weight(j);
  //  }
  //}
  xga = X * diagmat(weight);
  arma::mat U,V;
  arma::vec delta;
  svds(U,delta,V,sp_mat(xga),r);
  arma::mat vga(p,r);
  //for(int i=0;i<p;i++){
  //  for(int j=0;j<r;j++){
  //    vga(i,j) = V(i,j) * weight(i);
  //  }
  //}
  
  vga =diagmat(weight) * V  ;
  cout << size(V) << endl;
  arma::mat M = vga * diagmat(1/(pow(delta,2)+lambda)) * vga.t();
  //cout << M.row(1) << endl;
  arma::mat res = M * B;
  
  return(res);
}


// [[Rcpp::export]]
arma::mat ridregM4_c(arma::mat X, arma::vec B, double lambda, arma::vec weight, int r,double abstol){
  
  
  int p = X.n_cols;
  
  arma::mat xga = mat(r,p);
  //for(int i=0;i<r;i++){
  //  for(int j=0;j<p;j++){
  //    xga(i,j) = X(i,j) * weight(j);
  //  }
  //}
  xga = X * diagmat(weight);
  arma::mat M= xga * xga.t();
  //cout << xga(1,1) << endl;
  M.diag() = M.diag() + lambda;
  arma::mat M1= xga.t() * inv_sympd(M) * xga;
  arma::mat M2 = mat(p,p);
  //printf("%d\n",1);
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
  //arma::mat m15 = (eye(p,p)/lambda-M1/lambda/lambda) ;
  //M2 = diagmat(weight) * m15 * diagmat(weight);
  //cout << M2.row(1) << endl;
  arma::vec res = M2 * B /lambda;
  
  return(res);
}

// [[Rcpp::export]]
arma::mat test1(arma::mat A,arma::vec t){
  return(diagmat(t) * A * diagmat(t));
}

// [[Rcpp::export]]
arma::mat test2(arma::mat A,arma::vec t){
  int p = A.n_cols;
  arma::mat res = mat(p,p);
  int i,j;
  for(i=0;i<p;i++){
    for(j=0;j<=i;j++){
      res(i,j) = A(i,j) * t(j) * t(i);
      res(i,j) = res(j,i);
    }
  }
  return(res);
}
