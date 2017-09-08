
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utilfunc.h"
using namespace Rcpp;
using namespace arma; 


// sfun <- function(z,lambda){
//   res <- 0
//   if(z>lambda){
//     res <- z-lambda
//   } 
//   if(z < -1*lambda){
//     res <- z+lambda
//   }
//   res
// }

// [[Rcpp::export]]
double sfun_c(double z, double lambda){
  double res = 0;
  if(z > lambda){
    res = z-lambda;
  }
  // z < -lambda
  if(z + lambda < 0){
    res = z+lambda;
  }
  return(res);
}

// scadupdate <- function(z,lambda,gamq=3.7){
//   res <-0
//   if(abs(z)<=2*lambda){
//     res <- sfun(z,lambda)
//   }
//   if((abs(z)>2*lambda) & (abs(z)<= gamq*lambda)){
//     res <- sfun(z,lambda*gamq/(gamq-1))/(1-1/(gamq-1))
//   }
//   if(abs(z)>gamq*lambda){
//     res <- z
//   }
//   res
// }


// [[Rcpp::export]]
double scadupdate_c(double z, double lambda, double gam, double wjj){
  double res = 0;
  if(fabs(z) <= 2.0* lambda){
    res = sfun_c(z,lambda)/wjj;
  }
  if((fabs(z) > 2.0* lambda)&(fabs(z) <= gam* lambda)){
    res = sfun_c(z,lambda*gam / (gam-1.0))/(wjj-1.0/(gam-1.0));
  }
  if(fabs(z) > gam* lambda){
    res = z/wjj;
  }
  return(res);
  
}


// [[Rcpp::export]]
double elnetupdate_c(double z, double lambda, double alpha, double wjj){
  double res = 0;
  res = sfun_c(z,lambda*alpha)/(wjj+lambda*(1-alpha));
  return(res);
  
}


// scadpen <- function(X,Y,lambda){
//   b0 <- solve(t(X)%*% X) %*% (t(X) %*% Y)
//   b1 <- b0
//   n <- dim(X)[1]
//   p <- dim(X)[2]
//   z <- rep(0,p)
//   r <- Y-X %*% b0
//   i <- 0
//   while(i < 1000){
//     for(j in 1:p){
//       z[j] <- solve(t(X[,j])%*% X[,j]) %*% t(X[,j]) %*% r  + b0[j]
//       b1[j] <- scadupdate(z[j],lambda,3.7)
//       r <- r - (b1[j]-b0[j])*X[,j]
//     }
//     i <- i+1
// #print(t(b0))
//     if(abs(sum((Y-X%*%b0)^2)-sum((Y-X%*%b1)^2))<1e-5){
//       res <- b1
//       break
//     } else {
//       b0 <- b1
//     }
//   }
// #cat(i,"\n")
//   res
//     
// }


arma::vec elnetint_c(arma::mat A, arma::vec B, double lambda, arma::vec betao, arma::vec weight,  double alpha, double abstol, int* totiter){
  int p = A.n_cols;
  arma::vec b1= vec(p);
  arma::vec b0= betao;
  double z;
  arma::uvec actset =find(abs(b0)>0);
  if(actset.n_elem==0){actset = zeros<uvec>(1);}
  int iter = 0;
  arma::uvec temp(1);
  arma::vec lw = lambda*weight;
  //arma::vec temp;
  do
  {
    b1 = b0;
    for(int j=0; j < p; j++){
      temp(0) =j;
      z = B(j)-as_scalar(dot(A.submat(actset,temp) , b1(actset)))  +  A(j,j)* b1(j);
      b1(j) = elnetupdate_c(z, lw(j), alpha,A(j,j));
      if(fabs(b1(j))>0) {
        if(!any(actset == j)) {
          actset.insert_rows(actset.n_elem,temp);
          //cout << "Add!" << endl;
        }
      }
    }
    if(norm(b1-b0)<abstol){
      break;
    } else {
      b0 = b1;
      actset =find(abs(b0)>0);
      iter++;
    } 
  }while((iter<10000));
  if(iter==10000) cout << "ELNET fails to converge" << endl;
  *totiter = iter;
  return(b1);
}


arma::vec scadint_c(arma::mat A, arma::vec B, double lambda, arma::vec betao, double gam, double abstol, int* totiter){
  int p = A.n_cols;
  arma::vec b1= betao;
  arma::vec b0= betao;
  arma::vec bp = betao;
  double z;
  arma::uvec actset =find(abs(b0)>0);
  if(actset.n_elem==0){actset = zeros<uvec>(1);}
  int iter = 0;
  //int idx;
  arma::uvec temp(1);
  //arma::vec temp;
  do
  { //cout << "Outer loop" << iter << ": " << endl;
    //b1 = b0;
    //actset.print();
    for(int j=0; j < p; j++){
      //idx = actset(j);
      
      temp(0) =j;
      z = B(j)-as_scalar(dot(A.submat(actset,temp) , b1(actset)))  +  A(j,j)* b1(j);
      //cout << z << endl;
      b1(j) = scadupdate_c(z,lambda,gam,A(j,j));
      if(fabs(b1(j))>0) {
        if(!any(actset == j)) {
          actset.insert_rows(actset.n_elem,temp);
          //cout << "Inside " << j << ": " << endl;
          //b1.t().raw_print(cout,"b1: ");
          //actset.t().raw_print(cout,"actset: ");
          //cout << "Add!" << endl;
        }
      }
      
      
      //actset = unique(actset);
    }
    //for(int j=0; j < p; j++){
    //  b1(j) = (B(j)-as_scalar(dot(A.col(j) , b1))  +  A(j,j)* b0(j))*pow(b0(j),2)/(A(j,j)*pow(b0(j),2)+lambda);
    //}
    
    //b1.t().raw_print(cout,"b1: ");
    //actset.t().raw_print(cout,"actset: ");
    //cout << "Try to rescue" <<endl;
    //cout << iter << ":===================" << endl;
    //b1.t().raw_print(cout,"b1: ");
    //b0.t().raw_print(cout,"b0: ");
    //bp.t().raw_print(cout,"bp: ");
    if(norm(b1-b0)<abstol){
      break;
    } else {
      //if(norm(b1-bp)<1e-8){
      //  b1.fill(1);
      //}
      //if(iter>=2){
      //  bp=b0;
      //}
      b0 = b1;
      actset =find(abs(b1)>0);
      //b1.t().raw_print(cout);
      //actset.t().raw_print(cout,"=======Final actset: ");
      iter++;
    }
  }while((iter<1000));
  if(iter==1000) cout <<  "SCAD fails to converge" << endl;
  *totiter = iter;
  return(b1);
}

// // [[Rcpp::export]]
// arma::mat scadlseq_c(arma::mat X, arma::vec Y, arma::vec lambdav, double gam, double abstol){
//   mat res= mat(X.n_cols,lambdav.n_elem);
//   for (std::size_t i=0; i < lambdav.n_elem; i++){
//     res.col(i) = scadreg_c(X, Y, lambdav(i), gam, abstol);
//   }
//   return(res);
// }

// [[Rcpp::export]]
List penreg_c(arma::mat A, arma::vec B, int n, int p, arma::vec lambdav, arma::vec weight, double alpha, double gamma, std::string method, double abstol){
  int iter = 0;
  int j;
  List out(3);
  arma::vec beta0 = vec(A.n_cols);
  arma::mat res = mat(p,lambdav.n_elem);
  arma::rowvec totaliter = rowvec(lambdav.n_elem);
  if(!method.compare("scad")){
    int ifdiver =0;
    for(j=0;j < lambdav.n_elem; j++){
      iter=0;
      if(j==0){
        beta0.fill(0);
      } else {
        beta0 = res.col(j-1);
        //beta0.randu();
      }
      res.col(j) = scadint_c(A, B, lambdav(j), beta0, gamma, abstol,&iter);
      //if((ifdiver<=1)|(j<lambdav.n_elem/2)){
      //  res.col(j) = scadint_c(A, B, lambdav(j), beta0, gamma, abstol,&iter);
      //} else {
      //  res.col(j).zeros();
      //}
      //if(iter==10000) ifdiver=ifdiver+1;
      totaliter(j) = iter;
      //if(iter>=10000) cout << j << ": " << "SCAD update fails to converge" << endl;
      //if(){
      //  maxrun = j;
      //  break;
      //}
    }
  } else {
    for(j=0;j < lambdav.n_elem; j++){
      if(j==0){
        beta0.fill(0);
      } else {
        beta0 = res.col(j-1);
      }
      res.col(j) = elnetint_c(A, B, lambdav(j), beta0, weight, alpha,abstol,&iter);
      totaliter(j) = iter;
      //if(iter>=10000) cout << j << ": " << "ELNET update fails to converge" << endl;
    }
  }
  
  out[0]=res;
  out[1]=lambdav.t();
  out[2]=totaliter;
  return(out);
}

// [[Rcpp::export]]
List pencv_c(arma::mat A, arma::vec B, List trainL, List testL, int n, int p, arma::vec lambdav, arma::vec weight, int nfold, double alpha, double gamma, std::string method, double abstol){
  int iter = 0;
  int i,j;
  arma::mat beta = mat(p,lambdav.n_elem);
  arma::mat msemat = mat(lambdav.n_elem,nfold);
  //List out(5);
  List reseach(3);
  List orires(3);
  List oneres(3);
  List train(2);
  List test(2);
  if(!method.compare("scad")){
    orires = penreg_c(A, B, n, p, lambdav, weight,alpha,gamma,"scad",abstol);
    for(i=0; i<nfold; i++){
      train = Rcpp::as<Rcpp::List>(trainL[i]);
      test = Rcpp::as<Rcpp::List>(testL[i]);
      oneres = penreg_c(train[0], train[1], n-n/nfold, p, lambdav, weight,alpha,gamma,"scad",abstol);
      beta = Rcpp::as<arma::mat>(oneres[0]);
      //cout << i <<" :" << endl;
      //cout << "=================================" << endl;
      //beta.print();
      
      msemat.col(i) = rss_c(test[0],test[1],beta);
    }
  } else { 
    orires = penreg_c(A, B, n, p, lambdav, weight,alpha,gamma,"lasso",abstol);
    for(i=0; i<nfold; i++){
      train = Rcpp::as<Rcpp::List>(trainL[i]);
      test = Rcpp::as<Rcpp::List>(testL[i]);
      oneres = penreg_c(train[0], train[1], n-n/nfold, p, lambdav, weight,alpha,gamma,"lasso",abstol);
      beta = Rcpp::as<arma::mat>(oneres[0]);
      msemat.col(i) = rss_c(test[0],test[1],beta);
    }
  }
  arma::uvec status = uvec(lambdav.n_elem);
  status = find(all(abs(Rcpp::as<arma::mat>(orires[0]))<1e-10));
  arma::vec amse = mean(msemat,1);
  int idx;
  if(status.n_elem>0) {
    amse(status).fill(1e100);
    idx = amse.index_min();
  } else {
    idx = amse.index_min();
  }
  //arma::vec amse = mean(msemat,1);
  //int idx = amse.index_min();
  arma::vec betaopt = Rcpp::as<arma::mat>(orires[0]).col(idx);
  arma::vec paropt = Rcpp::as<arma::mat>(orires[1]).col(idx);
  
  return List::create(Named("betaopt") = betaopt,
                      Named("paropt") = paropt,
                      Named("barregpath") = orires,
                      Named("amse") = amse,
                      Named("msemat") = msemat,
                      Named("optidx") = idx);
  
}



