
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::mat zbar_c(arma::mat zmat, arma::vec Cv, double s){
  arma::vec indi(Cv.n_elem);
  int sumind = 0;
  for (int i =0; i < int(Cv.n_elem); i++){
    if (Cv[i]>=s){
      indi[i]=1;
      sumind++;
    } else {
      indi[i] = 0;
    }
  }
  //(zmat %*% (Cv>=s))/sum(Cv>=s)
  
  arma::mat zbarmat = (zmat * indi)/sumind;
  
  return(zbarmat) ;
}

// [[Rcpp::export]]
arma::mat zbarsum_c(arma::mat zmat, arma::vec Cv, arma::vec t){
  int l= t.n_elem;
  arma::mat zbarsum = arma::mat(zmat.n_rows,1,fill::zeros);
  for (int i =0; i < l; i++){
    zbarsum = zbarsum + zbar_c(zmat,Cv,t[i]);
  }
  //(zmat %*% (Cv>=s))/sum(Cv>=s)
  
  
  
  return(zbarsum) ;
}

//stzi <- function(z,s){
//  (z-zbar(s)) %*% t(z-zbar(s))
//}


// vinti <- function(zmat,Cv,i){
// 
//   Reduce("+",mapply("*", lapply(Cv[1:i],function(s) stzi(zmat[,i],s)) ,Cv[1:i]-c(0,Cv)[1:i],SIMPLIFY=F))
// }

// [[Rcpp::export]]
arma::mat V_c(arma::mat zmat, arma::vec Cv){
  int l = Cv.n_elem;
  //cout << l << endl;
  arma::mat vintsum = arma::mat(zmat.n_rows,zmat.n_rows,fill::zeros);
  for (int i=0; i<l;i++){
    arma::vec Zi = zmat.col(i);
    for(int j=0; j<i+1; j++){
      double tdiff;
      if (j==0){
        tdiff=Cv[0];
      } else {
        tdiff=Cv[j]-Cv[j-1];
      }
      arma::mat Zicen =Zi-zbar_c(zmat,Cv,Cv[j]);
      vintsum = vintsum + (Zicen * Zicen.t()) * tdiff;
    }
  }
  
  return(vintsum/l);
}


// [[Rcpp::export]]
arma::mat Omega_c(arma::mat zmat, List tlist, arma::vec Cv){
  arma::mat res= mat(zmat.n_rows,zmat.n_rows,fill::zeros);
  int n = tlist.size();
  double weight;
  for(int i=0;i < n; i++){
    for(int j =0; j <n; j++){
      arma::vec tv=tlist[j];
      int mj = tv.n_elem;
      for(int k=0; k < mj; k++){
        double csum =0;
        for(int l =0; l < n; l++){
          csum = csum + (Cv[l]>tv[k]);
        }
        weight = (Cv[i]>tv[k] ? 1 : 0)*1.0/csum;
        arma::mat Zicen = (zmat.col(i)-zbar_c(zmat,Cv,tv[k]))*tv[k];
        res = res + (Zicen * Zicen.t()) * weight;
      }
    }
  }
  return res;
}

// [[Rcpp::export]]
arma::mat P_c(arma::mat zmat, List tlist, List olist, arma::vec Cv){
  arma::mat res= mat(zmat.n_rows,1,fill::zeros);
  int n = tlist.size();
  
  for(int i=0;i < n; i++){
    double weight;
    arma::vec tv=tlist[i];
    arma::vec ov=olist[i];
    int mi = tv.n_elem;
    for(int j =0; j <mi; j++){
      weight = ov[j];
      arma::mat Zicen =(zmat.col(i)-zbar_c(zmat,Cv,tv[j]))*tv[j];
      res = res + Zicen * weight;
    }
  }
  return res;
}

// [[Rcpp::export]]
List reggenint_c(arma::vec Cv, arma::mat Zmat, arma::vec olist, arma::vec tlist, arma::uvec mind){
  List out(2);
  int n=Cv.n_elem;
  int p = Zmat.n_rows;
  
    arma::mat Ctemmat = repmat(Cv,1,tlist.n_elem);
    arma::mat Ttemmat = repmat(tlist.t(),n,1);
    arma::mat Indmat = arma::conv_to< arma::mat >::from(Ctemmat>=Ttemmat);
    arma::mat Indcolsum = sum(Indmat,0);
    
    //arma::vec stdInd = ;
    arma::mat Zbarmat = Zmat * Indmat;
    Zbarmat.each_row() /= Indcolsum;

    arma::mat Zextmat = reshape(repmat(Zmat,tlist.n_elem,1),p,tlist.n_elem*n);
    arma::mat Zbarextmat = Zextmat - repmat(Zbarmat,1,n);
    //Zdiff
    arma::mat Zdiff = Zbarextmat.each_row() % repmat(tlist.t(),1,n);
    arma::mat Zdifftw = Zdiff.each_row() % sqrt(vectorise(Indmat.each_row() / Indcolsum,1));
    arma::mat Omega = Zdifftw * Zdifftw.t();
    //length(tlist)*(rep(1:100,mlist)-1)+1:length(tlist)

    //mind.print();
    arma::vec Pv = (Zdiff.cols(mind)*olist);
    out[0] =Omega;
    out[1] = Pv;
    return(out);
}


// arma::mat diagmult_c(arma::vec diagv, arma::mat X){
//   int n = X.n_cols;
//   arma::mat res = mat(size(X));
//   int i,j;
//   for(i=0;i<n;i++){
//     for(j=0;j<i;j++){
//       res(i,j)=diagv(i) * diagv(j) * X(i,j);
//       res(j,i)=res(i,j);
//     }
//   }
//   for(i=0;i<n;i++){
//     res(i,i)=diagv(i) * diagv(i) * X(i,i);
//   }
//   return res;
// }
// 
// // [[Rcpp::export]]
// arma::vec barrdglupd_c(arma::mat A, arma::vec B, double lambda, arma::vec beta0){
//   arma::mat res = diagmult_c(beta0,inv(diagmult_c(beta0,A) + lambda * eye(size(A))))  * B;
//   return(res);
// }
// 
// // [[Rcpp::export]]
// arma::mat barrdghupd_c(arma::mat X, arma::vec B, double lambda, arma::vec weight){
//   int p = X.n_cols;
//   int r = X.n_rows;
//   arma::mat xga = mat(r,p);
//   xga = X * diagmat(weight);
//   arma::mat M= xga * xga.t();
//   M.diag() = M.diag() + lambda;
//   arma::mat M1= xga.t() * inv_sympd(M) * xga;
//   arma::mat M2 = mat(p,p);
//   double dd;
//   for(int i=0;i<p;i++){
//     for(int j=0;j<=i;j++){
//       if(j==i){
//         dd=1;
//       } else {
//         dd=0;
//       }
//       M2(i,j) = weight(i)*weight(j)*(dd-M1(i,j));
//       M2(j,i) = M2(i,j);
//     }
//   }
//   arma::vec res = M2 * B /lambda;
//   return(res);
// }

// updatefun <- function(X,Y,lambda,beta0)
// {
//   Gammat <- diag(c(beta0))
//   Xmat <- X %*% Gammat
//   Gammat %*% solve(t(Xmat) %*% Xmat + lambda * diag(rep(1,length(beta0)))) %*% t(Xmat) %*% Y
//   
// }

// // [[Rcpp::export]]
// arma::vec updatefun_c(arma::mat X, arma::vec Y, double lambda, arma::vec beta0){
//   arma::mat Gammat = diagmat(beta0);
//   arma::mat Xmat = X * Gammat;
//   arma::mat res = Gammat * inv(Xmat.t() * Xmat + lambda * eye(size(Gammat))) * Xmat.t() * Y;
//   return(res);
// }


// repeat{
//   
// #browser()
// #cat(c(beta0))
// #cat("\n")
//   beta1 <- updatefun(XX,YY,lambda,beta0)
//     beta1 <- ifelse(beta1<1e-10,0,beta1)
// #if(abs(beta1[1,1]-1.281336e-07)<1e-10) browser()
//     
//     if(mean(abs((beta1-beta0))) < abstol) break else beta0 <- beta1
//       iter <- iter +1 
// #print(c(iter,beta1))
// #resr <- rbind(resr,c(iter,beta1))
//     if(iter>1000) {
//       break
//     }
// #print(c(beta0))
// #cprint("\n")
// }


// // [[Rcpp::export]]
// arma::vec barrep_c(arma::mat X, arma::vec Y, double lambda, arma::vec betao, double abstol){
//   arma::vec beta1;
//   arma::vec beta0 = betao;
//   int iter = 0;
//   do
//   {
//     beta1 = updatefun_c(X,Y,lambda,beta0);
//     iter++;
//     //for(std::size_t i=0; i < betao.n_elem; i++){
//     //  if(beta1(i)<1e-10){
//     //    beta1(i) = 0;
//     //  }
//     //}
//     if(max(abs(beta1-beta0))<abstol){
//       break;
//     } else {
//       beta0 = beta1;
//     }
//   } while((iter<1000));
//   cout << "The iteration time is " << iter << endl;
//   if(iter>=1000) cout << lambda << ": " << "Fail to converge" << endl;
//   return(beta1);
// }
// 
// // [[Rcpp::export]]
// arma::mat barlseq_c(arma::mat X, arma::vec Y, arma::vec lambdav, arma::vec betao, double abstol){
//   arma::mat res= arma::mat(X.n_cols,lambdav.n_elem);
//   for (std::size_t i=0; i < lambdav.n_elem; i++){
//     res.col(i) = barrep_c(X, Y, lambdav(i), betao, abstol);
//   }
//   return(res);
// }
