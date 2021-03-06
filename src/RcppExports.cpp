// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// barrdglupd_c
arma::vec barrdglupd_c(arma::mat A, arma::vec B, double lambda, arma::vec beta0);
RcppExport SEXP bar_barrdglupd_c(SEXP ASEXP, SEXP BSEXP, SEXP lambdaSEXP, SEXP beta0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta0(beta0SEXP);
    rcpp_result_gen = Rcpp::wrap(barrdglupd_c(A, B, lambda, beta0));
    return rcpp_result_gen;
END_RCPP
}
// barrdghupd_c
arma::mat barrdghupd_c(arma::mat X, arma::vec B, double lambda, arma::vec weight);
RcppExport SEXP bar_barrdghupd_c(SEXP XSEXP, SEXP BSEXP, SEXP lambdaSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(barrdghupd_c(X, B, lambda, weight));
    return rcpp_result_gen;
END_RCPP
}
// barridge_c
List barridge_c(arma::mat A, arma::vec B, int n, int p, arma::vec xiv, arma::vec lambdav, double abstol, int stand);
RcppExport SEXP bar_barridge_c(SEXP ASEXP, SEXP BSEXP, SEXP nSEXP, SEXP pSEXP, SEXP xivSEXP, SEXP lambdavSEXP, SEXP abstolSEXP, SEXP standSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xiv(xivSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambdav(lambdavSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    Rcpp::traits::input_parameter< int >::type stand(standSEXP);
    rcpp_result_gen = Rcpp::wrap(barridge_c(A, B, n, p, xiv, lambdav, abstol, stand));
    return rcpp_result_gen;
END_RCPP
}
// barcd_c
List barcd_c(arma::mat A, arma::vec B, int n, int p, arma::vec xiv, arma::vec lambdav, double abstol);
RcppExport SEXP bar_barcd_c(SEXP ASEXP, SEXP BSEXP, SEXP nSEXP, SEXP pSEXP, SEXP xivSEXP, SEXP lambdavSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xiv(xivSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambdav(lambdavSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(barcd_c(A, B, n, p, xiv, lambdav, abstol));
    return rcpp_result_gen;
END_RCPP
}
// barcv_c
List barcv_c(arma::mat A, arma::vec B, List trainL, List testL, int n, int p, arma::vec xiv, arma::vec lambdav, int nfold, std::string method, double abstol, int stand);
RcppExport SEXP bar_barcv_c(SEXP ASEXP, SEXP BSEXP, SEXP trainLSEXP, SEXP testLSEXP, SEXP nSEXP, SEXP pSEXP, SEXP xivSEXP, SEXP lambdavSEXP, SEXP nfoldSEXP, SEXP methodSEXP, SEXP abstolSEXP, SEXP standSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< List >::type trainL(trainLSEXP);
    Rcpp::traits::input_parameter< List >::type testL(testLSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xiv(xivSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambdav(lambdavSEXP);
    Rcpp::traits::input_parameter< int >::type nfold(nfoldSEXP);
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    Rcpp::traits::input_parameter< int >::type stand(standSEXP);
    rcpp_result_gen = Rcpp::wrap(barcv_c(A, B, trainL, testL, n, p, xiv, lambdav, nfold, method, abstol, stand));
    return rcpp_result_gen;
END_RCPP
}
// zbar_c
arma::mat zbar_c(arma::mat zmat, arma::vec Cv, double s);
RcppExport SEXP bar_zbar_c(SEXP zmatSEXP, SEXP CvSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type zmat(zmatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Cv(CvSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(zbar_c(zmat, Cv, s));
    return rcpp_result_gen;
END_RCPP
}
// zbarsum_c
arma::mat zbarsum_c(arma::mat zmat, arma::vec Cv, arma::vec t);
RcppExport SEXP bar_zbarsum_c(SEXP zmatSEXP, SEXP CvSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type zmat(zmatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Cv(CvSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(zbarsum_c(zmat, Cv, t));
    return rcpp_result_gen;
END_RCPP
}
// V_c
arma::mat V_c(arma::mat zmat, arma::vec Cv);
RcppExport SEXP bar_V_c(SEXP zmatSEXP, SEXP CvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type zmat(zmatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Cv(CvSEXP);
    rcpp_result_gen = Rcpp::wrap(V_c(zmat, Cv));
    return rcpp_result_gen;
END_RCPP
}
// Omega_c
arma::mat Omega_c(arma::mat zmat, List tlist, arma::vec Cv);
RcppExport SEXP bar_Omega_c(SEXP zmatSEXP, SEXP tlistSEXP, SEXP CvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type zmat(zmatSEXP);
    Rcpp::traits::input_parameter< List >::type tlist(tlistSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Cv(CvSEXP);
    rcpp_result_gen = Rcpp::wrap(Omega_c(zmat, tlist, Cv));
    return rcpp_result_gen;
END_RCPP
}
// P_c
arma::mat P_c(arma::mat zmat, List tlist, List olist, arma::vec Cv);
RcppExport SEXP bar_P_c(SEXP zmatSEXP, SEXP tlistSEXP, SEXP olistSEXP, SEXP CvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type zmat(zmatSEXP);
    Rcpp::traits::input_parameter< List >::type tlist(tlistSEXP);
    Rcpp::traits::input_parameter< List >::type olist(olistSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Cv(CvSEXP);
    rcpp_result_gen = Rcpp::wrap(P_c(zmat, tlist, olist, Cv));
    return rcpp_result_gen;
END_RCPP
}
// reggenint_c
List reggenint_c(arma::vec Cv, arma::mat Zmat, arma::vec olist, arma::vec tlist, arma::uvec mind);
RcppExport SEXP bar_reggenint_c(SEXP CvSEXP, SEXP ZmatSEXP, SEXP olistSEXP, SEXP tlistSEXP, SEXP mindSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Cv(CvSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Zmat(ZmatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type olist(olistSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tlist(tlistSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type mind(mindSEXP);
    rcpp_result_gen = Rcpp::wrap(reggenint_c(Cv, Zmat, olist, tlist, mind));
    return rcpp_result_gen;
END_RCPP
}
// bar_rcpparma_hello_world
arma::mat bar_rcpparma_hello_world();
RcppExport SEXP bar_bar_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(bar_rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP bar_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP bar_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP bar_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}
// ridreg_c
arma::vec ridreg_c(arma::mat A, arma::vec B, arma::vec binital, double lambda, arma::vec weight, double abstol);
RcppExport SEXP bar_ridreg_c(SEXP ASEXP, SEXP BSEXP, SEXP binitalSEXP, SEXP lambdaSEXP, SEXP weightSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type binital(binitalSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(ridreg_c(A, B, binital, lambda, weight, abstol));
    return rcpp_result_gen;
END_RCPP
}
// ridreg0_c
arma::vec ridreg0_c(arma::mat X, arma::vec Y, arma::vec binital, double lambda, arma::vec weight, double abstol);
RcppExport SEXP bar_ridreg0_c(SEXP XSEXP, SEXP YSEXP, SEXP binitalSEXP, SEXP lambdaSEXP, SEXP weightSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type binital(binitalSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(ridreg0_c(X, Y, binital, lambda, weight, abstol));
    return rcpp_result_gen;
END_RCPP
}
// ridregi_c
arma::vec ridregi_c(arma::mat X, arma::vec Y, double lambda, arma::vec weight, double abstol);
RcppExport SEXP bar_ridregi_c(SEXP XSEXP, SEXP YSEXP, SEXP lambdaSEXP, SEXP weightSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(ridregi_c(X, Y, lambda, weight, abstol));
    return rcpp_result_gen;
END_RCPP
}
// ridreg1_c
arma::vec ridreg1_c(arma::mat X, arma::vec Y, double lambda, arma::vec weight, double abstol);
RcppExport SEXP bar_ridreg1_c(SEXP XSEXP, SEXP YSEXP, SEXP lambdaSEXP, SEXP weightSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(ridreg1_c(X, Y, lambda, weight, abstol));
    return rcpp_result_gen;
END_RCPP
}
// ridreg2_c
arma::vec ridreg2_c(arma::mat A, arma::vec B, double lambda, arma::vec weight, double abstol);
RcppExport SEXP bar_ridreg2_c(SEXP ASEXP, SEXP BSEXP, SEXP lambdaSEXP, SEXP weightSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(ridreg2_c(A, B, lambda, weight, abstol));
    return rcpp_result_gen;
END_RCPP
}
// ridreg3_c
arma::vec ridreg3_c(arma::mat A, arma::vec B, double lambda, arma::vec weight, double abstol);
RcppExport SEXP bar_ridreg3_c(SEXP ASEXP, SEXP BSEXP, SEXP lambdaSEXP, SEXP weightSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(ridreg3_c(A, B, lambda, weight, abstol));
    return rcpp_result_gen;
END_RCPP
}
// ridreg4_c
arma::vec ridreg4_c(arma::mat A, arma::vec B, double lambda, arma::vec weight, double abstol);
RcppExport SEXP bar_ridreg4_c(SEXP ASEXP, SEXP BSEXP, SEXP lambdaSEXP, SEXP weightSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(ridreg4_c(A, B, lambda, weight, abstol));
    return rcpp_result_gen;
END_RCPP
}
// ridregM1_c
arma::mat ridregM1_c(arma::mat A, arma::vec B, double lambda, arma::vec weight, double abstol);
RcppExport SEXP bar_ridregM1_c(SEXP ASEXP, SEXP BSEXP, SEXP lambdaSEXP, SEXP weightSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(ridregM1_c(A, B, lambda, weight, abstol));
    return rcpp_result_gen;
END_RCPP
}
// ridregM2_c
arma::vec ridregM2_c(arma::mat A, arma::vec B, double lambda, arma::vec weight, int r, double abstol);
RcppExport SEXP bar_ridregM2_c(SEXP ASEXP, SEXP BSEXP, SEXP lambdaSEXP, SEXP weightSEXP, SEXP rSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(ridregM2_c(A, B, lambda, weight, r, abstol));
    return rcpp_result_gen;
END_RCPP
}
// ridregM3_c
arma::vec ridregM3_c(arma::mat X, arma::vec B, double lambda, arma::vec weight, int r, double abstol);
RcppExport SEXP bar_ridregM3_c(SEXP XSEXP, SEXP BSEXP, SEXP lambdaSEXP, SEXP weightSEXP, SEXP rSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(ridregM3_c(X, B, lambda, weight, r, abstol));
    return rcpp_result_gen;
END_RCPP
}
// ridregM4_c
arma::mat ridregM4_c(arma::mat X, arma::vec B, double lambda, arma::vec weight, int r, double abstol);
RcppExport SEXP bar_ridregM4_c(SEXP XSEXP, SEXP BSEXP, SEXP lambdaSEXP, SEXP weightSEXP, SEXP rSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(ridregM4_c(X, B, lambda, weight, r, abstol));
    return rcpp_result_gen;
END_RCPP
}
// test1
arma::mat test1(arma::mat A, arma::vec t);
RcppExport SEXP bar_test1(SEXP ASEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(test1(A, t));
    return rcpp_result_gen;
END_RCPP
}
// test2
arma::mat test2(arma::mat A, arma::vec t);
RcppExport SEXP bar_test2(SEXP ASEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(test2(A, t));
    return rcpp_result_gen;
END_RCPP
}
// sfun_c
double sfun_c(double z, double lambda);
RcppExport SEXP bar_sfun_c(SEXP zSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(sfun_c(z, lambda));
    return rcpp_result_gen;
END_RCPP
}
// scadupdate_c
double scadupdate_c(double z, double lambda, double gam, double wjj);
RcppExport SEXP bar_scadupdate_c(SEXP zSEXP, SEXP lambdaSEXP, SEXP gamSEXP, SEXP wjjSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type gam(gamSEXP);
    Rcpp::traits::input_parameter< double >::type wjj(wjjSEXP);
    rcpp_result_gen = Rcpp::wrap(scadupdate_c(z, lambda, gam, wjj));
    return rcpp_result_gen;
END_RCPP
}
// elnetupdate_c
double elnetupdate_c(double z, double lambda, double alpha, double wjj);
RcppExport SEXP bar_elnetupdate_c(SEXP zSEXP, SEXP lambdaSEXP, SEXP alphaSEXP, SEXP wjjSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type wjj(wjjSEXP);
    rcpp_result_gen = Rcpp::wrap(elnetupdate_c(z, lambda, alpha, wjj));
    return rcpp_result_gen;
END_RCPP
}
// penreg_c
List penreg_c(arma::mat A, arma::vec B, int n, int p, arma::vec lambdav, arma::vec weight, double alpha, double gamma, std::string method, double abstol);
RcppExport SEXP bar_penreg_c(SEXP ASEXP, SEXP BSEXP, SEXP nSEXP, SEXP pSEXP, SEXP lambdavSEXP, SEXP weightSEXP, SEXP alphaSEXP, SEXP gammaSEXP, SEXP methodSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambdav(lambdavSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(penreg_c(A, B, n, p, lambdav, weight, alpha, gamma, method, abstol));
    return rcpp_result_gen;
END_RCPP
}
// pencv_c
List pencv_c(arma::mat A, arma::vec B, List trainL, List testL, int n, int p, arma::vec lambdav, arma::vec weight, int nfold, double alpha, double gamma, std::string method, double abstol);
RcppExport SEXP bar_pencv_c(SEXP ASEXP, SEXP BSEXP, SEXP trainLSEXP, SEXP testLSEXP, SEXP nSEXP, SEXP pSEXP, SEXP lambdavSEXP, SEXP weightSEXP, SEXP nfoldSEXP, SEXP alphaSEXP, SEXP gammaSEXP, SEXP methodSEXP, SEXP abstolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< List >::type trainL(trainLSEXP);
    Rcpp::traits::input_parameter< List >::type testL(testLSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambdav(lambdavSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< int >::type nfold(nfoldSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    rcpp_result_gen = Rcpp::wrap(pencv_c(A, B, trainL, testL, n, p, lambdav, weight, nfold, alpha, gamma, method, abstol));
    return rcpp_result_gen;
END_RCPP
}
// test
arma::vec test(arma::vec a, double b);
RcppExport SEXP bar_test(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(test(a, b));
    return rcpp_result_gen;
END_RCPP
}
// rss_c
arma::vec rss_c(arma::mat A, arma::vec B, arma::mat beta);
RcppExport SEXP bar_rss_c(SEXP ASEXP, SEXP BSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(rss_c(A, B, beta));
    return rcpp_result_gen;
END_RCPP
}
// rssf_c
arma::vec rssf_c(arma::mat A, arma::vec B, arma::mat beta, double tss);
RcppExport SEXP bar_rssf_c(SEXP ASEXP, SEXP BSEXP, SEXP betaSEXP, SEXP tssSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type tss(tssSEXP);
    rcpp_result_gen = Rcpp::wrap(rssf_c(A, B, beta, tss));
    return rcpp_result_gen;
END_RCPP
}
// standard_c
List standard_c(arma::mat A, arma::vec B, int n);
RcppExport SEXP bar_standard_c(SEXP ASEXP, SEXP BSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(standard_c(A, B, n));
    return rcpp_result_gen;
END_RCPP
}
