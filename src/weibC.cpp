// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// //' Weibull distribution negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each Weibull parameter
// //' @param X1 a design matrix for the Weibull log scale parameter
// //' @param X2 a design matrix for the Weibull log shape parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return weibd0 a scalar, the negative log-liklihood
// //' @return weibd12 a matrix, first then second derivatives w.r.t. Weibull parameters
// //' @return weibd34 a matrix, third then fourth derivatives w.r.t. Weibull parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double weibd0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, arma::vec yvec, const arma::uvec& dupid, int dcate)
{
    
arma::vec llambdavec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lkvec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = yvec.size();

if (dcate == 1) {
    llambdavec = llambdavec.elem(dupid);
    lkvec = lkvec.elem(dupid);
}

double y, llambda, lk, lambda, k;
double nllh=0.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
llambda = llambdavec[j];
lk = lkvec[j];
lambda = exp(llambda);
k = exp(lk);

nllh += -(lk - llambda + (k - 1.0) * log(y / lambda) - R_pow(y / lambda, k));

}

return(nllh);

}

// //' @rdname weibd0
// [[Rcpp::export]]
arma::mat weibd12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::vec yvec, const arma::uvec dupid, int dcate)
{
    
arma::vec llambdavec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lkvec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 5);

if (dcate == 1) {
    llambdavec = llambdavec.elem(dupid);
    lkvec = lkvec.elem(dupid);
}

double y, ll, lk;
double ee1, ee2, ee3, ee4, ee6, ee7, ee8, ee9; 

for (int j=0; j < nobs; j++) {

y = yvec[j];
ll = llambdavec[j];
lk = lkvec[j];

ee1 = exp(lk);
ee2 = exp(ll);
ee3 = y/ee2;
ee4 = ee1 - 1;
ee6 = log(y) - ll;
ee7 = R_pow(ee3, ee4);
ee8 = R_pow(ee3, ee1);
ee9 = ee1 * ee6;

out(j, 0) = (1 - y * ee7/ee2) * ee1;
out(j, 1) = -(1 + (1 - ee8) * ee1 * ee6);

out(j, 2)= ee4 * ((y * y)/(ee2 * ee2 * ee3 * ee3) - 1) + y * ee1 * (y * ee4 * R_pow(ee3, ee1 - 2)/ee2 + ee7)/ee2;
out(j, 3) = (1 - y * (ee9 * ee7 + ee7)/ee2) * ee1;
out(j, 4) = -((1 - (ee9 * ee8 + ee8)) * ee1 * ee6);

}

return out;

}

// //' @rdname weibd0
// [[Rcpp::export]]
arma::mat weibd34(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::vec yvec, const arma::uvec dupid, int dcate)
{
    
arma::vec llambdavec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lkvec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 9);

if (dcate == 1) {
    llambdavec = llambdavec.elem(dupid);
    lkvec = lkvec.elem(dupid);
}

double y, ll, lk;
double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee9; 
double ee10, ee11, ee12, ee13, ee14, ee16, ee17, ee18, ee19; 
double ee21, ee22, ee23, ee24, ee26, ee27, ee29; 
double ee33, ee34, ee35, ee40; 

for (int j=0; j < nobs; j++) {

y = yvec[j];
ll = llambdavec[j];
lk = lkvec[j];

ee1 = exp(lk);
ee2 = exp(ll);
ee3 = y/ee2;
ee4 = ee1 - 1;
ee5 = ee1 - 2;
ee6 = R_pow(ee3, ee4);
ee7 = R_pow(ee3, ee5);
ee9 = log(y) - ll;
ee10 = R_pow(ee3, ee1);
ee11 = ee1 - 3;
ee12 = R_pow(ee3, ee11);
ee13 = ee3 * ee3;
ee14 = ee1 * ee9;
ee16 = ee2 * ee2 * ee13;
ee17 = y * y;
ee18 = ee14 * ee6;
ee19 = ee18 + ee6;
ee21 = ee19 + ee6 + ee6;
ee22 = ee17/ee16;
ee23 = y * ee5;
ee24 = 1/ee13;
ee26 = 2 * ee22 - 3;
ee27 = ee4 * ee7;
ee29 = ee1 * ee21 * ee9;
ee33 = ee14 * ee10 + ee10 + ee10 + ee10;
ee34 = y * ee1;
ee35 = y * ee4;
ee40 = ee23 * ee12/ee2 + ee7 + ee7 + ee7;

// third derivatives
// 1=log(scale), 2=log(shape)
// order: 111, 112, 122, 222

out(j, 0) = (1 + ee17 * ee26/ee16) * ee4 - ee34 * (ee35 * 
    ee40/ee2 + ee6)/ee2;
out(j, 1) = ee1 * (y * (ee18 + y * 
    (ee24 + ee27 + ee1 * (ee4 * ee9 * ee7 + ee7))/ee2 + ee6)/ee2 - 1);
out(j, 2) = (1 - y * (ee29 + ee6)/ee2) * ee1;
out(j, 3) = -((1 - (ee1 * ee33 * ee9 + ee10)) * ee1 * ee9);

// fourth derivatives
// 1=log(scale), 2=log(shape)
// order: 1111, 1112, 1122, 1222, 2222

out(j, 4) = ee4 * 
    (ee17 * (7 + ee17 * (8 * ee22 - 14)/ee16)/ee16 - 1) + 
    ee34 * (ee35 * (ee23 * (y * ee11 * R_pow(ee3, ee1 - 4)/ee2 + 
    ee12 + ee12 + ee12 + ee12 + ee12 + ee12)/ee2 + ee7 + 
    ee7 + ee7 + ee7 + ee7 + ee7 + ee7)/ee2 + ee6)/ee2;
out(j, 5) = (1 + y * (y * (ee26/ee13 - (((2 * ee7 + 
    ee7) * ee4 * ee9 + y * (ee4 * (ee5 * ee9 * ee12 + 
    ee12) + ee5 * ee12)/ee2 + ee7 + ee7 + ee7) * ee1 + 
    ee4 * ee40))/ee2 - ee19)/ee2) * ee1;
out(j, 6) = -ee1 * 
    (y * (ee29 + y * (ee24 + ((2 * (ee1 * ee7) + ee4 * 
    (ee14 * ee7 + ee7 + ee7 + ee7)) * ee9 + ee7 + 
    ee7 + ee7) * ee1 + ee27)/ee2 + ee6)/ee2 - 1);
out(j, 7) = (1 - y * (ee1 * (ee1 * (ee21 + ee6 + ee6 + 
    ee6) * ee9 + ee6 + ee6 + ee6 + ee6 + ee6 + ee6 + 
    ee6) * ee9 + ee6)/ee2) * ee1;
out(j, 8) = -((1 - 
    (ee1 * (ee1 * (ee33 + ee10 + ee10 + ee10) * ee9 + 
    ee10 + ee10 + ee10 + ee10 + ee10 + ee10 + ee10) * 
    ee9 + ee10)) * ee1 * ee9);

}

return out;

}
