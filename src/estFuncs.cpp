# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
# include <cmath>
# include "dcUtils.h"

using namespace Rcpp;

NumericVector dcGrad_ (double x, NumericVector phi) {
  int k = phi.length(); 
  // Eq17 terms
  arma::mat invB(k+1, k+1);
  arma::mat c;
  arma::mat cDer(k, k+1);
  double tp;
  arma::mat res(k, 1);
  
  //Eq5 term
  arma::mat pk;
  
  c = cMat(k, phi);
  
  for (unsigned int i = 0; i < cDer.n_rows; i++) {
    tp = tan(phi[i]);
    
    for (unsigned int j = 0; j < cDer.n_cols; j++) {
      if (j == cDer.n_cols) {
        cDer(i,j) = -1 * c[j] * tp;
      } else if (j > i) {
        cDer(i,j) = -1 * c[j] * tp;
      } else if (j == i) {
        cDer(i,j) = c[j] * (1 / tp);
      } else if (j < i) {
        cDer(i,j) = 0;
      }
    }
  }
  
  invB = invBMat(k);
  
  res = (2 * cDer * invB.t() * expVec(x, k));
  
  // res, remains to be divided by P_k of Equation 5 in Woods & Lin, see Equation 17.
  // P_k is obtained as follows:
  pk = (invB * c).t() * expVec(x, k);
  
  // The returned value needs to be multiplied by N(theta_q) from E step. mirt provides the N(), hence it cannot done here...
  return Rcpp::wrap(res / pk[0]); 
}

//' Gradient of the log-likelihood of univariate Davidian curves
//'
//' Provides the gradient for use in estimation.
//'
//' @param x A vector of observations.
//' @param phi Davidian curve parameters. length(phi) < 11.
// [[Rcpp::export]]
NumericVector dcGrad (NumericVector x,  NumericVector phi) {
  NumericMatrix res(x.length(), phi.length());
  NumericMatrix::Row tmprow = res(1, _);
  
  for (int i = 0; i < x.length(); i++) {
    res.row(i) = dcGrad_(x[i], phi);
  }
  
  return res;
}