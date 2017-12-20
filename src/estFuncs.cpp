# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
# include <cmath>
# include "dcUtils.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dcGrad_ (double x, int k, NumericMatrix c, NumericVector phi) {
  arma::mat invB(k+1, k+1);
  arma::mat cDer(k, k+1);
  double tp;
  arma::mat tmp1(k, 1);

  for (int i = 0; i < cDer.n_rows; i++) {
    tp = tan(phi[i]);

    for (int j = 0; j < cDer.n_cols; j++) {
      if (j == cDer.n_cols) {
        cDer(i,j) = -1 * c[j] * tp;
      } else if (j > i) {
        cDer(i,j) = -1 * c[j] * tp;
      } else if (j == i) {
        cDer(i,j) = c[j] * (1 / tp);
      } else if (j < i) {
        cDer(i,j) = 0;
      } else {
        cDer(i,j) = -999;
      }
    }
  }

  invB = invBMat(k).t();

  tmp1 = (2* cDer * invB * expVec(x, k));

  return Rcpp::wrap(tmp1); // remains to be divided by P_k of (5) in Woods & Lin
}

// NumericVector dcGradC (NumericVector x, int k, NumericMatrix c, NumericVector phi) {
  //   NumericVector res()
  // }