# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
# include <cmath>
# include "dcUtils.h"

// using namespace Rcpp;

arma::mat expVec(double x, int deg) {
  arma::mat out(deg+1, 1);
  
  for (int i = 0; i < deg+1; i++) {
    out(i, 0) = pow(x, i);
  }
  
  return out;
}

arma::mat cMat (int k, NumericVector phi) {
  arma::mat prim_c(1,1);
  arma::mat c(1, 1);
  
  prim_c[1] = 1;
  
  // Fill primitive prim_c matrix
  if (k != 0) {
    prim_c = arma::mat(k+1, k);
    c = arma::mat(k+1, 1);
    
    for (int i = 0; i < prim_c.n_rows; i++) {
      if (i + 1 == prim_c.n_rows) { // if last row
        // write all cos's
        for (int j = 0; j < prim_c.n_cols; j++) {
        prim_c(i, j) = cos(phi[j]);
        }
      } else { // if not last row
        for (int j = 0; j < prim_c.n_cols; j++) {
        if (j > i) {
        prim_c(i, j) = 1;
        } else if (j == i) { // if last col
        prim_c(i, j) = sin(phi[j]);
        } else { // if not last col
        prim_c(i, j) = cos(phi[j]);
        }
        }
      }
  }
}
      
      // Collapse cols by multiplication to create c matrix
      for (int i = 0; i < prim_c.n_rows; i++) {
      double mult = 1;
      
      for (int j = 0; j < prim_c.n_cols; j++) {
      mult = mult * prim_c(i, j);
      }
      
      c(i, 0) = mult;
      }
      
      return c;
}
    
arma::mat invBMat (int k) {
    arma::mat B(k+1, k+1);
    arma::mat invB(k+1, k+1);
    NumericMatrix M(k+1, k+1);
    NumericVector vals = NumericVector::create(1);
    int MSize = (k+1) * (k+1);
    
    // M matrix values
    if (k == 1) {
    vals = NumericVector::create(1, 0, 0, 1);
    } else if (k == 2) {
    vals = NumericVector::create(1, 0, 1, 0, 1, 0, 1, 0, 3);
    } else if (k == 3) {
    vals = NumericVector::create(1, 0, 1, 0, 0, 1, 0, 3, 1, 0, 3, 0, 0, 3, 0, 15);
    }
    
    // Fill M
    for (int i = 0; i < MSize; i++) {
    M[i] = vals[i];
    }
    
    B = arma::chol(Rcpp::as<arma::mat>(M));
    invB = arma::inv(B);
    
    return invB;
    }