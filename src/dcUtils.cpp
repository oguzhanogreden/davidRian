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

arma::vec insZ (arma::vec vals) {
  arma::vec tmp(vals.n_elem * 2, arma::fill::zeros);
  
  for (int i = 0; i < vals.n_elem; i++) {
    tmp[i*2] = vals(i);
  }
  
  return tmp;
}

arma::mat fillM (arma::vec vals) {
  int k = vals.n_elem;
  arma::mat M(k, k, arma::fill::zeros);
  arma::rowvec valsZ;
  
  valsZ = insZ(vals).t();
  
  for (int i = 0; i < k; i++) {
    M.row(i) = valsZ.head(k);
    valsZ = arma::shift(valsZ, -1);
  }
  
  return M;
}

arma::mat invBMat (int k) {
    arma::mat B(k+1, k+1);
    arma::mat invB(k+1, k+1);
    arma::mat M(k+1, k+1);
    arma::vec vals(k+1);
    int MSize = (k+1) * (k+1);
    
    // M matrix values
    if (k == 1) {
    vals = arma::vec({1, 1});
    } else if (k == 2) {
    vals = arma::vec({1, 1, 3});
    } else if (k == 3) {
      vals = arma::vec({1, 1, 3, 15});
    } else if (k == 4) {
      vals = arma::vec({1, 1, 3, 15, 105});
    } else if (k == 5) {
      vals = arma::vec({1, 1, 3, 15, 105, 945});
    } else if (k == 6) {
      vals = arma::vec({1, 1, 3, 15, 105, 945, 10395});
    } else if (k == 7) {
      vals = arma::vec({1, 1, 3, 15, 105, 945, 10395, 135135});
    } else if (k == 8) {
      vals = arma::vec({1, 1, 3, 15, 105, 945, 10395, 135135, 2027025});
    }
    
    // Fill M
    M = fillM(vals);
    
    B = arma::chol(M);
    invB = arma::inv(B);
    
    return invB;
    }

