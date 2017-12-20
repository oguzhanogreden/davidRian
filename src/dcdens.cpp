# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
# include <cmath>

using namespace Rcpp;

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

double dcdens_ (double x, int k, double mean, double sd, NumericVector phi) {
  arma::mat invB(k+1, k+1);
  arma::mat c(k+1, 1);
  arma::mat res(k, k);
  
  double normdens;
  
  // invB
  invB = invBMat(k);
  
  // Create c
  c = cMat(k, phi);
  
  res = arma::pow((invB * c).t() * expVec(x, k), 2);
  normdens = dnorm(NumericVector::create(x), mean, sd, true)[0];
  
  return exp(log(res[0]) + normdens);
}

//' Density function for Davidian curves
//'
//' Returns the density for a vector of x.
//'
//' @param x An integer vector
//' @name dcdens
// [[Rcpp::export]]
NumericVector dcdensC (NumericVector x, int k, double mean, double sd, NumericVector phi) {
  NumericVector res(x.length());

  for (int i = 0; i < x.length(); i++) {
    res[i] = dcdens_(x[i], k, mean, sd, phi);
  }

  return res;
}

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

// [[Rcpp::export]]
NumericVector rdc (int n, int k, double mean, double sd, NumericVector phi) {
  NumericVector out(n);
  NumericVector c(1, 8.55);
  int accepted = 0;
  
  NumericVector y(1);
  NumericVector u(1);
  NumericVector ratio(1);
  
  while (accepted < n) {
    y = runif(1, mean-10, mean+10); // proposal density is runif, domain is defined as [mean-10, mean+10].
    u = runif(1);
    
    ratio = dcdensC(y, k, mean, sd, phi) / (c * dunif(y, mean-10, mean+10)); //  ;
    
    if (is_true(all(u < ratio))) {
      out[accepted] = y[0];
      accepted++;
    
    }
  }
  
  return out;
}
/*** R
rdcdens <- function(n, k, mean, sd, phi) {
  res <- c()
  c <- 8.55
  
  while (length(res) < n) {
    y <- runif(1, mean-10, mean+10) # proposal density is runif, domain is defined as [mean-10, mean+10].
    u <- runif(1)
    
    if (u < dcdensC(y, k, mean, sd, phi)/(c * dunif(y, mean-10, mean+10))) {
      res <- c(res, y)
    }
  }
  
  res
}

testcpp <- rdc(2000, 1, 0, 5, 1)
testr <- rdcdens(2000, 1, 0, 5, 1)

*/



