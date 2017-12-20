# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
# include <cmath>
# include "dcUtils.h"

using namespace Rcpp;

double ddc_ (double x, int k, double mean, double sd, NumericVector phi) {
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
// [[Rcpp::export]]
NumericVector ddc (NumericVector x, int k, double mean, double sd, NumericVector phi) {
  NumericVector res(x.length());

  for (int i = 0; i < x.length(); i++) {
    res[i] = ddc_(x[i], k, mean, sd, phi);
  }

  return res;
}

//' Sampling from a given Davidian curves
//'
//' Samples n realizations from the specified Davidian Curve.
//'
//' @param x An integer vector
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
    
    ratio = ddc(y, k, mean, sd, phi) / (c * dunif(y, mean-10, mean+10)); //  ;
    
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
    
    if (u < ddc(y, k, mean, sd, phi)/(c * dunif(y, mean-10, mean+10))) {
      res <- c(res, y)
    }
  }
  
  res
}

testcpp <- rdc(2000, 1, 0, .1, 1)
testr <- rdcdens(2000, 1, 0, 5, 1)
*/



