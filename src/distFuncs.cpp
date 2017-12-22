# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
# include <cmath>
# include "dcUtils.h"

using namespace Rcpp;

double ddc_ (double x, double mean, double sd, NumericVector phi) {
  int k = phi.length();
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

//' Density function for univariate Davidian Curves (DC)
//'
//' Returns the density for a vector of x.
//'
//' @param x A vector of observations.
//' @param mean Mean of the normal base of DC, see the package vignette.
//' @param sd SD of the normal base of DC, see the package vignette.
//' @param phi DC parameters as introduced in Woods & Lin.
//' 
// [[Rcpp::export]]
NumericVector ddc (NumericVector x, double mean, double sd, NumericVector phi) {
  
  
  if (phi.length() > 8) {
    stop("k > 8 is not supported.");
  }
  
  NumericVector res(x.length());

  for (int i = 0; i < x.length(); i++) {
    res[i] = ddc_(x[i], mean, sd, phi);
  }

  return res;
}

//' Random samples from a univariate Davidian Curve (DC)
//'
//' Returns n samples from a univariate DC.
//'
//' @param n Number of observations to be sampled.
//' @param mean Mean of the normal base of DC, see the package vignette.
//' @param sd SD of the normal base of DC, see the package vignette.
//' @param phi DC parameters as introduced in Woods & Lin.
//' 
// [[Rcpp::export]]
NumericVector rdc (int n, double mean, double sd, NumericVector phi) {
  
  if (phi.length() > 8) {
    stop("k > 8 is not supported.");
  }
  
  NumericVector out(n);
  NumericVector c(1, 8.55);
  int accepted = 0;
  
  NumericVector y(1);
  NumericVector u(1);
  NumericVector ratio(1);
  
  while (accepted < n) {
    y = runif(1, mean-10, mean+10); // proposal density is runif, domain is defined as [mean-10, mean+10].
    u = runif(1);
    
    ratio = ddc(y, mean, sd, phi) / (c * dunif(y, mean-10, mean+10)); //  ;
    
    if (is_true(all(u < ratio))) {
      out[accepted] = y[0];
      accepted++;
    
    }
  }
  
  return out;
}

