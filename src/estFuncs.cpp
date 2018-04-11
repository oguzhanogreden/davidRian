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
        if (tp == 0) {
          cDer(i,j) = 0; // it is not clear how to go about this in Eq (18), hence 0.
        } else {
          cDer(i,j) = c[j] * (1/tp);
        }
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
//' @param phi phi Davidian curve parameters.
//' A maximum of 10 parameters is allowed, all of which should be between -90 < phi <= 90.
//' 
//' @details Woods & Lin (2009) provide the gradient (Equations 17 and 18). Although the parameterization is
//' defined for -90 < phi < 90, there are elements which are not defined when phi = 0.0. A warning is
//' raised for this.
//' 
//' @references Woods, C. M., & Lin, N. (2009). Item response theory with estimation of the latent density using Davidian curves.
//' \emph{Applied Psychological Measurement, 33}(2), 102-117. \doi{10.1177/0146621608319512}
//' 
//' @examples
//' # The loglikelihood of a univariate Davidian curve is given by,
//' dc_LL <- function(phi, dat) {
//'   sum(log(ddc(dat, phi)))
//' }
//' 
//' # dc_grad can be used for obtaining the gradient of this loglikelihood as follows:
//' dc_LL_GR <- function(phi, dat) {
//'   colSums(dc_grad(dat, phi))
//' }
//'
//' # This can be verified by numerical approximation.
//' # For instance, using numDeriv package:
//' \dontrun{
//' phi <- c(-5, 2.5, 10)
//' d <- runif(10, -5, 5)
//' dc_LL_GR(phi, d)
//' numDeriv::grad(dc_LL, x = phi, dat = d)
//' }
//' 
// [[Rcpp::export]]
NumericVector dc_grad (NumericVector x,  NumericVector phi) {
  NumericMatrix res(x.length(), phi.length());
  NumericMatrix::Row tmprow = res(1, _);
  
  if (phi.length() > 10) {
    stop("length(phi) > 10 is not supported.");
  }
  
  if (is_true(any((phi <= -90) | (phi > 90)))) {
    stop("90 < phi <= 90 should hold for all phi.");
  }
  
  if(is_true(any(phi == 0))) {
    warning("Ensure that not all phi's are 0 when dc_grad() is called.");
  }
  
  for (int i = 0; i < x.length(); i++) {
    res.row(i) = dcGrad_(x[i], phi);
  }
  
  return res;
}