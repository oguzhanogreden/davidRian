
<!-- README.md is generated from README.Rmd. Please edit that file -->
Davidian curves in R
====================

A Davidian curve \[1, 2\] is a seminonparametric density representation, which defines flexible and easy to estimate densities. You can use Davidian curves for relaxing the distributional assumptions of random variables. For instance, \[1\] used this representation in linear mixed modeling context for relaxing the normal distribution assumption for random variables. In another case, \[3\] used Davidian curves in Item Response Theory setting. In this case, the Davidian curve replaced the standard normal density used in marginal maximum likelihood estimation.

This package provides a collection of function which are useful for making Davidian curves accessible for use in a context independent manner. The package makes three function available to users, following the parameterization provided by \[1\], illustrated by \[2\] and reproduced here:

1.  `ddc(x, phi)`, Davidian curve density function
2.  `rdc(n, phi)`, a random sampler for Davidian curves
3.  `dcGrad(x, phi)`, the gradient function of Davidian curves

The second argument `phi` of these functions is a vector, which contains the parameter(s) defining the shape of the curve. Following \[2\], this package provides support for Davidian curves up to 10 parameters. The standard normal case, where `length(phi) == 0` holds is not implemented.

References
==========

1.  Gallant, A. R., & Nychka, D. W. (1987). Semi-nonparametric maximum likelihood estimation. Econometrica: Journal of the Econometric Society, 363-390.
2.  Zhang, D., & Davidian, M. (2001). Linear mixed models with flexible distributions of random effects for longitudinal data. Biometrics, 57(3), 795-802.
3.  Woods, C. M., & Lin, N. (2009). Item response theory with estimation of the latent density using Davidian curves. Applied Psychological Measurement, 33(2), 102-117.
