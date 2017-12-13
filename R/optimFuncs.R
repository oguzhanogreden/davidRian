# DC estimation functions

dc.LL <- function(thetas, k, mean, sd, phi) {
  dens <- dcdensC(thetas, k, mean, sd, phi)

  sum(log(dens))
}

dc.LL.grad <- function (thetas, k, mean, sd, phi) {
  primGrad <- sapply(thetas, function(x) {
      dcGrad_(x, k, cMat(k, phi), phi)
    })
  
  primGrad <- t(primGrad)

  densDenom <- sqrt(dcdensC(thetas, k, mean, sd, phi) / dnorm(thetas, mean = mean, sd = sd))
  primGrad/densDenom
}
