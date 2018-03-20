library(cowplot)
library(davidRian)
library(plyr)
library(ggplot2)

# Functions ####

# LL and Gradient for Isolated Estimation

dc_LL <- function (phi, thetas) {
  dens <- ddc(thetas, mean=0, sd=1, phi)
  sum(log(dens))
}

dc_GR <- function (phi, thetas) {
  colSums(dcGrad(thetas, phi))
}

dc_E_LL <- function (phi, qp, Nq) {
  dens <- ddc(qp, mean=0, sd=1, phi)
  sum(Nq * log(dens))
}

dc_E_GR <- function (phi, qp, Nq) {
  colSums(Nq * dcGrad(qp, phi))
}
  
# Given a vector of thetas, number of Davidian curve parameters, LL and gradient function, return parameters
thetas_to_par <- function(thetas, k, ll, gr) {
  dc_init_phi <- rep(77, k)
  
  optim(dc_init_phi, fn = ll, gr = gr, thetas = thetas,
        control = list(fnscale = -1, trace = 0), method = "L-BFGS-B", lower = -89.99, upper = 90)
}

thetas_to_par_Etable <- function(qp, Nq, k, ll, gr) {
  dc_init_phi <- rep(77, k)
  optim(dc_init_phi, fn = ll, gr = gr, qp = qp, Nq = Nq,
        control = list(fnscale = -1, trace = 0), method = "L-BFGS-B", lower = -89.99, upper = 90)
}

# Density and sampler from bimodal in Woods & Lin (2009):
dens_bm <- function(x) {
  dnorm(x, mean = -.70, sd = .50) * .6 + dnorm(x, mean = 1.05, sd = .54) * .4
}

sample_bm <- function(N) {
  c(rnorm(N*.6, mean = -.70, sd = .50), rnorm(N*.4, mean = 1.05, sd = .54))
}

# Imitate E-table by density:
imitate_bm_e <- function(quadpts, N) {
  qp <- seq(-6, 6, length.out = quadpts)
  
  cbind(qp=qp, Nq=dens_bm(qp) / sum(dens_bm(qp)) * N)
}

# From mirt:
interpolateWoods <- function (point, std_point, nq1, nq2, delta) {
  term1 <- (point - std_point)/delta
  term2 <- nq2 - nq1
  
  term1 * term2 + nq1
}

extrapolateWoods <- function (point, std_point, nq1, nq2, delta, tail) {
  if (tail == "left") {
    ((nq1 / nq2) ^ ((std_point - point) / delta)) * nq1
  } else if (tail == "right") {
    ((nq2 / nq1) ^ ((point - std_point) / delta)) * nq2
  }
}

standardizeQuadrature <- function (qp, std_qp, nq) {
  # This gives almost identical results to approx().
  # Except for the tails, where this is a bit more smoothly decreasing.
  min_stdqp <- min(std_qp)
  max_stdqp <- max(std_qp)
  
  
  res <- matrix(NA, nrow = length(qp), ncol = 2)
  res[, 1] <- qp
  delta <- qp[2] - qp[1] # is the distance between any two q or, equivalently, between any two q*
  for (i in 1:length(qp)) {
    if (qp[i] <= min_stdqp) {
      res[i,2] <- extrapolateWoods(qp[i], min_stdqp, nq[1], nq[2], delta, "left")
    } else if (max_stdqp <= qp[i]) {
      res[i,2] <- extrapolateWoods(qp[i], max_stdqp, nq[length(nq)-1], nq[length(nq)], delta, "right")
    } else {
      std_ind <- max(which(qp[i] > std_qp))
      res[i,2] <- interpolateWoods(qp[i], std_qp[std_ind], nq[std_ind], nq[std_ind+1], delta)
    }
  }
  
  res[,2] <- res[,2] / sum(res[,2]) * sum(nq)
  res
}
