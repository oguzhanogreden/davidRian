source("../lib.R")
# 1. Estimate parameters of DC using simulated data ####

R <- 100 # Simulation runs
N <- c(121, 250, 500, 1000)
K <- c(2, 4, 6, 8)

sim_pars <- expand.grid(1:R, N, K)
sim_pars[,"Iter"] <- 1:nrow(sim_pars)

res <- apply(sim_pars, 1, function(x) {
  if (x["Iter"] %% 50 == 0) message(round(x["Iter"]/nrow(sim_pars), 2))
  
  seed <- runif(1, 1, 1e+9)
  set.seed(seed)
  
  n <- x["Var2"]
  k <- x["Var3"]
  dat <- sample_bm(n)
  
  est <- thetas_to_par(dat, k, dc_LL, dc_GR)
  
  c(seed = seed, N=n, K=k, par = c(est$pa, rep(NA, 8-k)), conv = est$convergence)
})

saveRDS(as.data.frame(t(res)), file="Test1.RDS")