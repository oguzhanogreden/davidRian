ie <- imitate_bm_e(121, 1000)
est <- thetas_to_par_Etable(ie[,"qp"], ie[,"Nq"], 3, dc_E_LL, dc_E_GR)

source("../lib.R")
# 2. Imitate E-table using the density

R <- 100 # Simulation runs
N <- 1000
K <- c(2, 4, 6, 8)
QP <- c(61, 121, 241)

sim_pars <- expand.grid(1:R, N, K, QP)
sim_pars[,"Iter"] <- 1:nrow(sim_pars)

res <- apply(sim_pars, 1, function(x) {
  if (x["Iter"] %% 50 == 0) message(round(x["Iter"]/nrow(sim_pars), 2))
  
  seed <- runif(1, 1, 1e+9)
  set.seed(seed)
  
  n <- x["Var2"]
  k <- x["Var3"]
  qp <- x["Var4"]
  etable <- imitate_bm_e(qp, n)
  
  est <- thetas_to_par(dat, k, dc_LL, dc_GR)
  
  c(seed = seed, N=n, K=k, par = c(est$pa, rep(NA, 8-k)), conv = est$convergence)
})

saveRDS(as.data.frame(t(res)), file="Test2.RDS")