source("../lib.R")
# 4.i Imitate E-table using the density ####

K <- c(1:8)

sim_pars <- expand.grid(K)

load("e-table.Rdata")

res_i <- apply(sim_pars, 1, function(x) {
  seed <- runif(1, 1, 1e+9)
  set.seed(seed)
  
  k <- x["Var1"]
  
  est <- thetas_to_par_Etable(as.vector(gTheta[[1]]), rr, k, dc_E_LL, dc_E_GR)
  
  c(K=k, par = c(est$pa, rep(NA, 8-k)), conv = est$convergence)
})

saveRDS(as.data.frame(t(res_i)), file="Test4i.RDS")

# 4.ii Imitate E-table using density, standardized

obs_m <- sum(gTheta[[1]]) / sum(rr)
obs_sd <- sqrt(sum(rr * (gTheta[[1]] - obs_m)^2) / sum(rr))
gTheta_std <- (gTheta[[1]] - obs_m) / obs_sd
std_qp <- standardizeQuadrature(gTheta[[1]], gTheta_std, rr)

res_ii <- apply(sim_pars, 1, function(x) {
  seed <- runif(1, 1, 1e+9)
  set.seed(seed)
  
  k <- x["Var1"]
  est <- thetas_to_par_Etable(std_qp[,1], std_qp[,2], k, dc_E_LL, dc_E_GR)
  
  c(K=k, par = c(est$pa, rep(NA, 8-k)), conv = est$convergence)
})

saveRDS(as.data.frame(t(res_ii)), file="Test4ii.RDS")
