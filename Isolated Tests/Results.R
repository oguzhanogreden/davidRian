source("lib.R")

res1 <- readRDS("1/Test1.RDS")
res2 <- readRDS("2/Test2.RDS")
load("4/e-table.Rdata")
res4i <- readRDS("4/Test4i.RDS")
res4ii <- readRDS("4/Test4ii.RDS")

# 1. Simulated data ####
K <- unique(res1$K.Var3)
N <- unique(res1$N.Var2)
sim_pars <- expand.grid(N, K)

plots1 <- apply(sim_pars[order(sim_pars$Var1),], 1, function(x) {
  
  match_conds <- (res1[,2] == x[1]) & (res1[,3] == x[2])
  res_subs <- res1[match_conds, c(1:(x[2]+3), 12)]
  
  qp <- seq(-6, 6, length.out = 121)
  
  base_plot <- ggplot(data.frame(X = qp), aes(x=X)) +
    stat_function(fun = dens_bm, colour = "red", size=1.1) +
    ggtitle(paste("N =", x[1], "|", "# DC Par =", x[2])) +
    theme_minimal()
  
  for (i in 1:nrow(res_subs)) {
    phi <- res_subs[i, 4:(x[2]+3)]
    base_plot <- base_plot + stat_function(fun = ddc, args = list(mean = 0, sd = 1, phi = unlist(phi)),
                                           alpha = .05)
  }
  
  base_plot
})

plot_out1 <- plot_grid(plotlist = plots1)
ggsave("plot1.PNG", plot_out1, width = 18, height = 12, dpi = 450)

# 2. Simulated DC ####
K <- unique(res2$K.Var3)
N <- unique(res2$N.Var2)
sim_pars <- expand.grid(N, K)

plots2 <- apply(sim_pars[order(sim_pars$Var1),], 1, function(x) {
  
  match_conds <- (res2[,2] == x[1]) & (res2[,3] == x[2])
  res_subs <- res2[match_conds, c(1:(x[2]+3), 12)]
  
  qp <- seq(-6, 6, length.out = 121)
  
  base_plot <- ggplot(data.frame(X = qp), aes(x=X)) +
    stat_function(fun = ddc, args = list(mean = 0, sd = 1, phi = c(77.32, 78.51, 76.33, 77.16)),
                  colour = "red", size=1.1) +
    ggtitle(paste("N =", x[1], "|", "# DC Par =", x[2])) +
    theme_minimal()
  
  for (i in 1:nrow(res_subs)) {
    phi <- res_subs[i, 4:(x[2]+3)]
    base_plot <- base_plot + stat_function(fun = ddc, args = list(mean = 0, sd = 1, phi = unlist(phi)),
                                           alpha = .05)
  }
  
  base_plot
})

plot_out2 <- plot_grid(plotlist = plots2, ncol = 2)
ggsave("plot2.PNG", plot_out2, width = 9, height = 12, dpi = 450)

# 4.1 E-table not standardized ####

K <- unique(res4i$K.Var1)
sim_pars <- expand.grid(K)

qp <- as.vector(gTheta[[1]])

plot4i_e <- ggplot(data.frame(X = qp, Y = rr), aes(x=X, y=Y)) + geom_line() +
  ggtitle("E-table")

plot4i <- apply(sim_pars, 1, function(x) {
  match_conds <- (res4i[,1] == x)
  phi <- res4i[match_conds, 2:(x[1]+1)]
  
  plot <- ggplot(data.frame(X = qp), aes(x=X)) +
    stat_function(fun = ddc, args = list(mean = 0, sd = 1, phi = unlist(phi)), alpha = 1) +
    ylim(0, .5) +
    ggtitle(paste("# DC Par =", x[1])) +
    theme_minimal()

  
  plot
})

plot4i[[9]] <- plot4_e

plot_out4i <- plot_grid(plotlist = plot4i[c(9, 1:8)], ncol = 3)
ggsave("plot4i.PNG", plot_out4i, width = 9, height = 12, dpi = 450)

# 4.2 E-table standardized ####

K <- unique(res4ii$K.Var1)
sim_pars <- expand.grid(K)

obs_m <- sum(gTheta[[1]]) / sum(rr)
obs_sd <- sqrt(sum(rr * (gTheta[[1]] - obs_m)^2) / sum(rr))
gTheta_std <- (gTheta[[1]] - obs_m) / obs_sd
std_qp <- standardizeQuadrature(gTheta[[1]], gTheta_std, rr)
std_qp[,2] <- std_qp[,2] / sum(std_qp[,2]) * sum(rr)

plot4ii_e <- ggplot(data.frame(X = std_qp[,1], Y = std_qp[,2]), aes(x=X, y=Y)) + geom_line() +
  ggtitle("E-table")

plot4ii <- apply(sim_pars, 1, function(x) {
  match_conds <- (res4ii[,1] == x)
  phi <- res4ii[match_conds, 2:(x[1]+1)]
  
  plot <- ggplot(data.frame(X = std_qp[,1]), aes(x=X)) +
    stat_function(fun = ddc, args = list(mean = 0, sd = 1, phi = unlist(phi)), alpha = 1) +
    ylim(0, .5) +
    ggtitle(paste("# DC Par =", x[1])) +
    theme_minimal()
  
  
  plot
})

plot4ii[[9]] <- plot4ii_e

plot_out4ii <- plot_grid(plotlist = plot4ii[c(9, 1:8)], ncol = 3)
ggsave("plot4ii.PNG", plot_out4ii, width = 9, height = 12, dpi = 450)
