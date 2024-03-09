### Sensitivity analysis
# Parameters to assess sensitivity: (as changes from default values)

source("main.R")

## Unhash these for regular run

init <- initialize(init_mst = T)
#init <- initialize()
mcmc <- init[[1]]
data <- init[[2]]

data$n_local = 10
data$sample_every = 10
data$n_global = 100000

## Unhash these for run on pre-computed data and initial MCMC

# load("data.RData")
# load("mcmc.RData")

a_gs <- c(4, 6) # Default = 5
lambda_gs <- a_gs/5 # Default = 1. Also update a_g to maintain mean of 5. Corresponds to variance of 5 (default), 25, 1
a_ss <- c(4, 6) # Default = 5
lambda_ss <- a_ss/5 # Default = 1. Also update a_s to maintain mean of 5. Corresponds to variance of 5 (default), 25, 1
rhos <- c(0.05, 0.2) # Default = 0.1. Also update psi = rho / (2.5 + rho)
psis <- 0.1 / (c(1.5, 3.5) + 0.1)
vs <- c(500, 2000)

# jth column is corresponds to each of the above params, in order
combos <- matrix(c(5, 1, 5, 1, 0.1, 0.1/(2.5 + 0.1), 1000), nrow = 15, ncol = 7, byrow = T)
combos[1:2, 1] <- a_gs
combos[3:4, 1] <- lambda_gs * 5
combos[3:4, 2] <- lambda_gs
combos[5:6, 3] <- a_ss
combos[7:8, 3] <- lambda_ss * 5
combos[7:8, 4] <- lambda_ss
combos[9:10, 5] <- rhos
combos[9:10, 6] <- rhos / (2.5 + rhos)
combos[11:12, 6] <- psis
combos[13:14, 7] <- vs


run_combo <- function(i, mcmc){
  mcmc$a_g <- combos[i, 1]
  mcmc$lambda_g <- combos[1, 2]
  mcmc$a_s <- combos[i, 3]
  mcmc$lambda_s <- combos[i, 4]
  mcmc$rho <- combos[i, 5]
  mcmc$psi <- combos[i, 6]
  mcmc$v <- combos[i, 7]

  out <- run_mcmc(mcmc, data)
  return(out)
}

outs <- mclapply(1:nrow(combos), run_combo, mcmc = mcmc, mc.cores = nrow(combos))

save(outs, file = "outs.RData")
#load("outs.RData")

sensitivity <- list()
all_liks <- list()
for (i in 1:nrow(combos)) {
  all_liks[[i]] <- outs[[i]][[1]]
  sensitivity[[i]] <- outs[[i]][[2]]

}

## Get adjacency matrix. Ancestor > data$n_obs doesn't count for anything. Indirect transmission OK
get_adj <- function(run, n_obs){
  out <- matrix(0, nrow = n_obs, ncol = n_obs)
  len <- length(run)
  for (i in (len*0.1 + 1):len) {
    h <- run[[i]]$h
    present <- which(!is.na(h) & h <= n_obs)
    present <- present[present <= n_obs]
    out[cbind(h[present], present)] <- out[cbind(h[present], present)] + 1
  }
  return(out / (0.9*len))
}

adjs <- lapply(sensitivity, get_adj, n_obs = data$n_obs)


# Run on default settings
default <- outs[[nrow(combos)]]
default_adj <- get_adj(default[[2]], data$n_obs)

# Histogram of change in adjacency matrix
network_change <- function(i){
  change <- as.vector(adjs[[i]] - default_adj)
  change <- change[abs(change) > 0]

  p <- ggplot(data.frame(x = change), aes(x = x)) +
    geom_histogram(aes(y=after_stat(density)), binwidth = 0.1, boundary = 0.1, color = "white", fill = "grey") +
    xlab("Change in Posterior Probability") +
    ylab("Probability Density") +
    #scale_y_continuous(trans='log2') +
    xlim(-1, 1) +
    theme_minimal()
  p
}

hists <- lapply(1:nrow(combos), network_change)
all_hists <- plot_grid(
  hists[[1]],
  hists[[2]],
  hists[[3]],
  hists[[4]],
  hists[[5]],
  hists[[6]],
  hists[[7]],
  hists[[8]],
  hists[[9]],
  hists[[10]],
  hists[[11]],
  hists[[12]],
  hists[[13]],
  hists[[14]],
  hists[[15]],
  ncol = 3,
  labels = "AUTO"
)

ggsave("./figs/hist.pdf", width = 7.8, height = 9)
ggsave("./figs/hist.png", width = 7.8, height = 9)

## Next up, density plots of mu
param_dens <- function(i, param){
  run <- sensitivity[[i]]
  out <- c()
  len <- length(run)
  for (j in (len*0.1 + 1):len) {
    out <- c(out, (run[[j]][[param]]))
  }
  p <- ggplot(data.frame(x = out), aes(x = x)) +
    geom_histogram(aes(y=after_stat(density)), color = "white", fill = "grey") +
    xlab(paste("Value of", param)) +
    ylab("Probability Density") +
    scale_x_continuous(breaks = signif(seq(min(out), max(out), by = (max(out) - min(out)) / 2) , 2)) +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    theme_minimal() +
    theme(plot.margin = margin(t = 0, r = 0.2, b = 0, l = 0, unit = "in"))
  p
}

output_plots <- function(param){
  dens <- lapply(1:nrow(combos), param_dens, param = param)
  all_dens <- plot_grid(
    dens[[1]],
    dens[[2]],
    dens[[3]],
    dens[[4]],
    dens[[5]],
    dens[[6]],
    dens[[7]],
    dens[[8]],
    dens[[9]],
    dens[[10]],
    dens[[11]],
    dens[[12]],
    dens[[13]],
    dens[[14]],
    dens[[15]],
    ncol = 3,
    labels = "AUTO"
  )

  ggsave(paste0("./figs/", param, ".pdf"), width = 7.8, height = 10.2)
  ggsave(paste0("./figs/", param, ".png"), width = 7.8, height = 10.2)

}

output_plots("mu")
output_plots("p")
#output_plots("lambda")
output_plots("b")




# cmd("gcloud storage cp ./figs/b.pdf gs://broad-viral/ispecht")
# cmd("gcloud storage cp ./figs/mu.pdf gs://broad-viral/ispecht")
# cmd("gcloud storage cp ./figs/p.pdf gs://broad-viral/ispecht")
# cmd("gcloud storage cp ./figs/hist.pdf gs://broad-viral/ispecht")
