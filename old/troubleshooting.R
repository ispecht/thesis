## Try to diagnose failure of MCMC

load("ispecht_mcmcs.RData")
load("data.RData")


datas <- list()
roots <- c()
for (i in 1:length(mcmcs)) {
  roots[i] <- mcmcs[[i]]$root
}


for (i in 1:length(mcmcs)) {
  joined <- c(mcmcs[[i]]$root, mcmcs[[i]]$cluster)
  datas[[i]] <- data
  datas[[i]]$s <- datas[[i]]$s[joined]
  datas[[i]]$n_obs <- sum(joined <= data$n_obs)
  datas[[i]]$snvs <- datas[[i]]$snvs[joined]
  datas[[i]]$frozen <- setdiff(which(joined %in% roots), 1)
}

set.seed(r)

for (m in seq(25, 145, 12)) {
  set.seed(r)
  all_res <- parallel::mclapply(
    157:167,
    function(i, mcmcs, datas){
      local_mcmc(mcmcs[[i]], datas[[i]])
    },
    mcmcs = mcmcs,
    datas = datas,
    mc.set.seed = F,
    mc.cores = 11
  )
  print(m)
}




set.seed(r)
mclapply(1:11, function(x){rnorm(1)}, mc.cores = 11, mc.set.seed = F)
set.seed(r)
rnorm(1)

Sys.time()
sapply(2:mcmc$n, g_lik, mcmc = mcmc, data = data)
Sys.time()
unlist(parallel::mclapply(2:mcmc$n, g_lik, mcmc = mcmc, data = data, mc.cores = 12))
Sys.time()
