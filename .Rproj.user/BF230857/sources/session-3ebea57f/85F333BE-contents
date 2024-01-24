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

# all_res <- parallel::mclapply(
#   1:length(mcmcs),
#   function(i, mcmcs, datas){
#     local_mcmc(mcmcs[[i]], datas[[i]])
#   },
#   mcmcs = mcmcs,
#   datas = datas,
#   mc.cores = 12
# )

all_res <- list()
for (j in 1:length(mcmcs)) {
  all_res[[j]] <- list()
  all_res[[j]][[1]] <- mcmcs[[j]]
}




for (i in 1:length(mcmcs)) {

  if(all(
    mcmcs[[i]]$g_lik[2:mcmcs[[i]]$n] == sapply(2:mcmcs[[i]]$n, g_lik, mcmc = mcmcs[[i]], data = datas[[i]])
  )){

  }else{
    print(i)
  }
}


i = 145
which(mcmcs[[i]]$g_lik[2:mcmcs[[i]]$n] != sapply(2:mcmcs[[i]]$n, g_lik, mcmc = mcmcs[[i]], data = datas[[i]]))


datas[[i]]$frozen

mcmcs[[i]]$anc_cluster
-492.786228
mcmcs[[166]]$g_lik[datas[[166]]$frozen]

double_rnorm <- function(x){
  s <- rnorm(1)
  t <- rnorm(1)
  c(s,t)
}

mclapply(1:5, double_rnorm, mc.preschedule = F, mc.set.seed = F)
rnorm(1, rnorm(1), 1)

s <- .Random.seed
set.seed(2)
#rnorm(1)
double_rnorm(1)
set.seed(2)
mclapply(1:10, double_rnorm, mc.cores = 10, mc.set.seed = F, mc.preschedule = F)



