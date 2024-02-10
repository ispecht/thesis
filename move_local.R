### Script to execute local mcmc on a subtree

source("likelihood.R")
source("moves.R")
source("prior.R")
source("subroutines.R")
source("local_mcmc.R")

set.seed(1)

i <- commandArgs(trailingOnly=TRUE)
load(paste0("./state/substate/tree_", i, "/mcmc.RData"))
load(paste0("./state/substate/tree_", i, "/data.RData"))

res <- local_mcmc(mcmc_tmp, data_tmp)

save(res, file = paste0("./state/substate/tree_", i, "/res.RData"))
