# Update global parameters

global_mcmc <- function(mcmc, data){

  mcmc <- moves$b(mcmc, data)
  #mcmc <- moves$a_g(mcmc, data)
  #mcmc <- moves$a_s(mcmc, data)
  mcmc <- moves$mu(mcmc, data)
  mcmc <- moves$p(mcmc, data)
  mcmc <- moves$v(mcmc, data)
  #mcmc <- moves$rho(mcmc, data)
  #mcmc <- moves$psi(mcmc, data)


  return(mcmc)
}

## After making global moves, chop up the tree, and make one copy of "mcmc" for each subtree
breakdown <- function(mcmc, data){

  subtrees <- chop(mcmc, data)
  mcmcs <- list()
  for (i in 1:data$n_subtrees) {
    mcmcs[[i]] <- mcmc
    mcmcs[[i]]$root <- subtrees[[1]][i]
    mcmcs[[i]]$cluster <- subtrees[[2]][[i]]
  }

  return(mcmcs)

}
