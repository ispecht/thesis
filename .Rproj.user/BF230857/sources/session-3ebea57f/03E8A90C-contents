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
  n_subtrees <- length(subtrees[[1]])
  mcmcs <- list()
  datas <- list()
  for (i in 1:n_subtrees) {
    mcmcs[[i]] <- mcmc
    mcmcs[[i]]$root <- subtrees[[1]][i]
    mcmcs[[i]]$cluster <- subtrees[[2]][[i]]
    if(mcmcs[[i]]$root == 1){
      mcmcs[[i]]$anc_cluster <- NA
    }else{
      for (j in 1:n_subtrees) {
        if(mcmcs[[i]]$root %in% c(subtrees[[1]][j], subtrees[[2]][[j]])){
          mcmcs[[i]]$anc_cluster <- j
        }
      }
    }


    joined <- c(mcmcs[[i]]$root, mcmcs[[i]]$cluster)

    # To save memory: extract only necessary components of MCMC and data
    mcmcs[[i]]$n <- length(joined)
    mcmcs[[i]]$h <- mcmcs[[i]]$h[joined]
    mcmcs[[i]]$h <- match(mcmcs[[i]]$h, joined)
    mcmcs[[i]]$w <- mcmcs[[i]]$w[joined]
    mcmcs[[i]]$t <- mcmcs[[i]]$t[joined]
    mcmcs[[i]]$m01 <- mcmcs[[i]]$m01[joined]
    mcmcs[[i]]$m10 <- mcmcs[[i]]$m10[joined]
    mcmcs[[i]]$m0y <- mcmcs[[i]]$m0y[joined]
    mcmcs[[i]]$m1y <- mcmcs[[i]]$m1y[joined]
    mcmcs[[i]]$mx0 <- mcmcs[[i]]$mx0[joined]
    mcmcs[[i]]$mx1 <- mcmcs[[i]]$mx1[joined]
    mcmcs[[i]]$mxy <- mcmcs[[i]]$mxy[joined]
    mcmcs[[i]]$d <- mcmcs[[i]]$d[joined]
    mcmcs[[i]]$g_lik <- mcmcs[[i]]$g_lik[joined]

    datas[[i]] <- data
    datas[[i]]$s <- datas[[i]]$s[joined]
    datas[[i]]$n_obs <- sum(joined <= data$n_obs)
    datas[[i]]$snvs <- datas[[i]]$snvs[joined]
    datas[[i]]$frozen <- setdiff(which(joined %in% subtrees[[1]]), 1)

    # Degree of roots of other trees needs to be set to 0
    mcmcs[[i]]$d[datas[[i]]$frozen] <- 0

    # e_lik needs to be re-computed based on the smaller set
    mcmcs[[i]]$e_lik <- e_lik(mcmcs[[i]], datas[[i]])

  }

  return(list(mcmcs, datas))

}
