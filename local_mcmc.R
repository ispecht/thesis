# Schedule of moves within a subtree

local_mcmc <- function(mcmc, data){

  res <- list()
#
#    set.seed(212)
#    mcmc <- mcmcs[[j]]
#    data$n_local <- 60

  for (r in 1:data$n_local) {
    mcmc <- moves$w(mcmc, data)
    mcmc <- moves$t(mcmc, data)
    mcmc <- moves$swap(mcmc, data, exchange_children = T)
    mcmc <- moves$h_step(mcmc, data, resample_t = T)
    mcmc <- moves$genotype(mcmc, data)
    #
    mcmc <- moves$h_global(mcmc, data)
    mcmc <- moves$create(mcmc, data)

    mcmc <- moves$w_t(mcmc, data)
    mcmc <- moves$h_step(mcmc, data)
    mcmc <- moves$swap(mcmc, data)
    mcmc <- moves$h_step(mcmc, data, resample_t = T, resample_w = T)

    if(r %% 100 == 0){
      res <- c(res, list(mcmc))
    }

  }

  return(res)

}

## Join together results calculated in parallel across subtrees
amalgamate <- function(all_res, mcmcs, data){

  # Number of samples for each subtree
  n_samples <- length(all_res[[1]])

  # Number of subtrees
  n_subtrees <- data$n_subtrees


  # If we didn't break up the tree, nothing to do here!
  if(n_subtrees == 1){
    return(all_res[[1]])
  }else{

    ## Loop through each sample and produce an amalgamated MCMC state:

    # Create a list to store the amalgamated results
    res <- list()
    for (i in 1:n_samples) {

      # Get the ancestral cluster of each cluster
      anc_clusters <- c()
      roots <- c()
      for (j in 1:n_subtrees) {
        roots[j] <- all_res[[j]][[i]]$root
        anc <- all_res[[j]][[i]]$h[all_res[[j]][[i]]$root]

        if(is.na(anc)){
          anc_clusters[j] <- NA
        }else{
          for (k in 1:n_subtrees) {
            if(anc %in% all_res[[k]][[i]]$cluster | anc == all_res[[k]][[i]]$root){
              anc_clusters[j] <- k
            }
          }
        }
      }

      # First determine who the unobserved hosts are in each cluster, so that they may be re-indexed
      unobs <- list()
      cls <- list() # Cluster including roots of upstream clusters, as
      for (j in 1:n_subtrees) {
        cls[[j]] <- c(roots[which(anc_clusters == j)], all_res[[j]][[i]]$cluster)

        unobs[[j]] <- cls[[j]][cls[[j]] > data$n_obs]
      }

      # Then create a mapping to the indexing of mcmc
      mappings <- list()
      displacement <- 0
      for (j in 1:n_subtrees) {
        # mappings[[j]] = new names in output mcmc
        # cls[[j]] = current names in all_res[[j]][[i]]

        mappings[[j]] <- cls[[j]]
        if(length(unobs[[j]]) > 0){
          new_unobs_index <- (data$n_obs + displacement + 1):(data$n_obs + displacement + length(unobs[[j]]))
          mappings[[j]][match(unobs[[j]], mappings[[j]])] <- new_unobs_index
        }

        displacement <- displacement + length(unobs[[j]])
      }

      ## Now, fill an amalgamated mcmc with info from the correct subtrees

      # Initialization doesn't really matter; we will initialize to the mcmc state of the first subtree
      mcmc <- all_res[[1]][[i]]
      # Correct lengths of entries, to avoid carrying over extraneous information
      mcmc$n <- data$n_obs + displacement
      mcmc$h <- mcmc$h[1:mcmc$n]
      mcmc$w <- mcmc$w[1:mcmc$n]
      mcmc$t <- mcmc$t[1:mcmc$n]
      mcmc$m01 <- mcmc$m01[1:mcmc$n]
      mcmc$m10 <- mcmc$m10[1:mcmc$n]
      mcmc$m0y <- mcmc$m0y[1:mcmc$n]
      mcmc$m1y <- mcmc$m1y[1:mcmc$n]
      mcmc$mx0 <- mcmc$mx0[1:mcmc$n]
      mcmc$mx1 <- mcmc$mx1[1:mcmc$n]
      mcmc$mxy <- mcmc$mxy[1:mcmc$n]
      mcmc$d <- mcmc$d[1:mcmc$n]
      mcmc$g_lik <- mcmc$g_lik[1:mcmc$n]
      mcmc$root <- NULL
      mcmc$cluster <- NULL


      # Update all entries of mcmc for each cluster (plus roots of upstream clusters)
      for (j in 1:n_subtrees) {

        # Ancestors to fill in
        ancs <- all_res[[j]][[i]]$h[cls[[j]]]

        # Which of these ancestors are in cls[[j]]?
        in_cl <- ancs %in% cls[[j]]

        # Update index for ancestors in cluster
        ancs[in_cl] <- mappings[[j]][match(ancs[in_cl], cls[[j]])]

        # If not in cluster, use the mapping of the ancestral cluster (if there is one)
        if(!is.na(anc_clusters[j])){
          ancs[!in_cl] <- mappings[[anc_clusters[[j]]]][match(ancs[!in_cl], cls[[anc_clusters[[j]]]])]
        }
        mcmc$h[mappings[[j]]] <- ancs

        # Update other components of mcmc (easier!)


        mcmc$w[mappings[[j]]] <- all_res[[j]][[i]]$w[cls[[j]]]
        mcmc$t[mappings[[j]]] <- all_res[[j]][[i]]$t[cls[[j]]]
        mcmc$m01[mappings[[j]]] <- all_res[[j]][[i]]$m01[cls[[j]]]
        mcmc$m10[mappings[[j]]] <- all_res[[j]][[i]]$m10[cls[[j]]]
        mcmc$m0y[mappings[[j]]] <- all_res[[j]][[i]]$m0y[cls[[j]]]
        mcmc$m1y[mappings[[j]]] <- all_res[[j]][[i]]$m1y[cls[[j]]]
        mcmc$mx0[mappings[[j]]] <- all_res[[j]][[i]]$mx0[cls[[j]]]
        mcmc$mx1[mappings[[j]]] <- all_res[[j]][[i]]$mx1[cls[[j]]]
        mcmc$mxy[mappings[[j]]] <- all_res[[j]][[i]]$mxy[cls[[j]]]
        mcmc$g_lik[mappings[[j]]] <- all_res[[j]][[i]]$g_lik[cls[[j]]]

        # For e_lik, compute as differences
        if(i == 1){
          if(j == 1){
            mcmc$e_lik <- all_res[[j]][[i]]$e_lik
          }else{
            mcmc$e_lik <- mcmc$e_lik + all_res[[j]][[i]]$e_lik - mcmcs[[j]]$e_lik
          }
        }else{
          if(j == 1){
            mcmc$e_lik <- res[[i-1]]$e_lik + all_res[[j]][[i]]$e_lik - all_res[[j]][[i-1]]$e_lik
          }else{
            mcmc$e_lik <- mcmc$e_lik + all_res[[j]][[i]]$e_lik - all_res[[j]][[i-1]]$e_lik
          }
        }
      }


      ## Correct node degrees
      # Test this; check not much changes
      mcmc$d <- sapply(1:mcmc$n, function(x){sum(mcmc$h[2:mcmc$n] == x)}) # Node degrees

      res[[i]] <- mcmc
    }

    return(res)

  }
}






