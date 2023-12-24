### MCMC moves

moves <- list()

## Update one of the w_i's by adding or subtracting rounded N(0,3)
moves$w <- function(mcmc, data){
  # Choose random host with ancestor
  i <- sample(2:mcmc$n, 1)

  # Proposal
  prop <- mcmc
  change <- round(rnorm(1, 0, 3))
  prop$w[i] <- mcmc$w[i] + change

  prop$e_lik <- e_lik(prop, data)
  prop$g_lik[i] <- g_lik(prop, data, i)
  prop$prior <- prior(prop)


  # Here the Hastings ratio is not 1, because we're imagining that the newly-added / deleted intermediate hosts are labeled.
  if(change > 0){
    # If added new hosts, P(new to old) is probability of correctly selecting the new hosts to delete
    # P(old to new) is probability choosing the correct new hosts to add from larger population
    hastings <-  -lchoose(prop$w[i], change) -
      (-lchoose(data$N, change))
  }else if(change < 0){
    hastings <- -lchoose(data$N, -change) -
      (-lchoose(mcmc$w[i], -change))
  }else{
    hastings <- 0
  }

  if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior + hastings){
    return(prop)
  }else{
    return(mcmc)
  }

}

## Update one of the t_i's using a N(0,1) proposal density if observed; N(0, 10) if not
moves$t <- function(mcmc, data){
  # Choose random host with ancestor
  i <- sample(2:mcmc$n, 1)
  # Proposal
  prop <- mcmc
  prop$t[i] <- rnorm(1, mcmc$t[i], ifelse(i <= data$n_obs, 1, 10))
  prop$e_lik <- e_lik(prop, data)
  update <- c(i, which(mcmc$h == i)) # For which hosts must we update the genomic likelihood?
  prop$g_lik[update] <- sapply(update, g_lik, mcmc = prop, data = data)
  prop$prior <- prior(prop)

  if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior){
    return(prop)
  }else{
    return(mcmc)
  }
}

## Update t_i and w_i and w_j's simultaneously, where j is the child of i
moves$w_t <- function(mcmc, data){
  # Choose random host with ancestor and children
  choices <- intersect(2:mcmc$n, which(mcmc$d > 0))
  if(length(choices) == 0){
    return(mcmc)
  }else{
    i <- ifelse(length(choices) == 1, choices, sample(choices, 1))
    js <- which(mcmc$h == i)
    # If multiple children, randomize the order
    if(length(js) > 1){
      js <- sample(js, length(js), replace = F)
    }
    h <- mcmc$h[i]

    # How big could the weight of i possibly be?
    max_dist <- mcmc$w[i] + min(mcmc$w[js])

    prop <- mcmc
    prop$w[i] <- sample(0:max_dist, 1)
    delta <- prop$w[i] - mcmc$w[i] # Change in weight
    prop$w[js] <- mcmc$w[js] - delta

    # Now sample the time i contracts the disease as a beta centered at its new position
    max_t <- min(mcmc$t[js])
    prop$t[i] <- mcmc$t[h] + (max_t - mcmc$t[h]) * rbeta(1, prop$w[i] + 1, max_dist - prop$w[i] + 1)

    ## Hastings ratio
    hastings <- 0
    if(length(js) > 1){
      if(delta > 0){
        # In this case, we lose weight from the edges leading into the js, so:
        hastings <- hastings - lchoose(data$N, delta) * (length(js) - 1) + # P(new -> old): add the correct nodes from the population onto each of the #[js] - 1 edges
          sum(lchoose(mcmc$w[js[2:length(js)]], delta)) # P(old -> new): delete the correct nodes from js[2], js[3], ...

      }
      if(delta < 0){
        hastings <- hastings - sum(lchoose(prop$w[js[2:length(js)]], -delta)) + # P(new -> old): delete the correct nodes from js[2], js[3], ...
          lchoose(data$N, -delta) * (length(js) - 1) # P(old -> new): add the correct nodes from the population to each of js[2], js[3], ...
      }
    }

    # Proposal density for the beta draw
    hastings <- hastings + dbeta((mcmc$t[i] - mcmc$t[h]) / (max_t - mcmc$t[h]), mcmc$w[i] + 1, max_dist - mcmc$w[i] + 1, log = T) - # P(new -> old)
      dbeta((prop$t[i] - prop$t[h]) / (max_t - prop$t[h]), prop$w[i] + 1, max_dist - prop$w[i] + 1, log = T) # P(old -> new)

    prop$e_lik <- e_lik(prop, data)
    update <- c(i, js) # For which hosts must we update the genomic likelihood?
    prop$g_lik[update] <- sapply(update, g_lik, mcmc = prop, data = data)
    prop$prior <- prior(prop)

    if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior + hastings){
      return(prop)
    }else{
      return(mcmc)
    }
  }
}

## Update b using a N(0,0.01) proposal density
moves$b <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$b <- rnorm(1, mcmc$b, 0.01)
  prop$e_lik <- e_lik(prop, data)
  prop$g_lik[2:mcmc$n] <- sapply(2:mcmc$n, g_lik, mcmc = prop, data = data)
  prop$prior <- prior(prop)

  if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior){
    return(prop)
  }else{
    return(mcmc)
  }
}

## Update a_g using a N(0,1) proposal density
moves$a_g <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$a_g <- rnorm(1, mcmc$a_g, 1)
  prop$e_lik <- e_lik(prop, data)
  prop$prior <- prior(prop)

  if(log(runif(1)) < prop$e_lik + prop$prior - mcmc$e_lik - mcmc$prior){
    return(prop)
  }else{
    return(mcmc)
  }
}

## Update a_s using a N(0,1) proposal density
moves$a_s <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$a_s <- rnorm(1, mcmc$a_s, 1)
  prop$e_lik <- e_lik(prop, data)
  prop$prior <- prior(prop)

  if(log(runif(1)) < prop$e_lik + prop$prior - mcmc$e_lik - mcmc$prior){
    return(prop)
  }else{
    return(mcmc)
  }
}

## Update mu using a N(0,1e-7) proposal density
moves$mu <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$mu <- rnorm(1, mcmc$mu, 1e-7)
  prop$e_lik <- e_lik(prop, data)
  prop$g_lik[2:mcmc$n] <- sapply(2:mcmc$n, g_lik, mcmc = prop, data = data)
  prop$prior <- prior(prop)

  if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior){
    return(prop)
  }else{
    return(mcmc)
  }
}

## Update p using a N(0,1e-7) proposal density
moves$p <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$p <- rnorm(1, mcmc$p, 1e-7)
  prop$e_lik <- e_lik(prop, data)
  prop$g_lik[2:mcmc$n] <- sapply(2:mcmc$n, g_lik, mcmc = prop, data = data)
  prop$prior <- prior(prop)

  if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior){
    return(prop)
  }else{
    return(mcmc)
  }
}

## Update v using a N(0,100) proposal density (rounded)
moves$v <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$v <- round(rnorm(1, mcmc$v, 100))
  prop$e_lik <- e_lik(prop, data)
  prop$g_lik[2:mcmc$n] <- sapply(2:mcmc$n, g_lik, mcmc = prop, data = data)
  prop$prior <- prior(prop)

  if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior){
    return(prop)
  }else{
    return(mcmc)
  }
}

## Update rho using a N(0,0.1) proposal density
moves$rho <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$rho <- rnorm(1, mcmc$rho, 0.1)
  prop$e_lik <- e_lik(prop, data)
  prop$prior <- prior(prop)

  if(log(runif(1)) < prop$e_lik + prop$prior - mcmc$e_lik - mcmc$prior){
    return(prop)
  }else{
    return(mcmc)
  }
}

## Update psi using a N(0,0.1) proposal density
moves$psi <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$psi <- rnorm(1, mcmc$psi, 0.1)
  prop$e_lik <- e_lik(prop, data)
  prop$prior <- prior(prop)

  if(log(runif(1)) < prop$e_lik + prop$prior - mcmc$e_lik - mcmc$prior){
    return(prop)
  }else{
    return(mcmc)
  }
}

## Update genotype at (a) missing sites in observed host, or (b) all sites in unobserved host

moves$genotype <- function(mcmc, data){

  # Choose random host with ancestor
  i <- sample(2:mcmc$n, 1)
  js <- which(mcmc$h == i) # Children
  # Let h denote the ancestor of i; never used in computations

  # Get a list of "SNVs of interest": sites that change going from h to i, or i to a child of i
  interest <- unique(c(
    mcmc$m01[[i]],
    mcmc$mx1[[i]],
    unlist(mcmc$m10[js]),
    unlist(mcmc$m1y[js]),
    mcmc$m10[[i]],
    mcmc$mx0[[i]],
    unlist(mcmc$m01[js]),
    unlist(mcmc$m0y[js])
  ))

  # If i is observed, we can only change sites with missing data
  if(i <= data$n_obs){
    interest <- intersect(interest, data$snvs[[i]]$missing$call)
  }

  if(length(interest) == 0){
    return(mcmc)
  }else{

    # Proposal
    prop <- mcmc

    # Pick one SNV to update. We switch whether it exists or not in i.
    snv <- ifelse(length(interest) == 1, interest, sample(interest, 1))

    prop <- flip_genotype(prop, mcmc, i, js, snv)

    prop$e_lik <- e_lik(prop, data)
    update <- c(i, js) # For which hosts must we update the genomic likelihood?
    prop$g_lik[update] <- sapply(update, g_lik, mcmc = prop, data = data)
    prop$prior <- prior(prop)

    if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior){
      #print(i)
      #print(snv)
      return(prop)
    }else{
      return(mcmc)
    }

  }
}

### Topological moves

## Move the ancestor of a node one step upstream (towards tips) or one step downstream (towards root) onto next/previous tracked host
moves$h_step <- function(mcmc, data, resample_t = FALSE){
  # Choose random host with ancestor
  i <- sample(2:mcmc$n, 1)
  h_old <- mcmc$h[i]

  # Proposal
  prop <- mcmc

  # Are we going upstream or downstream?
  upstream <- runif(1) < 1/2

  if(upstream){
    # Who are the other children of h_old?
    children <- setdiff(which(mcmc$h == h_old), i)

    # Which ones have a lesser edge weight?
    children <- children[mcmc$w[children] < mcmc$w[i]]

    # If no valid children, reject
    # Also reject if h_old is not observed and has <= 2 total children, because then we can't remove one
    if(length(children) == 0 | (h_old > data$n_obs & length(which(mcmc$h == h_old)) <= 2)){
      return(mcmc)
    }else{

      # Pick one
      h_new <- ifelse(length(children) == 1, children, sample(children, 1))

      prop <- shift_upstream(prop, i, h_old, h_new)

      # What's the change in edge weight for i? (Negative)
      change <- mcmc$w[h_new] + 1

      # Should we change the time of infection for i?
      if(resample_t){
        js <- which(prop$h == i)
        if(length(js) > 0){
          max_t <- min(prop$t[js])
        }else{
          max_t <- Inf
        }
        if(i <= data$n_obs){
          max_t <- min(max_t, data$s[i])
        }
        if(prop$t[h_new] > max_t){
          return(mcmc)
        }else{
          prop$t[i] <- runif(1, prop$t[h_new], max_t)
        }
      }


      prop$e_lik <- e_lik(prop, data)
      prop$g_lik[i] <- g_lik(prop, data, i)
      if(resample_t){
        if(length(js) > 0){
          prop$g_lik[js] <- sapply(js, g_lik, mcmc = prop, data = data)
        }
      }
      prop$prior <- prior(prop)

      hastings <- -lchoose(data$N, change) - # P(new -> old): choose [change] people to add back onto the edge
        (-lchoose(mcmc$w[i], change)) + # P(old -> new): choose [change] people to delete from the edge
        log(length(children)) # P(new -> old): 1; P(old -> new): choose from among #[children] people to be h_new

      if(resample_t){
        hastings <- hastings - log(max_t - prop$t[h_old]) + log(max_t - prop$t[h_new])
      }

      # Accept / reject
      if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior + hastings){
        return(prop)
      }else{
        return(mcmc)
      }

    }
  }else{
    if(h_old == 1 | (h_old > data$n_obs & length(which(mcmc$h == h_old)) <= 2)){
      # If no downstream move, reject
      # Also reject if h_old is not observed and has <= 2 total children, because then we can't remove one
      return(mcmc)
    }else{

      # New ancestor of i is ancestor's ancestor
      h_new <- mcmc$h[h_old]

      prop <- shift_downstream(prop, i, h_old, h_new)

      # What's the change in edge weight for i? (Positive)
      change <- mcmc$w[h_old] + 1

      # Should we change the time of infection for i?
      if(resample_t){
        js <- which(prop$h == i)
        if(length(js) > 0){
          max_t <- min(prop$t[js])
        }else{
          max_t <- Inf
        }
        if(i <= data$n_obs){
          max_t <- min(max_t, data$s[i])
        }
        if(prop$t[h_new] > max_t){
          return(mcmc)
        }else{
          prop$t[i] <- runif(1, prop$t[h_new], max_t)
        }
      }

      prop$e_lik <- e_lik(prop, data)
      prop$g_lik[i] <- g_lik(prop, data, i)
      if(resample_t){
        if(length(js) > 0){
          prop$g_lik[js] <- sapply(js, g_lik, mcmc = prop, data = data)
        }
      }
      prop$prior <- prior(prop)

      ## Compute the number of possible children who could be chosen by i in the new config

      # Who are the other children of h_old?
      children <- setdiff(which(prop$h == h_new), i)

      # Which ones have a lesser edge weight?
      children <- children[prop$w[children] < prop$w[i]]

      hastings <- -lchoose(prop$w[i], change) - # P(new -> old): choose [change] people to delete from the edge
        (-lchoose(data$N, change)) - # P(old -> new): choose [change] people to add back onto the edge
        log(length(children)) # P(new -> old): choose from among #[children] people to be h_new; P(old -> new): 1

      if(resample_t){
        hastings <- hastings - log(max_t - prop$t[h_old]) + log(max_t - prop$t[h_new])
      }

      # Accept / reject
      if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior + hastings){
        return(prop)
      }else{
        return(mcmc)
      }
    }
  }
}

## The swap
## Switch h -> i -> j to
## h -> j -> i
moves$swap <- function(mcmc, data){
  # Choose host with a child and a parent
  # If unobserved, i must have at least 3 children, because losing one
  choices <- union(which(mcmc$d >= 3), which(mcmc$d[1:data$n_obs] >= 1))
  choices <- setdiff(choices, 1)

  if(length(choices) == 0){
    return(mcmc)
  }else{
    # Pick i
    i <- ifelse(length(choices) == 1, choices, sample(choices, 1))

    # Pick h
    h <- mcmc$h[i]

    # Pick j
    children_i <- which(mcmc$h == i)
    j <- ifelse(length(children_i) == 1, children_i, sample(children_i, 1))

    children_i <- setdiff(children_i, j)
    children_j <- which(mcmc$h == j)

    # Update the state
    prop <- mcmc
    prop <- shift_downstream(prop, j, i, h) # Shift j from i onto h
    prop <- shift_upstream(prop, i, h, j) # Shift i from h onto j
    prop$w[j] <- mcmc$w[i] # Swapping edge weights
    prop$t[j] <- mcmc$t[i] # Swapping time of infection
    prop$w[i] <- mcmc$w[j]
    prop$t[i] <- mcmc$t[j]


    prop$e_lik <- e_lik(prop, data)
    update <- c(i, j, children_i, children_j) # For which hosts must we update the genomic likelihood?
    prop$g_lik[update] <- sapply(update, g_lik, mcmc = prop, data = data)
    prop$prior <- prior(prop)

    # Accept / reject
    if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior){
      print("baste")
      return(prop)
    }else{
      return(mcmc)
    }
  }
}


## Create / remove a node
moves$create <- function(mcmc, data){
  # Are we creating or deleting an unobserved node?
  if(runif(1) < 1/2){
    create <- T
  }else{
    create <- F
  }

  # Are we moving nodes onto the new node upstream or downstream?
  if(runif(1) < 1/2){
    upstream <- T
  }else{
    upstream <- F
  }

  if(create){
    # Pick any node with an ancestor. (Note, some choices impossible, but this is okay!)
    j1 <- sample(2:mcmc$n, 1)
    h <- mcmc$h[j1]

    if(mcmc$w[j1] == 0 | (upstream & mcmc$d[h] == 1) | (!upstream & mcmc$d[j1] == 0)){
      return(mcmc)
    }else{

      # Who else are we attaching to i?
      if(upstream){
        kids <- setdiff(which(mcmc$h == h), j1)
      }else{
        kids <- which(mcmc$h == j1)
      }

      j2s <- kids[runif(length(kids)) < data$p_move]

      # If moving nobody, or moving everyone upstream off an unobserved node, or moving all but 0 or 1 downstream off an unobserved node, reject
      if(length(j2s) == 0 | (upstream & h > data$n_obs & length(j2s) == length(kids)) | (!upstream & h > data$n_obs & length(j2s) >= length(kids) - 1)){
        return(mcmc)
      }else{
        js <- c(j1, j2s)

        # How far upstream from h is the new node?
        # Maximum is min(w[js]) - 1 to preserve sum of all edge weights
        if(upstream){
          max_dist <- min(mcmc$w[js]) - 1
        }else{
          max_dist <- mcmc$w[j1] - 1
        }


        # If min is 0, can't make the move
        if(max_dist < 0){
          return(mcmc)
        }else{

          # New edge weight coming into i, the new host
          dist <- sample(0:max_dist, 1)
          i <- mcmc$n + 1

          # Maximum time i could be infected
          max_t <- min(mcmc$t[js])

          # Proposal
          prop <- mcmc

          ## Stick i onto h
          prop$n <- mcmc$n + 1
          prop$h[i] <- h
          prop$w[i] <- dist
          prop$t[i] <- mcmc$t[h] + (max_t - mcmc$t[h]) * rbeta(1, dist + 1, max_dist - dist + 1) # Weighted average
          prop$d[i] <- 0
          prop$d[h] <- mcmc$d[h] + 1

          ## Initialize genotype for i. This is all changing, so we initialize as i loses all iSNVs to 0. Everything else stays the same
          prop$mx0[[i]] <- unique(c(
            mcmc$mx0[[j1]], mcmc$mxy[[j1]], mcmc$mx1[[j1]]
          ))
          prop$m01[[i]] <- character(0)
          prop$m10[[i]] <- character(0)
          prop$m0y[[i]] <- character(0)
          prop$m1y[[i]] <- character(0)
          prop$mx1[[i]] <- character(0)
          prop$mxy[[i]] <- character(0)

          ## Move all js onto i
          if(upstream){
            for (j in js) {
              prop <- shift_upstream(prop, j, h, i)
            }
          }else{
            prop <- shift_upstream(prop, j1, h, i)
            for (j2 in j2s) {
              prop <- shift_downstream(prop, j2, j1, i)
            }
          }

          ## Create new genotype for i
          geno <- genotype(prop, i, js, data$eps)
          prop <- geno[[1]]
          log_p <- geno[[2]]

          ## Hastings ratio
          hastings <- 0

          # Delta is change in number of hosts along edges into j2s
          if(upstream){
            delta <- dist + 1
            # In this case, we lose weight from the edges leading into the j2s, so:
            hastings <- hastings - lchoose(data$N, delta) * length(j2s) + # P(new -> old): add the correct nodes from the population onto each of the #[js] - 1 edges
              sum(lchoose(mcmc$w[j2s], delta)) # P(old -> new): delete the correct nodes from js[2], js[3], ...
          }else{
            delta <- prop$w[j1] + 1
            hastings <- hastings - sum(lchoose(prop$w[j2s], delta)) + # P(new -> old): delete the correct nodes from js[2], js[3], ...
              lchoose(data$N, delta) * length(j2s) # P(old -> new): add the correct nodes from the population to each of js[2], js[3], ...
          }

          hastings <- hastings -
            length(j2s) * log(data$p_move) - (length(kids) - length(j2s)) * log(1 - data$p_move) + # P(old -> new): Probability of choosing kids to move onto i
            log(max_dist + 1) - # P(old -> new): choose how far upstream
            dbeta((prop$t[i] - mcmc$t[h]) / (max_t - mcmc$t[h]), dist + 1, max_dist - dist + 1) - # P(old -> new): beta density for t_i
            log_p - # P(old -> new): probability of newly-created genotype for i
            log(sum(prop$h > data$n_obs, na.rm = T)) + # P(new -> old): pick a host with an unobserved ancestor
            log(mcmc$n - 1) # P(old -> new): pick a host with an ancestor

          update <- c(i, js)
          return(accept_or_reject(prop, mcmc, data, update, hastings))

        }
      }
    }
  }else{
    ## Delete a node by tucking it back inside its parent / child
    choices <- which(mcmc$h > data$n_obs)
    if(length(choices) == 0){
      return(mcmc)
    }else{
      j1 <- ifelse(length(choices) == 1, choices, sample(choices, 1))
      i <- mcmc$h[j1]
      h <- mcmc$h[i]
      js <- which(mcmc$h == i)
      j2s <- setdiff(js, j1)

      if(!upstream & any(mcmc$w[j2s] <= mcmc$w[j1])){
        return(mcmc)
      }else{

        prop <- mcmc

        # Put all children of i onto h or j1
        if(upstream){
          for (j in js) {
            prop <- shift_downstream(prop, j, i, h)
          }
        }else{
          for (j2 in j2s) {
            prop <- shift_upstream(prop, j2, i, j1)
          }
          # Put j1 onto h
          prop <- shift_downstream(prop, j1, i, h)
        }
        # Degree of h still needs to decrease by 1, because of deletion of i
        prop$d[h] <- prop$d[h] - 1

        ## For calculation of Hastings ratio:

        # How far upstream from h is the new node?
        # Maximum is min(w[js]) - 1 to preserve sum of all edge weights
        if(upstream){
          max_dist <- min(prop$w[js]) - 1
        }else{
          max_dist <- prop$w[j1] - 1
        }

        # Maximum time i could be infected
        max_t <- min(prop$t[js])

        # Who else are we attaching to i?
        if(upstream){
          kids <- setdiff(which(prop$h == h), j1)
        }else{
          kids <- which(prop$h == j1)
        }

        # Probability that genotype() returns the genotype of i in mcmc
        log_p <- genotype(mcmc, i, js, data$eps, comparison = T)



        ## Hastings ratio
        hastings <- 0

        # Delta is change in number of hosts along edges into j2s
        if(upstream){
          delta <- mcmc$w[i] + 1
          # In this case, we lose weight from the edges leading into the j2s, so:
          hastings <- hastings - sum(lchoose(prop$w[j2s], delta)) + # P(new -> old): delete the correct nodes from js[2], js[3], ...
            lchoose(data$N, delta) * length(j2s) # P(old -> new): add the correct nodes from the population onto each of the #[js] - 1 edges

        }else{
          delta <- mcmc$w[j1] + 1
          hastings <- hastings - lchoose(data$N, delta) * length(j2s) + # P(new -> old): add the correct nodes from the population to each of js[2], js[3], ...
            sum(lchoose(mcmc$w[j2s], delta)) # P(old -> new): delete the correct nodes from js[2], js[3], ...
        }

        hastings <- hastings +
          length(j2s) * log(data$p_move) + (length(kids) - length(j2s)) * log(1 - data$p_move) - # P(new -> old): Probability of choosing kids to move onto i
          log(max_dist + 1) + # P(new -> old): choose how far upstream
          dbeta((mcmc$t[i] - mcmc$t[h]) / (max_t - mcmc$t[h]), mcmc$w[i] + 1, max_dist - mcmc$w[i] + 1) + # P(new -> old): beta density for t_i
          log_p - # P(new -> old): probability newly-created genotype for i equals genotype for i in "mcmc"
          log(prop$n - 2) + # P(new -> old): pick a host with an ancestor (-2 because haven't yet updated n)
          log(sum(mcmc$h > data$n_obs, na.rm = T)) # P(old -> new): pick a host with an unobserved ancestor

        ## Re-indexing: everyone above i steps down 1, i gets deleted
        prop$h[which(prop$h > i)] <- prop$h[which(prop$h > i)] - 1
        prop$n <- prop$n - 1
        prop$h <- prop$h[-i]
        prop$w <- prop$w[-i]
        prop$t <- prop$t[-i]
        prop$m01 <- prop$m01[-i]
        prop$m10 <- prop$m10[-i]
        prop$m0y <- prop$m0y[-i]
        prop$m1y <- prop$m1y[-i]
        prop$mx0 <- prop$mx0[-i]
        prop$mx1 <- prop$mx1[-i]
        prop$mxy <- prop$mxy[-i]
        prop$d <- prop$d[-i]
        prop$g_lik <- prop$g_lik[-i]

        # Update the js
        js[js > i] <- js[js > i] - 1

        update <- js
        return(accept_or_reject(prop, mcmc, data, update, hastings))

      }
    }
  }
}



## Create a new unobserved node one step upstream
moves$create_upstream <- function(mcmc, data){

  # Are we creating or deleting an unobserved node?
  if(runif(1) < 1/2){
    create <- T
  }else{
    create <- F
  }

  if(create){
    # Select the ancestor of our newly created node.
    # Must be observed and have >= 2 children, or unobserved and >= 3 children
    choices <- union(
      which(mcmc$d[1:data$n_obs] >= 2),
      which(mcmc$d >= 3)
    )

    if(length(choices) == 0){
      return(mcmc)
    }else{
      h <- ifelse(length(choices) == 1, choices, sample(choices, 1))

      # Who will be the children of the new node? Assume it's always 2 children.
      js <- which(mcmc$h == h)
      if(length(js) > 2){
        js <- sample(js, 2, replace = F)
      }

      j1 <- js[1]
      j2 <- js[2]

      # How far upstream is the new node?
      # Maximum is min(w[j1], w[j2]) - 1 to preserve sum of all edge weights
      max_dist <- min(mcmc$w[j1], mcmc$w[j2]) - 1

      # If both are 0, can't make the move
      if(max_dist < 0){
        return(mcmc)
      }else{

        # New edge weight coming into i, the new host
        dist <- sample(0:max_dist, 1)
        i <- mcmc$n + 1

        # Maximum time i could be infected
        max_t <- min(mcmc$t[j1], mcmc$t[j2])

        # Proposal
        prop <- mcmc

        ## Update all entries in the Markov chain
        prop$n <- mcmc$n + 1
        prop$h[js] <- mcmc$n + 1 # Index of new host is n + 1
        prop$h[i] <- h
        prop$w[js] <- mcmc$w[js] - dist - 1
        prop$w[i] <- dist
        prop$t[i] <- mcmc$t[h] + (max_t - mcmc$t[h]) * rbeta(1, dist + 1, max_dist - dist + 1) # Weighted average
        prop$d[h] <- mcmc$d[h] - 1
        prop$d[i] <- 2

        # Create a genotype for the new host
        geno <- create_genotype(prop, mcmc, i, j1, j2)
        prop <- geno[[1]]
        ambiguous <- geno[[2]]

        ## Randomly swap the genotype at certain sites, to ensure reversibility
        # Probability of a swap is eps

        # Get a list of "SNVs of interest": sites that change going from h to i, or i to a child of i
        interest <- unique(c(
          prop$m01[[i]],
          prop$mx1[[i]],
          unlist(prop$m10[js]),
          unlist(prop$m1y[js]),
          prop$m10[[i]],
          prop$mx0[[i]],
          unlist(prop$m01[js]),
          unlist(prop$m0y[js])
        ))

        interest <- setdiff(interest, ambiguous) # No need to update ambiguous sites, because they're already 50/50

        swaps <- interest[runif(length(interest)) < data$eps]

        # Perform the flip
        for (snv in swaps) {
          prop <- flip_genotype(prop, prop, i, js, snv)
        }

        ## Okay, now we have the heinous task of computing the Hastings ratio.

        # How many nodes along edges were removed?
        change <- dist + 1 # (negative)

        # Assume, technically, that j2 moves onto the edge from h to j1

        hastings <- -lchoose(data$N, change) - # P(new -> old): choose [change] people to add back onto the edge
          (-lchoose(mcmc$w[j2], change)) - # P(old -> new): choose [change] people to delete from the edge
          log(sum(prop$d[(data$n_obs + 1):prop$n] == 2)) - log(2) + # P(new -> old): choose an unobserved host with exactly 2 children, then order them
          log(length(choices)) + lchoose(sum(mcmc$h == h, na.rm = T), 2) + log(2) + # P(old -> new): choose from among #[choices] people to be h, and then pick j1, j2
          log(max_dist + 1) - # P(old -> new): uniform density for max_dist
          dbeta((prop$t[i] - mcmc$t[h]) / (max_t - mcmc$t[h]), dist + 1, max_dist - dist + 1) + # P(old -> new): beta density for t_i
          length(ambiguous) * log(2) - # P(old -> new): make a choice of the genotype at each ambiguous site
          length(swaps) * log(data$eps) - # P(old -> new): each random swap has probability eps
          (length(interest) - length(swaps)) * log(1 - data$eps) # P(old -> new): each lack of random swap has probability 1 - eps

        prop$e_lik <- e_lik(prop, data)
        update <- c(i, js) # For which hosts must we update the genomic likelihood?
        prop$g_lik[update] <- sapply(update, g_lik, mcmc = prop, data = data)
        prop$prior <- prior(prop)

        # Accept / reject
        if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior + hastings){
          return(prop)
        }else{
          return(mcmc)
        }
      }
    }
  }else{
    ## Getting rid of an unobserved node

    # Choices: any unobserved node with degree exactly 2
    choices <- which(mcmc$d == 2)
    choices <- choices[choices > data$n_obs]

    if(length(choices) == 0){
      return(mcmc)
    }else{
      # Select one to delete
      i <- ifelse(length(choices) == 1, choices, sample(choices, 1))

      # Ancestor
      h <- mcmc$h[i]

      # Children
      js <- sample(which(mcmc$h == i), 2) # Randomize the order
      j1 <- js[1]
      j2 <- js[2]

      # Proposal
      prop <- mcmc

      # Update the state. Save re-indexing for last!
      prop$n <- mcmc$n - 1
      prop$h[js] <- h
      prop$w[js] <- mcmc$w[js] + mcmc$w[i] + 1
      prop$d[h] <- mcmc$d[h] + 1

      w_j2 <- prop$w[j2] # Note the edge weight of j2

      prop <- update_genetics_downstream(prop, prop, j1, i)
      prop <- update_genetics_downstream(prop, prop, j2, i)

      ## For computing Hastings ratio...

      # To calculate Hastings ratio, generate a new genotype for i from scratch based on j1, j2
      geno <- create_genotype(prop, prop, i, j1, j2)
      backprop <- geno[[1]] # Proposed return to state with unobserved host
      ambiguous <- geno[[2]]

      # Get a list of "SNVs of interest": sites that change going from h to i, or i to a child of i
      interest <- unique(c(
        backprop$m01[[i]],
        backprop$mx1[[i]],
        unlist(backprop$m10[js]),
        unlist(backprop$m1y[js]),
        backprop$m10[[i]],
        backprop$mx0[[i]],
        unlist(backprop$m01[js]),
        unlist(backprop$m0y[js])
      ))

      # Ambiguous sites don't matter
      interest <- setdiff(interest, ambiguous)

      # Now, create two vectors, noting the classification of each SNV in "interest" for mcmc and for backprop
      class0 <- rep(0, length(interest))
      class1 <- rep(0, length(interest))

      class0[interest %in% union(mcmc$m01[[i]], mcmc$mx1[[i]])] <- 1
      class1[interest %in% union(backprop$m01[[i]], backprop$mx1[[i]])] <- 1

      # How many discrepancies are there?
      n_discrepancies <- sum(class0 != class1)

      # Maximum time i could be infected
      max_t <- min(mcmc$t[j1], mcmc$t[j2])

      # Maximum distance upstream we could propose a new node on "prop"
      max_dist <- min(prop$w[j1], prop$w[j2]) - 1

      # How many nodes along edge were added?
      change <- mcmc$w[i] + 1 # positive

      # Note the degree of h in the new config
      deg_h <- prop$d[h]

      ## Re-indexing: everyone above i steps down 1, i gets deleted
      prop$h[which(prop$h > i)] <- prop$h[which(prop$h > i)] - 1

      prop$h <- prop$h[-i]
      prop$w <- prop$w[-i]
      prop$t <- prop$t[-i]
      prop$m01 <- prop$m01[-i]
      prop$m10 <- prop$m10[-i]
      prop$m0y <- prop$m0y[-i]
      prop$m1y <- prop$m1y[-i]
      prop$mx0 <- prop$mx0[-i]
      prop$mx1 <- prop$mx1[-i]
      prop$mxy <- prop$mxy[-i]
      prop$d <- prop$d[-i]
      prop$g_lik <- prop$g_lik[-i]

      # Update j1 and j2 themselves
      if(j1 > i){
        j1 <- j1 - 1
      }
      if(j2 > i){
        j2 <- j2 - 1
      }

      # How many choices to make the reverse move?
      back_choices <- union(
        which(prop$d[1:data$n_obs] >= 2),
        which(prop$d >= 3)
      )

      # Again, nightmarish Hastings ratio...
      # Assume, technically, that j2 off of the edge from h to j1

      hastings <- (-lchoose(w_j2, change)) - # P(new -> old): choose [change] people to delete from the edge
        -lchoose(data$N, change) - # P(old -> new): choose [change] people to add back onto the edge
        log(length(back_choices)) - lchoose(deg_h, 2) - log(2) + # P(new -> old): choose from among #[back_choices] people to be h, and then pick j1, j2
        log(length(choices)) + log(2) - # P(old -> new): choose an unobserved host with exactly 2 children, then order them
        log(max_dist + 1) + # P(new -> old): uniform density for max_dist
        dbeta((mcmc$t[i] - mcmc$t[h]) / (max_t - mcmc$t[h]), mcmc$w[i] + 1, max_dist - mcmc$w[i] + 1) - # P(new -> old): beta density for t_i
        length(ambiguous) * log(2) + # P(new -> old): make a choice of the genotype at each ambiguous site
        n_discrepancies * log(data$eps) + # P(new -> old): each random swap has probability eps
        (length(interest) - n_discrepancies) * log(1 - data$eps) # P(new -> old): each lack of random swap has probability 1 - eps

      prop$e_lik <- e_lik(prop, data)
      update <- c(j1, j2) # For which hosts must we update the genomic likelihood? (note, now j1, j2 have different values)
      prop$g_lik[update] <- sapply(update, g_lik, mcmc = prop, data = data)
      prop$prior <- prior(prop)

      # Accept / reject
      if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior + hastings){
        return(prop)
      }else{
        return(mcmc)
      }
    }
  }
}


###


## Create a new unobserved node one step downstream, i.e.
# h -> j1 -> j2 becomes
# h -> i -> j1, j2
moves$create_downstream <- function(mcmc, data){

  # Are we creating or deleting an unobserved node?
  if(runif(1) < 1/2){
    create <- T
  }else{
    create <- F
  }

  if(create){
    # Select j1
    # Must be observed and have ancestor and have >= 1 children, or unobserved and >= 3 children (because it will lose a child)
    choices <- union(
      which(mcmc$d[1:data$n_obs] >= 1),
      which(mcmc$d >= 3)
    )
    choices <- setdiff(choices, 1)

    if(length(choices) == 0){
      return(mcmc)
    }else{
      j1 <- ifelse(length(choices) == 1, choices, sample(choices, 1))
      h <- mcmc$h[j1] # Ancestor of j1

      # Pick a child of this node
      j2 <- which(mcmc$h == j1)
      if(length(j2) > 1){
        j2 <- sample(j2, 1, replace = F)
      }

      js <- c(j1, j2)

      # How far upstream from h is the new node?
      # Maximum is w[j1] - 1
      max_dist <- mcmc$w[j1] - 1

      # If negative distance, can't make the move
      if(max_dist < 0){
        return(mcmc)
      }else{

        # New edge weight coming into i, the new host
        dist <- sample(0:max_dist, 1)
        i <- mcmc$n + 1

        # Maximum time i could be infected
        max_t <- mcmc$t[j1]

        # Proposal
        prop <- mcmc

        ## Update all entries in the Markov chain
        prop$n <- mcmc$n + 1
        prop$h[js] <- i
        prop$h[i] <- h
        prop$w[j1] <- mcmc$w[j1] - dist - 1
        prop$w[j2] <- mcmc$w[j2] + prop$w[j1] + 1 # Distance increases by the weight of j1, plus 1 for j1 itself
        prop$w[i] <- dist
        prop$t[i] <- mcmc$t[h] + (max_t - mcmc$t[h]) * rbeta(1, dist + 1, max_dist - dist + 1) # Weighted average
        prop$d[j1] <- mcmc$d[j1] - 1
        prop$d[i] <- 2

        # Create a genotype for the new host
        geno <- create_genotype_downstream(prop, mcmc, i, j1, j2)
        prop <- geno[[1]]
        ambiguous <- geno[[2]]

        ## Randomly swap the genotype at certain sites, to ensure reversibility
        # Probability of a swap is eps

        # Get a list of "SNVs of interest": sites that change going from h to i, or i to a child of i
        interest <- unique(c(
          prop$m01[[i]],
          prop$mx1[[i]],
          unlist(prop$m10[js]),
          unlist(prop$m1y[js]),
          prop$m10[[i]],
          prop$mx0[[i]],
          unlist(prop$m01[js]),
          unlist(prop$m0y[js])
        ))

        interest <- setdiff(interest, ambiguous) # No need to update ambiguous sites, because they're already 50/50

        swaps <- interest[runif(length(interest)) < data$eps]

        # Perform the flip
        for (snv in swaps) {
          prop <- flip_genotype(prop, prop, i, js, snv)
        }

        ## Okay, now we have the heinous task of computing the Hastings ratio.

        # How many nodes along edges were ADDED?
        change <- prop$w[j1] + 1 # (positive)

        # Assume, technically, that j2 moves onto the edge from h to j1

        hastings <- -lchoose(prop$w[j2], change) - # P(new -> old): probability we choose #[change] people to delete from the edge
          (-lchoose(data$N, change)) - # P(old -> new): choose [change] people to add back onto the edge
          log(sum(prop$d[(data$n_obs + 1):prop$n] == 2)) - log(2) + # P(new -> old): choose an unobserved host with exactly 2 children, then order them
          log(length(choices)) + log(sum(mcmc$h == j1, na.rm = T)) + # P(old -> new): choose from among #[choices] people to be j1, and then pick j2
          log(max_dist + 1) - # P(old -> new): uniform density for max_dist
          dbeta((prop$t[i] - mcmc$t[h]) / (max_t - mcmc$t[h]), dist + 1, max_dist - dist + 1) + # P(old -> new): beta density for t_i
          length(ambiguous) * log(2) - # P(old -> new): make a choice of the genotype at each ambiguous site
          length(swaps) * log(data$eps) - # P(old -> new): each random swap has probability eps
          (length(interest) - length(swaps)) * log(1 - data$eps) # P(old -> new): each lack of random swap has probability 1 - eps

        prop$e_lik <- e_lik(prop, data)
        update <- c(i, js) # For which hosts must we update the genomic likelihood?
        prop$g_lik[update] <- sapply(update, g_lik, mcmc = prop, data = data)
        prop$prior <- prior(prop)

        # Accept / reject
        if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior + hastings){
          return(prop)
        }else{
          return(mcmc)
        }
      }
    }
  }else{
    ## Getting rid of an unobserved node

    # Choices: any unobserved node with degree exactly 2
    choices <- which(mcmc$d == 2)
    choices <- choices[choices > data$n_obs]

    if(length(choices) == 0){
      return(mcmc)
    }else{
      # Select one to delete
      i <- ifelse(length(choices) == 1, choices, sample(choices, 1))

      # Ancestor
      h <- mcmc$h[i]

      # Children
      js <- sample(which(mcmc$h == i), 2) # Randomize the order
      j1 <- js[1]
      j2 <- js[2]

      # If weight of j2 not big enough, can't make the move
      if(mcmc$w[j2] <= mcmc$w[j1]){
        return(mcmc)
      }else{
        # Proposal
        prop <- mcmc

        # Update the state. Save re-indexing for last!
        prop$n <- mcmc$n - 1
        prop$h[j1] <- h
        prop$h[j2] <- j1
        prop$w[j1] <- mcmc$w[j1] + mcmc$w[i] + 1
        prop$w[j2] <- mcmc$w[j2] - mcmc$w[j1] - 1
        prop$d[j1] <- mcmc$d[j1] + 1

        w_j1 <- prop$w[j1] # Note the edge weight of j1, for Hastings ratio

        # j2 moves upstream to become child of j1
        prop <- update_genetics_upstream(prop, prop, j2, j1)
        # j1 moves downstream to become child of h
        prop <- update_genetics_downstream(prop, prop, j1, i)

        ## For computing Hastings ratio...

        # To calculate Hastings ratio, generate a new genotype for i from scratch based on j1, j2
        geno <- create_genotype_downstream(prop, prop, i, j1, j2)
        backprop <- geno[[1]] # Proposed return to state with unobserved host
        ambiguous <- geno[[2]]

        # Get a list of "SNVs of interest": sites that change going from h to i, or i to a child of i
        interest <- unique(c(
          backprop$m01[[i]],
          backprop$mx1[[i]],
          unlist(backprop$m10[js]),
          unlist(backprop$m1y[js]),
          backprop$m10[[i]],
          backprop$mx0[[i]],
          unlist(backprop$m01[js]),
          unlist(backprop$m0y[js])
        ))

        # Ambiguous sites don't matter
        interest <- setdiff(interest, ambiguous)

        # Now, create two vectors, noting the classification of each SNV in "interest" for mcmc and for backprop
        class0 <- rep(0, length(interest))
        class1 <- rep(0, length(interest))

        class0[interest %in% union(mcmc$m01[[i]], mcmc$mx1[[i]])] <- 1
        class1[interest %in% union(backprop$m01[[i]], backprop$mx1[[i]])] <- 1

        # How many discrepancies are there?
        n_discrepancies <- sum(class0 != class1)

        # Maximum time i could be infected
        max_t <- mcmc$t[j1]

        # Maximum distance upstream we could propose a new node on "prop"
        max_dist <- prop$w[j1] - 1

        # How many nodes along edge were deleted?
        change <- mcmc$w[j1] + 1 # positive

        # Note the degree of j1 in the new config
        deg_j1 <- prop$d[j1]

        ## Re-indexing: everyone above i steps down 1, i gets deleted
        prop$h[which(prop$h > i)] <- prop$h[which(prop$h > i)] - 1

        prop$h <- prop$h[-i]
        prop$w <- prop$w[-i]
        prop$t <- prop$t[-i]
        prop$m01 <- prop$m01[-i]
        prop$m10 <- prop$m10[-i]
        prop$m0y <- prop$m0y[-i]
        prop$m1y <- prop$m1y[-i]
        prop$mx0 <- prop$mx0[-i]
        prop$mx1 <- prop$mx1[-i]
        prop$mxy <- prop$mxy[-i]
        prop$d <- prop$d[-i]
        prop$g_lik <- prop$g_lik[-i]

        # Update j1 and j2 themselves
        if(j1 > i){
          j1 <- j1 - 1
        }
        if(j2 > i){
          j2 <- j2 - 1
        }

        # How many choices to make the reverse move?
        back_choices <- union(
          which(prop$d[1:data$n_obs] >= 1),
          which(prop$d >= 3)
        )
        back_choices <- setdiff(back_choices, 1)

        # Again, nightmarish Hastings ratio...

        hastings <- (-lchoose(data$N, change)) - # P(new -> old): choose [change] people to add back onto the edge
          -lchoose(w_j1, change) - # P(old -> new): choose [change] people to delete from the edge
          log(length(back_choices)) - log(deg_j1) + # P(new -> old): choose from among #[back_choices] people to be j1, and then pick j2
          log(length(choices)) + log(2) - # P(old -> new): choose an unobserved host with exactly 2 children, then order them
          log(max_dist + 1) + # P(new -> old): uniform density for max_dist
          dbeta((mcmc$t[i] - mcmc$t[h]) / (max_t - mcmc$t[h]), mcmc$w[i] + 1, max_dist - mcmc$w[i] + 1) - # P(new -> old): beta density for t_i
          length(ambiguous) * log(2) + # P(new -> old): make a choice of the genotype at each ambiguous site
          n_discrepancies * log(data$eps) + # P(new -> old): each random swap has probability eps
          (length(interest) - n_discrepancies) * log(1 - data$eps) # P(new -> old): each lack of random swap has probability 1 - eps

        prop$e_lik <- e_lik(prop, data)
        update <- c(j1, j2) # For which hosts must we update the genomic likelihood? (note, now j1, j2 have different values)
        prop$g_lik[update] <- sapply(update, g_lik, mcmc = prop, data = data)
        prop$prior <- prior(prop)

        # Accept / reject
        if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior + hastings){
          return(prop)
        }else{
          return(mcmc)
        }
      }
    }
  }
}









