# MIT License
#
# Copyright (c) 2023 Ivan Specht
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

### MCMC moves

moves <- list()

## Update one of the w_i's by adding or subtracting either:
# rounded N(0, sqrt(delta_t * lambda_g / a_g))  (strategic) OR
# N(0, 3) (random)
moves$w <- function(mcmc, data){
  # Choose random host with ancestor
  i <- sample(2:mcmc$n, 1)
  h <- mcmc$h[i]

  delta_t <- mcmc$t[i] - mcmc$t[h]

  # Proposal
  prop <- mcmc

  if(delta_t > 50){
    change <- round(rnorm(1, 0, sqrt(delta_t * mcmc$lambda_g / mcmc$a_g)))
  }else{
    change <- round(rnorm(1, 0, 3))
  }

  prop$w[i] <- mcmc$w[i] + change

  prop$e_lik <- e_lik(prop, data)
  prop$g_lik[i] <- g_lik(prop, data, i)
  prop$prior <- prior(prop)

  if(data$pooled_coalescent){
    # With new coalescent:
    hastings <- 0
  }else{
    # Here the Hastings ratio is not 1, because we're imagining that the newly-added / deleted intermediate hosts are labeled.
    if(change > 0){
      # If added new hosts, P(new to old) is probability of correctly selecting the new hosts to delete
      # P(old to new) is probability choosing the correct new hosts to add from larger population
      hastings <- -
        (-lchoose(data$N, change))  # Choose from N, order them, and place them on the edge
    }else if(change < 0){
      hastings <- -lchoose(data$N, -change)
    }else{
      hastings <- 0
    }
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
  i <- sample(setdiff(2:mcmc$n, data$frozen), 1)
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
  choices <- (2:mcmc$n)[which(mcmc$d[2:mcmc$n] > 0)]
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
    if(data$pooled_coalescent){
      hastings <- 0
    }else{
      hastings <- 0
      if(length(js) > 1){
        if(delta > 0){
          # In this case, we lose weight from the edges leading into the js, so:
          hastings <- hastings - lchoose(data$N, delta) * (length(js) - 1)  # P(new -> old): add the correct nodes from the population onto each of the #[js] - 1 edges
             # P(old -> new): delete the correct nodes from js[2], js[3], ...

        }
        if(delta < 0){
          hastings <- hastings + # P(new -> old): delete the correct nodes from js[2], js[3], ...
            lchoose(data$N, -delta) * (length(js) - 1)  # P(old -> new): add the correct nodes from the population to each of js[2], js[3], ...
        }
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
  prop$b <- rnorm(1, mcmc$b, 0.1)
  prop$e_lik <- e_lik(prop, data)
  prop$g_lik[2:mcmc$n] <- sapply(2:mcmc$n, g_lik, mcmc = prop, data = data)
  prop$prior <- prior(prop)

  if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior){
    return(prop)
  }else{
    return(mcmc)
  }
}

## Update lambda using a N(0,0.5^2) proposal density
moves$lambda <- function(mcmc, data){
  # Proposal
  prop <- mcmc
  prop$lambda <- rnorm(1, mcmc$lambda, 0.5)
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
  prop$a_g <- rnorm(1, mcmc$a_g, 0.5)

  # Also update reproductive number and hence psi to maintain constant growth rate
  prop$psi <- prop$rho / (exp((prop$a_g / prop$lambda_g) * data$growth) + prop$rho) # second parameter, NBin offspring distribution (computed in terms of R0)

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
  if(data$virus == "SARS-CoV-2"){
    sd <- 5e-7
  }else if(data$virus == "H5N1"){
    sd <- 1e-6
  }
  prop$mu <- rnorm(1, mcmc$mu, sd)
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
  prop$p <- rnorm(1, mcmc$p, 5e-7)
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
  prop$v <- round(rnorm(1, mcmc$v, 1000))
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
#### May need to check hastings ratio here...
moves$genotype <- function(mcmc, data){

  # Choose random host with ancestor
  i <- sample(setdiff(2:mcmc$n, data$frozen), 1)
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
moves$h_step <- function(mcmc, data, upstream = TRUE, resample_t = FALSE, resample_w = FALSE){
  # Choose random host with ancestor
  if(resample_t){
    i <- sample(setdiff(2:mcmc$n, data$frozen), 1)
  }else{
    i <- sample(2:mcmc$n, 1)
  }

  h_old <- mcmc$h[i]

  # Proposal
  prop <- mcmc

  # Are we going upstream or downstream?
  #upstream <- runif(1) < 1/2

  if(upstream){
    # Who are the other children of h_old?
    children <- setdiff(which(mcmc$h == h_old), i)

    # What's the maximum time at which i can be infected?
    max_t <- get_max_t(mcmc, data, i)

    # Which ones have a compatible time of infection?
    children <- children[mcmc$t[children] < max_t]

    # Children not allowed to be frozen
    children <- setdiff(children, data$frozen)

    # If no valid children, reject
    # Also reject if h_old is not observed and has <= 2 total children, because then we can't remove one
    if(length(children) == 0 | (h_old > data$n_obs & length(which(mcmc$h == h_old)) <= 2)){
      return(mcmc)
    }else{

      # Pick one
      h_new <- ifelse(length(children) == 1, children, sample(children, 1))

      prop <- shift_upstream(prop, data, i, h_old, h_new, resample_t, resample_w)

      # What's the change in edge weight for i?
      change <- prop$w[i] - mcmc$w[i]

      update <- i
      # If updating t, need to also change genomic likelihood of children of i
      if(resample_t){
        update <- c(update, which(mcmc$h == i))
      }

      if(data$pooled_coalescent){
        hastings <- 0
      }else{
        hastings <- 0
        if(change < 0){
          hastings <- hastings - lchoose(data$N, -change)  # P(new -> old): choose [change] people to add back onto the edge
          # P(old -> new): choose [change] people to delete from the edge
        }
        if(change > 0){
          hastings <- hastings - # P(new -> old): choose [change] people to delete from the edge
            (-lchoose(data$N, change))   # P(old -> new): choose [change] people to add back onto the edge
        }
      }

      hastings <- hastings + log(length(children)) # P(new -> old): 1; P(old -> new): choose from among #[children] people to be h_new

      if(resample_t){
        hastings <- hastings - log(max_t - mcmc$t[h_old]) + # P(new -> old): uniform draw of time of infection
          log(max_t - prop$t[h_new]) # P(old -> new): uniform draw of time of infection
      }
      if(resample_w){ # Poisson draw for edge weights
        hastings <- hastings + dpois(mcmc$w[i], (mcmc$t[i] - mcmc$t[h_old]) * mcmc$lambda_g / mcmc$a_g, log = T) -
          dpois(prop$w[i], (prop$t[i] - prop$t[h_new]) * mcmc$lambda_g / mcmc$a_g, log = T)
      }

      return(accept_or_reject(prop, mcmc, data, update, hastings))

    }
  }else{
    if(h_old == 1 | (h_old > data$n_obs & length(which(mcmc$h == h_old)) <= 2)){
      # If no downstream move, reject
      # Also reject if h_old is not observed and has <= 2 total children, because then we can't remove one
      return(mcmc)
    }else{

      # New ancestor of i is ancestor's ancestor
      h_new <- mcmc$h[h_old]

      prop <- shift_downstream(prop, data, i, h_old, h_new, resample_t, resample_w)

      # What's the change in edge weight for i? (Positive)
      change <- mcmc$w[h_old] + 1

      ## Compute the number of possible children who could be chosen by i in the new config
      # Who are the other children of h_old?
      children <- setdiff(which(prop$h == h_new), i)

      # What's the maximum time at which i can be infected?
      max_t <- get_max_t(mcmc, data, i)

      # Which ones have a lesser time of infection than max_t?
      children <- children[prop$t[children] < max_t]

      # Children can't be frozen
      children <- setdiff(children, data$frozen)

      # What's the change in edge weight for i?
      change <- prop$w[i] - mcmc$w[i]

      update <- i
      # If updating t, need to also change genomic likelihood of children of i
      if(resample_t){
        update <- c(update, which(mcmc$h == i))
      }

      if(data$pooled_coalescent){
        hastings <- 0
      }else{
        hastings <- 0
        if(change < 0){
          hastings <- hastings - lchoose(data$N, -change)   # P(new -> old): choose [change] people to add back onto the edge
             # P(old -> new): choose [change] people to delete from the edge
        }
        if(change > 0){
          hastings <- hastings  - # P(new -> old): choose [change] people to delete from the edge
            (-lchoose(data$N, change)) # P(old -> new): choose [change] people to add back onto the edge
        }
      }

      hastings <- hastings - log(length(children)) # P(new -> old): choose from among #[children] people to be h_new; P(old -> new): 1

      if(resample_t){
        hastings <- hastings - log(max_t - mcmc$t[h_old]) + # P(new -> old): uniform draw of time of infection
          log(max_t - prop$t[h_new]) # P(old -> new): uniform draw of time of infection
      }
      if(resample_w){ # Poisson draw for edge weights
        hastings <- hastings + dpois(mcmc$w[i], (mcmc$t[i] - mcmc$t[h_old]) * mcmc$lambda_g / mcmc$a_g, log = T) -
          dpois(prop$w[i], (prop$t[i] - prop$t[h_new]) * mcmc$lambda_g / mcmc$a_g, log = T)
      }

      return(accept_or_reject(prop, mcmc, data, update, hastings))
    }
  }
}

## Global change in ancestor
# Importance sampling based on other nodes with similar additions / deletions
moves$h_global <- function(mcmc, data){
  # Sample any node with ancestor
  i <- sample(2:mcmc$n, 1)
  h_old <- mcmc$h[i]

  # Nodes which are infected earlier than i
  choices <- which(mcmc$t < mcmc$t[i])

  choices <- setdiff(choices, data$frozen)

  if(length(choices) == 0 | (h_old > data$n_obs & mcmc$d[h_old] <= 2)){
    return(mcmc)
  }else{

    # "Score" the choices: shared iSNV = +1
    scores <- softmax(sapply(choices, score, mcmc=mcmc, i=i), data$tau)

    h_new <- ifelse(length(choices) == 1, choices, sample(choices, 1, prob = scores))

    # Find the path from h_old to h_new
    route <- paths(mcmc$h, h_old, h_new)
    down <- route[[1]]
    up <- route[[2]]

    prop <- mcmc

    # If length of down < 2, don't need to do anything
    if(length(down) >= 2){
      for (j in 2:length(down)) {
        prop <- shift_downstream(prop, data, i, down[j-1], down[j])
      }
    }
    if(length(up) >= 2){
      for (j in 2:length(up)) {
        prop <- shift_upstream(prop, data, i, up[j-1], up[j])
      }
    }

    prop$e_lik <- e_lik(prop, data)
    update <- i
    prop$g_lik[update] <- sapply(update, g_lik, mcmc = prop, data = data)
    prop$prior <- prior(prop)

    rev_scores <- softmax(sapply(choices, score, mcmc=prop, i=i), data$tau)
    hastings <- log(rev_scores[which(choices == h_old)]) - log(scores[which(choices == h_new)])

    if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior + hastings){
      #print("way to go!")
      return(prop)
    }else{
      return(mcmc)
    }
  }
}

## The swap
## Switch h -> i -> j to
## h -> j -> i
moves$swap <- function(mcmc, data, exchange_children = FALSE){
  # Choose host with a parent and a grandparent
  choices <- which(mcmc$h != 1)
  choices <- setdiff(choices, data$frozen)

  if(length(choices) == 0){
    return(mcmc)
  }else{
    # Pick j
    j <- ifelse(length(choices) == 1, choices, sample(choices, 1))

    # Pick i
    i <- mcmc$h[j]

    # For exchange_children = FALSE
    # If i unobserved, i must have at least 3 children, because losing one
    # For exchange_children = TRUE
    # If i unobserved, j must have at least 2 children
    # If j unobserved, i must have at least 2 children (including j)

    if(
      (exchange_children == F & i > data$n_obs & mcmc$d[i] < 3) |
      (exchange_children == T & i > data$n_obs & mcmc$d[j] < 2) |
      (exchange_children == T & j > data$n_obs & mcmc$d[i] < 2)
    ){
      return(mcmc)
    }else{
      # Pick h
      h <- mcmc$h[i]

      # Children of each
      children_i <- setdiff(which(mcmc$h == i), j)
      children_j <- which(mcmc$h == j)

      # Update the state
      prop <- mcmc
      prop <- shift_downstream(prop, data, j, i, h) # Shift j from i onto h
      prop <- shift_upstream(prop, data, i, h, j) # Shift i from h onto j
      prop$w[j] <- mcmc$w[i] # Swapping edge weights
      prop$w[i] <- mcmc$w[j]
      prop$t[j] <- mcmc$t[i] # Swapping time of infection
      prop$t[i] <- mcmc$t[j]

      if(exchange_children){
        for (k in children_i) {
          prop <- shift_downstream(prop, data, k, i, j)
          prop$w[k] <- mcmc$w[k] # Keep edge weight the same
        }
        for (k in children_j) {
          prop <- shift_upstream(prop, data, k, j, i)
          prop$w[k] <- mcmc$w[k] # Keep edge weight the same
        }
      }


      prop$e_lik <- e_lik(prop, data)
      update <- c(i, j, children_i, children_j) # For which hosts must we update the genomic likelihood?
      prop$g_lik[update] <- sapply(update, g_lik, mcmc = prop, data = data)
      prop$prior <- prior(prop)

      # Accept / reject
      if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior){
        # if(exchange_children){
        #   print("BASED")
        # }
        return(prop)
      }else{
        return(mcmc)
      }
    }
  }
}


## Create / remove a node
moves$create <- function(mcmc, data, create = T, upstream = T){
  # Are we creating or deleting an unobserved node?
  # if(runif(1) < 1/2){
  #   create <- T
  # }else{
  #   create <- F
  # }
  #
  # # Are we moving nodes onto the new node upstream or downstream?
  # if(runif(1) < 1/2){
  #   upstream <- T
  # }else{
  #   upstream <- F
  # }

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
      if(length(j2s) == 0 | (upstream & h > data$n_obs & length(j2s) == length(kids)) | (!upstream & j1 > data$n_obs & length(j2s) >= length(kids) - 1)){
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
          #prop$cluster <- c(2:mcmc$n, i)

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
              prop <- shift_upstream(prop, data, j, h, i)
            }
          }else{
            prop <- shift_upstream(prop, data, j1, h, i)
            for (j2 in j2s) {
              prop <- shift_downstream(prop, data, j2, j1, i)
            }
          }

          ## Create new genotype for i
          geno <- genotype(prop, i, js, data$eps)
          prop <- geno[[1]]
          log_p <- geno[[2]]

          ## Hastings ratio
          if(data$pooled_coalescent){
            hastings <- 0
          }else{
            # Delta is change in number of hosts along edges into j2s
            hastings <- 0
            if(upstream){
              delta <- dist + 1
              # In this case, we lose weight from the edges leading into the j2s, so:
              hastings <- hastings - lchoose(data$N, delta) * length(j2s)  # P(new -> old): add the correct nodes from the population onto each of the #[js] - 1 edges
                 # P(old -> new): delete the correct nodes from js[2], js[3], ...
            }else{
              delta <- prop$w[j1] + 1
              hastings <- hastings  + # P(new -> old): delete the correct nodes from js[2], js[3], ...
                lchoose(data$N, delta) * length(j2s)  # P(old -> new): add the correct nodes from the population to each of js[2], js[3], ...
            }
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

      if((!upstream & any(mcmc$w[j2s] <= mcmc$w[j1])) | (!upstream & j1 %in% data$frozen)){
        return(mcmc)
      }else{

        prop <- mcmc

        # Put all children of i onto h or j1
        if(upstream){
          for (j in js) {
            prop <- shift_downstream(prop, data, j, i, h)
          }
        }else{
          for (j2 in j2s) {
            prop <- shift_upstream(prop, data, j2, i, j1)
          }
          # Put j1 onto h
          prop <- shift_downstream(prop, data, j1, i, h)
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
        if(data$pooled_coalescent){
          hastings <- 0
        }else{
          hastings <- 0
          #Delta is change in number of hosts along edges into j2s
          if(upstream){
            delta <- mcmc$w[i] + 1
            # In this case, we lose weight from the edges leading into the j2s, so:
            hastings <- hastings  + # P(new -> old): delete the correct nodes from js[2], js[3], ...
              lchoose(data$N, delta) * length(j2s)  # P(old -> new): add the correct nodes from the population onto each of the #[js] - 1 edges

          }else{
            delta <- mcmc$w[j1] + 1
            hastings <- hastings - lchoose(data$N, delta) * length(j2s)  # P(new -> old): add the correct nodes from the population to each of js[2], js[3], ...
            # P(old -> new): delete the correct nodes from js[2], js[3], ...
          }
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
        #prop$cluster <- setdiff(prop$cluster, i) # Delete i, then shift everyone bigger down by 1
        #prop$cluster[which(prop$cluster > i)] <-  prop$cluster[which(prop$cluster > i)] - 1

        # Update the js
        js[js > i] <- js[js > i] - 1

        update <- js
        return(accept_or_reject(prop, mcmc, data, update, hastings))

      }
    }
  }
}




