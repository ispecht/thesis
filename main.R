### Execute large-scale outbreak reconstruction algorithm
set.seed(213)
## Libraries
library(ape)
library(Rcpp)
library(igraph)
library(ggraph)
library(parallel)

source("likelihood.R")
source("moves.R")
source("prior.R")
source("subroutines.R")
source("initialize.R")
source("global_mcmc.R")
source("local_mcmc.R")

init <- initialize()

mcmc <- init[[1]]
data <- init[[2]]

#load("./input_data_huge/data.RData")
#load("./input_data_huge/mcmc.RData")

### M-H algo
output <- list()

liks <- c()

for (r in 1:data$n_global) {

  # Make global moves
  mcmc <- global_mcmc(mcmc, data)

  # Chop up the tree into pieces
  breakdowns <- breakdown(mcmc, data)
  mcmcs <- breakdowns[[1]]
  datas <- breakdowns[[2]]

  # Run MCMC in parallel over each subtree
  all_res <- parallel::mclapply(
    1:data$n_subtrees,
    function(i, mcmcs, datas){
      local_mcmc(mcmcs[[i]], datas[[i]])
    },
    mcmcs = mcmcs,
    datas = datas,
    mc.cores = data$n_subtrees
  )
  #...or run in series
  # all_res <- list()
  # for (j in 1:data$n_subtrees) {
  #   all_res[[j]] <- local_mcmc(mcmcs[[j]], datas[[j]])
  # }

  # Amalgamate results of parallel MCMC run
  amalgam <- amalgamate(all_res, mcmcs, datas, mcmc, data)

  # Record amalgamated results, filtering to parameters of interest
  for (i in 1:length(amalgam)) {
    output <- c(output, list(
      amalgam[[i]][data$record]
    ))
  }

  # "mcmc" is now the most recent result
  mcmc <- amalgam[[length(amalgam)]]

  #print(r)

  liks <- c(liks, mcmc$e_lik + sum(mcmc$g_lik[2:mcmc$n]) + mcmc$prior)


  print(paste(r * data$n_local * data$n_subtrees, "iterations complete. Log-likelihood =", round(liks[r], 2)))
  print(plot_current(mcmc$h, data$n_obs))
  #print(mcmc$w)

  if(r == 10){
    data$n_subtrees <- 3
  }

}

# for (i in 1:1000) {
#   print(output[[i]]$mu)
# }








#print(mcmc$g_lik)
plot(liks)
plot(liks[50000:100000])

## Diagnostics
diagnostics <- list()
diagnostics$n <- c()
diagnostics$n_mut <- c() # Total number of mutations
diagnostics$adj <- matrix(0, nrow = data$n_obs, ncol = data$n_obs)
for (i in 1:length(res)) {
  diagnostics$n[i] <- res[[i]]$n
  diagnostics$n_mut[i] <- length(c(
    unlist(res[[i]]$m01),
    unlist(res[[i]]$m0y),
    unlist(res[[i]]$m10),
    unlist(res[[i]]$m1y)
  ))

}
plot(diagnostics$n)
plot(diagnostics$n_mut)
plot(diagnostics$n_mut[300:1000])

## Idea for visualization / summary: for each unobserved node, get a list of which nodes are upstream.
# Then take the MAP of the vector of ancestors concatenated with the vector of lists of upstream nodes.
diagnostics$h <- list()
diagnostics$adj <- matrix(0, nrow = max(diagnostics$n), ncol = max(diagnostics$n))
for (i in 1:length(res)) {
  h <- res[[i]]$h
  if(res[[i]]$n > data$n_obs){
    unobserved_anc <- which(h > data$n_obs)
    unobs <- unique(h[unobserved_anc])
    # Rename "unobs" as "sort(unobs)"
    sorted <- sort(unobs)
    h[unobserved_anc] <- sorted[match(h[unobserved_anc], unobs)]
  }
  diagnostics$h[[i]] <- h

  if(i > 0){
    diagnostics$adj[cbind(h[2:res[[i]]$n], 2:res[[i]]$n)] <- diagnostics$adj[cbind(h[2:res[[i]]$n], 2:res[[i]]$n)] + 1
  }
}
diagnostics$adj <- diagnostics$adj / length(res)

adj <- diagnostics$adj
adj[adj < 0.05] <- 0

g <- graph_from_adjacency_matrix(adj, mode = "directed", weighted = T)
color <- rep("orange", length(V(g)))
color[1:data$n_obs] <- "blue"
plot(
  g,
  #vertex.label=NA,
  vertex.label.cex = 0.4,
  vertex.label.color = 'white',
  vertex.label.family = 'sans',
  vertex.size=4,
  vertex.color = color,
  vertex.frame.color = "#00000000",
  edge.width=E(g)$weight*3,
  edge.color = rgb(0,0,0,E(g)$weight),
  edge.arrow.size = 0.5
)

#save(res, file = "res_12-28.RData")









