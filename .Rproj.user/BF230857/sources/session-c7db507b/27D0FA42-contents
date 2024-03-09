### Visualize MCMC output

library(ggplot2)
library(ggraph)
library(igraph)

# Load in results
load("output.RData")

mus <- c()
ps <- c()
bs <- c()
vs <- c()
for (i in 1:10000) {
  mus[i] <- output[[i]]$mu
  ps[i] <- output[[i]]$p
  bs[i] <- output[[i]]$b
  vs[i] <- output[[i]]$v
}
