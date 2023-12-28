### Execute large-scale outbreak reconstruction algorithm
set.seed(213)
## Libraries
library(ape)
library(Rcpp)
library(igraph)
library(ggraph)
library(gganimate)

source("likelihood.R")
source("moves.R")
source("prior.R")
source("subroutines.R")

## Filters
filters <- list(
  af = 0.03,
  dp = 100,
  sb = 10
)
init_mst <- F

## Data Processing

# Load the reference sequence
ref_genome <- read.FASTA("./input_data/ref.fasta")

# Length of genome
n_bases <- length(ref_genome[[1]])

# Load the FASTA of sequences
fasta <- read.FASTA("./input_data/aligned.fasta")

# The first genome is itself the ref genome
fasta <- c(ref_genome, fasta)

# Number of samples
n <- length(fasta)

# Names of sequences
names <- names(fasta)

# VCF files present
vcfs <- list.files("./input_data/vcf/")

# Date
date <- read.csv("input_data/date.csv")
s <- c()
for (i in 1:n) {
  # Check if we have a VCF file
  included <- grepl(names[i], date[,1])
  if(sum(included) >= 2){
    stop(paste("Multiple sample collection dates found for sequence", names[i]))
  }else if(sum(included) == 1){
    s[i] <- date[,2][included]
  }else{
    s[i] <- NA
  }
}
s[1] <- 0

## List of SNVs present per sample
snvs <- list()
for (i in 1:n) {
  # Check if we have a VCF file
  included <- grepl(names[i], vcfs)
  if(sum(included) >= 2){
    stop(paste("Multiple VCF files found for sequence", names[i]))
  }else if(sum(included) == 1){
    vcf <- read.table(paste0("./input_data/vcf/", vcfs[included]))
    snvs[[i]] <- genetic_info(ref_genome[[1]], fasta[[i]], filters = filters, vcf = vcf)
  }else{
    snvs[[i]] <- genetic_info(ref_genome[[1]], fasta[[i]], filters = filters)
  }
}

# Compile vectors of all positions with SNVs (may have duplicates) and all SNVs (unique)
all_pos <- c()
all_snv <- c()
for (i in 1:n) {
  calls <- snvs[[i]]$snv$call
  new <- which(!(calls %in% all_snv))
  all_snv <- c(all_snv, calls[new])
  all_pos <- c(all_pos, snvs[[i]]$snv$pos[new])
  if(!is.null(snvs[[i]]$isnv)){
    calls <- snvs[[i]]$isnv$call
    new <- which(!(calls %in% all_snv))
    all_snv <- c(all_snv, calls[new])
    all_pos <- c(all_pos, snvs[[i]]$isnv$pos[new])
  }
}

# Remove irrelevant missing site info; add SNVs with missing data
for (i in 1:n) {
  # Which positions in all_pos are detected in the missing sites in the sample?
  present <- which(all_pos %in% snvs[[i]]$missing$pos)
  # Change missing$pos to be a vector of these sites (may have duplicates)
  snvs[[i]]$missing$pos <- all_pos[present]
  # Add a new vector of the SNVs for which there's no information
  snvs[[i]]$missing$call <- all_snv[present]
}

### Initialize MCMC and data

## For MCMC initialization: minimum spanning tree
if(init_mst){
  snv_dist <- ape::dist.dna(fasta, "N", as.matrix = T)
  tree <- ape::mst(snv_dist)
  init_h <- adj_to_anc(tree, 1)
}

data <- list()
data$s <- s
data$N <- 10000 #population size
data$n_obs <- n # number of observed hosts, plus 1 (index case)
data$n_bases <- n_bases
data$snvs <- snvs
data$eps <- 0.005 # Explore/exploit tradeoff for genotypes of new nodes
data$p_move <- 0.6
data$tau = 0.2

mcmc <- list()
mcmc$n <- n # number of tracked hosts
mcmc$h <- rep(1, n) # ancestors; initialized to index case
mcmc$w <- rep(0, n) # edge weights; initialized to 0
mcmc$w[1] <- 0 # For convenience
mcmc$h[1] <- NA
mcmc$t <- s - 5 # time of contracting
mcmc$m01 <- list() # fixed mutations added in each transmission link
mcmc$m10 <- list() # fixed mutations deleted in each transmission link
mcmc$m0y <- list() # 0% -> y%, 0 < y < 100
mcmc$m1y <- list() # 100% -> y%, 0 < y < 100
mcmc$mx0 <- list() # x% -> 0%, 0 < x < 100
mcmc$mx1 <- list() # x% -> 100%, 0 < x < 100
mcmc$mxy <- list() # x% -> y%, 0 < x < 100, 0 < y < 100
for (i in 1:n) {
  mcmc$m01[[i]] <- snvs[[i]]$snv$call
  mcmc$m10[[i]] <- character(0)
  mcmc$m0y[[i]] <- snvs[[i]]$isnv$call
  mcmc$m1y[[i]] <- character(0)
  mcmc$mx0[[i]] <- character(0)
  mcmc$mx1[[i]] <- character(0)
  mcmc$mxy[[i]] <- character(0)
}

if(init_mst){
  gens <- generations(init_h, 1)
  max_t <- min(s[2:n] - 5)
  for (g in 2:length(gens)) {
    for (i in gens[[g]]) {

      if(g >= 3){
        anc <- ancestry(init_h, i)
        for (j in 2:(length(anc) - 1)) {
          mcmc <- update_genetics_upstream(mcmc, mcmc, i, anc[j])
          mcmc$m01[[j]] <- setdiff(mcmc$m01[[j]], snvs[[j]]$missing$call) # Remove calls for missing positions
          mcmc$m10[[j]] <- setdiff(mcmc$m10[[j]], snvs[[j]]$missing$call)
        }
      }
      mcmc$t[i] <- max_t - 5*(length(gens) - g)
    }
  }
  mcmc$h <- init_h
}

mcmc$b <- 0.95 # Probability bottleneck has size 1
mcmc$a_g <- 5 # shape parameter of the generation interval
mcmc$lambda_g <- 1 # rate parameter of the generation interval. FOR NOW: fixing at 1.
mcmc$a_s <- 5 # shape parameter of the sojourn interval
mcmc$lambda_s <- 1 # rate parameter of the sojourn interval. FOR NOW: fixing at 1.
mcmc$mu <- 1e-6 # mutation rate, sites/day
mcmc$p <- 1e-6 # mutation rate, sites/cycle
mcmc$v <- 1000 # burst size
mcmc$rho <- 0.1 # first parameter, NBin offspring distribution (overdispersion param)
mcmc$psi <- 0.1 / (2.5 + 0.1) # second parameter, NBin offspring distribution (computed in terms of R0)

# Functions of MCMC params
mcmc$d <- sapply(1:n, function(x){sum(mcmc$h[2:n] == x)}) # Node degrees

# Also track the epidemiological and genomic likelihoods, and prior
# The genomic likelihood we will store on a per-person basis, for efficiency purposes
mcmc$e_lik <- e_lik(mcmc, data)
mcmc$g_lik <- c(NA, sapply(2:n, g_lik, mcmc = mcmc, data = data))
mcmc$prior <- prior(mcmc)

### M-H algo
liks <- c()
N_iters <- 1000
res <- list()
for (r in 1:N_iters) {
  mcmc <- moves$w(mcmc, data)
  mcmc <- moves$t(mcmc, data)
  mcmc <- moves$b(mcmc, data)
  #mcmc <- moves$a_g(mcmc, data)
  #mcmc <- moves$a_s(mcmc, data)
  mcmc <- moves$mu(mcmc, data)
  mcmc <- moves$p(mcmc, data)
  mcmc <- moves$v(mcmc, data)
  #mcmc <- moves$rho(mcmc, data)
  #mcmc <- moves$psi(mcmc, data)
  #print(mcmc$w)
  #print(mcmc$b)
  mcmc <- moves$swap(mcmc, data, exchange_children = T)
  mcmc <- moves$h_step(mcmc, data, resample_t = T)
  mcmc <- moves$genotype(mcmc, data)
  mcmc <- moves$h_global(mcmc, data)
  mcmc <- moves$create(mcmc, data)

  mcmc <- moves$w_t(mcmc, data)
  mcmc <- moves$h_step(mcmc, data)
  mcmc <- moves$swap(mcmc, data)
  mcmc <- moves$h_step(mcmc, data, resample_t = T, resample_w = T)

  if(mcmc$n > data$n_obs){
    if(any(mcmc$d[(data$n_obs + 1):mcmc$n] < 2)){
      print(r)
    }
  }


  # print(mcmc$h)
  # print(length(c(
  #   unlist(mcmc$m01),
  #   unlist(mcmc$m0y),
  #   unlist(mcmc$m10),
  #   unlist(mcmc$m1y)
  # )))
  #print(mcmc$b)


  liks <- c(liks, mcmc$e_lik + sum(mcmc$g_lik[2:mcmc$n]) + mcmc$prior)

  if(r %% 100 == 0){
    res <- c(res, list(mcmc))
    print(paste(r, "iterations complete. Log-likelihood =", round(liks[r], 2)))
    #print(plot_current(mcmc$h, data$n_obs))
  }


}

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



