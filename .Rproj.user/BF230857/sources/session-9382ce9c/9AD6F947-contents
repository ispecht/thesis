# Initialize mcmc (current state of the markov chain) and data
initialize <- function(){
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
    snv_dist <- ape::dist.dna(fasta, "N")
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
  data$n_cores <- parallel::detectCores()
  # Number of subtrees to chop into is n_cores, as long as each subtree has at least 100 people
  data$n_subtrees <- max(min(data$n_cores, floor(n / 100)), 1)
  data$n_global <- 50 # Number of global moves
  data$n_local <- 100 # Number of local moves per global move
  data$sample_every <- 100 # Per how many local moves do we draw one sample?
  #data$n_subtrees <- 3

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

  return(list(mcmc, data))
}
