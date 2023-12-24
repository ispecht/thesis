### Helper functions

# Functions to convert between raw format and nucleotide letters:
base_to_raw <- function(b){
  out <- rep(as.raw(04), length(b))
  out[b == "A"] <- as.raw(136)
  out[b == "C"] <- as.raw(40)
  out[b == "G"] <- as.raw(72)
  out[b == "T"] <- as.raw(24)
  return(out)
}

raw_to_base <- function(r){
  out <- rep("N", length(r))
  out[r == as.raw(136)] <- "A"
  out[r == as.raw(40)] <- "C"
  out[r == as.raw(72)] <- "G"
  out[r == as.raw(24)] <- "T"
  return(out)
}

genetic_info <- function(seq1, seq2, filters, vcf = NULL){
  # List output: list of SNVs, iSNVs, and positions with no information
  out <- list()

  ## Get iSNVs from VCF, if provided
  if(!is.null(vcf)){

    out$isnv <- list()

    # Position
    pos <- vcf$V2
    # Allele on reference genome
    ref <- vcf$V4
    # Alternate allele
    alt <- vcf$V5
    # Info column
    info <- vcf$V8

    # Read depth
    dp <- gsub(".*DP=", "", info)
    dp <- sub(";.*", "", dp)
    dp <- as.numeric(dp)

    # Allele frequency
    af <- gsub(".*AF=", "", info)
    af <- sub(";.*", "", af)
    af <- as.numeric(af)

    # Strand bias
    sb <- gsub(".*;SB=", "", info)
    sb <- sub(";.*", "", sb)
    sb <- as.numeric(sb)

    # Which sites pass the filters?
    keep <- which(dp >= filters$dp & sb < filters$sb & af >= filters$af & af <= 1 - filters$af)
    pos <- pos[keep]
    ref <- ref[keep]
    alt <- alt[keep]
    af <- af[keep]

    out$isnv$call <- paste0(ref, pos, alt)
    out$isnv$pos <- pos
    out$isnv$af <- af

  }

  out$snv <- list()

  ## Get SNVs from FASTA
  snv_pos <- which(
    seq1 != seq2 &
      seq1 %in% c(as.raw(136), as.raw(40), as.raw(72), as.raw(24)) &
      seq2 %in% c(as.raw(136), as.raw(40), as.raw(72), as.raw(24))
  )

  ## Get missing sites from FASTA in seq2
  missing_pos <- which(
    !(seq2 %in% c(as.raw(136), as.raw(40), as.raw(72), as.raw(24)))
  )

  # Remove positions already accounted for in VCF, if provided
  if(!is.null(vcf)){
    snv_pos <- setdiff(snv_pos, pos)
    missing_pos <- setdiff(missing_pos, pos)
  }

  old <- raw_to_base(seq1[snv_pos])
  new <- raw_to_base(seq2[snv_pos])

  out$snv$call <- paste0(old, snv_pos, new)
  out$snv$pos <- snv_pos
  out$missing <- list()
  out$missing$pos <- missing_pos

  return(out)

}

# Distribution of de novo iSNVs
denovo <- function(x, p, log = FALSE){
  k <- 1/sqrt(p)
  if(log){
    log(1-(1-x)^k * (1 + k*x)) - log(k) - 2*log(x)
  }else{
    (1-(1-x)^k * (1 + k*x)) / (k*x^2)
  }
}

# Update genetics for the following topological move:
# From g -> i, g -> h
# To g -> h -> i
update_genetics_upstream <- function(prop, mcmc, i, h){
  # Everything that doesn't stay the same in i
  all_i <- unique(c(
    mcmc$m01[[i]],
    mcmc$m10[[i]],
    mcmc$m0y[[i]],
    mcmc$m1y[[i]],
    mcmc$mx0[[i]],
    mcmc$mx1[[i]],
    mcmc$mxy[[i]]
  ))

  # Everything that doesn't stay the same in h
  all_h <- unique(c(
    mcmc$m01[[h]],
    mcmc$m10[[h]],
    mcmc$m0y[[h]],
    mcmc$m1y[[h]],
    mcmc$mx0[[h]],
    mcmc$mx1[[h]],
    mcmc$mxy[[h]]
  ))

  prop$m01[[i]] <- setdiff(mcmc$m01[[i]], all_h) # 00 in h, 01 in i
  prop$m01[[i]] <- union(prop$m01[[i]], intersect(mcmc$mx0[[h]], mcmc$mx1[[i]])) # x0 in h, x1 in i
  prop$m01[[i]] <- union(prop$m01[[i]], setdiff(mcmc$m10[[h]], all_i)) # 10 in h, 11 in i

  prop$m10[[i]] <- setdiff(mcmc$m10[[i]], all_h) # 11 in h, 10 in i
  prop$m10[[i]] <- union(prop$m10[[i]], intersect(mcmc$mx1[[h]], mcmc$mx0[[i]])) # x1 in h, x0 in i
  prop$m10[[i]] <- union(prop$m10[[i]], setdiff(mcmc$m01[[h]], all_i)) # 01 in h, 00 in i

  prop$m0y[[i]] <- setdiff(mcmc$m0y[[i]], all_h) # 00 in h, 0y in i
  prop$m0y[[i]] <- union(prop$m0y[[i]], intersect(mcmc$m10[[h]], mcmc$m1y[[i]])) # 10 in h, 1y in i
  prop$m0y[[i]] <- union(prop$m0y[[i]], intersect(mcmc$mx0[[h]], mcmc$mxy[[i]])) # x0 in h, xy in i

  prop$m1y[[i]] <- setdiff(mcmc$m1y[[i]], all_h) # 11 in h, 1y in i
  prop$m1y[[i]] <- union(prop$m1y[[i]], intersect(mcmc$m01[[h]], mcmc$m0y[[i]])) # 01 in h, 0y in i
  prop$m1y[[i]] <- union(prop$m1y[[i]], intersect(mcmc$mx1[[h]], mcmc$mxy[[i]])) # x1 in h, xy in i

  prop$mx0[[i]] <- intersect(mcmc$mxy[[h]], mcmc$mx0[[i]]) # xy in h, x0 in i
  prop$mx0[[i]] <- union(prop$mx0[[i]], setdiff(mcmc$m0y[[h]], all_i)) # 0y in h, 00 in i
  prop$mx0[[i]] <- union(prop$mx0[[i]], intersect(mcmc$m1y[[h]], mcmc$m10[[i]])) # 1y in h, 10 in i

  prop$mx1[[i]] <- intersect(mcmc$mxy[[h]], mcmc$mx1[[i]]) # xy in h, x1 in i
  prop$mx1[[i]] <- union(prop$mx1[[i]], setdiff(mcmc$m1y[[h]], all_i)) # 1y in h, 11 in i
  prop$mx1[[i]] <- union(prop$mx1[[i]], intersect(mcmc$m0y[[h]], mcmc$m01[[i]])) # 0y in h, 01 in i

  prop$mxy[[i]] <- intersect(mcmc$mxy[[h]], mcmc$mxy[[i]]) # xy in h, xy in i
  prop$mxy[[i]] <- union(prop$mxy[[i]], intersect(mcmc$m1y[[h]], mcmc$m1y[[i]])) # 1y in h, 1y in i
  prop$mxy[[i]] <- union(prop$mxy[[i]], intersect(mcmc$m0y[[h]], mcmc$m0y[[i]])) # 0y in h, 0y in i

  # Whew.
  return(prop)
}

# Wrap as a function: switch from
# h_old -> i, h_old -> h_new to
# h_old -> h_new -> i
shift_upstream <- function(mcmc, i, h_old, h_new){
  # Update all necessary components of MCMC
  mcmc$h[i] <- h_new # Update the ancestor
  mcmc$w[i] <- mcmc$w[i] - mcmc$w[h_new] - 1
  mcmc <- update_genetics_upstream(mcmc, mcmc, i, h_new) # Update genetics. i is inheriting from h_new.
  mcmc$d[h_old] <- mcmc$d[h_old] - 1
  mcmc$d[h_new] <- mcmc$d[h_new] + 1
  return(mcmc)
}

# Update genetics for the following topological move:
# From g -> h -> i
# To g -> i, g -> h
update_genetics_downstream <- function(prop, mcmc, i, h){
  # Everything that doesn't stay the same in i
  all_i <- unique(c(
    mcmc$m01[[i]],
    mcmc$m10[[i]],
    mcmc$m0y[[i]],
    mcmc$m1y[[i]],
    mcmc$mx0[[i]],
    mcmc$mx1[[i]],
    mcmc$mxy[[i]]
  ))

  # Everything that doesn't stay the same in h
  all_h <- unique(c(
    mcmc$m01[[h]],
    mcmc$m10[[h]],
    mcmc$m0y[[h]],
    mcmc$m1y[[h]],
    mcmc$mx0[[h]],
    mcmc$mx1[[h]],
    mcmc$mxy[[h]]
  ))

  prop$m01[[i]] <- setdiff(mcmc$m01[[h]], all_i) # 01 in h, 11 in i
  prop$m01[[i]] <- union(prop$m01[[i]], intersect(mcmc$m0y[[h]], mcmc$mx1[[i]])) # 0y in h, x1 in i
  prop$m01[[i]] <- union(prop$m01[[i]], setdiff(mcmc$m01[[i]], all_h)) # 00 in h, 01 in i

  prop$m10[[i]] <- setdiff(mcmc$m10[[h]], all_i) # 10 in h, 00 in i
  prop$m10[[i]] <- union(prop$m10[[i]], intersect(mcmc$m1y[[h]], mcmc$mx0[[i]])) # 1y in h, x0 in i
  prop$m10[[i]] <- union(prop$m10[[i]], setdiff(mcmc$m10[[i]], all_h)) # 11 in h, 10 in i

  prop$m0y[[i]] <- setdiff(mcmc$m0y[[i]], all_h) # 00 in h, 0y in i
  prop$m0y[[i]] <- union(prop$m0y[[i]], intersect(mcmc$m10[[h]], mcmc$m1y[[i]])) # 01 in h, 1y in i
  prop$m0y[[i]] <- union(prop$m0y[[i]], intersect(mcmc$m0y[[h]], mcmc$mxy[[i]])) # 0y in h, xy in i

  prop$m1y[[i]] <- setdiff(mcmc$m1y[[i]], all_h) # 11 in h, 1y in i
  prop$m1y[[i]] <- union(prop$m1y[[i]], intersect(mcmc$m10[[h]], mcmc$m0y[[i]])) # 10 in h, 0y in i
  prop$m1y[[i]] <- union(prop$m1y[[i]], intersect(mcmc$m1y[[h]], mcmc$mxy[[i]])) # 1y in h, xy in i

  prop$mx0[[i]] <- intersect(mcmc$mxy[[h]], mcmc$mx0[[i]]) # xy in h, x0 in i
  prop$mx0[[i]] <- union(prop$mx0[[i]], setdiff(mcmc$mx0[[h]], all_i)) # x0 in h, 00 in i
  prop$mx0[[i]] <- union(prop$mx0[[i]], intersect(mcmc$mx1[[h]], mcmc$m10[[i]])) # x1 in h, 10 in i

  prop$mx1[[i]] <- intersect(mcmc$mxy[[h]], mcmc$mx1[[i]]) # xy in h, x1 in i
  prop$mx1[[i]] <- union(prop$mx1[[i]], setdiff(mcmc$mx1[[h]], all_i)) # x1 in h, 11 in i
  prop$mx1[[i]] <- union(prop$mx1[[i]], intersect(mcmc$mx0[[h]], mcmc$m01[[i]])) # x0 in h, 01 in i

  prop$mxy[[i]] <- intersect(mcmc$mxy[[h]], mcmc$mxy[[i]]) # xy in h, xy in i
  prop$mxy[[i]] <- union(prop$mxy[[i]], intersect(mcmc$mx1[[h]], mcmc$m1y[[i]])) # x1 in h, 1y in i
  prop$mxy[[i]] <- union(prop$mxy[[i]], intersect(mcmc$mx0[[h]], mcmc$m0y[[i]])) # x0 in h, 0y in i

  # Whew.
  return(prop)
}

# Wrap as a function: switch from
# h_new -> h_old -> i
# h_old -> i, h_new -> i
shift_downstream <- function(mcmc, i, h_old, h_new){
  # Update all necessary components of MCMC
  mcmc$h[i] <- h_new # Update the ancestor
  mcmc$w[i] <- mcmc$w[i] + mcmc$w[h_old] + 1
  mcmc <- update_genetics_downstream(mcmc, mcmc, i, h_old) # Update genetics. i is inheriting from h_new, but compared to genetics of h_old
  mcmc$d[h_old] <- mcmc$d[h_old] - 1
  mcmc$d[h_new] <- mcmc$d[h_new] + 1
  return(mcmc)
}

# Flip the genotype for a SNV
flip_genotype <- function(prop, mcmc, i, js, snv){
  ## Run through cases of updating genetic info in i
  if(snv %in% mcmc$m01[[i]]){

    # Delete from 01 in i
    prop$m01[[i]] <- setdiff(mcmc$m01[[i]], snv)

    # Note that we're changing from 1 to 0 in i
    add <- FALSE
  }else if(snv %in% mcmc$mx1[[i]]){

    # Delete from x1 in i
    prop$mx1[[i]] <- setdiff(mcmc$mx1[[i]], snv)

    # Union to x0 in i
    prop$mx0[[i]] <- union(mcmc$mx0[[i]], snv)

    # Note that we're changing from 1 to 0 in i
    add <- FALSE
  }else if(snv %in% mcmc$m10[[i]]){

    # Delete from 01 in i
    prop$m10[[i]] <- setdiff(mcmc$m10[[i]], snv)

    # Note that we're changing from 0 to 1 in i
    add <- TRUE
  }else if(snv %in% mcmc$mx0[[i]]){

    # Delete from x1 in i
    prop$mx0[[i]] <- setdiff(mcmc$mx0[[i]], snv)

    # Union to x0 in i
    prop$mx1[[i]] <- union(mcmc$mx1[[i]], snv)

    # Note that we're changing from 0 to 1 in i
    add <- TRUE
  }else if(snv %in% c(unlist(mcmc$m10[js]), unlist(mcmc$m1y[js]))){ # 11 in i

    # Union to 10 in i
    prop$m10[[i]] <- union(mcmc$m10[[i]], snv)

    add <- FALSE
  }else if(snv %in% c(unlist(mcmc$m01[js]), unlist(mcmc$m0y[js]))){ # 00 in i

    # Union to 01 in i
    prop$m01[[i]] <- union(mcmc$m01[[i]], snv)

    add <- TRUE
  }else{
    print("Something is wrong.")
  }

  ## Now update genetic info for j in js
  if(add){
    for (j in js) {
      if(snv %in% mcmc$m01[[j]]){
        # Delete from 01 in j
        prop$m01[[j]] <- setdiff(mcmc$m01[[j]], snv)
      }else if(snv %in% mcmc$m0y[[j]]){
        # Delete from 0y in j
        prop$m0y[[j]] <- setdiff(mcmc$m0y[[j]], snv)
        # Add to 1y in j
        prop$m1y[[j]] <- union(mcmc$m1y[[j]], snv)
      }else{
        # Otherwise, it was in 00 in j
        prop$m10[[j]] <- union(mcmc$m10[[j]], snv)
      }
    }
  }else{
    for (j in js) {
      if(snv %in% mcmc$m10[[j]]){
        # Delete from 10 in j
        prop$m10[[j]] <- setdiff(mcmc$m10[[j]], snv)
      }else if(snv %in% mcmc$m0y[[j]]){
        # Delete from 1y in j
        prop$m1y[[j]] <- setdiff(mcmc$m1y[[j]], snv)
        # Add to 0y in j
        prop$m0y[[j]] <- union(mcmc$m0y[[j]], snv)
      }else{
        # Otherwise, it was in 11 in j
        prop$m01[[j]] <- union(mcmc$m01[[j]], snv)
      }
    }
  }

  return(prop)
}

## Resample the genotype for an unobserved host, or for an observed host with missing sites, based on (approximate) parsimony
genotype <- function(mcmc, i, js, eps, comparison = F){
  # Get all SNVs that might need to change
  snvs <- unique(c(
    mcmc$mx0[[i]],
    mcmc$m10[[i]],
    mcmc$m01[[i]],
    mcmc$mx1[[i]],
    unlist(mcmc$m01[js]),
    unlist(mcmc$m0y[js]),
    unlist(mcmc$m10[js]),
    unlist(mcmc$m1y[js])
  ))

  # Which ones go 0y or 1y in j, or x0 or x1 in i? (with multiplicity)
  isnvs <- c(
    mcmc$mx0[[i]],
    mcmc$mx1[[i]],
    unlist(mcmc$m0y[js]),
    unlist(mcmc$m1y[js])
  )
  tab_isnv <- table(isnvs)
  ind_isnv <- match(names(tab_isnv), snvs) # Indices of these iSNVs in "snvs"

  # And the others?
  non_isnvs <- c(
    mcmc$m10[[i]],
    mcmc$m01[[i]],
    unlist(mcmc$m01[js]),
    unlist(mcmc$m10[js])
  )
  tab_non_isnv <- table(non_isnvs)
  ind_non_isnv <- match(names(tab_non_isnv), snvs) # Indices of these SNVs in "snvs"

  # Probability of swapping
  probs <- rep(0, length(snvs))
  probs[ind_non_isnv] <- probs[ind_non_isnv] + unname(tab_non_isnv)
  probs[ind_isnv] <- probs[ind_isnv] + unname(tab_isnv)/2
  probs <- probs / (length(js) + 1)

  if(comparison){
    # Anything with a probability of 0 or 1 won't be considered when creating a new node i
    delete <- which((probs == 0 | probs == 1))
    probs <- probs[-delete]
  }

  # By parsimony, round probability
  probs[probs < 0.5] <- 0
  probs[probs > 0.5] <- 1

  # Add random noise
  probs <- (eps/2) + (1 - eps)*probs

  # If our goal is to compute the probability that a new host i has this genotype...
  if(comparison){
    # The probability of creating a new genotype at i equal to that in MCMC is the probability we don't swap anything in i
    return(sum(log(1-probs)))
  }else{
    # Which ones get swapped?
    which_swap <- which(runif(length(snvs)) < probs)
    swaps <- snvs[which_swap]

    # Log probability of picking this genotype
    log_p <- sum(log(probs[which_swap])) + sum(log(1 - probs[-which_swap]))

    for (snv in swaps) {
      mcmc <- flip_genotype(mcmc, mcmc, i, js, snv)
    }

    return(list(mcmc, log_p))
  }
}

# Accept / reject
accept_or_reject <- function(prop, mcmc, data, update, hastings){
  prop$e_lik <- e_lik(prop, data)
  prop$g_lik[update] <- sapply(update, g_lik, mcmc = prop, data = data)
  prop$prior <- prior(prop)

  # Accept / reject
  if(log(runif(1)) < prop$e_lik + sum(prop$g_lik[-1]) + prop$prior - mcmc$e_lik - sum(mcmc$g_lik[-1]) - mcmc$prior + hastings){
    print("yooooooo")
    return(prop)
  }else{
    return(mcmc)
  }
}







###### MAYBE DELETING

## Create the genotype of a newly added host
create_genotype <- function(prop, mcmc, i, j1, j2){
  ## Create lists of SNVs in category A in j1, B in j2
  s_01_01 <- intersect(mcmc$m01[[j1]], mcmc$m01[[j2]])
  s_10_10 <- intersect(mcmc$m10[[j1]], mcmc$m10[[j2]])

  s_x0_x0 <- intersect(mcmc$mx0[[j1]], mcmc$mx0[[j2]])
  s_x0_xy <- intersect(mcmc$mx0[[j1]], mcmc$mxy[[j2]])
  s_x0_x1 <- intersect(mcmc$mx0[[j1]], mcmc$mx1[[j2]])

  s_xy_x0 <- intersect(mcmc$mxy[[j1]], mcmc$mx0[[j2]])
  s_xy_xy <- intersect(mcmc$mxy[[j1]], mcmc$mxy[[j2]])
  s_xy_x1 <- intersect(mcmc$mxy[[j1]], mcmc$mx1[[j2]])

  s_x1_x0 <- intersect(mcmc$mx1[[j1]], mcmc$mx0[[j2]])
  s_x1_xy <- intersect(mcmc$mx1[[j1]], mcmc$mxy[[j2]])
  s_x1_x1 <- intersect(mcmc$mx1[[j1]], mcmc$mx1[[j2]])

  prop$m01[[i]] <- s_01_01
  prop$m10[[i]] <- s_10_10


  # Which SNVs go x0, x1 in i?
  prop$mx0[[i]] <- unique(c(s_x0_x0, s_x0_xy, s_xy_x0))
  prop$mx1[[i]] <- unique(c(s_xy_x1, s_x1_xy, s_x1_x1))

  ambiguous <- unique(c(s_x0_x1, s_xy_xy, s_x1_x0))
  picks <- runif(length(ambiguous)) < 1/2

  prop$mx0[[i]] <- union(prop$mx0[[i]], ambiguous[picks])
  prop$mx1[[i]] <- union(prop$mx1[[i]], ambiguous[!picks])

  prop$m0y[[i]] <- character(0)
  prop$mxy[[i]] <- character(0)
  prop$m1y[[i]] <- character(0)

  # Now update changes in j1, j2

  prop <- update_genetics_upstream(prop, prop, j1, i)
  prop <- update_genetics_upstream(prop, prop, j2, i)

  return(list(prop, ambiguous))
}

## Create the genotype of a newly added host DOWNSTREAM
# Recall this means:
# h -> j1 -> j2 goes to
# h-> i -> j1, j2

create_genotype_downstream <- function(prop, mcmc, i, j1, j2){

  # Everything that doesn't stay the same in h
  all_j2 <- unique(c(
    mcmc$m01[[j2]],
    mcmc$m10[[j2]],
    mcmc$m0y[[j2]],
    mcmc$m1y[[j2]],
    mcmc$mx0[[j2]],
    mcmc$mx1[[j2]],
    mcmc$mxy[[j2]]
  ))

  ## Create lists of SNVs in category A in j1, B in j2
  s_01_11 <- setdiff(mcmc$m01[[j1]], all_j2) # 01 in j1, 11 in j2
  s_10_00 <- setdiff(mcmc$m10[[j1]], all_j2) # 10 in j1, 00 in j2

  s_x0_00 <- setdiff(mcmc$mx0[[j1]], all_j2)
  s_x0_0y <- intersect(mcmc$mx0[[j1]], mcmc$m0y[[j2]])
  s_x0_01 <- intersect(mcmc$mx0[[j1]], mcmc$m01[[j2]])

  s_xy_x0 <- intersect(mcmc$mxy[[j1]], mcmc$mx0[[j2]])
  s_xy_xy <- intersect(mcmc$mxy[[j1]], mcmc$mxy[[j2]])
  s_xy_x1 <- intersect(mcmc$mxy[[j1]], mcmc$mx1[[j2]])

  s_x1_10 <- intersect(mcmc$mx1[[j1]], mcmc$m10[[j2]])
  s_x1_1y <- intersect(mcmc$mx1[[j1]], mcmc$m1y[[j2]])
  s_x1_11 <- setdiff(mcmc$mx1[[j1]],all_j2)

  prop$m01[[i]] <- s_01_11
  prop$m10[[i]] <- s_10_00


  # Which SNVs go x0, x1 in i?
  prop$mx0[[i]] <- unique(c(s_x0_00, s_x0_0y, s_xy_x0))
  prop$mx1[[i]] <- unique(c(s_xy_x1, s_x1_1y, s_x1_11))

  ambiguous <- unique(c(s_x0_01, s_xy_xy, s_x1_10))
  picks <- runif(length(ambiguous)) < 1/2

  prop$mx0[[i]] <- union(prop$mx0[[i]], ambiguous[picks])
  prop$mx1[[i]] <- union(prop$mx1[[i]], ambiguous[!picks])

  prop$m0y[[i]] <- character(0)
  prop$mxy[[i]] <- character(0)
  prop$m1y[[i]] <- character(0)

  # Now update changes in j1, j2

  prop <- update_genetics_upstream(prop, prop, j1, i) # j1 moves upstream from h onto i
  prop <- update_genetics_downstream(prop, prop, j2, j1) # j2 moves downstream from j1 onto i, given updates to j1

  return(list(prop, ambiguous))
}




# ## Varilly coalescent function
# cppFunction('double varilly(double r, double p, int d, int D, int N) {
#   double term = 1;
#   for(int k=1; k<=d; ++k){
#     term = term * (r - 1 + k) * (1-p) / (N - k + 1);
#   }
#   double sum = term;
#   for(int k=d+1; k<=N-D+d; ++k){
#     term = term * (r - 1 + k) * (1-p) * (N - D - k + d + 1) / (N - k + 1) / (k - d);
#     sum += term;
#   }
#   return sum;
# }')
#
# microbenchmark({
#   varilly(0.1, 0.01, 4, 10, 10000)
# })
#
# microbenchmark({
#   lgamma(10000)
# })




# cppFunction('double h2f1(double r, double p, int d, int D, int N) {
#   double p1 = 1;
#   int p2 = 1;
#   int p3 = 1;
#   int p4 = 1;
#   double p5 = 1;
#   double sum = 1;
#   for(int n=1; n<=N-D; ++n){
#     p1 *= d + r + n - 1;
#     p2 *= N - d - n + 1;
#     p3 *= N - D - n + 1;
#     p4 *= n;
#     p5 *= 1 - p;
#     sum += p3 * p1 * p5 / p4 / p2;
#   }
#   return sum;
# }')


# cppFunction('double h2f1(double r, double p, int d, int D, int N) {
#   double term = 1;
#   double sum = 1;
#   for(int n=1; n<=N-D; ++n){
#     term *= 1 - p;
#     term = term * (N - D - n + 1) / (N - d - n + 1);
#     term = term * (d + r + n - 1) / n;
#     sum += term;
#   }
#   return sum;
# }')


# prop$m10[[j1]] <- setdiff(mcmc$m10[[j1]], s_10_10)
# prop$m10[[j2]] <- setdiff(mcmc$m10[[j2]], s_10_10)

# prop$m01[[j1]] <- setdiff(mcmc$m01[[j1]], prop$m01[[i]])
# prop$m01[[j2]] <- setdiff(mcmc$m01[[j2]], prop$m01[[i]])

# for (j in c(j1, j2)) {
#   prop$mx0[[j]] <- character(0)
#   prop$mxy[[j]] <- character(0)
#   prop$mx1[[j]] <- character(0)
#
#   prop$m0y[[j]] <- union(prop$m0y[[j]], intersect(prop$mx0[[i]], mcmc$mxy[[j]])) # x0 in i, formerly xy in j
#   prop$m01[[j]] <- union(prop$m01[[j]], intersect(prop$mx0[[i]], mcmc$mx1[[j]])) # x0 in i, formerly x1 in j
#
#   prop$m10[[j]] <- union(prop$m10[[j]], intersect(prop$mx1[[i]], mcmc$mx0[[j]])) # x1 in i, formerly x0 in j
#   prop$m1y[[j]] <- union(prop$m1y[[j]], intersect(prop$mx1[[i]], mcmc$mxy[[j]])) # x1 in i, formerly xy in j
# }



