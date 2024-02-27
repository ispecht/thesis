# Recursive form of epi likelihood
coalescent <- function(mcmc, data){
  gens <- gen(mcmc$h, mcmc$w)

  G <- max(gens)
  qs <- get_qs(mcmc$n - 1 + sum(mcmc$w), data$N, G, mcmc$rho, mcmc$psi)

  out <- 0
  for (i in 1:mcmc$n) {
    out <- out + e_lik_node(mcmc$d[i], mcmc$w[i], qs[(gens[i] + 1 - mcmc$w[i]):(gens[i] + 1)], mcmc$rho, mcmc$psi)
  }

  out <- out + lchoose(data$N, mcmc$n - data$n_obs) + lfactorial(mcmc$n - data$n_obs) - (mcmc$n - 1)*log(data$N)
  return(out)
}


# Compute epidemiological log likelihood

e_lik <- function(mcmc, data){

  if(
    any(mcmc$w < 0) |
    mcmc$a_g < 0 |
    mcmc$lambda_g < 0 |
    mcmc$a_s < 0 |
    mcmc$lambda_s < 0 |
    mcmc$rho < 0 |
    mcmc$psi < 0
  ){
    return(-Inf)
  }else{
    return(
      # Generation intervals
      sum(dgamma(mcmc$t[2:mcmc$n] - mcmc$t[mcmc$h[2:mcmc$n]], shape = (mcmc$w[2:mcmc$n] + 1) * mcmc$a_g, rate = mcmc$lambda_g, log = T)) +

        # Sojourn intervals
        sum(dgamma(data$s[2:data$n_obs] - mcmc$t[2:data$n_obs], shape = mcmc$a_s, rate = mcmc$lambda_s, log = T)) +

        # Varilly Coalescent
        # sum(mcmc$d * log((1-mcmc$psi) / mcmc$psi) + lchoose(mcmc$d + mcmc$rho - 1, mcmc$d) - lchoose(data$N, mcmc$d)) +
        # sum(mcmc$w) * (log(mcmc$rho) + log((1-mcmc$psi) / mcmc$psi) - log(data$N))
        coalescent(mcmc, data)

    )
  }
}

# Compute genomic log likelihood for each person

g_lik <- function(mcmc, data, i){

  if(mcmc$v < 0 | mcmc$mu < 0 | mcmc$p < 0 | mcmc$b < 0 | mcmc$b > 1 | mcmc$w[i] < 0){
    -Inf
  }else{
    # Time of end of expo growth phase for ancestor of i, which is when we reach k = 1/sqrt(p) virions
    #g <- mcmc$t[mcmc$h[i]] - (mcmc$p * log(mcmc$p) * (-1 + 0.577215664901533 + log(mcmc$v)))/(2*mcmc$mu * log(mcmc$v))

    #g <- mcmc$t[mcmc$h[i]] - log(mcmc$p) / (2*log(mcmc$mu / mcmc$p))
    #g <- mcmc$t[mcmc$h[i]] + (mcmc$p *(-(log(mcmc$p)/2) + log(mcmc$v)))/(mcmc$mu* log(mcmc$v))
    #mcmc$p <- mcmc$mu / mcmc$lambda

    g <- mcmc$t[mcmc$h[i]] + log(1/sqrt(mcmc$p)) / mcmc$lambda

    # Evolutionary time
    delta_t <- mcmc$t[i] - g

    if(i > data$n_obs){
      delta_t <- mcmc$t[i] + log(1/sqrt(mcmc$p)) / mcmc$lambda - mcmc$t[mcmc$h[i]]
    }

    # Evolutionary time from first downstream host of h to infection time of i, approx
    delta_t_prime <- delta_t * mcmc$w[i] / (mcmc$w[i] + 1)

    if(delta_t <= 0){
      -Inf
    }else{
      # log probability of SPECIFIC iSNV in expo growth phase
      log_p_isnv <- log(
        (1 - (1 - mcmc$p)^(1/sqrt(mcmc$p))) * (1 - denovo_cdf(data$filters$af, mcmc$p)) / 3
      )

      # log probability of no iSNV in expo growth phase
      log_p_no_isnv <- log(
        (1 - mcmc$p)^(1/sqrt(mcmc$p)) +
          (1 - (1 - mcmc$p)^(1/sqrt(mcmc$p))) * denovo_cdf(data$filters$af, mcmc$p)
      )

      # Number of sites without a mutation
      no_mut <- data$n_bases -
        length(mcmc$m01[[i]]) -
        length(mcmc$m10[[i]]) -
        length(mcmc$m0y[[i]]) -
        length(mcmc$m1y[[i]]) -
        length(mcmc$mx0[[i]]) -
        length(mcmc$mx1[[i]]) -
        length(mcmc$mxy[[i]])

      # Likelihood from x = 0, y = 0 or x = 1, y = 1
      out <- no_mut * (log(1/4 + (3/4)*exp(-(4*mcmc$mu/3) * delta_t)) + ifelse(i<= data$n_obs, log_p_no_isnv, 0)) +

        # Likelihood from x = 0, y = 1 and x = 1, y = 0
        (length(mcmc$m01[[i]]) + length(mcmc$m10[[i]])) * (log(1/4 - (1/4)*exp(-(4*mcmc$mu/3) * delta_t)) + ifelse(i<= data$n_obs, log_p_no_isnv, 0))

      if(i <= data$n_obs){
        if(!is.null(data$snvs[[i]]$isnv)){

          # Frequencies of added iSNVs
          freq_0y <- data$snvs[[i]]$isnv$af[match(mcmc$m0y[[i]], data$snvs[[i]]$isnv$call)]

          # Frequencies of deleted iSNVs
          freq_1y <- data$snvs[[i]]$isnv$af[match(mcmc$m1y[[i]], data$snvs[[i]]$isnv$call)]

          out <- out +
            # Likelihood from x = 0, 0 < y < 1
            length(mcmc$m0y[[i]]) * log_p_isnv +
            sum(log(
              (1/4 + (3/4)*exp(-(4*mcmc$mu/3) * delta_t))*denovo_normed(freq_0y, mcmc$p, data$filters) + (1/4 - (1/4)*exp(-(4*mcmc$mu/3) * delta_t))*denovo_normed(1 - freq_0y, mcmc$p, data$filters)
            )) +

            # Likelihood from x = 1, 0 < y < 1
            length(mcmc$m1y[[i]]) * log_p_isnv +
            sum(log(
              (1/4 + (3/4)*exp(-(4*mcmc$mu/3) * delta_t))*denovo_normed(1 - freq_1y, mcmc$p, data$filters) + (1/4 - (1/4)*exp(-(4*mcmc$mu/3) * delta_t))*denovo_normed(freq_1y, mcmc$p, data$filters)
            ))

          if(mcmc$h[i] <= data$n_obs){
            if(!is.null(data$snvs[[mcmc$h[i]]]$isnv)){

              # Frequency of shared iSNV in ancestor of case i
              freq_xy_anc <- data$snvs[[mcmc$h[i]]]$isnv$af[match(mcmc$mxy[[i]], data$snvs[[mcmc$h[i]]]$isnv$call)]

              # Frequency of shared iSNV in case i
              freq_xy <- data$snvs[[i]]$isnv$af[match(mcmc$mxy[[i]], data$snvs[[i]]$isnv$call)]

              out <- out +
                # Likelihood from 0 < x < 1, 0 < y < 1
                length(mcmc$m1y[[i]]) * log_p_isnv +
                sum(log(
                  ((1/4 + (3/4)*exp(-(4*mcmc$mu/3) * delta_t_prime))*(1 - freq_xy_anc) + (1/4 - (1/4)*exp(-(4*mcmc$mu/3) * delta_t_prime))*freq_xy_anc) * denovo_normed(freq_xy, mcmc$p, data$filters) * (1 - (mcmc$b^(mcmc$w[i] + 1) * 2 * freq_xy_anc * (1 - freq_xy_anc) / 3^mcmc$w[i])) +
                    ((1/4 + (3/4)*exp(-(4*mcmc$mu/3) * delta_t_prime))*freq_xy_anc + (1/4 - (1/4)*exp(-(4*mcmc$mu/3) * delta_t_prime))*(1 - freq_xy_anc)) * denovo_normed(1 - freq_xy, mcmc$p, data$filters) * (1 - (mcmc$b^(mcmc$w[i] + 1) * 2 * freq_xy_anc * (1 - freq_xy_anc) / 3^mcmc$w[i])) +
                    (mcmc$b^(mcmc$w[i] + 1) * 2 * freq_xy_anc * (1 - freq_xy_anc) / 3^mcmc$w[i])
                ))
            }
          }
        }
      }

      if(mcmc$h[i] <= data$n_obs){
        if(!is.null(data$snvs[[mcmc$h[i]]]$isnv)){

          # Frequency of iSNV in ancestor of case i
          freq_x0_anc <- data$snvs[[mcmc$h[i]]]$isnv$af[match(mcmc$mx0[[i]], data$snvs[[mcmc$h[i]]]$isnv$call)]
          freq_x1_anc <- data$snvs[[mcmc$h[i]]]$isnv$af[match(mcmc$mx1[[i]], data$snvs[[mcmc$h[i]]]$isnv$call)]

          out <- out +
            # Likelihood from 0 < x < 1, y = 0
            length(mcmc$mx0[[i]]) * ifelse(i<= data$n_obs, log_p_no_isnv, 0) +
            sum(log(
              (1/4 + (3/4)*exp(-(4*mcmc$mu/3) * delta_t_prime))*(1 - freq_x0_anc) + (1/4 - (1/4)*exp(-(4*mcmc$mu/3) * delta_t_prime))*freq_x0_anc
            )) + sum(log(1 - (mcmc$b^(mcmc$w[i] + 1) * 2 * freq_x0_anc * (1 - freq_x0_anc) / 3^mcmc$w[i]))) + # probability we don't transmit successive split bottlenecks

            # Likelihood from 0 < x < 1, y = 1
            length(mcmc$mx1[[i]]) * ifelse(i<= data$n_obs, log_p_no_isnv, 0) +
            sum(log(
              (1/4 + (3/4)*exp(-(4*mcmc$mu/3) * delta_t_prime))*(freq_x1_anc) + (1/4 - (1/4)*exp(-(4*mcmc$mu/3) * delta_t_prime))*(1 - freq_x1_anc)
            )) + sum(log(1 - (mcmc$b^(mcmc$w[i] + 1) * 2 * freq_x1_anc * (1 - freq_x1_anc) / 3^mcmc$w[i]))) # probability we don't transmit successive split bottlenecks

        }
      }

      return(out)
    }
  }
}









