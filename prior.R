### Prior function and defaults

prior <- function(mcmc, priors = NULL){
  if(is.null(priors)){
    priors <- list()
  }
  if(is.null(priors$b)){
    priors$b <- function(x){
      # Unif(0,1)
      if(0 <= x & x <= 1){
        return(0)
      }else{
        return(-Inf)
      }
    }
  }

  if(is.null(priors$a_g)){
    priors$a_g <- function(x){
      # Unif(0,100)
      # if(0 <= x & x <= 100){
      #   return(-log(100))
      # }else{
      #   return(-Inf)
      # }
      return(dnorm(x, 5, 0.5, log = T))
    }
  }

  if(is.null(priors$a_s)){
    priors$a_s <- function(x){
      # Unif(0,100)
      if(0 <= x & x <= 100){
        return(-log(100))
      }else{
        return(-Inf)
      }
    }
  }

  if(is.null(priors$mu)){
    priors$mu <- function(x){
      # Unif(0,1)
      if(0 <= x & x <= 1){
        return(0)
      }else{
        return(-Inf)
      }
    }
  }

  if(is.null(priors$p)){
    priors$p <- function(x){
      # Unif(0,1)
      if(0 <= x & x <= 1){
        return(0)
      }else{
        return(-Inf)
      }
    }
  }

  if(is.null(priors$v)){
    priors$v <- function(x){
      # Unif(1,1e4)
      if(0 <= x & x <= 1e4){
        return(-log(1e4))
      }else{
        return(-Inf)
      }
    }
  }

  if(is.null(priors$lambda)){
    priors$lambda <- function(x){
      dnorm(x, 2, 1, log = T)
    }
  }

  if(is.null(priors$rho)){
    priors$rho <- function(x){
      # Unif(0,100)
      if(0 <= x & x <= 100){
        return(-log(100))
      }else{
        return(-Inf)
      }
    }
  }

  if(is.null(priors$psi)){
    priors$psi <- function(x){
      # Unif(0,100)
      if(0 <= x & x <= 100){
        return(-log(100))
      }else{
        return(-Inf)
      }
    }
  }

  priors$b(mcmc$b) +
    priors$a_g(mcmc$a_g) +
    priors$a_s(mcmc$a_s) +
    priors$mu(mcmc$mu) +
    priors$p(mcmc$p) +
    priors$v(mcmc$v) +
    priors$lambda(mcmc$lambda) +
    priors$rho(mcmc$rho) +
    priors$psi(mcmc$psi)

}






