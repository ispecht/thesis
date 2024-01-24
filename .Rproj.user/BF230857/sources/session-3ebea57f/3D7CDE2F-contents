
e_data <- data.frame()

for (i in 900:920) {
  p <- plot_current(res[[i]]$h, 155)
  es <- p[[1]][[1]]
  vert <- p[[1]][[2]]

  nd <- data.frame(frame = i, x = es$x, y = es$y, group = es$group)

  e_data <- rbind(e_data, nd)
}

fig <- e_data %>%
  plot_ly(
    x = ~x,
    y = ~y,
    frame = ~frame,
    type = 'scatter',
    mode = 'lines',
    transforms = list(
      list(
        type = 'groupby',
        groups = e_data$group
      )
    ),
    showlegend = F
  )

fig

ggplot(data= e_data)+
  geom_path(aes(x=x,y=y,group = group)) +
  transition_states(frame, transition_length = 2,
                    state_length = 1)

ggplot(data= vert)+
  geom_point(aes(x=x,y=y,group = group), size = 0.1) +
  coord_fixed()


# Chop tree into roughly equally sized pieces
chop_old <- function(mcmc, data){
  n <- length(mcmc$h)
  out <- list()
  roots <- c()
  deg <- total_degree(mcmc$h, mcmc$d)

  if(data$n_subtrees > 1){
    for (i in 1:(data$n_subtrees - 1)) {
      # On average, a branch should have n/data$n_subtrees nodes

      #prob <- log(deg)

      #prob <- exp(-100*(prob - log(n/data$n_subtrees))^2)

      prob <- rep(0,n)
      prob[deg > 0.5*(n - length(unlist(out)))/(data$n_subtrees - i + 1) & deg < 1.5*(n - length(unlist(out)))/(data$n_subtrees - i + 1)] <- 1 # Everything close to n/data$n_subtrees gets probability 1

      # Not allowed to choose "1", else the algo will end too early
      prob[1] <- 0

      # Also, something that makes our lives much easier: only allow observed hosts with observed parents to be roots
      if(n > data$n_obs){
        prob[(data$n_obs):n] <- 0
        prob[which(mcmc$h > data$n_obs)] <- 0
      }


      prob <- 0
      # If everything is 0, just pick the one closest to the target value
      if(all(prob == 0)){
        dists <- abs((n - length(unlist(out)))/(data$n_subtrees - i + 1) - deg)
        dists[1] <- Inf
        dists[deg == 0] <- Inf
        if(n > data$n_obs){
          dists[(data$n_obs):n] <- Inf
          dists[which(mcmc$h > data$n_obs)] <- Inf
        }

        dists <- deg - n/data$n_subtrees/2
        dists[dists < 0] <- Inf

        pick <- which.min(dists)


      }else{
        pick <- sample(1:n, 1, prob = prob)
      }

      anc <- ancestry(mcmc$h, pick)
      deg[anc] <- deg[anc] - deg[pick]

      ups <- get_upstream(mcmc$h, pick)
      ups <- setdiff(ups, unlist(out))

      deg[ups] <- 0
      out[[i]] <- sort(ups)


      # Set the ancestor of "ups" to a dummy, to save compute time on "get_upstream"
      mcmc$h[ups] <- 0
      roots[i] <- pick
      #print(i)
    }
  }


  out[[data$n_subtrees]] <- sort(setdiff(2:n, unlist(out)))
  roots[data$n_subtrees] <- 1

  return(list(roots, out))
}

