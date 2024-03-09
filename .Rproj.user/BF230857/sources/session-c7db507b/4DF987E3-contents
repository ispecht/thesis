### Compare empirical and approximate theoretical distrubutions for the RPW model
library(parallel)
library(extraDistr)
library(ggplot2)
library(tikzDevice)
library(cowplot)

## Simulate the RPW model, starting with 0 white balls and 1 black ball
# Assume p = p_W = 3p_B
# An efficient way of simulating this: select the steps at which the mutations happen as Geom(3p) draws
# At these steps, choose a ball. If white, return black. If black, return white with probability 1/3
rpw <- function(x, p, n){
  # Number of white balls
  w <- 0
  # Number of total balls
  tot <- 1
  while (tot < n) {
    # How many balls added to the urn before the transition?
    before <- rgeom(1, p)
    if(tot + before >= n){
      before <- n - tot
      w <- w + rbinom(1, before, rbeta(1, w, tot-w))
      break
    }
    if(w!=0){
      w <- w + rbinom(1, before, rbeta(1, w, tot-w))
    }
    tot <- tot + before
    # Did we draw a white ball?
    if(runif(1) < w/tot){
      # Add a black ball
      tot <- tot + 1
    }else{
      # Add a white ball with probability p_B = p/3
      if(runif(1) < 1/3){
        w <- w+1
      }
      tot <- tot + 1
    }
  }
  print(w)
  return(w/n)
}

## Theoretical CDF
t_cdf <- function(x, p, k){
  if(x == 0){
    (1-p/3)^k
  }else{
    # Integrate theoretical PDF
    (1-p/3)^k + (1-(1-p/3)^k)*(-1+(1-x)^(1+k)+x+k *x)/(k *x)
  }
}




# How many model runs?
n_iters <- 1e6

# Function to plot empirical and theoretical density
get_plot <- function(n = 1e6, p = 1e-6){

  # Set bin width
  bw <- sqrt(p)

  # Set max plot range
  max_range = 50*bw

  # Bin right edges
  bins <- seq(bw, max_range, bw)

  # Simulate RPW model
  sims <- unlist(mclapply(1:n_iters, rpw, p=p, n =n, mc.cores = 12))

  # Get proportion of iterations in each bin (bins being right edges)
  e_cdf <- unlist(mclapply(bins, function(x){mean(sims < x)}, mc.cores = 12))

  # Theoretical CDF
  # Take k = sqrt(1/p)
  t_vals <- sapply(bins, t_cdf, p=p, k = sqrt(1/p))

  df <- data.frame(x = rep(bins, 2), y = c(t_vals, e_cdf), CDF = c(rep("Theoretical", length(bins)), rep("Empirical", length(bins))))

  p <- ggplot(df, aes(x = x, y=y, group = CDF)) +
    geom_line(aes(color = CDF)) +
    ylab("Cumulative Density") +
    xlab("Proportion of White Balls") +
    scale_color_manual(values = c("black", "red")) +
    theme_minimal() +
    theme(legend.position = "top")

  return(p)
}

ns <- c(1e5, 1e6, 1e7)
ps <- c(1e-7, 1e-6, 1e-5)

plots <- list()
for (i in 1:length(ps)) {
  plots <- c(plots, list(get_plot(p = ps[i])))
  print(i)
}
for (i in 1:length(ns)) {
  plots <- c(plots, list(get_plot(n = ns[i])))
  print(i)
}

all_plots <- plot_grid(plots[[1]], plots[[4]], plots[[2]], plots[[5]], plots[[3]], plots[[6]], ncol = 2, labels = c("A","D","B","E","C","F"))

ggsave("./figs/rpw.pdf", width = 6.5, height = 7.5)
ggsave("./figs/rpw.png", width = 6.5, height = 7.5)

### Real data application

# Minimum alternate allele frequency
min_af_filter <- 0.03

# Maximum alternate allele frequency
max_af_filter <- 0.5

load("./maf.RData")

### Histogram and fitted density function

# PDF, CDF, and quantile funtions, using approximation in Methods
approx_pdf <- function(x){
  max_af_filter * min_af_filter / ((max_af_filter - min_af_filter) * x^2)
}

approx_cdf <- function(x){
  (max_af_filter * (min_af_filter - x))/((min_af_filter - max_af_filter) * x)
}

approx_qtile <- function(p){
  min_af_filter * max_af_filter / (max_af_filter + p*(min_af_filter - max_af_filter))
}

df <- data.frame(log(maf))

p1 <- ggplot(df, aes(x = (maf))) +
  geom_histogram(aes(y=after_stat(density)), binwidth = 0.01, boundary = 0.01, color = "white", fill = "grey") +
  geom_vline(xintercept=min_af_filter + 0.00001, linetype=3, color = "blue") +
  geom_function(fun = approx_pdf, color = "red", linewidth = 1, linetype = "dashed") +
  xlab("Minor Allele Frequency") +
  ylab("Probability Density") +
  theme_minimal()

#print(p1)

ggsave("./figs/hist.png", width = 6, height = 6)
ggsave("./figs/hist.pdf", width = 6, height = 6)

### Histogram with log-transform

round_label <- function(x){
  round(x, 2)
}

# Empirical log density
bw <- 0.01
rights <- seq(min_af_filter + bw, max_af_filter, bw)
lefts <- rights - bw
counts <- c()
for (i in 1:length(rights)) {
  counts[i] <- sum(maf < rights[i] & maf >= lefts[i])
}
counts <- counts/sum(counts)/bw

p2 <- ggplot(data.frame(x = lefts), aes(x=x)) +
  geom_line(aes(y = counts), stat = 'identity', linewidth = 1, color = "black") +
  geom_vline(xintercept=min_af_filter + 0.00001, linetype=3, color = "blue") +
  geom_function(fun = approx_pdf, color = "red", linewidth = 1, linetype = "dashed") +
  xlab("Log Minor Allele Frequency") +
  scale_x_continuous(trans = 'log', labels = round_label) +
  scale_y_continuous(trans = 'log', labels = round_label) +
  ylab("Log Probability Density") +
  theme_minimal()

#print(p2)

ggsave("./figs/loglog.png", width = 6, height = 6)
ggsave("./figs/loglog.pdf", width = 6, height = 6)



### Q-Q Plot

breaks <- approx_cdf(
  seq(min_af_filter, max_af_filter, 0.01)
)

# Rewrite CDF to correctly compute values below min_af_filter and above max_af_filter

full_cdf <- function(x){
  below <- which(x <= min_af_filter)
  above <- which(x  > max_af_filter)
  other <- which(x > min_af_filter & x <= max_af_filter)

  x[below] <- 0
  x[above] <- 1
  x[other] <- approx_cdf(x[other])

  x
}


ks.test(maf + rnorm(length(maf), 0, 0.000000001), full_cdf)


### CDF comparison
p3 <- ggplot(data.frame(x=maf), aes(x=x)) +
  stat_ecdf(geom = "line", color = "black", linewidth = 1, pad = F) +
  geom_function(fun = full_cdf, color = "red", linewidth = 1, linetype = "dashed", xlim = c(min_af_filter, max_af_filter)) +
  xlab("Minor Allele Frequency") +
  ylab("Cumulative Probability Density") +
  # xlim(c(0,0.01)) +
  # ylim(c(0,0.01)) +
  # scale_x_continuous(trans='log') +
  # scale_y_continuous(trans='log') +
  theme_minimal()

#print(p3)

ggsave("./figs/cdf.png", width = 6, height = 6)
ggsave("./figs/cdf.pdf", width = 6, height = 6)


fig2 <- plot_grid(p1, p3, labels = "AUTO", ncol = 2)
ggsave("./figs/fig2.pdf", width = 6.5, height = 3)
ggsave("./figs/fig2.png", width = 6.5, height = 3)





