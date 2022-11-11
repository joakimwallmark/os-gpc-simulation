library(gss)
# for empirical testgardner theta densities
tg_theta_density <- function(tg_thetas, alpha = 1.4) {
  n <- length(tg_thetas)
  theta_dens <- tg_thetas[0 < tg_thetas & tg_thetas < 100]
  dens <- ssden(~ theta_dens, alpha = alpha, domain = data.frame(theta_dens = c(0, 100)))
  p_max <- sum(tg_thetas == 100)/n
  p_min <- sum(tg_thetas == 0)/n
  
  df <- data.frame(x = seq(0, 100, by = 0.1), y = (1-p_max-p_min)*sapply(seq(0, 100, by = 0.1), function(x) { dssden(dens, x) } ))
  plot_tg <- ggplot(data.frame(tg_thetas), aes(x = tg_thetas)) +
    geom_histogram(aes(y = ..density..), colour = "black", fill = "white", 
                   bins = 14, boundary = 0) +
    geom_line(data = df, mapping = aes(x, y)) +
    geom_point(data = data.frame(x=c(0, 100), y=c(p_min, p_max)), mapping = aes(x, y), size = 3) +
    theme_bw() +
    xlab(expression(theta)) +
    ylab(expression(f(theta))) +
    scale_x_continuous(limits = c(0, 100), expand = c(0.02, 0.02)) +
    scale_y_continuous(expand = c(0, 0, 0.05, 0))
  return(list(density = list(p_min = p_min, p_max = p_max, dens_obj = dens), plot = plot_tg))
}

# get quantile for given p vector
tg_qss <- function(p, dens_list) {
  q <- vector("numeric", length = length(p))
  q[p < dens_list$p_min] <- 0
  q[p > (1-dens_list$p_max)] <- 100
  den_ind <- dens_list$p_min < p & p < (1-dens_list$p_max)
  if (sum(den_ind)) {
    p_dens <- p[den_ind]
    p_dens <- (p_dens-dens_list$p_min)/(1-dens_list$p_min-dens_list$p_max)
    q[den_ind] <- qssden(dens_list$dens_obj, p_dens)
  }
  return(q)
}


# get probabilities/densities for theta<=x for vector x
tg_pss <- function(x, dens_list) {
  p <- vector("numeric", length = length(x))
  p[x == 0] <- dens_list$p_min
  p[x == 100] <- 1
  den_ind <- 0 < x & x < 100
  if (sum(den_ind)) {
    p[den_ind] <- dens_list$p_min + (1-dens_list$p_min-dens_list$p_max)*pssden(dens_list$dens_obj, x[den_ind])
  }
  return(p)
}

# get probability/density for a given x
tg_dss <- function(x, dens_list) {
  d <- vector("numeric", length = length(x))
  d[x == 0] <- dens_list$p_min
  d[x == 100] <- dens_list$p_max
  den_ind <- 0 < x & x < 100
  if (sum(den_ind)) {
    d[den_ind] <- (1 - dens_list$p_min - dens_list$p_max)*dssden(dens_list$dens_obj, x[den_ind])
  }
  return(d)
}

# generate sample from tg theta distribution
tg_rss <- function(n, dens_list) {
  theta <- sample(
    c(0, 1, 100),
    n,
    replace = T,
    probs = c(
      dens_list$p_min,
      1 - dens_list$p_min - dens_list$p_max,
      dens_list$p_max
    )
  )
  theta[theta == 1] <- qssden(dens_list$dens_obj, runif(sum(theta == 1)))
  return(theta)
}