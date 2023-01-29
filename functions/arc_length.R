source("functions/tg_density.R")
library(pracma)
library(TestGardener)
library(mirt)

# stolen from mirt and rewritten to remove hessian and support surprisal deriv.
deriv_mirt_gpc <- function(x, theta, surp = TRUE) {
  theta <- as.matrix(theta)
  a <- mirt:::ExtractLambdas(x)
  d <- mirt:::ExtractZetas(x)
  ak <- 0:(x@ncat - 1L)
  p <- mirt:::P.nominal(c(a, ak, d), ncat = length(d), Theta = theta)
  num <- mirt:::P.nominal(c(a, ak, d),
    ncat = length(d), Theta = theta,
    returnNum = TRUE
  )
  den <- rowSums(num)
  grad <- vector("list", x@ncat)
  for (i in seq_len(x@ncat)) grad[[i]] <- matrix(0, nrow(theta), x@nfact)
  for (j in seq_len(x@nfact)) { # j are factors, i are item categories
    for (i in seq_len(x@ncat)) {
      grad[[i]][, j] <- ak[i] * a[j] * p[, i] - p[, i] *
        (num %*% (ak * a[j])) / den
      if (surp) {
        grad[[i]][, j] <- -grad[[i]][, j] / (p[, i] * log(2))
      }
    }
  }
  return(do.call(cbind, grad))
}

arc_length_surp <- function(mod, thetas, total_cat, max_surp = NULL, item_ids = NULL, mesh_accuracy = 1001) {
  derivs <- matrix(0, mesh_accuracy, total_cat)
  m2 <- 0
  if (class(mod) == "SingleGroupClass") { # check if mirt model
    if (is.null(item_ids)) item_ids <- 1:extract.mirt(mod, "nitems")
    # We want item category derivatives in columns and theta values in rows
    theta_mesh <- seq(min(thetas), max(thetas), length.out = mesh_accuracy)
    for (i in item_ids) {
      item <- extract.item(mod, i)
      m1 <- m2 + 1
      m2 <- m2 + item@ncat
      deriv_resc <- (max(theta_mesh) - min(theta_mesh)) / (length(theta_mesh) - 1)
      derivs[, m1:m2] <- deriv_mirt_gpc(item, theta_mesh, surp = TRUE) * deriv_resc
      if (!is.null(max_surp)) {
        capped_surp <- -log2(probtrace(item, theta_mesh)) > max_surp
        derivs[, m1:m2] <- (1 - capped_surp) * derivs[, m1:m2]
      }
    }
  } else { # assume TestGardener
    w_list <- mod$parList[[length(mod$parList)]]$WfdList
    if (is.null(item_ids)) item_ids <- seq_along(w_list)
    theta_mesh <- seq(0, 100, length.out = mesh_accuracy)
    for (i in item_ids) {
      ncat <- w_list[[i]]$M
      m1 <- m2 + 1
      m2 <- m2 + ncat
      deriv_resc <- (max(theta_mesh) - min(theta_mesh)) / (length(theta_mesh) - 1)
      derivs[, m1:m2] <- deriv_resc * eval.surp(theta_mesh, w_list[[i]]$Wfd, nderiv = 1) / log(2, ncat)
      if (!is.null(max_surp)) {
        capped_surp <- eval.surp(theta_mesh, w_list[[i]]$Wfd) / log(2, ncat) > max_surp
        derivs[, m1:m2] <- (1 - capped_surp) * derivs[, m1:m2]
      }
    }
  }
  al_mesh <- pracma::cumtrapz(sqrt(apply(derivs^2, 1, sum)))
  theta_al <- pracma::interp1(x = theta_mesh, y = as.numeric(al_mesh), xi = thetas)
  return(cbind(al = theta_al, theta = thetas))
}


arc_length_surp_mixed <-
  function(mirt_mod,
           mirt_theta,
           mirt_items,
           tg_mod,
           tg_theta,
           tg_items,
           total_cat,
           max_surp = NULL,
           mesh_accuracy = 51) {
    al_vec <- numeric(length(mirt_theta) + 1)
    w_list <- tg_mod$parList[[length(tg_mod$parList)]]$WfdList

    start_mirt <- min(mirt_theta)
    start_tg <- 0

    for (th_ind in seq_along(mirt_theta)) {
      derivs <- matrix(0, mesh_accuracy, total_cat)

      mirt_mesh <- seq(start_mirt, mirt_theta[th_ind], length.out = mesh_accuracy)
      os_mesh <- seq(start_tg, tg_theta[th_ind], length.out = mesh_accuracy)
      start_mirt <- mirt_theta[th_ind]
      start_tg <- tg_theta[th_ind]

      m2 <- 0
      for (i in mirt_items) {
        item <- extract.item(mirt_mod, i)
        m1 <- m2 + 1
        m2 <- m2 + item@ncat
        deriv_resc <- (max(mirt_mesh) - min(mirt_mesh)) / (length(mirt_mesh) - 1)
        derivs[, m1:m2] <- deriv_mirt_gpc(item, mirt_mesh, surp = TRUE) * deriv_resc
        if (!is.null(max_surp)) {
          capped_surp <- -log2(probtrace(item, mirt_mesh)) > max_surp
          derivs[, m1:m2] <- (1 - capped_surp) * derivs[, m1:m2]
        }
      }

      for (i in tg_items) {
        ncat <- w_list[[i]]$M
        m1 <- m2 + 1
        m2 <- m2 + ncat
        deriv_resc <- (max(os_mesh) - min(os_mesh)) / (length(os_mesh) - 1)
        derivs[, m1:m2] <- deriv_resc * eval.surp(os_mesh, w_list[[i]]$Wfd, nderiv = 1) / log(2, ncat)
        if (!is.null(max_surp)) {
          capped_surp <- eval.surp(os_mesh, w_list[[i]]$Wfd) / log(2, ncat) > max_surp
          derivs[, m1:m2] <- (1 - capped_surp) * derivs[, m1:m2]
        }
      }

      als <- pracma::cumtrapz(sqrt(apply(derivs^2, 1, sum)))
      al_vec[th_ind + 1] <- al_vec[th_ind] + als[length(als)]
    }

    return(cbind(al = al_vec[-1], mirt_theta = mirt_theta, tg_theta = tg_theta))
  }
