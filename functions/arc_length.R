source("functions/tg_density.R")
library(pracma)
library(TestGardener)
library(mirt)

# stolen from mirt and rewritten to remove hessian and support surprisal deriv.
deriv_mirt_gpc <- function(x, Theta, surp = T) 
{
  Theta <- as.matrix(Theta)
  a <- mirt:::ExtractLambdas(x)
  d <- mirt:::ExtractZetas(x)
  ak <- 0:(x@ncat - 1L)
  P <- mirt:::P.nominal(c(a, ak, d), ncat = length(d), Theta = Theta)
  Num <- mirt:::P.nominal(c(a, ak, d), ncat = length(d), Theta = Theta, 
                          returnNum = TRUE)
  Den <- rowSums(Num)
  # grad <- hess <- vector("list", x@ncat)
  # for (i in seq_len(x@ncat)) grad[[i]] <- hess[[i]] <- matrix(0, nrow(Theta), x@nfact)
  grad <- vector("list", x@ncat)
  for (i in seq_len(x@ncat)) grad[[i]] <- matrix(0, nrow(Theta), x@nfact)
  for (j in seq_len(x@nfact)) { # j are factors, i are item categories
    for (i in seq_len(x@ncat)) {
      grad[[i]][, j] <- ak[i] * a[j] * P[, i] - P[, i] * 
        (Num %*% (ak * a[j]))/Den
      # hess[[i]][, j] <- ak[i]^2 * a[j]^2 * P[, i] - 2 * 
      #   ak[i] * a[j] * P[, i] * (Num %*% (ak * a[j]))/Den + 
      #   2 * P[, i] * ((Num %*% (ak * a[j]))/Den)^2 - 
      #   P[, i] * ((Num %*% (ak^2 * a[j]^2))/Den)
      if (surp) { 
        grad[[i]][, j] <- -grad[[i]][, j]/(P[, i]*log(2))
      }
    }
  }
  # return(list(grad = grad, hess = hess))
  return(do.call(cbind, grad))
}

arc_length_surp <- function(mod, thetas, total_cat, max_surp = NULL, item_ids = NULL, mesh_accuracy = 1001) {
  no_thetas <- length(thetas)
  al <- vector("numeric", no_thetas)
  derivs <- matrix(0, mesh_accuracy, total_cat)
  m2 <- 0
  if (class(mod)=="SingleGroupClass") { # check if mirt model
    if (is.null(item_ids)) item_ids <- 1:extract.mirt(mod, "nitems")
    # We want item category derivatives in columns and theta values in rows
    theta_mesh <- seq(min(thetas), max(thetas), length.out = mesh_accuracy)
    for (i in item_ids) {
      item <- extract.item(mod, i)
      m1 <- m2 + 1
      m2 <- m2 + item@ncat
      deriv_resc <- (max(theta_mesh)-min(theta_mesh))/(length(theta_mesh)-1)
      derivs[, m1:m2] <- deriv_mirt_gpc(item, theta_mesh, surp = T)*deriv_resc
      if (!is.null(max_surp)) {
        capped_surp <- -log2(probtrace(item, theta_mesh))>max_surp
        derivs[, m1:m2] <- (1-capped_surp)*derivs[, m1:m2]
      }
    }
  }
  else { # assume TestGardener
    w_list <- mod$parList[[length(mod$parList)]]$WfdList
    if (is.null(item_ids)) item_ids <- 1:length(w_list)
    theta_mesh <- seq(0, 100, length.out = mesh_accuracy)
    for (i in item_ids) {
      ncat <- w_list[[i]]$M
      m1 = m2 + 1
      m2 = m2 + ncat
      deriv_resc <- (max(theta_mesh)-min(theta_mesh))/(length(theta_mesh)-1)
      derivs[, m1:m2] <- deriv_resc*eval.surp(theta_mesh, w_list[[i]]$Wfd, nderiv = 1)/log(2, ncat)
      if (!is.null(max_surp)) {
        capped_surp <- eval.surp(theta_mesh, w_list[[i]]$Wfd)/log(2, ncat)>max_surp
        derivs[, m1:m2] <- (1-capped_surp)*derivs[, m1:m2]
      }
    }
  }
  al_mesh <- pracma::cumtrapz(sqrt(apply(derivs^2, 1, sum)))
  theta_al <- pracma::interp1(x = theta_mesh, y = as.numeric(al_mesh), xi = thetas)
  return(cbind(al = theta_al, theta = thetas))
}


arc_length_surp_mixed <- function(mirt_mod, mirt_theta, mirt_items, tg_mod, tg_theta, tg_items, total_cat, max_surp = NULL, mesh_accuracy = 51) {
  al_vec <- numeric(length(mirt_theta)+1)
  w_list <- tg_mod$parList[[length(tg_mod$parList)]]$WfdList
  
  start_mirt <- min(mirt_theta)
  start_tg <- 0
  
  for (th_ind in 1:length(mirt_theta)) {
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
      deriv_resc <- (max(mirt_mesh)-min(mirt_mesh))/(length(mirt_mesh)-1)
      derivs[, m1:m2] <- deriv_mirt_gpc(item, mirt_mesh, surp = T)*deriv_resc
      if (!is.null(max_surp)) {
        capped_surp <- -log2(probtrace(item, mirt_mesh))>max_surp
        derivs[, m1:m2] <- (1-capped_surp)*derivs[, m1:m2]
      }
    }
    
    for (i in tg_items) {
      ncat <- w_list[[i]]$M
      m1 = m2 + 1
      m2 = m2 + ncat
      deriv_resc <- (max(os_mesh)-min(os_mesh))/(length(os_mesh)-1)
      derivs[, m1:m2] <- deriv_resc*eval.surp(os_mesh, w_list[[i]]$Wfd, nderiv = 1)/log(2, ncat)
      if (!is.null(max_surp)) {
        capped_surp <- eval.surp(os_mesh, w_list[[i]]$Wfd)/log(2, ncat)>max_surp
        derivs[, m1:m2] <- (1-capped_surp)*derivs[, m1:m2]
      }
    }
    
    als <- pracma::cumtrapz(sqrt(apply(derivs^2, 1, sum)))
    al_vec[th_ind+1] <- al_vec[th_ind]+als[length(als)]
    # theta_al <- pracma::interp1(x = theta_mesh, y = as.numeric(al_mesh), xi = thetas)
  }
  
  return(cbind(al = al_vec[-1], mirt_theta = mirt_theta, tg_theta = tg_theta))
}

# prob_mat_to_entropy <- function(p, base = 2) {
#   rowSums(-p*log(p, base = base))
# }
# 
# surp_mat_to_entropy <- function(s, base = 2) {
#   rowSums((1/base^s)*s)
# }
# 
# # assumes thetas are sorted
# arc_length_all_mixed_mod <- function(mirt_thetas, tg_dens, mirt_mod, mirt_items, tg_mod, tg_items, mode = "surp", surp_limit = 5) {
#   no_thetas <- length(mirt_thetas)
#   al <- vector("numeric", no_thetas)
#   tol = 5e-05
#   if (mode == "surp") {
#     manifold_func <- mixed_surps
#     tol = 5e-04 # higher tol. This is slow otherwise
#   }
#   else if (mode == "prob") {
#     manifold_func <- mixed_probs
#   }
#   else if (mode == "entropy") {
#     manifold_func <- mixed_entropy
#   }
#   al[1] <- 0 # first test taker is where we start, thus 0 arc length
#   for (i in 2:no_thetas) {
#     al[i] <- al[i-1] + arclength(manifold_func, mirt_thetas[i-1], mirt_thetas[i],
#                                  tg_dens = tg_dens,
#                                  mirt_mod = mirt_mod,
#                                  mirt_items = mirt_items,
#                                  tg_mod = tg_mod,
#                                  tg_items = tg_items,
#                                  surp_limit = surp_limit,
#                                  tol = tol)$length
#   }
#   return(al)
# }
# 
# mixed_surps <- function(mirt_theta, tg_dens, mirt_mod, mirt_items, tg_mod, tg_items, surp_limit = 5) {
#   # get mirt surps
#   surps <- vector("numeric", 0)
#   for (item in mirt_items) {
#     extr <- extract.item(mirt_mod, item)
#     probs <- probtrace(extr, mirt_theta)
#     surps <-  c(surps, -log(probs, 2)) # always use 2bit
#   }
#   surps[surps > surp_limit] <- surp_limit # limit surprisals for unlikely events
#   
#   # get corresponding tg_theta
#   tg_theta <- tg_qss(pnorm(mirt_theta), tg_dens)
#   # get testgardener surp
#   iter <- length(tg_mod$parList)
#   for (item in tg_items) {
#     WListi <- tg_mod$parList[[iter]]$WfdList[[item]]
#     Wfdi <- WListi$Wfd
#     surps <- c(surps, eval.surp(tg_theta, Wfdi)/log(2, WListi$M) )# always use 2bit
#   }
#   
#   return(surps)
# }

# arc_length_all <- function(mod, start, thetas, mode = "surp", tol = 5e-04) {
#   no_thetas <- length(thetas)
#   al <- vector("numeric", no_thetas)
#   if (class(mod)=="SingleGroupClass") { # check if mirt model
#     if (mode == "surp") {
#       manifold_func <- mirt_surps
#     }
#     else if (mode == "prob") {
#       manifold_func <- mirt_probs
#     }
#     else if (mode == "entropy") {
#       manifold_func <- mirt_entropy
#     }
#   }
#   else { # assume testGardener
#     if (mode == "surp") {
#       manifold_func <- tg_surps
#     }
#     else if (mode == "prob") {
#       manifold_func <- tg_probs
#     }
#     else if (mode == "entropy") {
#       manifold_func <- tg_entropy
#     }
#   }
#   
#   al[1] <- arclength(manifold_func, start, thetas[1], mod = mod, tol = tol)$length
#   for (i in 2:no_thetas) {
#     al[i] <- al[i-1] + arclength(manifold_func, thetas[i-1], thetas[i], mod = mod, tol = tol)$length
#   }
#   return(cbind(theta = thetas, al = al))
# }
# 
# 
# arc_length <- function(mod, from, to, mode = "surp") {
#   tol = 5e-05
#   if (class(mod)=="SingleGroupClass") { # check if mirt model
#     if (mode == "surp") {
#       manifold_func <- mirt_surps
#       tol = 5e-04 # higher tol. This is slow otherwise
#     }
#     else if (mode == "prob") {
#       manifold_func <- mirt_probs
#     }
#     else if (mode == "entropy") {
#       manifold_func <- mirt_entropy
#     }
#   } 
#   else {
#     if (mode == "surp") {
#       manifold_func <- tg_surps
#       tol = 5e-03 # higher tol. This is slow otherwise
#     }
#     else if (mode == "prob") {
#       manifold_func <- tg_probs
#     }
#     else if (mode == "entropy") {
#       manifold_func <- tg_entropy
#     }
#   }
#   arclength(manifold_func, from, to, mod = mod, tol = tol)
# }
# 
# tg_probs <- function(theta, mod) {
#   iter <- length(mod$parList)
#   no_items <- length(mod$parList[[iter]]$WfdList)
#   probs <- vector("numeric", 0)
#   for (item in 1:no_items) {
#     WListi <- mod$parList[[iter]]$WfdList[[item]]
#     Wfdi <- WListi$Wfd
#     # logMi <- log(WListi$M)
#     # surpisal <- eval.surp(theta, Wfdi)
#     # probs <- c(probs, exp(-surpisal * logMi))
#     surpisal <- eval.surp(theta, Wfdi)
#     probs <- c(probs, WListi$M^(-surpisal))
#   }
#   probs
# }
# 
# tg_surps <- function(theta, mod) {
#   iter <- length(mod$parList)
#   no_items <- length(mod$parList[[iter]]$WfdList)
#   surps <- vector("numeric", 0)
#   for (item in 1:no_items) {
#     WListi <- mod$parList[[iter]]$WfdList[[item]]
#     Wfdi <- WListi$Wfd
#     # surps <- c(surps, eval.surp(theta, Wfdi))
#     surps <- c(surps, eval.surp(theta, Wfdi)/log(2, WListi$M) )# always use 2bit
#   }
#   surps
# }
# 
# tg_entropy <- function(theta, mod) {
#   iter <- length(mod$parList)
#   no_items <- length(mod$parList[[iter]]$WfdList)
#   entropies <- vector("numeric", 0)
#   for (item in 1:no_items) {
#     WListi <- mod$parList[[iter]]$WfdList[[item]]
#     Wfdi <- WListi$Wfd
#     surps <- eval.surp(theta, Wfdi)/log(2, WListi$M) # always use 2bit
#     entropies <- c(entropies, surp_mat_to_entropy(surps))
#   }
#   entropies
# }
# 
# mirt_probs <- function(theta, mod) {
#   c(probtrace(mod, theta))
# }
# 
# mirt_surps <- function(theta, mod, surp_limit = 5) {
#   probs <- c(probtrace(mod, theta))
#   surps <- -log(probs, 2) # always use 2bit
#   surps[surps > surp_limit] <- surp_limit # limit surprisals for unlikely events
#   surps
# }
# 
# mirt_entropy <- function(theta, mod) {
#   entropies <- vector("numeric", 0)
#   for (item in 1:extract.mirt(mirt_mod, "nitems")) {
#     extr <- extract.item(mirt_mod, item)
#     probs <- probtrace(extr, theta)
#     entropies <- c(entropies, prob_mat_to_entropy(probs))
#   }
#   entropies
# }
# 

# 
# mixed_probs <- function(mirt_theta, tg_dens, mirt_mod, mirt_items, tg_mod, tg_items) {
#   # get mirt probs
#   probs <- vector("numeric", 0)
#   for (item in mirt_items) {
#     extr <- extract.item(mirt_mod, item)
#     probs <- c(probs, probtrace(extr, mirt_theta))
#   }
#   
#   # get corresponding tg_theta
#   tg_theta <- tg_qss(pnorm(mirt_theta), tg_dens)
#   # get testgardener prob
#   iter <- length(tg_mod$parList)
#   for (item in tg_items) {
#     WListi <- tg_mod$parList[[iter]]$WfdList[[item]]
#     Wfdi <- WListi$Wfd
#     surps <- eval.surp(tg_theta, Wfdi) 
#     probs <- c(probs, WListi$M^(-surps))
#   }
#   
#   return(probs)
# }
# 
# mixed_entropy <- function(mirt_theta, tg_dens, mirt_mod, mirt_items, tg_mod, tg_items) {
#   entropies <- vector("numeric", 0)
#   for (item in mirt_items) {
#     extr <- extract.item(mirt_mod, item)
#     probs <- probtrace(extr, mirt_theta)
#     entropies <- c(entropies, prob_mat_to_entropy(probs))
#   }
#   
#   # get corresponding tg theta
#   tg_theta <- tg_qss(pnorm(mirt_theta), tg_dens)
#   iter <- length(tg_mod$parList)
#   for (item in tg_items) {
#     WListi <- tg_mod$parList[[iter]]$WfdList[[item]]
#     Wfdi <- WListi$Wfd
#     surps <- eval.surp(tg_theta, Wfdi)/log(2, WListi$M)
#     entropies <- c(entropies, surp_mat_to_entropy(surps))
#   }
#   
#   return(entropies)
# }

