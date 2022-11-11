prob_mat_to_entropy <- function(p, base = 2) {
  rowSums(-p*log(p, base = base))
}

surp_mat_to_entropy <- function(s, base = 2) {
  rowSums((1/base^s)*s)
}

# Returns sum mapping matrix for thetas number of evenly spread theta values
# modes are surp, prob and entropy
get_dist_travelled_tg <- function(mod, items = 1:length(mod$parList[[iter]]$WfdList), theta = seq(0, 100, length.out = 10000), 
                                  mode = "surp", remove0 = F) {
  iter <- length(mod$parList)
  dist <- vector("numeric", length = length(theta))
  for (item in items) {
    surps <- eval.surp(theta, mod$parList[[iter]]$WfdList[[item]]$Wfd)
    categories <- ncol(surps)
    surps <- surps/log(2, categories) # always use 2bit
    total_item_dist <- 0
    if (mode == "entropy") {
      entr <- surp_mat_to_entropy(surps)
      for (i in 2:length(entr)) {
        # add up total item distance to theta corresponding to surp[i, ]
        total_item_dist <- total_item_dist + sum(abs(entr[i]-entr[i-1]))
        # add this to the total for the same theta
        dist[i] <- dist[i] + total_item_dist
      }
    }
    if (mode == "surp") {
      if (mode == remove0) { # ignore 0 event. We want bits for more learning only
        surps <- surps[, -1, drop = F]
      }
      for (i in 2:nrow(surps)) {
        total_item_dist <- total_item_dist + sum(abs(surps[i, ]-surps[i-1, ]))
        dist[i] <- dist[i] + total_item_dist
      }
    }
    if (mode == "prob") {
      probs <- 1/2^surps
      if (mode == remove0) { # ignore 0 event. We want bits for more learning only
        surps <- surps[, -1, drop = F]
      }
      for (i in 2:nrow(surps)) {
        total_item_dist <- total_item_dist + sum(abs(probs[i, ]-probs[i-1, ]))
        dist[i] <- dist[i] + total_item_dist
      }
    }
  }
  return(cbind(theta = theta, distance = dist))
}


get_dist_travelled_mirt <- function(mod, items=1:extract.mirt(mod, "nitems"), theta = c(-8, seq(-4, 4, length.out = 99998), 8), 
                                    mode = "surp", surp_limit = 5, remove0 = F) {
  dist <- vector("numeric", length = length(theta))
  for (item in items) {
    extr <- extract.item(mod, item)
    probs <- probtrace(extr, theta)
    total_item_dist <- 0
    if (mode == "entropy") {
      entr <- prob_mat_to_entropy(probs)
      for (i in 2:length(entr)) {
        # add up total item distance to theta corresponding to surp[i, ]
        total_item_dist <- total_item_dist + sum(abs(entr[i]-entr[i-1]))
        # add this to the total for the same theta
        dist[i] <- dist[i] + total_item_dist
      }
    }
    if (mode == "surp") {
      surps <-  -log(probs, base = 2)
      surps[surps > surp_limit] <- surp_limit # limit surprisals 
      if (mode == remove0) { # ignore 0 event. We want bits for more learning only
        surps <- surps[, !grepl(".P.0", colnames(surps))]
      }
      for (i in 2:nrow(surps)) {
        total_item_dist <- total_item_dist + sum(abs(surps[i, ]-surps[i-1, ]))
        # add this to the total for the same theta
        dist[i] <- dist[i] + total_item_dist
      }
    }
    if (mode == "prob") {
      if (mode == remove0) { # ignore 0 event. We want bits for more learning only
        probs <- probs[, !grepl(".P.0", colnames(probs))]
      }
      for (i in 2:nrow(probs)) {
        total_item_dist <- total_item_dist + sum(abs(probs[i, ]-probs[i-1, ]))
        # add this to the total for the same theta
        dist[i] <- dist[i] + total_item_dist
      }
    }
  }
  return(cbind(theta = theta, distance = dist))
}




# # Returns sum mapping matrix for thetas number of evenly spread theta values
# # modes are surpsum, entopy and remove0
# get_dist_travelled_tg2 <- function(mod, from, to, thetas = 10000, 
#                                    mode = "surpsum") {
#   theta <- seq(from, to, length.out = thetas)
#   iter <- length(mod$parList)
#   no_items <- length(mod$parList[[iter]]$WfdList)
#   dist <- vector("numeric", length = thetas)
#   for (item in 1:no_items) {
#     wfd <- mod$parList[[iter]]$WfdList[[item]]$Wfd
#     surps <- eval.surp(theta, wfd)
#     total_item_dist <- 0
#     if (mode == "entropy") {
#       entr <- surp_mat_to_entropy(surps)
#       for (i in 2:length(entr)) {
#         # add up total item distance to theta corresponding to surp[i, ]
#         total_item_dist <- total_item_dist + sum(abs(entr[i]-entr[i-1]))
#         # add this to the total for the same theta
#         dist[i] <- dist[i] + total_item_dist
#       }
#     }
#     if (mode == "surpsum" | mode == "remove0") {
#       if (mode == "remove0") { # ignore 0 event. We want bits for more learning only
#         surps <- surps[, -1, drop = F]
#       }
#       for (i in 2:nrow(surps)) {
#         total_item_dist <- total_item_dist + sum(abs(surps[i, ]-surps[i-1, ]))
#         dist[i] <- dist[i] + total_item_dist
#       }
#     }
#   }
#   cbind(theta = theta, distance = dist)
# }
# 
# get_dist_travelled_mirt2 <- function(mod, from, to, thetas = 1000, 
#                                         mode = "surpsum") {
#   theta <- matrix(seq(from, to, length.out = thetas))
#   dist <- vector("numeric", length = thetas)
#   if (mode == "entropy") {
#     no_items <- length(coef(mod))-1
#     for (item in 1:no_items) {
#       total_item_dist <- 0
#       extr <- extract.item(mod, item)
#       entr <- prob_mat_to_entropy(probtrace(extr, theta))
#       for (i in 2:length(entr)) {
#         # add up total item distance to theta corresponding to surp[i, ]
#         total_item_dist <- total_item_dist + sum(abs(entr[i]-entr[i-1]))
#         # add this to the total for the same theta
#         dist[i] <- dist[i] + total_item_dist
#       }
#     }
#   }
#   
#   if (mode == "surpsum" | mode == "remove0") {
#     surps <-  -log(probtrace(mod, theta), base = 2)
#     if (mode == "remove0") { # ignore 0 event. We want bits for more learning only
#       surps <- surps[, !grepl(".P.0", colnames(surps))]
#     }
#     for (i in 2:nrow(surps)) {
#       dist <- dist + sum(abs(surps[i, ]-surps[i-1, ]))
#       # add this to the total for the same theta
#       dist[i] <- dist[i-1] + sum(abs(surps[i, ]-surps[i-1, ]))
#     }
#   }
#   dist
# }
# 

# 
# get_dist_travelled_test_takers <- function(mod, thetas, from, acc = 1000, 
#                                mode = "surpsum", remove_0 = F) {
#   sorted <- c(from, sort(thetas)) # add 0 for loop, remove later
#   dists <- vector("numeric", length = length(sorted))
#   max <- max(sorted)
#   for (i in 2:length(dists)) {
#     if (sorted[i] == sorted[i-1]) {
#       dists[i] <- dists[i-1]
#     }
#     else {
#       new_acc <- max(2, acc*(sorted[i]-sorted[i-1])/max) # set new accuracy
#       # new_acc <- max(sorted[i]-sorted[i-1], max*acc/(sorted[i]-sorted[i-1])) # set new accuracy
#       if (class(mod)=="SingleGroupClass") { # check if mirt model
#         dists[i] <- dists[i-1]+get_dist_travelled_mirt_old(mod, sorted[i-1], sorted[i], acc = new_acc, mode = mode)
#       }
#       else { # assume testGardener
#         dists[i] <- dists[i-1]+get_dist_travelled_tg_old(mod, sorted[i-1], sorted[i], acc = new_acc, mode = mode)
#       }
#     }
#   }
#   dists <- dists[-1]
#   dists[match(thetas, sorted[-1])] # sort distances in original test taker order
# }
# 
# # modes are surpsum, entopy and remove0
# get_dist_travelled_mirt_old <- function(mod, from, to, acc = 1000, 
#                                     mode = "surpsum") {
#   theta <- matrix(seq(from, to, length.out = acc))
#   dist <- 0
#   if (mode == "entropy") {
#     no_items <- length(coef(mod))-1
#     for (item in 1:no_items) {
#       extr <- extract.item(mod, item)
#       entr <- prob_mat_to_entropy(probtrace(extr, theta))
#       for (i in 2:length(entr)) {
#         dist <- dist + sum(abs(entr[i]-entr[i-1]))
#       }
#     }
#   }
#   
#   if (mode == "surpsum" | mode == "remove0") {
#     surps <-  -log(probtrace(mod, theta), base = 2)
#     if (mode == "remove0") { # ignore 0 event. We want bits for more learning only
#       surps <- surps[, !grepl(".P.0", colnames(surps))]
#     }
#     for (i in 2:nrow(surps)) {
#       dist <- dist + sum(abs(surps[i, ]-surps[i-1, ]))
#     }
#   }
#   dist
# }
# 
# # modes are surpsum, entopy and remove0
# get_dist_travelled_tg_old <- function(mod, from, to, acc = 500, 
#                                     mode = "surpsum") {
#   theta <- matrix(seq(from, to, length.out = acc))
#   iter <- length(mod$parList)
#   no_items <- length(mod$parList[[iter]]$WfdList)
#   dist <- 0
#   for (item in 1:no_items) {
#     wfd <- mod$parList[[iter]]$WfdList[[item]]$Wfd
#     surps <- eval.surp(theta, wfd)
#     if (mode == "entropy") {
#       entr <- surp_mat_to_entropy(surps)
#       for (i in 2:length(entr)) {
#         dist <- dist + sum(abs(entr[i]-entr[i-1]))
#       }
#     }
#     if (mode == "surpsum" | mode == "remove0") {
#       if (mode == "remove0") { # ignore 0 event. We want bits for more learning only
#         surps <- surps[, -1, drop = F]
#       }
#       for (i in 2:nrow(surps)) {
#         dist <- dist + sum(abs(surps[i, ]-surps[i-1, ]))
#       }
#     }
#   }
#   dist
# }
