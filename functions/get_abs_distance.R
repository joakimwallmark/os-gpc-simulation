prob_mat_to_entropy <- function(p, base = 2) {
  rowSums(-p * log(p, base = base))
}

surp_mat_to_entropy <- function(s, base = 2) {
  rowSums((1 / base^s) * s)
}

# Returns sum mapping matrix for thetas number of evenly spread theta values
# modes are surp, prob and entropy
get_dist_travelled_tg <-
  function(mod,
           items = seq_along(mod$parList[[iter]]$WfdList),
           theta = seq(0, 100, length.out = 10000),
           mode = "surp",
           remove0 = FALSE) {
    iter <- length(mod$parList)
    dist <- vector("numeric", length = length(theta))
    for (item in items) {
      surps <- eval.surp(theta, mod$parList[[iter]]$WfdList[[item]]$Wfd)
      categories <- ncol(surps)
      surps <- surps / log(2, categories) # always use 2bit
      total_item_dist <- 0
      if (mode == "entropy") {
        entr <- surp_mat_to_entropy(surps)
        for (i in 2:length(entr)) {
          # add up total item distance to theta corresponding to surp[i, ]
          total_item_dist <- total_item_dist + sum(abs(entr[i] - entr[i - 1]))
          # add this to the total for the same theta
          dist[i] <- dist[i] + total_item_dist
        }
      }
      if (mode == "surp") {
        if (mode == remove0) { # ignore 0 event. We want bits for more learning only
          surps <- surps[, -1, drop = FALSE]
        }
        for (i in 2:nrow(surps)) {
          total_item_dist <- total_item_dist + sum(abs(surps[i, ] - surps[i - 1, ]))
          dist[i] <- dist[i] + total_item_dist
        }
      }
      if (mode == "prob") {
        probs <- 1 / 2^surps
        if (mode == remove0) { # ignore 0 event. We want bits for more learning only
          surps <- surps[, -1, drop = FALSE]
        }
        for (i in 2:nrow(surps)) {
          total_item_dist <- total_item_dist + sum(abs(probs[i, ] - probs[i - 1, ]))
          dist[i] <- dist[i] + total_item_dist
        }
      }
    }
    return(cbind(theta = theta, distance = dist))
  }


get_dist_travelled_mirt <-
  function(mod,
           items = 1:extract.mirt(mod, "nitems"),
           theta = c(-8, seq(-4, 4, length.out = 99998), 8),
           mode = "surp",
           surp_limit = 5,
           remove0 = FALSE) {
    dist <- vector("numeric", length = length(theta))
    for (item in items) {
      extr <- extract.item(mod, item)
      probs <- probtrace(extr, theta)
      total_item_dist <- 0
      if (mode == "entropy") {
        entr <- prob_mat_to_entropy(probs)
        for (i in 2:length(entr)) {
          # add up total item distance to theta corresponding to surp[i, ]
          total_item_dist <- total_item_dist + sum(abs(entr[i] - entr[i - 1]))
          # add this to the total for the same theta
          dist[i] <- dist[i] + total_item_dist
        }
      }
      if (mode == "surp") {
        surps <- -log(probs, base = 2)
        surps[surps > surp_limit] <- surp_limit # limit surprisals
        if (mode == remove0) { # ignore 0 event. We want bits for more learning only
          surps <- surps[, !grepl(".P.0", colnames(surps))]
        }
        for (i in 2:nrow(surps)) {
          total_item_dist <- total_item_dist + sum(abs(surps[i, ] - surps[i - 1, ]))
          # add this to the total for the same theta
          dist[i] <- dist[i] + total_item_dist
        }
      }
      if (mode == "prob") {
        if (mode == remove0) { # ignore 0 event. We want bits for more learning only
          probs <- probs[, !grepl(".P.0", colnames(probs))]
        }
        for (i in 2:nrow(probs)) {
          total_item_dist <- total_item_dist + sum(abs(probs[i, ] - probs[i - 1, ]))
          # add this to the total for the same theta
          dist[i] <- dist[i] + total_item_dist
        }
      }
    }
    return(cbind(theta = theta, distance = dist))
  }
