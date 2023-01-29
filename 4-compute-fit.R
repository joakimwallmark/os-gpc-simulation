source("1-sim-setup.R")
source("functions/sim_measures.R")

print("Calculating fit MSE/KL...")
cl <- makeCluster(min(detectCores(logical = FALSE), 4)) # Too many cores needs too much RAM memory
registerDoParallel(cl)
res <- foreach(sn = seq_along(scenarios), .packages = c("mirt", "TestGardener")) %dopar% {
  scen <- scenarios[[sn]]
  no_thetas <- nrow(scen$test_taker_thetas)
  no_items <- length(scen$mirt_items) + length(scen$os_items)
  item_ind <- c(scen$mirt_items, scen$os_items)

  # Get true thetas for scenario
  true_probs <- list()
  # true mirt items
  if (length(scen$mirt_items) > 0) {
    for (item in seq_along(scen$mirt_items)) {
      t_extr <- extract.item(mirt_mod, scen$mirt_items[item])
      true_probs[[item]] <- probtrace(t_extr, scen$test_taker_thetas[, "mirt"])
    }
  }
  # true tg items
  if (length(scen$os_items) > 0) {
    for (item in (length(scen$mirt_items) + 1):no_items) {
      w_list_i <- os_mod$parList[[os_cycles]]$WfdList[[scen$os_items[item - length(scen$mirt_items)]]]
      surpisal <- eval.surp(scen$test_taker_thetas[, "tg"], w_list_i$Wfd)
      true_probs[[item]] <- w_list_i$M^(-surpisal)
    }
  }

  item_fit <- list(
    probs_mirt = list(), probs_tg = list(),
    kl_mirt = list(), kl_tg = list()
  )
  for (item in 1:no_items) {
    item_fit$kl_mirt[[item]] <- matrix(0, nrow = iterations, ncol = no_thetas)
    item_fit$kl_tg[[item]] <- matrix(0, nrow = iterations, ncol = no_thetas)
    # 1 array index is iteration, second thetas and third is item category
    item_fit$probs_mirt[[item]] <- array(dim = c(iterations, no_thetas, response_cats[item_ind[item]]))
    item_fit$probs_tg[[item]] <- array(dim = c(iterations, no_thetas, response_cats[item_ind[item]]))
  }

  for (r in 1:iterations) {
    print(r)
    file_name <- paste("simulation-data/models/r", r, "-", scen$name, ".RData", sep = "")
    load(file = file_name)
    sim_mirt_thetas <- fscores(sim_mirt_mod)
    sim_os_thetas <- sim_os_mod$parList[[os_cycles]]$theta

    for (item in 1:no_items) {
      # mirt model probs
      extr <- extract.item(sim_mirt_mod, item)
      sim_mirt_probs <- probtrace(extr, sim_mirt_thetas)
      true_p <- true_probs[[item]]
      if (ncol(sim_mirt_probs) != ncol(true_p)) { # add missing category if it was all 0 in data
        temp <- matrix(0.00001, nrow = nrow(sim_mirt_probs), ncol = ncol(true_p)) # we can have 0 surp. set small prob
        ind <- 1
        for (categ in 0:(ncol(true_p) - 1)) {
          if (categ %in% unique(sim_mirt_mod@Data$data[, item])) {
            temp[, categ + 1] <- sim_mirt_probs[, ind]
            ind <- ind + 1
          }
        }
        sim_mirt_probs <- temp
      }
      kl_mirt <- rowSums(true_p * log(true_p, base = 2) - true_p * log(sim_mirt_probs, base = 2))
      # os model probs
      w_list_i <- sim_os_mod$parList[[os_cycles]]$WfdList[[item]]
      surpisal <- eval.surp(sim_os_thetas, w_list_i$Wfd)
      sim_tg_probs <- w_list_i$M^(-surpisal)
      true_p <- true_probs[[item]]
      kl_tg <- rowSums(true_p * log(true_p, base = 2) - true_p * log(sim_tg_probs, base = 2))

      item_fit$probs_mirt[[item]][r, , ] <- sim_mirt_probs
      item_fit$probs_tg[[item]][r, , ] <- sim_tg_probs
      item_fit$kl_mirt[[item]][r, ] <- kl_mirt
      item_fit$kl_tg[[item]][r, ] <- kl_tg
    }
  }
  all_items_measures <- get_simulation_measures(
    item_fit,
    true_probs,
    seq_along(true_probs),
    item_ind,
    response_cats
  )
  os_items_measures <- NULL
  mirt_items_measures <- NULL
  if (length(scen$mirt_items) > 0) {
    mirt_items_measures <- get_simulation_measures(
      item_fit,
      true_probs, seq_along(scen$mirt_items),
      item_ind,
      response_cats
    )
  }
  if (length(scen$os_items) > 0) {
    item_start_ind <- length(scen$mirt_items) + 1
    os_items_measures <- get_simulation_measures(
      item_fit,
      true_probs,
      item_start_ind:(item_start_ind + length(scen$os_items) - 1),
      item_ind,
      response_cats
    )
  }

  file_name <- paste("simulation-data/measures/", scen$name, ".RData", sep = "")
  save(all_items_measures, mirt_items_measures, os_items_measures, file = file_name)
  NULL
}
stopCluster(cl)
