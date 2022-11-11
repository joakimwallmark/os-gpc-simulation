get_simulation_measures <- function(item_fit, true_probs, items, item_ind, response_cats) {
  # Matrices with measures. rows thetas and cols are items 
  avg_mirt_kl <- avg_mirt_kl_see <- avg_os_kl <- avg_os_kl_see <- 
    matrix(0, nrow = nrow(scen$test_taker_thetas), ncol = length(items))
  avg_mirt_bias <- avg_mirt_see <- avg_mirt_rmse <- avg_os_bias <- avg_os_see <- 
    avg_os_rmse <- list()
  tab_item_i <- 1
  for (item in items) {
    avg_mirt_kl[, tab_item_i] <- colMeans(item_fit$kl_mirt[[item]])
    avg_mirt_kl_see[, tab_item_i] <- apply(item_fit$kl_mirt[[item]], 2, sd)
    avg_os_kl[, tab_item_i] <- colMeans(item_fit$kl_tg[[item]])
    avg_os_kl_see[, tab_item_i] <- apply(item_fit$kl_tg[[item]], 2, sd)
    
    avg_mirt_bias[[tab_item_i]] <- avg_mirt_see[[tab_item_i]] <- avg_mirt_rmse[[tab_item_i]] <- 
      avg_os_bias[[tab_item_i]] <- avg_os_see[[tab_item_i]] <- avg_os_rmse[[tab_item_i]] <- 
      matrix(0, nrow = nrow(scen$test_taker_thetas), ncol = response_cats[item_ind[item]])
    for (cat in 1:response_cats[item_ind[item]]) {
      avg_mirt_bias[[tab_item_i]][, cat] <- 
        rowMeans(t(item_fit$probs_mirt[[item]][, , cat])-true_probs[[item]][, cat])
      avg_os_bias[[tab_item_i]][, cat] <- 
        rowMeans(t(item_fit$probs_tg[[item]][, , cat])-true_probs[[item]][, cat])
      
      avg_mirt_see[[tab_item_i]][, cat] <- 
        apply(item_fit$probs_mirt[[item]][, , cat], 2, sd)
      avg_os_see[[tab_item_i]][, cat] <- 
        apply(item_fit$probs_tg[[item]][, , cat], 2, sd)
      
      avg_mirt_rmse[[tab_item_i]][, cat] <- 
        sqrt(rowMeans((t(item_fit$probs_mirt[[item]][, , cat])-true_probs[[item]][, cat])^2))
      avg_os_rmse[[tab_item_i]][, cat] <- 
        sqrt(rowMeans((t(item_fit$probs_tg[[item]][, , cat])-true_probs[[item]][, cat])^2))
    }
    tab_item_i <- tab_item_i + 1
  }
  list(
    avg_mirt_kl = avg_mirt_kl,
    avg_os_kl = avg_os_kl,
    avg_mirt_kl_see = avg_mirt_kl_see,
    avg_os_kl_see = avg_os_kl_see,
    avg_mirt_bias = avg_mirt_bias,
    avg_os_bias = avg_os_bias,
    avg_mirt_see = avg_mirt_see,
    avg_os_see = avg_os_see,
    avg_mirt_rmse = avg_mirt_rmse,
    avg_os_rmse = avg_os_rmse
  )
}