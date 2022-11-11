sim_sum_table <- function(table, cols = c(1:4, 9:14, 5:6)) {
  table <- table[, cols]
  table[, 1] <- table[, 1] == "skewed"
  table[, 5:10] <- table[, 5:10]*100
  table[, -(1:10)] <- round(table[, -(1:10)], 3)
  table[, 5:10] <- round(table[, 5:10], 2)
  
  table <- table[
    with(table, order(table[, 1], table[, 2], table[, 3], table[, 4])),
  ]
  # Bold the smallest value in each column
  for (c in seq(5, ncol(table), by = 2)) {
    table[, c] = ifelse(table[, c] <= table[, c+1], 
                        paste0('\\textbf{', table[, c],'}'), 
                        table[, c])
    table[, c+1] = ifelse(table[, c+1] <= table[, c], 
                          paste0('\\textbf{', table[, c+1],'}'), 
                          table[, c+1])
  }
  print(xtable(cbind("", table[, -1]), digits = c(rep(0, 4), 1, rep(2, 6), 3, 3)), 
              type="latex", comment=F, include.rownames=FALSE, 
              sanitize.text.function = function(x) x, table.placement = "!ht")
}

sim_item_table <- function(measure_list) {
  item_table <- data.frame("", 
                           1:ncol(measure_list$avg_mirt_kl),
                           sapply(measure_list$avg_mirt_bias, function(x) { mean(abs(x)) } ),
                           sapply(measure_list$avg_os_bias, function(x) { mean(abs(x)) } ),
                           sapply(measure_list$avg_mirt_see, function(x) { mean(x) } ),
                           sapply(measure_list$avg_os_see, function(x) { mean(x) } ),
                           sapply(measure_list$avg_mirt_rmse, function(x) { mean(x) } ),
                           sapply(measure_list$avg_os_rmse, function(x) { mean(x) } ),
                           colMeans(measure_list$avg_mirt_kl),
                           colMeans(measure_list$avg_os_kl))
  colnames(item_table) <- c("True item type", "item", 
                            "mirt bias", "os bias", 
                            "mirt SE", "os SE", 
                            "mirt RMSE", "os RMSE", 
                            "mirt KL div.", "os KL div.")  
  
  item_table[, -(1:2)] <- round(item_table[, -(1:2)], 3)*100
  
  # Bold the smallest value in each column
  for (c in seq(3, ncol(item_table), by = 2)) {
    item_table[, c] = ifelse(item_table[, c] <= item_table[, c+1], 
                             paste0('\\textbf{', item_table[, c],'}'), 
                             item_table[, c])
    item_table[, c+1] = ifelse(item_table[, c+1] <= item_table[, c], 
                               paste0('\\textbf{', item_table[, c+1],'}'), 
                               item_table[, c+1])
  }
  
  print(xtable(item_table), 
        type="latex", comment=F, include.rownames=FALSE, 
        sanitize.text.function = function(x) x, table.placement = "!ht") 
}