source("1-sim-setup.R")

# Generate data -----------------------------------------------------------
print("Generating data...")
set.seed(1235)
for (sn in 1:length(scenarios)) {
  print(paste("scenario:", sn))
  scen <- scenarios[[sn]]
  cl <- makeCluster(detectCores(logical = F))
  registerDoParallel(cl)
  foreach(r = 1:iterations, .packages = c("mirt", "TestGardener")) %dopar% {
    os_sample <- Usimulate(scen$test_taker_thetas[, "tg"],
                           os_mod$parList[[length(os_mod$parList)]]$WfdList)
    data <- os_sample[, scen$os_items]-1
    if (length(scen$mirt_items) != 0) {
      mirt_sample <- simdata(
        model = mirt_mod,
        Theta = scen$test_taker_thetas[, "mirt", drop = F],
        which.items = scen$mirt_items
      )
      data <- cbind(mirt_sample, data)
    }
    
    if (length(scen$mirt_items) == 0) 
    {
      colnames(data) <- paste("I", scen$os_items, sep = "")
    } else if (length(scen$os_items) == 0) 
    {
      colnames(data) <- paste("I", scen$mirt_items, sep = "")
    } else {
      colnames(data) <- c(paste("I", scen$mirt_items, sep = ""), 
                          paste("I", scen$os_items, sep = ""))
    }
    
    file_name <- paste("simulation-data/datasets/r", r,"-", scen$name, ".RData", sep = "")
    save(data, file = file_name)
    NULL
  }
  stopCluster(cl)
}
