source("1-sim-setup.R")

# Fit models --------------------------------------------------------------
item_score_start <- 0
print("Fitting models...")
for (sn in 1:length(scenarios)) {
  print(paste("scenario:", sn))
  scen <- scenarios[[sn]]
  opt_scr <-
    lapply(response_cats[c(scen$mirt_items, scen$os_items)], function(x) {
      item_score_start:(x-(1-item_score_start))
    })
  opt_list <- list(itemLab=NULL, optLab=NULL, optScr=opt_scr)
  
  cl <- makeCluster(detectCores(logical = F))
  registerDoParallel(cl)
  foreach(r = 1:iterations, .packages = c("mirt", "TestGardener")) %dopar% {
    file_name <- paste("simulation-data/datasets/r", r, "-", scen$name, ".RData", sep = "")
    load(file = file_name)
    
    sim_mirt_mod <- mirt(data = data, model=1, itemtype = "gpcm", SE=F)
    
    dl <- make.dataList(data+1, key = NULL, opt_list, nbin = scen$bins,  
                        NumBasis = scen$basis, Wnorder = scen$spline_order)
    sim_os_mod <- Analyze(dl$percntrnk, dl$thetaQnt, dl, ncycle = os_cycles, verbose = F)
    
    file_name <- paste("simulation-data/models/r", r, "-", scen$name, ".RData", sep = "")
    save(sim_mirt_mod, sim_os_mod, dl, file = file_name)
    NULL
  }
  stopCluster(cl)
}
