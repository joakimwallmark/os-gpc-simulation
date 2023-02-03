source("functions/arc-length.R")

sim_generate_data <- function(mirt_mod, tg_mod, mirt_thetas, tg_thetas,
                              mirt_items, tg_items) {
  # mirt model data ---------------------------------------------------------
  sample_mirt <-
    simdata(
      model = mirt_mod,
      Theta = matrix(mirt_thetas, ncol = 1),
      which.items = mirt_items,
      equal.K = FALSE
    )
  # testgardener model data -------------------------------------------------
  sample_tg <- Usimulate(tg_thetas, tg_mod$parList[[length(tg_mod$parList)]]$WfdList)
  sample_tg <- sample_tg[, tg_items] # remove items not used
  sample_tg <- sample_tg - 1

  return(cbind(sample_mirt, sample_tg))
}

sim_fit_models <- function(sample, tg_items, tg_bins, tg_basis,
                           tg_opt_score, cyc = 10, opt_scr = NULL) {
  # fit mirt model to data -----------------------------------------
  mirt_mod <- mirt(data = sample, model = 1, itemtype = "gpcm")
  # fit test gardener model to data -----------------------------------------
  dl <- make_dataList(
    U = sample + 1,
    key = NULL,
    optList = list(itemLab = NULL, optLab = NULL, optScr = tg_opt_score),
    NumBasis = tg_basis,
    nbin = tg_bins
  )
  tg_mod <- Analyze(dl$percntrnk,
    dl$thetaQnt,
    dl,
    ncycle = cyc,
    verbose = FALSE
  )
  return(list(mirt_mod, tg_mod))
}
