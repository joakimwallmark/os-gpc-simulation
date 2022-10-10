options(scipen = 999)
# source("functions/fit_tg_models.R")
# source("functions/mod_item_fits.R")
# source("functions/mean_fit_table.R")
# source("functions/latex_mean_table.R")
library(patchwork)
library(ggplot2)

source("functions/sim_functions.R")
library(TestGardener)
library(mirt)

load("swesat14_mirt_mod.RData")
load("swesat14_tg_mod.RData")

# get thetas for both models
# mirt_thetas_ml <- c(fscores(mirt_mod, method = "ML"))
mirt_thetas <- c(fscores(mirt_mod))
tg_thetas <- c(best_tg_mod$parList[[10]]$theta)

# remove actual data
# load("swesat14_mirt_mod.RData")
# mirt_mod@Data$data[mirt_mod@Data$data > 0] <- 0
# save(mirt_thetas, mirt_thetas_ml, mirt_mod, file = "swesat14_mirt_mod.RData")

mirt_item_sample_ids <- sample(1:90, 45)
mirt_items <- colnames(mirt_mod@Data$data[, mirt_item_sample_ids])
tg_items <- (1:90)[!1:90 %in% mirt_item_sample_ids]

sim_generate_data(mirt_mod, best_tg_mod, 
                  c(0, 0, 2, 2), c(10, 10, 90, 90), 
                  mirt_item_sample_ids, tg_items)
sim_generate_data(mirt_mod, best_tg_mod, 
                  mirt_thetas, tg_thetas, 
                  mirt_item_sample_ids, tg_items)

sample_mirt <-
  simdata(
    model = mirt_mod,
    Theta = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), ncol = 1)
  )
