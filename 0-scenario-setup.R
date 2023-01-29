library(patchwork)
library(mirt)
library(TestGardener)
library(sn)
source("functions/tg_density.R")
source("functions/arc_length.R")
source("functions/get_abs_distance.R")

# skew normal parameters
xi <- 1.5 # location
omega <- 1.2 # scale (~variance)
alpha <- -4 # shape
max_surp <- 20 # needs to match sim-setup.R file...

# NATMAT ------------------------------------------------------------------
data_string <- "natmat18"
test_lengths <- c(14, 28)
props_os_items <- c(1, 0.5, 0)
nbins <- 10 # Match original model
nbasis <- 6
spline_order <- 4
# sample test takers ------------------------------------------------------
set.seed(12335)
test_takers <- list()
test_takers[[1]] <- sort(rnorm(1000))
test_takers[[2]] <- sort(rsn(1000, xi, omega, alpha))
test_takers[[3]] <- sort(rnorm(5000))
test_takers[[4]] <- sort(rsn(5000, xi, omega, alpha))
test_takers_pop <- c("normal", "skewed", "normal", "skewed")

# SWESAT ------------------------------------------------------------------
data_string <- "swesat14"
test_lengths <- c(30, 60)
props_os_items <- c(1, 0.5, 0)
nbins <- 7 # Match original model
nbasis <- 4
spline_order <- 4
# sample test takers ------------------------------------------------------
set.seed(12335)
test_takers <- list()
test_takers[[1]] <- sort(rnorm(1000))
test_takers[[2]] <- sort(rsn(1000, xi, omega, alpha))
test_takers[[3]] <- sort(rnorm(5000))
test_takers[[4]] <- sort(rsn(5000, xi, omega, alpha))
test_takers_pop <- c("normal", "skewed", "normal", "skewed", "normal", "skewed")

# Load models -------------------------------------------------------------
load(paste("simulation-data/sim_", data_string, "_mirt_mod.RData", sep = ""))
load(paste("simulation-data/sim_", data_string, "_tg_mod.RData", sep = ""))
load(paste("simulation-data/sim_", data_string, "_item_data.RData", sep = ""))

# OS population theta density ---------------------------------------------
set.seed(131) # donno why density spline is random...
data_dens <- tg_theta_density(c(os_mod$parList[[10]]$theta), alpha = 1.5)
data_dens$plot

# SETUP SCENARIOS ---------------------------------------------------------
items <- seq_along(response_cats)
resp_cats <- sort(as.numeric(unique(names(table(response_cats)))), decreasing = TRUE)
scenarios <- list() # every unique scenario is stored in this list
scenario_ind <- 1
set.seed(2422)
for (length in test_lengths) {
  for (prop_os_items in props_os_items) {
    # sample items
    if (prop_os_items == 1) {
      item_sample <- sample(response_cats, length)
      os_items <- which(names(response_cats) %in% names(item_sample))
      mirt_items <- c()
    } else if (prop_os_items == 0) {
      os_items <- c()
      item_sample <- sample(response_cats, length)
      mirt_items <- which(names(response_cats) %in% names(item_sample))
    } else {
      item_sample <- sample(response_cats, length)
      os_sample <- sample(seq_along(item_sample), prop_os_items * length)
      os_items <- which(names(response_cats) %in% names(item_sample[os_sample]))
      mirt_items <- which(names(response_cats) %in% names(item_sample[-os_sample]))
    }
    for (sample in seq_along(test_takers)) {
      mirt_samp <- test_takers[[sample]]
      n <- length(mirt_samp)
      # convert to corresponding tg thetas
      sample_tg <- tg_qss(pnorm(mirt_samp), data_dens$density)

      # compute different measures of ability...
      true_ability_list <- list()
      if (is.null(mirt_items)) {
        total_cat <- sum(response_cats[c(os_items)])
        true_ability_list$surp_al <-
          arc_length_surp(
            os_mod,
            sample_tg,
            total_cat,
            max_surp,
            item_ids = os_items,
            mesh_accuracy = 10001
          )[, 1]
      } else if (is.null(os_items)) {
        total_cat <- sum(response_cats[c(mirt_items)])
        true_ability_list$surp_al <-
          arc_length_surp(
            mirt_mod,
            mirt_samp,
            total_cat,
            max_surp,
            item_ids = mirt_items,
            mesh_accuracy = 10001
          )[, 1]
      } else {
        true_ability_list$surp_al <-
          arc_length_surp_mixed(
            mirt_mod,
            mirt_samp,
            mirt_items,
            os_mod,
            sample_tg,
            os_items,
            total_cat,
            100,
            mesh_accuracy = 101
          )[, 1]
      }

      surp_rand_mirt <-
        get_dist_travelled_mirt(mirt_mod, mirt_items, sort(c(
          mirt_samp, seq(min(mirt_samp), max(mirt_samp), length.out = 10000)
        )))
      surp_rand_tg <-
        get_dist_travelled_tg(os_mod, os_items, sort(c(sample_tg, seq(0, 100, length.out = 10000))))

      prob_rand_mirt <-
        get_dist_travelled_mirt(mirt_mod, mirt_items, sort(c(
          mirt_samp, seq(min(mirt_samp), max(mirt_samp), length.out = 10000)
        )), mode = "prob")
      prob_rand_tg <-
        get_dist_travelled_tg(os_mod, os_items, sort(c(sample_tg, seq(0, 100, length.out = 10000))), mode = "prob")

      entr_rand_mirt <-
        get_dist_travelled_mirt(mirt_mod, mirt_items, sort(c(
          mirt_samp, seq(min(mirt_samp), max(mirt_samp), length.out = 10000)
        )), mode = "entropy")
      entr_rand_tg <-
        get_dist_travelled_tg(os_mod, os_items, sort(c(sample_tg, seq(0, 100, length.out = 10000))), mode = "entropy")
      # pick out the thetas for the sample. Add same person thetas from both models
      true_ability_list$surp_dist <-
        surp_rand_tg[match(sample_tg, surp_rand_tg[, 1]), 2] + surp_rand_mirt[match(mirt_samp, surp_rand_mirt[, 1]), 2]
      true_ability_list$prob_dist <-
        prob_rand_tg[match(sample_tg, prob_rand_tg[, 1]), 2] + prob_rand_mirt[match(mirt_samp, prob_rand_mirt[, 1]), 2]
      true_ability_list$entr_dist <-
        entr_rand_tg[match(sample_tg, entr_rand_tg[, 1]), 2] + entr_rand_mirt[match(mirt_samp, entr_rand_mirt[, 1]), 2]


      scnario_name <-
        paste(
          data_string,
          "-n",
          n,
          "-i",
          length,
          "-os-prop",
          round(prop_os_items, 2),
          "-",
          test_takers_pop[sample],
          sep = ""
        )
      scenarios[[scenario_ind]] <- list(
        name = scnario_name,
        pop = test_takers_pop[sample],
        n = n,
        prop_os = prop_os_items,
        os_items = os_items,
        mirt_items = mirt_items,
        tg_mod = data_dens$density,
        test_taker_thetas = cbind(mirt = mirt_samp, tg = sample_tg),
        true_ability_list = true_ability_list,
        bins = nbins,
        basis = nbasis,
        spline_order = spline_order
      )
      print(paste("scenario", scenario_ind))
      scenario_ind <- scenario_ind + 1
    }
  }
}

save(scenarios, file = paste("simulation-data/sim-", data_string, "-scenarios.RData", sep = ""))
load(file = paste("simulation-data/sim-", data_string, "-scenarios.RData", sep = ""))
