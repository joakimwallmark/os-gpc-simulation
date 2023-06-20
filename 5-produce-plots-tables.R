source("1-sim-setup.R")
source("functions/simulation-tables.R")
source("functions/simulation-plots.R")
library(xtable)

print("Producing measure plots and tables")
len <- length(scenarios)
summary <- sum_mirt_items <- sum_os_items <- data.frame(
  character(len),
  numeric(len),
  numeric(len),
  numeric(len),
  numeric(len),
  numeric(len),
  numeric(len),
  numeric(len),
  numeric(len),
  numeric(len),
  numeric(len),
  numeric(len),
  numeric(len),
  numeric(len)
)
colnames(summary) <-
  colnames(sum_mirt_items) <-
  colnames(sum_os_items) <-
  c(
    "pop",
    "n",
    "items",
    "OS item prop.",
    "kl GPC",
    "kl OS",
    "kl SEE GPC",
    "kl SEE GPC",
    "Bias GPC",
    "Bias OS",
    "SEE GPC",
    "SEE OS",
    "RMSE GPC",
    "RMSE OS"
  )
for (sn in seq_along(scenarios)) {
  print(paste("scenario:", sn))
  scen <- scenarios[[sn]]
  file_name <- paste("simulation-data/measures/", scen$name, ".RData", sep = "")
  load(file = file_name)

  x <- scen$true_ability_list$surp_al

  p1_kldi <- sim_fit_plot_kl(x, all_items_measures, 1:2, "Arc length", "KL divergence")
  p1_bias <- sim_fit_plot_bias(x, all_items_measures, "Arc length", "Probability bias")
  p1_rmse <- sim_fit_plot_rmse(x, all_items_measures, "Arc length", "Probability RMSE")
  summary[sn, ] <- c(
    scen$pop,
    scen$n,
    length(scen$os_items) + length(scen$mirt_items),
    scen$prop_os,
    mean(all_items_measures$avg_mirt_kl),
    mean(all_items_measures$avg_os_kl),
    mean(all_items_measures$avg_mirt_kl_see),
    mean(all_items_measures$avg_os_kl_see),
    mean(sapply(all_items_measures$avg_mirt_bias, function(x) rowMeans(abs(x)))),
    mean(sapply(all_items_measures$avg_os_bias, function(x) rowMeans(abs(x)))),
    mean(sapply(all_items_measures$avg_mirt_see, function(x) rowMeans(x))),
    mean(sapply(all_items_measures$avg_os_see, function(x) rowMeans(x))),
    mean(sapply(all_items_measures$avg_mirt_rmse, function(x) rowMeans(x))),
    mean(sapply(all_items_measures$avg_os_rmse, function(x) rowMeans(x)))
  )
  if (is.null(mirt_items_measures) || is.null(os_items_measures)) { # if we only have 1 item type
    file_name <- paste("simulation-data/plots/plot-", scen$name, ".png", sep = "")
    ggsave(file_name,
      plot = p1_kldi / p1_bias / p1_rmse + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
      width = 21 - 2 * 2.54, height = 20, units = "cm"
    )
  } else {
    p1_kldi <- p1_kldi + ggtitle("All items")
    p1_bias <- p1_bias + ggtitle("All items")
    p1_rmse <- p1_rmse + ggtitle("All items")
    p2 <- sim_fit_plot_kl(x, mirt_items_measures, 1:2, "Arc length", "KL divergence", "GPC items")
    p3 <- sim_fit_plot_kl(x, os_items_measures, 1:2, "Arc length", "KL divergence", "OS items")
    file_name <- paste("simulation-data/plots/kl-", scen$name, ".png", sep = "")
    ggsave(file_name,
      plot = p1_kldi / p2 / p3 + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
      width = 21 - 2 * 2.54, height = 20, units = "cm"
    )

    p2 <- sim_fit_plot_bias(x, mirt_items_measures, "Arc length", "Probability bias", "GPC items")
    p3 <- sim_fit_plot_bias(x, os_items_measures, "Arc length", "Probability bias", "OS items")
    file_name <- paste("simulation-data/plots/bias-", scen$name, ".png", sep = "")
    ggsave(file_name,
      plot = p1_bias / p2 / p3 + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
      width = 21 - 2 * 2.54, height = 20, units = "cm"
    )

    p2 <- sim_fit_plot_rmse(x, mirt_items_measures, "Arc length", "Probability RMSE", "GPC items")
    p3 <- sim_fit_plot_rmse(x, os_items_measures, "Arc length", "Probability RMSE", "OS items")
    file_name <- paste("simulation-data/plots/rmse-", scen$name, ".png", sep = "")
    ggsave(file_name,
      plot = p1_rmse / p2 / p3 + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
      width = 21 - 2 * 2.54, height = 20, units = "cm"
    )

    sum_mirt_items[sn, ] <- c(
      scen$pop,
      scen$n,
      length(scen$os_items) + length(scen$mirt_items),
      scen$prop_os,
      mean(mirt_items_measures$avg_mirt_kl),
      mean(mirt_items_measures$avg_os_kl),
      mean(mirt_items_measures$avg_mirt_kl_see),
      mean(mirt_items_measures$avg_os_kl_see),
      mean(sapply(mirt_items_measures$avg_mirt_bias, function(x) rowMeans(abs(x)))),
      mean(sapply(mirt_items_measures$avg_os_bias, function(x) rowMeans(abs(x)))),
      mean(sapply(mirt_items_measures$avg_mirt_see, function(x) rowMeans(x))),
      mean(sapply(mirt_items_measures$avg_os_see, function(x) rowMeans(x))),
      mean(sapply(mirt_items_measures$avg_mirt_rmse, function(x) rowMeans(x))),
      mean(sapply(mirt_items_measures$avg_os_rmse, function(x) rowMeans(x)))
    )

    sum_os_items[sn, ] <- c(
      scen$pop,
      scen$n,
      length(scen$os_items) + length(scen$mirt_items),
      scen$prop_os,
      mean(os_items_measures$avg_mirt_kl),
      mean(os_items_measures$avg_os_kl),
      mean(os_items_measures$avg_mirt_kl_see),
      mean(os_items_measures$avg_os_kl_see),
      mean(sapply(os_items_measures$avg_mirt_bias, function(x) rowMeans(abs(x)))),
      mean(sapply(os_items_measures$avg_os_bias, function(x) rowMeans(abs(x)))),
      mean(sapply(os_items_measures$avg_mirt_see, function(x) rowMeans(x))),
      mean(sapply(os_items_measures$avg_os_see, function(x) rowMeans(x))),
      mean(sapply(os_items_measures$avg_mirt_rmse, function(x) rowMeans(x))),
      mean(sapply(os_items_measures$avg_os_rmse, function(x) rowMeans(x)))
    )
  }
}

summary[2:14] <- apply(summary[2:14], 2, as.numeric)
sum_mirt_items[2:14] <- apply(sum_mirt_items[2:14], 2, as.numeric)
sum_os_items[2:14] <- apply(sum_os_items[2:14], 2, as.numeric)
file_name <- paste("simulation-data/tables/table-", data_string, ".RData", sep = "")
save(summary, sum_mirt_items, sum_os_items, file = file_name)
load(file = file_name)


sim_sum_table(summary)
sim_sum_table(sum_mirt_items)
sim_sum_table(sum_os_items)


# Specific scenario item measures -----------------------------------------
sapply(scenarios, function(x) x$name)
scen_ind <- 15
scen_ind <- 3
scen_ind <- 1
scenarios[[scen_ind]]$name
scenar <- scenarios[[scen_ind]]
file_name <- paste("simulation-data/measures/", scenar$name, ".RData", sep = "")
load(file = file_name)

sim_item_table(all_items_measures)

# Specific item plot ------------------------------------------------------
x <- 100 * scenar$true_ability_list$surp_al / max(scenar$true_ability_list$surp_al)
plot_item <- 9
plot_item <- 17
plot_item <- 24
item_plot <-
  sim_item_fit_plot_bias(x,
    all_items_measures,
    2,
    item = plot_item,
    "% arc length",
    "Probability bias",
    "GPC items"
  )
item_plot
item_plot2 <-
  sim_item_fit_plot_bias_col(x,
    all_items_measures,
    item = plot_item,
    "% arc length",
    "Probability bias",
    "GPC items"
  )
item_plot2


# Average bias over arc length --------------------------------------------
sapply(scenarios, function(x) x$name)
scen_ind <- 25
scen_ind <- 22
scen_ind <- 21
scen_ind <- 16
scen_ind <- 15
scenarios[[scen_ind]]$name
scenar <- scenarios[[scen_ind]]
file_name <- paste("simulation-data/measures/", scenar$name, ".RData", sep = "")
load(file = file_name)
x <- 100 * scenar$true_ability_list$surp_al / max(scenar$true_ability_list$surp_al)

mirt_bias <- rowMeans(sapply(all_items_measures$avg_mirt_bias, function(x) {
  rowMeans(abs(x))
}))
os_bias <- rowMeans(sapply(all_items_measures$avg_os_bias, function(x) {
  rowMeans(abs(x))
}))

mirt_se <- rowMeans(sapply(all_items_measures$avg_mirt_see, function(x) {
  rowMeans(x)
}))
os_se <- rowMeans(sapply(all_items_measures$avg_os_see, function(x) {
  rowMeans(x)
}))

mirt_rmse <- rowMeans(sapply(all_items_measures$avg_mirt_rmse, function(x) {
  rowMeans(x)
}))
os_rmse <- rowMeans(sapply(all_items_measures$avg_os_rmse, function(x) {
  rowMeans(x)
}))

biases <- c(mirt_bias, os_bias)
ses <- c(mirt_se, os_se)
rmses <- c(mirt_rmse, os_rmse)
model <- c(rep("GPC", length(x)), rep("Optimal scoring", length(x)))
df <- data.frame(x = x, biases = biases, ses = ses, rmses = rmses, Model = model)

plot <- ggplot(df, aes(x = x, y = biases, color = Model)) +
  geom_point(size = 0.5) +
  theme_bw() +
  theme(legend.margin = margin(t = -0.2, b = -0.2, unit = "cm")) + # smaller margins legend
  xlab("% Arc length") +
  ylab("Avg. abs. probability bias") +
  scale_x_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0.025, 0.025), add = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, 0.09),
    expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
  ) +
  scale_color_grey(start = 0.1, end = 0.4, name = "")
plot
plot2 <- ggplot(df, aes(x = x, y = ses, color = Model)) +
  geom_point(size = 0.5) +
  theme_bw() +
  theme(legend.margin = margin(t = -0.2, b = -0.2, unit = "cm")) + # smaller margins legend
  xlab("% Arc length") +
  ylab("Avg. probability SE") +
  scale_x_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0.025, 0.025), add = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, 0.09),
    expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
  ) +
  scale_color_grey(start = 0.1, end = 0.4, name = "")
plot2
plot3 <- ggplot(df, aes(x = x, y = rmses, color = Model)) +
  geom_point(size = 0.5) +
  theme_bw() +
  theme(legend.margin = margin(t = -0.2, b = -0.2, unit = "cm")) + # smaller margins legend
  xlab("% Arc length") +
  ylab("Avg. probability RMSE") +
  scale_x_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0.025, 0.025), add = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, 0.09),
    expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
  ) +
  scale_color_grey(start = 0.1, end = 0.4, name = "")
plot3
plot / plot3
plot_norm <- plot
plot_norm3 <- plot3
ggsave(
  paste("plots/sim-item-bias-rmse-", scenarios[[scen_ind]]$name, ".eps", sep = ""),
  plot = plot / plot3,
  width = 21 - 2 * 2.54,
  height = 11.5,
  units = "cm"
)

plot_norm <-
  plot_norm + ylab("Bias") + ggtitle("Dataset population") + theme(
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    axis.title = element_text(size = 9)
  )
plot <-
  plot + ylab("Bias") + ggtitle("Skewed population") + theme(
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    axis.title = element_text(size = 9)
  )
plot_norm3 <-
  plot_norm3 + ylab("RMSE") + ggtitle("Skewed population") + theme(
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    axis.title = element_text(size = 9)
  )
plot3 <-
  plot3 + ylab("RMSE") + ggtitle("Dataset population") + theme(
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    axis.title = element_text(size = 9)
  )
plot_norm <- plot_norm + guides(color = guide_legend(override.aes = list(size = 2)))
plot <- plot + guides(color = guide_legend(override.aes = list(size = 2)))
plot_norm3 <- plot_norm3 + guides(color = guide_legend(override.aes = list(size = 2)))
plot3 <- plot3 + guides(color = guide_legend(override.aes = list(size = 2)))
plot_norm / plot / plot_norm3 / plot3
plot_norm+guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(paste("plots/sim-item-bias-rmse-comp.eps", sep = ""),
  plot = plot_norm / plot / plot_norm3 / plot3, width = 21 - 2 * 2.54, height = 25 - 2 * 2.54, units = "cm"
)
ggsave(paste("plots/sim-item-bias-rmse-comp.jpg", sep = ""),
  plot = plot_norm / plot / plot_norm3 / plot3, width = 21 - 2 * 2.54, height = 25 - 2 * 2.54, units = "cm"
)
