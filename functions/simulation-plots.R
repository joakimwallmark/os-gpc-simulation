sim_fit_plot_bias <- function(x, measures, cat, x_lab = "Arc length", y_lab = "Bias", title = "") {
  mirt_bias <- measures$avg_mirt_bias[[item]][, cat]
  os_bias <- measures$avg_os_bias[[item]][, cat]

  ys <- c(mirt_bias, os_bias)
  model <- c(rep("GPC", length(x)), rep("Optimal score", length(x)))
  df <- data.frame(x = x, y = ys, Model = model)
  plot <- ggplot(df, aes(x = x, y = y, linetype = Model)) +
    geom_line() +
    theme_bw() +
    theme(legend.margin = margin(t = -0.2, b = -0.2, unit = "cm")) + # smaller margins legend
    xlab(x_lab) +
    ylab(y_lab) +
    scale_x_continuous(
      limits = c(NA, NA),
      expand = expansion(mult = c(0.025, 0.025), add = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(NA, NA),
      expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
    ) +
    scale_linetype_discrete(name = "")

  if (title != "") {
    return(plot + ggtitle(title))
  } else {
    return(plot)
  }
}

sim_item_fit_plot_bias <- function(x, measures, item, cat, x_lab = "Arc length", y_lab = "Bias", title = "") {
  mirt_bias <- measures$avg_mirt_bias[[item]][, cat]
  os_bias <- measures$avg_os_bias[[item]][, cat]

  ys <- c(mirt_bias, os_bias)
  model <- c(rep("GPC", length(x)), rep("Optimal score", length(x)))
  df <- data.frame(x = x, y = ys, Model = model)
  plot <- ggplot(df, aes(x = x, y = y, linetype = Model)) +
    geom_line() +
    theme_bw() +
    theme(legend.margin = margin(t = -0.2, b = -0.2, unit = "cm")) + # smaller margins legend
    xlab(x_lab) +
    ylab(y_lab) +
    scale_x_continuous(
      limits = c(NA, NA),
      expand = expansion(mult = c(0.025, 0.025), add = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(NA, NA),
      expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
    ) +
    scale_linetype_discrete(name = "")

  if (title != "") {
    return(plot + ggtitle(title))
  } else {
    return(plot)
  }
}

sim_item_fit_plot_bias_col <- function(x, measures, item, x_lab = "Arc length", y_lab = "Bias", title = "") {
  mirt_bias <- measures$avg_mirt_bias[[item]]
  os_bias <- measures$avg_os_bias[[item]]

  ys <- c(mirt_bias, os_bias)
  categories <- 0:(ncol(os_bias) - 1)
  score <- as.factor(c(sapply(categories, rep, length(x))))
  model <- c(rep("GPC", length(categories) * length(x)), rep("Optimal score", length(categories) * length(x)))
  df <- data.frame(x = x, y = ys, Model = model, Score = score)
  plot <- ggplot(df, aes(x = x, y = y, linetype = Score, color = Model)) +
    geom_line() +
    theme_bw() +
    theme(legend.margin = margin(t = -0.2, b = -0.2, unit = "cm")) + # smaller margins legend
    xlab(x_lab) +
    ylab(y_lab) +
    scale_x_continuous(
      limits = c(NA, NA),
      expand = expansion(mult = c(0.025, 0.025), add = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(NA, NA),
      expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
    ) +
    scale_linetype_discrete(name = "")

  if (title != "") {
    return(plot + ggtitle(title))
  } else {
    return(plot)
  }
}

sim_fit_plot_rmse <- function(x, measures, x_lab = "Arc length", y_lab = "RMSE", title = "") {
  ys1 <- rowMeans(sapply(measures$avg_mirt_rmse, function(x) rowMeans(x)))
  ys2 <- rowMeans(sapply(measures$avg_os_rmse, function(x) rowMeans(x)))
  ys <- c(ys1, ys2)
  df <- data.frame(x = x, y = ys, model <- c(rep("GPC", length(x)), rep("Optimal score", length(x))))
  plot <- ggplot(df, aes(x = x, y = y, shape = model, linetype = model)) +
    geom_line() +
    theme_bw() +
    theme(legend.margin = margin(t = -0.2, b = -0.2, unit = "cm")) + # smaller margins legend
    xlab(x_lab) +
    ylab(y_lab) +
    scale_x_continuous(
      limits = c(NA, NA),
      expand = expansion(mult = c(0.025, 0.025), add = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(NA, NA),
      expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
    ) +
    scale_linetype_discrete(name = "")

  if (title != "") {
    return(plot + ggtitle(title))
  } else {
    return(plot)
  }
}

sim_fit_plot_bias <- function(x, measures, x_lab = "Arc length", y_lab = "Bias", title = "") {
  ys1 <- rowMeans(sapply(measures$avg_mirt_bias, function(x) rowMeans(abs(x))))
  ys2 <- rowMeans(sapply(measures$avg_os_bias, function(x) rowMeans(abs(x))))
  ys <- c(ys1, ys2)
  df <- data.frame(x = x, y = ys, model <- c(rep("GPC", length(x)), rep("Optimal score", length(x))))
  plot <- ggplot(df, aes(x = x, y = y, shape = model, linetype = model)) +
    geom_line() +
    theme_bw() +
    theme(legend.margin = margin(t = -0.2, b = -0.2, unit = "cm")) + # smaller margins legend
    xlab(x_lab) +
    ylab(y_lab) +
    scale_x_continuous(
      limits = c(NA, NA),
      expand = expansion(mult = c(0.025, 0.025), add = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(NA, NA),
      expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
    ) +
    scale_linetype_discrete(name = "")

  if (title != "") {
    return(plot + ggtitle(title))
  } else {
    return(plot)
  }
}

sim_fit_plot_kl <- function(x, measures, ind, x_lab = "Arc length", y_lab = "KL divergence", title = "") {
  ys <- c(rowMeans(measures[[ind[1]]]), rowMeans(measures[[ind[2]]]))
  df <- data.frame(x = x, y = ys, model <- c(rep("GPC", length(x)), rep("Optimal score", length(x))))
  plot <- ggplot(df, aes(x = x, y = y, shape = model, linetype = model)) +
    geom_line() +
    theme_bw() +
    theme(legend.margin = margin(t = -0.2, b = -0.2, unit = "cm")) + # smaller margins legend
    xlab(x_lab) +
    ylab(y_lab) +
    scale_x_continuous(
      limits = c(NA, NA),
      expand = expansion(mult = c(0.025, 0.025), add = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(NA, NA),
      expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
    ) +
    scale_linetype_discrete(name = "")

  if (title != "") {
    return(plot + ggtitle(title))
  } else {
    return(plot)
  }
}



# takes simulation tables and makes plots
plots_mse_kl_cdf_perc <-
  function(scen_name,
           x,
           mirt_rmse,
           tg_rmse,
           mirt_kl,
           tg_kl,
           file_name_mse,
           file_name_kl) {
    df_mse <-
      data.frame(
        x = x,
        y = c(rowMeans(mirt_rmse), rowMeans(tg_rmse)),
        model <- c(rep("GPC", length(x)), rep("Optimal score", length(x)))
      )
    mse_plot <- ggplot(df_mse, aes(x = x, y = y, linetype = model)) +
      geom_line() +
      theme_bw() +
      ggtitle("Average probability RMSE") +
      xlab(expression(theta ~ "distribution cdf percentage")) +
      ylab("RMSE") +
      scale_x_continuous(
        limits = c(NA, NA),
        expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
      ) +
      scale_y_continuous(
        limits = c(NA, NA),
        expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
      )
    df_mse_mirt_items <-
      data.frame(
        x = x,
        y = c(
          rowMeans(mirt_rmse[, seq_along(scen$mirt_items)]),
          rowMeans(tg_rmse[, seq_along(scen$mirt_items)])
        ),
        model <- c(rep("GPC", length(x)), rep("Optimal score", length(x)))
      )
    mse_plot_mirt <- ggplot(df_mse_mirt_items, aes(x = x, y = y, linetype = model)) +
      geom_line() +
      theme_bw() +
      ggtitle("Average probability RMSE. GPC items") +
      xlab(expression(theta ~ "distribution cdf percentage")) +
      ylab("RMSE") +
      scale_x_continuous(
        limits = c(NA, NA),
        expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
      ) +
      scale_y_continuous(
        limits = c(NA, NA),
        expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
      )
    df_mse_tg_items <-
      data.frame(
        x = x,
        y = c(
          rowMeans(mirt_rmse[, -(seq_along(scen$mirt_items))]),
          rowMeans(tg_rmse[, -(seq_along(scen$mirt_items))])
        ),
        model <-
          c(rep("GPC", length(x)), rep("Optimal score", length(x)))
      )
    mse_plot_tg <- ggplot(df_mse_tg_items, aes(x = x, y = y, linetype = model)) +
      geom_line() +
      theme_bw() +
      ggtitle("Average probability RMSE. OS items") +
      xlab(expression(theta ~ "distribution cdf percentage")) +
      ylab("RMSE") +
      scale_x_continuous(
        limits = c(NA, NA),
        expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
      ) +
      scale_y_continuous(
        limits = c(NA, NA),
        expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
      )
    df_kl <-
      data.frame(
        x = x,
        y = c(rowMeans(mirt_kl), rowMeans(tg_kl)),
        model <- c(rep("GPC", length(x)), rep("Optimal score", length(x)))
      )
    kl_plot <- ggplot(df_kl, aes(x = x, y = y, linetype = model)) +
      geom_line() +
      theme_bw() +
      ggtitle("Average KL diverence") +
      xlab(expression(theta ~ "distribution cdf percentage")) +
      ylab("KL divergence") +
      scale_x_continuous(
        limits = c(NA, NA),
        expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
      ) +
      scale_y_continuous(
        limits = c(NA, NA),
        expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
      )

    df_kl_mirt_items <-
      data.frame(
        x = x,
        y = c(rowMeans(mirt_kl[, seq_along(scen$mirt_items)]), rowMeans(tg_kl[, seq_along(scen$mirt_items)])),
        model <- c(rep("GPC", length(x)), rep("Optimal score", length(x)))
      )
    kl_plot_mirt <- ggplot(df_kl_mirt_items, aes(x = x, y = y, linetype = model)) +
      geom_line() +
      theme_bw() +
      ggtitle("Average KL diverence. GPC items") +
      xlab(expression(theta ~ "distribution cdf percentage")) +
      ylab("KL divergence") +
      scale_x_continuous(
        limits = c(NA, NA),
        expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
      ) +
      scale_y_continuous(
        limits = c(NA, NA),
        expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
      )

    df_kl_tg_items <-
      data.frame(
        x = x,
        y = c(
          rowMeans(mirt_kl[, -(seq_along(scen$mirt_items))]),
          rowMeans(tg_kl[, -(seq_along(scen$mirt_items))])
        ),
        model <-
          c(rep("GPC", length(x)), rep("Optimal score", length(x)))
      )
    kl_plot_tg <- ggplot(df_kl_tg_items, aes(x = x, y = y, linetype = model)) +
      geom_line() +
      theme_bw() +
      ggtitle("Average KL diverence. OS items") +
      xlab(expression(theta ~ "distribution cdf percentage")) +
      ylab("KL divergence") +
      scale_x_continuous(
        limits = c(NA, NA),
        expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
      ) +
      scale_y_continuous(
        limits = c(NA, NA),
        expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
      )
    ggsave(
      file_name_mse,
      plot = mse_plot / mse_plot_mirt / mse_plot_tg + plot_annotation(tag_levels = "A"),
      width = 21 - 2 * 2.54,
      height = 20,
      units = "cm"
    )
    ggsave(
      file_name_kl,
      plot = kl_plot / kl_plot_mirt / kl_plot_tg + plot_annotation(tag_levels = "A"),
      width = 21 - 2 * 2.54,
      height = 20,
      units = "cm"
    )
    return(list(mse_plot, mse_plot_mirt, mse_plot_tg, kl_plot, kl_plot_mirt, kl_plot_tg))
  }
