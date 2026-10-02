# Replication of the main-text figures of the GP-MF paper.
# Run with replication/ as working directory: Rscript out_script.R
#
# Raw mode starts from the collected source files and writes the intermediate
# CSVs in data/; CSV mode produces the same figures from these CSVs alone.
# Both modes read the CSVs, so all figures come from the same code path.

library(dplyr)
library(tidyr)
library(reshape2)
library(lubridate)

library(ggplot2)
library(ggh4x)
library(scales)
library(cowplot)

## --------- Settings ---------------------------- #
use.raw <- NA # NA = use raw files if both exist, TRUE = require them, FALSE = CSV mode
raw.eval <- "raw/eval_out.rda"
raw.mcs <- "raw/mcs_out.rda"

data.dir <- "data/"
save.dir <- "figures/"

has.raw <- file.exists(raw.eval) && file.exists(raw.mcs)
if (is.na(use.raw)) use.raw <- has.raw
if (use.raw && !has.raw) stop("raw mode requested, but ", raw.eval, " or ", raw.mcs, " is missing")
cat("Mode:", ifelse(use.raw, "raw", "CSV"), "\n")

# colours and chart settings
col_blue <- "navy"
col_red <- "firebrick"
tile_alpha <- 0.2

alpha_bars <- 1
col_ind_best <- "black"
lwd_bar_overall <- 0.2

mean_ord <- c("GP", "BLR", "BART")

outsig <- function(x, digs, lab) {
  paste0(format(round(x, digs), nsmall = digs), lab)
}

sl.endo <- c("GDPC1", "GDPCTPI")

# dates (labelled by the first day of the last month of the quarter)
sl.dates <- as.character(seq(1990, 2023.25, by = 0.25))
datelabs_q <- as.character(seq(as.Date("1990-03-01"), as.Date("2023-06-01"), by = "quarter"))
datelabs_m <- as.character(seq(as.Date("1990-01-01"), as.Date("2023-09-01"), by = "month"))

sl.precovid <- as.character(seq(as.Date("1980-01-01"), as.Date("2019-12-01"), by = "month"))
sl.nocovid <- c(
  as.character(seq(as.Date("1980-01-01"), as.Date("2020-03-01"), by = "month")),
  as.character(seq(as.Date("2020-10-01"), as.Date("2023-09-01"), by = "month"))
)
sl.postcovid <- datelabs_m[!(datelabs_m %in% sl.precovid)]

sl.dates.subs <- list(
  "Full" = datelabs_m,
  "Full (*)" = datelabs_m[datelabs_m %in% sl.nocovid],
  "Pre 2020" = sl.precovid,
  "Post 2020 (*)" = sl.postcovid[sl.postcovid %in% sl.nocovid]
)
date_sets <- names(sl.dates.subs)

# labels (raw -> display)
midas_raw <- c(
  "br", "xalm",
  "alm3", "ber3", "leg3", "fou3", "ber5", "leg5",
  "avg", "u", "f"
)
midas_lab <- c(
  "eq", "xalm",
  "alm3", "ber3", "leg3", "fou3", "ber5", "leg5",
  "br", "blu", "svd"
)

mean_raw <- c("lin-arp", "lin-mfvar", "lin-hs", "bart-250", "gp-sqexp")
mean_lab <- c("AR(P)", "MFVAR", "BLR", "BART", "GP")

size_raw <- c("s", "f", "b")
size_lab <- c("s", "f", "b")

month_d_lab <- c("5/3", "4/3", "1", "2/3", "1/3", "0", "-1/3")
mcs_splits <- c("full", "full-outl", "pre2020", "post2020", "post2020-outl")
splits_lab_map <- c(
  "Full" = "full", "Full (*)" = "full-outl",
  "Pre 2020" = "pre2020", "Post 2020 (*)" = "post2020-outl"
)

bench.mdl <- c(
  "size" = "s",
  "midas" = "br",
  "mean" = "BLR",
  "var" = "hom"
) # Benchmark model

# writes doubles with 17 significant digits (base write.csv keeps only 15)
write_csv17 <- function(df, file) {
  num <- vapply(df, is.double, logical(1))
  df[num] <- lapply(df[num], function(x) sprintf("%.17g", x))
  write.csv(df, file = file, row.names = FALSE, quote = which(!num))
}

dir.create(data.dir, showWarnings = FALSE)
dir.create(save.dir, showWarnings = FALSE)

## --------- Raw branch: source files -> intermediate CSVs ---------------------------- #
if (use.raw) {
  load(raw.eval) # eval.obj, missings

  MSE <- melt(eval.obj$MSE[, sl.dates, , ]) %>%
    mutate("metric" = "mse") %>%
    select(metric, Var1, Var2, Var3, Var4, value)
  CRPS <- melt(eval.obj$QWCRPS["none", , sl.dates, , ]) %>%
    mutate("metric" = "none") %>%
    select(metric, Var1, Var2, Var3, Var4, value)
  rm(eval.obj, missings)
  gc()

  colnames(MSE) <- colnames(CRPS) <- c("metric", "target", "date", "horizon", "model", "value")
  metrics_df <- bind_rows(MSE, CRPS) %>%
    mutate(date = as.factor(date))
  rm(MSE, CRPS)
  levels(metrics_df$date) <- datelabs_q
  levels(metrics_df$model) <- gsub("fou_", "fou", gsub("ber_", "ber", gsub(
    "leg_", "leg",
    gsub("_alm_", "_alm", levels(metrics_df$model))
  )))

  metrics_df <- metrics_df %>%
    mutate(date = as.character(date)) %>%
    separate(horizon, into = c("horizon", "month")) %>%
    mutate(
      horizon = as.numeric(ifelse(horizon == "nc", "0", "1")),
      month = as.numeric(gsub("m", "", month))
    ) %>%
    separate(model, into = c("size", "midas", "mean", "var"), sep = "_")

  # forecasts (horizon 1) are dated by the target quarter
  metrics_h0 <- metrics_df %>%
    subset(horizon %in% 0) %>%
    mutate(
      mbq = month,
      month = (month - 3) / 3
    )
  metrics_h1 <- metrics_df %>%
    subset(horizon %in% 1) %>%
    mutate(date = as.character(as.Date(date) %m+% months(3))) %>%
    mutate(
      mbq = month,
      month = ((month - 6)) / 3
    )
  metrics_df <- bind_rows(metrics_h0, metrics_h1)
  rm(metrics_h0, metrics_h1)

  metrics_df$var <- factor(metrics_df$var)
  levels(metrics_df$var) <- c("hom", "sv")
  metrics_df$metric <- factor(metrics_df$metric, levels = c("none", "mse"))
  levels(metrics_df$metric) <- c("CRPS", "MSE")
  metrics_df$size <- factor(metrics_df$size, levels = size_raw)
  levels(metrics_df$size) <- size_lab
  metrics_df$midas <- factor(metrics_df$midas, levels = midas_raw)
  levels(metrics_df$midas) <- midas_lab
  metrics_df$mean <- factor(metrics_df$mean, levels = mean_raw)
  levels(metrics_df$mean) <- mean_lab

  # attach the benchmark
  metrics_bench <- metrics_df %>%
    subset(size %in% bench.mdl["size"]) %>%
    subset(midas %in% bench.mdl["midas"]) %>%
    subset(mean %in% bench.mdl["mean"]) %>%
    subset(var %in% bench.mdl["var"]) %>%
    rename("bench" = "value") %>%
    ungroup() %>%
    select(metric, target, date, horizon, month, bench)
  metrics_df <- left_join(metrics_df,
    metrics_bench,
    by = c("metric", "target", "date", "horizon", "month")
  ) %>%
    select(metric, target, date, horizon, mbq, month, size, midas, mean, var, value, bench)
  rm(metrics_bench)

  metrics_df$horizon <- factor(metrics_df$horizon)
  metrics_df$mbq <- factor(metrics_df$mbq)

  # subsample averages (RMSE = square root of the mean squared error)
  split_ls <- list()
  for (subs in date_sets) {
    metrics_sub <- metrics_df %>% subset(date %in% sl.dates.subs[[subs]])

    mse_df <- metrics_sub %>%
      ungroup() %>%
      subset(metric %in% "MSE") %>%
      group_by(metric, target, horizon, mbq, month, size, midas, mean, var) %>%
      summarise(
        value = sqrt(mean(value)),
        bench = sqrt(mean(bench)),
        .groups = "drop"
      )
    other_df <- metrics_sub %>%
      ungroup() %>%
      subset(!(metric %in% "MSE")) %>%
      group_by(metric, target, horizon, mbq, month, size, midas, mean, var) %>%
      summarise(
        value = mean(value),
        bench = mean(bench),
        .groups = "drop"
      )
    metrics_avg <- bind_rows(other_df, mse_df)
    rm(mse_df, other_df)
    levels(metrics_avg$metric)[levels(metrics_avg$metric) == "MSE"] <- "RMSE"

    split_ls[[subs]] <- metrics_avg %>%
      mutate(split = subs)
  }
  rm(metrics_df, metrics_sub, metrics_avg)

  losses_avg <- do.call(bind_rows, split_ls) %>%
    subset(mean != "AR(P)")
  rm(split_ls)
  losses_avg$month_d <- factor(losses_avg$month)
  levels(losses_avg$month_d) <- month_d_lab

  losses_avg <- losses_avg %>%
    mutate(across(c(metric, target, horizon, mbq, month_d, size, midas, mean, var), as.character)) %>%
    select(metric, target, split, horizon, mbq, month_d, size, midas, mean, var, value, bench)
  write_csv17(as.data.frame(losses_avg), paste0(data.dir, "losses_avg.csv"))
  rm(losses_avg)

  # model confidence set indicators
  load(raw.mcs) # mcs_df
  levels(mcs_df$month_d) <- paste0(gsub("h", "", levels(mcs_df$month_d)), "/3")
  levels(mcs_df$month_d)[levels(mcs_df$month_d) == "3/3"] <- "1"
  levels(mcs_df$month_d)[levels(mcs_df$month_d) == "0/3"] <- "0"
  levels(mcs_df$metric)[levels(mcs_df$metric) == "MSE"] <- "RMSE"

  mcs_out <- mcs_df %>%
    subset(metric %in% "CRPS") %>%
    subset(split %in% splits_lab_map) %>%
    mutate(across(c(metric, target, month_d, split), as.character))
  write_csv17(as.data.frame(mcs_out), paste0(data.dir, "mcs.csv"))
  rm(mcs_df, mcs_out)
}

## --------- Common entry point: intermediate CSVs ---------------------------- #
splits_out <- read.csv(paste0(data.dir, "losses_avg.csv"),
  colClasses = c(rep("character", 10), "numeric", "numeric")
)
splits_out <- splits_out %>%
  mutate(
    metric = factor(metric, levels = c("CRPS", "RMSE")),
    target = factor(target, levels = sl.endo),
    split = factor(split, levels = date_sets),
    horizon = factor(horizon, levels = c("0", "1")),
    mbq = factor(mbq, levels = c("1", "2", "3", "4")),
    month_d = factor(month_d, levels = month_d_lab),
    size = factor(size, levels = size_lab),
    midas = factor(midas, levels = midas_lab),
    mean = factor(mean, levels = setdiff(mean_lab, "AR(P)")),
    var = factor(var, levels = c("hom", "sv")),
    rel = value / bench,
    pred_type = factor(ifelse(mbq == "4", "BC", ifelse(horizon == "0", "NC", "FC")),
      levels = c("BC", "NC", "FC")
    )
  )

mcs_df <- read.csv(paste0(data.dir, "mcs.csv"),
  colClasses = c(rep("character", 8), rep("numeric", 3))
)
mcs_df <- mcs_df %>%
  mutate(
    metric = factor(metric, levels = c("CRPS", "RMSE")),
    target = factor(target, levels = sl.endo),
    month_d = factor(month_d, levels = month_d_lab),
    split = factor(split, levels = mcs_splits)
  )

## --------- Best-model tables (CRPS) ---------------------------- #
set_digits <- 2
set.max <- 1 # maximum number of top models
set.max_gp <- 1 # maximum number of top GP model

row_h <- 0.25 # inches per tile row
overhead <- 0.75 # fixed overhead per plot (title, strips, axes, margins)
min_plot_h <- 1.4 # minimum plot height to avoid squishing low-row panels

splits_sets <- list("full" = c("Full", "Full (*)", "Pre 2020", "Post 2020 (*)"))

horz <- levels(splits_out$month_d)

plot.dir <- "tabs"
dir.create(paste0(save.dir, plot.dir), showWarnings = FALSE)

mcs_label <- function(rel, mcs, digs) {
  ifelse(is.na(mcs) | mcs == 0, outsig(rel, digs, " "),
    ifelse(mcs == 1, outsig(rel, digs, "'"),
      ifelse(mcs == 2, outsig(rel, digs, "°"), NA)))
}

for (tt in sl.endo) {
  for (mm in c("CRPS")) {
    for (ii in 1:length(splits_sets)) {
      splits <- splits_sets[[ii]]
      pp_ls <- list()
      n_rows_ls <- list()
      for (ss in splits) {
        splits_lab <- splits_lab_map[[ss]]

        # best GP and best non-GP specification per horizon
        best_vec <- NULL
        best_multi <- NULL
        for (hh in horz) {
          sl_best <- splits_out %>%
            ungroup() %>%
            subset(target %in% tt) %>%
            subset(split %in% ss) %>%
            subset(metric %in% mm) %>%
            subset(month_d %in% hh) %>%
            mutate(uid = paste0(mean, "-", var, "-", midas, "-", size)) %>%
            ungroup() %>%
            group_by(mean) %>%
            arrange(value) %>%
            ungroup()
          sl_best1 <- sl_best %>%
            subset(mean %in% "GP") %>%
            group_by(mean) %>%
            slice_head(n = set.max_gp)
          sl_best2 <- sl_best %>%
            subset(mean != "GP") %>%
            slice_head(n = set.max)
          sl_best3 <- sl_best %>%
            subset(mean == "MFVAR") %>%
            slice_head(n = set.max)
          sl_best <- bind_rows(sl_best1, sl_best2)
          best_vec <- c(best_vec, sl_best$uid)
          best_multi <- c(best_multi, sl_best3$uid)
        }

        # select best models across horizons
        sl_best <- unique(best_vec)
        splits_best <- splits_out %>%
          subset(target %in% tt) %>%
          subset(split %in% ss) %>%
          subset(metric %in% mm) %>%
          mutate(uid = paste0(mean, "-", var, "-", midas, "-", size)) %>%
          subset(uid %in% sl_best)

        splits_best <- splits_best %>%
          mutate(split_merge = splits_lab) %>%
          left_join(mcs_df, by = c("metric", "target", "month_d", "mean", "var", "midas", "size", "split_merge" = "split"))

        # one row with the best horizon-specific loss of the multivariate model
        sl_best_multi <- unique(best_multi)
        splits_best_multi <- splits_out %>%
          subset(target %in% tt) %>%
          subset(split %in% ss) %>%
          subset(metric %in% mm) %>%
          mutate(uid = paste0(mean, "-", var, "-", midas, "-", size)) %>%
          subset(uid %in% sl_best_multi) %>%
          mutate(split_merge = splits_lab) %>%
          left_join(mcs_df, by = c("metric", "target", "month_d", "mean", "var", "midas", "size", "split_merge" = "split")) %>%
          group_by(month_d) %>%
          filter(value == min(value)) %>%
          # MCS indication
          mutate(value_print = mcs_label(rel, mcs, set_digits)) %>%
          mutate(best = 0, value_bold = "", uid2 = "best per h")
        splits_best_multi$mean <- "MV"
        splits_best_multi$month_d <- factor(splits_best_multi$month_d, levels = rev(levels(splits_best_multi$month_d)))

        # adjust with indicators for plotting
        tab_out <- splits_best %>%
          # MCS indication
          mutate(value_print = mcs_label(rel, mcs, set_digits)) %>%
          ungroup() %>%
          group_by(metric, target, split, month_d) %>%
          mutate(best = ifelse(value == min(value), 1, 0)) %>%
          mutate(
            value_bold = ifelse(best == 1, value_print, ""),
            value_print = ifelse(best == 0, value_print, "")
          )

        tab_out$month_d <- factor(tab_out$month_d, levels = rev(levels(tab_out$month_d)))

        # reshape table for plotting
        pp_tmp <- tab_out %>%
          mutate(uid2 = paste0(var, "-", midas, "-", size)) %>%
          subset(mean != "MFVAR")
        pp_tmp <- bind_rows(pp_tmp, splits_best_multi)
        n_rows_ls[[ss]] <- pp_tmp %>%
          group_by(mean) %>%
          summarise(n = n_distinct(uid2), .groups = "drop") %>%
          pull(n) %>%
          sum()

        pp_ls[[ss]] <- pp_tmp %>%
          ggplot() +
          geom_tile(aes(x = month_d, y = uid2, fill = rel), alpha = tile_alpha) +
          geom_text(aes(x = month_d, y = uid2, label = value_bold), size = 3, fontface = "bold") +
          geom_text(aes(x = month_d, y = uid2, label = value_print), size = 3) +
          labs(title = ss, x = "h", y = "") +
          facet_nested(mean ~ pred_type,
            switch = "y", scales = "free", space = "free"
          ) +
          scale_fill_gradient2(
            midpoint = 1, low = col_blue, high = col_red, mid = "grey98", na.value = "grey80",
            limits = c(0.8, 1.2), oob = squish
          ) +
          scale_color_gradient2(
            midpoint = 1, low = col_blue, high = col_red, mid = "grey98", na.value = "grey80",
            limits = c(0.8, 1.2), oob = squish
          ) +
          coord_cartesian(expand = FALSE) +
          theme_cowplot() +
          theme(
            strip.placement = "outside",
            strip.background = element_blank(),
            strip.text.y.left = element_text(angle = 0, hjust = 0),
            legend.position = "none", legend.key.width = unit(1.5, "cm"), legend.key.height = unit(0.1, "cm"),
            axis.ticks = element_blank(),
            axis.text.y = element_text(hjust = 1),
            axis.line = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0, size = 14, margin = margin(l = 8, unit = "mm")),
            plot.title.position = "plot",
            plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
          )
      }
      n_rows_vec <- unlist(n_rows_ls[splits])
      plot_heights <- pmax(n_rows_vec * row_h + overhead, min_plot_h)
      idx1 <- seq(1, length(splits), by = 2)
      idx2 <- seq(2, length(splits), by = 2)
      col1_height <- sum(plot_heights[idx1])
      col2_height <- sum(plot_heights[idx2])
      col1_plot <- plot_grid(plotlist = pp_ls[splits[idx1]], ncol = 1,
                             align = "v", axis = "tb", rel_heights = plot_heights[idx1])
      col2_plot <- plot_grid(plotlist = pp_ls[splits[idx2]], ncol = 1,
                             align = "v", axis = "tb", rel_heights = plot_heights[idx2])
      if (col1_height < col2_height) {
        col1_plot <- plot_grid(col1_plot, NULL, ncol = 1,
                               rel_heights = c(col1_height, col2_height - col1_height))
      } else if (col2_height < col1_height) {
        col2_plot <- plot_grid(col2_plot, NULL, ncol = 1,
                               rel_heights = c(col2_height, col1_height - col2_height))
      }
      pdf_height <- max(col1_height, col2_height)
      pdf(file = paste0(save.dir, plot.dir, "/", "bestmods_", tt, "_", mm, "_", names(splits_sets)[ii], ".pdf"), width = 10, height = pdf_height)
      print(plot_grid(col1_plot, col2_plot, ncol = 2))
      dev.off()
    }
  }
}

## --------- Mean comparison (levels, Full (*)) ---------------------------- #
plot.dir <- "expl_lvl"
dir.create(paste0(save.dir, plot.dir, "/mean"), recursive = TRUE, showWarnings = FALSE)

plot_vv <- splits_out %>%
  ungroup() %>%
  mutate(plot_value = value)

for (tt in sl.endo) {
  for (mm in c("CRPS", "RMSE")) {
    for (ss in c("Full (*)")) {
      plot_out <- plot_vv %>%
        subset(target %in% tt) %>%
        subset(split %in% ss) %>%
        subset(metric %in% mm) %>%
        subset(mean != "MFVAR")

      plot_out$mean <- factor(plot_out$mean, levels = mean_ord)
      y_breaks <- seq(floor(10 * min(plot_out$plot_value)) / 10, ceiling(10 * max(plot_out$plot_value)) / 10, length.out = 6)

      ss_lab <- gsub(" ", "", gsub("\\)", "", gsub("\\*", "_outl", gsub(" \\(", "", ss))))
      plot_out$month_d <- factor(plot_out$month_d, levels = rev(levels(plot_out$month_d)))

      pp_means <- plot_out %>%
        ggplot(aes(x = mean, y = plot_value)) +
        geom_hline(yintercept = y_breaks, color = "grey90", linewidth = 0.5) +
        # full range (min–max)
        stat_summary(
          fun.data = function(x) data.frame(ymin = min(x), ymax = max(x), y = median(x)),
          geom = "crossbar", width = 0.8, linewidth = 0.2, alpha = alpha_bars, fill = "white"
        ) +
        # inner range (10th–90th percentile)
        stat_summary(
          fun.data = function(x) data.frame(
            ymin = quantile(x, 0.1), ymax = quantile(x, 0.9), y = median(x)
          ),
          geom = "crossbar", width = 0.8, linewidth = 0.2, alpha = alpha_bars, fill = "grey80"
        ) +
        # best (minimum) model
        stat_summary(
          fun = "min", fun.min = "min", fun.max = "min",
          geom = "point", shape = 16, size = 1.2, color = col_ind_best
        ) +
        stat_summary(
          fun = "min", fun.min = "min", fun.max = "min",
          geom = "crossbar", linewidth = lwd_bar_overall, width = 0.83, color = col_ind_best
        ) +
        scale_y_continuous(breaks = y_breaks) +
        facet_nested(. ~ pred_type + month_d,
          scales = "free", space = "free",
          nest_line = element_line(), solo_line = TRUE
        ) +
        coord_cartesian(expand = FALSE, clip = "off") +
        xlab("Mean") + ylab(mm) +
        theme_minimal() +
        theme(
          # axis
          axis.text.x  = element_text(hjust = 0, vjust = 0.5, angle = 90),
          axis.ticks.x = element_blank(),
          # grid
          panel.grid.major.x = element_line(linewidth = 4, color = "white"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          # strips
          strip.background = element_blank(),
          strip.text       = element_text(size = 9, margin = margin(2, 0, 2, 0)),
          strip.placement  = "outside"
        )

      pdf(
        file = paste0(save.dir, plot.dir, "/mean/expl_mean_", tt, "_", mm, "_", ss_lab, ".pdf"),
        width = 3.8, height = 2.3
      )
      print(pp_means)
      dev.off()
    }
  }
}

## --------- MCS shares (25 percent, CRPS, all horizons) ---------------------------- #
plot.dir <- "mcs"
dir.create(paste0(save.dir, plot.dir), showWarnings = FALSE)

# compression order (display labels with the polynomial degree stripped)
midas_ord <- c("eq", "xalm", "alm", "ber", "leg", "fou", "br", "blu", "svd")

for (conflvl in c(25)) {
  for (spl in c("full", "full-outl")) {
    for (sl_metric in c("CRPS")) {
      lab <- "all"
      sl_horz <- c("-1/3", "0", "1/3", "2/3", "1", "4/3", "5/3")

      avg <- mcs_df %>%
        subset(mean != "MFVAR") %>%
        mutate(midas = gsub("5", "", gsub("3", "", midas))) %>%
        subset(metric %in% sl_metric) %>%
        subset(split %in% spl) %>%
        subset(month_d %in% sl_horz) %>%
        group_by(target, mean, midas) %>%
        summarize(
          conf25_perc = mean(conf25, na.rm = TRUE),
          conf10_perc = mean(conf10, na.rm = TRUE)
        )
      avg$midas <- factor(avg$midas, levels = midas_ord)

      # panel average: unweighted mean of the bars
      avg_summary <- avg %>%
        dplyr::group_by(target, mean) %>%
        dplyr::summarise(panel_avg = mean(conf25_perc), .groups = "drop")
      avg_plot <- avg %>%
        ggplot(aes(y = conf25_perc))

      pp_shares <- avg_plot +
        geom_hline(
          yintercept = seq(0, 1, by = 0.25),
          color = "grey85",
          linewidth = 0.5,
          linetype = "solid"
        ) +
        geom_bar(aes(x = midas), stat = "identity", position = position_dodge()) +
        geom_hline(yintercept = c(0, 1)) +
        geom_hline(
          data = avg_summary,
          aes(yintercept = panel_avg),
          linewidth = 1,
          color = "black"
        ) +
        facet_grid(target ~ mean, scales = "free", space = "free", switch = "y") +
        coord_cartesian(expand = FALSE, clip = "off", ylim = c(0, 1)) +
        ylab("Specifications in SMS (%)") + xlab("Compression") +
        theme_minimal() +
        theme(
          legend.position = "bottom",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.background = element_blank(),
          strip.placement = "outside",
          panel.grid = element_blank(),
          panel.spacing.x = unit(0.25, "cm"),
          panel.spacing.y = unit(0.5, "cm"),
          panel.grid.major.x = element_line(linewidth = 3.4, color = "grey95"),
          strip.text = element_text(size = 9, margin = margin(2, 0, 2, 0)),
          plot.title = element_text(size = 9),
          legend.position.inside = c(0, 0),
          legend.justification = c(0, 0),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.margin = margin(0, 0, 0, 0),
          legend.spacing.x = unit(0, "pt"),
          legend.spacing.y = unit(0, "pt")
        )

      pdf(file = paste0(save.dir, plot.dir, "/mcs-shares_", conflvl, "_", spl, "_", sl_metric, "_", gsub("/", "", lab), ".pdf"), width = 4, height = 2.5)
      print(pp_shares)
      dev.off()
    }
  }
}
