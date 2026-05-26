find_project_root <- function(start_dir) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in 1:20) {
    if (file.exists(file.path(cur, "SSA-in-Image-processing.Rproj"))) {
      return(cur)
    }
    next_dir <- dirname(cur)
    if (identical(next_dir, cur)) {
      break
    }
    cur <- next_dir
  }
  stop("Could not find project root.")
}

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile, winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/")
}

script_dir <- get_script_dir()
repo_root <- find_project_root(script_dir)
text_dir <- file.path(repo_root, "VKR", "VKR_TEXT")
tables_dir <- file.path(text_dir, "chapters", "chapter3_experiments", "tables")
data_dir <- file.path(tables_dir, "data")
images_dir <- file.path(text_dir, "assets", "images", "chapter3")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(images_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(repo_root, "ssa-based methods", "cssa-transform.r"))

make_unit_row <- function(n, signal_index, intensity = 1) {
  x <- numeric(n)
  x[signal_index] <- intensity
  x
}

cssa_rank1_row <- function(row_vec) {
  row_vec |>
    matrix(nrow = 1L) |>
    dft() |>
    cssa.row(num.line = 1L) |>
    idft.row() |>
    Re() |>
    as.numeric()
}

row_metrics <- function(est, signal_index) {
  truth <- make_unit_row(length(est), signal_index)
  pred_index <- which.max(est)
  est_pos <- pmax(est, 0)
  pos_sum <- sum(est_pos)
  q <- if (pos_sum > 0) est_pos / pos_sum else rep(NA_real_, length(est))
  peak_distance <- abs(seq_along(est) - pred_index)
  near_idx <- max(1L, signal_index - 1L):min(length(est), signal_index + 1L)

  data.frame(
    pred_index = pred_index,
    argmax_error = abs(pred_index - signal_index),
    mse = mean((est - truth)^2),
    near_mass_share = if (pos_sum > 0) sum(est_pos[near_idx]) / pos_sum else NA_real_,
    peak_spread_l1 = if (pos_sum > 0) sum(peak_distance * q) else NA_real_,
    peak_spread_l2_sq = if (pos_sum > 0) sum(peak_distance^2 * q) else NA_real_,
    active_share_0001 = mean(est > 0.001),
    active_share_01 = mean(est > 0.1),
    max_value = max(est),
    true_index_value = est[signal_index],
    stringsAsFactors = FALSE
  )
}

run_frequency_rep <- function(rep_i, rep_seed, n, sigma, signal_indices) {
  set.seed(rep_seed)
  noise_vec <- rnorm(n, sd = sigma)

  rows <- lapply(signal_indices, function(signal_index) {
    clean <- make_unit_row(n, signal_index)
    noisy <- clean + noise_vec
    est <- cssa_rank1_row(noisy)

    omega_full <- round((signal_index - 1L) / n, 12)
    omega_wrapped <- round(atan2(sin(2 * pi * omega_full), cos(2 * pi * omega_full)) / (2 * pi), 12)
    omega <- round(abs(omega_wrapped), 12)

    data.frame(
      rep = rep_i,
      signal_index = signal_index,
      omega_full = omega_full,
      omega = omega,
      omega_wrapped = omega_wrapped,
      row_metrics(est, signal_index),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

split_indices <- function(x, n_groups) {
  split(x, rep(seq_len(n_groups), length.out = length(x)))
}

run_frequency_experiment_parallel <- function(n, n_rep, sigma, seed, n_workers, signal_indices) {
  set.seed(seed)
  rep_seeds <- sample.int(.Machine$integer.max, n_rep)
  rep_ids <- seq_len(n_rep)

  n_workers <- max(1L, min(as.integer(n_workers), n_rep))
  if (n_workers == 1L) {
    return(do.call(rbind, lapply(rep_ids, function(rep_i) {
      run_frequency_rep(rep_i, rep_seeds[rep_i], n, sigma, signal_indices)
    })))
  }

  chunks <- split_indices(rep_ids, n_workers)
  cl <- parallel::makeCluster(n_workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterExport(cl, varlist = "repo_root_worker", envir = list2env(list(repo_root_worker = repo_root)))
  parallel::clusterEvalQ(cl, {
    source(file.path(repo_root_worker, "ssa-based methods", "cssa-transform.r"))
    NULL
  })
  parallel::clusterExport(
    cl,
    varlist = c(
      "repo_root",
      "make_unit_row",
      "cssa_rank1_row",
      "row_metrics",
      "run_frequency_rep",
      "rep_seeds",
      "n",
      "sigma",
      "signal_indices"
    ),
    envir = environment()
  )

  out <- parallel::parLapply(cl, chunks, function(rep_chunk) {
    do.call(rbind, lapply(rep_chunk, function(rep_i) {
      run_frequency_rep(rep_i, rep_seeds[rep_i], n, sigma, signal_indices)
    }))
  })

  do.call(rbind, out)
}

mean_or_na <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

add_peak_spread_l2 <- function(df) {
  if ("peak_spread_l2_sq" %in% names(df)) {
    df$peak_spread_l2 <- sqrt(df$peak_spread_l2_sq)
  }
  df
}

summarise_by_index <- function(samples) {
  out <- aggregate(
    cbind(
      pred_index,
      argmax_error,
      mse,
      near_mass_share,
      peak_spread_l1,
      peak_spread_l2_sq,
      active_share_0001,
      active_share_01,
      max_value,
      true_index_value
    ) ~ signal_index + omega_full + omega + omega_wrapped,
    data = samples,
    FUN = mean_or_na
  )
  out <- add_peak_spread_l2(out)
  out[order(out$signal_index), ]
}

summarise_by_frequency <- function(samples) {
  summary <- aggregate(
    cbind(
      argmax_error,
      mse,
      near_mass_share,
      peak_spread_l1,
      peak_spread_l2_sq,
      active_share_0001,
      active_share_01,
      max_value,
      true_index_value
    ) ~ omega,
    data = samples,
    FUN = mean_or_na
  )
  counts <- aggregate(rep ~ omega, data = samples, FUN = length)
  names(counts)[names(counts) == "rep"] <- "n_samples"
  out <- merge(summary, counts, by = "omega", all.x = TRUE)
  out <- add_peak_spread_l2(out)
  out[order(out$omega), ]
}

write_df_tex <- function(df, file, digits = 3L) {
  df_out <- df
  for (j in seq_along(df_out)) {
    if (is.numeric(df_out[[j]])) {
      df_out[[j]] <- format(round(df_out[[j]], digits), nsmall = digits, trim = TRUE)
    } else {
      df_out[[j]] <- as.character(df_out[[j]])
      is_math <- grepl("^\\$.*\\$$", df_out[[j]])
      df_out[[j]][!is_math] <- gsub("_", "\\\\_", df_out[[j]][!is_math])
    }
  }

  align <- paste(rep("l", ncol(df_out)), collapse = "|")
  lines <- c(
    paste0("\\begin{tabular}{|", align, "|}"),
    "\\hline",
    paste0(paste(names(df_out), collapse = " & "), " \\\\"),
    "\\hline"
  )
  for (i in seq_len(nrow(df_out))) {
    lines <- c(lines, paste0(paste(df_out[i, ], collapse = " & "), " \\\\"), "\\hline")
  }
  lines <- c(lines, "\\end{tabular}")
  writeLines(lines, file, useBytes = TRUE)
}

draw_frequency_plot <- function(summary_df, file, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::pdf(file, width = 6.4, height = 4.0)
  } else {
    grDevices::png(file, width = 1300, height = 820, res = 180)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(mar = c(5.1, 6.8, 4.1, 2.1))
  on.exit(graphics::par(old_par), add = TRUE)

  plot(
    summary_df$omega,
    summary_df$argmax_error,
    type = "b",
    pch = 19,
    col = "#1f4e79",
    xlab = expression(omega[phantom(.) * j]),
    ylab = expression(E~group("|", hat(j) - j, "|")),
    main = "Mean localization error"
  )
  grid(col = "gray85")
  graphics::lines(
    stats::lowess(summary_df$omega, summary_df$argmax_error, f = 0.35),
    col = "#c00000",
    lwd = 2
  )
}

draw_full_plot <- function(summary_df, file, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::pdf(file, width = 6.4, height = 4.0)
  } else {
    grDevices::png(file, width = 1300, height = 820, res = 180)
  }
  on.exit(grDevices::dev.off(), add = TRUE)

  plot(
    summary_df$omega,
    summary_df$argmax_error,
    type = "b",
    pch = 19,
    col = "#1f4e79",
    xlab = expression(omega[j] == (j - 1) / N),
    ylab = "Mean |j_hat - j|",
    main = "Mean localization error"
  )
  grid(col = "gray85")
  graphics::lines(
    stats::lowess(summary_df$omega, summary_df$argmax_error, f = 0.35),
    col = "#c00000",
    lwd = 2
  )
}

draw_peak_spread_plot <- function(summary_df, file, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::pdf(file, width = 6.4, height = 4.0)
  } else {
    grDevices::png(file, width = 1300, height = 820, res = 180)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(mar = c(5.1, 6.8, 4.1, 2.1))
  on.exit(graphics::par(old_par), add = TRUE)

  y_values <- c(summary_df$peak_spread_l1, summary_df$peak_spread_l2)
  plot(
    summary_df$omega,
    summary_df$peak_spread_l2,
    type = "b",
    pch = 19,
    lwd = 2,
    col = "#1f4e79",
    ylim = range(y_values, na.rm = TRUE),
    xlab = expression(omega[phantom(.) * j]),
    ylab = expression(S[p](omega)),
    main = "Mean spread around estimated peak"
  )
  grid(col = "gray85")
  graphics::lines(
    summary_df$omega,
    summary_df$peak_spread_l2,
    type = "b",
    pch = 19,
    lwd = 2,
    col = "#1f4e79"
  )
  graphics::lines(
    summary_df$omega,
    summary_df$peak_spread_l1,
    type = "b",
    pch = 17,
    lwd = 2,
    col = "#c00000"
  )
  legend(
    "topright",
    legend = c("S2", "S1"),
    col = c("#1f4e79", "#c00000"),
    lwd = 2,
    pch = c(19, 17),
    bty = "n"
  )
}

parse_indices <- function(value, n) {
  if (!nzchar(value)) {
    return(seq_len(n))
  }
  idx <- as.integer(strsplit(value, "[,;[:space:]]+", perl = TRUE)[[1]])
  idx <- idx[is.finite(idx)]
  idx <- unique(idx[idx >= 1L & idx <= n])
  if (length(idx) == 0L) {
    stop("FREQ_INDICES did not contain valid indices.")
  }
  sort(idx)
}

n <- as.integer(Sys.getenv("CH3_FREQ_N", "100"))
n_rep <- as.integer(Sys.getenv("CH3_FREQ_N_REP", "500"))
sigma <- as.numeric(Sys.getenv("CH3_FREQ_SIGMA", "0.2"))
seed <- as.integer(Sys.getenv("CH3_FREQ_SEED", "1111"))
n_workers <- as.integer(Sys.getenv("CH3_FREQ_WORKERS", as.character(min(8L, parallel::detectCores()))))
signal_indices <- parse_indices(Sys.getenv("CH3_FREQ_INDICES", ""), n)

cat(
  "Running chapter 3 frequency experiment: N=", n,
  ", n_rep=", n_rep,
  ", sigma=", sigma,
  ", workers=", n_workers,
  ", indices=", length(signal_indices),
  "\n",
  sep = ""
)

samples <- run_frequency_experiment_parallel(
  n = n,
  n_rep = n_rep,
  sigma = sigma,
  seed = seed,
  n_workers = n_workers,
  signal_indices = signal_indices
)

summary_by_index <- summarise_by_index(samples)
summary_frequency <- summarise_by_frequency(samples)

write.csv(samples, file.path(data_dir, "frequency_error_samples.csv"), row.names = FALSE)
write.csv(summary_by_index, file.path(data_dir, "frequency_error_summary_by_index.csv"), row.names = FALSE)
write.csv(summary_frequency, file.path(data_dir, "frequency_error_summary.csv"), row.names = FALSE)

params <- data.frame(
  "Параметр" = c(
    "$N$",
    "$\\sigma$",
    "Число повторов",
    "Число ядер",
    "Частоты",
    "Объем выборки для внутренней точки"
  ),
  "Значение" = c(
    as.character(n),
    format(sigma, trim = TRUE),
    as.character(n_rep),
    as.character(min(n_workers, n_rep)),
    "$0\\le\\omega\\le1/2$",
    paste0("$2\\cdot", n_rep, "$")
  ),
  check.names = FALSE
)
write_df_tex(params, file.path(tables_dir, "frequency_experiment_params.tex"), digits = 3L)

draw_frequency_plot(
  summary_frequency,
  file.path(images_dir, "frequency_localization.pdf"),
  device = "pdf"
)
draw_frequency_plot(
  summary_frequency,
  file.path(images_dir, "frequency_localization.png"),
  device = "png"
)
draw_full_plot(
  summary_by_index,
  file.path(images_dir, "frequency_localization_full.pdf"),
  device = "pdf"
)
draw_peak_spread_plot(
  summary_frequency,
  file.path(images_dir, "frequency_peak_spread.pdf"),
  device = "pdf"
)
draw_peak_spread_plot(
  summary_frequency,
  file.path(images_dir, "frequency_peak_spread.png"),
  device = "png"
)

cat("Done: frequency_localization.pdf, frequency_peak_spread.pdf and frequency_experiment_params.tex\n")
