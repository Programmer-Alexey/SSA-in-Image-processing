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
  stop("РќРµ СѓРґР°Р»РѕСЃСЊ РЅР°Р№С‚Рё РєРѕСЂРµРЅСЊ РїСЂРѕРµРєС‚Р°.")
}

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[1])
} else {
  tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
}
script_dir <- if (!is.null(script_path)) {
  dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
} else {
  getwd()
}

project_root <- find_project_root(script_dir)
setwd(project_root)

source(file.path(project_root, "ssa-based methods", "cssa-transform.r"))

clip01 <- function(m) {
  m <- as.matrix(m)
  m[!is.finite(m)] <- 0
  m[m < 0.001] <- 0
  m[m > 0.99] <- 1
  m
}

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
  est_pos <- pmax(est, 0)
  pos_sum <- sum(est_pos)
  near_idx <- max(1L, signal_index - 1L):min(length(est), signal_index + 1L)

  data.frame(
    pred_index = which.max(est),
    argmax_error = abs(which.max(est) - signal_index),
    mse = mean((est - truth)^2),
    near_mass_share = if (pos_sum > 0) sum(est_pos[near_idx]) / pos_sum else NA_real_,
    active_share_0001 = mean(est > 0.001),
    active_share_01 = mean(est > 0.1),
    max_value = max(est),
    true_index_value = est[signal_index],
    stringsAsFactors = FALSE
  )
}

run_unit_frequency_experiment <- function(n = 100L, n_rep = 1L, sigma = 0.2, seed = 1111L) {
  rows <- vector("list", n * n_rep)
  out_i <- 1L

  set.seed(seed)
  for (rep_i in seq_len(n_rep)) {
    for (signal_index in seq_len(n)) {
      clean <- make_unit_row(n, signal_index)
      noisy <- clean + rnorm(n, sd = sigma)
      est <- cssa_rank1_row(noisy)

      omega <- -2 * pi * (signal_index - 1L) / n
      omega_wrapped <- atan2(sin(omega), cos(omega))
      metrics <- row_metrics(est, signal_index)

      rows[[out_i]] <- data.frame(
        rep = rep_i,
        signal_index = signal_index,
        omega = omega,
        omega_wrapped = omega_wrapped,
        metrics,
        stringsAsFactors = FALSE
      )
      out_i <- out_i + 1L
    }
  }

  do.call(rbind, rows)
}

summarise_by_frequency <- function(samples) {
  aggregate(
    cbind(
      pred_index,
      argmax_error,
      mse,
      near_mass_share,
      active_share_0001,
      active_share_01,
      max_value,
      true_index_value
    ) ~ signal_index + omega + omega_wrapped,
    data = samples,
    FUN = mean
  )
}

draw_frequency_error_plot <- function(summary_df, file, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::pdf(file, width = 6.2, height = 3.8)
  } else {
    grDevices::png(file, width = 1300, height = 800, res = 180)
  }
  on.exit(grDevices::dev.off(), add = TRUE)

  plot(
    summary_df$omega,
    summary_df$argmax_error,
    type = "b",
    pch = 19,
    col = "#1f4e79",
    xlab = expression(omega[j] == -2*pi*(j-1)/N),
    ylab = expression(abs(hat(j) - j)),
    main = "РћС€РёР±РєР° Р»РѕРєР°Р»РёР·Р°С†РёРё РґР»СЏ РµРґРёРЅРёС‡РЅРѕРіРѕ РѕСЂС‚Р°"
  )
  grid(col = "gray85")
  abline(stats::lm(argmax_error ~ omega, data = summary_df), col = "red", lwd = 2)
}

slides_dir <- file.path(project_root, "VKR", "VKR_SLIDES")
tables_dir <- file.path(slides_dir, "tables")
images_dir <- file.path(slides_dir, "images")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(images_dir, recursive = TRUE, showWarnings = FALSE)

n <- as.integer(Sys.getenv("UNIT_FREQUENCY_N", "100"))
n_rep <- as.integer(Sys.getenv("UNIT_FREQUENCY_N_REP", "1"))
sigma <- as.numeric(Sys.getenv("UNIT_FREQUENCY_SIGMA", "0.2"))
seed <- as.integer(Sys.getenv("UNIT_FREQUENCY_SEED", "1111"))

samples <- run_unit_frequency_experiment(
  n = n,
  n_rep = n_rep,
  sigma = sigma,
  seed = seed
)
summary_df <- summarise_by_frequency(samples)
summary_df <- summary_df[order(summary_df$signal_index), ]

write.csv(
  samples,
  file.path(tables_dir, "unit_ort_frequency_error_samples.csv"),
  row.names = FALSE
)
write.csv(
  summary_df,
  file.path(tables_dir, "unit_ort_frequency_error_summary.csv"),
  row.names = FALSE
)

draw_frequency_error_plot(
  summary_df,
  file.path(images_dir, "unit_ort_frequency_localization.pdf"),
  device = "pdf"
)
draw_frequency_error_plot(
  summary_df,
  file.path(images_dir, "unit_ort_frequency_localization.png"),
  device = "png"
)

cat(
  "Р“РѕС‚РѕРІРѕ: РїРѕСЃС‚СЂРѕРµРЅ РіСЂР°С„РёРє РѕС€РёР±РєРё РѕС‚ omega_j РґР»СЏ j=1..", n,
  ", n_rep=", n_rep,
  ", sigma=", sigma,
  "\n",
  sep = ""
)

