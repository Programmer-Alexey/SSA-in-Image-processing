find_project_root <- function(start_dir) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in 1:20) {
    if (file.exists(file.path(cur, "SSA-in-Image-processing.Rproj"))) {
      return(cur)
    }
    next_dir <- dirname(cur)
    if (identical(next_dir, cur)) break
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

project_lib <- normalizePath(".r-local-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

required_pkgs <- c("Rcpp", "Rssa", "dplyr", "imager", "EBImage", "reticulate")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("РћС‚СЃСѓС‚СЃС‚РІСѓСЋС‚ РїР°РєРµС‚С‹: ", paste(missing_pkgs, collapse = ", "))
}

source(file.path(project_root, "ssa-based methods", "cssa-transform.r"))
source(file.path(project_root, "hough transform", "hough_transform.r"))

library(dplyr)

slides_dir <- file.path(project_root, "VKR", "VKR_SLIDES", "SSA Slides")
tables_dir <- file.path(slides_dir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

clip01 <- function(m) {
  m <- as.matrix(m)
  m[!is.finite(m)] <- 0
  m[m < 0.001] <- 0
  m[m > 0.99] <- 1
  m
}

threshold01 <- function(m, threshold = 0.1) {
  out <- clip01(m)
  out[out < threshold] <- 0
  out
}

make_line_image <- function(n_row, n_col, a, b, intensity = 1) {
  matrix(0, nrow = n_row, ncol = n_col) |>
    add.line(a = a, b = b, intensity = intensity)
}

cssa_denoise <- function(m, num_of_lines = 1L, method = "row.row") {
  cleaned <- switch(
    method,
    "row.row" = m |> dft() |> cssa.row(num.line = num_of_lines) |> idft.row() |> Re(),
    stop("РќРµРёР·РІРµСЃС‚РЅС‹Р№ CSSA-РјРµС‚РѕРґ: ", method)
  )
  clip01(cleaned)
}

bm3d_available <- reticulate::py_module_available("bm3d")
bm3d_lib <- NULL

bm3d_denoise <- function(m, sigma_noise) {
  if (!bm3d_available) {
    stop("Python-РјРѕРґСѓР»СЊ 'bm3d' РЅРµРґРѕСЃС‚СѓРїРµРЅ.")
  }
  if (is.null(bm3d_lib)) {
    bm3d_lib <<- reticulate::import("bm3d")
  }
  x <- as.matrix(m)
  x[!is.finite(x)] <- 0
  clip01(bm3d_lib$bm3d(x, sigma_psd = as.numeric(sigma_noise)))
}

median_denoise <- function(m, n = 3L) {
  clip01(as.matrix(imager::medianblur(imager::as.cimg(clip01(m)), n = as.integer(n))))
}

wiener_denoise <- function(m, ksize = 5L) {
  x <- as.matrix(m)
  ksize <- as.integer(ksize)
  if (ksize %% 2L == 0L) {
    ksize <- ksize + 1L
  }

  kernel <- matrix(1 / (ksize * ksize), nrow = ksize, ncol = ksize)
  local_mean <- EBImage::filter2(x, kernel, boundary = "replicate")
  local_var <- EBImage::filter2(x^2, kernel, boundary = "replicate") - local_mean^2
  local_var[local_var < 0] <- 0
  noise_var <- mean(local_var, na.rm = TRUE)

  gain <- pmax(local_var - noise_var, 0) / pmax(local_var, 1e-8)
  clip01(local_mean + gain * (x - local_mean))
}

quantile_denoise <- function(m, q = 0.9) {
  x <- clip01(m)
  thr <- as.numeric(stats::quantile(x, probs = q, na.rm = TRUE))
  clip01(ifelse(x >= thr, x, 0))
}

convert_ab_to_rho_theta <- function(a, b) {
  theta <- atan2(1, -a)
  rho <- b / sqrt(a^2 + 1)
  if (theta < 0) {
    theta <- theta + pi
    rho <- -rho
  }
  c(rho = rho, theta = theta)
}

find_k_max <- function(acc, k, qrho, qtheta) {
  acc_copy <- acc
  maxima <- vector("list", k)

  for (i in seq_len(k)) {
    idx <- which(acc_copy == max(acc_copy), arr.ind = TRUE)[1, ]
    maxima[[i]] <- c(rho = qrho[idx[1]], theta = qtheta[idx[2]])
    acc_copy[idx[1], idx[2]] <- 0
  }

  out <- do.call(rbind, maxima)
  colnames(out) <- c("rho", "theta")
  out
}

line_from_rho_theta <- function(rho, theta, n_row, n_col, eps = 1e-8) {
  out <- matrix(0, nrow = n_row, ncol = n_col)

  if (abs(sin(theta)) > eps) {
    return(add.line(
      out,
      a = -cos(theta) / sin(theta),
      b = rho / sin(theta),
      intensity = 1
    ))
  }

  x0 <- round(rho / cos(theta))
  if (is.finite(x0) && x0 >= 1 && x0 <= n_col) {
    out[, x0] <- 1
  }
  out
}

distance_metrics_to_true_line <- function(line_matrix, a_true, b_true) {
  pts <- which(as.matrix(line_matrix) > 0, arr.ind = TRUE)
  if (nrow(pts) == 0) {
    return(c(mean_distance = NA_real_, sum_distance = 0, active_pixels = 0))
  }

  y <- pts[, "row"]
  x <- pts[, "col"]
  distances <- abs(a_true * x - y + b_true) / sqrt(a_true^2 + 1)

  c(
    mean_distance = mean(distances),
    sum_distance = sum(distances),
    active_pixels = nrow(pts)
  )
}

run_one_ht <- function(input_matrix, processed_matrix, cfg, rho_step, theta_step) {
  ht_result <- make_accumulator(
    input_matrix,
    detector = function(m) processed_matrix,
    rho_step = rho_step,
    theta_step = theta_step
  )
  pred_line <- find_k_max(ht_result$accumulator, k = 1L, qrho = ht_result$rho, qtheta = ht_result$theta)
  reconstructed_line <- line_from_rho_theta(
    rho = pred_line[1, 1],
    theta = pred_line[1, 2],
    n_row = cfg$n_row,
    n_col = cfg$n_col
  )
  dist <- distance_metrics_to_true_line(reconstructed_line, a_true = cfg$a, b_true = cfg$b)
  true_line <- convert_ab_to_rho_theta(cfg$a, cfg$b)

  data.frame(
    pred_rho = pred_line[1, 1],
    pred_theta = pred_line[1, 2],
    rho_mse = (pred_line[1, 1] - true_line["rho"])^2,
    theta_mse = (pred_line[1, 2] - true_line["theta"])^2,
    reconstructed_mean_distance = as.numeric(dist["mean_distance"]),
    processed_active_pixels = sum(processed_matrix > 0),
    stringsAsFactors = FALSE
  )
}

escape_tex <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([_%&#$])", "\\\\\\1", x, perl = TRUE)
  x
}

write_df_tex <- function(df, file, digits = 4L) {
  df_out <- df
  for (j in seq_along(df_out)) {
    if (is.numeric(df_out[[j]])) {
      df_out[[j]] <- format(round(df_out[[j]], digits), nsmall = digits, trim = TRUE, scientific = FALSE)
    } else {
      df_out[[j]] <- escape_tex(df_out[[j]])
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

wide_metric_table <- function(summary_df, metric, config_levels, method_levels) {
  wide <- xtabs(summary_df[[metric]] ~ summary_df$method_label + summary_df$config_label)
  wide <- wide[method_levels, config_levels, drop = FALSE]
  data.frame(
    "РњРµС‚РѕРґ" = rownames(wide),
    as.data.frame.matrix(wide),
    check.names = FALSE,
    row.names = NULL
  )
}

sigma_tag <- function(sigma_value) {
  gsub("\\.", "", sprintf("%.2f", sigma_value))
}

grid_metric_table <- function(summary_df, metric, method_levels) {
  grid_levels <- c("$2; 0.02$", "$1; 0.01$", "$0.5; 0.005$")
  grid_wide <- xtabs(summary_df[[metric]] ~ summary_df$method_label + summary_df$grid_label)
  grid_wide <- grid_wide[method_levels, grid_levels, drop = FALSE]
  data.frame(
    "РњРµС‚РѕРґ" = rownames(grid_wide),
    as.data.frame.matrix(grid_wide),
    check.names = FALSE,
    row.names = NULL
  )
}

configs <- data.frame(
  config_id = c("diag_pos_full", "diag_neg_full", "steep_full", "steep_shifted"),
  config_label = c("РџРѕР»РѕР¶РёС‚РµР»СЊРЅР°СЏ РґРёР°РіРѕРЅР°Р»СЊ", "РћС‚СЂРёС†Р°С‚РµР»СЊРЅР°СЏ РґРёР°РіРѕРЅР°Р»СЊ", "РљСЂСѓС‚Р°СЏ РїСЂСЏРјР°СЏ", "РљСЂСѓС‚Р°СЏ СЃРјРµС‰РµРЅРЅР°СЏ"),
  n_col = rep(100L, 4),
  n_row = rep(100L, 4),
  a = c(1, -1, 2, 2),
  b = c(0, 101, -1, -40),
  stringsAsFactors = FALSE
)

method_labels <- c(
  cssa_row_row = "CSSA ROW-ROW",
  bm3d = "BM3D",
  median = "Median",
  wiener = "Wiener",
  quantile_09 = "Quantile 0.9"
)

detector_q <- 0.9
detectors <- list(
  cssa_row_row = function(m, sigma_noise) cssa_denoise(m, num_of_lines = 1L, method = "row.row"),
  median = function(m, sigma_noise) median_denoise(m, n = 3L),
  wiener = function(m, sigma_noise) wiener_denoise(m, ksize = 5L),
  quantile_09 = function(m, sigma_noise) quantile_denoise(m, q = detector_q)
)

if (bm3d_available) {
  detectors <- append(
    detectors,
    list(bm3d = function(m, sigma_noise) bm3d_denoise(m, sigma_noise = sigma_noise)),
    after = 1L
  )
} else {
  warning("BM3D РїСЂРѕРїСѓС‰РµРЅ: Python-РјРѕРґСѓР»СЊ 'bm3d' РЅРµ РЅР°Р№РґРµРЅ.")
}

method_order <- names(method_labels)
method_order_current <- method_order[method_order %in% names(detectors)]
method_levels_current <- unname(method_labels[method_order_current])

sigma_noise <- as.numeric(Sys.getenv("THRESHOLD_GRID_SIGMA", "0.2"))
threshold_n_rep <- as.integer(Sys.getenv("THRESHOLD_01_N_REP", "100"))
grid_n_rep <- as.integer(Sys.getenv("GRID_STEP_N_REP", "30"))
set.seed(as.integer(Sys.getenv("THRESHOLD_GRID_SEED", "1111")))
tag <- sigma_tag(sigma_noise)

run_threshold_01 <- function(n_rep) {
  rows <- vector("list", length = n_rep * nrow(configs) * length(detectors))
  out_i <- 1L

  for (cfg_i in seq_len(nrow(configs))) {
    cfg <- configs[cfg_i, ]
    base_matrix <- make_line_image(cfg$n_row, cfg$n_col, cfg$a, cfg$b)

    for (rep_i in seq_len(n_rep)) {
      noisy_matrix <- add.noise(base_matrix, sigma = sigma_noise)

      for (method_name in names(detectors)) {
        processed_matrix <- detectors[[method_name]](noisy_matrix, sigma_noise = sigma_noise)
        processed_matrix <- threshold01(processed_matrix, threshold = 0.1)
        metrics <- run_one_ht(noisy_matrix, processed_matrix, cfg, rho_step = 1, theta_step = 0.01)

        rows[[out_i]] <- cbind(
          data.frame(
            config_id = cfg$config_id,
            config_label = cfg$config_label,
            sigma = sigma_noise,
            rep = rep_i,
            method = method_name,
            method_label = method_labels[[method_name]],
            stringsAsFactors = FALSE
          ),
          metrics
        )
        out_i <- out_i + 1L
      }
    }
  }

  dplyr::bind_rows(rows)
}

run_grid_step <- function(n_rep) {
  grid_specs <- data.frame(
    grid_label = c("$2; 0.02$", "$1; 0.01$", "$0.5; 0.005$"),
    rho_step = c(2, 1, 0.5),
    theta_step = c(0.02, 0.01, 0.005),
    stringsAsFactors = FALSE
  )
  rows <- vector("list", length = n_rep * nrow(configs) * length(detectors) * nrow(grid_specs))
  out_i <- 1L

  for (cfg_i in seq_len(nrow(configs))) {
    cfg <- configs[cfg_i, ]
    base_matrix <- make_line_image(cfg$n_row, cfg$n_col, cfg$a, cfg$b)

    for (rep_i in seq_len(n_rep)) {
      noisy_matrix <- add.noise(base_matrix, sigma = sigma_noise)

      for (method_name in names(detectors)) {
        processed_matrix <- detectors[[method_name]](noisy_matrix, sigma_noise = sigma_noise)

        for (grid_i in seq_len(nrow(grid_specs))) {
          grid <- grid_specs[grid_i, ]
          metrics <- run_one_ht(
            noisy_matrix,
            processed_matrix,
            cfg,
            rho_step = grid$rho_step,
            theta_step = grid$theta_step
          )

          rows[[out_i]] <- cbind(
            data.frame(
              config_id = cfg$config_id,
              config_label = cfg$config_label,
              sigma = sigma_noise,
              rep = rep_i,
              method = method_name,
              method_label = method_labels[[method_name]],
              grid_label = grid$grid_label,
              rho_step = grid$rho_step,
              theta_step = grid$theta_step,
              stringsAsFactors = FALSE
            ),
            metrics
          )
          out_i <- out_i + 1L
        }
      }
    }
  }

  dplyr::bind_rows(rows)
}

threshold_samples_path <- file.path(tables_dir, paste0("threshold_01_samples_sigma_", tag, ".csv"))
threshold_summary_path <- file.path(tables_dir, paste0("threshold_01_summary_sigma_", tag, ".csv"))

if (threshold_n_rep > 0L) {
  threshold_samples <- run_threshold_01(threshold_n_rep)
  write.csv(threshold_samples, threshold_samples_path, row.names = FALSE)
} else if (file.exists(threshold_samples_path)) {
  threshold_samples <- read.csv(threshold_samples_path, stringsAsFactors = FALSE)
} else {
  threshold_samples <- NULL
}

if (!is.null(threshold_samples)) {
  threshold_summary <- threshold_samples |>
    dplyr::group_by(config_id, config_label, method, method_label, sigma) |>
    dplyr::summarise(
      mean_distance = mean(reconstructed_mean_distance, na.rm = TRUE),
      mean_rho_mse = mean(rho_mse, na.rm = TRUE),
      mean_theta_mse = mean(theta_mse, na.rm = TRUE),
      mean_active_pixels = mean(processed_active_pixels, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      config_label = factor(config_label, levels = configs$config_label),
      method_label = factor(method_label, levels = method_levels_current)
    ) |>
    dplyr::arrange(config_label, method_label)

  write.csv(threshold_summary, threshold_summary_path, row.names = FALSE)

  for (metric_name in c("mean_distance", "mean_rho_mse", "mean_theta_mse")) {
    threshold_df <- wide_metric_table(
      threshold_summary,
      metric = metric_name,
      config_levels = configs$config_label,
      method_levels = method_levels_current
    )
    out_name <- sub("^mean_", "threshold_01_mean_", metric_name)
    write_df_tex(
      threshold_df,
      file.path(tables_dir, paste0(out_name, "_sigma_", tag, ".tex")),
      digits = if (metric_name == "mean_theta_mse") 6L else 4L
    )
  }
}

grid_samples_path <- file.path(tables_dir, paste0("grid_step_samples_sigma_", tag, ".csv"))
grid_summary_path <- file.path(tables_dir, paste0("grid_step_summary_sigma_", tag, ".csv"))

if (grid_n_rep > 0L) {
  grid_samples <- run_grid_step(grid_n_rep)
  write.csv(grid_samples, grid_samples_path, row.names = FALSE)
} else if (file.exists(grid_samples_path)) {
  grid_samples <- read.csv(grid_samples_path, stringsAsFactors = FALSE)
} else {
  grid_samples <- NULL
}

if (!is.null(grid_samples)) {
  grid_summary <- grid_samples |>
    dplyr::group_by(method, method_label, grid_label, rho_step, theta_step, sigma) |>
    dplyr::summarise(
      mean_distance = mean(reconstructed_mean_distance, na.rm = TRUE),
      mean_rho_mse = mean(rho_mse, na.rm = TRUE),
      mean_theta_mse = mean(theta_mse, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(method_label = factor(method_label, levels = method_levels_current)) |>
    dplyr::arrange(method_label, rho_step)

  write.csv(grid_summary, grid_summary_path, row.names = FALSE)

  for (metric_name in c("mean_distance", "mean_rho_mse", "mean_theta_mse")) {
    grid_df <- grid_metric_table(
      grid_summary,
      metric = metric_name,
      method_levels = method_levels_current
    )
    out_name <- sub("^mean_", "grid_step_mean_", metric_name)
    write_df_tex(
      grid_df,
      file.path(tables_dir, paste0(out_name, "_sigma_", tag, ".tex")),
      digits = if (metric_name == "mean_theta_mse") 6L else 4L
    )
  }
}

cat("Р“РѕС‚РѕРІРѕ: С‚Р°Р±Р»РёС†С‹ threshold_01 Рё grid_step РґР»СЏ sigma=", sigma_noise, "\n", sep = "")

