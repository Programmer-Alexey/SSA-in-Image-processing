## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  cache = FALSE
)

set.seed(1111)


## ----libraries, include=FALSE-------------------------------------------------
project_lib <- normalizePath(".r-local-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

required_pkgs <- c(
  "Rcpp",
  "Rssa",
  "lattice",
  "gridExtra",
  "dplyr",
  "imager",
  "EBImage",
  "reticulate"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop("РћС‚СЃСѓС‚СЃС‚РІСѓСЋС‚ РїР°РєРµС‚С‹: ", paste(missing_pkgs, collapse = ", "))
}

source("ssa-based methods/cssa-transform.r")
source("hough transform/hough_transform.r")

library(dplyr)
library(reticulate)

bm3d_available <- reticulate::py_module_available("bm3d")
bm3d_lib <- NULL
line_experiments_dir <- file.path(project_root, "line_experiments")


## ----helpers------------------------------------------------------------------
clip01 <- function(m) {
  m <- as.matrix(m)
  m[!is.finite(m)] <- 0
  m[m < 0.001] <- 0
  m[m > 0.99] <- 1
  m
}

normalize_xtab <- function(tab) {
  t(apply(tab, 1, function(x) {
    row_min <- min(x, na.rm = TRUE)
    if (!is.finite(row_min)) {
      return(rep(NA_real_, length(x)))
    }
    if (row_min == 0) {
      ifelse(x == 0, 1, Inf)
    } else {
      x / row_min
    }
  }))
}

make_line_image <- function(n_row, n_col, a, b, intensity = 1) {
  matrix(0, nrow = n_row, ncol = n_col) |>
    add.line(a = a, b = b, intensity = intensity)
}

cssa_denoise <- function(m, num_of_lines = 1L, method = "row.row") {
  cleaned <- switch(
    method,
    "col.col" = m |> dft() |> cssa.col(num.line = num_of_lines) |> idft.col() |> Re(),
    "col.row" = m |> dft() |> cssa.col(num.line = num_of_lines) |> idft.row() |> Re(),
    "row.col" = m |> dft() |> cssa.row(num.line = num_of_lines) |> idft.col() |> Re(),
    "row.row" = m |> dft() |> cssa.row(num.line = num_of_lines) |> idft.row() |> Re(),
    stop("РќРµРёР·РІРµСЃС‚РЅС‹Р№ CSSA-РјРµС‚РѕРґ: ", method)
  )
  clip01(cleaned)
}

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
  x <- m
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
  out <- local_mean + gain * (x - local_mean)
  clip01(out)
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

find_k_max <- function(acc, k, qrho, qtheta, suppress = FALSE, window = 6L) {
  acc_copy <- acc
  maxima <- vector("list", k)

  for (i in seq_len(k)) {
    idx <- which(acc_copy == max(acc_copy), arr.ind = TRUE)[1, ]
    maxima[[i]] <- c(rho = qrho[idx[1]], theta = qtheta[idx[2]])

    if (suppress) {
      r_range <- max(1, idx[1] - window):min(nrow(acc_copy), idx[1] + window)
      c_range <- max(1, idx[2] - window):min(ncol(acc_copy), idx[2] + window)
      acc_copy[r_range, c_range] <- 0
    } else {
      acc_copy[idx[1], idx[2]] <- 0
    }
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
    return(c(
      mean_distance = NA_real_,
      sum_distance = 0,
      active_pixels = 0
    ))
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
      df_out[[j]] <- format(round(df_out[[j]], digits), nsmall = digits, trim = TRUE)
    } else {
      df_out[[j]] <- escape_tex(df_out[[j]])
    }
  }

  align <- paste(rep("l", ncol(df_out)), collapse = "|")
  header <- paste(names(df_out), collapse = " & ")
  body <- apply(df_out, 1, function(row) paste(row, collapse = " & "))

  lines <- c(
    paste0("\\begin{tabular}{|", align, "|}"),
    "\\hline",
    paste0(header, " \\\\"),
    "\\hline"
  )
  for (line in body) {
    lines <- c(lines, paste0(line, " \\\\"), "\\hline")
  }
  lines <- c(lines, "\\end{tabular}")
  writeLines(lines, con = file, useBytes = TRUE)
}

wide_metric_table <- function(summary_df, sigma_value, metric, config_levels, method_levels) {
  df <- summary_df |>
    dplyr::filter(sigma == sigma_value)

  wide <- xtabs(df[[metric]] ~ df$method_label + df$config_label)
  wide <- wide[method_levels, config_levels, drop = FALSE]

  data.frame(
    "РњРµС‚РѕРґ" = rownames(wide),
    as.data.frame.matrix(wide),
    check.names = FALSE,
    row.names = NULL
  )
}

wide_normalized_metric_table <- function(summary_df, sigma_value, metric, config_levels, method_levels) {
  df <- summary_df |>
    dplyr::filter(sigma == sigma_value)

  tab <- xtabs(df[[metric]] ~ df$config_label + df$method_label)
  tab <- tab[config_levels, method_levels, drop = FALSE]
  tab_norm <- normalize_xtab(tab)
  wide <- t(tab_norm)

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

method_label <- function(method) {
  out <- unname(method_labels[method])
  ifelse(is.na(out), method, out)
}

add_distance_metrics_to_parameter_samples <- function(parameter_samples) {
  required_cols <- c("config_id", "description", "n_col", "n_row", "a", "b",
                     "sigma", "rep", "method", "true_rho", "true_theta",
                     "pred_rho", "pred_theta")
  missing_cols <- setdiff(required_cols, names(parameter_samples))
  if (length(missing_cols) > 0) {
    stop("Р’ С„Р°Р№Р»Рµ СЃ РїР°СЂР°РјРµС‚СЂР°РјРё РЅРµ С…РІР°С‚Р°РµС‚ СЃС‚РѕР»Р±С†РѕРІ: ", paste(missing_cols, collapse = ", "))
  }

  rows <- vector("list", nrow(parameter_samples))

  for (i in seq_len(nrow(parameter_samples))) {
    cur <- parameter_samples[i, ]
    reconstructed_line <- line_from_rho_theta(
      rho = cur$pred_rho,
      theta = cur$pred_theta,
      n_row = cur$n_row,
      n_col = cur$n_col
    )

    distance_metrics <- distance_metrics_to_true_line(
      reconstructed_line,
      a_true = cur$a,
      b_true = cur$b
    )

    rows[[i]] <- data.frame(
      config_id = cur$config_id,
      description = cur$description,
      n_col = cur$n_col,
      n_row = cur$n_row,
      a = cur$a,
      b = cur$b,
      sigma = cur$sigma,
      rep = cur$rep,
      method = cur$method,
      method_label = method_label(cur$method),
      true_rho = cur$true_rho,
      true_theta = cur$true_theta,
      pred_rho = cur$pred_rho,
      pred_theta = cur$pred_theta,
      processed_active_pixels = if ("active_pixels" %in% names(cur)) cur$active_pixels else NA_real_,
      reconstructed_active_pixels = as.numeric(distance_metrics["active_pixels"]),
      reconstructed_sum_distance = as.numeric(distance_metrics["sum_distance"]),
      reconstructed_mean_distance = as.numeric(distance_metrics["mean_distance"]),
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(rows)
}

load_existing_parameter_samples <- function(sigma_values) {
  known_files <- data.frame(
    sigma = c(0.2, 0.05),
    file = c(
      file.path(line_experiments_dir, "line_experiment3_ht_error_samples.csv"),
      file.path(line_experiments_dir, "line_experiment3_sigma_005_ht_error_samples.csv")
    ),
    stringsAsFactors = FALSE
  )

  requested <- known_files[known_files$sigma %in% sigma_values, , drop = FALSE]
  if (nrow(requested) != length(sigma_values) || any(!file.exists(requested$file))) {
    return(NULL)
  }

  dplyr::bind_rows(lapply(requested$file, read.csv, stringsAsFactors = FALSE))
}


## ----experiment-parameters----------------------------------------------------
num_of_lines <- 1L
detector_q <- 0.9
rho_step_ht <- 1
theta_step_ht <- 0.01
n_rep <- as.integer(Sys.getenv("LINE_DISTANCE_N_REP", "100"))
sigma_values <- as.numeric(strsplit(Sys.getenv("LINE_DISTANCE_SIGMA", "0.2,0.05"), ",")[[1]])
use_nms <- FALSE

configs <- data.frame(
  config_id = c(
    "diag_pos_full",
    "diag_neg_full",
    "steep_full",
    "steep_shifted"
  ),
  config_label = c(
    "РџРѕР»РѕР¶РёС‚РµР»СЊРЅР°СЏ РґРёР°РіРѕРЅР°Р»СЊ",
    "РћС‚СЂРёС†Р°С‚РµР»СЊРЅР°СЏ РґРёР°РіРѕРЅР°Р»СЊ",
    "РљСЂСѓС‚Р°СЏ РїСЂСЏРјР°СЏ",
    "РљСЂСѓС‚Р°СЏ СЃРјРµС‰РµРЅРЅР°СЏ"
  ),
  n_col = rep(100L, 4),
  n_row = rep(100L, 4),
  a = c(1, -1, 2, 2),
  b = c(0, 101, -1, -40),
  description = c(
    "Positive full diagonal",
    "Negative full diagonal",
    "Steep line crossing most of the image",
    "Steep shifted line with shorter visible segment"
  ),
  stringsAsFactors = FALSE
)

configs


## ----detector-list------------------------------------------------------------
detectors <- list(
  cssa_row_row = function(m, sigma_noise) {
    cssa_denoise(m, num_of_lines = num_of_lines, method = "row.row")
  },
  median = function(m, sigma_noise) {
    median_denoise(m, n = 3L)
  },
  wiener = function(m, sigma_noise) {
    wiener_denoise(m, ksize = 5L)
  },
  quantile_09 = function(m, sigma_noise) {
    quantile_denoise(m, q = detector_q)
  }
)

method_labels <- c(
  cssa_row_row = "CSSA ROW-ROW",
  bm3d = "BM3D",
  median = "Median",
  wiener = "Wiener",
  quantile_09 = "Quantile 0.9"
)

method_order <- names(method_labels)

if (bm3d_available) {
  detectors <- append(
    detectors,
    list(
      bm3d = function(m, sigma_noise) {
        bm3d_denoise(m, sigma_noise = sigma_noise)
      }
    ),
    after = 1L
  )
} else {
  warning("BM3D РїСЂРѕРїСѓС‰РµРЅ: Python-РјРѕРґСѓР»СЊ 'bm3d' РЅРµ РЅР°Р№РґРµРЅ.")
}

method_levels <- unname(method_labels[names(detectors)])
names(detectors)


## ----experiment-function------------------------------------------------------
run_rt_distance_experiment <- function(
  configs,
  detectors,
  sigma_values,
  num_of_lines,
  rho_step_ht,
  theta_step_ht,
  n_rep,
  use_nms = FALSE
) {
  total_rows <- length(sigma_values) * nrow(configs) * n_rep * length(detectors)
  results_list <- vector("list", length = total_rows)
  idx_out <- 1L

  for (sigma_noise in sigma_values) {
    for (cfg_i in seq_len(nrow(configs))) {
      cfg <- configs[cfg_i, ]
      base_matrix <- make_line_image(
        n_row = cfg$n_row,
        n_col = cfg$n_col,
        a = cfg$a,
        b = cfg$b
      )
      true_line <- matrix(convert_ab_to_rho_theta(cfg$a, cfg$b), nrow = 1)

      for (rep_i in seq_len(n_rep)) {
        noisy_matrix <- add.noise(base_matrix, sigma = sigma_noise)

        for (method_name in names(detectors)) {
          processed_matrix <- clip01(detectors[[method_name]](
            noisy_matrix,
            sigma_noise = sigma_noise
          ))

          ht_result <- make_accumulator(
            noisy_matrix,
            detector = function(m) processed_matrix,
            rho_step = rho_step_ht,
            theta_step = theta_step_ht
          )

          pred_line <- find_k_max(
            ht_result$accumulator,
            k = num_of_lines,
            qrho = ht_result$rho,
            qtheta = ht_result$theta,
            suppress = use_nms,
            window = 6L
          )

          reconstructed_line <- line_from_rho_theta(
            rho = pred_line[1, 1],
            theta = pred_line[1, 2],
            n_row = cfg$n_row,
            n_col = cfg$n_col
          )

          distance_metrics <- distance_metrics_to_true_line(
            reconstructed_line,
            a_true = cfg$a,
            b_true = cfg$b
          )

          results_list[[idx_out]] <- data.frame(
            config_id = cfg$config_id,
            config_label = cfg$config_label,
            description = cfg$description,
            n_col = cfg$n_col,
            n_row = cfg$n_row,
            a = cfg$a,
            b = cfg$b,
            sigma = sigma_noise,
            rep = rep_i,
            method = method_name,
            method_label = method_labels[[method_name]],
            true_rho = true_line[1, 1],
            true_theta = true_line[1, 2],
            pred_rho = pred_line[1, 1],
            pred_theta = pred_line[1, 2],
            processed_active_pixels = sum(processed_matrix > 0),
            reconstructed_active_pixels = as.numeric(distance_metrics["active_pixels"]),
            reconstructed_sum_distance = as.numeric(distance_metrics["sum_distance"]),
            reconstructed_mean_distance = as.numeric(distance_metrics["mean_distance"]),
            stringsAsFactors = FALSE
          )

          idx_out <- idx_out + 1L
        }
      }
    }
  }

  do.call(rbind, results_list)
}


## ----run-experiment-----------------------------------------------------------
samples_csv <- file.path(line_experiments_dir, "line_experiment_rt_distance_samples.csv")
force_recalculate <- tolower(Sys.getenv("LINE_DISTANCE_FORCE_RECALCULATE", "false")) == "true"

recalculate <- TRUE
if (file.exists(samples_csv)) {
  distance_samples <- read.csv(samples_csv, stringsAsFactors = FALSE)
  recalculate <- !(
    setequal(unique(distance_samples$sigma), sigma_values) &&
      setequal(unique(distance_samples$config_id), configs$config_id) &&
      "reconstructed_mean_distance" %in% names(distance_samples)
  )
}

if (recalculate && !force_recalculate) {
  parameter_samples <- load_existing_parameter_samples(sigma_values)
  if (!is.null(parameter_samples)) {
    distance_samples <- add_distance_metrics_to_parameter_samples(parameter_samples)
    write.csv(distance_samples, samples_csv, row.names = FALSE)
    recalculate <- FALSE
  }
}

if (recalculate) {
  distance_samples <- run_rt_distance_experiment(
    configs = configs,
    detectors = detectors,
    sigma_values = sigma_values,
    num_of_lines = num_of_lines,
    rho_step_ht = rho_step_ht,
    theta_step_ht = theta_step_ht,
    n_rep = n_rep,
    use_nms = use_nms
  )

  write.csv(distance_samples, samples_csv, row.names = FALSE)
}

distance_samples$config_label <- configs$config_label[
  match(distance_samples$config_id, configs$config_id)
]

unknown_methods <- setdiff(unique(distance_samples$method), method_order)
method_order_current <- c(method_order[method_order %in% unique(distance_samples$method)], unknown_methods)
method_levels_current <- method_label(method_order_current)

head(distance_samples, 12)


## ----summary------------------------------------------------------------------
distance_summary <- distance_samples |>
  dplyr::group_by(config_id, config_label, description, sigma, method, method_label) |>
  dplyr::summarise(
    mean_distance = mean(reconstructed_mean_distance, na.rm = TRUE),
    median_distance = median(reconstructed_mean_distance, na.rm = TRUE),
    point_weighted_mean_distance = sum(reconstructed_sum_distance, na.rm = TRUE) /
      sum(reconstructed_active_pixels, na.rm = TRUE),
    mean_processed_active_pixels = mean(processed_active_pixels, na.rm = TRUE),
    mean_reconstructed_active_pixels = mean(reconstructed_active_pixels, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    config_id = factor(config_id, levels = configs$config_id),
    config_label = factor(config_label, levels = configs$config_label),
    method = factor(method, levels = method_order_current),
    method_label = factor(method_label, levels = method_levels_current)
  ) |>
  dplyr::arrange(sigma, config_id, method)

distance_summary


## ----xtabs-mean-distance------------------------------------------------------
mean_distance_tables <- lapply(
  split(distance_summary, distance_summary$sigma),
  function(df) xtabs(mean_distance ~ config_id + method, data = df)
)
mean_distance_tables

mean_distance_tables_norm <- lapply(mean_distance_tables, normalize_xtab)
mean_distance_tables_norm


## ----xtabs-median-distance----------------------------------------------------
median_distance_tables <- lapply(
  split(distance_summary, distance_summary$sigma),
  function(df) xtabs(median_distance ~ config_id + method, data = df)
)
median_distance_tables

median_distance_tables_norm <- lapply(median_distance_tables, normalize_xtab)
median_distance_tables_norm


## ----save-tables--------------------------------------------------------------
slides_tables_dir <- file.path(project_root, "VKR", "VKR_SLIDES", "tables")
dir.create(slides_tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(line_experiments_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  distance_samples,
  samples_csv,
  row.names = FALSE
)

write.csv(
  distance_summary,
  file.path(line_experiments_dir, "line_experiment_rt_distance_summary.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    params = list(
      sigma_values = sigma_values,
      n_rep = n_rep,
      rho_step_ht = rho_step_ht,
      theta_step_ht = theta_step_ht,
      use_nms = use_nms,
      detectors = names(detectors)
    ),
    distance_samples = distance_samples,
    distance_summary = distance_summary
  ),
  file.path(line_experiments_dir, "line_experiment_rt_distance_results.rds")
)

write.csv(
  distance_samples,
  file.path(slides_tables_dir, "line_distance_samples.csv"),
  row.names = FALSE
)

write.csv(
  distance_summary,
  file.path(slides_tables_dir, "line_distance_summary.csv"),
  row.names = FALSE
)

config_levels <- configs$config_label

for (sigma_value in sigma_values) {
  tag <- sigma_tag(sigma_value)

  mean_df <- wide_metric_table(
    distance_summary,
    sigma_value = sigma_value,
    metric = "mean_distance",
    config_levels = config_levels,
    method_levels = method_levels_current
  )

  median_df <- wide_metric_table(
    distance_summary,
    sigma_value = sigma_value,
    metric = "median_distance",
    config_levels = config_levels,
    method_levels = method_levels_current
  )

  point_weighted_mean_df <- wide_metric_table(
    distance_summary,
    sigma_value = sigma_value,
    metric = "point_weighted_mean_distance",
    config_levels = config_levels,
    method_levels = method_levels_current
  )

  mean_norm_df <- wide_normalized_metric_table(
    distance_summary,
    sigma_value = sigma_value,
    metric = "mean_distance",
    config_levels = config_levels,
    method_levels = method_levels_current
  )

  median_norm_df <- wide_normalized_metric_table(
    distance_summary,
    sigma_value = sigma_value,
    metric = "median_distance",
    config_levels = config_levels,
    method_levels = method_levels_current
  )

  write_df_tex(
    mean_df,
    file.path(slides_tables_dir, paste0("line_distance_mean_sigma_", tag, ".tex")),
    digits = 4L
  )

  write_df_tex(
    median_df,
    file.path(slides_tables_dir, paste0("line_distance_median_sigma_", tag, ".tex")),
    digits = 4L
  )

  write_df_tex(
    point_weighted_mean_df,
    file.path(slides_tables_dir, paste0("line_distance_point_weighted_mean_sigma_", tag, ".tex")),
    digits = 4L
  )

  write_df_tex(
    mean_norm_df,
    file.path(slides_tables_dir, paste0("line_distance_mean_norm_sigma_", tag, ".tex")),
    digits = 4L
  )

  write_df_tex(
    median_norm_df,
    file.path(slides_tables_dir, paste0("line_distance_median_norm_sigma_", tag, ".tex")),
    digits = 4L
  )
}

list.files(slides_tables_dir, pattern = "^line_distance", full.names = TRUE)


