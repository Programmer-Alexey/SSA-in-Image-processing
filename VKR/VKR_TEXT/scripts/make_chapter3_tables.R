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
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(project_root, ".."), winslash = "/", mustWork = TRUE)
setwd(project_root)

project_lib <- normalizePath(file.path(repo_root, ".r-local-lib"),
                             winslash = "/", mustWork = FALSE)
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

required_pkgs <- c("dplyr")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop("РћС‚СЃСѓС‚СЃС‚РІСѓСЋС‚ РїР°РєРµС‚С‹: ", paste(missing_pkgs, collapse = ", "))
}

if (!requireNamespace("reticulate", quietly = TRUE)) {
  install.packages("reticulate", repos = "https://cloud.r-project.org")
}

source(file.path(repo_root, "ssa-based methods", "cssa-transform.r"))
source(file.path(repo_root, "hough transform", "hough_transform.r"))

library(dplyr)
library(reticulate)

clip01 <- function(m) {
  m <- as.matrix(m)
  m[!is.finite(m)] <- 0
  m[m < 0.001] <- 0
  m[m > 0.99] <- 1
  m
}

make_line_image <- function(n_row, n_col, a, b, intensity = 1) {
  matrix(0, nrow = n_row, ncol = n_col) |>
    add.line(a = a, b = b, intensity = intensity)
}

binary_by_quantile <- function(m, q = 0.9, use_abs = FALSE) {
  x <- as.matrix(m)
  x[!is.finite(x)] <- 0
  score <- if (use_abs) abs(x) else x
  thr <- as.numeric(stats::quantile(score, probs = q, na.rm = TRUE))
  out <- ifelse(score >= thr, clip01(x), 0)
  clip01(out)
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

bm3d_available <- reticulate::py_module_available("bm3d")
if (!bm3d_available) {
  stop("Python-РјРѕРґСѓР»СЊ 'bm3d' РЅРµРґРѕСЃС‚СѓРїРµРЅ, Р° РѕРЅ РЅСѓР¶РµРЅ РґР»СЏ РіР»Р°РІС‹ 3.")
}
bm3d_lib <- reticulate::import("bm3d")

bm3d_denoise <- function(m, sigma_noise) {
  x <- as.matrix(m)
  x[!is.finite(x)] <- 0
  clip01(bm3d_lib$bm3d(x, sigma_psd = as.numeric(sigma_noise)))
}

pad_replicate <- function(x, pad) {
  nr <- nrow(x)
  nc <- ncol(x)
  out <- matrix(0, nrow = nr + 2 * pad, ncol = nc + 2 * pad)
  row_idx <- pmin(pmax(seq_len(nr + 2 * pad) - pad, 1L), nr)
  col_idx <- pmin(pmax(seq_len(nc + 2 * pad) - pad, 1L), nc)
  out[,] <- x[row_idx, col_idx]
  out
}

median_denoise <- function(m, n = 3L) {
  x <- clip01(m)
  n <- as.integer(n)
  if (n %% 2L == 0L) {
    n <- n + 1L
  }
  pad <- (n - 1L) %/% 2L
  padded <- pad_replicate(x, pad)
  out <- matrix(0, nrow = nrow(x), ncol = ncol(x))

  for (i in seq_len(nrow(x))) {
    for (j in seq_len(ncol(x))) {
      block <- padded[i:(i + 2L * pad), j:(j + 2L * pad)]
      out[i, j] <- stats::median(block)
    }
  }

  clip01(out)
}

wiener_denoise <- function(m, ksize = 5L) {
  x <- clip01(m)
  ksize <- as.integer(ksize)
  if (ksize %% 2L == 0L) {
    ksize <- ksize + 1L
  }

  pad <- (ksize - 1L) %/% 2L
  padded <- pad_replicate(x, pad)
  local_mean <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  local_var <- matrix(0, nrow = nrow(x), ncol = ncol(x))

  for (i in seq_len(nrow(x))) {
    for (j in seq_len(ncol(x))) {
      block <- padded[i:(i + 2L * pad), j:(j + 2L * pad)]
      local_mean[i, j] <- mean(block)
      local_var[i, j] <- stats::var(as.vector(block))
    }
  }

  local_var[!is.finite(local_var)] <- 0
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

convert_rho_theta_to_ab <- function(rho, theta, eps = 1e-8) {
  if (abs(sin(theta)) < eps) {
    return(c(a = NA_real_, b = NA_real_))
  }
  c(
    a = -cos(theta) / sin(theta),
    b = rho / sin(theta)
  )
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

wide_metric_table <- function(summary_df, metric, config_levels, method_levels) {
  wide <- xtabs(summary_df[[metric]] ~ summary_df$config_label + summary_df$method_label)
  wide <- wide[config_levels, method_levels, drop = FALSE]
  out <- data.frame(
    "РљРѕРЅС„РёРіСѓСЂР°С†РёСЏ" = rownames(wide),
    as.data.frame.matrix(wide),
    check.names = FALSE,
    row.names = NULL
  )
  out
}

dir.create("tables/chapter3", recursive = TRUE, showWarnings = FALSE)
dir.create("images/chapter3", recursive = TRUE, showWarnings = FALSE)

file.copy(
  file.path(repo_root, "line_experiment3_files", "figure-html", "configuration-preview-1.png"),
  "images/chapter3/configuration-preview-1.png",
  overwrite = TRUE
)

sigma_noise <- 0.2
num_of_lines <- 1L
detector_q <- 0.9
rho_step_ht <- 1
theta_step_ht <- 0.01
n_rep <- 100L
use_nms <- FALSE

configs <- data.frame(
  config_id = c(
    "diag_pos_full",
    "diag_neg_full",
    "steep_full",
    "steep_shifted"
  ),
  config_label = c(
    "diag\\_pos\\_full",
    "diag\\_neg\\_full",
    "steep\\_full",
    "steep\\_shifted"
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

detectors <- list(
  cssa_row_row = function(m, sigma_noise) {
    cssa_denoise(m, num_of_lines = num_of_lines, method = "row.row")
  },
  cssa_row_col = function(m, sigma_noise) {
    cssa_denoise(m, num_of_lines = num_of_lines, method = "row.col")
  },
  cssa_col_row = function(m, sigma_noise) {
    cssa_denoise(m, num_of_lines = num_of_lines, method = "col.row")
  },
  cssa_col_col = function(m, sigma_noise) {
    cssa_denoise(m, num_of_lines = num_of_lines, method = "col.col")
  },
  median = function(m, sigma_noise) {
    median_denoise(m, n = 3L)
  },
  wiener = function(m, sigma_noise) {
    wiener_denoise(m, ksize = 5L)
  },
  bm3d = function(m, sigma_noise) {
    bm3d_denoise(m, sigma_noise = sigma_noise)
  },
  quantile_09 = function(m, sigma_noise) {
    quantile_denoise(m, q = detector_q)
  }
)

method_labels <- c(
  cssa_row_row = "CSSA row.row",
  cssa_row_col = "CSSA row.col",
  cssa_col_row = "CSSA col.row",
  cssa_col_col = "CSSA col.col",
  median = "Median",
  wiener = "Wiener",
  bm3d = "BM3D",
  quantile_09 = "Quantile 0.9"
)

result_csv <- "tables/chapter3/chapter3_ab_error_samples.csv"
recalculate <- !file.exists(result_csv)

if (recalculate) {
  set.seed(1111)
  results_list <- vector("list", length = nrow(configs) * n_rep * length(detectors))
  idx_out <- 1L

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

        pred_ab <- convert_rho_theta_to_ab(pred_line[1, 1], pred_line[1, 2])

        results_list[[idx_out]] <- data.frame(
          config_id = cfg$config_id,
          config_label = cfg$config_label,
          description = cfg$description,
          n_col = cfg$n_col,
          n_row = cfg$n_row,
          a_true = cfg$a,
          b_true = cfg$b,
          sigma = sigma_noise,
          rep = rep_i,
          method = method_name,
          method_label = method_labels[[method_name]],
          true_rho = true_line[1, 1],
          true_theta = true_line[1, 2],
          pred_rho = pred_line[1, 1],
          pred_theta = pred_line[1, 2],
          pred_a = pred_ab["a"],
          pred_b = pred_ab["b"],
          abs_err_a = abs(pred_ab["a"] - cfg$a),
          abs_err_b = abs(pred_ab["b"] - cfg$b),
          active_pixels = sum(processed_matrix > 0),
          stringsAsFactors = FALSE
        )

        idx_out <- idx_out + 1L
      }
    }
  }

  ab_error_samples <- do.call(rbind, results_list)
  write.csv(ab_error_samples, result_csv, row.names = FALSE)
} else {
  ab_error_samples <- read.csv(result_csv, stringsAsFactors = FALSE)
}

ab_error_summary <- ab_error_samples |>
  dplyr::group_by(config_id, config_label, description, method, method_label) |>
  dplyr::summarise(
    mean_abs_a = mean(abs_err_a, na.rm = TRUE),
    median_abs_a = median(abs_err_a, na.rm = TRUE),
    mean_abs_b = mean(abs_err_b, na.rm = TRUE),
    median_abs_b = median(abs_err_b, na.rm = TRUE),
    mean_active_pixels = mean(active_pixels, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  ab_error_summary,
  "tables/chapter3/chapter3_ab_error_summary.csv",
  row.names = FALSE
)

config_levels <- configs$config_label
method_levels <- unname(method_labels)

mean_a_df <- wide_metric_table(ab_error_summary, "mean_abs_a", config_levels, method_levels)
median_a_df <- wide_metric_table(ab_error_summary, "median_abs_a", config_levels, method_levels)
mean_b_df <- wide_metric_table(ab_error_summary, "mean_abs_b", config_levels, method_levels)
median_b_df <- wide_metric_table(ab_error_summary, "median_abs_b", config_levels, method_levels)

write_df_tex(mean_a_df, "tables/chapter3/mean_abs_a.tex", digits = 4L)
write_df_tex(median_a_df, "tables/chapter3/median_abs_a.tex", digits = 4L)
write_df_tex(mean_b_df, "tables/chapter3/mean_abs_b.tex", digits = 4L)
write_df_tex(median_b_df, "tables/chapter3/median_abs_b.tex", digits = 4L)

experiment_params_df <- data.frame(
  "РџР°СЂР°РјРµС‚СЂ" = c(
    "Р Р°Р·РјРµСЂ РёР·РѕР±СЂР°Р¶РµРЅРёСЏ",
    "Р§РёСЃР»Рѕ РїСЂСЏРјС‹С…",
    "РЈСЂРѕРІРµРЅСЊ С€СѓРјР° sigma",
    "Р§РёСЃР»Рѕ РїРѕРІС‚РѕСЂРѕРІ",
    "РЁР°Рі rho",
    "РЁР°Рі theta",
    "NMS"
  ),
  "Р—РЅР°С‡РµРЅРёРµ" = c(
    "100 x 100",
    "1",
    "0.2",
    as.character(n_rep),
    as.character(rho_step_ht),
    as.character(theta_step_ht),
    if (use_nms) "Р”Р°" else "РќРµС‚"
  ),
  check.names = FALSE
)

method_params_df <- data.frame(
  "РњРµС‚РѕРґ" = c(
    "CSSA row.row",
    "CSSA row.col",
    "CSSA col.row",
    "CSSA col.col",
    "Median",
    "Wiener",
    "BM3D",
    "Quantile 0.9"
  ),
  "РџР°СЂР°РјРµС‚СЂС‹" = c(
    "1 РєРѕРјРїРѕРЅРµРЅС‚Р°; DFT РїРѕ СЃС‚СЂРѕРєР°Рј, CSSA РїРѕ СЃС‚СЂРѕРєР°Рј, РѕР±СЂР°С‚РЅРѕРµ DFT РїРѕ СЃС‚СЂРѕРєР°Рј",
    "1 РєРѕРјРїРѕРЅРµРЅС‚Р°; DFT РїРѕ СЃС‚СЂРѕРєР°Рј, CSSA РїРѕ СЃС‚СЂРѕРєР°Рј, РѕР±СЂР°С‚РЅРѕРµ DFT РїРѕ СЃС‚РѕР»Р±С†Р°Рј",
    "1 РєРѕРјРїРѕРЅРµРЅС‚Р°; DFT РїРѕ СЃС‚СЂРѕРєР°Рј, CSSA РїРѕ СЃС‚РѕР»Р±С†Р°Рј, РѕР±СЂР°С‚РЅРѕРµ DFT РїРѕ СЃС‚СЂРѕРєР°Рј",
    "1 РєРѕРјРїРѕРЅРµРЅС‚Р°; DFT РїРѕ СЃС‚СЂРѕРєР°Рј, CSSA РїРѕ СЃС‚РѕР»Р±С†Р°Рј, РѕР±СЂР°С‚РЅРѕРµ DFT РїРѕ СЃС‚РѕР»Р±С†Р°Рј",
    "РћРєРЅРѕ 3 x 3",
    "Р›РѕРєР°Р»СЊРЅРѕРµ РѕРєРЅРѕ 5 x 5",
    "sigma_psd = 0.2",
    "РџРѕСЂРѕРі РїРѕ РєРІР°РЅС‚РёР»СЋ 0.9"
  ),
  check.names = FALSE
)

write_df_tex(experiment_params_df, "tables/chapter3/experiment_params.tex", digits = 4L)
write_df_tex(method_params_df, "tables/chapter3/method_params.tex", digits = 4L)

