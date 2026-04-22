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

source(file.path(project_root, "ssa-based methods", "cssa-transform.r"))
source(file.path(project_root, "hough transform", "hough_transform.r"))

slides_dir <- file.path(project_root, "VKR", "VKR_SLIDES", "SSA Slides")
tables_dir <- file.path(slides_dir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

make_line_image <- function(n_row, n_col, a, b, intensity = 1) {
  matrix(0, nrow = n_row, ncol = n_col) |>
    add.line(a = a, b = b, intensity = intensity)
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

distance_to_true_line <- function(line_matrix, a_true, b_true) {
  pts <- which(as.matrix(line_matrix) > 0, arr.ind = TRUE)
  if (nrow(pts) == 0) {
    return(NA_real_)
  }
  y <- pts[, "row"]
  x <- pts[, "col"]
  mean(abs(a_true * x - y + b_true) / sqrt(a_true^2 + 1))
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

run_ideal_ht <- function(clean_matrix, cfg, rho_step, theta_step) {
  ht_result <- make_accumulator(
    clean_matrix,
    detector = function(m) clean_matrix,
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
  distance_to_true_line(reconstructed_line, a_true = cfg$a, b_true = cfg$b)
}

run_ideal_param_errors <- function(clean_matrix, cfg, rho_step, theta_step) {
  ht_result <- make_accumulator(
    clean_matrix,
    detector = function(m) clean_matrix,
    rho_step = rho_step,
    theta_step = theta_step
  )
  pred_line <- find_k_max(ht_result$accumulator, k = 1L, qrho = ht_result$rho, qtheta = ht_result$theta)
  true_line <- convert_ab_to_rho_theta(cfg$a, cfg$b)

  c(
    rho_mse = unname((pred_line[1, 1] - true_line["rho"])^2),
    theta_mse = unname((pred_line[1, 2] - true_line["theta"])^2)
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

configs <- data.frame(
  config_id = c("diag_pos_full", "diag_neg_full", "steep_full", "steep_shifted"),
  config_label = c("РџРѕР»РѕР¶РёС‚РµР»СЊРЅР°СЏ РґРёР°РіРѕРЅР°Р»СЊ", "РћС‚СЂРёС†Р°С‚РµР»СЊРЅР°СЏ РґРёР°РіРѕРЅР°Р»СЊ", "РљСЂСѓС‚Р°СЏ РїСЂСЏРјР°СЏ", "РљСЂСѓС‚Р°СЏ СЃРјРµС‰РµРЅРЅР°СЏ"),
  n_col = rep(100L, 4),
  n_row = rep(100L, 4),
  a = c(1, -1, 2, 2),
  b = c(0, 101, -1, -40),
  stringsAsFactors = FALSE
)

grid_specs <- data.frame(
  grid_label = c("$2; 0.02$", "$1; 0.01$", "$0.5; 0.005$"),
  rho_step = c(2, 1, 0.5),
  theta_step = c(0.02, 0.01, 0.005),
  stringsAsFactors = FALSE
)

ideal_grid_rows <- vector("list", nrow(configs) * nrow(grid_specs))
idx <- 1L
ideal_threshold_rows <- vector("list", nrow(configs))
ideal_parameter_rows <- vector("list", nrow(configs))

for (cfg_i in seq_len(nrow(configs))) {
  cfg <- configs[cfg_i, ]
  clean_matrix <- make_line_image(cfg$n_row, cfg$n_col, cfg$a, cfg$b)
  param_errors <- run_ideal_param_errors(clean_matrix, cfg, rho_step = 1, theta_step = 0.01)

  ideal_threshold_rows[[cfg_i]] <- data.frame(
    config_id = cfg$config_id,
    config_label = cfg$config_label,
    mean_distance = run_ideal_ht(clean_matrix, cfg, rho_step = 1, theta_step = 0.01),
    stringsAsFactors = FALSE
  )

  ideal_parameter_rows[[cfg_i]] <- data.frame(
    config_id = cfg$config_id,
    config_label = cfg$config_label,
    rho_mse = as.numeric(param_errors["rho_mse"]),
    theta_mse = as.numeric(param_errors["theta_mse"]),
    stringsAsFactors = FALSE
  )

  for (grid_i in seq_len(nrow(grid_specs))) {
    grid <- grid_specs[grid_i, ]
    ideal_grid_rows[[idx]] <- data.frame(
      config_id = cfg$config_id,
      config_label = cfg$config_label,
      grid_label = grid$grid_label,
      rho_step = grid$rho_step,
      theta_step = grid$theta_step,
      mean_distance = run_ideal_ht(clean_matrix, cfg, grid$rho_step, grid$theta_step),
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
}

ideal_grid <- do.call(rbind, ideal_grid_rows)
ideal_threshold <- do.call(rbind, ideal_threshold_rows)
ideal_parameters <- do.call(rbind, ideal_parameter_rows)

write.csv(ideal_grid, file.path(tables_dir, "ideal_grid_step_mean_distance.csv"), row.names = FALSE)
write.csv(ideal_threshold, file.path(tables_dir, "ideal_default_grid_mean_distance.csv"), row.names = FALSE)
write.csv(ideal_parameters, file.path(tables_dir, "ideal_default_grid_parameter_errors.csv"), row.names = FALSE)

ideal_grid_summary <- aggregate(mean_distance ~ grid_label, data = ideal_grid, FUN = mean)
ideal_grid_summary <- ideal_grid_summary[match(grid_specs$grid_label, ideal_grid_summary$grid_label), ]
ideal_grid_df <- data.frame(
  "РЎР»СѓС‡Р°Р№" = "Р±РµР· С€СѓРјР°",
  t(as.matrix(ideal_grid_summary$mean_distance)),
  check.names = FALSE
)
names(ideal_grid_df) <- c("РЎР»СѓС‡Р°Р№", grid_specs$grid_label)
write_df_tex(ideal_grid_df, file.path(tables_dir, "ideal_grid_step_mean_distance.tex"), digits = 4L)

ideal_threshold_wide <- data.frame(
  "РЎР»СѓС‡Р°Р№" = "Р±РµР· С€СѓРјР°",
  t(as.matrix(ideal_threshold$mean_distance)),
  check.names = FALSE
)
names(ideal_threshold_wide) <- c("РЎР»СѓС‡Р°Р№", configs$config_label)
write_df_tex(ideal_threshold_wide, file.path(tables_dir, "ideal_default_grid_mean_distance.tex"), digits = 4L)

ideal_rho_wide <- data.frame(
  "РЎР»СѓС‡Р°Р№" = "Р±РµР· С€СѓРјР°",
  t(as.matrix(ideal_parameters$rho_mse)),
  check.names = FALSE
)
names(ideal_rho_wide) <- c("РЎР»СѓС‡Р°Р№", configs$config_label)
write_df_tex(ideal_rho_wide, file.path(tables_dir, "ideal_default_grid_rho_mse.tex"), digits = 4L)

ideal_theta_wide <- data.frame(
  "РЎР»СѓС‡Р°Р№" = "Р±РµР· С€СѓРјР°",
  t(as.matrix(ideal_parameters$theta_mse)),
  check.names = FALSE
)
names(ideal_theta_wide) <- c("РЎР»СѓС‡Р°Р№", configs$config_label)
write_df_tex(ideal_theta_wide, file.path(tables_dir, "ideal_default_grid_theta_mse.tex"), digits = 6L)

cat("Р“РѕС‚РѕРІРѕ: ideal_* С‚Р°Р±Р»РёС†С‹ РґР»СЏ РіРµРѕРјРµС‚СЂРёС‡РµСЃРєРѕР№ Рё РїР°СЂР°РјРµС‚СЂРёС‡РµСЃРєРѕР№ РѕС€РёР±РєРё\n")

