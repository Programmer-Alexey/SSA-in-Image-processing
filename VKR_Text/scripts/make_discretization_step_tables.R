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

project_lib <- normalizePath(
  file.path(repo_root, ".r-local-lib"),
  winslash = "/",
  mustWork = FALSE
)
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

required_pkgs <- c("Rssa", "Rcpp", "dplyr")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop("Отсутствуют пакеты: ", paste(missing_pkgs, collapse = ", "))
}

source(file.path(repo_root, "ssa-based methods", "cssa-transform.r"))
source(file.path(repo_root, "hough transform", "hough_transform.r"))

library(dplyr)

clip01 <- function(m) {
  m <- as.matrix(m)
  m[!is.finite(m)] <- 0
  m[m < 0.1] <- 0
  m[m > 0.99] <- 1
  m
}

make_line_image <- function(n_row, n_col, a, b, intensity = 1) {
  matrix(0, nrow = n_row, ncol = n_col) |>
    add.line(a = a, b = b, intensity = intensity)
}

cssa_denoise <- function(m, num_of_lines = 1L, method = "row.row") {
  cleaned <- switch(
    method,
    "row.row" = m |> dft() |> cssa.row(num.line = num_of_lines) |> idft.row() |> Re(),
    stop("Поддерживается только метод row.row")
  )
  clip01(cleaned)
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

compute_err <- function(true, pred) {
  true <- as.matrix(true)
  pred <- as.matrix(pred)

  pair_err <- function(tr, tt, pr, pt) {
    pt_alt <- (pt + pi) %% (2 * pi)
    dtheta <- atan2(sin(tt - pt), cos(tt - pt))^2
    dtheta_alt <- atan2(sin(tt - pt_alt), cos(tt - pt_alt))^2

    loss_1 <- (tr - pr)^2 + dtheta
    loss_2 <- (tr + pr)^2 + dtheta_alt

    if (loss_1 <= loss_2) {
      c(dr = (tr - pr)^2, dtheta = dtheta, total = loss_1)
    } else {
      c(dr = (tr + pr)^2, dtheta = dtheta_alt, total = loss_2)
    }
  }

  out <- pair_err(true[1, 1], true[1, 2], pred[1, 1], pred[1, 2])
  c(dr = unname(out[1]), dtheta = unname(out[2]))
}

format_num <- function(x, digits = 6L) {
  format(round(x, digits), nsmall = digits, scientific = FALSE, trim = TRUE)
}

underline_min <- function(x, digits = 6L) {
  out <- format_num(x, digits)
  min_x <- min(x, na.rm = TRUE)
  idx <- which(abs(x - min_x) < 1e-12)
  out[idx] <- paste0("\\underline{", out[idx], "}")
  out
}

write_df_tex_raw <- function(df, file, align = NULL) {
  if (is.null(align)) {
    align <- paste(rep("c", ncol(df)), collapse = "|")
  }
  header <- paste(names(df), collapse = " & ")
  body <- apply(df, 1, function(row) paste(row, collapse = " & "))

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

dir.create("tables", recursive = TRUE, showWarnings = FALSE)

sigma_noise <- 0.05
n_rep <- 100L
n_row <- 100L
n_col <- 100L
a_true <- 2
b_true <- -1

grid_df <- data.frame(
  grid_id = c("coarse", "medium", "fine", "very_fine"),
  grid_label = c("Крупный шаг", "Средний шаг", "Мелкий шаг", "Очень мелкий шаг"),
  rho_step = c(5, 1, 0.5, 0.25),
  theta_step = c(1, 0.01, 0.005, 0.0025),
  stringsAsFactors = FALSE
)

set.seed(1111)
base_matrix <- make_line_image(n_row = n_row, n_col = n_col, a = a_true, b = b_true)
true_line <- matrix(convert_ab_to_rho_theta(a_true, b_true), nrow = 1)

results_list <- vector("list", length = n_rep * nrow(grid_df))
idx_out <- 1L

for (rep_i in seq_len(n_rep)) {
  noisy_matrix <- add.noise(base_matrix, sigma = sigma_noise)
  processed_matrix <- cssa_denoise(noisy_matrix, num_of_lines = 1L, method = "row.row")

  for (grid_i in seq_len(nrow(grid_df))) {
    grid <- grid_df[grid_i, ]

    ht_result <- make_accumulator(
      noisy_matrix,
      detector = function(m) processed_matrix,
      rho_step = grid$rho_step,
      theta_step = grid$theta_step
    )

    pred_line <- find_k_max(
      ht_result$accumulator,
      k = 1L,
      qrho = ht_result$rho,
      qtheta = ht_result$theta,
      suppress = FALSE,
      window = 0L
    )

    err <- compute_err(true_line, pred_line)

    results_list[[idx_out]] <- data.frame(
      rep = rep_i,
      grid_id = grid$grid_id,
      grid_label = grid$grid_label,
      rho_step = grid$rho_step,
      theta_step = grid$theta_step,
      pred_rho = pred_line[1, 1],
      pred_theta = pred_line[1, 2],
      dr = as.numeric(err["dr"]),
      dtheta = as.numeric(err["dtheta"]),
      stringsAsFactors = FALSE
    )

    idx_out <- idx_out + 1L
  }
}

samples_df <- do.call(rbind, results_list)
write.csv(
  samples_df,
  "tables/discretization_step_samples.csv",
  row.names = FALSE
)

summary_df <- samples_df |>
  dplyr::group_by(grid_id, grid_label, rho_step, theta_step) |>
  dplyr::summarise(
    mean_dr = mean(dr),
    median_dr = median(dr),
    mean_dtheta = mean(dtheta),
    median_dtheta = median(dtheta),
    mean_pred_rho = mean(pred_rho),
    mean_pred_theta = mean(pred_theta),
    .groups = "drop"
  ) |>
  dplyr::arrange(match(grid_id, grid_df$grid_id))

write.csv(
  summary_df,
  "tables/discretization_step_summary.csv",
  row.names = FALSE
)

params_df <- data.frame(
  "Параметр" = c(
    "Изображение",
    "Прямая",
    "Уровень шума",
    "Предобработка",
    "Число повторов"
  ),
  "Значение" = c(
    "100 x 100",
    "$y = 2x - 1$",
    "$\\sigma = 0.05$",
    "CSSA row.row, 1 компонента",
    as.character(n_rep)
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

summary_tex_df <- data.frame(
  "Сетка" = summary_df$grid_label,
  "$h_\\rho$" = format_num(summary_df$rho_step, digits = 4L),
  "$h_\\theta$" = format_num(summary_df$theta_step, digits = 4L),
  "Средняя $\\Delta_\\rho$" = underline_min(summary_df$mean_dr, digits = 6L),
  "Медианная $\\Delta_\\rho$" = underline_min(summary_df$median_dr, digits = 6L),
  "Средняя $\\Delta_\\theta$" = underline_min(summary_df$mean_dtheta, digits = 6L),
  "Медианная $\\Delta_\\theta$" = underline_min(summary_df$median_dtheta, digits = 6L),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write_df_tex_raw(
  params_df,
  "tables/discretization_step_params.tex",
  align = "l|c"
)

write_df_tex_raw(
  summary_tex_df,
  "tables/discretization_step_summary.tex",
  align = "l|c|c|c|c|c|c"
)

cat("Готово. Таблицы сохранены в папку tables.\n")
