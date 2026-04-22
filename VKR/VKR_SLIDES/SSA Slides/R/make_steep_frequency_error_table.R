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

if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("РќСѓР¶РµРЅ РїР°РєРµС‚ dplyr.")
}

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

cssa_rank1_row <- function(row_vec) {
  row_vec |>
    matrix(nrow = 1L) |>
    dft() |>
    cssa.row(num.line = 1L) |>
    idft.row() |>
    Re() |>
    as.numeric()
}

row_metrics <- function(est, true_col) {
  clean <- numeric(length(est))
  clean[true_col] <- 1
  pos <- pmax(est, 0)
  near_idx <- max(1L, true_col - 1L):min(length(est), true_col + 1L)
  pos_sum <- sum(pos)

  data.frame(
    argmax_error = abs(which.max(est) - true_col),
    mse = mean((est - clean)^2),
    near_mass_share = if (pos_sum > 0) sum(pos[near_idx]) / pos_sum else NA_real_,
    active_share_0001 = mean(est > 0.001),
    active_share_01 = mean(est > 0.1),
    stringsAsFactors = FALSE
  )
}

run_config <- function(config_id, a, b, sigma = 0.2, n_rep = 100L, n = 100L) {
  signal_matrix <- make_line_image(n, n, a = a, b = b)
  signal_points <- which(signal_matrix > 0, arr.ind = TRUE)
  # Р”Р»СЏ РІС‹Р±СЂР°РЅРЅС‹С… РєСЂСѓС‚С‹С… РїСЂСЏРјС‹С… РІ РєР°Р¶РґРѕР№ СЃС‚СЂРѕРєРµ СЂРѕРІРЅРѕ РѕРґРёРЅ РїРёРєСЃРµР»СЊ СЃРёРіРЅР°Р»Р°.
  signal_rows <- split(signal_points[, "col"], signal_points[, "row"])

  rows <- vector("list", length(signal_rows) * n_rep)
  out_i <- 1L

  set.seed(1111)
  for (rep_i in seq_len(n_rep)) {
    noisy_matrix <- add.noise(signal_matrix, sigma = sigma)

    for (row_name in names(signal_rows)) {
      row_id <- as.integer(row_name)
      true_col <- as.integer(signal_rows[[row_name]][1])
      est <- cssa_rank1_row(noisy_matrix[row_id, ])
      metrics <- row_metrics(est, true_col)
      omega <- -2 * pi * (true_col - 1L) / n
      omega_wrapped <- atan2(sin(omega), cos(omega))

      rows[[out_i]] <- data.frame(
        config_id = config_id,
        rep = rep_i,
        row_id = row_id,
        true_col = true_col,
        omega = omega,
        omega_wrapped = omega_wrapped,
        metrics,
        stringsAsFactors = FALSE
      )
      out_i <- out_i + 1L
    }
  }

  dplyr::bind_rows(rows)
}

format_range <- function(x) {
  sprintf("%d--%d", min(x), max(x))
}

write_df_tex <- function(df, file, digits = 4L) {
  df_out <- df
  for (j in seq_along(df_out)) {
    if (is.numeric(df_out[[j]])) {
      df_out[[j]] <- format(round(df_out[[j]], digits), nsmall = digits, trim = TRUE)
    } else {
      df_out[[j]] <- as.character(df_out[[j]])
      df_out[[j]] <- gsub("_", "\\\\_", df_out[[j]])
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

slides_dir <- file.path(project_root, "VKR", "VKR_SLIDES", "SSA Slides")
dir.create(file.path(slides_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

n_rep <- as.integer(Sys.getenv("STEEP_FREQUENCY_N_REP", "100"))
sigma <- as.numeric(Sys.getenv("STEEP_FREQUENCY_SIGMA", "0.2"))

samples <- dplyr::bind_rows(
  run_config("steep_full", a = 2, b = -1, sigma = sigma, n_rep = n_rep),
  run_config("steep_shifted", a = 2, b = -40, sigma = sigma, n_rep = n_rep)
)

samples <- samples |>
  dplyr::mutate(
    frequency_band = cut(
      true_col,
      breaks = seq(0, 100, by = 10),
      labels = paste(seq(1, 91, by = 10), seq(10, 100, by = 10), sep = "--"),
      include.lowest = TRUE
    )
  )

summary <- samples |>
  dplyr::group_by(config_id, frequency_band) |>
  dplyr::summarise(
    true_col_range = format_range(true_col),
    omega_wrapped_min = min(omega_wrapped),
    omega_wrapped_max = max(omega_wrapped),
    n_rows = dplyr::n_distinct(row_id),
    mean_argmax_error = mean(argmax_error),
    median_argmax_error = median(argmax_error),
    mean_mse = mean(mse),
    mean_near_mass_share = mean(near_mass_share, na.rm = TRUE),
    mean_active_share_0001 = mean(active_share_0001),
    mean_active_share_01 = mean(active_share_01),
    .groups = "drop"
  )

write.csv(samples, file.path(slides_dir, "tables", "steep_frequency_error_samples.csv"), row.names = FALSE)
write.csv(summary, file.path(slides_dir, "tables", "steep_frequency_error_summary.csv"), row.names = FALSE)

slide_table <- summary |>
  dplyr::mutate(
    config_id = dplyr::recode(
      config_id,
      steep_full = "РєСЂСѓС‚Р°СЏ",
      steep_shifted = "РєСЂСѓС‚Р°СЏ СЃРјРµС‰."
    )
  ) |>
  dplyr::select(
    config_id,
    true_col_range,
    omega_wrapped_min,
    omega_wrapped_max,
    mean_argmax_error,
    mean_near_mass_share,
    mean_active_share_0001,
    mean_active_share_01
  )
names(slide_table) <- c(
  "РљРѕРЅС„.",
  "$j$",
  "$\\omega_{\\min}$",
  "$\\omega_{\\max}$",
  "$|\\hat j-j|$",
  "РґРѕР»СЏ РјР°СЃСЃС‹ $j\\pm1$",
  "$\\widehat x>10^{-3}$",
  "$\\widehat x>0.1$"
)

write_df_tex(slide_table, file.path(slides_dir, "tables", "steep_frequency_error_summary.tex"), digits = 3L)

cat("Р“РѕС‚РѕРІРѕ: tables/steep_frequency_error_summary.csv Рё .tex\n")

