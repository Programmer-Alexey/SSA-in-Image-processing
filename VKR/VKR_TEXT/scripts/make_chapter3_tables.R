script_dir_for_source <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

source(file.path(script_dir_for_source(), "00_common.R"))

ensure_report_dirs()
load_project_code(load_hough = TRUE)

comparison_n_row <- as.integer(Sys.getenv("VKR_N_ROW", "100"))
comparison_n_col <- as.integer(Sys.getenv("VKR_N_COL", "100"))
comparison_n_rep <- as.integer(Sys.getenv("VKR_N_REP", "100"))
comparison_sigmas <- as.numeric(strsplit(Sys.getenv("VKR_SIGMAS", "0.05,0.20"), ",")[[1L]])
comparison_threshold_multipliers <- as.numeric(strsplit(Sys.getenv("VKR_THRESHOLDS", "0,1,2"), ",")[[1L]])
comparison_seed <- as.integer(Sys.getenv("VKR_SEED", "2026"))
comparison_line_method <- Sys.getenv("VKR_LINE_METHOD", "default")
comparison_default_intensity <- 0.8
comparison_table_scope <- Sys.getenv("VKR_TABLE_SCOPE", "report")

one_line_rho_step <- as.numeric(Sys.getenv("VKR_ONE_RHO_STEP", "0.05"))
one_line_theta_step <- as.numeric(Sys.getenv("VKR_ONE_THETA_STEP", "0.005"))
two_line_rho_step <- as.numeric(Sys.getenv("VKR_TWO_RHO_STEP", "0.02"))
two_line_theta_step <- as.numeric(Sys.getenv("VKR_TWO_THETA_STEP", "0.002"))

one_line_configs <- list(
  diag_pos_full = data.frame(a = 1, b = 0, intensity = comparison_default_intensity),
  diag_neg_full = data.frame(a = -1, b = 101, intensity = comparison_default_intensity),
  steep_full = data.frame(a = 2, b = -1, intensity = comparison_default_intensity),
  steep_shifted = data.frame(a = 2, b = -40, intensity = comparison_default_intensity)
)

two_line_configs <- list(
  two_lines = data.frame(
    a = c(2, -1),
    b = c(-1, 101),
    intensity = c(0.8, 0.5)
  )
)

method_code_order <- c(
  "row.row",
  "col.row",
  "max.row",
  "median",
  "wiener",
  "row.row.esprit"
)

method_display <- c(
  row.row = "\\texttt{ROW.ROW}",
  col.row = "\\texttt{COL.ROW}",
  max.row = "\\texttt{MAX.ROW}",
  median = "Median",
  wiener = "Wiener",
  row.row.esprit = "\\texttt{ROW.ROW}+ESPRIT"
)

metric_order <- c("mean_rho_mse", "mean_theta_mse", "mean_active_pixels")
metric_labels <- c(
  mean_rho_mse = "MSE по $\\rho$",
  mean_theta_mse = "MSE по $\\theta$",
  mean_active_pixels = "Активные пиксели"
)

resolve_positive_count <- function(value, default) {
  value <- suppressWarnings(as.integer(value)[1L])
  if (!is.finite(value) || is.na(value)) {
    return(max(1L, as.integer(default)[1L]))
  }
  max(1L, value)
}

method_settings <- function(lines, n_row, n_col, threshold_multiplier) {
  lines <- normalize_config_lines(lines)
  row_rank <- resolve_method_rank("row.row", lines, n_row, n_col)
  col_rank <- resolve_method_rank("col.row", lines, n_row, n_col)
  max_k <- nrow(lines)

  data.frame(
    method_code = method_code_order,
    method = c(
      "CSSA ROW.ROW",
      "CSSA COL.ROW",
      "MAX.ROW",
      "Median",
      "Wiener",
      "CSSA ROW.ROW + ESPRIT"
    ),
    threshold_multiplier = threshold_multiplier,
    use_noise_cut = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
    m_coef = c(threshold_multiplier, threshold_multiplier, NA, threshold_multiplier, threshold_multiplier, threshold_multiplier),
    num_components = c(row_rank$value, col_rank$value, max_k, 0L, 0L, row_rank$value),
    cssa_rank_source = c(row_rank$source, col_rank$source, "k", "filter", "filter", row_rank$source),
    cssa_auto_rank = c(row_rank$auto, col_rank$auto, max_k, 0L, 0L, row_rank$auto),
    frequency_grid_multiplier = c(1L, esprit_grid_multiplier_from_lines(lines), 1L, 1L, 1L, 1L),
    k = c(1L, 1L, max_k, NA_integer_, NA_integer_, 1L),
    stringsAsFactors = FALSE
  )
}

preprocess_binary <- function(noisy, spec, sigma_noise) {
  switch(
    spec$method_code,
    row.row = cssa_method_binary(
      noisy,
      method = "row.row",
      noisy = noisy,
      use_noise_cut = spec$use_noise_cut,
      m_coef = spec$m_coef,
      num_components = spec$num_components
    ),
    col.row = cssa_method_binary(
      noisy,
      method = "col.row",
      noisy = noisy,
      use_noise_cut = spec$use_noise_cut,
      m_coef = spec$m_coef,
      num_components = spec$num_components
    ),
    max.row = row_max_binary(noisy, k = spec$k),
    median = threshold_binary(
      denoised = median_denoise(noisy, n = 3L),
      noisy = noisy,
      use_noise_cut = spec$use_noise_cut,
      m_coef = spec$m_coef
    ),
    wiener = threshold_binary(
      denoised = wiener_denoise(noisy, ksize = 5L),
      noisy = noisy,
      use_noise_cut = spec$use_noise_cut,
      m_coef = spec$m_coef
    ),
    row.row.esprit = threshold_binary(
      denoised = cssa_esprit_denoise(
        noisy,
        method = "row.row",
        num_components = spec$num_components,
        frequency_grid_multiplier = spec$frequency_grid_multiplier,
        clip = TRUE
      ),
      noisy = noisy,
      use_noise_cut = spec$use_noise_cut,
      m_coef = spec$m_coef
    ),
    stop("Unknown method: ", spec$method_code)
  )
}

run_single_method <- function(noisy,
                              lines,
                              spec,
                              sigma,
                              rho_step,
                              theta_step) {
  processed <- preprocess_binary(noisy, spec, sigma_noise = sigma)
  ht <- make_accumulator(
    processed$binary,
    detector = function(m) m,
    rho_step = rho_step,
    theta_step = theta_step
  )

  true_lines <- config_true_lines(lines)
  pred_lines <- find_k_max(
    ht$accumulator,
    k = nrow(true_lines),
    qrho = ht$rho,
    qtheta = ht$theta,
    suppress = nrow(true_lines) > 1L
  )
  err <- parameter_line_errors_by_line(true_lines, pred_lines)
  labels <- line_config_labels(lines)

  data.frame(
    config_id = labels[err$line_index],
    line_index = err$line_index,
    pred_index = err$pred_index,
    method_code = spec$method_code,
    method = spec$method,
    sigma = sigma,
    threshold_multiplier = spec$threshold_multiplier,
    num_components = spec$num_components,
    cssa_rank_source = spec$cssa_rank_source,
    cssa_auto_rank = spec$cssa_auto_rank,
    frequency_grid_multiplier = spec$frequency_grid_multiplier,
    k = spec$k,
    threshold_value = processed$threshold_value,
    noise_sd = processed$noise_sd,
    true_rho = err$true_rho,
    true_theta = err$true_theta,
    pred_rho = err$pred_rho,
    pred_theta = err$pred_theta,
    rho_mse = err$rho_mse,
    theta_mse = err$theta_mse,
    active_pixels = sum(processed$binary > 0),
    row.names = NULL
  )
}

run_experiment_family <- function(configs,
                                  family,
                                  rho_step,
                                  theta_step,
                                  n_rep = comparison_n_rep) {
  out <- list()
  counter <- 1L

  for (sigma in comparison_sigmas) {
    set.seed(comparison_seed + as.integer(round(1000 * sigma)) + if (family == "two") 100000L else 0L)
    for (config_name in names(configs)) {
      lines <- normalize_config_lines(configs[[config_name]])
      clean <- make_line_image(
        n_row = comparison_n_row,
        n_col = comparison_n_col,
        lines = lines,
        line_method = comparison_line_method
      )

      for (rep_i in seq_len(n_rep)) {
        noisy <- add.noise(clean, sigma = sigma)
        for (threshold_multiplier in comparison_threshold_multipliers) {
          settings <- method_settings(lines, comparison_n_row, comparison_n_col, threshold_multiplier)
          for (i in seq_len(nrow(settings))) {
            cur <- run_single_method(
              noisy = noisy,
              lines = lines,
              spec = settings[i, ],
              sigma = sigma,
              rho_step = rho_step,
              theta_step = theta_step
            )
            cur$family <- family
            cur$parent_config_id <- config_name
            cur$rep <- rep_i
            cur$rho_step <- rho_step
            cur$theta_step <- theta_step
            out[[counter]] <- cur
            counter <- counter + 1L
          }
        }
      }
    }
  }

  do.call(rbind, out)
}

mean_or_na <- function(v) {
  if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
}

summarise_samples <- function(samples) {
  out <- aggregate(
    cbind(rho_mse, theta_mse, active_pixels, threshold_value, noise_sd) ~
      family + parent_config_id + config_id + line_index + sigma +
      threshold_multiplier + method_code + method + num_components +
      cssa_rank_source + cssa_auto_rank + frequency_grid_multiplier +
      rho_step + theta_step,
    data = samples,
    FUN = mean_or_na,
    na.action = na.pass
  )
  names(out)[names(out) == "rho_mse"] <- "mean_rho_mse"
  names(out)[names(out) == "theta_mse"] <- "mean_theta_mse"
  names(out)[names(out) == "active_pixels"] <- "mean_active_pixels"
  names(out)[names(out) == "threshold_value"] <- "mean_threshold"
  names(out)[names(out) == "noise_sd"] <- "mean_noise_sd"
  out$n_rep <- length(unique(samples$rep))
  out[order(
    out$family,
    out$sigma,
    out$threshold_multiplier,
    out$parent_config_id,
    out$line_index,
    match(out$method_code, method_code_order)
  ), ]
}

sigma_suffix <- function(sigma) {
  sprintf("%03d", as.integer(round(100 * as.numeric(sigma))))
}

threshold_suffix <- function(threshold_multiplier) {
  if (threshold_multiplier <= 0) "nothr" else paste0("thr", as.integer(threshold_multiplier))
}

threshold_caption <- function(threshold_multiplier) {
  if (threshold_multiplier <= 0) {
    "без порога"
  } else {
    paste0("порог $", threshold_multiplier, "\\hat\\sigma$")
  }
}

make_latex_table <- function(summary_df, family_value, sigma_value, threshold_value, file_id, caption, label) {
  cur <- subset(
    summary_df,
    family == family_value &
      abs(sigma - sigma_value) < 1e-12 &
      abs(threshold_multiplier - threshold_value) < 1e-12
  )

  config_order <- unique(cur$config_id[order(cur$parent_config_id, cur$line_index)])
  rows <- list()
  for (cfg in config_order) {
    for (metric in metric_order) {
      vals <- vapply(method_code_order, function(code) {
        x <- cur[cur$config_id == cfg & cur$method_code == code, metric, drop = TRUE]
        if (length(x) == 0L) return("")
        format_metric_value(x[1L], metric)
      }, character(1L))
      rows[[length(rows) + 1L]] <- c(
        if (metric == metric_order[1L]) paste0("$", format_config_label(cfg), "$") else "",
        metric_labels[[metric]],
        vals
      )
    }
  }

  df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(df) <- c(
    "Конфигурация",
    "Метрика",
    unname(method_display[method_code_order])
  )

  align <- paste0("|l|l|", paste(rep("c", length(method_code_order)), collapse = "|"), "|")
  header <- paste(names(df), collapse = " & ")
  body <- unlist(lapply(seq_len(nrow(df)), function(i) {
    c(paste(df[i, ], collapse = " & "), "\\\\", "\\hline")
  }))

  lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\resizebox{\\textwidth}{!}{%",
    sprintf("\\begin{tabular}{%s}", align),
    "\\hline",
    header,
    "\\\\",
    "\\hline",
    body,
    "\\end{tabular}%",
    "}",
    "\\end{table}"
  )
  write_lines_utf8(lines, file.path(chapter3_table_dir, paste0(file_id, ".tex")))
}

write_parameter_tables <- function() {
  params_03 <- data.frame(
    "$N_{\\text{row}}$" = comparison_n_row,
    "$N_{\\text{col}}$" = comparison_n_col,
    "$n$" = comparison_n_rep,
    "$\\sigma$" = "$\\{0.05,\\,0.20\\}$",
    "пороги" = "без порога, $1\\hat\\sigma$, $2\\hat\\sigma$",
    "$h_\\rho$" = format_plain_num(one_line_rho_step, 4L),
    "$h_\\theta$" = format_plain_num(one_line_theta_step, 4L),
    check.names = FALSE
  )
  params_04 <- data.frame(
    "$N_{\\text{row}}$" = comparison_n_row,
    "$N_{\\text{col}}$" = comparison_n_col,
    "$n$" = comparison_n_rep,
    "$\\sigma$" = "$\\{0.05,\\,0.20\\}$",
    "пороги" = "без порога, $1\\hat\\sigma$, $2\\hat\\sigma$",
    "$h_\\rho$" = format_plain_num(two_line_rho_step, 4L),
    "$h_\\theta$" = format_plain_num(two_line_theta_step, 4L),
    check.names = FALSE
  )
  write_lines_utf8(latex_table_lines(params_03, resize = NULL), file.path(chapter3_table_dir, "ch3_exp_params_03.tex"))
  write_lines_utf8(latex_table_lines(params_04, resize = NULL), file.path(chapter3_table_dir, "ch3_exp_params_04.tex"))
}

write_all_result_tables <- function(summary_df) {
  if (identical(comparison_table_scope, "all")) {
    specs <- expand.grid(
      family = c("one", "two"),
      sigma = comparison_sigmas,
      threshold = comparison_threshold_multipliers,
      stringsAsFactors = FALSE
    )
  } else {
    specs <- data.frame(
      family = c("one", "one", "two", "two", "one", "one", "one", "two", "two"),
      sigma = c(0.05, 0.20, 0.05, 0.20, 0.05, 0.20, 0.20, 0.05, 0.20),
      threshold = c(0, 0, 0, 0, 1, 1, 2, 1, 1),
      stringsAsFactors = FALSE
    )
    specs <- specs[
      specs$sigma %in% comparison_sigmas &
        specs$threshold %in% comparison_threshold_multipliers,
      ,
      drop = FALSE
    ]
  }

  for (i in seq_len(nrow(specs))) {
    family <- specs$family[i]
    sigma <- specs$sigma[i]
    threshold <- specs$threshold[i]
    prefix <- if (family == "one") "ch3_one" else "ch3_two"
    title <- if (family == "one") "Одна прямая" else "Две прямые"
    file_id <- paste(prefix, sigma_suffix(sigma), threshold_suffix(threshold), sep = "_")
    caption <- paste0(
      title,
      ": $\\sigma=",
      formatC(sigma, format = "f", digits = 2L),
      "$, ",
      threshold_caption(threshold)
    )
    make_latex_table(
      summary_df = summary_df,
      family_value = family,
      sigma_value = sigma,
      threshold_value = threshold,
      file_id = file_id,
      caption = caption,
      label = paste0("tab:", file_id)
    )
  }
}

write_metadata <- function(summary_df) {
  one_settings <- do.call(rbind, lapply(names(one_line_configs), function(name) {
    lines <- normalize_config_lines(one_line_configs[[name]])
    cbind(
      family = "one",
      parent_config_id = name,
      method_settings(lines, comparison_n_row, comparison_n_col, threshold_multiplier = 1)
    )
  }))
  two_settings <- do.call(rbind, lapply(names(two_line_configs), function(name) {
    lines <- normalize_config_lines(two_line_configs[[name]])
    cbind(
      family = "two",
      parent_config_id = name,
      method_settings(lines, comparison_n_row, comparison_n_col, threshold_multiplier = 1)
    )
  }))
  settings <- rbind(one_settings, two_settings)

  write_csv_utf8(settings, file.path(chapter3_data_dir, "method_settings.csv"))
  write_csv_utf8(
    data.frame(
      parameter = c(
        "n_row",
        "n_col",
        "n_rep",
        "sigmas",
        "threshold_multipliers",
        "one_line_rho_step",
        "one_line_theta_step",
        "two_line_rho_step",
        "two_line_theta_step",
        "line_method",
        "seed",
        "table_scope"
      ),
      value = c(
        comparison_n_row,
        comparison_n_col,
        comparison_n_rep,
        paste(comparison_sigmas, collapse = ","),
        paste(comparison_threshold_multipliers, collapse = ","),
        one_line_rho_step,
        one_line_theta_step,
        two_line_rho_step,
        two_line_theta_step,
        comparison_line_method,
        comparison_seed,
        comparison_table_scope
      )
    ),
    file.path(chapter3_data_dir, "experiment_parameters.csv")
  )

  invisible(settings)
}

cat("Running one-line comparisons...\n")
one_samples <- run_experiment_family(one_line_configs, "one", one_line_rho_step, one_line_theta_step)
cat("Running two-line comparisons...\n")
two_samples <- run_experiment_family(two_line_configs, "two", two_line_rho_step, two_line_theta_step)

samples <- rbind(one_samples, two_samples)
summary <- summarise_samples(samples)

write_csv_utf8(samples, file.path(chapter3_data_dir, "chapter3_raw_samples.csv"))
write_csv_utf8(summary, file.path(chapter3_data_dir, "chapter3_summary.csv"))
write_metadata(summary)
write_parameter_tables()
write_all_result_tables(summary)

cat("Done. LaTeX tables are in: ", chapter3_table_dir, "\n", sep = "")
cat("Raw samples and summary CSV are in: ", chapter3_data_dir, "\n", sep = "")
