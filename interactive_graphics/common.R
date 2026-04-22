find_project_root <- function(start_dir = getwd()) {
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
  stop("Не удалось найти корень проекта.")
}

common_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
common_dir <- if (!is.null(common_file)) {
  dirname(normalizePath(common_file, winslash = "/", mustWork = FALSE))
} else {
  getwd()
}
project_root <- find_project_root(common_dir)

required_pkgs <- c(
  "shiny",
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
  stop("Отсутствуют пакеты: ", paste(missing_pkgs, collapse = ", "))
}

source(file.path(project_root, "ssa-based methods", "cssa-transform.r"))
hough_dir <- file.path(project_root, "hough transform")
rcpp_cache_dir <- file.path(project_root, ".rcpp-cache")
dir.create(rcpp_cache_dir, showWarnings = FALSE, recursive = TRUE)
Rcpp::sourceCpp(file.path(hough_dir, "hough_fast.cpp"), cacheDir = rcpp_cache_dir)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

bm3d_available <- reticulate::py_module_available("bm3d")
bm3d_lib <- if (bm3d_available) reticulate::import("bm3d") else NULL

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
    stop("Неизвестный CSSA-метод: ", method)
  )
  clip01(cleaned)
}

bm3d_denoise <- function(m, sigma_noise) {
  if (!bm3d_available) {
    stop("Python-модуль 'bm3d' недоступен.")
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
  out <- local_mean + gain * (x - local_mean)
  clip01(out)
}

available_detectors <- function(num_of_lines, include_quantile = FALSE) {
  detectors <- list(
    cssa_row_row = function(m, sigma_noise) {
      cssa_denoise(m, num_of_lines = num_of_lines, method = "row.row")
    },
    median = function(m, sigma_noise) {
      median_denoise(m, n = 3L)
    },
    wiener = function(m, sigma_noise) {
      wiener_denoise(m, ksize = 5L)
    }
  )

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
  }

  if (include_quantile) {
    detectors$quantile_09 <- function(m, sigma_noise) {
      x <- clip01(m)
      thr <- as.numeric(stats::quantile(x, probs = 0.9, na.rm = TRUE))
      clip01(ifelse(x >= thr, x, 0))
    }
  }

  detectors
}

detector_labels <- function(detectors) {
  labels <- c(
    cssa_row_row = "CSSA row.row",
    bm3d = "BM3D",
    median = "Median",
    wiener = "Wiener",
    quantile_09 = "Quantile 0.9"
  )
  labels[names(detectors)]
}

convert_ab_to_rho_theta <- function(a, b) {
  out <- convert_ab_to_rt_cpp(a, b)
  unname(out[1, ])
}

convert_rho_theta_to_ab <- function(rho, theta) {
  out <- convert_rt_to_ab_cpp(rho, theta)
  unname(out[1, ])
}

find_k_max <- function(acc, k, qrho, qtheta, suppress = FALSE, window = 6L) {
  find_k_max_cpp(
    acc = as.matrix(acc),
    k = as.integer(k),
    qrho = as.numeric(qrho),
    qtheta = as.numeric(qtheta),
    suppress = isTRUE(suppress),
    window = as.integer(window)
  )
}

compute_err <- function(true, pred) {
  return(compute_err_cpp(as.matrix(true), as.matrix(pred)))

  true <- as.matrix(true)
  pred <- as.matrix(pred)

  if (ncol(true) != 2 || ncol(pred) != 2) {
    stop("true и pred должны иметь по два столбца: rho и theta")
  }
  if (nrow(true) != nrow(pred)) {
    stop("true и pred должны содержать одинаковое число прямых")
  }

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

  n_lines <- nrow(true)
  if (n_lines == 1L) {
    out <- pair_err(true[1, 1], true[1, 2], pred[1, 1], pred[1, 2])
    return(c(dr = unname(out[1]), dtheta = unname(out[2])))
  }

  cost <- matrix(0, nrow = n_lines, ncol = n_lines)
  for (i in seq_len(n_lines)) {
    for (j in seq_len(n_lines)) {
      cost[i, j] <- pair_err(
        true[i, 1], true[i, 2],
        pred[j, 1], pred[j, 2]
      )[3]
    }
  }

  all_perms <- function(v) {
    if (length(v) == 1L) {
      return(matrix(v, nrow = 1L))
    }
    out <- lapply(seq_along(v), function(i) {
      rest <- all_perms(v[-i])
      cbind(v[i], rest)
    })
    do.call(rbind, out)
  }

  perm_mat <- all_perms(seq_len(n_lines))
  best_perm <- NULL
  best_cost <- Inf
  for (row_i in seq_len(nrow(perm_mat))) {
    cur_perm <- perm_mat[row_i, ]
    cur_cost <- sum(cost[cbind(seq_len(n_lines), cur_perm)])
    if (cur_cost < best_cost) {
      best_cost <- cur_cost
      best_perm <- cur_perm
    }
  }

  dr_total <- 0
  dtheta_total <- 0
  for (i in seq_len(n_lines)) {
    out <- pair_err(
      true[i, 1], true[i, 2],
      pred[best_perm[i], 1], pred[best_perm[i], 2]
    )
    dr_total <- dr_total + out[1]
    dtheta_total <- dtheta_total + out[2]
  }

  c(dr = unname(dr_total), dtheta = unname(dtheta_total))
}

parse_line_text <- function(text, expected_n = NULL) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0L) {
    stop("Список прямых пуст.")
  }

  parsed <- lapply(lines, function(line) {
    parts <- unlist(strsplit(line, "[,;[:space:]]+", perl = TRUE))
    parts <- parts[nzchar(parts)]
    nums <- as.numeric(parts)
    if (any(!is.finite(nums))) {
      stop("Не удалось распознать строку: ", line)
    }
    if (!(length(nums) %in% c(2L, 3L))) {
      stop("Каждая строка должна содержать a,b или a,b,intensity. Ошибка в строке: ", line)
    }
    data.frame(
      a = nums[1],
      b = nums[2],
      intensity = if (length(nums) == 3L) nums[3] else 1,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, parsed)
  if (!is.null(expected_n) && nrow(out) != expected_n) {
    stop("Ожидалось ", expected_n, " прямых, а получено ", nrow(out), ".")
  }
  out
}

standard_single_configs <- function(n_row, n_col) {
  defs <- data.frame(
    config_id = c("diag_pos_full", "diag_neg_full", "steep_full", "steep_shifted"),
    description = c(
      "Positive full diagonal",
      "Negative full diagonal",
      "Steep line crossing most of the image",
      "Steep shifted line with shorter visible segment"
    ),
    a = c(1, -1, 2, 2),
    b = c(0, 101, -1, -40),
    stringsAsFactors = FALSE
  )

  lapply(seq_len(nrow(defs)), function(i) {
    list(
      config_id = defs$config_id[i],
      description = defs$description[i],
      n_row = n_row,
      n_col = n_col,
      lines = defs[i, c("a", "b")],
      intensity = 1
    )
  }) |>
    lapply(function(cfg) {
      cfg$lines$intensity <- 1
      cfg
    })
}

make_custom_config <- function(lines_df, n_row, n_col, config_id = "custom", description = "Custom configuration") {
  list(
    config_id = config_id,
    description = description,
    n_row = n_row,
    n_col = n_col,
    lines = lines_df
  )
}

build_config_list <- function(num_lines, n_row, n_col, use_standard_single, line_text) {
  num_lines <- as.integer(num_lines)
  n_row <- as.integer(n_row)
  n_col <- as.integer(n_col)

  if (num_lines == 1L && isTRUE(use_standard_single)) {
    return(standard_single_configs(n_row = n_row, n_col = n_col))
  }

  lines_df <- parse_line_text(line_text, expected_n = num_lines)
  list(make_custom_config(lines_df, n_row = n_row, n_col = n_col))
}

config_to_matrix <- function(cfg) {
  out <- matrix(0, nrow = cfg$n_row, ncol = cfg$n_col)
  for (i in seq_len(nrow(cfg$lines))) {
    out <- add.line(
      out,
      a = cfg$lines$a[i],
      b = cfg$lines$b[i],
      intensity = cfg$lines$intensity[i]
    )
  }
  out
}

config_true_params <- function(cfg) {
  out <- t(vapply(seq_len(nrow(cfg$lines)), function(i) {
    convert_ab_to_rho_theta(cfg$lines$a[i], cfg$lines$b[i])
  }, numeric(2)))
  colnames(out) <- c("rho", "theta")
  out
}

draw_config_preview <- function(config_list) {
  mats <- lapply(config_list, config_to_matrix)
  labels <- vapply(config_list, function(cfg) {
    if (nrow(cfg$lines) == 1L) {
      sprintf("%s\n(a = %.2f, b = %.2f)", cfg$config_id, cfg$lines$a[1], cfg$lines$b[1])
    } else {
      sprintf("%s\n%d lines", cfg$config_id, nrow(cfg$lines))
    }
  }, character(1))

  nplots <- if (length(mats) <= 2L) length(mats) else 2L
  plot.matrix(mats, from.0.to.1 = TRUE, labels = labels, nplots = nplots)
}

threshold_processed_matrix <- function(
  processed_matrix,
  noisy_matrix,
  num_of_lines,
  threshold_mode = "auto_rowrow_sd",
  threshold_value = 0.1,
  threshold_multiplier = 1,
  residual_sd = NULL
) {
  processed_matrix <- clip01(processed_matrix)

  if (threshold_mode == "none") {
    return(list(
      processed = processed_matrix,
      threshold_value = 0,
      residual_sd = NA_real_
    ))
  }

  if (threshold_mode == "manual") {
    thr <- as.numeric(threshold_value)
    out <- processed_matrix
    out[out < thr] <- 0
    return(list(
      processed = out,
      threshold_value = thr,
      residual_sd = NA_real_
    ))
  }

  if (is.null(residual_sd)) {
    rowrow_signal <- cssa_denoise(noisy_matrix, num_of_lines = num_of_lines, method = "row.row")
    residual_sd <- stats::sd(as.vector(noisy_matrix - rowrow_signal), na.rm = TRUE)
  }
  thr <- as.numeric(threshold_multiplier) * residual_sd
  out <- processed_matrix
  out[out < thr] <- 0
  list(
    processed = out,
    threshold_value = thr,
    residual_sd = residual_sd
  )
}

make_accumulator_generic <- function(processed_matrix, rho_step = 1, theta_step = 0.01, weighted = FALSE) {
  bound_matrix <- as.matrix(processed_matrix)
  bound_matrix[!is.finite(bound_matrix)] <- 0
  if (is.complex(bound_matrix)) {
    bound_matrix <- Re(bound_matrix)
  }

  return(make_accumulator_processed(
    bound_matrix,
    rho_step = as.numeric(rho_step),
    theta_step = as.numeric(theta_step),
    weighted = isTRUE(weighted)
  ))

  bound_matrix <- as.matrix(processed_matrix)
  bound_matrix[!is.finite(bound_matrix)] <- 0
  if (is.complex(bound_matrix)) {
    bound_matrix <- Re(bound_matrix)
  }

  n_row <- nrow(bound_matrix)
  n_col <- ncol(bound_matrix)
  points <- which(bound_matrix > 0, arr.ind = TRUE)

  n_theta <- floor(pi / theta_step) + 1L
  theta <- seq.int(0L, n_theta - 1L) * theta_step

  rho_max <- ceiling(sqrt(n_row * n_row + n_col * n_col))
  n_rho <- floor((2 * rho_max) / rho_step) + 1L
  rho <- -rho_max + seq.int(0L, n_rho - 1L) * rho_step

  accumulator <- matrix(0, nrow = n_rho, ncol = n_theta)
  if (nrow(points) == 0L) {
    return(list(
      accumulator = accumulator,
      rho = rho,
      theta = theta,
      active_pixels = 0L,
      active_weight = 0
    ))
  }

  x <- points[, "col"]
  y <- points[, "row"]
  weights <- if (weighted) bound_matrix[points] else rep(1, nrow(points))

  for (j in seq_along(theta)) {
    rho_values <- x * cos(theta[j]) + y * sin(theta[j])
    rho_idx <- round((rho_values - rho[1]) / rho_step) + 1L
    keep <- rho_idx >= 1L & rho_idx <= n_rho
    cur_sum <- tapply(weights[keep], rho_idx[keep], sum)
    accumulator[as.integer(names(cur_sum)), j] <- as.numeric(cur_sum)
  }

  list(
    accumulator = accumulator,
    rho = rho,
    theta = theta,
    active_pixels = length(weights),
    active_weight = sum(weights)
  )
}

run_ht_experiment_app <- function(
  config_list,
  method_names,
  sigma_noise,
  num_of_lines,
  rho_step_ht,
  theta_step_ht,
  n_rep,
  ht_type,
  threshold_mode,
  threshold_value,
  threshold_multiplier,
  progress = NULL
) {
  detectors <- available_detectors(num_of_lines = num_of_lines)
  method_names <- intersect(method_names, names(detectors))
  if (length(method_names) == 0L) {
    stop("Не выбрано ни одного доступного метода.")
  }

  results_list <- vector("list", length(config_list) * n_rep * length(method_names))
  idx_out <- 1L
  weighted <- identical(ht_type, "weighted")
  total_steps <- length(config_list) * n_rep * length(method_names)
  step_i <- 0L

  for (cfg in config_list) {
    base_matrix <- config_to_matrix(cfg)
    true_line <- config_true_params(cfg)

    for (rep_i in seq_len(n_rep)) {
      noisy_matrix <- add.noise(base_matrix, sigma = sigma_noise)
      residual_sd <- NULL
      if (identical(threshold_mode, "auto_rowrow_sd")) {
        rowrow_signal <- cssa_denoise(noisy_matrix, num_of_lines = num_of_lines, method = "row.row")
        residual_sd <- stats::sd(as.vector(noisy_matrix - rowrow_signal), na.rm = TRUE)
      }

      for (method_name in method_names) {
        processed_matrix <- detectors[[method_name]](noisy_matrix, sigma_noise = sigma_noise)
        threshold_info <- threshold_processed_matrix(
          processed_matrix = processed_matrix,
          noisy_matrix = noisy_matrix,
          num_of_lines = num_of_lines,
          threshold_mode = threshold_mode,
          threshold_value = threshold_value,
          threshold_multiplier = threshold_multiplier,
          residual_sd = residual_sd
        )

        ht_result <- make_accumulator_generic(
          threshold_info$processed,
          rho_step = rho_step_ht,
          theta_step = theta_step_ht,
          weighted = weighted
        )

        pred_line <- find_k_max(
          ht_result$accumulator,
          k = num_of_lines,
          qrho = ht_result$rho,
          qtheta = ht_result$theta,
          suppress = num_of_lines > 1L,
          window = 6L
        )

        err <- compute_err(true_line, pred_line)

        results_list[[idx_out]] <- data.frame(
          config_id = cfg$config_id,
          description = cfg$description,
          n_col = cfg$n_col,
          n_row = cfg$n_row,
          sigma = sigma_noise,
          rep = rep_i,
          method = method_name,
          true_rho = paste(round(true_line[, 1], 6), collapse = "; "),
          true_theta = paste(round(true_line[, 2], 6), collapse = "; "),
          pred_rho = paste(round(pred_line[, 1], 6), collapse = "; "),
          pred_theta = paste(round(pred_line[, 2], 6), collapse = "; "),
          active_pixels = ht_result$active_pixels,
          active_weight = ht_result$active_weight,
          threshold_value = threshold_info$threshold_value,
          dr = as.numeric(err["dr"]),
          dtheta = as.numeric(err["dtheta"]),
          stringsAsFactors = FALSE
        )
        idx_out <- idx_out + 1L

        step_i <- step_i + 1L
        if (!is.null(progress)) {
          progress(step_i / total_steps, detail = sprintf("%s / rep %d / %s", cfg$config_id, rep_i, method_name))
        }
      }
    }
  }

  dplyr::bind_rows(results_list)
}

summarise_ht_errors_app <- function(ht_error_samples) {
  ht_error_samples |>
    dplyr::group_by(config_id, description, method) |>
    dplyr::summarise(
      mean_dr = mean(dr),
      median_dr = median(dr),
      mean_dtheta = mean(dtheta),
      median_dtheta = median(dtheta),
      mean_active_pixels = mean(active_pixels),
      mean_active_weight = mean(active_weight),
      mean_threshold = mean(threshold_value),
      .groups = "drop"
    ) |>
    dplyr::arrange(config_id, mean_dr, mean_dtheta)
}

matrix_to_df <- function(tab, row_name = "config_id") {
  out <- as.data.frame.matrix(tab, stringsAsFactors = FALSE)
  out[[row_name]] <- rownames(out)
  out <- out[, c(row_name, setdiff(names(out), row_name)), drop = FALSE]
  rownames(out) <- NULL
  out
}

ideal_discretization_error <- function(true_line, rho_step_ht, theta_step_ht) {
  return(ideal_discretization_error_cpp(
    as.matrix(true_line),
    rho_step_ht = as.numeric(rho_step_ht),
    theta_step_ht = as.numeric(theta_step_ht)
  ))

  pred <- cbind(
    rho = round(true_line[, 1] / rho_step_ht) * rho_step_ht,
    theta = round(true_line[, 2] / theta_step_ht) * theta_step_ht
  )
  pred[, 2] <- pmin(pmax(pred[, 2], 0), pi)

  err <- compute_err(true_line, pred)
  baseline_dr <- max(as.numeric(err["dr"]), (rho_step_ht / 2)^2)
  baseline_dtheta <- max(as.numeric(err["dtheta"]), (theta_step_ht / 2)^2)

  list(
    pred = pred,
    dr = as.numeric(err["dr"]),
    dtheta = as.numeric(err["dtheta"]),
    baseline_dr = baseline_dr,
    baseline_dtheta = baseline_dtheta
  )
}

find_big_error_cases <- function(
  cfg,
  method_name,
  sigma_noise,
  num_of_lines,
  num_maxima = 1L,
  rho_step_ht,
  theta_step_ht,
  ht_type,
  threshold_mode,
  threshold_value,
  threshold_multiplier,
  n_search,
  factor_threshold,
  max_cases,
  progress = NULL
) {
  detectors <- available_detectors(num_of_lines = num_of_lines)
  if (!method_name %in% names(detectors)) {
    stop("Метод ", method_name, " недоступен.")
  }

  num_maxima <- max(as.integer(num_maxima), as.integer(num_of_lines), 1L)
  base_matrix <- config_to_matrix(cfg)
  true_line <- config_true_params(cfg)
  ideal <- ideal_discretization_error(true_line, rho_step_ht, theta_step_ht)
  weighted <- identical(ht_type, "weighted")
  cases <- list()

  for (rep_i in seq_len(n_search)) {
    if (length(cases) >= max_cases) {
      break
    }

    noisy_matrix <- add.noise(base_matrix, sigma = sigma_noise)
    residual_sd <- NULL
    if (identical(threshold_mode, "auto_rowrow_sd")) {
      rowrow_signal <- cssa_denoise(noisy_matrix, num_of_lines = num_of_lines, method = "row.row")
      residual_sd <- stats::sd(as.vector(noisy_matrix - rowrow_signal), na.rm = TRUE)
    }
    processed_matrix <- detectors[[method_name]](noisy_matrix, sigma_noise = sigma_noise)
    threshold_info <- threshold_processed_matrix(
      processed_matrix = processed_matrix,
      noisy_matrix = noisy_matrix,
      num_of_lines = num_of_lines,
      threshold_mode = threshold_mode,
      threshold_value = threshold_value,
      threshold_multiplier = threshold_multiplier,
      residual_sd = residual_sd
    )

    ht_result <- make_accumulator_generic(
      threshold_info$processed,
      rho_step = rho_step_ht,
      theta_step = theta_step_ht,
      weighted = weighted
    )

    pred_line_all <- find_k_max(
      ht_result$accumulator,
      k = num_maxima,
      qrho = ht_result$rho,
      qtheta = ht_result$theta,
      suppress = num_maxima > 1L,
      window = 6L
    )

    pred_line <- pred_line_all[seq_len(num_of_lines), , drop = FALSE]
    err <- compute_err(true_line, pred_line)
    dr_ratio <- as.numeric(err["dr"]) / max(as.numeric(ideal$dr), 1e-12)
    dtheta_ratio <- as.numeric(err["dtheta"]) / max(as.numeric(ideal$dtheta), 1e-12)

    if (dr_ratio >= factor_threshold || dtheta_ratio >= factor_threshold) {
      cases[[length(cases) + 1L]] <- list(
        case_id = sprintf("rep_%03d", rep_i),
        rep = rep_i,
        dr = as.numeric(err["dr"]),
        dtheta = as.numeric(err["dtheta"]),
        dr_ratio = dr_ratio,
        dtheta_ratio = dtheta_ratio,
        threshold_value = threshold_info$threshold_value,
        noisy_matrix = noisy_matrix,
        processed_matrix = threshold_info$processed,
        ht_result = ht_result,
        pred_line = pred_line,
        pred_line_all = pred_line_all,
        true_line = true_line
      )
    }

    if (!is.null(progress)) {
      progress(rep_i / n_search, detail = sprintf("rep %d / %d", rep_i, n_search))
    }
  }

  list(ideal = ideal, cases = cases)
}

draw_case_matrix_triplet <- function(case) {
  rgb.palette <- grDevices::colorRampPalette(c("white", "black"), space = "rgb")
  make_panel <- function(mat, title) {
    mat <- as.matrix(mat)
    if (is.complex(mat)) {
      mat <- Re(mat)
    }
    mat[mat < 0] <- 0
    mat[mat > 1] <- 1

    lattice::levelplot(
      t(mat),
      xlab = "x",
      ylab = "y",
      main = title,
      col.regions = rgb.palette,
      at = seq(0, 1, 0.01),
      colorkey = FALSE,
      col = "transparent",
      border = NA,
      cuts = 255,
      scales = list(
        draw = TRUE,
        x = list(at = pretty(seq_len(ncol(mat)), n = 5)),
        y = list(at = pretty(seq_len(nrow(mat)), n = 5))
      )
    )
  }

  return(gridExtra::grid.arrange(
    make_panel(clip01(case$noisy_matrix), "Шумная матрица"),
    make_panel(case$processed_matrix, "После обработки"),
    make_panel(1 * (case$processed_matrix > 0), "Бинаризация > 0"),
    ncol = 3
  ))

  plot.matrix(
    list(
      "Шумная матрица" = clip01(case$noisy_matrix),
      "После обработки" = case$processed_matrix,
      "Бинаризация > 0" = 1 * (case$processed_matrix > 0)
    ),
    from.0.to.1 = TRUE,
    labels = c("Шумная матрица", "После обработки", "Бинаризация > 0"),
    nplots = 3
  )
}

draw_accumulator <- function(case) {
  ht_result <- case$ht_result
  pred_lines <- case$pred_line_all %||% case$pred_line
  image(
    x = ht_result$theta,
    y = ht_result$rho,
    z = t(ht_result$accumulator),
    col = gray.colors(256, start = 1, end = 0),
    xlab = expression(theta),
    ylab = expression(rho),
    main = "Аккумуляторный массив",
    useRaster = TRUE,
    axes = FALSE
  )
  axis(1, at = pretty(ht_result$theta, n = 6))
  axis(2, at = pretty(ht_result$rho, n = 6))
  box()
  points(case$true_line[, "theta"], case$true_line[, "rho"], pch = 19, col = "#1f77b4", cex = 1.2)
  points(pred_lines[, "theta"], pred_lines[, "rho"], pch = 4, col = "#d62728", cex = 1.5, lwd = 2)
  legend(
    "topright",
    legend = c("Истинные параметры", "Оцененные параметры"),
    pch = c(19, 4),
    col = c("#1f77b4", "#d62728"),
    pt.cex = c(1.2, 1.5),
    bty = "n"
  )
  return(invisible(NULL))

  image(
    x = ht_result$theta,
    y = ht_result$rho,
    z = t(ht_result$accumulator),
    col = gray.colors(256, start = 1, end = 0),
    xlab = expression(theta),
    ylab = expression(rho),
    main = "Аккумуляторный массив",
    useRaster = TRUE
  )
  points(case$true_line[, "theta"], case$true_line[, "rho"], pch = 19, col = "#1f77b4", cex = 1.2)
  points(case$pred_line[, "theta"], case$pred_line[, "rho"], pch = 4, col = "#d62728", cex = 1.5, lwd = 2)
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
    noise_vec <- rnorm(n, sd = sigma)

    for (signal_index in seq_len(n)) {
      clean <- make_unit_row(n, signal_index)
      noisy <- clean + noise_vec
      est <- cssa_rank1_row(noisy)

      omega <- 2 * pi * (signal_index - 1L) / n
      omega_wrapped <- atan2(sin(omega), cos(omega))

      rows[[out_i]] <- data.frame(
        rep = rep_i,
        signal_index = signal_index,
        omega = omega,
        omega_wrapped = omega_wrapped,
        row_metrics(est, signal_index),
        stringsAsFactors = FALSE
      )
      out_i <- out_i + 1L
    }
  }

  dplyr::bind_rows(rows)
}

summarise_unit_frequency <- function(samples) {
  out <- aggregate(
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
  out[order(out$signal_index), ]
}
