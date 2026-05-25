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
if (is.null(common_file)) {
  common_file <- getOption("interactive_graphics.common_path", NULL)
}
common_dir <- if (!is.null(common_file)) {
  dirname(normalizePath(common_file, winslash = "/", mustWork = FALSE))
} else {
  getwd()
}
common_candidates <- unique(stats::na.omit(c(
  common_file,
  file.path(common_dir, "common.R"),
  file.path(dirname(common_dir), "common.R"),
  file.path(common_dir, "interactive_graphics", "common.R")
)))
common_path <- common_candidates[file.exists(common_candidates)][1]
if (is.na(common_path) || !nzchar(common_path)) {
  stop("Не удалось определить путь к interactive_graphics/common.R.")
}
common_path <- normalizePath(common_path, winslash = "/", mustWork = TRUE)
common_dir <- dirname(common_path)
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

available_worker_count <- function(reserve = 1L, max_workers = 8L) {
  cores <- parallel::detectCores(logical = TRUE)
  if (!is.finite(cores) || is.na(cores)) {
    cores <- 1L
  }
  cores <- max(1L, as.integer(cores) - as.integer(reserve))
  max_workers <- max(1L, as.integer(max_workers))
  min(cores, max_workers)
}

normalize_worker_count <- function(n_workers, n_tasks = Inf) {
  cores <- parallel::detectCores(logical = TRUE)
  if (!is.finite(cores) || is.na(cores)) {
    cores <- 1L
  }
  n_workers <- suppressWarnings(as.integer(n_workers %||% 1L))
  if (!is.finite(n_workers) || is.na(n_workers)) {
    n_workers <- 1L
  }
  n_workers <- max(1L, min(n_workers, as.integer(cores)))
  if (is.finite(n_tasks)) {
    n_workers <- min(n_workers, max(1L, as.integer(n_tasks)))
  }
  n_workers
}

make_app_cluster <- function(n_workers) {
  cl <- parallel::makeCluster(as.integer(n_workers))
  initialized <- FALSE
  on.exit({
    if (!initialized) {
      try(parallel::stopCluster(cl), silent = TRUE)
    }
  }, add = TRUE)

  init_worker <- function(path) {
    options(interactive_graphics.common_path = path)
    source(path, local = globalenv())
    NULL
  }
  for (worker_i in seq_along(cl)) {
    parallel::clusterCall(cl[worker_i], init_worker, common_path)
  }

  initialized <- TRUE
  cl
}

make_rep_seeds <- function(n) {
  sample.int(.Machine$integer.max, as.integer(n))
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

normalize_line_method <- function(line_method = "bresenham") {
  if (length(line_method) == 0L || is.na(line_method[1L])) {
    return("bresenham")
  }

  line_method <- as.character(line_method[1L])
  if (!line_method %in% c("bresenham", "default")) {
    return("bresenham")
  }

  line_method
}

make_line_image <- function(n_row, n_col, a, b, intensity = 1, line_method = "bresenham") {
  matrix(0, nrow = n_row, ncol = n_col) |>
    add.line(a = a, b = b, method = normalize_line_method(line_method), intensity = intensity)
}

cssa_component_limit <- function(n, L = (n + 1L) %/% 2L) {
  n <- as.integer(n)[1L]
  L <- as.integer(L)[1L]
  if (!is.finite(n) || is.na(n) || n < 1L) {
    return(1L)
  }
  if (!is.finite(L) || is.na(L) || L < 1L) {
    L <- (n + 1L) %/% 2L
  }
  K <- n - L + 1L
  max(1L, min(L, K))
}

cssa_component_limit_for_matrix <- function(m, method) {
  method <- match.arg(method, c("col.col", "col.row", "row.col", "row.row"))
  if (method %in% c("col.col", "col.row")) {
    return(cssa_component_limit(nrow(m)))
  }
  cssa_component_limit(ncol(m))
}

cssa_denoise <- function(m, num_of_lines = 1L, method = "row.row") {
  m <- as.matrix(m)
  num_of_lines <- suppressWarnings(as.integer(num_of_lines)[1L])
  if (!is.finite(num_of_lines) || is.na(num_of_lines) || num_of_lines < 1L) {
    num_of_lines <- 1L
  }
  num_of_lines <- min(num_of_lines, cssa_component_limit_for_matrix(m, method))

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

esprit_project_cssa_corrected <- function(x,
                                          L = (length(x) + 1L) %/% 2L,
                                          num_components = 1L,
                                          frequency_grid_multiplier = 1) {
  x <- as.vector(x)
  N <- length(x)
  num_components <- suppressWarnings(as.integer(num_components)[1L])
  if (!is.finite(num_components) || is.na(num_components) || num_components < 1L) {
    num_components <- 1L
  }
  num_components <- min(num_components, cssa_component_limit(N, L = L))
  frequency_grid_multiplier <- as.numeric(frequency_grid_multiplier)[1L]
  if (!is.finite(frequency_grid_multiplier) || frequency_grid_multiplier <= 0) {
    frequency_grid_multiplier <- 1
  }

  tryCatch({
    fit <- Rssa::ssa(x, L = L, kind = "cssa", svd.method = "svd")
    esprit_res <- Rssa::parestimate(
      fit,
      groups = list(signal = seq_len(num_components)),
      method = "esprit",
      normalize.roots = FALSE
    )

    num_components <- min(num_components, length(esprit_res$frequencies))
    if (num_components < 1L) {
      return(x)
    }

    signal_rec <- Rssa::reconstruct(
      fit,
      groups = list(signal = seq_len(num_components))
    )
    x_signal <- as.vector(signal_rec$signal)
    if (length(x_signal) != N) {
      x_signal <- x
    }

    grid_size <- frequency_grid_multiplier * N
    omega <- round(esprit_res$frequencies[seq_len(num_components)] * grid_size) / grid_size
    omega <- unique(omega)
    mu <- exp(2 * pi * omega * 1i)
    basis <- outer(0:(N - 1L), mu, function(i, z) z ^ i)
    amplitudes <- as.vector(qr.solve(basis, x_signal))

    as.vector(basis %*% amplitudes)
  }, error = function(e) {
    x
  })
}

esprit_cssa_denoise <- function(m,
                                num_of_lines = 1L,
                                method = "row.row",
                                frequency_grid_multiplier = 1) {
  method <- match.arg(method, c("row.row", "col.row"))
  m <- as.matrix(m)
  z <- dft(m)

  z_hat <- switch(
    method,
    "row.row" = {
      L <- (ncol(z) + 1L) %/% 2L
      t(vapply(seq_len(nrow(z)), function(i) {
        as.complex(esprit_project_cssa_corrected(
          z[i, ],
          L = L,
          num_components = num_of_lines,
          frequency_grid_multiplier = frequency_grid_multiplier
        ))
      }, complex(ncol(z))))
    },
    "col.row" = {
      L <- (nrow(z) + 1L) %/% 2L
      vapply(seq_len(ncol(z)), function(j) {
        as.complex(esprit_project_cssa_corrected(
          z[, j],
          L = L,
          num_components = num_of_lines,
          frequency_grid_multiplier = frequency_grid_multiplier
        ))
      }, complex(nrow(z)))
    }
  )

  clip01(Re(idft.row(z_hat)))
}

max_row_denoise <- function(m, k = 1L) {
  m <- as.matrix(m)
  k <- suppressWarnings(as.integer(k)[1L])
  if (!is.finite(k) || is.na(k)) {
    k <- 1L
  }
  k <- max(1L, min(k, ncol(m)))
  out <- matrix(0, nrow = nrow(m), ncol = ncol(m))

  for (i in seq_len(nrow(m))) {
    ind <- order(m[i, ], decreasing = TRUE)[seq_len(k)]
    out[i, ind] <- m[i, ind]
  }

  clip01(out)
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

available_detectors <- function(num_of_lines,
                                include_quantile = FALSE,
                                include_esprit = TRUE,
                                max_row_k = 1L,
                                cssa_components = num_of_lines,
                                esprit_components = cssa_components,
                                frequency_grid_multiplier = 1,
                                cfg = NULL,
                                line_method = "bresenham") {
  cssa_components <- suppressWarnings(as.integer(cssa_components)[1L])
  if (!is.finite(cssa_components) || is.na(cssa_components)) {
    cssa_components <- num_of_lines
  }
  cssa_components <- max(1L, cssa_components)
  esprit_components <- suppressWarnings(as.integer(esprit_components)[1L])
  if (!is.finite(esprit_components) || is.na(esprit_components)) {
    esprit_components <- cssa_components
  }
  esprit_components <- max(1L, esprit_components)

  detectors <- list(
    cssa_row_row = function(m, sigma_noise) {
      cssa_denoise(m, num_of_lines = cssa_components, method = "row.row")
    },
    cssa_col_row = function(m, sigma_noise) {
      cssa_denoise(m, num_of_lines = cssa_components, method = "col.row")
    },
    max_row = function(m, sigma_noise) {
      max_row_denoise(m, k = max_row_k)
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

  if (isTRUE(include_esprit)) {
    detectors <- append(
      detectors,
      list(
        cssa_row_row_esprit = function(m, sigma_noise) {
          esprit_cssa_denoise(
            m,
            num_of_lines = esprit_components,
            method = "row.row",
            frequency_grid_multiplier = frequency_grid_multiplier
          )
        },
        cssa_col_row_esprit = function(m, sigma_noise) {
          esprit_cssa_denoise(
            m,
            num_of_lines = esprit_components,
            method = "col.row",
            frequency_grid_multiplier = frequency_grid_multiplier
          )
        }
      ),
      after = 2L
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
    cssa_col_row = "CSSA col.row",
    cssa_row_row_esprit = "CSSA row.row + ESPRIT",
    cssa_col_row_esprit = "CSSA col.row + ESPRIT",
    max_row = "MAX.row",
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

compute_err_by_line <- function(true, pred) {
  true <- as.matrix(true)
  pred <- as.matrix(pred)

  if (ncol(true) != 2 || ncol(pred) != 2) {
    stop("true and pred must have two columns: rho and theta")
  }
  if (nrow(true) != nrow(pred)) {
    stop("true and pred must contain the same number of lines")
  }

  pair_err <- function(tr, tt, pr, pt) {
    tr <- unname(tr)
    tt <- unname(tt)
    pr <- unname(pr)
    pt <- unname(pt)

    pt_alt <- (pt + pi) %% (2 * pi)
    dtheta <- atan2(sin(tt - pt), cos(tt - pt))^2
    dtheta_alt <- atan2(sin(tt - pt_alt), cos(tt - pt_alt))^2

    loss_1 <- (tr - pr)^2 + dtheta
    loss_2 <- (tr + pr)^2 + dtheta_alt

    if (loss_2 < loss_1) {
      c(dr = (tr + pr)^2, dtheta = dtheta_alt, total = loss_2)
    } else {
      c(dr = (tr - pr)^2, dtheta = dtheta, total = loss_1)
    }
  }

  all_perms <- function(v) {
    if (length(v) == 1L) {
      return(matrix(v, nrow = 1L))
    }
    out <- lapply(seq_along(v), function(i) {
      cbind(v[i], all_perms(v[-i]))
    })
    do.call(rbind, out)
  }

  n_lines <- nrow(true)
  pair_errors <- vector("list", n_lines * n_lines)
  cost <- matrix(0, nrow = n_lines, ncol = n_lines)
  for (i in seq_len(n_lines)) {
    for (j in seq_len(n_lines)) {
      cur <- pair_err(true[i, 1], true[i, 2], pred[j, 1], pred[j, 2])
      pair_errors[[(i - 1L) * n_lines + j]] <- cur
      cost[i, j] <- cur["total"]
    }
  }

  if (n_lines == 1L) {
    best_perm <- 1L
  } else {
    perm_mat <- all_perms(seq_len(n_lines))
    costs <- apply(perm_mat, 1L, function(cur_perm) {
      sum(cost[cbind(seq_len(n_lines), cur_perm)])
    })
    best_perm <- perm_mat[which.min(costs), ]
  }

  rows <- lapply(seq_len(n_lines), function(i) {
    cur <- pair_errors[[(i - 1L) * n_lines + best_perm[i]]]
    data.frame(
      line_index = i,
      pred_index = best_perm[i],
      true_rho = true[i, 1],
      true_theta = true[i, 2],
      pred_rho = pred[best_perm[i], 1],
      pred_theta = pred[best_perm[i], 2],
      dr = unname(cur["dr"]),
      dtheta = unname(cur["dtheta"]),
      stringsAsFactors = FALSE
    )
  })

  out <- dplyr::bind_rows(rows)
  rownames(out) <- NULL
  out
}

compute_err_with_missing <- function(true, pred, rho_step_ht = 1, theta_step_ht = 0.01) {
  true <- as.matrix(true)
  pred <- as.matrix(pred)
  if (ncol(true) != 2 || ncol(pred) != 2) {
    stop("true and pred must have two columns: rho and theta")
  }

  n_true <- nrow(true)
  n_pred <- nrow(pred)
  if (n_pred == n_true) {
    return(compute_err(true, pred))
  }

  if (n_pred > n_true) {
    best <- c(dr = Inf, dtheta = Inf)
    pred_sets <- utils::combn(seq_len(n_pred), n_true, simplify = FALSE)
    for (pred_idx in pred_sets) {
      cur <- compute_err(true, pred[pred_idx, , drop = FALSE])
      if (sum(cur) < sum(best)) {
        best <- cur
      }
    }
    return(best)
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

  all_perms <- function(v) {
    if (length(v) == 1L) {
      return(matrix(v, nrow = 1L))
    }
    out <- lapply(seq_along(v), function(i) {
      cbind(v[i], all_perms(v[-i]))
    })
    do.call(rbind, out)
  }

  err <- c(dr = 0, dtheta = 0)
  if (n_pred > 0L) {
    best <- c(dr = Inf, dtheta = Inf, total = Inf)
    true_sets <- utils::combn(seq_len(n_true), n_pred, simplify = FALSE)
    pred_perms <- all_perms(seq_len(n_pred))
    for (true_idx in true_sets) {
      for (row_i in seq_len(nrow(pred_perms))) {
        pred_idx <- pred_perms[row_i, ]
        cur <- c(dr = 0, dtheta = 0, total = 0)
        for (i in seq_along(true_idx)) {
          pair <- pair_err(
            true[true_idx[i], 1], true[true_idx[i], 2],
            pred[pred_idx[i], 1], pred[pred_idx[i], 2]
          )
          cur <- cur + pair
        }
        if (cur["total"] < best["total"]) {
          best <- cur
        }
      }
    }
    err <- c(dr = unname(best["dr"]), dtheta = unname(best["dtheta"]))
  }

  missing_n <- n_true - n_pred
  c(
    dr = unname(err["dr"]) + missing_n * max(as.numeric(rho_step_ht), 1)^2,
    dtheta = unname(err["dtheta"]) + missing_n * max(as.numeric(theta_step_ht), 0.01)^2
  )
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

split_custom_config_blocks <- function(text) {
  raw_lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  blocks <- list()
  current <- character(0)

  for (line in raw_lines) {
    line <- trimws(line)
    if (!nzchar(line)) {
      if (length(current) > 0L) {
        blocks[[length(blocks) + 1L]] <- current
        current <- character(0)
      }
    } else {
      current <- c(current, line)
    }
  }

  if (length(current) > 0L) {
    blocks[[length(blocks) + 1L]] <- current
  }
  blocks
}

line_config_ids <- function(lines_df) {
  lines_df <- as.data.frame(lines_df)
  if (!"intensity" %in% names(lines_df)) {
    lines_df$intensity <- 1
  }
  sprintf("(%g, %g, %g)", lines_df$a, lines_df$b, lines_df$intensity)
}

line_display_ids <- function(lines_df) {
  ids <- line_config_ids(lines_df)
  sprintf("line %d: %s", seq_along(ids), ids)
}

integer_lcm <- function(x) {
  x <- unique(abs(as.integer(x)))
  x <- x[is.finite(x) & !is.na(x) & x > 0L]
  if (length(x) == 0L) {
    return(1L)
  }

  gcd_int <- function(a, b) {
    a <- abs(as.integer(a))
    b <- abs(as.integer(b))
    while (b != 0L) {
      r <- a %% b
      a <- b
      b <- r
    }
    a
  }

  as.integer(Reduce(function(a, b) abs(a * b) / gcd_int(a, b), x))
}

esprit_grid_multiplier_from_lines <- function(lines_df) {
  lines_df <- as.data.frame(lines_df)
  if (!"a" %in% names(lines_df)) {
    return(1L)
  }

  a <- suppressWarnings(as.numeric(lines_df$a))
  a <- a[is.finite(a) & !is.na(a)]
  if (length(a) == 0L) {
    return(1L)
  }

  a_int <- round(abs(a))
  a_int[a_int < 1L] <- 1L
  integer_lcm(a_int)
}

cssa_colrow_safe_rank_from_lines <- function(lines_df, n_row, n_col) {
  lines_df <- as.data.frame(lines_df)
  if (!all(c("a", "b") %in% names(lines_df))) {
    return(1L)
  }

  a <- suppressWarnings(as.numeric(lines_df$a))
  b <- suppressWarnings(as.numeric(lines_df$b))
  keep <- is.finite(a) & !is.na(a) & is.finite(b) & !is.na(b)
  if (!any(keep)) {
    return(1L)
  }

  n_row <- as.integer(n_row)[1L]
  n_col <- as.integer(n_col)[1L]
  if (!is.finite(n_row) || is.na(n_row) || n_row < 1L) {
    n_row <- 1L
  }
  if (!is.finite(n_col) || is.na(n_col) || n_col < 1L) {
    n_col <- 1L
  }

  L <- (n_row + 1L) %/% 2L
  K <- n_row - L + 1L
  ranks <- vapply(which(keep), function(i) {
    cur_a <- a[i]
    cur_b <- b[i]
    alpha <- abs(cur_a)
    if (alpha < 1) {
      alpha <- 1
    }

    if (cur_a > 0) {
      y1 <- cur_a + cur_b
      yN <- cur_a * n_col + cur_b
    } else {
      y1 <- n_row + 1 - (cur_a + cur_b)
      yN <- n_row + 1 - (cur_a * n_col + cur_b)
    }

    max(1L, as.integer(floor(min(
      L,
      K,
      max(yN, 0),
      max(n_row + 1 - y1, 0),
      max(alpha, y1) + max(n_row + 1 - yN - alpha, 0)
    ))))
  }, integer(1))

  max(1L, as.integer(sum(ranks)))
}

cssa_rowrow_auto_rank_from_lines <- function(lines_df) {
  max(1L, nrow(as.data.frame(lines_df)))
}

esprit_grid_multiplier_for_method <- function(lines_df, method_name) {
  if (grepl("^cssa_row_row", method_name)) {
    return(1L)
  }
  esprit_grid_multiplier_from_lines(lines_df)
}

cssa_auto_rank_for_method <- function(cfg, method_name, line_method = "bresenham") {
  if (is.null(cfg) || is.null(cfg$lines)) {
    return(1L)
  }

  if (grepl("^cssa_col_row", method_name)) {
    return(cssa_colrow_safe_rank_from_lines(cfg$lines, cfg$n_row, cfg$n_col))
  }
  if (grepl("^cssa_row_row", method_name)) {
    return(cssa_rowrow_auto_rank_from_lines(cfg$lines))
  }
  cssa_rowrow_auto_rank_from_lines(cfg$lines)
}

cssa_rank_limit_for_config <- function(cfg, method_name) {
  if (is.null(cfg)) {
    return(1L)
  }
  n_row <- as.integer(cfg$n_row %||% 1L)[1L]
  n_col <- as.integer(cfg$n_col %||% 1L)[1L]
  if (grepl("^cssa_col_row|^cssa_col_col", method_name)) {
    return(cssa_component_limit(n_row))
  }
  if (grepl("^cssa_row_row|^cssa_row_col", method_name)) {
    return(cssa_component_limit(n_col))
  }
  max(cssa_component_limit(n_row), cssa_component_limit(n_col))
}

cssa_auto_rank_pair <- function(cfg, line_method = "bresenham") {
  if (is.null(cfg) || is.null(cfg$lines)) {
    return(c(row_row = 1L, col_row = 1L))
  }
  c(
    row_row = min(cssa_rowrow_auto_rank_from_lines(cfg$lines), cssa_component_limit(cfg$n_col)),
    col_row = min(cssa_colrow_safe_rank_from_lines(cfg$lines, cfg$n_row, cfg$n_col), cssa_component_limit(cfg$n_row))
  )
}

parse_cssa_rank_pair <- function(cssa_components) {
  if (is.null(cssa_components) || length(cssa_components) == 0L) {
    return(c(row_row = NA_integer_, col_row = NA_integer_))
  }

  text <- paste(as.character(cssa_components), collapse = ",")
  parts <- unlist(strsplit(text, "[,;[:space:]]+", perl = TRUE))
  parts <- parts[nzchar(parts)]
  values <- suppressWarnings(as.integer(parts))
  values <- values[is.finite(values) & !is.na(values) & values >= 1L]
  if (length(values) == 0L) {
    return(c(row_row = NA_integer_, col_row = NA_integer_))
  }
  if (length(values) == 1L) {
    values <- rep(values, 2L)
  }
  c(row_row = values[[1L]], col_row = values[[2L]])
}

resolve_cssa_components <- function(cssa_components, cfg, method_name, line_method = "bresenham") {
  auto <- cssa_auto_rank_for_method(cfg, method_name, line_method = line_method)
  limit <- cssa_rank_limit_for_config(cfg, method_name)
  auto <- min(auto, limit)
  manual_pair <- parse_cssa_rank_pair(cssa_components)
  manual <- if (grepl("^cssa_col_row", method_name)) {
    manual_pair[["col_row"]]
  } else if (grepl("^cssa_row_row", method_name)) {
    manual_pair[["row_row"]]
  } else {
    manual_pair[["row_row"]]
  }

  if (!is.finite(manual) || is.na(manual) || manual < 1L) {
    return(list(value = auto, source = "auto", auto = auto))
  }
  value <- min(max(1L, manual), limit)
  source <- if (value < manual) "manual_clamped" else "manual"
  list(value = value, source = source, auto = auto)
}

custom_config_id <- function(lines_df) {
  paste(line_config_ids(lines_df), collapse = "; ")
}

parse_custom_config_list <- function(text, num_lines, n_row, n_col) {
  num_lines <- max(1L, as.integer(num_lines)[1L])
  blocks <- split_custom_config_blocks(text)
  if (length(blocks) == 0L) {
    stop("РЎРїРёСЃРѕРє РїСЂСЏРјС‹С… РїСѓСЃС‚.")
  }

  if (length(blocks) == 1L) {
    lines_df <- parse_line_text(paste(blocks[[1L]], collapse = "\n"))
    if (nrow(lines_df) %% num_lines != 0L) {
      stop(
        "Р§РёСЃР»Рѕ СЃС‚СЂРѕРє РІ custom-РєРѕРЅС„РёРіСѓСЂР°С†РёРё РґРѕР»Р¶РЅРѕ РґРµР»РёС‚СЊСЃСЏ РЅР° РєРѕР»РёС‡РµСЃС‚РІРѕ РїСЂСЏРјС‹С…: ",
        num_lines,
        "."
      )
    }
    idx <- split(seq_len(nrow(lines_df)), ceiling(seq_len(nrow(lines_df)) / num_lines))
    blocks <- lapply(idx, function(i) lines_df[i, , drop = FALSE])
  } else {
    blocks <- lapply(blocks, function(block) {
      parse_line_text(paste(block, collapse = "\n"), expected_n = num_lines)
    })
  }

  n_configs <- length(blocks)
  lapply(seq_along(blocks), function(i) {
    config_label <- custom_config_id(blocks[[i]])
    make_custom_config(
      blocks[[i]],
      n_row = n_row,
      n_col = n_col,
      config_id = config_label,
      description = config_label
    )
  })
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
      intensity = 0.8
    )
  }) |>
    lapply(function(cfg) {
      cfg$lines$intensity <- 0.8
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

  parse_custom_config_list(
    text = line_text,
    num_lines = num_lines,
    n_row = n_row,
    n_col = n_col
  )
}

config_to_matrix <- function(cfg, line_method = "bresenham") {
  line_method <- normalize_line_method(line_method)
  out <- matrix(0, nrow = cfg$n_row, ncol = cfg$n_col)
  for (i in seq_len(nrow(cfg$lines))) {
    out <- add.line(
      out,
      a = cfg$lines$a[i],
      b = cfg$lines$b[i],
      method = line_method,
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

draw_config_preview <- function(config_list, line_method = "bresenham") {
  mats <- lapply(config_list, config_to_matrix, line_method = line_method)
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
  threshold_quantile_p = 0.9,
  threshold_multiplier = 1,
  residual_sd = NULL,
  cssa_components = num_of_lines
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
    thr <- suppressWarnings(as.numeric(threshold_value)[1L])
    if (!is.finite(thr) || is.na(thr)) {
      thr <- 0
    }
    thr <- max(0, thr)
    out <- processed_matrix
    out[out < thr] <- 0
    return(list(
      processed = out,
      threshold_value = thr,
      residual_sd = NA_real_
    ))
  }

  if (threshold_mode == "quantile") {
    p <- as.numeric(threshold_quantile_p %||% 0.9)
    if (!is.finite(p)) {
      p <- 0.9
    }
    p <- min(max(p, 0), 1)
    thr <- as.numeric(stats::quantile(processed_matrix, probs = p, na.rm = TRUE, names = FALSE))
    out <- processed_matrix
    out[out < thr] <- 0
    return(list(
      processed = out,
      threshold_value = thr,
      residual_sd = NA_real_
    ))
  }

  if (is.null(residual_sd)) {
    rowrow_signal <- cssa_denoise(noisy_matrix, num_of_lines = cssa_components, method = "row.row")
    residual_sd <- stats::sd(as.vector(noisy_matrix - rowrow_signal), na.rm = TRUE)
  }
  threshold_multiplier <- suppressWarnings(as.numeric(threshold_multiplier)[1L])
  if (!is.finite(threshold_multiplier) || is.na(threshold_multiplier)) {
    threshold_multiplier <- 1
  }
  threshold_multiplier <- max(0, threshold_multiplier)
  if (!is.finite(residual_sd) || is.na(residual_sd)) {
    residual_sd <- 0
  }
  thr <- threshold_multiplier * residual_sd
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

run_ht_repetition_app <- function(
  cfg,
  rep_i,
  seed,
  method_names,
  sigma_noise,
  num_of_lines,
  rho_step_ht,
  theta_step_ht,
  ht_type,
  threshold_mode,
  threshold_value,
  threshold_quantile_p,
  threshold_multiplier,
  threshold_method_names,
  max_row_k = 1L,
  cssa_components = NULL,
  line_method = "bresenham"
) {
  set.seed(seed)
  num_of_lines <- nrow(cfg$lines)
  base_matrix <- config_to_matrix(cfg, line_method = line_method)
  true_line <- config_true_params(cfg)
  noisy_matrix <- add.noise(base_matrix, sigma = sigma_noise)
  threshold_enabled <- !identical(threshold_mode, "none") && length(threshold_method_names) > 0L
  residual_sd <- NULL
  if (threshold_enabled && identical(threshold_mode, "auto_rowrow_sd")) {
    residual_rank <- resolve_cssa_components(cssa_components, cfg, "cssa_row_row", line_method = line_method)
    rowrow_signal <- cssa_denoise(noisy_matrix, num_of_lines = residual_rank$value, method = "row.row")
    residual_sd <- stats::sd(as.vector(noisy_matrix - rowrow_signal), na.rm = TRUE)
  }

  weighted <- identical(ht_type, "weighted")
  line_ids <- line_display_ids(cfg$lines)
  rows <- vector("list", length(method_names))
  for (method_i in seq_along(method_names)) {
    method_name <- method_names[[method_i]]
    rank_info <- resolve_cssa_components(cssa_components, cfg, method_name, line_method = line_method)
    method_frequency_grid_multiplier <- esprit_grid_multiplier_for_method(cfg$lines, method_name)
    detectors <- available_detectors(
      num_of_lines = num_of_lines,
      max_row_k = max_row_k,
      cssa_components = rank_info$value,
      esprit_components = rank_info$value,
      frequency_grid_multiplier = method_frequency_grid_multiplier,
      cfg = cfg,
      line_method = line_method
    )
    processed_matrix <- detectors[[method_name]](noisy_matrix, sigma_noise = sigma_noise)
    cur_threshold_mode <- if (method_name %in% threshold_method_names) threshold_mode else "none"
    threshold_info <- threshold_processed_matrix(
      processed_matrix = processed_matrix,
      noisy_matrix = noisy_matrix,
      num_of_lines = num_of_lines,
      threshold_mode = cur_threshold_mode,
      threshold_value = threshold_value,
      threshold_quantile_p = threshold_quantile_p,
      threshold_multiplier = threshold_multiplier,
      residual_sd = residual_sd,
      cssa_components = rank_info$value
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

    err_by_line <- compute_err_by_line(true_line, pred_line)
    cur_line_ids <- line_ids[err_by_line$line_index]

    rows[[method_i]] <- data.frame(
      config_id = cur_line_ids,
      parent_config_id = cfg$config_id,
      description = cur_line_ids,
      line_index = err_by_line$line_index,
      pred_index = err_by_line$pred_index,
      num_lines = num_of_lines,
      cssa_components = rank_info$value,
      cssa_rank_source = rank_info$source,
      cssa_auto_rank = rank_info$auto,
      frequency_grid_multiplier = method_frequency_grid_multiplier,
      n_col = cfg$n_col,
      n_row = cfg$n_row,
      sigma = sigma_noise,
      rep = rep_i,
      method = method_name,
      true_rho = err_by_line$true_rho,
      true_theta = err_by_line$true_theta,
      pred_rho = err_by_line$pred_rho,
      pred_theta = err_by_line$pred_theta,
      active_pixels = ht_result$active_pixels,
      active_weight = ht_result$active_weight,
      threshold_value = threshold_info$threshold_value,
      dr = err_by_line$dr,
      dtheta = err_by_line$dtheta,
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(rows)
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
  threshold_quantile_p = 0.9,
  threshold_multiplier = 1,
  threshold_method_names = NULL,
  progress = NULL,
  n_workers = 1L,
  max_row_k = 1L,
  cssa_components = num_of_lines,
  line_method = "bresenham"
) {
  detectors <- available_detectors(
    num_of_lines = num_of_lines,
    max_row_k = max_row_k
  )
  method_names <- intersect(method_names, names(detectors))
  if (length(method_names) == 0L) {
    stop("Не выбрано ни одного доступного метода.")
  }
  if (is.null(threshold_method_names)) {
    threshold_method_names <- method_names
  } else {
    threshold_method_names <- intersect(threshold_method_names, method_names)
  }

  tasks <- unlist(lapply(seq_along(config_list), function(cfg_i) {
    lapply(seq_len(n_rep), function(rep_i) {
      list(cfg = config_list[[cfg_i]], rep_i = rep_i)
    })
  }), recursive = FALSE)
  seeds <- make_rep_seeds(length(tasks))
  for (i in seq_along(tasks)) {
    tasks[[i]]$seed <- seeds[[i]]
  }

  n_workers <- normalize_worker_count(n_workers, length(tasks))
  worker_args <- list(
    method_names = method_names,
    sigma_noise = sigma_noise,
    num_of_lines = num_of_lines,
    rho_step_ht = rho_step_ht,
    theta_step_ht = theta_step_ht,
    ht_type = ht_type,
    threshold_mode = threshold_mode,
    threshold_value = threshold_value,
    threshold_quantile_p = threshold_quantile_p,
    threshold_multiplier = threshold_multiplier,
    threshold_method_names = threshold_method_names,
    max_row_k = max_row_k,
    cssa_components = cssa_components,
    line_method = line_method
  )

  if (n_workers > 1L) {
    if (!is.null(progress)) {
      progress(0.05, detail = sprintf("Параллельный запуск: %d воркеров", n_workers))
    }
    cl <- make_app_cluster(n_workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    results_list <- parallel::parLapplyLB(cl, tasks, function(task, args) {
      do.call(
        run_ht_repetition_app,
        c(
          list(cfg = task$cfg, rep_i = task$rep_i, seed = task$seed),
          args
        )
      )
    }, worker_args)
    if (!is.null(progress)) {
      progress(1, detail = "Параллельный расчет завершен")
    }
    return(dplyr::bind_rows(results_list))
  }

  results_list <- vector("list", length(tasks))
  step_i <- 0L

  for (task_i in seq_along(tasks)) {
    task <- tasks[[task_i]]
    results_list[[task_i]] <- do.call(
      run_ht_repetition_app,
      c(
        list(cfg = task$cfg, rep_i = task$rep_i, seed = task$seed),
        worker_args
      )
    )
    step_i <- step_i + 1L
    if (!is.null(progress)) {
      progress(step_i / length(tasks), detail = sprintf("%s / rep %d", task$cfg$config_id, task$rep_i))
    }
  }

  dplyr::bind_rows(results_list)
}

summarise_ht_errors_app <- function(ht_error_samples) {
  ht_error_samples |>
    dplyr::group_by(config_id, description, method) |>
    dplyr::summarise(
      parent_config_id = paste(unique(parent_config_id), collapse = " | "),
      line_index = dplyr::first(line_index),
      num_lines = dplyr::first(num_lines),
      cssa_components = dplyr::first(cssa_components),
      cssa_rank_source = dplyr::first(cssa_rank_source),
      cssa_auto_rank = dplyr::first(cssa_auto_rank),
      frequency_grid_multiplier = dplyr::first(frequency_grid_multiplier),
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

ideal_discretization_table_app <- function(config_list, rho_step_ht, theta_step_ht) {
  rows <- lapply(config_list, function(cfg) {
    true_line <- config_true_params(cfg)
    pred_line <- cbind(
      rho = round(true_line[, "rho"] / rho_step_ht) * rho_step_ht,
      theta = round(true_line[, "theta"] / theta_step_ht) * theta_step_ht
    )
    pred_line[, "theta"] <- pmin(pmax(pred_line[, "theta"], 0), pi)

    err_by_line <- compute_err_by_line(true_line, pred_line)
    line_ids <- line_display_ids(cfg$lines)

    data.frame(
      config_id = line_ids[err_by_line$line_index],
      parent_config_id = cfg$config_id,
      line_index = err_by_line$line_index,
      true_rho = err_by_line$true_rho,
      true_theta = err_by_line$true_theta,
      ideal_rho = err_by_line$pred_rho,
      ideal_theta = err_by_line$pred_theta,
      ideal_dr = err_by_line$dr,
      ideal_dtheta = err_by_line$dtheta,
      baseline_dr = pmax(err_by_line$dr, (rho_step_ht / 2)^2),
      baseline_dtheta = pmax(err_by_line$dtheta, (theta_step_ht / 2)^2),
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(rows)
}

find_big_error_rep_app <- function(
  cfg,
  rep_i,
  seed,
  method_name,
  sigma_noise,
  num_of_lines,
  num_maxima,
  rho_step_ht,
  theta_step_ht,
  ht_type,
  threshold_mode,
  threshold_value,
  threshold_quantile_p,
  threshold_multiplier,
  factor_threshold,
  ideal,
  max_row_k = 1L,
  cssa_components = num_of_lines,
  line_method = "bresenham"
) {
  set.seed(seed)
  base_matrix <- config_to_matrix(cfg, line_method = line_method)
  true_line <- config_true_params(cfg)
  weighted <- identical(ht_type, "weighted")

  noisy_matrix <- add.noise(base_matrix, sigma = sigma_noise)
  rank_info <- resolve_cssa_components(cssa_components, cfg, method_name, line_method = line_method)
  residual_sd <- NULL
  if (identical(threshold_mode, "auto_rowrow_sd")) {
    residual_rank <- resolve_cssa_components(cssa_components, cfg, "cssa_row_row", line_method = line_method)
    rowrow_signal <- cssa_denoise(noisy_matrix, num_of_lines = residual_rank$value, method = "row.row")
    residual_sd <- stats::sd(as.vector(noisy_matrix - rowrow_signal), na.rm = TRUE)
  }
  frequency_grid_multiplier <- esprit_grid_multiplier_for_method(cfg$lines, method_name)
  detectors <- available_detectors(
    num_of_lines = num_of_lines,
    max_row_k = max_row_k,
    cssa_components = rank_info$value,
    esprit_components = rank_info$value,
    frequency_grid_multiplier = frequency_grid_multiplier,
    cfg = cfg,
    line_method = line_method
  )
  processed_matrix <- detectors[[method_name]](noisy_matrix, sigma_noise = sigma_noise)
  threshold_info <- threshold_processed_matrix(
    processed_matrix = processed_matrix,
    noisy_matrix = noisy_matrix,
    num_of_lines = num_of_lines,
    threshold_mode = threshold_mode,
    threshold_value = threshold_value,
    threshold_quantile_p = threshold_quantile_p,
    threshold_multiplier = threshold_multiplier,
    residual_sd = residual_sd,
    cssa_components = rank_info$value
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

  err <- compute_err_with_missing(
    true_line,
    pred_line_all,
    rho_step_ht = rho_step_ht,
    theta_step_ht = theta_step_ht
  )
  dr_ratio <- as.numeric(err["dr"]) / max(as.numeric(ideal$dr), 1e-12)
  dtheta_ratio <- as.numeric(err["dtheta"]) / max(as.numeric(ideal$dtheta), 1e-12)

  if (dr_ratio < factor_threshold && dtheta_ratio < factor_threshold) {
    return(NULL)
  }

  list(
    case_id = sprintf("rep_%03d", rep_i),
    rep = rep_i,
    dr = as.numeric(err["dr"]),
    dtheta = as.numeric(err["dtheta"]),
    dr_ratio = dr_ratio,
    dtheta_ratio = dtheta_ratio,
    num_maxima = num_maxima,
    detected_lines = nrow(pred_line_all),
    cssa_components = rank_info$value,
    cssa_rank_source = rank_info$source,
    cssa_auto_rank = rank_info$auto,
    frequency_grid_multiplier = frequency_grid_multiplier,
    threshold_value = threshold_info$threshold_value,
    noisy_matrix = noisy_matrix,
    processed_matrix = threshold_info$processed,
    ht_result = ht_result,
    pred_line = pred_line_all,
    pred_line_all = pred_line_all,
    true_line = true_line
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
  threshold_quantile_p = 0.9,
  threshold_multiplier = 1,
  n_search,
  factor_threshold,
  max_cases,
  progress = NULL,
  n_workers = 1L,
  max_row_k = 1L,
  cssa_components = num_of_lines,
  line_method = "bresenham"
) {
  detectors <- available_detectors(
    num_of_lines = num_of_lines,
    max_row_k = max_row_k
  )
  if (!method_name %in% names(detectors)) {
    stop("Метод ", method_name, " недоступен.")
  }

  num_maxima <- suppressWarnings(as.integer(num_maxima)[1L])
  if (!is.finite(num_maxima) || is.na(num_maxima)) {
    num_maxima <- 1L
  }
  num_maxima <- max(num_maxima, 1L)
  base_matrix <- config_to_matrix(cfg, line_method = line_method)
  true_line <- config_true_params(cfg)
  ideal <- ideal_discretization_error(true_line, rho_step_ht, theta_step_ht)
  cases <- list()
  rep_ids <- seq_len(n_search)
  seeds <- make_rep_seeds(length(rep_ids))
  n_workers <- normalize_worker_count(n_workers, length(rep_ids))

  worker_args <- list(
    cfg = cfg,
    method_name = method_name,
    sigma_noise = sigma_noise,
    num_of_lines = num_of_lines,
    num_maxima = num_maxima,
    rho_step_ht = rho_step_ht,
    theta_step_ht = theta_step_ht,
    ht_type = ht_type,
    threshold_mode = threshold_mode,
    threshold_value = threshold_value,
    threshold_quantile_p = threshold_quantile_p,
    threshold_multiplier = threshold_multiplier,
    factor_threshold = factor_threshold,
    ideal = ideal,
    max_row_k = max_row_k,
    cssa_components = cssa_components,
    line_method = line_method
  )

  if (n_workers > 1L) {
    if (!is.null(progress)) {
      progress(0.05, detail = sprintf("Параллельный поиск: %d воркеров", n_workers))
    }
    tasks <- Map(function(rep_i, seed) list(rep_i = rep_i, seed = seed), rep_ids, seeds)
    parallel_cases <- tryCatch({
      cl <- make_app_cluster(n_workers)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      found <- parallel::parLapplyLB(cl, tasks, function(task, args) {
        do.call(
          find_big_error_rep_app,
          c(
            list(rep_i = task$rep_i, seed = task$seed),
            args
          )
        )
      }, worker_args)
      found_cases <- found[!vapply(found, is.null, logical(1))]
      if (length(found_cases) > 0L) {
        case_reps <- vapply(found_cases, `[[`, integer(1), "rep")
        found_cases <- found_cases[order(case_reps)]
        found_cases <- found_cases[seq_len(min(length(found_cases), max_cases))]
      }
      found_cases
    }, error = function(e) {
      if (!is.null(progress)) {
        progress(0.05, detail = paste("Parallel search failed; falling back to one worker:", conditionMessage(e)))
      }
      NULL
    })
    if (!is.null(parallel_cases)) {
      if (!is.null(progress)) {
        progress(1, detail = sprintf("Найдено случаев: %d", length(parallel_cases)))
      }
      return(list(ideal = ideal, cases = parallel_cases))
    }
  }

  for (rep_i in seq_len(n_search)) {
    if (length(cases) >= max_cases) {
      break
    }

    case <- do.call(
      find_big_error_rep_app,
      c(
        list(rep_i = rep_i, seed = seeds[[rep_i]]),
        worker_args
      )
    )
    if (!is.null(case)) {
      cases[[length(cases) + 1L]] <- case
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

run_unit_frequency_rep <- function(rep_i, rep_seed, n, sigma) {
  set.seed(rep_seed)
  noise_vec <- rnorm(n, sd = sigma)
  rows <- vector("list", n)

  for (signal_index in seq_len(n)) {
    clean <- make_unit_row(n, signal_index)
    noisy <- clean + noise_vec
    est <- cssa_rank1_row(noisy)

    omega <- 2 * pi * (signal_index - 1L) / n
    omega_wrapped <- atan2(sin(omega), cos(omega))

    rows[[signal_index]] <- data.frame(
      rep = rep_i,
      signal_index = signal_index,
      omega = omega,
      omega_wrapped = omega_wrapped,
      row_metrics(est, signal_index),
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(rows)
}

run_unit_frequency_experiment <- function(
  n = 100L,
  n_rep = 1L,
  sigma = 0.2,
  seed = 1111L,
  n_workers = 1L,
  progress = NULL
) {
  set.seed(seed)
  rep_seeds <- make_rep_seeds(n_rep)
  rep_ids <- seq_len(n_rep)
  n_workers <- normalize_worker_count(n_workers, n_rep)

  if (n_workers > 1L) {
    if (!is.null(progress)) {
      progress(0.05, detail = sprintf("Параллельный расчет: %d воркеров", n_workers))
    }
    tasks <- Map(function(rep_i, rep_seed) list(rep_i = rep_i, rep_seed = rep_seed), rep_ids, rep_seeds)
    cl <- make_app_cluster(n_workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    rows <- parallel::parLapplyLB(cl, tasks, function(task, n, sigma) {
      run_unit_frequency_rep(
        rep_i = task$rep_i,
        rep_seed = task$rep_seed,
        n = n,
        sigma = sigma
      )
    }, n = n, sigma = sigma)
    if (!is.null(progress)) {
      progress(1, detail = "Параллельный расчет завершен")
    }
    return(dplyr::bind_rows(rows))
  }

  rows <- vector("list", n_rep)
  for (rep_i in rep_ids) {
    rows[[rep_i]] <- run_unit_frequency_rep(
      rep_i = rep_i,
      rep_seed = rep_seeds[[rep_i]],
      n = n,
      sigma = sigma
    )
    if (!is.null(progress)) {
      progress(rep_i / n_rep, detail = sprintf("rep %d / %d", rep_i, n_rep))
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
