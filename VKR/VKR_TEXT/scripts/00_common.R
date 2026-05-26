find_repo_root <- function(start_dir = getwd()) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in seq_len(30L)) {
    if (file.exists(file.path(cur, "SSA-in-Image-processing.Rproj")) ||
        file.exists(file.path(cur, "ssa-based methods", "cssa-transform.r"))) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) {
      break
    }
    cur <- parent
  }
  stop("Cannot find repository root.")
}

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/")))
  }
  if (!is.null(sys.frames()[[1L]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1L]]$ofile, winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- get_script_dir()
repo_root <- find_repo_root(script_dir)
report_root <- file.path(repo_root, "VKR", "VKR_TEXT")

chapter1_image_dir <- file.path(report_root, "assets", "images", "chapter1")
chapter2_image_dir <- file.path(report_root, "assets", "images", "chapter2")
chapter3_image_dir <- file.path(report_root, "assets", "images", "chapter3")
chapter3_table_dir <- file.path(report_root, "chapters", "chapter3_experiments", "tables")
chapter3_data_dir <- file.path(chapter3_table_dir, "data")

ensure_report_dirs <- function() {
  dirs <- c(
    chapter1_image_dir,
    chapter2_image_dir,
    chapter3_image_dir,
    chapter3_table_dir,
    chapter3_data_dir
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

load_project_code <- function(load_hough = TRUE) {
  required <- c("Rssa", "Rcpp")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing R packages: ", paste(missing, collapse = ", "))
  }

  source(file.path(repo_root, "ssa-based methods", "cssa-transform.r"))
  if (isTRUE(load_hough)) {
    source(file.path(repo_root, "hough transform", "hough_transform.r"))
  }
  invisible(TRUE)
}

clip01 <- function(m, lower = 0.001, upper = 0.99) {
  x <- as.matrix(m)
  x[!is.finite(x)] <- 0
  x[x < lower] <- 0
  x[x > upper] <- 1
  x
}

make_line_image <- function(n_row = 100L,
                            n_col = 100L,
                            lines,
                            line_method = "default") {
  lines <- normalize_config_lines(lines)
  out <- matrix(0, nrow = n_row, ncol = n_col)
  for (i in seq_len(nrow(lines))) {
    out <- add.line(
      out,
      a = lines$a[i],
      b = lines$b[i],
      method = line_method,
      intensity = lines$intensity[i]
    )
  }
  out
}

normalize_config_lines <- function(lines) {
  if (is.null(lines)) {
    stop("lines must not be NULL")
  }
  lines <- as.data.frame(lines)
  if (!"a" %in% names(lines) || !"b" %in% names(lines)) {
    stop("lines must contain columns a and b")
  }
  if (!"intensity" %in% names(lines)) {
    lines$intensity <- 1
  }
  lines$a <- as.numeric(lines$a)
  lines$b <- as.numeric(lines$b)
  lines$intensity <- as.numeric(lines$intensity)
  lines
}

line_config_labels <- function(lines) {
  lines <- normalize_config_lines(lines)
  sprintf("(%g, %g, %g)", lines$a, lines$b, lines$intensity)
}

line_config_label <- function(lines) {
  paste(line_config_labels(lines), collapse = "; ")
}

format_config_label <- function(config_id) {
  compact <- gsub("\\s+", "", config_id)
  parts <- strsplit(gsub("[()]", "", compact), ",", fixed = FALSE)[[1L]]
  if (length(parts) >= 3L) {
    paste0("(", parts[1L], ",", parts[2L], ";c=", parts[3L], ")")
  } else {
    compact
  }
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

config_true_lines <- function(lines) {
  lines <- normalize_config_lines(lines)
  out <- t(vapply(seq_len(nrow(lines)), function(i) {
    convert_ab_to_rho_theta(lines$a[i], lines$b[i])
  }, numeric(2L)))
  colnames(out) <- c("rho", "theta")
  out
}

cssa_component_limit <- function(n, L = (n + 1L) %/% 2L) {
  n <- as.integer(n)[1L]
  L <- as.integer(L)[1L]
  K <- n - L + 1L
  max(1L, min(L, K))
}

gcd_int <- function(a, b) {
  a <- abs(as.integer(a))
  b <- abs(as.integer(b))
  while (b != 0L) {
    tmp <- b
    b <- a %% b
    a <- tmp
  }
  a
}

integer_lcm <- function(x) {
  x <- abs(as.integer(x))
  x <- x[is.finite(x) & !is.na(x) & x > 0L]
  if (length(x) == 0L) {
    return(1L)
  }
  as.integer(Reduce(function(a, b) abs(a * b) / gcd_int(a, b), x))
}

esprit_grid_multiplier_from_lines <- function(lines) {
  lines <- normalize_config_lines(lines)
  a <- round(abs(lines$a))
  a[a < 1L] <- 1L
  integer_lcm(a)
}

cssa_rowrow_auto_rank_from_lines <- function(lines) {
  max(1L, nrow(normalize_config_lines(lines)))
}

cssa_colrow_safe_rank_from_lines <- function(lines, n_row, n_col) {
  lines <- normalize_config_lines(lines)
  n_row <- max(1L, as.integer(n_row)[1L])
  n_col <- max(1L, as.integer(n_col)[1L])
  L <- (n_row + 1L) %/% 2L
  K <- n_row - L + 1L

  ranks <- vapply(seq_len(nrow(lines)), function(i) {
    cur_a <- lines$a[i]
    cur_b <- lines$b[i]
    alpha <- max(1, abs(cur_a))

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
  }, integer(1L))

  max(1L, as.integer(sum(ranks)))
}

parse_component_pair <- function(x) {
  if (length(x) == 0L || all(is.na(x))) {
    return(c(row.row = NA_integer_, col.row = NA_integer_))
  }
  parts <- unlist(strsplit(paste(as.character(x), collapse = ","), "[,;[:space:]]+"))
  parts <- parts[nzchar(parts)]
  values <- suppressWarnings(as.integer(parts))
  values <- values[is.finite(values) & !is.na(values) & values >= 1L]
  if (length(values) == 0L) {
    return(c(row.row = NA_integer_, col.row = NA_integer_))
  }
  if (length(values) == 1L) {
    values <- rep(values, 2L)
  }
  c(row.row = values[1L], col.row = values[2L])
}

auto_rank_pair <- function(lines, n_row, n_col) {
  c(
    row.row = min(
      cssa_rowrow_auto_rank_from_lines(lines),
      cssa_component_limit(n_col)
    ),
    col.row = min(
      cssa_colrow_safe_rank_from_lines(lines, n_row, n_col),
      cssa_component_limit(n_row)
    )
  )
}

resolve_method_rank <- function(method, lines, n_row, n_col, cssa_components = NA) {
  auto <- auto_rank_pair(lines, n_row, n_col)
  manual <- parse_component_pair(cssa_components)
  key <- if (grepl("col\\.row", method)) "col.row" else "row.row"
  limit <- if (identical(key, "col.row")) {
    cssa_component_limit(n_row)
  } else {
    cssa_component_limit(n_col)
  }

  if (!is.finite(manual[[key]]) || is.na(manual[[key]])) {
    value <- auto[[key]]
    source <- "auto"
  } else {
    value <- min(max(1L, manual[[key]]), limit)
    source <- if (value < manual[[key]]) "manual_clamped" else "manual"
  }

  list(value = value, source = source, auto = auto[[key]])
}

threshold_binary <- function(denoised,
                             noisy,
                             threshold = 0.1,
                             use_noise_cut = TRUE,
                             m_coef = 1) {
  denoised <- as.matrix(denoised)
  noisy <- as.matrix(noisy)
  noise_sd <- stats::sd(as.numeric(noisy - denoised), na.rm = TRUE)
  threshold_value <- if (isTRUE(use_noise_cut)) {
    coef <- as.numeric(m_coef)[1L]
    if (!is.finite(coef)) {
      coef <- 1
    }
    coef * noise_sd
  } else {
    thr <- as.numeric(threshold)[1L]
    if (!is.finite(thr)) {
      thr <- 0
    }
    thr
  }

  active <- if (threshold_value <= 0) {
    denoised > 0
  } else {
    denoised >= threshold_value
  }

  list(
    denoised = denoised,
    binary = ifelse(active, 1, 0),
    threshold_value = threshold_value,
    noise_sd = noise_sd
  )
}

cssa_denoise <- function(m,
                         method = c("row.row", "col.row", "row.col", "col.col"),
                         num_components = 1L,
                         clip = TRUE) {
  method <- match.arg(method)
  num_components <- max(1L, as.integer(num_components)[1L])
  z <- dft(as.matrix(m))
  cleaned <- switch(
    method,
    "row.row" = z |> cssa.row(num.line = num_components) |> idft.row() |> Re(),
    "col.row" = z |> cssa.col(num.line = num_components) |> idft.row() |> Re(),
    "row.col" = z |> cssa.row(num.line = num_components) |> idft.col() |> Re(),
    "col.col" = z |> cssa.col(num.line = num_components) |> idft.col() |> Re()
  )
  if (isTRUE(clip)) clip01(cleaned) else cleaned
}

cssa_method_binary <- function(m,
                               method = "row.row",
                               noisy = m,
                               threshold = 0.1,
                               use_noise_cut = TRUE,
                               m_coef = 1,
                               num_components = 1L) {
  threshold_binary(
    denoised = cssa_denoise(
      m,
      method = method,
      num_components = num_components,
      clip = TRUE
    ),
    noisy = noisy,
    threshold = threshold,
    use_noise_cut = use_noise_cut,
    m_coef = m_coef
  )
}

row_max_binary <- function(m, k = 1L) {
  m <- as.matrix(m)
  k <- max(1L, min(as.integer(k)[1L], ncol(m)))
  binary <- matrix(0, nrow = nrow(m), ncol = ncol(m))
  for (i in seq_len(nrow(m))) {
    ind <- order(m[i, ], decreasing = TRUE)[seq_len(k)]
    binary[i, ind] <- m[i, ind]
  }
  list(
    denoised = binary,
    binary = binary,
    threshold_value = NA_real_,
    noise_sd = NA_real_
  )
}

pad_replicate <- function(x, pad) {
  nr <- nrow(x)
  nc <- ncol(x)
  row_idx <- pmin(pmax(seq_len(nr + 2L * pad) - pad, 1L), nr)
  col_idx <- pmin(pmax(seq_len(nc + 2L * pad) - pad, 1L), nc)
  x[row_idx, col_idx]
}

median_denoise <- function(m, n = 3L) {
  x <- clip01(m)
  n <- max(1L, as.integer(n)[1L])
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

line_preserving_median_denoise <- function(m,
                                           n = 3L,
                                           seed_sigma = 1.5,
                                           support_sigma = 1.0,
                                           max_abs_step = 2L,
                                           clip = FALSE) {
  m <- as.matrix(m)
  nr <- nrow(m)
  nc <- ncol(m)
  base <- median_denoise(m, n = n)
  center <- stats::median(as.numeric(m), na.rm = TRUE)
  sigma_hat <- stats::mad(as.numeric(m), center = center, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(sigma_hat) || sigma_hat <= 0) {
    sigma_hat <- stats::sd(as.numeric(m), na.rm = TRUE)
  }
  if (!is.finite(sigma_hat) || sigma_hat <= 0) {
    sigma_hat <- 1e-12
  }

  seed_level <- center + seed_sigma * sigma_hat
  support_level <- center + support_sigma * sigma_hat
  seed <- m > seed_level
  support <- m > support_level

  make_directions <- function(max_abs_step) {
    dirs <- list()
    keys <- character()
    for (dr in seq.int(-max_abs_step, max_abs_step)) {
      for (dc in seq.int(-max_abs_step, max_abs_step)) {
        if (dr == 0L && dc == 0L) next
        g <- gcd_int(dr, dc)
        dr0 <- as.integer(dr / g)
        dc0 <- as.integer(dc / g)
        if (dc0 < 0L || (dc0 == 0L && dr0 < 0L)) {
          dr0 <- -dr0
          dc0 <- -dc0
        }
        key <- paste(dr0, dc0, sep = ":")
        if (!key %in% keys) {
          keys <- c(keys, key)
          dirs[[length(dirs) + 1L]] <- c(dr = dr0, dc = dc0)
        }
      }
    }
    dirs
  }

  shift_matrix <- function(x, dr, dc, fill = FALSE) {
    out <- matrix(fill, nrow = nrow(x), ncol = ncol(x))
    dst_r <- seq_len(nrow(x))
    dst_c <- seq_len(ncol(x))
    src_r <- dst_r + dr
    src_c <- dst_c + dc
    ok_r <- src_r >= 1L & src_r <= nrow(x)
    ok_c <- src_c >= 1L & src_c <= ncol(x)
    if (any(ok_r) && any(ok_c)) {
      out[ok_r, ok_c] <- x[src_r[ok_r], src_c[ok_c], drop = FALSE]
    }
    out
  }

  dirs <- make_directions(as.integer(max_abs_step)[1L])
  line_supported <- matrix(FALSE, nrow = nr, ncol = nc)
  for (d in dirs) {
    before <- shift_matrix(support, -d[["dr"]], -d[["dc"]], fill = FALSE)
    after <- shift_matrix(support, d[["dr"]], d[["dc"]], fill = FALSE)
    line_supported <- line_supported | (seed & before & after)
  }

  gap_value <- matrix(NA_real_, nrow = nr, ncol = nc)
  gap_supported <- matrix(FALSE, nrow = nr, ncol = nc)
  for (d in dirs) {
    before_support <- shift_matrix(support, -d[["dr"]], -d[["dc"]], fill = FALSE)
    after_support <- shift_matrix(support, d[["dr"]], d[["dc"]], fill = FALSE)
    before_value <- shift_matrix(m, -d[["dr"]], -d[["dc"]], fill = NA_real_)
    after_value <- shift_matrix(m, d[["dr"]], d[["dc"]], fill = NA_real_)
    cur_supported <- before_support & after_support
    cur_value <- apply(
      cbind(as.numeric(before_value), as.numeric(m), as.numeric(after_value)),
      1L,
      stats::median,
      na.rm = TRUE
    )
    cur_value <- matrix(cur_value, nrow = nr, ncol = nc)
    take <- cur_supported & (is.na(gap_value) | cur_value > gap_value)
    gap_value[take] <- cur_value[take]
    gap_supported <- gap_supported | cur_supported
  }

  out <- base
  out[line_supported] <- pmax(out[line_supported], m[line_supported])
  fill <- gap_supported & is.finite(gap_value) & gap_value > seed_level
  out[fill] <- pmax(out[fill], gap_value[fill])
  if (isTRUE(clip)) {
    out[out < 0] <- 0
    out[out > 1] <- 1
  }
  out
}

wiener_denoise <- function(m, ksize = 5L) {
  x <- clip01(m)
  ksize <- max(1L, as.integer(ksize)[1L])
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
  clip01(local_mean + gain * (x - local_mean))
}

esprit_project_cssa_corrected <- function(x,
                                          L = (length(x) + 1L) %/% 2L,
                                          num_components = 1L,
                                          frequency_grid_multiplier = 1) {
  x <- as.vector(x)
  N <- length(x)
  num_components <- max(1L, as.integer(num_components)[1L])
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

    rec <- Rssa::reconstruct(fit, groups = list(signal = seq_len(num_components)))
    x_signal <- as.vector(rec$signal)
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

cssa_esprit_denoise <- function(m,
                               method = c("row.row", "col.row"),
                               num_components = 1L,
                               frequency_grid_multiplier = 1,
                               clip = TRUE) {
  method <- match.arg(method)
  z <- dft(as.matrix(m))
  z_hat <- switch(
    method,
    "row.row" = {
      L <- (ncol(z) + 1L) %/% 2L
      t(vapply(seq_len(nrow(z)), function(i) {
        as.complex(esprit_project_cssa_corrected(
          z[i, ],
          L = L,
          num_components = num_components,
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
          num_components = num_components,
          frequency_grid_multiplier = frequency_grid_multiplier
        ))
      }, complex(nrow(z)))
    }
  )
  out <- Re(idft.row(z_hat))
  if (isTRUE(clip)) clip01(out) else out
}

line_parameter_error <- function(true_line, pred_line) {
  theta_error <- function(x, y) atan2(sin(x - y), cos(x - y))^2
  rho <- true_line[["rho"]]
  theta <- true_line[["theta"]]
  pred_rho <- pred_line[["rho"]]
  pred_theta <- pred_line[["theta"]]

  direct <- c(
    rho_mse = unname((rho - pred_rho)^2),
    theta_mse = unname(theta_error(theta, pred_theta))
  )
  flipped <- c(
    rho_mse = unname((rho + pred_rho)^2),
    theta_mse = unname(theta_error(theta, pred_theta + pi))
  )
  if (sum(direct) <= sum(flipped)) direct else flipped
}

parameter_line_errors_by_line <- function(true_lines, pred_lines) {
  true_lines <- as.matrix(true_lines)
  pred_lines <- as.matrix(pred_lines)
  if (nrow(true_lines) != nrow(pred_lines)) {
    stop("true_lines and pred_lines must have the same number of rows")
  }

  pair_error <- function(true_line, pred_line) {
    line_parameter_error(
      c(rho = unname(true_line[1L]), theta = unname(true_line[2L])),
      c(rho = unname(pred_line[1L]), theta = unname(pred_line[2L]))
    )
  }

  n_lines <- nrow(true_lines)
  if (n_lines == 1L) {
    err <- pair_error(true_lines[1L, ], pred_lines[1L, ])
    return(data.frame(
      line_index = 1L,
      pred_index = 1L,
      true_rho = true_lines[1L, "rho"],
      true_theta = true_lines[1L, "theta"],
      pred_rho = pred_lines[1L, "rho"],
      pred_theta = pred_lines[1L, "theta"],
      rho_mse = unname(err["rho_mse"]),
      theta_mse = unname(err["theta_mse"]),
      row.names = NULL
    ))
  }

  all_permutations <- function(v) {
    if (length(v) == 1L) {
      return(matrix(v, nrow = 1L))
    }
    do.call(rbind, lapply(seq_along(v), function(i) {
      cbind(v[i], all_permutations(v[-i]))
    }))
  }

  perms <- all_permutations(seq_len(n_lines))
  best_total <- Inf
  best_perm <- seq_len(n_lines)
  for (i in seq_len(nrow(perms))) {
    cur <- t(vapply(seq_len(n_lines), function(j) {
      pair_error(true_lines[j, ], pred_lines[perms[i, j], ])
    }, numeric(2L)))
    cur_total <- sum(colMeans(cur))
    if (cur_total < best_total) {
      best_total <- cur_total
      best_perm <- perms[i, ]
    }
  }

  out <- do.call(rbind, lapply(seq_len(n_lines), function(i) {
    err <- pair_error(true_lines[i, ], pred_lines[best_perm[i], ])
    data.frame(
      line_index = i,
      pred_index = best_perm[i],
      true_rho = true_lines[i, "rho"],
      true_theta = true_lines[i, "theta"],
      pred_rho = pred_lines[best_perm[i], "rho"],
      pred_theta = pred_lines[best_perm[i], "theta"],
      rho_mse = unname(err["rho_mse"]),
      theta_mse = unname(err["theta_mse"]),
      row.names = NULL
    )
  }))
  rownames(out) <- NULL
  out
}

find_k_max <- function(acc, k, qrho, qtheta, suppress = FALSE, window = 6L) {
  acc_copy <- acc
  maxima <- vector("list", k)
  for (i in seq_len(k)) {
    idx <- which(acc_copy == max(acc_copy), arr.ind = TRUE)[1L, ]
    maxima[[i]] <- c(rho = qrho[idx[1L]], theta = qtheta[idx[2L]])
    if (isTRUE(suppress)) {
      rows <- max(1L, idx[1L] - window):min(nrow(acc_copy), idx[1L] + window)
      cols <- max(1L, idx[2L] - window):min(ncol(acc_copy), idx[2L] + window)
      acc_copy[rows, cols] <- 0
    } else {
      acc_copy[idx[1L], idx[2L]] <- 0
    }
  }
  out <- do.call(rbind, maxima)
  colnames(out) <- c("rho", "theta")
  out
}

project_to_hough_grid <- function(true_line, rho_step, theta_step) {
  pred_line <- cbind(
    rho = round(true_line[, "rho"] / rho_step) * rho_step,
    theta = round(true_line[, "theta"] / theta_step) * theta_step
  )
  pred_line[, "theta"] <- pmin(pmax(pred_line[, "theta"], 0), pi)
  pred_line
}

format_metric_value <- function(x, metric) {
  if (is.na(x)) {
    return("")
  }
  if (metric == "mean_active_pixels") {
    out <- formatC(x, format = "f", digits = 1L)
    return(sub("\\.0$", "", out))
  }
  if (x == 0) {
    return("0")
  }
  abs_x <- abs(x)
  if (metric == "mean_theta_mse") {
    if (abs_x < 1e-3) {
      return(sprintf("%.2e", x))
    }
    return(sprintf("%.4f", x))
  }
  if (abs_x < 0.01) {
    return(sprintf("%.2e", x))
  }
  if (abs_x >= 100) {
    return(sprintf("%.1f", x))
  }
  sprintf("%.4f", x)
}

format_plain_num <- function(x, digits = 4L) {
  ifelse(
    abs(x) > 0 & abs(x) < 1e-3,
    sprintf("%.2e", x),
    format(round(x, digits), nsmall = digits, scientific = FALSE, trim = TRUE)
  )
}

write_csv_utf8 <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(file)
}

write_lines_utf8 <- function(lines, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, con = file, useBytes = TRUE)
  invisible(file)
}

latex_table_lines <- function(df, caption = NULL, label = NULL, resize = "\\textwidth") {
  n_cols <- ncol(df)
  header <- paste(colnames(df), collapse = " & ")
  align <- paste0("|", paste(rep("c", n_cols), collapse = "|"), "|")
  body <- unlist(lapply(seq_len(nrow(df)), function(i) {
    c(paste(df[i, ], collapse = " & "), "\\\\", "\\hline")
  }))
  lines <- c()
  if (!is.null(caption) || !is.null(label)) {
    lines <- c(lines, "\\begin{table}[H]", "\\centering")
    if (!is.null(caption)) lines <- c(lines, sprintf("\\caption{%s}", caption))
    if (!is.null(label)) lines <- c(lines, sprintf("\\label{%s}", label))
    if (!is.null(resize)) lines <- c(lines, sprintf("\\resizebox{%s}{!}{%%", resize))
  }
  lines <- c(
    lines,
    sprintf("\\begin{tabular}{%s}", align),
    "\\hline",
    header,
    "\\\\",
    "\\hline",
    body,
    "\\end{tabular}"
  )
  if (!is.null(caption) || !is.null(label)) {
    if (!is.null(resize)) lines <- c(lines, "}")
    lines <- c(lines, "\\end{table}")
  }
  lines
}

save_matrix_image <- function(m,
                              file,
                              title = "",
                              mode = c("bw", "blue_black", "blue_red", "accumulator"),
                              from_0_to_1 = FALSE,
                              width = 1200,
                              height = 1000,
                              res = 180,
                              xlab = "",
                              ylab = "") {
  mode <- match.arg(mode)
  mat <- as.matrix(m)
  if (is.complex(mat)) mat <- Re(mat)
  mat[!is.finite(mat)] <- 0
  if (isTRUE(from_0_to_1)) {
    mat[mat < 0] <- 0
    mat[mat > 1] <- 1
  }

  ext <- tolower(tools::file_ext(file))
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  if (ext == "pdf") {
    grDevices::pdf(file, width = width / res, height = height / res)
  } else {
    grDevices::png(file, width = width, height = height, res = res)
  }
  on.exit(grDevices::dev.off(), add = TRUE)

  if (mode == "blue_black") {
    max_abs <- max(abs(mat), na.rm = TRUE)
    if (!is.finite(max_abs) || max_abs == 0) {
      max_abs <- 1
    }
    pos <- pmin(pmax(mat / max_abs, 0), 1)
    neg <- pmin(pmax(-mat / max_abs, 0), 1)
    strength <- pmax(pos, neg)
    red <- 1 - strength
    green <- 1 - strength
    blue <- 1 - strength
    blue[neg > 0] <- 1
    colors <- grDevices::rgb(red, green, blue)
    plot_colors <- unique(as.vector(colors))
    color_index <- matrix(match(as.vector(colors), plot_colors), nrow = nrow(mat))

    graphics::par(mar = c(3, 3, 3, 1))
    graphics::image(
      x = seq_len(ncol(mat)),
      y = seq_len(nrow(mat)),
      z = t(color_index),
      col = plot_colors,
      breaks = seq(0.5, length(plot_colors) + 0.5, by = 1),
      xlab = xlab,
      ylab = ylab,
      main = title,
      useRaster = TRUE
    )
    graphics::box()
    return(invisible(file))
  }

  if (mode == "bw") {
    pal <- grDevices::colorRampPalette(c("white", "black"))(256L)
    rng <- if (isTRUE(from_0_to_1)) c(0, 1) else range(mat)
  } else if (mode == "blue_red") {
    pal <- grDevices::colorRampPalette(c("navy", "white", "firebrick"))(256L)
    rng <- range(mat)
  } else {
    pal <- grDevices::colorRampPalette(c("white", "yellow", "red"))(256L)
    mat <- pmax(mat, 0)
    rng <- c(0, max(mat))
  }
  if (!is.finite(rng[1L]) || !is.finite(rng[2L]) || rng[1L] == rng[2L]) {
    rng <- rng + c(-0.5, 0.5)
  }

  graphics::par(mar = c(3, 3, 3, 1))
  graphics::image(
    x = seq_len(ncol(mat)),
    y = seq_len(nrow(mat)),
    z = t(mat),
    col = pal,
    breaks = seq(rng[1L], rng[2L], length.out = length(pal) + 1L),
    xlab = xlab,
    ylab = ylab,
    main = title,
    useRaster = TRUE
  )
  graphics::box()
  invisible(file)
}
