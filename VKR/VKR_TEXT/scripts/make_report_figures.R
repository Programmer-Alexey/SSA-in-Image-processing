script_dir_for_source <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/")))
  }

  source_frames <- Filter(
    function(frame) !is.null(frame$ofile),
    sys.frames()
  )
  if (length(source_frames) > 0L) {
    return(dirname(normalizePath(source_frames[[1L]]$ofile, winslash = "/")))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    editor_path <- tryCatch(
      rstudioapi::getSourceEditorContext()$path,
      error = function(...) ""
    )
    if (nzchar(editor_path)) {
      return(dirname(normalizePath(editor_path, winslash = "/")))
    }
  }

  candidates <- c(
    getwd(),
    file.path(getwd(), "scripts"),
    file.path(getwd(), "VKR", "VKR_TEXT", "scripts")
  )
  common_path <- file.path(candidates, "00_common.R")
  match <- which(file.exists(common_path))[1L]
  if (!is.na(match)) {
    return(normalizePath(candidates[match], winslash = "/"))
  }

  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

source(file.path(script_dir_for_source(), "00_common.R"))

ensure_report_dirs()
load_project_code(load_hough = TRUE)

draw_rt_params <- function(file, a = -2, b = 4) {
  xlim <- c(-1, 5)
  ylim <- c(-1, 7)

  margin_in <- c(bottom = 0.65, left = 0.65, top = 0.12, right = 0.12)
  plot_unit <- 0.55
  plot_width <- diff(xlim) * plot_unit
  plot_height <- diff(ylim) * plot_unit

  grDevices::pdf(
    file,
    width = plot_width + margin_in["left"] + margin_in["right"],
    height = plot_height + margin_in["bottom"] + margin_in["top"]
  )
  on.exit(grDevices::dev.off(), add = TRUE)

  line_candidates <- rbind(
    c(xlim[1L], a * xlim[1L] + b),
    c(xlim[2L], a * xlim[2L] + b),
    c((ylim[1L] - b) / a, ylim[1L]),
    c((ylim[2L] - b) / a, ylim[2L])
  )
  inside <-
    line_candidates[, 1L] >= xlim[1L] &
    line_candidates[, 1L] <= xlim[2L] &
    line_candidates[, 2L] >= ylim[1L] &
    line_candidates[, 2L] <= ylim[2L]
  line_points <- unique(round(line_candidates[inside, , drop = FALSE], 10), MARGIN = 1)
  if (nrow(line_points) < 2L) {
    stop("Line y = a x + b does not intersect the plotting window in two distinct points.")
  }

  theta <- atan2(1, -a)
  if (theta < 0) {
    theta <- theta + pi
  }
  rho <- b * sin(theta)
  foot <- c(rho * cos(theta), rho * sin(theta))
  line_label <- sprintf("y=%sx%+g", formatC(a, format = "fg", digits = 3), b)
  line_angle <- atan(a) * 180 / pi
  dash_step <- 0.33
  dash_gap <- 0.20
  dash_unit <- foot / sqrt(sum(foot^2))
  dash_start <- seq(0, rho, by = dash_step + dash_gap)
  dash_end <- pmin(dash_start + dash_step, rho)

  graphics::par(mai = margin_in, xaxs = "i", yaxs = "i", lend = "butt")
  graphics::plot(
    NA,
    xlim = xlim,
    ylim = ylim,
    xlab = "x",
    ylab = "",
    asp = 1,
    axes = FALSE,
    main = ""
  )
  graphics::axis(1, at = seq(xlim[1L], xlim[2L], by = 1))
  graphics::axis(2, at = seq(ylim[1L], ylim[2L], by = 2), las = 1)
  graphics::abline(h = 0, v = 0, col = "grey35", lwd = 0.9)

  graphics::segments(line_points[1L, 1L], line_points[1L, 2L], line_points[2L, 1L], line_points[2L, 2L], lwd = 2.2, col = "black")
  graphics::segments(
    dash_start * dash_unit[1L],
    dash_start * dash_unit[2L],
    dash_end * dash_unit[1L],
    dash_end * dash_unit[2L],
    lwd = 1.8,
    col = "grey25"
  )
  graphics::text(
    mean(line_points[, 1L]) + 0.28,
    mean(line_points[, 2L]),
    labels = line_label,
    srt = line_angle,
    adj = c(0.5, -0.2),
    cex = 1.05
  )
  graphics::text(foot[1] / 2, foot[2] / 2 + 0.35, expression(rho), cex = 1.0)

  ang <- seq(0, theta, length.out = 50)
  graphics::lines(1.15 * cos(ang), 1.15 * sin(ang), lwd = 1.2)
  graphics::text(1.45 * cos(theta / 2), 1.45 * sin(theta / 2), expression(theta), pos = 4, cex = 0.95)
  graphics::box()
}
#draw_rt_params(file.path(chapter1_image_dir, "rt_params.pdf"))

save_profile_plot <- function(est,
                              signal_index,
                              file,
                              title = "",
                              xlim = NULL) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(file, width = 1600, height = 850, res = 200)
  on.exit(grDevices::dev.off(), add = TRUE)

  est <- as.numeric(est)
  x <- seq_along(est)
  max_index <- which.max(est)
  if (is.null(xlim)) {
    left <- max(1L, min(signal_index, max_index) - 18L)
    right <- min(length(est), max(signal_index, max_index) + 18L)
    xlim <- c(left, right)
  }
  show <- x >= xlim[1L] & x <= xlim[2L]
  ylim <- range(c(est[show], 0, est[signal_index], est[max_index]), finite = TRUE)
  span <- diff(ylim)
  if (!is.finite(span) || span == 0) {
    span <- 1
  }
  ylim <- ylim + c(-0.10, 0.14) * span

  graphics::par(mar = c(4.1, 4.2, 2.4, 0.8), xaxs = "i", yaxs = "i")
  graphics::plot(
    x[show],
    est[show],
    type = "h",
    lwd = 1.7,
    col = "grey35",
    xlab = "Позиция",
    ylab = "Восстановленное значение",
    main = title,
    xlim = xlim,
    ylim = ylim
  )
  graphics::abline(h = 0, col = "grey70", lwd = 0.5)
  graphics::abline(v = signal_index, col = "firebrick", lty = 2, lwd = 1.4)
  graphics::points(max_index, est[max_index], pch = 19, col = "black", cex = 1.15)
  graphics::legend(
    "topright",
    legend = c("Истинный индекс", "Максимум"),
    col = c("firebrick", "black"),
    lty = c(2, NA),
    lwd = c(1.4, NA),
    pch = c(NA, 19),
    pt.bg = c(NA, "black"),
    pt.cex = c(1.0, 1.0),
    bty = "n"
  )
}

save_accumulator_plot <- function(ht, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(file, width = 1200, height = 950, res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)

  acc <- as.matrix(ht$accumulator)
  acc[!is.finite(acc)] <- 0
  display <- pmax(acc, 0)
  max_display <- max(display, na.rm = TRUE)
  pal <- grDevices::colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(256L)

  graphics::par(mar = c(4.0, 4.4, 0.8, 0.8), xaxs = "i", yaxs = "i")
  graphics::plot(
    NA,
    xlim = range(ht$theta),
    ylim = range(ht$rho),
    xlab = expression(theta),
    ylab = expression(rho),
    axes = FALSE
  )
  graphics::axis(1)
  graphics::axis(2, las = 1)
  pos <- which(display > 0, arr.ind = TRUE)
  if (nrow(pos) > 0L) {
    values <- display[pos]
    ord <- order(values)
    pos <- pos[ord, , drop = FALSE]
    values <- values[ord]
    scaled <- sqrt(values / max_display)
    color_index <- pmax(1L, pmin(256L, as.integer(ceiling(scaled * 255)) + 1L))
    graphics::points(
      ht$theta[pos[, 2L]],
      ht$rho[pos[, 1L]],
      pch = 15,
      cex = 0.08 + 0.22 * scaled,
      col = pal[color_index]
    )
    high <- values >= max(2, 0.35 * max_display)
    if (any(high)) {
      graphics::points(
        ht$theta[pos[high, 2L]],
        ht$rho[pos[high, 1L]],
        pch = 15,
        cex = 0.60,
        col = "red3"
      )
    }
  }
  graphics::box()
}

one_row_estimate <- function(signal_index, noise_row, ncomp = 1L) {
  clean <- matrix(0, nrow = 1L, ncol = ncol(noise_row))
  clean[1, signal_index] <- 1
  noisy <- clean + noise_row
  est <- noisy |>
    dft() |>
    cssa.row(num.line = ncomp) |>
    idft.row() |>
    Re()
  list(signal_index = signal_index, est = as.numeric(est[1, ]))
}

find_shift_case <- function(noise_row, search_idx) {
  best <- NULL
  best_dist <- Inf
  best_value <- -Inf
  for (idx in search_idx) {
    cur <- one_row_estimate(idx, noise_row, ncomp = 1L)
    max_index <- which.max(cur$est)
    dist <- abs(max_index - idx)
    if (dist > 0 && (dist < best_dist || (dist == best_dist && cur$est[max_index] > best_value))) {
      best <- cur
      best_dist <- dist
      best_value <- cur$est[max_index]
    }
  }
  if (is.null(best)) {
    best <- one_row_estimate(search_idx[1L], noise_row, ncomp = 1L)
  }
  best
}

find_smear_case <- function(noise_row, search_idx) {
  best <- NULL
  best_score <- -Inf
  for (idx in search_idx) {
    cur <- one_row_estimate(idx, noise_row, ncomp = 1L)
    pos <- pmax(cur$est, 0)
    window <- pmax(1, idx - 4):pmin(length(pos), idx + 4)
    score <- sum(pos[window]) - max(pos)
    if (which.max(cur$est) == idx && score > best_score) {
      best <- cur
      best_score <- score
    }
  }
  if (is.null(best)) {
    best <- one_row_estimate(search_idx[1L], noise_row, ncomp = 1L)
  }
  best
}

find_first_shift_seed <- function(start_seed = 1L,
                                  max_seed = 10000L,
                                  search_idx = 40:360,
                                  sigma_noise = 0.2) {
  for (seed in seq.int(start_seed, max_seed)) {
    set.seed(seed)
    noise_row <- add.noise(matrix(0, nrow = 1L, ncol = 400L), sigma = sigma_noise)
    cur <- find_shift_case(noise_row, search_idx)
    max_index <- which.max(cur$est)
    if (max_index != cur$signal_index) {
      return(list(
        seed = seed,
        noise_row = noise_row,
        shift_case = cur,
        max_index = max_index,
        shift = abs(max_index - cur$signal_index)
      ))
    }
  }
  stop("No shifted maximum found in the requested seed range.")
}

make_shift_profile_case <- function(seed,
                                    search_idx = 40:360,
                                    sigma_noise = 0.2,
                                    low_noise_scale = 0.25) {
  set.seed(seed)
  noise_row <- add.noise(matrix(0, nrow = 1L, ncol = 400L), sigma = sigma_noise)
  shift_case <- find_shift_case(noise_row, search_idx)
  shift_case_low <- one_row_estimate(shift_case$signal_index, noise_row * low_noise_scale, ncomp = 1L)
  list(seed = seed, noise_row = noise_row, shift_case = shift_case, shift_case_low = shift_case_low)
}

make_smear_profile_case <- function(seed,
                                    search_idx = 40:360,
                                    sigma_noise = 0.2) {
  set.seed(seed)
  noise_row <- add.noise(matrix(0, nrow = 1L, ncol = 400L), sigma = sigma_noise)
  list(seed = seed, noise_row = noise_row, smear_case = find_smear_case(noise_row, search_idx))
}

draw_singular_values <- function(sigma_noise, file) {
  set.seed(1111)
  n <- 100L
  clean <- matrix(0, nrow = 1L, ncol = n)
  clean[1, 10] <- 1
  noisy <- add.noise(clean, sigma = sigma_noise)
  clean_fit <- Rssa::ssa(dft(clean)[1, ], kind = "cssa", svd.method = "svd")
  noisy_fit <- Rssa::ssa(dft(noisy)[1, ], kind = "cssa", svd.method = "svd")

  grDevices::png(file, width = 1600, height = 1200, res = 200)
  on.exit(grDevices::dev.off(), add = TRUE)

  k <- 1:12
  noisy_sigma <- noisy_fit$sigma[k]
  clean_sigma <- clean_fit$sigma[1L]
  y_max <- max(c(clean_sigma, noisy_sigma), na.rm = TRUE)

  graphics::par(mar = c(4.1, 4.2, 2.3, 0.8), xaxs = "i", yaxs = "i")
  graphics::plot(
    k,
    noisy_sigma,
    type = "b",
    pch = 19,
    lty = 1,
    col = "firebrick",
    xlab = "Номер компоненты",
    ylab = "Сингулярное значение",
    main = paste0("sigma = ", sigma_noise),
    xlim = c(0, max(k) + 0.5),
    ylim = c(0, 1.08 * y_max),
    xaxt = "n"
  )
  graphics::axis(1, at = 0:max(k))
  graphics::points(1, clean_sigma, pch = 19, col = "black", cex = 1.15)
  graphics::legend(
    "topright",
    legend = c("Чистая строка", "Зашумленная строка"),
    col = c("black", "firebrick"),
    lty = 1,
    pch = 19,
    bty = "n"
  )
}

threshold_by_residual <- function(processed, noisy, multiplier = 1) {
  sigma_hat <- stats::sd(as.numeric(noisy - processed), na.rm = TRUE)
  ifelse(processed >= multiplier * sigma_hat, 1, 0)
}

save_preprocessing_method <- function(noisy,
                                      method,
                                      file,
                                      lines,
                                      threshold = FALSE,
                                      for_demo_median = FALSE) {
  row_rank <- resolve_method_rank("row.row", lines, nrow(noisy), ncol(noisy))$value
  col_rank <- resolve_method_rank("col.row", lines, nrow(noisy), ncol(noisy))$value
  max_k <- nrow(normalize_config_lines(lines))

  processed <- switch(
    method,
    row.row = cssa_denoise(noisy, method = "row.row", num_components = row_rank, clip = FALSE),
    col.row = cssa_denoise(noisy, method = "col.row", num_components = col_rank, clip = FALSE),
    max.row = row_max_binary(noisy, k = max_k)$denoised,
    median = {
      if (isTRUE(for_demo_median)) {
        line_preserving_median_denoise(noisy, n = 3L, clip = FALSE)
      } else {
        median_denoise(noisy, n = 3L)
      }
    },
    wiener = wiener_denoise(noisy, ksize = 5L),
    stop("Unknown figure method: ", method)
  )

  shown <- if (isTRUE(threshold) && method %in% c("row.row", "col.row")) {
    threshold_by_residual(processed, noisy, multiplier = 1)
  } else {
    processed
  }

  save_matrix_image(
    shown,
    file = file,
    title = "",
    mode = if (isTRUE(threshold) || method == "max.row") "bw" else "blue_black",
    from_0_to_1 = isTRUE(threshold) || method == "max.row",
    width = 1000,
    height = 1000,
    res = 180
  )
}

cat("Generating chapter 1 figures...\n")
draw_rt_params(file.path(chapter1_image_dir, "rt_params.pdf"))

one_line <- data.frame(a = 1, b = -10, intensity = 1)
one_clean <- make_line_image(100, 100, one_line, line_method = "default")
save_matrix_image(
  one_clean,
  file.path(chapter1_image_dir, "chapter1_one_line_clean.png"),
  mode = "bw",
  from_0_to_1 = TRUE,
  width = 900,
  height = 900,
  res = 180
)

two_lines <- data.frame(a = c(2, -1), b = c(-1, 101), intensity = c(0.8, 0.5))
two_clean <- make_line_image(100, 100, two_lines, line_method = "default")
save_matrix_image(
  two_clean,
  file.path(chapter1_image_dir, "chapter1_two_lines_clean.png"),
  mode = "bw",
  from_0_to_1 = TRUE,
  width = 1000,
  height = 900,
  res = 180,
  xlab = "x",
  ylab = "y"
)
two_ht <- make_accumulator(
  two_clean,
  detector = function(m) ifelse(m > 0, 1, 0),
  rho_step = 0.02,
  theta_step = 0.002
)
save_accumulator_plot(
  two_ht,
  file.path(chapter1_image_dir, "chapter1_two_lines_accumulator.png")
)

set.seed(322)
demo_lines <- data.frame(a = 1, b = -10, intensity = 0.8)
demo_clean <- make_line_image(100, 100, demo_lines, line_method = "default")
demo_noisy <- add.noise(demo_clean, sigma = 0.05)

save_preprocessing_method(demo_noisy, "row.row", file.path(chapter1_image_dir, "chapter1_cssa_row_row_processed.png"), demo_lines)
save_preprocessing_method(demo_noisy, "col.row", file.path(chapter1_image_dir, "chapter1_cssa_col_row_processed.png"), demo_lines)
save_preprocessing_method(demo_noisy, "max.row", file.path(chapter1_image_dir, "chapter1_max_row_processed.png"), demo_lines)
save_preprocessing_method(demo_noisy, "median", file.path(chapter1_image_dir, "chapter1_median_processed.png"), demo_lines, for_demo_median = TRUE)
save_preprocessing_method(demo_noisy, "wiener", file.path(chapter1_image_dir, "chapter1_wiener_processed.png"), demo_lines)

cat("Generating chapter 2 figures...\n")
draw_singular_values(0.2, file.path(chapter2_image_dir, "singular-values-02.png"))
draw_singular_values(0.05, file.path(chapter2_image_dir, "singular-values-005.png"))

profile_seed <- 1211L
shift_profile <- make_shift_profile_case(profile_seed)
shift_case <- shift_profile$shift_case
shift_case_low <- shift_profile$shift_case_low

smear_profile <- make_smear_profile_case(profile_seed)
smear_case <- smear_profile$smear_case

save_profile_plot(
  shift_case$est,
  shift_case$signal_index,
  file.path(chapter2_image_dir, "shift_profiling.png"),
  "Смещение максимума"
)
save_profile_plot(
  shift_case_low$est,
  shift_case_low$signal_index,
  file.path(chapter2_image_dir, "shift_profiling_sigma_005.png"),
  "Смещение максимума"
)
save_profile_plot(
  smear_case$est,
  smear_case$signal_index,
  file.path(chapter2_image_dir, "smear_profiling.png"),
  "Растекание"
)

set.seed(2027)
threshold_clean <- make_line_image(100, 100, two_lines, line_method = "default")
threshold_noisy <- add.noise(threshold_clean, sigma = 0.2)
save_preprocessing_method(threshold_noisy, "row.row", file.path(chapter2_image_dir, "row_row_before_threshold.png"), two_lines)
save_preprocessing_method(threshold_noisy, "col.row", file.path(chapter2_image_dir, "col_row_before_threshold.png"), two_lines)
save_preprocessing_method(threshold_noisy, "row.row", file.path(chapter2_image_dir, "row_row_after_threshold.png"), two_lines, threshold = TRUE)
save_preprocessing_method(threshold_noisy, "col.row", file.path(chapter2_image_dir, "col_row_after_threshold.png"), two_lines, threshold = TRUE)

cat("Done. Figures are in: ", file.path(report_root, "assets", "images"), "\n", sep = "")
