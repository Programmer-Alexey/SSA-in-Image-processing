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

source(file.path(repo_root, "ssa-based methods", "cssa-transform.r"))

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

cssa_row_row <- function(m, num_of_lines = 1L) {
  m |>
    dft() |>
    cssa.row(num.line = num_of_lines) |>
    idft.row() |>
    Re() |>
    clip01()
}

binary_active <- function(m) {
  ifelse(as.matrix(m) > 0, 1, 0)
}

row_components <- function(row_vec, ncomp = 3L) {
  z <- fft(row_vec)
  s <- ssa(z, kind = "cssa", svd.method = "svd")
  comps <- vector("list", ncomp)

  for (j in seq_len(ncomp)) {
    r <- reconstruct(s, groups = list(Seasonality = j))
    comps[[j]] <- Re(fft(r$Seasonality, inverse = TRUE) / length(r$Seasonality))
  }

  comps
}

build_one_row_cssa_example <- function(n = 100L, signal_ind = 10L, sigma_noise = 0.2) {
  clean_matrix <- matrix(0, nrow = 1, ncol = n)
  clean_matrix[1, signal_ind] <- 1
  noisy_matrix <- add.noise(clean_matrix, sigma = sigma_noise)

  clean_dft <- dft(clean_matrix)[1, ]
  noisy_dft <- dft(noisy_matrix)[1, ]

  list(
    clean_fit = ssa(clean_dft, kind = "cssa", svd.method = "svd"),
    noisy_fit = ssa(noisy_dft, kind = "cssa", svd.method = "svd"),
    sigma_noise = sigma_noise
  )
}

draw_singular_values_compare <- function() {
  example_high <- build_one_row_cssa_example(n = 100L, signal_ind = 10L, sigma_noise = 0.2)
  example_low <- build_one_row_cssa_example(n = 100L, signal_ind = 10L, sigma_noise = low_noise_sigma)

  spectra <- list(example_high, example_low)
  y_values <- unlist(lapply(
    spectra,
    function(ex) c(ex$clean_fit$sigma[1:12], ex$noisy_fit$sigma[1:12])
  ))
  y_values <- y_values[is.finite(y_values) & y_values > 0]
  global_ylim <- range(y_values)

  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

  for (ex in spectra) {
    sigma_df <- data.frame(
      component = 1:12,
      clean = round(ex$clean_fit$sigma[1:12], 6),
      noisy = round(ex$noisy_fit$sigma[1:12], 6)
    )

    matplot(
      sigma_df$component,
      sigma_df[, c("clean", "noisy")],
      type = "b",
      pch = 19,
      lty = 1,
      log = "y",
      col = c("black", "firebrick"),
      xlab = "Номер компоненты",
      ylab = "Сингулярное значение",
      main = bquote(sigma == .(ex$sigma_noise)),
      ylim = global_ylim
    )

    legend(
      "topright",
      legend = c("Чистая строка", "Шумная строка"),
      col = c("black", "firebrick"),
      lty = 1,
      pch = 19,
      bty = "n"
    )
  }
}

collect_row_diagnostics <- function(noisy_matrix, signal_matrix, ncomp = 3L) {
  true_rows <- which(rowSums(signal_matrix) > 0)
  out <- vector("list", length(true_rows))

  for (idx in seq_along(true_rows)) {
    row_id <- true_rows[idx]
    true_col <- which(signal_matrix[row_id, ] > 0)[1]
    comps <- row_components(noisy_matrix[row_id, ], ncomp = ncomp)

    comp1 <- pmax(comps[[1]], 0)
    comp2 <- pmax(comps[[2]], 0)
    support <- max(1, true_col - 1):min(ncol(signal_matrix), true_col + 1)
    mass1 <- sum(comp1[support]) / max(sum(comp1), 1e-12)

    out[[idx]] <- data.frame(
      row_id = row_id,
      true_col = true_col,
      freq = (true_col - 1) / ncol(signal_matrix),
      argmax1 = which.max(comp1),
      argmax2 = which.max(comp2),
      mass1 = mass1,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out)
}

pick_example_rows <- function(diag_high, diag_low) {
  low_reduced <- diag_low[, c("row_id", "argmax1", "mass1")]
  names(low_reduced) <- c("row_id", "argmax1_low", "mass1_low")
  merged <- merge(diag_high, low_reduced, by = "row_id")

  shift_fixed_pool <- subset(
    merged,
    argmax1 != true_col & argmax1_low == true_col
  )
  smear_pool <- subset(
    merged,
    argmax1 == true_col & argmax1_low == true_col & mass1 < 0.55 & mass1_low > mass1
  )

  if (nrow(shift_fixed_pool) == 0 || nrow(smear_pool) == 0) {
    stop("Не удалось автоматически подобрать строки для обновленного рисунка 2.4.")
  }

  list(
    wrong = shift_fixed_pool$row_id[
      which.max(abs(shift_fixed_pool$argmax1 - shift_fixed_pool$true_col))
    ],
    smear = smear_pool$row_id[which.min(smear_pool$mass1)]
  )
}

plot_component_panel <- function(noisy_matrix, signal_matrix, row_id, main_title) {
  true_col <- which(signal_matrix[row_id, ] > 0)[1]
  comps <- row_components(noisy_matrix[row_id, ], ncomp = 1L)
  comp1 <- pmax(comps[[1]], 0)
  x <- seq_along(comp1)
  ylim_max <- max(comp1, na.rm = TRUE)

  plot(
    x,
    comp1,
    type = "h",
    col = "black",
    lwd = 2,
    xlab = "Номер столбца",
    ylab = "Интенсивность",
    ylim = c(0, max(ylim_max, 1e-6)),
    main = main_title
  )
  abline(v = true_col, col = "black", lty = 3)
}

plot_complex_parts <- function(z, main_text) {
  plot(
    Re(z),
    type = "l",
    col = "steelblue",
    lwd = 2,
    xlab = "Индекс",
    ylab = "Значение",
    main = main_text
  )
  lines(Im(z), col = "firebrick", lwd = 2, lty = 2)
  legend(
    "topright",
    legend = c("Re", "Im"),
    col = c("steelblue", "firebrick"),
    lty = c(1, 2),
    bty = "n"
  )
}

dir.create("images/chapter2", recursive = TRUE, showWarnings = FALSE)

low_noise_sigma <- 0.05

set.seed(1111)
signal_matrix <- make_line_image(100, 100, a = 2, b = -1)
noisy_02 <- add.noise(signal_matrix, sigma = 0.2)
processed_02 <- cssa_row_row(noisy_02, num_of_lines = 1L)
binary_02 <- binary_active(processed_02)

set.seed(1111)
noisy_low <- add.noise(signal_matrix, sigma = low_noise_sigma)
processed_low <- cssa_row_row(noisy_low, num_of_lines = 1L)
binary_low <- binary_active(processed_low)

diag_02 <- collect_row_diagnostics(noisy_02, signal_matrix, ncomp = 3L)
diag_low <- collect_row_diagnostics(noisy_low, signal_matrix, ncomp = 3L)
examples <- pick_example_rows(diag_02, diag_low)

set.seed(1211)
fixed_noise_matrix <- add.noise(matrix(0, nrow = 1, ncol = 1000), sigma = 0.2)

build_fixed_noise_case <- function(signal_ind, noise_matrix) {
  clean_matrix <- matrix(0, nrow = 1, ncol = ncol(noise_matrix))
  clean_matrix[1, signal_ind] <- 1
  noisy_matrix <- clean_matrix + noise_matrix
  noisy_dft_matrix <- dft(noisy_matrix)
  noisy_fit <- ssa(noisy_dft_matrix[1, ], kind = "cssa", svd.method = "svd")
  est <- noisy_dft_matrix |>
    cssa.row(num.line = 1L) |>
    idft.row() |>
    Re()

  list(
    signal_ind = signal_ind,
    fit = noisy_fit,
    est = as.numeric(est[1, ])
  )
}

bad_case <- NULL
for (candidate_signal_ind in 40:ncol(fixed_noise_matrix)) {
  candidate_case <- build_fixed_noise_case(candidate_signal_ind, fixed_noise_matrix)
  if (which.max(candidate_case$est) != candidate_case$signal_ind) {
    bad_case <- candidate_case
    break
  }
}

if (is.null(bad_case)) {
  stop("Не удалось построить плохой случай для однострочного эксперимента.")
}

noise_only_dft_matrix <- dft(fixed_noise_matrix)
noise_only_fit <- ssa(noise_only_dft_matrix[1, ], kind = "cssa", svd.method = "svd")
noise_only_est <- noise_only_dft_matrix |>
  cssa.row(num.line = 1L) |>
  idft.row() |>
  Re()
noise_only_est <- as.numeric(noise_only_est[1, ])

reconstruct_single_component <- function(fit, comp_idx) {
  rec <- Rssa::reconstruct(
    fit,
    groups = list(Seasonality = comp_idx)
  )$Seasonality

  row <- rec |>
    matrix(nrow = 1) |>
    idft.row() |>
    Re()

  as.numeric(row[1, ])
}

bad_case_row5 <- reconstruct_single_component(bad_case$fit, 5L)
noise_only_row5 <- reconstruct_single_component(noise_only_fit, 5L)

png("images/chapter2/row_row_binary_compare_noise.png", width = 1600, height = 800, res = 200)
plot.matrix(
  list(
    "sigma = 0.2" = binary_02,
    "sigma = low" = binary_low
  ),
  from.0.to.1 = TRUE,
  labels = c(expression(sigma == 0.2), bquote(sigma == .(low_noise_sigma))),
  nplots = 2
)
dev.off()

pdf("images/chapter2/one_row_singular_values.pdf", width = 10, height = 4.5)
draw_singular_values_compare()
dev.off()

png("images/chapter2/one_row_singular_values.png", width = 2000, height = 900, res = 180)
draw_singular_values_compare()
dev.off()

draw_component_examples <- function() {
  layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE),
         heights = c(1, 1, 0.18))
  par(oma = c(0, 0, 2, 0))
  par(mar = c(4, 4, 3, 1))
  plot_component_panel(
    noisy_02,
    signal_matrix,
    examples$wrong,
    sprintf("Смещение максимума, строка %d", examples$wrong)
  )
  plot_component_panel(
    noisy_low,
    signal_matrix,
    examples$wrong,
    sprintf("Та же строка при sigma = %.2f, строка %d", low_noise_sigma, examples$wrong)
  )
  plot_component_panel(
    noisy_02,
    signal_matrix,
    examples$smear,
    sprintf("Растекание около верного пикселя, строка %d", examples$smear)
  )
  plot_component_panel(
    noisy_low,
    signal_matrix,
    examples$smear,
    sprintf("Та же строка при sigma = %.2f, строка %d", low_noise_sigma, examples$smear)
  )
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend(
    "center",
    horiz = TRUE,
    bty = "n",
    legend = c("1-я компонента", "истинный пиксель"),
    col = c("black", "black"),
    lwd = c(2, 1),
    lty = c(1, 3),
    cex = 0.9
  )
  mtext("Первая компонента для отдельных строк при CSSA row.row", outer = TRUE, line = 0.5)
}

draw_bad_case_vectors <- function() {
  par(mfrow = c(2, 2), mar = c(3, 4, 3, 1))
  plot(
    bad_case_row5,
    type = "h",
    lwd = 2,
    col = "grey20",
    xlab = "Индекс",
    ylab = "Амплитуда",
    main = paste("Плохой случай:", bad_case$signal_ind, "->", which.max(bad_case_row5))
  )
  abline(v = bad_case$signal_ind, col = "firebrick", lty = 2)
  abline(v = which.max(bad_case_row5), col = "steelblue", lty = 2)

  plot(
    noise_only_row5,
    type = "h",
    lwd = 2,
    col = "grey20",
    xlab = "Индекс",
    ylab = "Амплитуда",
    main = paste("Тот же шум без сигнала:", which.max(noise_only_row5))
  )
  abline(v = which.max(noise_only_row5), col = "steelblue", lty = 2)

  plot_complex_parts(
    bad_case$fit$U[, 5],
    "bad case: U[, 5]"
  )
  plot_complex_parts(
    noise_only_fit$U[, 5],
    "noise only: U[, 5]"
  )
}

draw_frequency_profile <- function() {
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  plot(
    diag_02$row_id,
    diag_02$freq,
    type = "b",
    pch = 16,
    cex = 0.7,
    xlab = "Номер строки",
    ylab = "Нормированная частота",
    main = "Изменение частоты по строкам"
  )
  plot(
    diag_02$row_id,
    diag_02$mass1,
    type = "n",
    xlab = "Номер строки",
    ylab = "Доля массы 1-й компоненты\nв окрестности истинного пикселя",
    main = "Качество выделения 1-й компоненты"
  )
  points(diag_02$row_id[diag_02$argmax1 == diag_02$true_col],
         diag_02$mass1[diag_02$argmax1 == diag_02$true_col],
         pch = 16, cex = 0.7)
  points(diag_02$row_id[diag_02$argmax1 != diag_02$true_col],
         diag_02$mass1[diag_02$argmax1 != diag_02$true_col],
         pch = 1, cex = 0.9)
  abline(h = 0.55, lty = 3)
  legend(
    "bottomright",
    bty = "n",
    legend = c("argmax 1-й компоненты совпадает", "argmax 1-й компоненты не совпадает"),
    pch = c(16, 1),
    cex = 0.8
  )
}

pdf("images/chapter2/row_row_component_examples.pdf", width = 11, height = 8)
draw_component_examples()
dev.off()

png("images/chapter2/row_row_component_examples.png", width = 2200, height = 1600, res = 200)
draw_component_examples()
dev.off()

pdf("images/chapter2/one_row_bad_case_vectors.pdf", width = 10, height = 4.5)
draw_bad_case_vectors()
dev.off()

png("images/chapter2/one_row_bad_case_vectors.png", width = 2000, height = 900, res = 180)
draw_bad_case_vectors()
dev.off()

pdf("images/chapter2/row_row_frequency_profile.pdf", width = 10, height = 4.5)
draw_frequency_profile()
dev.off()

png("images/chapter2/row_row_frequency_profile.png", width = 2000, height = 900, res = 180)
draw_frequency_profile()
dev.off()
