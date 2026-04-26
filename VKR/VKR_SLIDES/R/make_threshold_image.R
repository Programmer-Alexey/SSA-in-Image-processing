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

threshold_positive <- function(m, threshold) {
  x <- as.matrix(m)
  x[!is.finite(x)] <- 0
  ifelse(x >= threshold, 1, 0)
}

slides_dir <- file.path(project_root, "VKR", "VKR_SLIDES")
dir.create(file.path(slides_dir, "images"), recursive = TRUE, showWarnings = FALSE)

set.seed(1111)
signal_matrix <- make_line_image(100, 100, a = 2, b = -1)
noisy_matrix <- add.noise(signal_matrix, sigma = 0.2)
processed_matrix <- cssa_row_row(noisy_matrix, num_of_lines = 1L)
thresholded_matrix <- threshold_positive(processed_matrix, threshold = 0.1)

pdf(file.path(slides_dir, "images", "row_row_threshold_01_hough_input.pdf"), width = 5.5, height = 5.2)
plot.matrix(
  thresholded_matrix,
  from.0.to.1 = TRUE,
  labels = expression(hat(M)[ij] >= 0.1),
  nplots = 1
)
dev.off()

png(file.path(slides_dir, "images", "row_row_threshold_01_hough_input.png"), width = 1100, height = 1000, res = 180)
plot.matrix(
  thresholded_matrix,
  from.0.to.1 = TRUE,
  labels = expression(hat(M)[ij] >= 0.1),
  nplots = 1
)
dev.off()

write.csv(
  data.frame(
    threshold = 0.1,
    active_pixels_before = sum(processed_matrix > 0),
    active_pixels_after = sum(thresholded_matrix > 0)
  ),
  file.path(slides_dir, "tables", "row_row_threshold_01_summary.csv"),
  row.names = FALSE
)

cat("Р“РѕС‚РѕРІРѕ: images/row_row_threshold_01_hough_input.pdf\n")

