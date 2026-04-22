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
  stop("Could not find project root.")
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
slides_dir <- file.path(project_root, "VKR", "VKR_SLIDES", "SSA Slides")
images_dir <- file.path(slides_dir, "images")
dir.create(images_dir, recursive = TRUE, showWarnings = FALSE)

n_col <- 200
n_row <- 240

# Peaks are read from the accumulator slide using the grid (Delta rho, Delta theta) = (1, 0.01).
peaks <- data.frame(
  rho = c(0, 85),
  theta = c(2.03, 0.79)
)
peaks$a <- -tan(peaks$theta)
peaks$b <- peaks$rho / cos(peaks$theta)

draw_segment <- function(a, b, col) {
  x <- seq(1, n_col, length.out = 1000)
  y <- a * x + b
  keep <- y >= 1 & y <= n_row
  if (any(keep)) {
    lines(x[keep], y[keep], col = col, lwd = 4)
  }
}

draw_recovered_lines <- function() {
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot(
    NA,
    xlim = c(1, n_col),
    ylim = c(n_row, 1),
    asp = 1,
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  rect(1, 1, n_col, n_row, border = "gray40", lwd = 1.2)
  draw_segment(peaks$a[1], peaks$b[1], "#D43F3A")
  draw_segment(peaks$a[2], peaks$b[2], "#1F77B4")
}

pdf(file.path(images_dir, "accumulator_recovered_lines.pdf"), width = 3.5, height = 4.0)
draw_recovered_lines()
dev.off()

png(file.path(images_dir, "accumulator_recovered_lines.png"), width = 700, height = 800, res = 180)
draw_recovered_lines()
dev.off()

cat("Generated accumulator_recovered_lines.pdf and accumulator_recovered_lines.png\n")

