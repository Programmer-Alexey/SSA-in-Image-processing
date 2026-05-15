library("Rcpp")

.find_project_root <- function(start_dir) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in 1:20) {
    if (file.exists(file.path(cur, "SSA-in-Image-processing.Rproj"))) return(cur)
    next_dir <- dirname(cur)
    if (identical(next_dir, cur)) break
    cur <- next_dir
  }
  normalizePath(start_dir, winslash = "/", mustWork = FALSE)
}

.project_root_option <- getOption("ssa_image_project_root", NULL)
.project_root <- if (!is.null(.project_root_option) && dir.exists(.project_root_option)) {
  normalizePath(.project_root_option, winslash = "/", mustWork = TRUE)
} else {
  .this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  .start_dir <- if (!is.null(.this_file)) dirname(normalizePath(.this_file)) else getwd()
  .find_project_root(.start_dir)
}
.hough_dir <- file.path(.project_root, "hough transform")

source(file.path(.project_root, "ssa-based methods", "cssa-transform.r"))
sourceCpp(file.path(.hough_dir, "convolution.cpp"))

intensity_detector_parameter <- 0.8
gradient_detector_parameter <- 2
laplace_detector_parameter <- 2

# CSSA Detector
cssa_detector <- function(matrix, num_of_lines = 2, method = "row.row") {
  if (method == "col.col") {
    cleaned_matrix <- matrix |> dft() |> cssa.col(num.line = num_of_lines) |> idft.col() |> Re()
  } else if (method == "row.row") {
    cleaned_matrix <- matrix |> dft() |> cssa.row(num.line = num_of_lines) |> idft.row() |> Re()
  } else if (method == "col.row") {
    cleaned_matrix <- matrix |> dft() |> cssa.col(num.line = num_of_lines) |> idft.row() |> Re()
  } else if (method == "row.col") {
    cleaned_matrix <- matrix |> dft() |> cssa.row(num.line = num_of_lines) |> idft.col() |> Re()
  }
  
  
  cleaned_matrix <- ifelse(cleaned_matrix < 0.01, 0, ifelse(cleaned_matrix > 0.99, 1, cleaned_matrix))
  threshold_value <- quantile(abs(cleaned_matrix), 0.9)
  result_matrix <- ifelse(abs(cleaned_matrix) > threshold_value, cleaned_matrix, 0)
  result_matrix <- ifelse(result_matrix < 0.01, 0, ifelse(result_matrix > 0.99, 1, result_matrix))
  
  return(result_matrix)
}




