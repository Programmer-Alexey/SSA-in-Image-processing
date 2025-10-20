library("Rcpp")
source("ssa-based methods/cssa-transform.r")
sourceCpp("hough transform/convolution.cpp")

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




