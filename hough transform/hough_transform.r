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

source(file.path(.hough_dir, "convolution.r"))
library("Rcpp")
sourceCpp(file.path(.hough_dir, "make_accumulator.cpp"))


make_accumulator0 <- function(matrix, detector, theta_step = 1, rho_step = 1) {
  m <- detector(matrix)
  edge_points <- which(m > 0, arr.ind=T)
  
  x <- edge_points[, 1]
  y <- edge_points[, 2]
  
  theta <- seq(0, pi, length.out = 180 / theta_step)
  rho_max <- ceiling(sqrt(nrow(matrix)^2 + ncol(matrix)^2))
  rho <- seq(-rho_max, rho_max, by = rho_step)
  
  accumulator <- matrix(0, nrow = length(rho), ncol = length(theta))
  
  for (i in seq_along(x)) {
    for (j in seq_along(theta)) {
      rho_value <- x[i] * cos(theta[j]) + y[i] * sin(theta[j])
      rho_idx <- round((rho_value - min(rho)) / rho_step) + 1
      if (rho_idx >= 1 && rho_idx <= length(rho)) {
        accumulator[rho_idx, j] <- accumulator[rho_idx, j] + 1
      }
    }
  }
  
  list(accumulator = accumulator, rho = rho, theta = theta)
}


sourceCpp(file.path(.hough_dir, "detect_lines.cpp"))
detect_lines0 <- function(accumulator, quant_rho, quant_theta, N, ncol=100, nrow=100){
  
  eps <- 1/max(ncol, nrow)
  output <- matrix(0, ncol=ncol, nrow=nrow)
  
  candidates_ind <- sort(accumulator, decreasing = TRUE, index.return = TRUE)$ix
  sliced_ind <- candidates_ind[1:N]
  
  # indicies of candidates
  indices <- arrayInd(sliced_ind, dim(accumulator))
  rows <- indices[, 1]
  cols <- indices[, 2]
  
  for( i in 1:(length(rows))){
    rho <- quant_rho[rows[i]]
    theta <- quant_theta[cols[i]]
    
    #print(-tan(theta))
    #print(rho/cos(theta))
    
    if (abs(sin(theta)) > eps) {
      a_value <-1/tan(theta)
      
      output <- output |> add.line(a=-tan(theta), b=rho/cos(theta))
      
    }
  }
  output
  
}

detect_lines_nms <- function(accumulator, quant_rho, quant_theta, N, ncol=100, nrow=100, suppression_window=6){
  eps <- 1/max(ncol, nrow)
  output <- matrix(0, ncol=ncol, nrow=nrow)
  
  accumulator_suppressed <- accumulator
  detected_lines <- list()
  
  for (k in 1:N) {
    # Локальные максимумы
    max_idx <- which(accumulator_suppressed == max(accumulator_suppressed), arr.ind=TRUE)
    if (nrow(max_idx) == 0) break
    
    row_idx <- max_idx[1, 1]
    col_idx <- max_idx[1, 2]
    
    rho <- quant_rho[row_idx]
    theta <- quant_theta[col_idx]
    
    #print(-tan(theta))
    #print(rho/cos(theta))
    
    # Добавляем линию, если она не слишком близка к уже найденным
    output <- output |> add.line(a=-tan(theta), b=rho/cos(theta))
    detected_lines <- c(detected_lines, list(rho=rho, theta=theta))
    
    # Подавляем окрестности найденной точки
    row_min <- max(1, row_idx - suppression_window)
    row_max <- min(nrow(accumulator), row_idx + suppression_window)
    col_min <- max(1, col_idx - suppression_window)
    col_max <- min(ncol(accumulator), col_idx + suppression_window)
    
    accumulator_suppressed[row_min:row_max, col_min:col_max] <- 0
  }
  output
}
