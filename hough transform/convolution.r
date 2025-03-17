convolution <- function(matr, kernel) { # this is a version with padding to avoid incorrect size of result matrix
  matr.ncol <- ncol(matr)
  matr.nrow <- nrow(matr)
  
  kernel.ncol <- ncol(kernel)
  kernel.nrow <- nrow(kernel)
  
  pad_x <- floor(kernel.ncol / 2)
  pad_y <- floor(kernel.nrow / 2)
  
  # Make a large matrix
  padded_matr <- matrix(0, nrow = matr.nrow + 2 * pad_y, ncol = matr.ncol + 2 * pad_x)
  padded_matr[(pad_y + 1):(pad_y + matr.nrow), (pad_x + 1):(pad_x + matr.ncol)] <- matr
  
  # Result matrix
  result <- matrix(0, nrow = matr.nrow, ncol = matr.ncol)
  
  # standart convolution
  for (row in 1:matr.nrow) {
    for (col in 1:matr.ncol) {
      submatrix <- padded_matr[row:(row + kernel.nrow - 1), col:(col + kernel.ncol - 1)]
      result[row, col] <- sum(submatrix * kernel)
    }
  }
  
  return(result)
}


# Intensity
intensity_detector_parameter <- 0.8
intensity_detector <- function(matrix, threshold=intensity_detector_parameter){
  matrix <- ifelse(matrix > quantile(matrix, threshold), matrix, 0)
  matrix
}

# Gradient
gradient_detector_parameter <- 2
intensity_gradient_detector <- function(matrix, threshold=gradient_detector_parameter){
  sobel_x <- matrix(c(-1, 0, 1, 
                      -2, 0, 2, 
                      -1, 0, 1), 
                    nrow = 3, byrow = TRUE)
  
  G_x <- convolution(matrix, sobel_x)
  G_y <- convolution(matrix, t(sobel_x))
  matrix <- sqrt(G_x*G_x + G_y*G_y)
  matrix <- ifelse(matrix > threshold, matrix, 0)
}

# Laplacian
zero_crossings <- function(matrix) {
  rows <- nrow(matrix)
  cols <- ncol(matrix)
  edges <- matrix(0, nrow = rows, ncol = cols)
  
  for (i in 2:(rows-1)) {
    for (j in 2:(cols-1)) {
      neg <- sum(matrix[(i-1):(i+1), (j-1):(j+1)] < 0)
      pos <- sum(matrix[(i-1):(i+1), (j-1):(j+1)] > 0)
      if (neg > 0 && pos > 0) {
        edges[i, j] <- 1
      }
    }
  }
  return(edges)
}

gaussian_kernel <- function(n, sigma){
  # size of kernel is odd
  
  kernel <- matrix(0, ncol=n, nrow=n)
  center_x <- n %/% 2 + 1
  center_y <- center_x
  
  for(col in 1:n){
    for(row in 1:n){
      kernel[row, col] = exp(  -( (col-center_y)*(col-center_y) + (row-center_x)*(row-center_x)) / (2*sigma*sigma)  ) / (2*pi*sigma*sigma)
    }
  }
  
  return(kernel)
}

laplace_detector_parameter <- 2
laplace_detector<- function(matrix, sigma=laplace_detector_parameter){
  n <- as.integer(6*sigma) - 1
  gauss.kernel <- gaussian_kernel(n, sigma) 
  laplace_operator <- matrix(c(1, 4, 1,
                               4, -20, 4,
                               1, 4, 1), ncol=3, byrow=T)
  
  return(matrix |> convolution(gauss.kernel) |> convolution(laplace_operator) |> zero_crossings())
}