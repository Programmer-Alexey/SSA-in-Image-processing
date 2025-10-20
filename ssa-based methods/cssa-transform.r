library("Rssa")
library("lattice")
library(gridExtra)

# Create a zero matrix
new.matrix <- function(N_t = 100, N_c = 100){
  matrix(0, nrow = N_c, ncol = N_t)
}

# Adding lines in matrix
add.line <- function(m, a = 1, b = 0, method = "bresenham") {
  # method <- match.arg(method)
  N_c <- nrow(m)
  N_t <- ncol(m)
  
  if (method == "default") {
    for (k in 1:N_t) {
      y <- a * k + b
      if (y >= 1 && y <= N_c) {
        m[round(y), k] <- 1
      }
    }
    
  } else if (method == "bresenham") {
    # Алгоритм Брезенхема
    # Преобразуем уравнение y = a*x + b в координаты начала и конца
    x0 <- 1
    y0 <- round(a * x0 + b)
    x1 <- N_t
    y1 <- round(a * x1 + b)
    
    dx <- abs(x1 - x0)
    dy <- abs(y1 - y0)
    sx <- if (x0 < x1) 1 else -1
    sy <- if (y0 < y1) 1 else -1
    err <- dx - dy
    
    repeat {
      if (x0 >= 1 && x0 <= N_t && y0 >= 1 && y0 <= N_c) {
        m[y0, x0] <- 1
      }
      if (x0 == x1 && y0 == y1) break
      e2 <- 2 * err
      if (e2 > -dy) {
        err <- err - dy
        x0 <- x0 + sx
      }
      if (e2 < dx) {
        err <- err + dx
        y0 <- y0 + sy
      }
    }
  }
  
  return(m)
}


add.line.rand <- function(m, a = 1 , b = 0, sd=0.2){
  N_c <- nrow(m)
  N_t <- ncol(m)
  k <- 1
  for (k in 1:N_t){
    if (a*k + b >= 1 && a*k + b <= N_c){
      m[a*k + b, k] <- max(min(rnorm(1, sd=sd), 1), 0)
      k <- k+1
    }
  }
  m
}


# Adding noise in matrix
add.noise <- function(m, sigma = 0.2){
  N_c <- nrow(m)
  N_t <- ncol(m)
  m <- m + rnorm(N_c*N_t, sd=sigma)
  m
}

# Step DFT
dft <- function(m) {
  t(mvfft(t(m)))
}

# Step Complex SSA by rows
# You can use the Threshold: step = TRUE
# eps -- value threshold
# functional -- used functional
cssa.row <- function(m, num.line = 1, step = FALSE, eps = 0.001, 
                     functional = angle.fun){
  N_c <- nrow(m)
  N_t <- ncol(m)
  for (i in (1:N_c)){
    s <- ssa(m[i,], kind = "cssa", , svd.method = "svd")
    if (!step){
      r <- reconstruct(s, groups = list(Seasonality = 1:num.line))
      m[i,] <- r$Seasonality
    } 
    else{
      v <- numeric(N_t)
      for (j in 1:num.line){
        if (functional(s$U[,j]) < eps){  
          r <- reconstruct(s, groups = list(Seasonality = j))
          v <- v + r$Seasonality
        }
        m[i,] <- v
      }
    }
  }
  m
}

# Step Complex SSA by columns
cssa.col <- function(m, num.line = 1){
  N_c <- nrow(m)
  N_t <- ncol(m)
  for (i in (1:N_t)){
    s <- ssa(m[,i], kind = "cssa", svd.method = "svd")
    r <- reconstruct(s, groups = list(Seasonality = 1:num.line))
    m[,i] <- r$Seasonality
  }
  m
}



idft.row <- function(m) {
  t(mvfft(t(m), inverse = TRUE)) / ncol(m)
}

idft.col <- function(m) {
  mvfft(m, inverse = TRUE) / nrow(m)
}

# Functional for Complex SSA by rows
angle.fun <- function(V) {
  P <- Re(V); Q <- Im(V)
  P1 <- P[-length(P)]; P2 <- P[-1L]
  Q1 <- Q[-length(Q)]; Q2 <- Q[-1L]
  num <- P1*P2 + Q1*Q2
  den <- sqrt(P1*P1 + Q1*Q1) * sqrt(P2*P2 + Q2*Q2)
  cang <- pmax(-1, pmin(1, num / den))
  stats::var(acos(cang))
}

length.fun <- function(V) {
  P <- Re(V); Q <- Im(V)
  stats::var(P*P + Q*Q)
}


# Drawing matrix
# If from.0.to.1 = TRUE, all values less than 0, replaced by 0,
# all values greater than 1, are replaced by 1

plot.matrix <- function(m, from.0.to.1 = FALSE, labels = NULL, nplots = NULL){   
  rgb.palette <- colorRampPalette(c("white", "black"), space = "rgb")
  
  # одиночная матрица
  if (!is.list(m)) {
    m <- t(m)
    if (!from.0.to.1){
      return(levelplot(m, xlab="", ylab="", 
                       col.regions=rgb.palette,
                       scales=list(draw=FALSE), auto.key=FALSE, 
                       colorkey=FALSE,
                       col="transparent", border=NA, cuts=255))
    } else {
      m[m > 1] <- 1
      return(levelplot(m, xlab="", ylab="", 
                       col.regions=rgb.palette,
                       scales=list(draw=FALSE), auto.key=FALSE, 
                       at=seq(0,1,0.01), colorkey=FALSE,
                       col="transparent", border=NA, cuts=255))
    }
  }
  
  # список матриц
  plots <- list()
  for (i in seq_along(m)) {
    mat <- t(m[[i]])
    if (from.0.to.1) mat[mat > 1] <- 1
    
    p <- levelplot(mat, xlab="", ylab="", 
                   col.regions=rgb.palette,
                   scales=list(draw=FALSE), auto.key=FALSE,
                   at = if(from.0.to.1) seq(0,1,0.01) else NULL,
                   main = if(!is.null(labels)) labels[i] else " ",
                   colorkey=FALSE,
                   col="transparent", border=NA, cuts=255)
    plots[[i]] <- p
  }
  
  # сетка
  do.call(grid.arrange, c(plots, ncol = nplots %||% length(plots)))
}

# Rotation matrix - for the correct display of matrices
rotate <- function(x) t(apply(x, 2, rev))
