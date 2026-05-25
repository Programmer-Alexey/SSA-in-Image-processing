library("Rssa")
library("lattice")
library(gridExtra)

# Create a zero matrix
new.matrix <- function(N_t = 100, N_c = 100){
  matrix(0, nrow = N_c, ncol = N_t)
}

# Adding lines in matrix
add.line <- function(m, a = 1, b = 0, method = "bresenham", intensity = 1) {
  # method <- match.arg(method)
  N_c <- nrow(m)
  N_t <- ncol(m)
  
  if (method == "default") {
    for (k in 1:N_t) {
      y <- a * k + b
      if (y >= 1 && y <= N_c) {
        m[round(y), k] <- m[round(y), k] + intensity
      }
    }
    
  } else if (method == "bresenham") {
    # РђР»РіРѕСЂРёС‚Рј Р‘СЂРµР·РµРЅС…РµРјР°
    # РџСЂРµРѕР±СЂР°Р·СѓРµРј СѓСЂР°РІРЅРµРЅРёРµ y = a*x + b РІ РєРѕРѕСЂРґРёРЅР°С‚С‹ РЅР°С‡Р°Р»Р° Рё РєРѕРЅС†Р°
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
        m[y0, x0] <- m[y0, x0] + intensity
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

plot.matrix <- function(m, from.0.to.1 = FALSE, heatmap = "none", labels = NULL, nplots = NULL){
  rgb.palette <- colorRampPalette(c("white", "black"), space = "rgb")
  sign.palette <- colorRampPalette(c("blue", "white", "black"), space = "rgb")
  diverging.palette <- colorRampPalette(c("navy", "white", "firebrick"), space = "rgb")
  accumulator.palette <- colorRampPalette(c("white", "yellow", "red"), space = "rgb")

  as_plot_matrix <- function(x) {
    if (is.null(x)) {
      stop("Cannot plot NULL matrix")
    }
    x <- as.matrix(x)
    if (is.complex(x)) x <- Re(x)
    x
  }

  normalize_heatmap <- function(value) {
    if (is.logical(value)) {
      return(if (isTRUE(value)) "yellow_red" else "none")
    }
    if (is.null(value)) {
      return("none")
    }
    value <- as.character(value)[1L]
    aliases <- c(
      "false" = "none",
      "no" = "none",
      "none" = "none",
      "blue-black" = "blue_black",
      "blue_black" = "blue_black",
      "blue.red" = "blue_red",
      "blue-red" = "blue_red",
      "blue_red" = "blue_red",
      "yellow.red" = "yellow_red",
      "yellow-red" = "yellow_red",
      "yellow_red" = "yellow_red",
      "accumulator" = "yellow_red",
      "true" = "yellow_red"
    )
    key <- tolower(value)
    if (!key %in% names(aliases)) {
      stop("Unknown heatmap mode: ", value)
    }
    aliases[[key]]
  }

  heatmap_mode <- normalize_heatmap(heatmap)

  make_levelplot <- function(mat, main = NULL) {
    if (from.0.to.1) {
      mat[mat < 0] <- 0
      mat[mat > 1] <- 1
    }

    if (heatmap_mode == "none") {
      return(levelplot(mat, xlab="", ylab="",
                       col.regions=rgb.palette,
                       scales=list(draw=FALSE), auto.key=FALSE,
                       at = if(from.0.to.1) seq(0,1,0.01) else NULL,
                       main = main,
                       colorkey=FALSE,
                       col="transparent", border=NA, cuts=255))
    }

    finite_values <- mat[is.finite(mat)]
    if (length(finite_values) == 0L) {
      finite_values <- 0
      mat[] <- 0
    }
    mat[!is.finite(mat)] <- 0

    if (heatmap_mode == "blue_red") {
      value_range <- range(finite_values)
      if (value_range[1L] == value_range[2L]) {
        value_range <- value_range + c(-0.5, 0.5)
      }
      return(levelplot(mat, xlab="", ylab="",
                       col.regions=diverging.palette,
                       scales=list(draw=FALSE), auto.key=FALSE,
                       at=seq(value_range[1L], value_range[2L], length.out=257),
                       main = main,
                       colorkey=FALSE,
                       col="transparent", border=NA, cuts=255))
    }

    if (heatmap_mode == "yellow_red") {
      mat <- pmax(mat, 0)
      value_max <- max(mat[is.finite(mat)])
      if (!is.finite(value_max) || value_max == 0) {
        value_max <- 1
      }
      return(levelplot(mat, xlab="", ylab="",
                       col.regions=accumulator.palette,
                       scales=list(draw=FALSE), auto.key=FALSE,
                       at=seq(0, value_max, length.out=257),
                       main = main,
                       colorkey=FALSE,
                       col="transparent", border=NA, cuts=255))
    }

    max_abs <- if (from.0.to.1) 1 else max(abs(finite_values))
    if (!is.finite(max_abs) || max_abs == 0) {
      max_abs <- 1
    }

    levelplot(mat, xlab="", ylab="",
              col.regions=sign.palette,
              scales=list(draw=FALSE), auto.key=FALSE,
              at=seq(-max_abs, max_abs, length.out=257),
              main = main,
              colorkey=FALSE,
              col="transparent", border=NA, cuts=255)
  }
  
  # РѕРґРёРЅРѕС‡РЅР°СЏ РјР°С‚СЂРёС†Р°
  if (!is.list(m)) {
    m <- t(as_plot_matrix(m))
    return(make_levelplot(m))
  }
  
  # СЃРїРёСЃРѕРє РјР°С‚СЂРёС†
  plots <- list()
  for (i in seq_along(m)) {
    if (is.null(m[[i]])) {
      item_name <- names(m)[i]
      if (is.null(item_name) || item_name == "") {
        item_name <- paste0("#", i)
      }
      stop("Cannot plot NULL matrix in list item: ", item_name)
    }
    mat <- t(as_plot_matrix(m[[i]]))
    p <- make_levelplot(mat, main = if(!is.null(labels)) labels[i] else " ")
    plots[[i]] <- p
  }
  
  # СЃРµС‚РєР°
  do.call(grid.arrange, c(plots, ncol = nplots %||% length(plots)))
}

# Rotation matrix - for the correct display of matrices
rotate <- function(x) t(apply(x, 2, rev))
