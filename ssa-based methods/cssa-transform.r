library("Rssa")
library("lattice")

# Create a zero matrix
new.matrix <- function(N_t = 100, N_c = 100){
  matrix(0, nrow = N_c, ncol = N_t)
}

# Adding lines in matrix
add.line <- function(m, a = 1 , b = 0){
  N_c <- nrow(m)
  N_t <- ncol(m)
  k <- 1
  for (k in 1:N_t){
    if (a*k + b >= 1 && a*k + b <= N_c){
      m[a*k + b, k] <- 1
      k <- k+1
    }
  }
  m
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
dft <- function(m){
  N_c <- nrow(m)
  N_t <- ncol(m)
  for (i in (1:N_c)){
    m[i,] <- fft(m[i,])
  }
  m
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
    s <- ssa(m[i,], kind = "cssa", svd.method = "svd")
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

# Step IDFT by rows
idft.row <- function(m){
  N_c <- nrow(m)
  N_t <- ncol(m)
  for (i in (1:N_c)){
    m[i,] <- fft(m[i,], inverse=TRUE)/length(m[i,])
  }
  m
}

# Step IDFT by columns
idft.col <- function(m){
  N_c <- nrow(m)
  N_t <- ncol(m)
  for (i in (1:N_t)){
    m[,i] <- fft(m[,i], inverse=TRUE)/length(m[,i])
  }
  m
}

# Functional for Complex SSA by rows
angle.fun <- function(V){
  P <- Re(V)
  Q <- Im(Q)
  angle <- function(P1,P2,Q1,Q2){
    acos((P1*P2 + Q1*Q2)/sqrt(P1^2+Q1^2)/sqrt(P2^2+Q2^2))
  }
  var(angle(P[-length(P)],P[-1],Q[-length(Q)],Q[-1]))
}

length.fun <- function(V){
  P <- Re(V)
  Q <- Im(Q)
  sum(var(P^2 + Q^2)) 
}


# Drawing matrix
# If from.0.to.1 = TRUE, all values less than 0, replaced by 0,
# all values greater than 1, are replaced by 1

plot.matrix <- function(m, from.0.to.1 = FALSE){
  m <- t(m) # for graph-like image
  
  
  rgb.palette <- colorRampPalette(c("white", "black"), space = "rgb")
  if (!from.0.to.1){ 
    levelplot(m, xlab="", ylab="", col.regions=rgb.palette,  
              scales = list(draw=FALSE), auto.key= FALSE, 
              colorkey = FALSE)
  }
  else {
    for (i in 1:nrow(m)){
      for (j in 1:ncol(m)){
        if (m[i,j] > 1) {
          m[i,j] <- 1
        }
      }
    }
    levelplot(m, xlab="", ylab="", col.regions=rgb.palette,  
              scales = list(draw=FALSE), auto.key= FALSE, 
              at=seq(0,1,0.01), colorkey =FALSE)
  }
}

# Rotation matrix - for the correct display of matrices
rotate <- function(x) t(apply(x, 2, rev))