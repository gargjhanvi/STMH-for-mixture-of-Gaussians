# Load libraries
library(ggplot2)
library(dplyr)
library(cubature)

set.seed(123)

# 2D mixture of Gaussians
p <- function(D, x1, x2) {
  mu <- D / (2 * sqrt(2))
  0.5 * dnorm(x1, mean = mu, sd = 1) * dnorm(x2, mean = mu, sd = 1) +
    0.5 * dnorm(x1, mean = -mu, sd = 1) * dnorm(x2, mean = -mu, sd = 1)
}

# Negative log-density
f <- function(D, x1, x2) {
  -log(p(D, x1, x2))
}

# Inverse temperature schedule
beta <- function(D) {
  beta1 <- 1 / 400
  vec <- c(beta1)
  repeat {
    a <- 1
    b <- -10 / 4 * beta1
    c <- 5 / 4 * beta1^2
    new <- (-b + sqrt(b^2 - 4 * a * c)) / (2 * a)
    if (new >= 1) break
    vec <- c(vec, new)
    beta1 <- new
  }
  c(vec, 1)
}

# Compute normalization constant Z(D, beta)
Z <- function(D, beta_val) {
  integrand <- function(x) {
    p(D, x[1], x[2])^beta_val
  }
  center <- D / (2 * sqrt(2))
  support_radius <- 4
  lower <- c(-center - support_radius, -center - support_radius)
  upper <- c(center + support_radius, center + support_radius)
  
  result <- cubintegrate(
    f = integrand,
    lower = lower,
    upper = upper,
    relTol = 1e-8,
    absTol = 1e-12
  )
  
  return(result$integral)
}

# ----- Simulation Parameters -----
lambda <- 0.5
iterations <- 1000
N <- 5000
eta <- 400
D <- 20

# Compute beta schedule and normalization constants
beta_D <- beta(D)
L <- length(beta_D)
Z_vals <- sapply(beta_D, function(b) Z(D, b))

# Initialize matrices to store x1 and x2 samples
samples_matrix_x1 <- matrix(nrow = iterations, ncol = N)
samples_matrix_x2 <- matrix(nrow = iterations, ncol = N)

# ----- Run STMH Chain -----
for (iter in 1:iterations) {
  x1 <- 10
  x2 <- 10
  i <- 1
  n <- 1
  
  while (n <= N) {
    b <- rbinom(1, 1, lambda)
    
    if (b == 1) {
      # Propose temperature swap
      inew <- i + sample(c(-1, 1), size = 1)
      if (inew >= 1 && inew <= L) {
        fx <- f(D, x1, x2)
        if (is.finite(fx)) {
          swap_ratio <- (Z_vals[i] * exp(-beta_D[inew] * fx)) /
            (Z_vals[inew] * exp(-beta_D[i] * fx))
        } else {
          swap_ratio <- 0
        }
        if (runif(1) <= min(1, swap_ratio)) {
          i <- inew
        }
      }
    } else {
      # Propose move via MH 
      y1 <- rnorm(1, mean = x1, sd = sqrt(eta))
      y2 <- rnorm(1, mean = x2, sd = sqrt(eta))
      fx <- f(D, x1, x2)
      fy <- f(D, y1, y2)
      if (runif(1) <= min(1, exp(-beta_D[i] * (fy - fx)))) {
        x1 <- y1
        x2 <- y2
      }
    }
    
    # Store x1 and x2 samples
    samples_matrix_x1[iter, n] <- x1
    samples_matrix_x2[iter, n] <- x2
    n <- n + 1
  }
}

# ----- Convergence Analysis -----
mean_x1 <- colMeans(samples_matrix_x1)
mean_x2 <- colMeans(samples_matrix_x2)

cumulative_avg_x1 <- cumsum(mean_x1) / (1:N)
cumulative_avg_x2 <- cumsum(mean_x2) / (1:N)

# Plot for x1
plot(1:N, cumulative_avg_x1, type = "l", col = "blue", lwd = 2,
     xlab = "Number of Steps (N)",
     ylab = expression(hat(mu)^N[1]))
abline(h = 0, col = "red", lty = 2)

# Plot for x2
plot(1:N, cumulative_avg_x2, type = "l", col = "blue", lwd = 2,
     xlab = "Number of Steps (N)",
     ylab = expression(hat(mu)^N[2]))
abline(h = 0, col = "red", lty = 2)
