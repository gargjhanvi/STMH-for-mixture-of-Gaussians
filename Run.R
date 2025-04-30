# Load libraries
library(ggplot2)
library(dplyr)
library(cubature)

set.seed(123)

# 2D mixture of Gaussians density
p <- function(D, x1, x2) {
  mu <- D / (2 * sqrt(2))
  0.5 * dnorm(x1, mean = mu, sd = 1) * dnorm(x2, mean = mu, sd = 1) +
  0.5 * dnorm(x1, mean = -mu, sd = 1) * dnorm(x2, mean = -mu, sd = 1)
}

# Negative log-density
f <- function(D, x1, x2) {
  -log(p(D, x1, x2))
}

# Generate inverse temperature schedule
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

# ----- Simulation parameters -----
D <- 4
lambda <- 0.5
iterations <- 1000
N <- 100000
eta <- 400

# Compute beta schedule and normalization constants
beta_D <- beta(D)
L <- length(beta_D)
Z_vals <- sapply(beta_D, function(b) Z(D, b))

# Initialze matrix to store samples
samples_matrix <- matrix(nrow = iterations, ncol = N)

# ----- Run STMH chain -----
for (iter in 1:iterations) {
  x1 <- rnorm(1, mean = 10 / sqrt(2), sd = 1)
  x2 <- rnorm(1, mean = 10 / sqrt(2), sd = 1)
  i <- 1
  n <- 1
  
  while (n <= N) {
    b <- rbinom(1, 1, lambda)
    
    if (b == 1) {
      # Propose temperature swap
      inew <- i + sample(c(-1, 1), size = 1)
      if (inew >= 1 && inew <= L) {
        fx <- f(D, x1, x2)
        u <- runif(1)
        if (is.finite(fx)) {
          swap_ratio <- (Z_vals[i] * exp(-beta_D[inew] * fx)) /
            (Z_vals[inew] * exp(-beta_D[i] * fx))
        } else {
          swap_ratio <- 0
        }
        if (u <= min(1, swap_ratio)) {
          i <- inew
        }
      }
    } else {
      # Propose spatial move via MH
      y1 <- rnorm(1, mean = x1, sd = sqrt(eta))
      y2 <- rnorm(1, mean = x2, sd = sqrt(eta))
      fx <- f(D, x1, x2)
      fy <- f(D, y1, y2)
      u <- runif(1)
      if (u <= min(1, exp(-beta_D[i] * (fy - fx)))) {
        x1 <- y1
        x2 <- y2
      }
    }
    
    # Save sample (e.g., x1 coordinate)
    samples_matrix[iter, n] <- sqrt(x1^2 + x2^2)
    n <- n + 1
  }
}

# Analyze convergence 
mean <- colMeans(samples_matrix)
plot(mean, type = "l", col = "blue",
     xlab = "Iteration", ylab = "Mean")
