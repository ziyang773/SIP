# moving block bootstrap
library(boot)

# Romano & Thombs (1996) Moving Block Bootstrap for autocorrelations
mbb_rt96 <- function(x, lag = 1, b = 10, R = 1000, conf_level = 0.95) {
  n <- length(x)
  k <- ceiling(n / b)               # Number of blocks needed
  q <- n - b + 1                     # Total overlapping blocks
  blocks <- lapply(1:q, function(i) x[i:(i + b - 1)])

  # 1. Compute observed autocorrelation
  obs_rho <- cor(x[1:(n-lag)], x[(lag+1):n])

  # 2. Generate R bootstrap replicates
  boot_rhos <- replicate(R, {
    sampled_indices <- sample(1:q, size = k, replace = TRUE)
    x_star <- unlist(blocks[sampled_indices])[1:n]
    cor(x_star[1:(n-lag)], x_star[(lag+1):n])
  })

  # 3. Center around observed rho and scale by sqrt(n)
  scaled_boot <- sqrt(n) * abs(boot_rhos - obs_rho)

  # 4. Construct symmetric confidence interval
  alpha <- 1 - conf_level
  q_upper <- quantile(scaled_boot, 1 - alpha/2)
  ci_lower <- obs_rho - q_upper / sqrt(n)
  ci_upper <- obs_rho + q_upper / sqrt(n)

  # 5. p-value for H0: rho = 0
  obs_scaled <- sqrt(n) * abs(obs_rho - 0)       # distance from null
  p_val <- mean(scaled_boot >= obs_scaled)

  list(
    obs_rho = obs_rho,
    scaled_boot = scaled_boot,
    ci = c(ci_lower, ci_upper),
    p_value = p_val
  )
}

# -----------------------
# RT96-style MBB for multiple lags
mbb_multi_lag <- function(x, lags = 1:5, b = 10, R = 1000) {
  n <- length(x)
  k <- ceiling(n / b)
  q <- n - b + 1
  blocks <- lapply(1:q, function(i) x[i:(i + b - 1)])

  # 1. Observed autocorrelations
  obs_rhos <- sapply(lags, function(lag) cor(x[1:(n-lag)], x[(lag+1):n]))

  # 2. Bootstrap replicates
  boot_rhos <- replicate(R, {
    sampled_indices <- sample(1:q, size = k, replace = TRUE)
    x_star <- unlist(blocks[sampled_indices])[1:n]
    sapply(lags, function(lag) cor(x_star[1:(n-lag)], x_star[(lag+1):n]))
  })

  # 3. Compute max scaled deviation across lags (joint test)
  scaled_boot <- apply(boot_rhos, 2, function(rho_star) sqrt(n) * abs(rho_star - obs_rhos))
  max_scaled_boot <- apply(scaled_boot, 2, max)  # max across lags for each bootstrap

  obs_scaled <- sqrt(n) * abs(obs_rhos)           # distance from null for each lag
  obs_stat <- max(obs_scaled)                     # max across lags

  # 4. p-value for joint test
  p_val <- mean(max_scaled_boot >= obs_stat)

  list(obs_rhos = obs_rhos, p_value = p_val)
}

dwb_dependence_test <- function(x, L = 5, B = 1000, bandwidth = 20){
  n <- length(x)
  mu <- mean(x)
  x_centered <- x - mu
  var_x <- sum(x_centered^2) / n

  # 1. Compute observed max-scaled autocorrelation
  # use the formula: sum((x_t-mu)*(x_{t+k}-mu)) / (n * var_x)
  obs_rho <- sapply(1:L, function(k) {
    sum(x_centered[1:(n-k)] * x_centered[(k+1):n]) / (n * var_x)
  })
  T_obs <- max(sqrt(n) * abs(obs_rho))

  # 2. Optimized DWB Multiplier Generator
  # Using a moving average of normals to create the Bartlett kernel dependence
  generate_dwb <- function(n, l){
    # Using a window of size l
    w <- rnorm(n + l)
    xi <- numeric(n)
    # The kernel weights for a Bartlett-like dependence
    weights <- (1:l) / sqrt(sum((1:l)^2))
    for(t in 1:n){
      xi[t] <- sum(w[t:(t+l-1)] * weights)
    }
    return(xi)
  }

  # 3. Bootstrap replicates of the statistic
  T_star <- numeric(B)
  for(b in 1:B){
    xi <- generate_dwb(n, bandwidth)

    # Core DWB step: multiply centered series by dependent weights
    # Then calculate the "bootstrap" autocorrelation
    rho_star <- sapply(1:L, function(k) {
      z_k <- x_centered[1:(n-k)] * x_centered[(k+1):n]
      # Weighted sum divided by original variance
      sum(z_k * xi[1:(n-k)]) / (n * var_x)
    })

    T_star[b] <- max(sqrt(n) * abs(rho_star))
  }

  p_val <- mean(T_star >= T_obs)

  return(list(p_value = p_val, T_obs = T_obs, obs_rho = obs_rho))
}
