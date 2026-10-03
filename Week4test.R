# ==============================================================================
# Tutorial 4 - Exercise 4: Monte Carlo Simulation (No Custom Functions)
# ==============================================================================

set.seed(4008)

# Parameters specified in the exercise
N <- 1000       # Number of replications
c_true <- 1.0   # True parameter
c0 <- 1.0       # Null hypothesis value: H0: c = 1
phi <- 0.5      # MA(1) coefficient
crit_val <- 1.959964  # Standard normal 5% critical value (two-sided)

# ==============================================================================
# Simulation 1: T = 100
# ==============================================================================
T_val <- 100
m <- floor(4 * (T_val / 100)^(2/9)) # Newey-West bandwidth parameter

# Counters for false rejections
reject_short_100 <- 0
reject_nw_100    <- 0
reject_exact_100 <- 0

for (i in 1:N) {
  # 1. Generate innovations and error process
  eps <- rnorm(T_val + 1, mean = 0, sd = 1)
  u <- eps[2:(T_val + 1)] + phi * eps[1:T_val]
  y <- c_true + u
  
  # 2. Estimate parameter and residuals
  c_hat <- mean(y)
  u_hat <- y - c_hat
  
  # 3. Variance Estimators
  # (i) Short-run variance
  s2_short <- mean(u_hat^2)
  
  # (ii) Newey-West estimator with Bartlett weights
  gamma_0 <- mean(u_hat^2)
  nw_term <- 0
  for (j in 1:m) {
    gamma_j <- (1 / T_val) * sum(u_hat[(j + 1):T_val] * u_hat[1:(T_val - j)])
    weight  <- 1 - (j / (m + 1))
    nw_term <- nw_term + 2 * weight * gamma_j
  }
  s2_nw <- max(gamma_0 + nw_term, 1e-8)
  
  # (iii) Exact closed-form MA(1) estimator
  gamma_1 <- (1 / (T_val - 1)) * sum(u_hat[2:T_val] * u_hat[1:(T_val - 1)])
  s2_exact <- max(gamma_0 + 2 * (1 - 1 / T_val) * gamma_1, 1e-8)
  
  # 4. Compute t-statistics and check rejection
  t_short <- (c_hat - c0) / sqrt(s2_short / T_val)
  t_nw    <- (c_hat - c0) / sqrt(s2_nw / T_val)
  t_exact <- (c_hat - c0) / sqrt(s2_exact / T_val)
  
  if (abs(t_short) > crit_val) reject_short_100 <- reject_short_100 + 1
  if (abs(t_nw) > crit_val)    reject_nw_100    <- reject_nw_100 + 1
  if (abs(t_exact) > crit_val) reject_exact_100 <- reject_exact_100 + 1
}


# ==============================================================================
# Simulation 2: T = 1000
# ==============================================================================
T_val <- 1000
m <- floor(4 * (T_val / 100)^(2/9)) # Updated bandwidth for T = 1000

# Counters for false rejections
reject_short_1000 <- 0
reject_nw_1000    <- 0
reject_exact_1000 <- 0

for (i in 1:N) {
  # 1. Generate innovations and error process
  eps <- rnorm(T_val + 1, mean = 0, sd = 1)
  u <- eps[2:(T_val + 1)] + phi * eps[1:T_val]
  y <- c_true + u
  
  # 2. Estimate parameter and residuals
  c_hat <- mean(y)
  u_hat <- y - c_hat
  
  # 3. Variance Estimators
  # (i) Short-run variance
  s2_short <- mean(u_hat^2)
  
  # (ii) Newey-West estimator with Bartlett weights
  gamma_0 <- mean(u_hat^2)
  nw_term <- 0
  for (j in 1:m) {
    gamma_j <- (1 / T_val) * sum(u_hat[(j + 1):T_val] * u_hat[1:(T_val - j)])
    weight  <- 1 - (j / (m + 1))
    nw_term <- nw_term + 2 * weight * gamma_j
  }
  s2_nw <- max(gamma_0 + nw_term, 1e-8)
  
  # (iii) Exact closed-form MA(1) estimator
  gamma_1 <- (1 / (T_val - 1)) * sum(u_hat[2:T_val] * u_hat[1:(T_val - 1)])
  s2_exact <- max(gamma_0 + 2 * (1 - 1 / T_val) * gamma_1, 1e-8)
  
  # 4. Compute t-statistics and check rejection
  t_short <- (c_hat - c0) / sqrt(s2_short / T_val)
  t_nw    <- (c_hat - c0) / sqrt(s2_nw / T_val)
  t_exact <- (c_hat - c0) / sqrt(s2_exact / T_val)
  
  if (abs(t_short) > crit_val) reject_short_1000 <- reject_short_1000 + 1
  if (abs(t_nw) > crit_val)    reject_nw_1000    <- reject_nw_1000 + 1
  if (abs(t_exact) > crit_val) reject_exact_1000 <- reject_exact_1000 + 1
}


# ==============================================================================
# Print Results Table
# ==============================================================================
results <- data.frame(
  Estimator = c("Short-run (s^2 = hat{sigma}_u^2)", 
                "Newey-West (s^2 = hat{lambda}_1^2)", 
                "MA(1) closed-form (s^2 = hat{lambda}_2^2)"),
  Size_T100 = c(reject_short_100 / N, reject_nw_100 / N, reject_exact_100 / N),
  Size_T1000 = c(reject_short_1000 / N, reject_nw_1000 / N, reject_exact_1000 / N)
)

print(results)