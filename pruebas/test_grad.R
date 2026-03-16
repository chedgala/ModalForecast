library(numDeriv)

# simulate data
set.seed(42)
y <- arima.sim(n = 200, list(ar = 0.5, ma = 0.5))
p <- 1; q <- 1; d <- 0; n <- length(y)
mean_y <- mean(y)

neg_log_lik <- function(params) {
  c_mu <- params[1]
  phi <- if (p > 0) params[2:(p+1)] else numeric(0)
  theta <- if (q > 0) params[(p+2):(1+p+q)] else numeric(0)
  log_sigma <- params[length(params) - 1]
  log_gamma <- params[length(params)]
  sigma <- exp(log_sigma)
  gamma <- exp(log_gamma)

  mu_t <- numeric(n)
  eps <- numeric(n)
  
  for (t in 1:n) {
    ar_term <- 0
    if (p > 0) {
      for (i in seq_len(p)) {
        if (t - i > 0) ar_term <- ar_term + phi[i] * y[t - i]
        else ar_term <- ar_term + phi[i] * mean_y
      }
    }
    ma_term <- 0
    if (q > 0) {
      for (j in seq_len(q)) if (t - j > 0) ma_term <- ma_term + theta[j] * eps[t - j]
    }
    mu_t[t] <- c_mu + ar_term + ma_term
    eps[t] <- y[t] - mu_t[t]
  }

  z <- (y - mu_t) / sigma
  log_f0 <- function(val) -0.5 * log(2 * pi) - 0.5 * val^2

  idx_ge <- which(y >= mu_t)
  idx_lt <- which(y < mu_t)

  log_lik <- numeric(n)
  if (length(idx_ge) > 0) log_lik[idx_ge] <- log(2) - log_sigma - log(gamma + 1/gamma) + log_f0(z[idx_ge] / gamma)
  if (length(idx_lt) > 0) log_lik[idx_lt] <- log(2) - log_sigma - log(gamma + 1/gamma) + log_f0(z[idx_lt] * gamma)

  sum(-log_lik)
}

grad_log_lik <- function(params) {
  c_mu <- params[1]
  phi <- if (p > 0) params[2:(p+1)] else numeric(0)
  theta <- if (q > 0) params[(p+2):(1+p+q)] else numeric(0)
  log_sigma <- params[length(params) - 1]
  log_gamma <- params[length(params)]
  sigma <- exp(log_sigma)
  gamma <- exp(log_gamma)

  mu_t <- numeric(n)
  eps <- numeric(n)

  # For gradient
  num_beta <- 1 + p + q
  V <- matrix(0, nrow = n, ncol = num_beta) # V_t
  
  for (t in 1:n) {
    ar_term <- 0
    if (p > 0) {
      for (i in seq_len(p)) {
        if (t - i > 0) ar_term <- ar_term + phi[i] * y[t - i]
        else ar_term <- ar_term + phi[i] * mean_y
      }
    }
    ma_term <- 0
    if (q > 0) {
      for (j in seq_len(q)) if (t - j > 0) ma_term <- ma_term + theta[j] * eps[t - j]
    }
    mu_t[t] <- c_mu + ar_term + ma_term
    eps[t] <- y[t] - mu_t[t]
    
    # Compute Inputs for V_t
    Z_t <- numeric(num_beta)
    Z_t[1] <- 1
    if (p > 0) {
      for (i in 1:p) Z_t[1 + i] <- if (t - i > 0) y[t - i] else mean_y
    }
    if (q > 0) {
      for (j in 1:q) Z_t[1 + p + j] <- if (t - j > 0) eps[t - j] else 0
    }
    
    # Add MA terms to V_t
    V_ma <- rep(0, num_beta)
    if (q > 0) {
        for (j in 1:q) {
            if (t - j > 0) V_ma <- V_ma + theta[j] * V[t - j, ]
        }
    }
    V[t, ] <- Z_t - V_ma
  }

  G <- numeric(length(params))
  
  for (t in 1:n) {
    e_t <- eps[t]
    if (e_t >= 0) {
        W_t <- 1 / (sigma^2 * gamma^2)
        d_log_gamma <- - (gamma^2 - 1)/(gamma^2 + 1) + W_t * e_t^2
    } else {
        W_t <- gamma^2 / sigma^2
        d_log_gamma <- - (gamma^2 - 1)/(gamma^2 + 1) - W_t * e_t^2
    }
    
    # Deriv of neg log lik w.r.t beta
    G[1:num_beta] <- G[1:num_beta] - W_t * e_t * V[t, ]
    # Deriv w.r.t u_sigma
    G[length(params) - 1] <- G[length(params) - 1] + (1 - W_t * e_t^2)
    # Deriv w.r.t u_gamma
    G[length(params)] <- G[length(params)] - d_log_gamma
  }
  
  return(G)
}

init_params <- c(0.1, 0.4, 0.4, log(1), log(1))
num_g <- grad(neg_log_lik, init_params)
ana_g <- grad_log_lik(init_params)

print(num_g)
print(ana_g)
print(max(abs(num_g - ana_g)))
