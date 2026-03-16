#' @export
predict.modal_arima <- function(object, n.ahead = 10, ...) {
  y_orig <- object$y
  order <- object$order
  coefs <- object$coefficients

  p <- order[1]
  d <- order[2]
  q <- order[3]

  y_diff <- y_orig
  if (d > 0) {
    y_diff <- diff(y_orig, differences = d)
  }
  n <- length(y_diff)

  c_mu <- coefs["intercept"]
  phi <- if (p > 0) coefs[paste0("ar", 1:p)] else numeric(0)
  theta <- if (q > 0) coefs[paste0("ma", 1:q)] else numeric(0)

  # Reconstruct residuals for historical period
  mu_t <- numeric(n)
  eps <- numeric(n)
  mean_y <- mean(y_diff)
  for (t in 1:n) {
    ar_term <- 0
    if (p > 0) {
        for (i in 1:p) {
            ar_term <- ar_term + phi[i] * (if (t-i > 0) y_diff[t-i] else mean_y)
        }
    }
    ma_term <- 0
    if (q > 0) {
        for (j in 1:q) {
            ma_term <- ma_term + theta[j] * (if (t-j > 0) eps[t-j] else 0)
        }
    }
    mu_t[t] <- c_mu + ar_term + ma_term
    eps[t] <- y_diff[t] - mu_t[t]
  }

  # Forecast mode diffs out of sample
  pred_diff <- numeric(n.ahead)
  for (h in 1:n.ahead) {
    ar_term <- 0
    if (p > 0) {
        for (i in 1:p) {
            val <- if (h - i > 0) pred_diff[h - i] else y_diff[n + h - i]
            ar_term <- ar_term + phi[i] * val
        }
    }
    ma_term <- 0
    if (q > 0) {
        for (j in 1:q) {
            val <- if (h - j > 0) 0 else eps[n + h - j]
            ma_term <- ma_term + theta[j] * val
        }
    }
    pred_diff[h] <- c_mu + ar_term + ma_term
  }

  # Undifference if d > 0
  if (d > 0) {
    last_y <- utils::tail(y_orig, d)
    pred <- stats::diffinv(pred_diff, differences = d, xi = last_y)[-(1:d)]
  } else {
    pred <- pred_diff
  }

  return(pred)
}

#' @importFrom forecast forecast
#' @export
forecast::forecast

#' Forecast methodology for Modal ARIMA
#'
#' @param object A modal_arima object.
#' @param h The forecast horizon.
#' @param level Confidence level for prediction intervals.
#' @param interval Method for computing prediction intervals ("asymptotic" or "bootstrap"). Defaults to "asymptotic".
#' @param npaths Number of simulated paths for bootstrap intervals. Defaults to 1000.
#' @param ... Additional arguments.
#' @return A forecast object.
#' @export
#'
#' @examples
#' # 1. Simulate an asymmetric AR(1) time series
#' set.seed(123)
#' y <- arima.sim(n = 200, list(ar = 0.5))
#' 
#' # 2. Fit a Modal ARIMA model
#' mod <- fit_modal_arima(y, order = c(1, 0, 0))
#' 
#' # 3. Model Information and AIC/BIC
#' summary(mod)
#' AIC(mod)
#' BIC(mod)
#' 
#' # 4. Run residual diagnostics with hypothesis tests
#' diagnostics(mod)
#' 
#' # 5. Produce forecasts with Asymptotic & Bootstrap prediction bands
#' pred_asymp <- forecast(mod, h = 10, level = c(80, 95), interval = "asymptotic")
#' pred_boot <- forecast(mod, h = 10, level = c(80, 95), interval = "bootstrap", npaths = 500)
#' 
#' # The output contains the upper and lower bounds for the respective levels:
#' print(pred_asymp$lower)
forecast.modal_arima <- function(object, h = 10, level = c(80, 95), interval = c("asymptotic", "bootstrap"), npaths = 1000, ...) {
  interval <- match.arg(interval)
  
  if (length(object$y) <= 50 && interval == "asymptotic") {
    message("Notice: The sample size is small (n <= 50). It is recommended to use Bootstrap prediction intervals (interval = 'bootstrap') for more reliable results.")
  }

  pred <- predict(object, n.ahead = h)

  coefs <- object$coefficients
  p <- object$order[1]
  d <- object$order[2]
  q <- object$order[3]
  ar_coef <- if (p > 0) coefs[paste0("ar", 1:p)] else numeric(0)
  ma_coef <- if (q > 0) coefs[paste0("ma", 1:q)] else numeric(0)
  sigma <- coefs["sigma"]
  gamma <- coefs["gamma"]
  c_mu <- coefs["intercept"]

  lower <- matrix(NA, nrow = h, ncol = length(level))
  upper <- matrix(NA, nrow = h, ncol = length(level))

  if (interval == "asymptotic") {
    # Handle integration (d > 0) to get effective AR polynomial
    poly_diff <- c(1)
    if (d > 0) {
      for (i in seq_len(d)) {
        poly_diff <- stats::convolve(poly_diff, rev(c(1, -1)), type = "open")
      }
    }
    ar_poly <- c(1, -ar_coef)
    total_ar <- stats::convolve(ar_poly, rev(poly_diff), type = "open")
    total_ar_coef <- -total_ar[-1]
    
    total_ar_coef <- total_ar_coef[total_ar_coef != 0]

    if (h > 1) {
      psi <- stats::ARMAtoMA(ar = total_ar_coef, ma = ma_coef, lag.max = h - 1)
      psi_weights <- c(1, psi)
    } else {
      psi_weights <- c(1)
    }

    sigma_h <- sigma * sqrt(cumsum(psi_weights^2))

    # Skew-Normal Quantile Function
    qsn <- function(prob, mu = 0, sigma = 1, gamma = 1) {
      p_mu <- 1 / (gamma^2 + 1)
      res <- numeric(length(prob))
      for (i in seq_along(prob)) {
        if (prob[i] < p_mu) {
          z <- (1 / gamma) * stats::qnorm(prob[i] * (1 + gamma^2) / 2)
        } else {
          val <- (prob[i] * (1 + gamma^2) + gamma^2 - 1) / (2 * gamma^2)
          z <- gamma * stats::qnorm(val)
        }
        res[i] <- mu + sigma[i] * z
      }
      return(res)
    }

    for (i in seq_along(level)) {
      alpha <- 1 - level[i] / 100
      p_lower <- alpha / 2
      p_upper <- 1 - alpha / 2

      q_low <- qsn(rep(p_lower, h), mu = 0, sigma = sigma_h, gamma = gamma)
      q_up <- qsn(rep(p_upper, h), mu = 0, sigma = sigma_h, gamma = gamma)

      lower[, i] <- pred + q_low
      upper[, i] <- pred + q_up
    }

  } else if (interval == "bootstrap") {
    
    # Skew-Normal Random Generation
    rsn <- function(n, mu = 0, sigma = 1, gamma = 1) {
      u <- stats::runif(n)
      p_mu <- 1 / (gamma^2 + 1)
      z <- numeric(n)
      idx_lower <- which(u < p_mu)
      idx_upper <- which(u >= p_mu)
      
      if (length(idx_lower) > 0) {
        z[idx_lower] <- (1 / gamma) * stats::qnorm(u[idx_lower] * (1 + gamma^2) / 2)
      }
      if (length(idx_upper) > 0) {
        val <- (u[idx_upper] * (1 + gamma^2) + gamma^2 - 1) / (2 * gamma^2)
        z[idx_upper] <- gamma * stats::qnorm(val)
      }
      return(mu + sigma * z)
    }

    y_orig <- object$y
    y_diff <- y_orig
    if (d > 0) y_diff <- diff(y_orig, differences = d)
    n_hist <- length(y_diff)

    # Reconstruct historical residuals to bootstrap MA components
    mean_y <- mean(y_diff)
    eps_hist <- numeric(n_hist)
    for (t in 1:n_hist) {
      ar_term <- 0
      if (p > 0) {
        for (i in seq_len(p)) ar_term <- ar_term + ar_coef[i] * (if (t-i > 0) y_diff[t-i] else mean_y)
      }
      ma_term <- 0
      if (q > 0) {
        for (j in seq_len(q)) ma_term <- ma_term + ma_coef[j] * (if (t-j > 0) eps_hist[t-j] else 0)
      }
      eps_hist[t] <- y_diff[t] - (c_mu + ar_term + ma_term)
    }
    
    sim_paths <- matrix(NA, nrow = h, ncol = npaths)

    for (b in 1:npaths) {
      yt_sim <- numeric(h)
      eps_sim <- rsn(h, mu = 0, sigma = sigma, gamma = gamma)
      
      for (t in 1:h) {
        ar_term <- 0
        if (p > 0) {
          for (i in 1:p) {
            val <- if (t - i > 0) yt_sim[t - i] else y_diff[n_hist + t - i]
            ar_term <- ar_term + ar_coef[i] * val
          }
        }
        ma_term <- 0
        if (q > 0) {
          for (j in 1:q) {
            val <- if (t - j > 0) eps_sim[t - j] else eps_hist[n_hist + t - j]
            ma_term <- ma_term + ma_coef[j] * val
          }
        }
        yt_sim[t] <- c_mu + ar_term + ma_term + eps_sim[t]
      }
      
      # Undifference simulated path
      if (d > 0) {
        last_y <- utils::tail(y_orig, d)
        path <- stats::diffinv(yt_sim, differences = d, xi = last_y)[-(1:d)]
      } else {
        path <- yt_sim
      }
      
      sim_paths[, b] <- path
    }

    for (i in seq_along(level)) {
      alpha <- 1 - level[i] / 100
      p_lower <- alpha / 2
      p_upper <- 1 - alpha / 2
      
      lower[, i] <- apply(sim_paths, 1, stats::quantile, probs = p_lower)
      upper[, i] <- apply(sim_paths, 1, stats::quantile, probs = p_upper)
    }
  }

  res <- list(
    mean   = pred,
    lower  = lower,
    upper  = upper,
    level  = level,
    method = paste0("Modal ARIMA(", p, ",", d, ",", q, ")"),
    x      = object$y,
    model  = object
  )

  class(res) <- "forecast"
  return(res)
}
