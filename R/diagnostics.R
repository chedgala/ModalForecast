#' Diagnostic Plots for Modal ARIMA Models
#'
#' @description
#' Provides visual and statistical diagnostics for the residuals of a fitted modal ARIMA model.
#' Produces a comprehensive diagnostic panel including time series plot with
#' fitted modes, ACF/PACF of residuals, QQ-plot for normality, histogram of
#' residuals, and Ljung-Box p-values, all implemented using ggplot2.
#'
#' @param object An object of class \code{modal_arima}.
#' @param ... Additional arguments (unused).
#' @return A list of ggplot objects (invisibly) and draws the panel.
#' @import ggplot2
#' @importFrom stats acf pacf shapiro.test Box.test qqnorm qqline predict is.ts ts start frequency dnorm qnorm pnorm time residuals
#' @importFrom scales comma_format
#' @importFrom grid textGrob gpar
#' @importFrom gridExtra grid.arrange
#' @export
#'
#' @examples
#' \donttest{
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
#' }
diagnostics <- function(object, ...) {
  UseMethod("diagnostics")
}

#' @rdname diagnostics
#' @export
diagnostics.modal_arima <- function(object, ...) {
  # Global variable hack for ggplot2
  Time <- Observed <- Fitted <- Lag <- ACF <- PACF <- Residuals <- Pval <- density <- NULL

  y_orig <- object$y
  order <- object$order
  coefs <- object$coefficients

  p <- order[1]
  d <- order[2]
  q <- order[3]

  y_diff <- y_orig
  if (d > 0) y_diff <- diff(y_orig, differences = d)
  n <- length(y_diff)

  c_mu <- coefs["intercept"]
  phi <- if (p > 0) coefs[paste0("ar", 1:p)] else numeric(0)
  theta <- if (q > 0) coefs[paste0("ma", 1:q)] else numeric(0)
  sigma_hat <- coefs["sigma"]
  gamma_hat <- coefs["gamma"]

  # Reconstruct fitted modes and residuals
  mu_t <- numeric(n)
  eps <- numeric(n)
  mean_y <- mean(y_diff)
  for (t in 1:n) {
    ar_term <- 0
    if (p > 0) {
      for (i in 1:p) {
        ar_term <- ar_term + phi[i] * (if (t - i > 0) y_diff[t - i] else mean_y)
      }
    }
    ma_term <- 0
    if (q > 0) {
      for (j in 1:q) {
        ma_term <- ma_term + theta[j] * (if (t - j > 0) eps[t - j] else 0)
      }
    }
    mu_t[t] <- c_mu + ar_term + ma_term
    eps[t] <- y_diff[t] - mu_t[t]
  }

  # Randomized Quantile Residuals via CDF of the skew-normal
  p_skewnorm <- function(y, mu, sigma, gamma) {
    z <- (y - mu) / sigma
    p_thr <- 1 / (gamma^2 + 1)
    cdf <- numeric(length(y))
    idx_lt <- which(y < mu)
    if (length(idx_lt) > 0) cdf[idx_lt] <- 2 * p_thr * stats::pnorm(z[idx_lt] * gamma)
    idx_ge <- which(y >= mu)
    if (length(idx_ge) > 0) cdf[idx_ge] <- p_thr + 2 * (1 - p_thr) * (stats::pnorm(z[idx_ge] / gamma) - 0.5)
    return(cdf)
  }

  u <- p_skewnorm(y_diff, mu_t, sigma_hat, gamma_hat)
  u <- pmax(1e-7, pmin(1 - 1e-7, u))
  rqr <- stats::qnorm(u)

  # Prepare data for ggplot
  time_idx <- if(stats::is.ts(y_diff)) as.numeric(stats::time(y_diff)) else 1:n
  df_fits <- data.frame(
    Time = time_idx,
    Observed = as.numeric(y_diff),
    Fitted = as.numeric(mu_t)
  )

  # 1. Series & Fitted Modes
  p1 <- ggplot(df_fits, aes(x = Time)) +
    geom_line(aes(y = Observed, color = "Observed"), alpha = 0.6) +
    geom_line(aes(y = Fitted, color = "Fitted Mode"), linewidth = 0.75) +
    scale_color_manual(values = c("Observed" = "gray40", "Fitted Mode" = "blue")) +
    labs(title = "Series & Fitted Modes", y = "Value", x = "Time", color = NULL) +
    theme_minimal(base_size = 7) +
    theme(legend.position = "top")

  # 2. ACF of RQR
  bacf <- stats::acf(rqr, plot = FALSE, lag.max = 20)
  df_acf <- data.frame(Lag = as.numeric(bacf$lag[-1]), ACF = as.numeric(bacf$acf[-1]))
  p2 <- ggplot(df_acf, aes(x = Lag, y = ACF)) +
    geom_bar(stat = "identity", width = 0.1, fill = "blue") +
    geom_hline(yintercept = c(-1.96/sqrt(n), 1.96/sqrt(n)), linetype = "dashed", color = "red") +
    labs(title = "ACF of Residuals (RQR)") +
    theme_minimal(base_size = 7)

  # 3. PACF of RQR
  bpacf <- stats::pacf(rqr, plot = FALSE, lag.max = 20)
  df_pacf <- data.frame(Lag = as.numeric(bpacf$lag), PACF = as.numeric(bpacf$acf))
  p3 <- ggplot(df_pacf, aes(x = Lag, y = PACF)) +
    geom_bar(stat = "identity", width = 0.1, fill = "blue") +
    geom_hline(yintercept = c(-1.96/sqrt(n), 1.96/sqrt(n)), linetype = "dashed", color = "red") +
    labs(title = "PACF of Residuals (RQR)") +
    theme_minimal(base_size = 7)

  # 4. QQ-Plot
  df_rqr <- data.frame(Residuals = rqr)
  p4 <- ggplot(df_rqr, aes(sample = Residuals)) +
    stat_qq(color = "steelblue", alpha = 0.6) +
    stat_qq_line(color = "red", linewidth = 0.8) +
    labs(title = "Normal QQ-Plot (RQR)") +
    theme_minimal(base_size = 7)

  # 5. Histogram
  p5 <- ggplot(df_rqr, aes(x = Residuals)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20, fill = "steelblue", alpha = 0.4, color = "white") +
    stat_function(fun = stats::dnorm, color = "red", linewidth = 0.8) +
    labs(title = "Histogram (RQR)") +
    theme_minimal(base_size = 7)

  # 6. Ljung-Box p-values
  max_lag <- min(20, n - 1)
  lb_pvals <- sapply(1:max_lag, function(lag) {
    stats::Box.test(rqr, lag = lag, type = "Ljung-Box")$p.value
  })
  df_lb <- data.frame(Lag = 1:max_lag, Pval = lb_pvals)
  p6 <- ggplot(df_lb, aes(x = Lag, y = Pval)) +
    geom_point(color = "blue", size = 2) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.8) +
    ylim(0, 1) +
    labs(title = "Ljung-Box p-values", y = "p-value") +
    theme_minimal(base_size = 7)

  # Combine plots using gridExtra
  title_grob <- grid::textGrob("Modal ARIMA Diagnostic Panel",
                               gp = grid::gpar(fontsize = 10, fontface = "bold"))
  gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3, top = title_grob)

  # Print formal test results
  cat("\n=== Diagnostic Tests ===\n")
  sw <- stats::shapiro.test(rqr)
  cat(sprintf("Shapiro-Wilk Normality Test: W = %.4f, p-value = %.4f\n", sw$statistic, sw$p.value))
  if (sw$p.value > 0.05) {
    cat("  -> Conclusion (alpha = 0.05): Do not reject the assumption of normality.\n")
    cat("  -> The Randomized Quantile Residuals can be considered normally distributed.\n")
  } else {
    cat("  -> Conclusion (alpha = 0.05): Reject the assumption of normality.\n")
    cat("  -> The residuals deviate significantly from a normal distribution.\n")
  }

  lb <- stats::Box.test(rqr, lag = 10, type = "Ljung-Box")
  cat(sprintf("\nLjung-Box Test (lag=10):     X2 = %.4f, p-value = %.4f\n", lb$statistic, lb$p.value))
  if (lb$p.value > 0.05) {
    cat("  -> Conclusion (alpha = 0.05): Do not reject the assumption of independence.\n")
    cat("  -> The sample residuals can be considered independent and identically distributed.\n")
  } else {
    cat("  -> Conclusion (alpha = 0.05): Reject the assumption of independence.\n")
    cat("  -> The sample residuals exhibit significant autocorrelation and structure.\n")
  }

  cat(sprintf("\nEstimated gamma (skewness):  %.4f\n", gamma_hat))
  gamma_se <- tryCatch(sqrt(diag(solve(object$hessian))[length(object$coefficients)]), error=function(e) NA)
  
  if (!is.na(gamma_se)) {
      # Test H0: gamma = 1
      z_gamma <- (gamma_hat - 1) / gamma_se
      p_gamma <- 2 * (1 - stats::pnorm(abs(z_gamma)))
  } else {
      p_gamma <- NA
  }
  
  if (!is.na(p_gamma) && p_gamma > 0.05) {
      cat("  -> Conclusion (alpha = 0.05): Do not reject the assumption of symmetry (gamma = 1).\n")
      cat("  -> The distribution can be considered symmetric (Gaussian behavior).\n\n")
  } else if (!is.na(p_gamma)) {
      cat("  -> Conclusion (alpha = 0.05): Reject the assumption of symmetry.\n")
      if (gamma_hat > 1) {
          cat("  -> The distribution is significantly right-skewed.\n\n")
      } else {
          cat("  -> The distribution is significantly left-skewed.\n\n")
      }
  } else {
      cat("  -> Unable to compute standard error for gamma. Cannot test symmetry formally.\n\n")
  }

  invisible(list(p1, p2, p3, p4, p5, p6))
}

#' Extract Residuals from Modal ARIMA models
#'
#' @param object An object of class \code{modal_arima}.
#' @param ... Additional arguments.
#' @return A numeric vector of residuals.
#' @export
residuals.modal_arima <- function(object, ...) {
  # For internal use, we just need residuals
  y_orig <- object$y
  order <- object$order
  coefs <- object$coefficients

  p <- order[1]
  d <- order[2]
  q <- order[3]

  y_diff <- y_orig
  if (d > 0) y_diff <- diff(y_orig, differences = d)
  n <- length(y_diff)

  c_mu <- coefs["intercept"]
  phi <- if (p > 0) coefs[paste0("ar", 1:p)] else numeric(0)
  theta <- if (q > 0) coefs[paste0("ma", 1:q)] else numeric(0)

  mu_t <- numeric(n)
  eps <- numeric(n)
  mean_y <- mean(y_diff)
  for (t in 1:n) {
    ar_term <- 0
    if (p > 0) {
      for (i in 1:p) {
        ar_term <- ar_term + phi[i] * (if (t - i > 0) y_diff[t - i] else mean_y)
      }
    }
    ma_term <- 0
    if (q > 0) {
      for (j in 1:q) {
        ma_term <- ma_term + theta[j] * (if (t - j > 0) eps[t - j] else 0)
      }
    }
    mu_t[t] <- c_mu + ar_term + ma_term
    eps[t] <- y_diff[t] - mu_t[t]
  }
  return(eps)
}
