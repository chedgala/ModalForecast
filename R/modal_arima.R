#' Fit a Parametric Modal ARIMA Model using the Skew Normal Distribution
#'
#' This framework is based on the Fernandez-Steel (1998) skew-normal 
#' distribution, which allows for a flexible error structure where the 
#' mode is explicitly modeled.
#'
#' @references
#' Fernandez, C. and Steel, M. F. J. (1998). On Bayesian Modeling of Fat Tails 
#' and Skewness. Journal of the American Statistical Association, 93(441), 359-371.
#'
#' @seealso \code{\link{auto.modal.arima}}
#'
#' @note GitHub repository: \url{https://github.com/chedgala/ModalForecast}
#'
#' @param y numeric vector or time series of observations
#' @param order A specification of the non-seasonal part of the ARIMA
#'   model: the three components (p, d, q) are the AR order, the
#'   degree of differencing, and the MA order.
#'
#' @return An object of class \code{modal_arima}.
#' @importFrom stats arima optim qnorm runif dnorm pnorm coef AIC BIC logLik printCoefmat
#' @importFrom graphics plot
#' @export
#'
#' @examples
#' set.seed(123)
#' y <- arima.sim(n = 200, list(ar = 0.5))
#' mod <- fit_modal_arima(y, order = c(1, 0, 0))
#' summary(mod)
fit_modal_arima <- function(y, order = c(1, 0, 0)) {
  if (length(order) != 3) stop("'order' must have length 3 (p, d, q)")
  
  p <- order[1]
  d <- order[2]
  q <- order[3]
  
  y_orig <- y
  if (d > 0) {
    y <- diff(y, differences = d)
  }
  
  n <- length(y)
  
  # Log-likelihood function for Fernandez-Steel Skew Normal
  # params: c, phi (p), theta (q), log_sigma, log_gamma
  neg_log_lik <- function(params) {
    c_mu <- params[1]
    phi <- if (p > 0) params[2:(p+1)] else numeric(0)
    theta <- if (q > 0) params[(p+2):(1+p+q)] else numeric(0)
    
    log_sigma <- params[length(params) - 1]
    log_gamma <- params[length(params)]
    sigma <- exp(log_sigma)
    gamma <- exp(log_gamma)
    
    # Stationarity and Invertibility checks
    if (p > 0) {
      ar_roots <- polyroot(c(1, -phi))
      if (any(abs(ar_roots) <= 1.001)) return(1e10)
    }
    if (q > 0) {
      ma_roots <- polyroot(c(1, theta))
      if (any(abs(ma_roots) <= 1.001)) return(1e10)
    }
    
    # Initialize recursive values
    mu_t <- numeric(n)
    eps <- numeric(n)
    
    # Simple unconditional mean approximation for burn-in
    mean_y <- mean(y)
    
    for (t in 1:n) {
      # Compute AR part
      ar_term <- 0
      for (i in seq_len(p)) {
        if (t - i > 0) {
          ar_term <- ar_term + phi[i] * y[t - i]
        } else {
          ar_term <- ar_term + phi[i] * mean_y # startup approximation
        }
      }
      
      # Compute MA part
      ma_term <- 0
      for (j in seq_len(q)) {
        if (t - j > 0) {
          ma_term <- ma_term + theta[j] * eps[t - j]
        }
      }
      
      mu_t[t] <- c_mu + ar_term + ma_term
      eps[t] <- y[t] - mu_t[t]
    }
    
    # Likelihood computation for Split / Skew Normal where mu is the mode
    # f(y) = 2 / (sigma * (gamma + 1/gamma)) * (phi((y-mu)/(sigma * gamma)) I(y >= mu) + phi((y-mu)*gamma/sigma) I(y < mu))
    
    # Standard normal density phi(z) = (2 * pi)^(-1/2) * exp(-0.5 * z^2)
    # Log density:
    
    z <- (y - mu_t) / sigma
    
    log_f0 <- function(val) {
      -0.5 * log(2 * pi) - 0.5 * val^2
    }
    
    idx_ge <- which(y >= mu_t)
    idx_lt <- which(y < mu_t)
    
    log_lik <- numeric(n)
    if (length(idx_ge) > 0) {
      log_lik[idx_ge] <- log(2) - log_sigma - log(gamma + 1/gamma) + log_f0(z[idx_ge] / gamma)
    }
    if (length(idx_lt) > 0) {
      log_lik[idx_lt] <- log(2) - log_sigma - log(gamma + 1/gamma) + log_f0(z[idx_lt] * gamma)
    }
    
    ans <- -sum(log_lik)
    if (is.na(ans) || is.infinite(ans)) {
      ans <- 1e10
    }
    return(ans)
  }
  
  # Initial values based on standard ARIMA
  init_fit <- suppressWarnings(stats::arima(y, order = c(p, 0, q), method = "CSS"))
  init_c <- if ("intercept" %in% names(init_fit$coef)) init_fit$coef["intercept"] else 0
  init_phi <- if (p > 0) init_fit$coef[paste0("ar", 1:p)] else numeric(0)
  init_theta <- if (q > 0) init_fit$coef[paste0("ma", 1:q)] else numeric(0)
  
  init_sigma <- sqrt(init_fit$sigma2)
  if (init_sigma < 1e-4) init_sigma <- 1
  
  init_params <- c(init_c, init_phi, init_theta, log(init_sigma), 0)
  
  opt <- optim(par = init_params, fn = neg_log_lik, method = "BFGS", hessian = TRUE)
  
  est_c <- opt$par[1]
  est_phi <- if (p > 0) opt$par[2:(p+1)] else numeric(0)
  est_theta <- if (q > 0) opt$par[(p+2):(1+p+q)] else numeric(0)
  est_sigma <- exp(opt$par[length(opt$par) - 1])
  est_gamma <- exp(opt$par[length(opt$par)])
  
  coef_names <- c("intercept")
  if (p > 0) coef_names <- c(coef_names, paste0("ar", 1:p))
  if (q > 0) coef_names <- c(coef_names, paste0("ma", 1:q))
  coef_names <- c(coef_names, "sigma", "gamma")
  
  coefficients <- c(est_c, est_phi, est_theta, est_sigma, est_gamma)
  names(coefficients) <- coef_names
  
  # Return fitted object
  out <- list(
    y = y_orig,
    order = order,
    coefficients = coefficients,
    loglik = -opt$value,
    hessian = opt$hessian,
    convergence = opt$convergence
  )
  
  class(out) <- "modal_arima"
  return(out)
}
