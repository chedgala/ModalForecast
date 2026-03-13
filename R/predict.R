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

#' Forecast methodology for Modal ARIMA
#' @param object A modal_arima object.
#' @param h The forecast horizon.
#' @param ... Additional arguments.
#' @return A forecast object.
#' @importFrom forecast forecast
#' @export
forecast.modal_arima <- function(object, h = 10, ...) {
    pred <- predict(object, n.ahead = h)
    res <- list(mean = pred, method = paste0("Modal ARIMA(", object$order[1], ",", object$order[2], ",", object$order[3], ")"), x = object$y, model = object)
    class(res) <- "forecast"
    return(res)
}
