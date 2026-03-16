#' @export
print.modal_arima <- function(x, ...) {
  dist <- if (!is.null(x$dist)) x$dist else "normal"
  cat("\nCall:\nfit_modal_arima(order = c(", paste(x$order, collapse=", "), "), dist = \"", dist, "\")\n", sep="")
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\nLog-likelihood:", round(x$loglik, 2), "\n")
  invisible(x)
}

#' @export
coef.modal_arima <- function(object, ...) {
  object$coefficients
}

#' @export
logLik.modal_arima <- function(object, ...) {
  structure(object$loglik, df = length(object$coefficients), nobs = length(object$y), class = "logLik")
}

#' @export
AIC.modal_arima <- function(object, ..., k = 2) {
  -2 * object$loglik + k * length(object$coefficients)
}

#' @export
BIC.modal_arima <- function(object, ...) {
  -2 * object$loglik + log(length(object$y)) * length(object$coefficients)
}

#' @export
summary.modal_arima <- function(object, ...) {
  dist <- if (!is.null(object$dist)) object$dist else "normal"

  # Simple asymptotic standard errors from hessian
  se <- tryCatch(sqrt(diag(solve(object$hessian))), error = function(e) rep(NA, length(object$coefficients)))
  z <- object$coefficients / se
  pval <- 2 * (1 - pnorm(abs(z)))

  res <- cbind(Estimate = object$coefficients, `Std. Error` = se, `z value` = z, `Pr(>|z|)` = pval)

  dist_label <- switch(dist, "normal"="Skew-Normal", "t"="Skewed Student-t", "laplace"="Skewed Laplace")
  cat("\nModal ARIMA(", paste(object$order, collapse=","), ") Model [", dist_label, "]\n", sep="")
  cat("====================================================\n")
  printCoefmat(res)
  cat("---\n")
  cat(sprintf("Log-likelihood: %.2f   AIC: %.2f   BIC: %.2f\n", object$loglik, AIC(object), BIC(object)))
}

#' @export
plot.modal_arima <- function(x, ...) {
  Time <- Value <- NULL # Hack for ggplot2 notes
  dist <- if (!is.null(x$dist)) x$dist else "SKN"
  df <- data.frame(Time = 1:length(x$y), Value = as.numeric(x$y))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Time, y = Value)) +
    ggplot2::geom_line(color = "gray40") +
    ggplot2::labs(title = paste0("Modal ARIMA(", paste(x$order, collapse=","), ") [", dist, "] Fit"),
                  x = "Time", y = "Value") +
    ggplot2::theme_minimal()
  print(p)
}
