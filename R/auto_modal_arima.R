#' Automatic selection of Parametric Modal ARIMA model
#'
#' @details
#' This function performs a grid search over AR and MA orders to find the
#' optimal Modal ARIMA model based on the selected information criterion (AIC or BIC).
#'
#' @references
#' Fernandez, C. and Steel, M. F. J. (1998). On Bayesian Modeling of Fat Tails
#' and Skewness. Journal of the American Statistical Association, 93(441), 359-371.
#'
#' Galarza, C. E., Lachos, V. H., Cabral, C. R. B., and Castro, L. M. (2017).
#' Robust quantile regression using a generalized class of skewed distributions.
#' Stat, 6(1), 113-130.
#'
#' @seealso \code{\link{fit_modal_arima}}
#'
#' @note GitHub repository: \url{https://github.com/chedgala/ModalForecast}
#'
#' @param y numeric vector or time series of observations
#' @param d Integer, degree of differencing. If NA, it's determined automatically.
#' @param max.p Maximum AR order
#' @param max.q Maximum MA order
#' @param ic Information criterion to be used in model selection ("aic", "bic")
#' @param dist Character string specifying the error distribution.
#'   \code{"normal"} (default) for Skew-Normal, \code{"t"} for Skewed Student-t,
#'   \code{"laplace"} for Skewed Laplace.
#'
#' @return An object of class \code{modal_arima}.
#' @importFrom stats AIC BIC
#' @importFrom forecast ndiffs
#' @export
#'
#' @examples
#' \donttest{
#' # 1. Simulate an asymmetric AR(1) time series
#' set.seed(123)
#' y <- arima.sim(n = 200, list(ar = 0.5))
#'
#' # 2. Fit Modal ARIMA models with different distributions
#' mod_n <- fit_modal_arima(y, order = c(1, 0, 0), dist = "normal")
#' mod_t <- fit_modal_arima(y, order = c(1, 0, 0), dist = "t")
#' mod_l <- fit_modal_arima(y, order = c(1, 0, 0), dist = "laplace")
#'
#' # 3. Compare models
#' summary(mod_n)
#' summary(mod_l)
#' AIC(mod_n); AIC(mod_t); AIC(mod_l)
#'
#' # 4. Run residual diagnostics
#' diagnostics(mod_n)
#'
#' # 5. Produce forecasts with prediction bands
#' pred <- forecast(mod_n, h = 10, level = c(80, 95), interval = "asymptotic")
#' print(pred$lower)
#' }
auto.modal_arima <- function(y, d = NA, max.p = 5, max.q = 5,
                              ic = c("aic", "bic"),
                              dist = c("normal", "t", "laplace")) {
  ic <- match.arg(ic)
  dist <- match.arg(dist)

  if (is.na(d)) {
    if (!requireNamespace("forecast", quietly = TRUE)) {
      stop("Package 'forecast' is needed to automatically determine 'd'. Please install it or specify 'd'.")
    }
    d <- forecast::ndiffs(y)
  }

  best_ic <- Inf
  best_mod <- NULL

  for (p in 0:max.p) {
    for (q in 0:max.q) {
      if (p == 0 && q == 0) next
      mod <- tryCatch({
        fit_modal_arima(y, order = c(p, d, q), dist = dist)
      }, error = function(e) NULL)

      if (!is.null(mod) && mod$convergence == 0) {
        current_ic <- if (ic == "aic") AIC(mod) else BIC(mod)
        if (current_ic < best_ic) {
          best_ic <- current_ic
          best_mod <- mod
        }
      }
    }
  }

  if (is.null(best_mod)) {
    best_mod <- fit_modal_arima(y, order = c(0, d, 0), dist = dist)
  }

  return(best_mod)
}
