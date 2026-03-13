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
#' @seealso \code{\link{fit_modal_arima}}
#'
#' @note GitHub repository: \url{https://github.com/chedgala/ModalForecast}
#'
#' @param y numeric vector or time series of observations
#' @param d Integer, degree of differencing. If NA, it's determined automatically.
#' @param max.p Maximum AR order
#' @param max.q Maximum MA order
#' @param ic Information criterion to be used in model selection ("aic", "bic")
#'
#' @return An object of class \code{modal_arima}.
#' @importFrom stats AIC BIC
#' @importFrom forecast ndiffs
#' @export
#'
#' @examples
#' set.seed(123)
#' y <- arima.sim(n = 200, list(ar = 0.5))
#' mod <- auto.modal.arima(y, d = 0, max.p = 2, max.q = 2)
#' summary(mod)
auto.modal.arima <- function(y, d = NA, max.p = 5, max.q = 5, ic = c("aic", "bic")) {
  ic <- match.arg(ic)
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
        fit_modal_arima(y, order = c(p, d, q))
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
    # Fallback to white noise modal representation if no ARIMA fits
    best_mod <- fit_modal_arima(y, order = c(0, d, 0))
  }
  
  return(best_mod)
}
