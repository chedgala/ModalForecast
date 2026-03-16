library(ModalForecast)

# Simulate an AR(1) time series
set.seed(123)
y <- arima.sim(n = 200, list(ar = 0.5))

# Fit a modal ARIMA(1,0,0) model
mod <- fit_modal_arima(y, order = c(1, 0, 0))
summary(mod)

# Automatically select the best modal ARIMA model
mod <- auto.modal_arima(y, max.p = 2, max.q = 2)
summary(mod)

# Extract model metrics
AIC(mod); BIC(mod)
# Extract model residuals
residuals(mod)
# Run residual diagnostics
diagnostics(mod)
# Produce forecasts
forecast(mod)


# Instala si no lo tienes
# install.packages("devtools")

devtools::document(roclets = c('rd', 'collate', 'namespace'))
devtools::check(args = "--as-cran", manual = TRUE, vignettes = TRUE)


# Servidores Windows (win-builder)
devtools::check_win_devel()   # R en desarrollo (El más importante para CRAN)
devtools::check_win_release() # R versión actual (Recomendado)

# Servidores macOS (mac-builder)
devtools::check_mac_release()
