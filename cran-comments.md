## Submission summary

This is an update of ModalForecast from version 0.1.0 to 0.2.0. It adds seasonal
(SARIMA) models, marginal modal forecasts and exact prediction intervals, and
fixes several bugs in the estimation of the Skewed Student-t and Skewed Laplace
models (see NEWS.md). The package title and description were updated to
mention the seasonal models, and two contributors were added.

## Test environments

* local Windows 11 x64, R 4.5.0
* win-builder, R-devel (2026-09-21 r90579 ucrt)
* win-builder, R-release (R 4.6.1)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Possibly misspelled words in DESCRIPTION: SARIMA

  "SARIMA" is the standard acronym for seasonal ARIMA and is spelled
  correctly; the Description introduces it as "seasonal ARIMA (SARIMA)".

## Reverse dependencies

There are no reverse dependencies.
