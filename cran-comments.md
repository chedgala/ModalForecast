## Submission summary
This is a new submission. ModalForecast implements modal ARIMA models under the SKD family of distributions.

## Test environments
* local Windows 11 x64, R 4.5.0
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There were no NOTEs.

## Resubmission
This is a resubmission. In this version I have addressed the following points raised by the CRAN maintainer:
* Explained the acronym ARIMA in the Description text.
* Added the reference describing the methods in the Description field (Galarza et al., 2017) using the requested format `<doi:...>`.
* Unwrapped all examples by removing `\donttest{}` and simplified the examples (reduced maximum model dimensions and bootstrap iterations) so they can be executed quickly during testing.
* Removed the `LICENSE` file and its reference from the `DESCRIPTION` file, keeping only `License: GPL-3`, as there are no additional restrictions.
* Fixed previous NOTEs regarding non-standard files (`Rplots.pdf`, `submit_mock.R`) at the top level.
* The possibly misspelled word "SKD" in the DESCRIPTION file is an acronym for "Skewed Distribution" and is spelled correctly.
