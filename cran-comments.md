## Submission summary
This is a new submission. ModalForecast implements modal ARIMA models under the SKD family of distributions.

## Test environments
* local Windows 11 x64, R 4.5.0
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There were no NOTEs.

## Resubmission
This is a resubmission. In this version I have:
* Fixed the NOTE regarding the non-standard file `Rplots.pdf` found at top level by removing it and adding it to `.Rbuildignore` and `.gitignore`.
* Fixed the NOTE regarding the non-standard file `submit_mock.R` found at top level by adding it to `.Rbuildignore`.
* The possibly misspelled word "SKD" in the DESCRIPTION file is an acronym for "Skewed Distribution" and is spelled correctly.
