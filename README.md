# Bayesian predictive stacking for spatially-temporally misaligned outcome and exposure data

>[!IMPORTANT]
>Install the R package **spStackCOS** from the GitHub repository [SPan-18/spStackCOS-dev](https://github.com/SPan-18/spStackCOS-dev).

***

This repository contains code to implement different analyses, as it appears in the manuscript:

Soumyakanti Pan, and Sudipto Banerjee. 2025. _Bayesian inference for spatially-temporally misaligned data using preditive stacking._

Implemented models include Bayesian spatial-temporal regression models under the change of support, and the modifiable areal unit problem.
The article also analyzes the effect of ozone on annual age-adjusted asthma-related emergency department (ED) visit rates among residents of California during the study period 2015-22. 
A key challenge is that the outcome and the exposure data are misaligned both spatially and temporally, with the asthma-related ED visit rates recorded annually and at the county-level, while
the ozone concentrations are measured hourly by 200 monitoring stations scattered across different locations inside California.

## Authors

| Name   | Email       |              |
|:------ |:----------- | :----------- |
| Soumyakanti Pan | span18@ucla.edu | PhD candidate, UCLA Biostatistics |

## Instructions

To reproduce results or to try out the example code, issue
```bash
git clone https://github.com/SPan-18/AsthmaOzoneCA
```
on your personal computer or, remote server to clone this repository.

Install the R package `spStackCOS` package from GitHub.
```r
devtools::install_github("SPan-18/spStackCOS-dev")
```

Knit the markdown document `vignette.Rmd`.


