rm(list = ls())

pkgs_req <- c("devtools", "ggplot2", "dplyr", "ggpubr", "geoR", "MBA", "knitr",
              "kableExtra")

pkgs_missing <- pkgs_req[!pkgs_req %in% installed.packages()[,"Package"]]

if (length(pkgs_missing) > 0) {
  install.packages(pkgs_missing)
}

devtools::install_github("SPan-18/spStackCOS-dev")
