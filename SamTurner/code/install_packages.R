## Script: install_packages.R
## Author: Sam Turner sat19@ic.ac.uk
## About: Installs required R packages if they are missing


required_packages <- c("matrixStats",
                       "dplyr",
                       "minpack.lm",
                       "ggplot2",
                       "lme4",
                       "gridExtra",
                       "toOrdinal",
                       "reshape2"
                       )


new.packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)