
required_packages <- c("matrixStats",
                       "dplyr",
                       "minpack.lm",
                       "ggplot2",
                       "lme4",
                       "gridExtra",
                       "toOrdinal",
                       "reshape2"
                       "runner"
                       )


new.packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)