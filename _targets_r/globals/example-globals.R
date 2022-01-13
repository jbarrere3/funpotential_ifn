options(tidyverse.quiet = TRUE)
options(stringsAsFactors=FALSE)
source("R/functions_data.R")
source("R/functions_plot.R")
source("R/functions_analysis.R")
tar_option_set(packages = c("dplyr", "ggplot2", "targets", "tidyr", "RColorBrewer", "lme4",
          "data.table", "knitr", "stringr", "measurements", "sf", "raster",
          "rgdal", "exactextractr", "rgeos", "rnaturalearth", "rnaturalearthdata", 
          "ggspatial", "cowplot", "rjags", "coda", "R2jags", "MASS", "ggmcmc"))
