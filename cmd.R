#' cmd.R
#' @description Script to source the code compilation rapidly from command line 
#'              To compile, go to 'funpotential_ifn' directory and run source('cmd.R')
Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc")
rmarkdown::render("index.Rmd")