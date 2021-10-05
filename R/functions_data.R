#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_data.R  
#' @description R script containing all functions relative to data
#               importation and formatting
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#### Section 1 - Read French NFI ####
#' @description Group of functions used to read French NFI data
#' @author Georges Kunstler (INRAE - LESSEM)

#' Read IFN data for each year
#'
#' @description This is the base function to process the different types of IFN
#' data. The IFN data must already be downloaded and can be zipped or unzipped.
#' @export
#' @param variable Which variable to read
#' @param path Path to data folder
#' @param years Vector with years for which data are to be read
#' @param addyear boolean; add year to dataframe
#' @param zipped boolean if the files are zipped
#' @param ... Further arguments to @seealso [data.table::fread()].
#' @details The variable must be either: "placettes_foret", "arbres_foret", "plots_forest_5"
#'"couverts_foret", "arbres_morts_foret", "ecologie", "documentation"
#' @return a data.table object
read_ifn_data <- function(variable, path, years,
                          addyear = TRUE, zipped = TRUE, ...) {
  if (zipped) {
    zipfiles <- list.files(path, pattern = "\\.zip$")
    zipfiles <- file.path(path, zipfiles[sapply(years, grep, x = zipfiles)])
    files <- file.path(paste0(years, "-fr"), paste0(variable, "_", years, ".csv"))
    connections <- mapply(unz, zipfiles, files, SIMPLIFY = FALSE)
  } else {
    connections <- file.path(path, years,
                             paste0(variable, "_", years, ".csv"))
  }
  data_list <- lapply(connections, utils::read.table,
                      header = TRUE, sep = ";", dec = ".",
                      encoding = "UTF-8",
                      stringsAsFactors = FALSE)
  
  if (addyear) {
    for (i in seq_along(data_list)) {
      names(data_list[[i]]) <- gsub("\\.", "", gsub("X", "", names(data_list[[i]])))
      data_list[[i]]$year = years[[i]]
    }
  }
  data.table::rbindlist(data_list, use.names = TRUE, fill = TRUE)
}

#' Read IFN data remeasure for each year
#'
#' @description This is the base function to process the different types of IFN only for remeasure
#' data. The IFN data must already be downloaded and can be zipped or unzipped.
#' @export
#' @param variable Which variable to read
#' @param path Path to data folder
#' @param years Vector with years for which data are to be read in the format 2005 then we will read remeasured data for the year and the year + 5 such as 2005-2010
#' @param addyear boolean; add year to dataframe
#' @param zipped boolean if the files are zipped
#' @param ... Further arguments to @seealso [data.table::fread()].
#' @details The variable must be either:  "plots_forest_5"
#' @return a data.table object
read_ifn_remeasure_data <- function(variable, path, years,
                                    addyear = TRUE, zipped = TRUE, ...) {
  if (zipped) {
    zipfiles <- list.files(path, pattern = "\\.zip$")
    zipfiles <- file.path(path, zipfiles[sapply(
      paste(years, years+5, sep = "-"), 
      grep, x = zipfiles)])
    files <- file.path(paste0(paste(years, years+5, sep = "-"), 
                              "-en"), 
                       paste0(variable, "_", years, ".csv"))
    connections <- mapply(unz, zipfiles, files, SIMPLIFY = FALSE)
  } else {
    connections <- file.path(path, years,
                             paste0(variable, "_", years, ".csv"))
  }
  data_list <- lapply(connections, utils::read.table,
                      header = TRUE, sep = ";", dec = ".",
                      encoding = "UTF-8",
                      stringsAsFactors = FALSE)
  
  if (addyear) {
    for (i in seq_along(data_list)) {
      names(data_list[[i]]) <- gsub("\\.", "", gsub("X", "", names(data_list[[i]])))
      data_list[[i]]$year = years[[i]]
    }
  }
  data.table::rbindlist(data_list, use.names = TRUE, fill = TRUE)
}
# read_ifn_remeasure_data(variable = "plots_forest_5", path = "data/FrenchNFI_remeasure", years= 2005:2014)


#' Read plot data for desired years
#'
#' @description Function to call the base calling function for the IFN plot
#' data
#' @export
#' @param path Path to data folder
#' @param years Vector with years for which data are to be read; defaults
#' are 2008 to 2013
#' @return a data.table object
read_plot <- function(path = "data", years = 2008:2014) {
  read_ifn_data("placettes_foret", path, years,
                
                # colClasses currently ignored. Currently specified in format
                # required for data.table. Overrides boolean stringsAsFactors
                colClasses = list(character = c("dep", "ser",
                                                "uta1", "uta2", "dc",
                                                "dcespar1", "dcespar2",
                                                "tplant",
                                                "tpespar1", "tpespar2")))
}

#' Read tree level data for each year
#'
#' @description Function to call the base calling function for the IFN tree
#' level data for each year.
#' @export
#' @inheritParams read_plot
#' @return a data.table object
read_tree <- function(path = "data", years = 2008:2014) {
  read_ifn_data("arbres_foret", path, years)
}

#' Read cover data for each year
#'
#' @description Function to call the base calling function for the IFN cover
#' data for each year.
#' @export
#' @inheritParams read_plot
#' @return a data.table object
read_couv <- function(path = "data", years = 2008:2014) {
  read_ifn_data("couverts_foret", path, years)
}

#' Read dead tree data for each year
#'
#' @description Function to call the base calling function for the IFN dead
#' dead tree data for each year.
#' @export
#' @inheritParams read_plot
#' @return a data.table object
read_dead_tree <- function(path = "data", years = 2008:2014) {
  read_ifn_data("arbres_morts_foret", path, years)
}

#' Read documentation data for each year
#'
#' @description Function to call the base calling function for the IFN
#' documentation data for each year.
#' @export
#' @inheritParams read_plot
#' @return a data.table object
read_doc <- function(path = "data", years = 2008:2014) {
  data_doc <- read_ifn_data("documentation", path, years)
  data_doc$code <- ifelse(data_doc$code=="2,50E+04", "25E3", data_doc$code)
  data_doc$code <- ifelse(data_doc$code=="2,50E+4", "25E3", data_doc$code)
  data_doc$code <- ifelse(data_doc$code=="2,50E+06", "25E5", data_doc$code)
  data_doc$code <- ifelse(data_doc$code=="2,50E+6", "25E5", data_doc$code)
  return(data_doc)
  
}

#' Read ecological data (soil and abiotic variables) for each year
#'
#' @description Function to call the base calling function for the IFN
#' ecological data (soil and abiotic variables) for each year.
#' @export
#' @inheritParams read_plot
#' @return a data.table object
read_ecological_data <- function(path = "data", years = 2008:2014) {
  read_ifn_data("ecologie", path, years)
}

#' Read plot remeasure data for each year
#'
#' @description Function to call the base calling function for the IFN remeasure
#' plot remeasure data for each year.
#' @export
#' @inheritParams read_plot
#' @return a data.table object
read_reameasure_plot <- function(path = "data", years = 2008:2014) {
  read_ifn_remeasure_data("plots_forest_5", path, years)
}


#' Merge data on alive and dead trees
#'
#' @description Merge the data.tables returned by read_tree and read_dead_tree
#' into a single data table creating a new binary variable for alive/dead: 'dead'
#' @export
#' @param data_tree Data on living trees
#' @param data_tree_dead Data on dead trees
#' @return a data.table object
merge_dead_alive <- function(data_tree, data_tree_dead) {
  data.table::set(data_tree, j = "dead", value = 0)
  data.table::set(data_tree_dead, j = "dead", value = 1)
  data.table::rbindlist(list(data_tree, data_tree_dead),
                        use.names = TRUE, fill = TRUE)
}

#' Map plot data
#'
#' @description Plots the IFN data location for the data passed over the
#' French map
#' @export
#' @param data_plot IFN data to plot in the form returned by this
#' package functions
map_plot <- function(data_plot) {
  ### plot of map of plots over FRANCE
  x <- cbind(data_plot$xl93, data_plot$yl93)
  data <- as.data.frame(x)
  names(data) <- c("X", "Y")
  colors <- grDevices::densCols(x)
  graphics::plot(x, col = colors, pch = 20, cex = 0.25,
                 xlab = "X", ylab = "Y", asp = 1)
}



# Merge NFI data
#' @details function imported from the scrip "Merge_NFI_data.R", unlike the other 
#           functions of this first section that were imported from 'READ_FRENCH_NFI.R'

merge_NFI <- function(NFI_tree, NFI_tree_alive_remeasure, NFI_tree_dead_remeasure){
  
  # NFI_tree_alive_remeasure that were alive at first measurement
  # NFI_tree_dead_remeasure that were dead at first measurement
  NFI_tree_dead_remeasure$c135 <- NA
  NFI_tree_remeasure <- rbind(NFI_tree_alive_remeasure, NFI_tree_dead_remeasure)
  library(dplyr)
  NFI_tree_j <- left_join(NFI_tree, NFI_tree_remeasure , 
                          by = c("idp" = "idp", "a" = "a"))
  return(NFI_tree_j)
}
