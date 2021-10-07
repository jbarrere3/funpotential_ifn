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
#' @authors Georges Kunstler - Patrick Vallet - Björn Reineking (INRAE - LESSEM)

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
read_plot <- function(path = "data", years = 2008:2019) {
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
read_tree <- function(path = "data", years = 2008:2019) {
  read_ifn_data("arbres_foret", path, years)
}

#' Read cover data for each year
#'
#' @description Function to call the base calling function for the IFN cover
#' data for each year.
#' @export
#' @inheritParams read_plot
#' @return a data.table object
read_couv <- function(path = "data", years = 2008:2019) {
  read_ifn_data("couverts_foret", path, years)
}

#' Read dead tree data for each year
#'
#' @description Function to call the base calling function for the IFN dead
#' dead tree data for each year.
#' @export
#' @inheritParams read_plot
#' @return a data.table object
read_dead_tree <- function(path = "data", years = 2008:2019) {
  read_ifn_data("arbres_morts_foret", path, years)
}

#' Read documentation data for each year
#'
#' @description Function to call the base calling function for the IFN
#' documentation data for each year.
#' @export
#' @inheritParams read_plot
#' @return a data.table object
read_doc <- function(path = "data", years = 2008:2019) {
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
read_ecological_data <- function(path = "data", years = 2008:2019) {
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



# Merge NFI data for trees alive at first measurement
#' @details function imported from the scrip "Merge_NFI_data.R", unlike the other 
#           functions of this first section that were imported from 'READ_FRENCH_NFI.R'

merge_NFI_alive <- function(NFI_tree, NFI_tree_alive_remeasure){
  
  # NFI_tree_alive_remeasure that were alive at first measurement
  # NFI_tree_dead_remeasure that were dead at first measurement
  library(dplyr)
  NFI_tree_j <- left_join(NFI_tree, NFI_tree_alive_remeasure, 
                          by = c("idp" = "idp", "a" = "a"))
  return(NFI_tree_j)
}

#' Merge NFI remeasured dead and alive trees
#' @details Function that merges the two tables of remeasured trees (dead and alive) 
#' @param NFI_tree_alive_remeasure
#' @param NFI_tree_dead_remeasure
#' @return A data.table object

merge_NFI_remeasured <- function(NFI_tree_alive_remeasure, NFI_tree_dead_remeasure){
  NFI_tree_dead_remeasure$c135 <- NA
  NFI_tree_dead_remeasure$dead <- 1
  NFI_tree_alive_remeasure %>%
    mutate(dead = case_when(veget5 == "0" ~ 0, 
                            TRUE ~ 1)) %>%
    rbind(NFI_tree_dead_remeasure)
}



#### Section 2 - Format NFI for TreeMort ####
#' @details Group of functions used to format NFI data for the TreeMort project. 
#' @author Julien Barrere


#' Format Species TreeMort
#' @details Function to format the NFI species table to fit into TreeMort template
#' @param NFI_species Table containing NFI species data (code, french and latin name)
#' @param NFI_genus_species_correspondence Table linking genus of NFI_species to their family
#' @return a data.table object

Format_species_TreeMort <- function(NFI_species, NFI_genus_species_correspondence){
  NFI_species %>% 
    separate(Latin_name, into = c("genus", "species"), sep = "_") %>%
    mutate(species = case_when(species == "sp" ~ NA_character_, 
                               TRUE ~ species)) %>%
    merge(NFI_genus_species_correspondence, by = "genus", all.x = T, all.y = F) %>%
    rename(species.id = code) %>%
    select(species.id, species, genus, family) %>%
    distinct()
} 


#' Format Trees 1st Census TreeMort
#' @details Function to format the NFI tree data (1st measurement) of dead and alive trees
#'          to fit into TreeMort template
#' @param NFI_tree NFI data at tree level for both dead and alive trees
#' @param TreeMort_species NFI species data formatted for TreeMort
#' @return A data.table object

Format_trees_census1_TreeMort <- function(NFI_tree, TreeMort_species){
  NFI_tree %>%
    mutate(tree.id = paste(idp, a, sep = "_"), 
           census.id = paste(idp, year, sep = "_"), 
           d = c13*10/pi, 
           pom = 130, 
           ba = (c13*10)^2/(4*pi), 
           mode.death = case_when(veget == "5" ~ "1s", 
                                  veget %in% c("C", "A") ~ "1f"), 
           mode.death.other = NA_character_, 
           canopy.position = case_when(lib == 0 ~ 0, 
                                       lib %in% c(1, 2) ~ 1), 
           multistem = case_when(tige %in% c(1,7) ~ 0, 
                                 tige %in% c(5,6) ~ 1)) %>%
    rename(plot.id = idp, census.date = year, height = htot, tree.status = dead) %>%
    merge(TreeMort_species, by.x = "espar", by.y = "species.id", all.x = T, all.y = F) %>%
    select(tree.id, plot.id, census.id, census.date, species, genus, family, d, 
           pom, height, ba, tree.status, mode.death, mode.death.other,
           canopy.position, multistem)
}



#' Format Remeasured Trees (2nd census) TreeMort
#' @details Function to format the NFI tree data (1st measurement) of dead and alive trees
#'          to fit into TreeMort template
#' @param NFI_tree_remeasured NFI data at tree level for both dead and alive remeasured trees
#' @param TreeMort_tree_census1 NFI tree census1 data formatted for TreeMort
#' @return A data.table object

Format_trees_census2_TreeMort <- function(NFI_tree_remeasured, TreeMort_tree_census1){
  NFI_tree_remeasured %>%
    mutate(tree.id = paste(idp, a, sep = "_")) %>%
    merge(TreeMort_tree_census1, by = "tree.id", all.x = T, all.y = F) %>%
    mutate(census.id = paste(plot.id, (census.date +5), sep = "_"), 
           census.date = census.date + 5, 
           d = round(c135*1000/pi, digits = 0), 
           ba = round((c135*1000)^2/(4*pi), digits = 0), 
           height = NA_real_, 
           tree.status = dead, 
           mode.death = case_when(veget5 == "M" ~ "1s", 
                                  veget5 %in% c("A", "1", "2") ~ "1f", 
                                  veget5 %in% c("6", "7") ~ "2", 
                                  veget5 %in% c("N", "T") ~ "3"), 
           canopy.position = NA_real_) %>%
    select(tree.id, plot.id, census.id, census.date, species, genus, family, d, 
           pom, height, ba, tree.status, mode.death, mode.death.other,
           canopy.position, multistem)
}