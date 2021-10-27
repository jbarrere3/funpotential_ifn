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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 1 - Read French NFI ####
#' @description Group of functions used to read French NFI data
#' @authors Georges Kunstler - Patrick Vallet - Björn Reineking (INRAE - LESSEM)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 2 - Format NFI for TreeMort ####
#'
#' @details Group of functions used to format NFI data for the TreeMort project. 
#' @author Julien Barrere
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
           census.n = 1,
           d = round(c13*10/pi, digits = 0), 
           pom = 130, 
           ba = round((c13*10)^2/(4*pi), digits = 0), 
           mode.death = case_when(veget == "5" ~ "1s", 
                                  veget %in% c("C", "A") ~ "1f"), 
           mode.death.other = NA_real_, 
           canopy.position = case_when(lib == 0 ~ 0, 
                                       lib %in% c(1, 2) ~ 1), 
           multistem = case_when(tige %in% c(1,7) ~ 0, 
                                 tige %in% c(5,6) ~ 1), 
           statistical.weight = w) %>%
    rename(plot.id = idp, census.date = year, height = htot, tree.status = dead) %>%
    merge(TreeMort_species, by.x = "espar", by.y = "species.id", all.x = T, all.y = F) %>%
    select(tree.id, plot.id, census.id, census.date, census.n, species, genus, family, d, 
           pom, height, statistical.weight, ba, tree.status, mode.death, mode.death.other,
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
           census.n = 2,
           d = round(c135*1000/pi, digits = 0), 
           ba = round((c135*1000)^2/(4*pi), digits = 0), 
           height = NA_real_, 
           tree.status = dead, 
           mode.death = case_when(veget5 == "M" ~ "1s", 
                                  veget5 %in% c("A", "1", "2") ~ "1f", 
                                  veget5 %in% c("6", "7") ~ "2", 
                                  veget5 %in% c("N", "T") ~ "3"), 
           mode.death.other = NA_real_,
           canopy.position = NA_real_, 
           statistical.weight = NA_real_) %>%
    select(tree.id, plot.id, census.id, census.date, census.n, species, genus, family, d, 
           pom, height, statistical.weight, ba, tree.status, mode.death, mode.death.other,
           canopy.position, multistem)
}




#' Convert L93 to degree
#' @details Function to convert the L93 plot coordinates of the NFI plots into degrees
#' @param l93_data dataframe with two columns: x and y coordinates in L93
#' @return A dataframe with two columns, lon and lat. 

l93_to_degrees <- function(l93_data){
  colnames(l93_data) <- c("lon", "lat")
  l93_data <- l93_data %>%
    st_as_sf(coords = c("lon", "lat")) %>%
    st_set_crs(2154) %>% 
    st_transform(4326) %>%
    as.data.frame %>%
    dplyr::mutate(geometry = as.character(geometry))
  out <- data.frame(longitude = array(NA, dim = dim(l93_data)[1]), 
                    latitude = array(NA, dim = dim(l93_data)[1]))
  for(i in 1:dim(l93_data)[1]){
    out[i, ] <- strsplit(substring(as.character(l93_data[i, 1]), 3, nchar(as.character(l93_data[i, 1])) -1),"\\, ")[[1]]
  }  
  out                  
}


#' Compute NFI Plot Stand Age
#' @details Function to the stand age of NFI plots from the NFI tree table
#' @param NFI_tree Table containing NFI tree data from the first measurement (dead and alive trees)
#' @return a data.table object

Compute_NFI_stand_age <- function(NFI_tree){
  NFI_tree %>%
    group_by(idp) %>%
    summarise(n = sum(!is.na(age)), 
              mean = mean(age, na.rm = TRUE)) %>%
    mutate(stand.age = case_when(n > 0 ~ mean, 
                                 TRUE ~ NA_real_)) %>%
    select(idp, stand.age)
}


#' Format Plots TreeMort
#' @details Function to format the NFI plot table to fit into TreeMort template
#' @param NFI_plot Table containing NFI plot data (code, french and latin name)
#' @return a data.table object

Format_plots_TreeMort <- function(NFI_plot, NFI_ecological_data, NFI_plot_elevation, NFI_stand_age){
  cbind(NFI_plot, l93_to_degrees(NFI_plot[, c("xl93", "yl93")])) %>%
    merge(NFI_ecological_data, by = "idp") %>%
    merge(NFI_plot_elevation %>% mutate(idp = as.integer(as.character(idp))), by = "idp", all.x = T, all.y = F) %>%
    merge(NFI_stand_age, by = "idp", all.x = T, all.y = F) %>%
    mutate(plot.id = idp, 
           cluster = NA_character_,
           country = "France", 
           active = 0, 
           slope = round(atan(pent2/100)/0.0174515207, digits = 2),
           aspect = case_when(expo %in% c(c(0:49), c(350:400)) ~ "N", 
                              expo %in% c(50:149) ~ "E", 
                              expo %in% c(150:249) ~ "S", 
                              expo %in% c(250:349) ~ "W"),
           elevation = zp,
           soil.depth = prof2, 
           soil.depth.accuracy = case_when(obsprof == 1 ~ 1, 
                                           obsprof %in% c(2,3) ~ 2, 
                                           obsprof == 4 ~ 3),
           soil.type = tsol,
           change.protocol = 0,
           d.threshold = 75, 
           plot.contact = "Julien BARRERE", 
           plot.contact.email = "julien.barrere@inrae.fr") %>%
    select(plot.id, latitude, longitude, cluster, country, active, slope, 
           aspect, elevation, soil.depth, soil.depth.accuracy, soil.type, 
           stand.age, change.protocol, d.threshold, plot.contact, plot.contact.email)
}


#' Compute management per census
#' @details Function to associate a "management" value to each census
#' @param NFI_plot Table containing NFI plot data 
#' @param NFI_plot_remeasure Table containing NFI plot data for remeasured plots
#' @return a data.table object

compute_management_census <- function(NFI_plot, NFI_plot_remeasure){
  rbind((NFI_plot %>% select(idp, dc, year)), 
        (NFI_plot_remeasure %>% select(idp, dc5, year) %>% rename(dc = dc5))) %>%
    mutate(management = case_when(dc == 0 ~ 0, 
                                  dc > 0 ~ 1), 
           census.id = paste(idp, year, sep = "_")) %>%
    select(census.id, management)
}




#' Format Census TreeMort
#' @details Function to create a census table to fit into TreeMort template
#' @param TreeMort_tree Table containing NFI tree data formatted for treeMort
#' @param NFI_census_management Table containing the variable management for each census 
#' @return a data.table object

Format_census_TreeMort <- function(TreeMort_tree, NFI_census_management){
  TreeMort_tree %>%
    merge(NFI_census_management, by = "census.id", all.x = T, all.y = F) %>%
    mutate(plot.area = NA_real_, 
           lianas = 0, 
           census.contact = "Julien BARRERE", 
           census.contact.email = "julien.barrere@inrae.fr") %>%
    select(plot.id, census.id, census.date, census.n, plot.area, 
           management, lianas, census.contact, census.contact.email) %>%
    distinct()
}




#' Format Meta-data TreeMort
#' @details Function to create a metadata table containing important information on TreeMort formatting
#' @return a data.table object 

Meta_data_TreeMort <- function(){
  data.frame(File = c("tree_data.csv", "plot_data.csv", "plot_data.csv", "plot_data.csv", "plot_data.csv", "census_data.csv"), 
             Variable = c("statistical.weight", "soil.depth", "soil.depth.accuracy", "stand.age", "soil.type", "management"), 
             Observation = c("Number of trees that this sample represent per ha (.ha-1).", 
                             "Available by class. 0=0-4cm; 1=5-14cm; 2=15-24cm; […]; 8=75-84cm; 9>84cm ", 
                             "Soil depth is 1=accurate, 2=likely understimated, 3=likely overestimated", 
                             "Age at 130cm of the two oldest trees in the stand", 
                             "For correspondance, see the table on the two last pages of https://inventaire-forestier.ign.fr/IMG/pdf/IFN_campagne2008_documentation_donnees-brutes_pointforet.pdf", 
                             "Occurrence of a cut (=1) or not (=0) during the 5 years preceeding the census"))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 3 - TreeMort to FUNDIV ####
#'
#' @details Group of functions used to format and merge French NFI remeasured data to FUNDIV dataset. 
#' @author Julien Barrere
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Format tree from to FUNDIV
#' @details Function to format tree data from TreeMort format to FUNDIV template
#' @param TreeMort_tree Tree table of NFI data formatted for TreeMort
#' @param TreeMort_plot Plot table of NFI data formatted for TreeMort
#' @param FUNDIV_plot Plot tables already formatted to FUNDIV template
#' @return a tree table of French NFI remeasured data formatted for FUNDIV

Format_trees_TreeMort_to_FUNDIV <- function(TreeMort_tree, TreeMort_plot, FUNDIV_plot){
  # 1st step: format FNDIV_plot to merge it with TreeMort_tree table (to get climatic data)
  FUNDIV_plot.in <- FUNDIV_plot %>%
    filter(country == "FG") %>%
    separate(plotcode, into = c("plot.id", "country2"), sep = "_") %>%
    mutate(plot.id = as.integer(plot.id)) %>%
    select(plot.id, map.wc, wai.wc, mat.wc, mat, sgdd, spring_frosts, map, 
           wai, p_pet_yr, sws, p_pet_summer, spei_min, spei_mean)
  
  # 2nd step: Create FUNDIV tree from 1st census TreeMort table
  subset(TreeMort_tree, census.n == 1) %>%
    # Add longitude and latitude
    merge((TreeMort_plot %>% select(plot.id, latitude, longitude)), 
          by = "plot.id", all.x = T, all.y = F) %>%
    rename(dbh1 = d, start_year = census.date, n_ha2 = statistical.weight) %>%
    # Only keep trees that were alive during 1st census
    filter(tree.status == 0) %>%
    select(-tree.status, -mode.death) %>%
    # Add remeasured trees
    merge((subset(TreeMort_tree, census.n == 2) %>% select(tree.id, d, tree.status, mode.death) %>% rename(dbh2 = d)), 
          by = "tree.id", all.x = T, all.y = F) %>%
    mutate(plotcode = as.character(plot.id), 
           sp = paste(genus, species, sep = " "), 
           yearsbetweensurveys = as.integer(5), 
           cluster = NA_integer_, 
           country = "FR", 
           ba_ha2 = NA_real_, 
           management = 0, 
           surveydate2 = NA_character_, 
           BATOTcomp = NA_real_, 
           treestatus_th = case_when((dbh1 >= 100 & tree.status == 0) ~ 2, # Alive and dbh1 >10: "survivor"
                                     (dbh1 < 100 & dbh2 < 100) ~ 99, # Still lower than 10cm : To remove ?
                                     mode.death == "2" ~ 3, # Dead harvested : "Harvested"
                                     mode.death == "1s" ~ 4, # Dead standing : "dead + stem present" (natural mortality)
                                     mode.death %in% c("1f", "3") ~ 5, # Dead fallen or no information about death: "dead + stem absent"
                                     (dbh1 < 100 & dbh2 >= 100) ~ 1)) %>% # Grew above dbh100 : "ingrowth"
    filter(treestatus_th != 99) %>%
    merge(FUNDIV_plot.in, by = "plot.id", all.x = T, all.y = F) %>% # Add climatic variables
    # Correct inapropriate formats
    mutate(treecode2 = tree.id,
           plotcode = paste0(plotcode, "_FR"),
           latitude = as.numeric(latitude), 
           longitude = as.numeric(longitude), 
           treestatus_th = as.integer(treestatus_th),
           start_year = as.integer(start_year)) %>%
    select(plotcode, treestatus_th, dbh1, dbh2, ba_ha2, country, sp, cluster, 
           longitude, latitude, yearsbetweensurveys, management, surveydate2, 
           start_year, map.wc, wai.wc, mat.wc, mat, sgdd, spring_frosts, map, 
           wai, p_pet_yr, sws, p_pet_summer, spei_min, spei_mean, n_ha2,
           BATOTcomp, treecode2)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 4 - R. Seidl & C. Senf disturbance data       ####
#' @description Functions used to read disturbance spatial data
#' @authors Julien BARRERE
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(raster)
library(rgdal)
library(sf)


#' Get disturbance per plot
#' @description Function to get the disturbanc etype and year from C. Senf for each plot
#' @param country character: english name of the country, in lower cases (e.g., "france", "spain")
#' @param buffer numeric: radius of the buffer aroutn plot center to get disturbance value
#' @param FUNDIV_tree Tree table of FUNDIV data
#' @param dir Directory where disturbance data are stored

Get_disturbance_per_plot <- function(country, buffer, FUNDIV_tree, dir){
  print(paste0("Getting disturbance value for ", country))
  
  ## Step 1 - Get raster data
  print("-- Getting disturbance raster")
  disturbance_type_raster.in <- raster(paste0(dir, "/storm_fire_other_classification_", 
                                              country, ".tif"))
  disturbance_year_raster.in <- raster(paste0(dir, "/disturbance_year_1986-2020_", 
                                              country, ".tif"))
  
  ## Step 2 - Create a SPDF with a buffer around each FUNDIV plots
  print("-- Creating SPDF for FUNDIV plots")
  # Extract the country code in FUNDIV dataset
  countrycodes.in <- data.frame(code = c("DE", "ES", "FI", "FR", "SW", "WA"), 
                                country = c("germany", "spain", "finland", 
                                            "france", "sweden", "belgium"))
  this_countrycode.in <- countrycodes.in$code[which(countrycodes.in$country == country)]
  # Create Spatial Polygon from plot coordinates
  FUNDIV_plots_polygon.in <-  FUNDIV_tree %>%
    rename(lon = longitude, lat = latitude) %>%
    filter(country == this_countrycode.in) %>%
    dplyr::select(plotcode, lon, lat) %>%
    distinct 
  coordinates(FUNDIV_plots_polygon.in) = ~lon+lat
  proj4string(FUNDIV_plots_polygon.in) <- CRS("+proj=longlat +datum=WGS84")
  FUNDIV_plots_polygon.in <- spTransform(FUNDIV_plots_polygon.in, 
                                         crs(disturbance_type_raster.in))
  # Create a buffer around each plot
  FUNDIV_plots_polygon.in <- gBuffer(FUNDIV_plots_polygon.in, 
                                     width = buffer, byid = T)
  
  ## Step 3 : Extract raster values within each plot
  print("-- Extracting raster value within each plot")
  out <- data.frame(plotcode = FUNDIV_plots_polygon.in@data$plotcode)
  out$disturbance.type <- exact_extract(disturbance_type_raster.in, 
                                        FUNDIV_plots_polygon.in, 
                                        fun = 'majority')
  out$disturbance.year <- exact_extract(disturbance_year_raster.in, 
                                        FUNDIV_plots_polygon.in, 
                                        fun = 'majority')
  out <- out %>%
    mutate(disturbance.year = case_when(disturbance.type > 0 ~ disturbance.year, 
                                        TRUE ~ NA_real_)) %>% data.frame
  return(out)
}
 