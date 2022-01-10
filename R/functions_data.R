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
        (NFI_plot_remeasure %>% select(idp, dc5, year) %>% rename(dc = dc5) %>% mutate(year = year + 5))) %>%
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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 3 - FUNDIV importation and formatting ####
#'
#' @details Group of functions used to format and merge French NFI remeasured data to FUNDIV dataset. 
#' @author Julien Barrere
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Format tree from TreeMort to FUNDIV
#' @details Function to format tree data from TreeMort format to FUNDIV template
#' @param TreeMort_tree Tree table of NFI data formatted for TreeMort
#' @param FUNDIV_species Species table already formatted to FUNDIV template
#' @return a tree table of French NFI remeasured data formatted for FUNDIV
Format_trees_TreeMort_to_FUNDIV <- function(TreeMort_tree, FUNDIV_species){
  out <- subset(TreeMort_tree, census.n == 1) %>%
    rename(dbh1 = d, height1 = height, weight1 = statistical.weight) %>%
    # Only keep trees that were alive during 1st census
    filter(tree.status == 0) %>%
    select(-tree.status, -mode.death) %>%
    # Add remeasured trees
    merge((subset(TreeMort_tree, census.n == 2) %>% select(tree.id, d, tree.status, mode.death) %>% rename(dbh2 = d)), 
          by = "tree.id", all.x = T, all.y = F) %>%
    mutate(plotcode = paste0(plot.id, "_FG"), 
           sp = paste(genus, species, sep = " "), 
           treestatus_th = case_when((dbh1 >= 100 & tree.status == 0) ~ 2, # Alive and dbh1 >10: "survivor"
                                     (dbh1 < 100 & dbh2 < 100) ~ 99, # Still lower than 10cm : To remove ?
                                     mode.death == "2" ~ 3, # Dead harvested : "Harvested"
                                     mode.death == "1s" ~ 4, # Dead standing : "dead + stem present" (natural mortality)
                                     mode.death %in% c("1f", "3") ~ 5, # Dead fallen or no information about death: "dead + stem absent"
                                     (dbh1 < 100 & dbh2 >= 100) ~ 1), # Grew above dbh100 : "ingrowth"
           height2 = NA_real_, 
           ba1 = (pi*((dbh1/1000)/2)^2),
           ba2 = (pi*((dbh2/1000)/2)^2),
           dbh1_mod = NA_real_,
           ba1_mod = NA_real_,
           bachange_ha_yr_mod = NA_real_,
           weight2 = NA_real_, 
           country = "FG") %>% 
    mutate(ba_ha1 = ba1*weight1,
           ba_ha2 = ba2*weight1,
           bachange_ha_yr = weight1*(ba2 - ba1)/5) %>%
    merge((FUNDIV_species %>% mutate(sp = paste(genus, species, sep = " ")) %>% dplyr::select(sp, id)), 
          by = "sp", all.x = T, all.y = F) %>%
    rename(speciesid = id) %>%
    filter(treestatus_th != 99) %>%
    # Correct inapropriate formats
    mutate(treecode = tree.id,
           treestatus_th = as.integer(treestatus_th)) %>%
    select(treecode, plotcode, speciesid, treestatus_th, dbh1, dbh2, height1, height2, ba1, ba_ha1, 
           ba2, ba_ha2, bachange_ha_yr, dbh1_mod, ba1_mod, bachange_ha_yr_mod, weight1, weight2, country)
  return(out)
}


#' Format plots for French data from TreeMort to FUNDIV template
#' @details Function to format tree data from TreeMort format to FUNDIV template
#' @param FUNDIV_tree_FR French NFI tree table formatted for FUNDIV with function Format_trees_TreeMort_to_FUNDIV
#' @param NFI_ecological_data Table containing NFI plot ecological data
#' @param TreeMort_plot Plot table of NFI data formatted for TreeMort
#' @return a plot table of French NFI remeasured data formatted for FUNDIV
Format_plots_TreeMort_to_FUNDIV <- function(NFI_ecological_data, TreeMort_plot, FUNDIV_tree_FR){
  # Get the exact date of the first census
  out.surveydate <- NFI_ecological_data %>%
    mutate(plotcode = paste0(idp, "_FG"), 
           surveydate1 = as.character(as.Date(dateeco, format = "%d/%m/%Y"))) %>%
    dplyr::select(plotcode, surveydate1)
  
  # Get the coordinates
  out.coordinates <- TreeMort_plot %>% 
    mutate(plotcode = paste0(plot.id, "_FG")) %>%
    dplyr::select(plotcode, latitude, longitude)
  
  # Format final dataset
  out <- FUNDIV_tree_FR %>%
    group_by(plotcode, country) %>%
    summarise(ba_ha1 = sum(ba_ha1, na.rm = T), 
              ba_ha2 = sum(ba_ha2, na.rm = T)) %>%
    merge(out.surveydate, by = "plotcode", all.x = T, all.y = F) %>%
    merge(out.coordinates, by = "plotcode", all.x = T, all.y = F) %>%
    mutate(yearsbetweensurveys = 5, 
           surveydate2 = as.character(as.numeric(substr(surveydate1, 1, 4)) + 5), 
           biome = NA_real_, 
           management = NA_real_, 
           cluster = NA_real_) %>%
    dplyr::select(plotcode, cluster, country, longitude, latitude, yearsbetweensurveys, surveydate1, 
                  surveydate2, biome, ba_ha1, ba_ha2, management)
  return(out)
}

#' Correct the weight in FUNDIV raw data
#' @param FUNDIV_tree FUNDIV tree table
correct_weight_FUNDIV <- function(FUNDIV_tree){
  FUNDIV_tree %>%
    mutate(weight1 = case_when(country == "DE" ~ weight1, 
                               country %in% c("ES", "FG", "FI", "SW", "WA") ~ 10000/(pi*weight1^2)), 
           weight2 = case_when(country == "DE" ~ weight2, 
                               country %in% c("ES", "FG", "FI", "SW", "WA") ~ 10000/(pi*weight2^2))) %>% 
    mutate_if(is.numeric, list(~na_if(., Inf)))
}


#' Format FUNDIV data for the IPM
#' @author G. Kunstler (slightly modified by J. Barrere)
#' @param FUNDIV_tree FUNDIV tree table
#' @param FUNDIV_plots FUNDIV plot table
#' @param FUNDIV_species FUNDIV species table
#' @param FUNDIV_climate climate index for each FUNDIV plots
#' @param FUNDIV_management Occurrence of harvesting in the different FUNDIV plots
#' @param remove_harv Boolean to specify whether plots where harvesting occurred should be removed. 
read_FUNDIV_tree_data <- function(FUNDIV_tree, FUNDIV_plots, FUNDIV_species, 
                                  FUNDIV_climate, FUNDIV_management, remove_harv = TRUE){
  # Add species in tree table, and correct anomalies
  FUNDIV_tree.in <- FUNDIV_tree
  FUNDIV_tree.in$speciesid[FUNDIV_tree.in$speciesid %in% c(46,47)] <- 48
  FUNDIV_species.in <- FUNDIV_species
  FUNDIV_species.in$species[FUNDIV_species.in$id ==277] <- "pubescens"
  FUNDIV_species.in <- FUNDIV_species.in %>% mutate(sp = paste(genus, species)) %>%
    dplyr::select(c(id, sp))
  FUNDIV_species.in[FUNDIV_species.in$id == 277, "sp"] <- "Quercus pubescens"
  FUNDIV_species.in[FUNDIV_species.in$id == 48, "sp"] <- "Betula"
  data <- left_join(FUNDIV_tree.in, FUNDIV_species.in, by = c("speciesid" = "id"))
  
  # Add data at the plot level
  FUNDIV_plots.in <- FUNDIV_plots %>% 
    mutate(start_year = as.numeric(substr(surveydate1, 1, 4)), 
           end_year = as.numeric(substr(surveydate2, 1, 4)), 
           longitude = as.numeric(longitude), 
           latitude = as.numeric(latitude)) %>%
    dplyr::select(-country, -ba_ha1, -ba_ha2)
  FUNDIV_climate.in <- FUNDIV_climate %>% 
    filter(country != "FR") %>% # French climatic data are obsolete
    dplyr::select(-longitude, -latitude, -country,-surveydate1, -surveydate2, -start_year, -end_year)
  
  # merging data
  data <- data %>%
    merge(FUNDIV_plots.in, by = "plotcode", all.x = T, all.y = F) %>%
    merge(FUNDIV_climate.in, by = "plotcode", all.x = T, all.y = F) %>% 
    group_by(plotcode) %>%
    mutate(BATOT_ha1=sum(ba_ha1))
  
  # remove plot with harvesting if required
  if (remove_harv){
    data <- group_by(data, plotcode) %>%
      mutate(N_harv = sum(treestatus_th == 3)) %>%
      filter(N_harv <1)
    # remove French plots where management has been recorded, we don't have details
    # on individually harvested trees, the management column contains a 1 for those French
    # plots in which management has been recorded
    data$management[is.na(data$management)] <- 0
    data <- filter(data, management ==0)
    # remove plots with harvesting based on plot code
    FUNDIV_management.in <- FUNDIV_management
    plots_with_harv <- FUNDIV_management.in$plotcode[FUNDIV_management.in$management2 >0 & !is.na(FUNDIV_management.in$management2)]
    data <- filter(data, ! plotcode %in% plots_with_harv)
  }
  
  
  # Compute additional variables for competition
  data <- data %>%
    filter(latitude >30) %>% # remove plots in the Canary Islands
    # calculate the number per hectare from the weight
    mutate(n_ha1 = weight1, 
           n_ha2 = case_when(country == "FG" ~ ba_ha2/ba2, 
                             TRUE ~ weight2)) %>% 
    mutate(BATOTcomp = BATOT_ha1 - ba_ha1, # Compute basal area of competitors
           country=replace(country, country=='FG', 'FR')) %>%
    filter(!is.na(BATOT_ha1) & !is.na(surveydate2)& !is.na(surveydate1))
  data$treecode2 <- paste0("code_", seq_len(length.out = nrow(data)))
  
  # remove unused variables
  if (remove_harv){
    data <- select(data, -c(treecode, speciesid, height1, height2, ba1,
                            ba_ha1, ba2, bachange_ha_yr, dbh1_mod,
                            ba1_mod, bachange_ha_yr_mod, weight1, weight2,
                            biome, surveydate1, end_year, tile,
                            bio1, BATOT_ha1, N_harv, n_ha1))
  }else{
    data <- select(data, -c(treecode, speciesid, height1, height2, ba1,
                            ba_ha1, ba2, bachange_ha_yr, dbh1_mod,
                            ba1_mod, bachange_ha_yr_mod, weight1, weight2,
                            biome, surveydate1, end_year, tile,
                            bio1, BATOT_ha1, n_ha1))
  }
  return(data)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 4 - Disturbance data                          ####
#' @description Functions used to read disturbance spatial data
#' @authors Julien BARRERE
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(raster)
library(rgdal)
library(sf)


#' Get disturbance area and year per plot for a specific country
#' @description Function to get the type, area and year of each disturbance intercepting a FUNDIV plot
#' @param country character: english name of the country, in lower cases (e.g., "france", "spain")
#' @param buffer numeric: radius of the buffer aroutn plot center to get disturbance value
#' @param FUNDIV_tree Tree table of FUNDIV data
#' @param dir Directory where disturbance data are stored

Get_disturbance_per_plot_per_country <- function(country, buffer, FUNDIV_tree, dir){
  ## Step 1 - Get raster data
  print("-- Getting disturbance raster")
  disturbance_type_raster.in <- raster(paste0(dir, "/storm_fire_other_classification_", 
                                              country, ".tif"))
  disturbance_year_raster.in <- raster(paste0(dir, "/disturbance_year_1986-2020_", 
                                              country, ".tif"))
  disturbance_severity_raster.in <- raster(paste0(dir, "/disturbance_severity_1986-2016_", 
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
  print("----------Raster 1: disturbance type")
  disturbance_type_raster.in_extract <- exact_extract(
    disturbance_type_raster.in, FUNDIV_plots_polygon.in, coverage_area = T)
  print("----------Raster 2: disturbance year")
  disturbance_year_raster.in_extract <- exact_extract(
    disturbance_year_raster.in, FUNDIV_plots_polygon.in, coverage_area = T)
  print("----------Raster 3: disturbance severity")
  disturbance_severity_raster.in_extract <- exact_extract(
    disturbance_severity_raster.in, FUNDIV_plots_polygon.in, coverage_area = T)
  
  ## Step 4 : Format the dataset
  print("-- Formating the data")
  # Inclusion of type
  out <- data.frame(t(sapply(disturbance_type_raster.in_extract, c))) %>%
    mutate(V = substr(value, start = 3, stop = nchar(value) - 1), 
           A = substr(coverage_area, start = 3, 
                      stop = nchar(coverage_area) - 1)) %>%
    separate(V, into = paste0("type.", c(1:max(sapply(disturbance_type_raster.in_extract, nrow)))),
             sep = ", ", fill = "right") %>%
    separate(A, into = paste0("area.", c(1:max(sapply(disturbance_type_raster.in_extract, nrow)))), 
             sep = ", ", fill = "right") %>%
    dplyr::select(-value, -coverage_area) %>%
    # Inclusion of year
    cbind(data.frame(t(sapply(disturbance_year_raster.in_extract, c)))) %>%
    mutate(value = substr(value, start = 3, stop = nchar(value) - 1), 
           plotcode = FUNDIV_plots_polygon.in@data$plotcode) %>%
    separate(value, into = paste0("year.", c(1:max(sapply(disturbance_type_raster.in_extract, nrow)))),
             sep = ", ", fill = "right") %>%
    dplyr::select(-coverage_area) %>%
    # Inclusion of severity
    cbind(data.frame(t(sapply(disturbance_severity_raster.in_extract, c)))) %>%
    mutate(value = substr(value, start = 3, stop = nchar(value) - 1)) %>%
    separate(value, into = paste0("seve.", c(1:max(sapply(disturbance_severity_raster.in_extract, nrow)))),
             sep = ", ", fill = "right") %>%
    dplyr::select(-coverage_area) %>%
    # Finalize formatting
    gather(key, value, -plotcode) %>%
    mutate(info = sub("\\..+", "", key), 
           cell = as.numeric(sub("....\\.", "", key))) %>%
    dplyr::select(-key) %>%
    spread(info, value) %>%
    filter(type %in% as.character(c(1:3))) %>%
    mutate(type = as.numeric(type), 
           year = as.numeric(year), 
           severity = as.numeric(seve)/100, 
           coverage = as.numeric(area)/(pi*(buffer^2))) %>%
    dplyr::select(plotcode, type, year, severity, coverage)
                 
  return(out)
}



#' Get disturbance area and year per plot for all countries
#' @description Useful if the same buffer is applied for each country 
#' @param buffer numeric: radius of the buffer aroutn plot center to get disturbance value
#' @param FUNDIV_tree Tree table of FUNDIV data
#' @param dir Directory where disturbance data are stored

Get_disturbance_per_plot <- function(buffer, FUNDIV_tree, dir){
  countries.in = c("germany", "spain", "finland", 
                   "france", "sweden", "belgium")
  for(i in 1:length(countries.in)){
    print(paste0("Getting disturbance data for ", countries.in[i]))
    if(i == 1) out <- Get_disturbance_per_plot_per_country(countries.in[i], buffer, FUNDIV_tree, dir)
    else out <- rbind.data.frame(out, Get_disturbance_per_plot_per_country(countries.in[i], buffer, FUNDIV_tree, dir))
  }
  return(out)
}




#' Add disturbance to FUNDIV
#' @description Function to sum all disturbances that occured in FUNDIV plots between the two censuses
#'              Add 3 additional columns to the entry dataset, cumulative % of the area affected by 
#'              each of the 3 types of disturbance
#' @param FUNDIV_tree Tree table of FUNDIV data
#' @param disturbance_per_plot dataframe containging the type, area and year of each 
#'                                       disturbance intercepting a FUNDIV plot
#'                                       (computed with the previous function)

add_disturbance_to_FUNDIV <- function(FUNDIV_tree, disturbance_per_plot){
  
  # Extract area of each type of disturbance that occured between census 1 and 2
  FUNDIV.disturbance.in <- FUNDIV_tree %>%
    dplyr::select(plotcode, start_year, yearsbetweensurveys) %>%
    mutate(end_year = as.numeric(start_year + yearsbetweensurveys), 
           start_year = as.numeric(start_year)) %>%
    distinct() %>%
    merge(disturbance_per_plot, by = "plotcode", all.y = T) %>%
    #filter(year >= start_year & year <= end_year) %>%
    mutate(disturbance = case_when(type == 1 ~ "disturbance.other", 
                                   type == 2 ~ "disturbance.storm", 
                                   type == 3 ~ "disturbance.fire")) %>%
    dplyr::select(plotcode, disturbance, coverage, severity) %>%
    mutate(coverage2 = coverage*severity) %>%
    group_by(plotcode, disturbance) %>%
    summarize(Area = sum(coverage2, na.rm = T)) %>%
    spread(key = disturbance, value = Area, fill = NA) 
  
  # Add to the FUNDIV dataset and organize columns to ensure uniformity
  out <- FUNDIV_tree %>%
    merge(FUNDIV.disturbance.in, by = "plotcode", all.x = T, all.y = F) %>%
    mutate(disturbance.other = if("disturbance.other" %in% colnames(.)) disturbance.other else NA_real_, 
           disturbance.storm = if("disturbance.storm" %in% colnames(.)) disturbance.storm else NA_real_,
           disturbance.fire = if("disturbance.fire" %in% colnames(.)) disturbance.fire else NA_real_) %>%
    dplyr::select(c(colnames(FUNDIV_tree), "disturbance.storm", "disturbance.fire", "disturbance.other"))
  return(out)
}



#' Get Annual prevalence of each disturbance per country
#' @description Function that compute the percentage of all plot buffer area 
#'              affected by each disturbance per country and per year
#' @param disturbance_per_plot dataframe containging the type, area and year of each 
#'                                       disturbance intercepting a FUNDIV plot buffer
#' @param FUNDIV_tree Tree table of FUNDIV data
#' @return A list with two df: @data.plot contains the prevalence (in %) per country, 
#'                                        year and disturbance
#'                             @country.date contains the year of the oldest 1st census, 
#'                                           and most recent 2nd census per country
get_annualPrevalence <- function(disturbance_per_plot, FUNDIV_tree){
  disturbance.in <- disturbance_per_plot %>%
    merge((FUNDIV_tree %>% dplyr::select(plotcode, country) %>% distinct()), 
          by = "plotcode", all.x = T, all.y = F) %>%
    mutate(disturbance.type = case_when(type == 3 ~ "Fire", 
                                        type == 2 ~ "Storm", 
                                        type == 1 ~ "Other")) %>%
    group_by(country, year, disturbance.type) %>%
    summarise(Area = sum(coverage)) %>%
    filter(year %in% c(1986:2020)) %>%
    mutate(id = paste(country, year, disturbance.type, sep = "_")) %>%
    ungroup() %>%
    dplyr::select(id, Area)
  
  countryLevelInfo.in <- FUNDIV_tree %>%
    mutate(end_year = start_year + yearsbetweensurveys) %>%
    dplyr::select(plotcode, country, start_year, end_year) %>%
    distinct() %>%
    group_by(country) %>%
    summarise(n.plot = n(), 
              min.year = min(start_year, na.rm = T), 
              max.year = max(end_year, na.rm = T))
  
  out.data <- expand.grid(country = c("DE", "ES", "FI", "FR", "SW", "WA"), 
                          year = c(1986:2020), 
                          disturbance.type = c("Fire", "Other", "Storm")) %>%
    mutate(id = paste(country, year, disturbance.type, sep = "_")) %>%
    merge(disturbance.in, by = "id", all.x = T) %>%
    replace_na(list(Area = 0)) %>%
    merge((countryLevelInfo.in %>% dplyr::select(country, n.plot)), 
          by = "country", all.x = T) %>%
    mutate(Prevalence = round(Area/n.plot*100, digits = 4), 
           Country = case_when(country == "DE" ~ "Germany", 
                               country == "ES" ~ "Spain", 
                               country == "FI" ~ "Finland", 
                               country == "FR" ~ "France", 
                               country == "SW" ~ "Sweden", 
                               country == "WA" ~ "Belgium")) %>%
    dplyr::select(Country, year, disturbance.type, Prevalence)
  
  countryLevelInfo.in <- countryLevelInfo.in %>%
    mutate(Country = case_when(country == "DE" ~ "Germany", 
                               country == "ES" ~ "Spain", 
                               country == "FI" ~ "Finland", 
                               country == "FR" ~ "France", 
                               country == "SW" ~ "Sweden", 
                               country == "WA" ~ "Belgium")) %>%
    gather(key = "type", value = "Year", "min.year", "max.year") %>%
    dplyr::select(Country, Year)
  
  out <- list()
  out$data.plot <- out.data
  out$country.date <- countryLevelInfo.in
  return(out)
}


# Import and format AGRESTE salvage logging data
#' @description Import and format AGRESTE data on salvage logging at national level
#' @param path.in path to the csv file containing agreste data
#' @return A dataframe formated in english

import_agreste <- function(path.in){
  as.data.frame(t(read.csv(path.in, header = F, sep = ";"))) %>%
    `colnames<-`(.[1, ]) %>%
    .[-1, ] %>%
    mutate_if(is.character, as.numeric) %>%
    rename('year' = 'Recolte', 'Log' = 'grumes', 'Total' = 'total', 
           'Industry' = 'industrie', 'Energy' = 'energie') %>%
    gather('Total', 'Industry', 'Log', 'Energy', 
           key = 'Use', value = 'Volume')
}






#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 5 - Format to Finnish NFI ####
#'
#' @details Group of functions used to 
#'          format French NFI remeasured 
#'          data to Finnish NFI format. 
#' @author Julien Barrere
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Format tree level table from FUNDIV to Finnish NFI
#' @param TreeMort_tree Tree table of NFI data formatted for TreeMort
#' @param NFI_tree NFI data at tree level for both dead and alive trees
#' @param NFI_plot_elevation Table linking plot id with elevation asl (m)
#' @param NFI_plot_remeasure Table containing NFI plot data for remeasured plots
#' @param NFI_ecological_data Table containing NFI plot ecological data (to get soil information)
#' @return a tree table of French NFI remeasured data formatted for Finnish NFI

Format_trees_FUNDIV_to_FinnishNFI <- function(FUNDIV_FrenchNFI_tree, NFI_tree, 
                                              NFI_plot_elevation, NFI_plot_remeasure, 
                                              NFI_ecological_data){
  # extract information contained only in the original NFI dataset
  original.NFI.in <- NFI_tree %>%
    mutate(tree_id = paste(idp, a, sep = "_")) %>%
    filter(tree_id %in% FUNDIV_FrenchNFI_tree$treecode2) %>%
    dplyr::select(tree_id, lib, htot, age)
  
  # create a finnish dataset for the first French NFI census
  out.census1 <- FUNDIV_FrenchNFI_tree %>%
    rename(tree_species = sp, 
           sampling_year = start_year, 
           dbh = dbh1,
           tree_id = treecode2, 
           Years_between_samplings = yearsbetweensurveys) %>%
    mutate(area_identifier1 = sub("_..", "", plotcode)) %>%
    merge(original.NFI.in, by = "tree_id", all.x = T) %>%
    merge((NFI_plot_elevation %>% mutate(idp = as.character(idp))), 
          by.x = "area_identifier1", by.y = "idp", all.x = T, all.y = F) %>%
    merge((NFI_plot_remeasure %>% mutate(idp = as.character(idp))), 
          by.x = "area_identifier1", by.y = "idp", all.x = T, all.y = F) %>%
    merge((NFI_ecological_data %>% mutate(idp = as.character(idp)) %>% dplyr::select(idp, prof2)), 
          by.x = "area_identifier1", by.y = "idp", all.x = T, all.y = F) %>%
    mutate(alive_now = 1, 
           alive_in_the_next_inventory = case_when(treestatus_th %in% c(1, 2) ~ 1, 
                                                   TRUE ~ 0), 
           NFI_identifier = "census1", 
           harvesting_done = case_when((dc5 %in% c(1:4) | prelev5 == 1) ~ 1, 
                                       TRUE ~ 0), 
           occurence_of_disturbance = case_when(nincid5 %in% c(1:5) ~ 1, 
                                                nincid5 == 0 ~ 0), 
           disturbance_agent = case_when(nincid5 %in% c(1:5) ~ nincid5, 
                                         TRUE ~ NA_integer_)) %>%
    rename(canopy_layer = lib, 
           height = htot, 
           lat = latitude, 
           lon = longitude, 
           elevation = zp, 
           severity_of_disturbance = incid5,
           org_layer_thickness = prof2, 
           statistical_weight = n_ha2, 
           status = treestatus_th) %>%
    mutate(area_identifier2 = NA_integer_, area_identifier3 = NA_integer_, 
           study_plot = NA_integer_, study_plot_id = NA_integer_, 
           study_plot_part = NA_integer_, height_to_living_canopy = NA_real_, 
           direction = NA_real_, distance = NA_real_, height_to_dry_branches = NA_real_, 
           study_plot_type = NA_character_, sampling_day = NA_real_, sampling_month = NA_real_, 
           sampling_group= NA_character_, site_type = NA_integer_, org_layer_quality = NA_integer_, 
           soil_type = NA_integer_, grain_size = NA_real_, time_since_harvesting = NA_real_, 
           tree_cohort = NA_integer_, dominat_tree_species = NA_integer_, 
           BA_in_plot_surroundings = NA_real_, appearance_of_disturbance = NA_integer_, 
           time_since_latest_disturbance = NA_real_, x_in_plot = NA_real_, y_in_plot = NA_real_, 
           bal0 = NA_real_, bal4 = NA_real_, bal0_dbh = NA_real_, bal4_dbh = NA_real_, 
           bal0_dbh_sq = NA_real_, bal4_dbh_sq = NA_real_, bal0_dbh_sq_rank = NA_real_) %>%
    dplyr::select(area_identifier1, area_identifier2, area_identifier3, study_plot, study_plot_id, 
                  study_plot_part, tree_id, alive_now, alive_in_the_next_inventory, tree_species, 
                  dbh, canopy_layer, height_to_living_canopy, height, age, direction, distance, 
                  height_to_dry_branches, NFI_identifier, study_plot_type, lat, lon, elevation, 
                  sampling_day, sampling_month, sampling_year, Years_between_samplings, sampling_group, 
                  site_type, org_layer_quality, org_layer_thickness, soil_type, grain_size, 
                  harvesting_done, time_since_harvesting, tree_cohort, dominat_tree_species, 
                  BA_in_plot_surroundings, occurence_of_disturbance, severity_of_disturbance, 
                  appearance_of_disturbance, time_since_latest_disturbance, disturbance_agent, 
                  x_in_plot, y_in_plot, bal0, bal4, bal0_dbh, bal4_dbh, bal0_dbh_sq, 
                  bal4_dbh_sq, bal0_dbh_sq_rank, statistical_weight, status)
  
  out <- FUNDIV_FrenchNFI_tree %>%
    dplyr::select(treecode2, dbh2) %>%
    rename(tree_id = treecode2) %>%
    merge(out.census1, by = "tree_id") %>%
    mutate(dbh = dbh2, 
           alive_now = alive_in_the_next_inventory, 
           alive_in_the_next_inventory = NA_real_, 
           sampling_year = sampling_year + 5, 
           height = NA_real_, 
           age = age+5, 
           canopy_layer = NA_real_, 
           NFI_identifier = "census2", 
           statistical_weight = NA_real_, 
           status = NA_integer_) %>%
    dplyr::select(colnames(out.census1)) %>%
    rbind.data.frame(out.census1) %>%
    mutate(occurence_of_disturbance = case_when(NFI_identifier == "census1" ~ NA_real_, 
                                                NFI_identifier == "census2" ~ occurence_of_disturbance), 
           severity_of_disturbance = case_when(NFI_identifier == "census1" ~ NA_integer_, 
                                               NFI_identifier == "census2" ~ severity_of_disturbance)) %>%
    arrange(area_identifier1, tree_id, NFI_identifier)
  
}






#' Format Meta-data Finnish NFI
#' @details Function to create a metadata table containing important information on Finnish NFI formating
#' @return a data.table object 

Meta_data_FinnishNFI <- function(){
  data.frame(Variable = c("canopy_layer", "NFI_identifier", "lat", "lon", "harvesting_done", "severity_of_disturbance",
                          "disturbance_agent", "org_layer_thickness", "statistical.weight", "status"), 
             Observation = c("Indicates the proportion of the crown that has direct access to light. 0=0%; 1=1%-66%; 2=67%-100%", 
                             "Indicates whether the tree was sampled during the first (census1) or second (census2) inventory round", 
                             "latitude in WGS projection", 
                             "longitude in WGS projection",
                             "Did a harvest occured between the two inventory rounds ?",
                             "% of damage caused by the disturbance, available by class. 0=0%, 1=1-25%, 2=25-50%, 3=50-75%, 4=75-100%", 
                             "1=Fire, 2=Natural mortality, 3=Landslide, 4=Storm, 5=Other incident",
                             "Available by class. 0=0-4cm; 1=5-14cm; 2=15-24cm; […]; 8=75-84cm; 9>84cm ", 
                             "NEW VARIABLE - Number of trees that this sample represents per ha (.ha-1)", 
                             "NEW VARIABLE - Detailed status of the tree. 1=ingrowth, 2=survivor, 3=dead(harvested), 4=dead(stem present), 5=dead(stem absent)"))
}

