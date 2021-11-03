#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_plot.R  
#' @description R script containing all functions relative to data
#               visualisation
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 1 - Data exploration - FUNDIV / Disturbances ####
#' @description Functions to explore the relation between tree
#'              mortality (FUNDIV) and disturbances (Senf 2021)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Plot survival rate disturbance
#'
#' @description Function to plot the survival rate at plot level for each
#'              type of disturbance (fire, storm or other)
#' @param FUNDIV_tree_withdisturbance FUNDIV tree table with disturbance data

Plot_survival_perdisturbance <- function(FUNDIV_tree_withdisturbance){
  FUNDIV_tree_withdisturbance.in <- FUNDIV_tree_withdisturbance %>%
    mutate(Disturbance = case_when(disturbance.type == 0 ~ "No disturbance", 
                                   disturbance.type == 1 ~ "Other disturbance", 
                                   disturbance.type == 2 ~ "Storm disturbance", 
                                   disturbance.type == 3 ~ "Fire disturbance"), 
           tree.alive = case_when(treestatus_th %in% c(1, 2) ~ 1, 
                                  TRUE ~ 0)) %>%
    mutate(Disturbance = case_when((disturbance.year >= start_year & disturbance.year <= start_year+yearsbetweensurveys) ~ Disturbance, 
                                   TRUE ~ "No disturbance")) %>%
    mutate(Disturbance = factor(Disturbance, 
                                levels = c("No disturbance", "Other disturbance", 
                                           "Storm disturbance", "Fire disturbance"))) %>%
    group_by(plotcode, Disturbance) %>%
    summarize(proportion.tree.alive = sum(tree.alive)/n()*100)
  
  label.in <- FUNDIV_tree_withdisturbance.in %>%
    group_by(Disturbance) %>%
    summarize(n = n()) %>%
    mutate(label = paste0("n = ", n), 
           x = 50, y = 0.8)
  
  FUNDIV_tree_withdisturbance.in %>%
    ggplot(aes(x = proportion.tree.alive)) + 
    geom_histogram(aes(y = stat(density) * 5), 
                   binwidth = 5, fill = '#AE2012', color = 'black') + 
    facet_wrap(~ Disturbance) + 
    xlab("Survival percentage") + ylab("Frequency") + 
    theme(panel.background = element_rect(color = 'black', fill = 'white'), 
          panel.grid = element_blank(),
          strip.background = element_blank()) + 
    geom_text(data = label.in, 
              mapping = aes(x = x, y = y, label = label), 
              size = 3.5, fontface = 'italic')
}

#' Map Disturbance FUNDIV plots
#' @description Create 4 maps of FUNDIV plots, one per type fo disturbance
#' @param FUNDIV_tree_withdisturbance FUNDIV tree table with disturbance data

Map_FUNDIVplots_perDisturbance2 <- function(FUNDIV_tree_withdisturbance){
  # Create sf object from FUNDIV_tree_withdisturbance
  FUNDIV_tree_withdisturbance.in <- FUNDIV_tree_withdisturbance %>%
  mutate(Disturbance = case_when(disturbance.type == 0 ~ "No disturbance", 
                                 disturbance.type == 1 ~ "Other disturbance", 
                                 disturbance.type == 2 ~ "Storm disturbance", 
                                 disturbance.type == 3 ~ "Fire disturbance")) %>%
  mutate(Disturbance = factor(Disturbance, 
                              levels = c("No disturbance", "Storm disturbance", 
                                         "Fire disturbance", "Other disturbance"))) %>%
  select(plotcode, longitude, latitude, Disturbance) %>%
  distinct %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
  
  # Map the generated sf object
  ne_countries(scale = "medium", returnclass = "sf") %>%
    ggplot() +
    geom_sf(fill = "#E9ECEF", show.legend = F) +
    geom_sf(data = FUNDIV_tree_withdisturbance.in, 
            size = 0.2, shape = 20, aes(color = Disturbance), 
            show.legend = "point")+
    scale_color_manual(values = c("#ADB5BD", "#0091AD", "#D00000", "#38B000")) +
    coord_sf(expand = FALSE, xlim = c(-25, 40), ylim = c(35, 72)) + 
    annotation_scale(location = "bl", width_hint = 0.13) +
    annotation_north_arrow(location = "tl", which_north = "true", 
                           pad_x = unit(0.35, "in"), pad_y = unit(0.1, "in"),
                           style = north_arrow_fancy_orienteering) +
    facet_wrap(~ Disturbance) +
    theme(panel.background = element_rect(color = 'black', fill = 'white'), 
          panel.grid = element_line(colour = 'lightgray', linetype = "dashed"),
          legend.position = "none", 
          strip.background = element_blank(), 
          strip.text = element_text(size = 12, face = "bold"))
}


#' Compare disturbance NFI Senf
#' @description Function that map the disturbance that occured in France, between the two surveys
#'              with two different sources: Senf and field estimation by NFI agents
#' @param FUNDIV_tree_withdisturbance FUNDIV tree table with disturbance data
#' @param NFI_plot_remeasure Table containing NFI plot data for remeasured plots

Map_compare_NFI_Senf_disturbance <- function(FUNDIV_tree_withdisturbance, NFI_plot_remeasure){
  # Format data
  data.in <- FUNDIV_tree_withdisturbance %>%
    filter(country == "FR") %>%
    mutate(Senf = case_when(disturbance.type == 0 ~ "No disturbance", 
                            disturbance.type == 1 ~ "Other disturbance", 
                            disturbance.type == 2 ~ "Storm disturbance", 
                            disturbance.type == 3 ~ "Fire disturbance")) %>%
    mutate(Senf = case_when((disturbance.year >= start_year & disturbance.year <= start_year+yearsbetweensurveys) ~ Senf, 
                            TRUE ~ "No disturbance")) %>%
    select(plotcode, longitude, latitude, Senf) %>%
    distinct %>%
    merge((NFI_plot_remeasure %>% mutate(plotcode = paste0(idp, "_FR"))), 
          by = "plotcode", all.x = T, all.y = F) %>%
    mutate(FrenchNFI = case_when(nincid5 == 1 ~ "Fire disturbance", 
                                 nincid5 == 3 ~ "Landslide disturbance", 
                                 nincid5 == 4 ~ "Storm disturbance", 
                                 nincid5 %in% c(2, 5) ~ "Other disturbance", 
                                 TRUE ~ "No disturbance")) %>%
    pivot_longer(c(FrenchNFI, Senf), names_to = "Source", values_to = "Disturbance") %>%
    filter(Disturbance != "No disturbance") %>%
    mutate(Disturbance = factor(Disturbance, 
                                levels = c("Storm disturbance", "Fire disturbance", 
                                           "Other disturbance", "Landslide disturbance"))) %>%
    select(plotcode, longitude, latitude, Source, Disturbance) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
  
  # Plot data
  ne_countries(scale = "medium", returnclass = "sf") %>%
    ggplot() +
    geom_sf(fill = "white", show.legend = F) +
    geom_sf(data = data.in, 
            size = 2, shape = 21, aes(fill = Disturbance), 
            color = 'black', show.legend = "point")+
    scale_fill_manual(values = c("#0091AD", "#D00000", "#38B000", "#B5179E")) +
    coord_sf(expand = FALSE, xlim = c(-6, 9), ylim = c(42, 52)) + 
    annotation_scale(location = "br", width_hint = 0.2) +
    annotation_north_arrow(location = "tl", which_north = "true", 
                           pad_x = unit(0.22, "in"), pad_y = unit(0.1, "in"),
                           style = north_arrow_fancy_orienteering) +
    facet_wrap(~ Source) +
    theme(panel.background = element_rect(color = 'black', fill = 'aliceblue'), 
          panel.grid = element_line(colour = 'lightgray', linetype = "dashed"),
          legend.title = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 12, face = "bold"))
}

