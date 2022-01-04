list(
  tar_target(FUNDIV_tree_original_NoFR, read.csv(file.path("data/FUNDIV", "FunDiv_trees_Nadja.csv"),
                                                 stringsAsFactors=FALSE) %>% dplyr::filter(country != "FG")),
  tar_target(FUNDIV_species, read.csv(file.path("data/FUNDIV", "FunDiv_species_Nadja.csv"),
                                     stringsAsFactors=FALSE, fileEncoding = "cp1252")), 
  tar_target(FUNDIV_plots_original_noFR, read.csv(file.path("data/FUNDIV", "FunDiv_plots_Nadja.csv"),
                                                  stringsAsFactors=FALSE) %>% dplyr::filter(country != "FG")),
  tar_target(FUNDIV_climate, read.csv(file.path("data/FUNDIV", "FunDiv_plots_climate_sapropos.csv"),
                                      stringsAsFactors=FALSE)),
  tar_target(FUNDIV_management, read.csv(file.path("data/FUNDIV", "FunDivEUROPE_plot_management.csv"),
                                      stringsAsFactors=FALSE))
)
