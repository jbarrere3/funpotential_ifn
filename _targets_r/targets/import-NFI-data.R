list(
  tar_target(NFI_plot, read_plot(path = "data/FrenchNFI")), 
  tar_target(NFI_plot_elevation, read.csv("data/elevation_NFIplots.csv", sep = ";")),
  tar_target(NFI_ecological_data, read_ecological_data(path = "data/FrenchNFI")),
  tar_target(NFI_tree_alive, read_tree(path = "data/FrenchNFI")), 
  tar_target(NFI_tree_dead, read_dead_tree(path = "data/FrenchNFI")), 
  tar_target(NFI_plot_remeasure, read_reameasure_plot(path = "data/FrenchNFI_remeasure")), 
  tar_target(NFI_tree_alive_remeasure, read.csv("data/donnees_arbres_retour/foret_c135.csv", sep = ";")), 
  tar_target(NFI_tree_dead_remeasure, read.csv("data/donnees_arbres_retour/foret_morts.csv", sep = ";"))
)

