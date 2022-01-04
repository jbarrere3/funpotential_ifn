list(
  tar_target(FUNDIV_tree_original, 
             correct_weight_FUNDIV(rbind(FUNDIV_tree_original_NoFR, FUNDIV_tree_original_FR))), 
  tar_target(FUNDIV_plots_original, 
             rbind(FUNDIV_plots_original_noFR, FUNDIV_plots_original_FR))
)
