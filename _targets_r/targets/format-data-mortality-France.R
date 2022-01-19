list(
  tar_target(data_jags_Aalba, 
           format_FrenchNFI_mortality(FUNDIV_tree, FUNDIV_climate, 
                                      disturbance_per_plot_200m, 
                                      NFI_plot_remeasure, "Abies alba")), 
  tar_target(data_jags_Aalba_2, 
           format_FrenchNFI_mortality_2(FUNDIV_tree, FUNDIV_climate, 
                                      disturbance_per_plot_200m, 
                                      NFI_plot_remeasure, "Abies alba"))
)
