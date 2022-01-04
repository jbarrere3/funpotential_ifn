tar_target(FinnishNFI_French_tree, 
           Format_trees_FUNDIV_to_FinnishNFI(subset(FUNDIV_tree, country == "FR"), NFI_tree, 
                                              NFI_plot_elevation, NFI_plot_remeasure, 
                                              NFI_ecological_data))
