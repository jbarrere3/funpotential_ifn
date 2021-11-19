list(tar_target(FUNDIV_tree_disturbance15m, 
           add_disturbance_to_FUNDIV(FUNDIV_tree, disturbance_per_plot_15m)), 
     tar_target(FUNDIV_tree_disturbance100m, 
           add_disturbance_to_FUNDIV(FUNDIV_tree, disturbance_per_plot_100m)), 
     tar_target(FUNDIV_tree_disturbance200m, 
           add_disturbance_to_FUNDIV(FUNDIV_tree, disturbance_per_plot_200m)))
