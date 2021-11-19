list(tar_target(disturbance_per_plot_15m, 
           Get_disturbance_per_plot(15, FUNDIV_tree, "data/Disturbance")), 
     tar_target(disturbance_per_plot_100m, 
           Get_disturbance_per_plot(100, FUNDIV_tree, "data/Disturbance")), 
     tar_target(disturbance_per_plot_200m, 
           Get_disturbance_per_plot(200, FUNDIV_tree, "data/Disturbance")))
