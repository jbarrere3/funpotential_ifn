list(
  tar_target(plot_disturbedArea_perMortalityRate_15m, 
           plot_areaDisturbance_perMortalityRate(FUNDIV_tree_disturbance15m, 
                                                 "Buffer = 15m radius", 
                                                 death.in = "natural mortality")), 
  tar_target(plot_disturbedArea_perMortalityRate_100m, 
           plot_areaDisturbance_perMortalityRate(FUNDIV_tree_disturbance100m, 
                                                 "Buffer = 100m radius", 
                                                 death.in = "natural mortality")), 
  tar_target(plot_disturbedArea_perMortalityRate_200m, 
           plot_areaDisturbance_perMortalityRate(FUNDIV_tree_disturbance200m, 
                                                 "Buffer = 200m radius", 
                                                 death.in = "natural mortality"))
)
