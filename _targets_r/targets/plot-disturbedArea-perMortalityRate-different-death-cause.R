list(
  tar_target(plot_disturbedArea_perMortalityRate_200m_harvest, 
           plot_areaDisturbance_perMortalityRate(FUNDIV_tree_disturbance200m, 
                                                 "Buffer = 200m radius", 
                                                 death.in = "harvested")), 
  tar_target(plot_disturbedArea_perMortalityRate_200m_unknown, 
           plot_areaDisturbance_perMortalityRate(FUNDIV_tree_disturbance200m, 
                                                 "Buffer = 200m radius", 
                                                 death.in = "unknown cause")), 
  tar_target(plot_disturbedArea_perMortalityRate_200m_alldeath, 
           plot_areaDisturbance_perMortalityRate(FUNDIV_tree_disturbance200m, 
                                                 "Buffer = 200m radius", 
                                                 death.in = c("natural mortality", "harvested", "unknown cause")))
)
