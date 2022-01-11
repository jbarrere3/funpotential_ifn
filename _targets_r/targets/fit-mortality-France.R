list(
  tar_target(fit.mortality.FR_A.alba, 
             fit_mortality_FR(FUNDIV_tree, FUNDIV_climate, disturbance_per_plot_200m, 
                             NFI_plot_remeasure, species.in = "Abies alba", 
                             n.chains = 3, n.iter = 500, n.burn = 100, n.thin = 1)), 
  tar_target(fit.mortality.FR_P.sylvestris, 
             fit_mortality_FR(FUNDIV_tree, FUNDIV_climate, disturbance_per_plot_200m, 
                             NFI_plot_remeasure, species.in = "Pinus sylvestris", 
                             n.chains = 3, n.iter = 500, n.burn = 100, n.thin = 1))
)

