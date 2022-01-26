list(
  tar_target(jags_simulated, 
             fit_mortality(data_jags_generated$data, n.chains = 3, n.iter = 5000, n.burn = 1000, n.thin = 1)), 
  tar_target(jags_simulated_2, 
             fit_mortality_2(data_jags_generated_2$data, n.chains = 3, n.iter = 5000, n.burn = 3000, n.thin = 1)))

