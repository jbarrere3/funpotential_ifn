list(
  tar_target(data_jags_generated, 
           generate_data_jags(10, 1000, data_jags_Aalba)), 
  tar_target(data_jags_generated_2, 
           generate_data_jags_2(10, 1000, data_jags_Aalba_2))
)
