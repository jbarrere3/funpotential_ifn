list(
  tar_target(FUNDIV_tree_original, read_FUNDIV_tree_data(data_path = "data/FUNDIV", remove_harv = TRUE)), 
  tar_target(FUNDIV_plot, read_FUNDIV_plot_data(data_path = "data/FUNDIV"))
)
