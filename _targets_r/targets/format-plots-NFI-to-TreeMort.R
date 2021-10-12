list(
  tar_target(NFI_stand_age, Compute_NFI_stand_age(NFI_tree)), 
  tar_target(TreeMort_plot, Format_plots_TreeMort(NFI_plot, NFI_ecological_data, NFI_plot_elevation, NFI_stand_age))
)
