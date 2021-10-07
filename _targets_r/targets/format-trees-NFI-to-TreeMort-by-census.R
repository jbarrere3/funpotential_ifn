list(
  tar_target(TreeMort_tree_census1, Format_trees_census1_TreeMort(NFI_tree, TreeMort_species)),
  tar_target(TreeMort_tree_census2, Format_trees_census2_TreeMort(NFI_tree_remeasured, TreeMort_tree_census1))
)
