list(
  tar_target(NFI_tree, merge_dead_alive(NFI_tree_alive, NFI_tree_dead)),
  tar_target(NFI_tree_remeasured, 
             merge_NFI_remeasured(NFI_tree_alive_remeasure, NFI_tree_dead_remeasure))
)
