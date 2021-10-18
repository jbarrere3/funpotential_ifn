tar_target(TreeMort_tree, 
           rbind(subset(TreeMort_tree_census1, 
                        tree.id %in% TreeMort_tree_census2$tree.id),
                 TreeMort_tree_census2))
