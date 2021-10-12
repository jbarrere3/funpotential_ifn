tar_target(TreeMort_tree, 
           rbind(subset(TreeMort_tree_census1, 
                        plot.id %in% TreeMort_tree_census2$plot.id),
                 TreeMort_tree_census2))
