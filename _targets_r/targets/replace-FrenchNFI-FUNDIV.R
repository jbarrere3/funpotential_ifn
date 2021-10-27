tar_target(FUNDIV_tree, 
           (FUNDIV_tree_original %>% 
             filter(country != "FR") %>% 
             rbind(FUNDIV_FrenchNFI_tree)))
