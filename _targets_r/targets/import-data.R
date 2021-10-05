tar_target(data_imported, 
           read_ifn_data(variable = "arbres_morts_foret", 
                         path = "data/FrenchNFI", 
                         years = c(2012, 2013),
                         addyear = TRUE, 
                         zipped = TRUE))
