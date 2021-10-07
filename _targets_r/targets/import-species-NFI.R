list(
  tar_target(NFI_species, read.csv("data/FrenchNFI/species.csv", sep = ",")), 
  tar_target(NFI_genus_species_correspondence, read.csv("data/TreeMort/correspondanceGenusFamilyNFI.csv", sep = ";"))
)
