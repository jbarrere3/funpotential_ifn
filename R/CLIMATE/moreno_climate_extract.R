#####################################################################################
# Script to extract the climate variables for the given country and survey years
#####################################################################################

# MAT
print("MAT")
for(i in 1:nrow(plots_country)){
  print(i)
  x <- plots_country$latitude[i];x
  y <- plots_country$longitude[i];y
  
  tile <- tile_extents$tile[(y>tile_extents$ymin & y<tile_extents$ymax) & (x>tile_extents$xmin & x<tile_extents$xmax)]
  print(tile)
  if(length(tile)>0){
    # load the climate stack for the tile
    load(file=paste(country,"/Tmean/yearly/mat_", tile,".RData",sep=""))
    plots_country$tile[i] <- tile
    plots_country$mat[i] <- extract_climate_variable_from_coords(plots_country[i,c('longitude','latitude')], 
                                                                 mat_stack,
                                                                 surveys_start, surveys_end, 
                                                                 plots_country$start_year[i],
                                                                 plots_country$end_year[i])
  }
}
save(plots_country, file=save_file_name)

# MAP
print("MAP")
for(a in 1:length(tiles)){
  print(tiles[a])
  # use the saved tile value
  tile_plots <- plots_country[plots_country$tile==tiles[a],]
  # load the climate stack for the tile
  load(file=paste(country,"/Prec/map/yearly_", tiles[a],".RData",sep=""))
  
  print(nrow(tile_plots))
  for(i in 1:nrow(tile_plots)){
    plotcode <- tile_plots$plotcode[i]
    if(!is.na(plotcode)){
      map <- extract_climate_variable_from_coords(tile_plots[i,c('longitude','latitude')], 
                                                  map_stack,
                                                  surveys_start, surveys_end, 
                                                  tile_plots$start_year[i],
                                                  tile_plots$end_year[i])
      plots_country$map[plots_country$plotcode==plotcode] <- map
    }
  }
  save(plots_country, file=save_file_name)
}

# SGDD
print("SGDD")
for(a in 1:length(tiles)){
  print(tiles[a])
  tile_plots <- plots_country[plots_country$tile==tiles[a],]
  # load the climate stack for the tile
  load(file=paste(country,"/Tmean/sgdd/yearly_", tiles[a],".RData",sep=""))
  
  print(nrow(tile_plots))
  for(i in 1:nrow(tile_plots)){
    plotcode <- tile_plots$plotcode[i]
    if(!is.na(plotcode)){
      sgdd <- extract_climate_variable_from_coords(tile_plots[i,c('longitude','latitude')], 
                                                   sgdd_stack,
                                                   surveys_start, surveys_end, 
                                                   tile_plots$start_year[i],
                                                   tile_plots$end_year[i])
      plots_country$sgdd[plots_country$plotcode==plotcode] <- sgdd
    }
  }
  save(plots_country, file=save_file_name)
}

# Spring frosts
print("Spring frosts")
for(a in 1:length(tiles)){
  print(tiles[a])
  tile_plots <- plots_country[plots_country$tile==tiles[a],]
  # load the climate stack for the tile
  load(file=paste(country,"/Frost/yearly_", tiles[a],".RData",sep=""))
  
  print(nrow(tile_plots))
  for(i in 1:nrow(tile_plots)){
    plotcode <- tile_plots$plotcode[i]
    if(!is.na(plotcode)){
      sf <- extract_climate_variable_from_coords(tile_plots[i,c('longitude','latitude')], 
                                                 sf_stack,
                                                 surveys_start, surveys_end, 
                                                 tile_plots$start_year[i],
                                                 tile_plots$end_year[i])
      plots_country$spring_frosts[plots_country$plotcode==plotcode] <- sf
    }
  }
  save(plots_country, file=save_file_name)
}

# WAI (yearly)
print("WAI")
for(a in 1:length(tiles)){
  print(tiles[a])
  tile_plots <- plots_country[plots_country$tile==tiles[a],]
  load(file=paste(country,"/Prec/map/yearly_", tiles[a],".RData",sep=""))
  
  print(nrow(tile_plots))
  for(i in 1:nrow(tile_plots)){
    plotcode <- tile_plots$plotcode[i]
    if(!is.na(plotcode)){
      if(is.na(tile_plots$wai_yearly[i])){
        wai <- extract_yearly_water_availablity_from_coords(tile_plots[i,c('longitude','latitude')], 
                                                            pet.country, map_stack,
                                                            surveys_start, surveys_end,
                                                            tile_plots$start_year[i]-2, 
                                                            tile_plots$end_year[i])
        plots_country$wai_yearly[plots_country$plotcode==plotcode] <- wai[1]
        plots_country$p_pet_yearly[plots_country$plotcode==plotcode] <- wai[2]
      }
    }
  }
}
save(plots_country, file=save_file_name)

# Summer water stress (yearly)
print("SWS")
for(a in 1:length(tiles)){
  print(tiles[a])
  tile_plots <- plots_country[plots_country$tile==tiles[a],]
  load(file=paste(country,"/Prec/map/yearly_summer_", tiles[a],".RData",sep=""))
  
  print(nrow(tile_plots))
  for(i in 1:nrow(tile_plots)){
    plotcode <- tile_plots$plotcode[i]
    if(!is.na(plotcode)){
      
      wai <- extract_summer_water_stress_from_coords(tile_plots[i,c('longitude','latitude')], 
                                                     pet.summer.country, 
                                                     summer_p_stack,
                                                     surveys_start, 
                                                     surveys_end,
                                                     tile_plots$start_year[i]-2, 
                                                     tile_plots$end_year[i])
      plots_country$sws_yearly[plots_country$plotcode==plotcode] <- wai[1]
      plots_country$p_pet_summer[plots_country$plotcode==plotcode] <- wai[2]
      
    }
  }
}
save(plots_country, file=save_file_name)

