####################################################################
# Download the tiff and headers files for Tmin, Tmax and Prec
####################################################################

# Create the directories
dirclimate <- paste0(getwd(), "/data/CLIMATE/")
if(!dir.exists(paste(dirclimate, country, sep = "/"))) dir.create(paste(dirclimate, country, sep = "/"))
for(a in 1:length(tiles)){
  if(!dir.exists(paste(dirclimate, country, tiles[a], sep = "/"))){
    dir.create(paste(dirclimate, country, tiles[a], sep = "/"))
  }
}

# Download data
print("Downloading data")
for(a in 1:length(tiles)){
  print(tiles[a])
  for(y in survey_dates){
    print(y)
    url <- paste("ftp://palantir.boku.ac.at/Public/ClimateData/TiledClimateData/", gsub("_", "/", tiles[a]),"/Tmax",y,"_",tiles[a],".tif", sep="")
    output_file <- paste(getwd(), "/data/CLIMATE/", country,"/", tiles[a],"/Tmax",y,"_",tiles[a],".tif", sep="")
    try(GET(url, authenticate('guest', ""), write_disk(output_file, overwrite = TRUE)))
    
    url <- paste("ftp://palantir.boku.ac.at/Public/ClimateData/TiledClimateData/", gsub("_", "/", tiles[a]),"/Tmin",y,"_",tiles[a],".tif", sep="")
    output_file <- paste(getwd(), "/data/CLIMATE/", country,"/", tiles[a],"/Tmin",y,"_",tiles[a],".tif", sep="")
    try(GET(url, authenticate('guest', ""), write_disk(output_file, overwrite = TRUE)))
    
    url <- paste("ftp://palantir.boku.ac.at/Public/ClimateData/TiledClimateData/", gsub("_", "/", tiles[a]),"/Tmax",y,"_",tiles[a],".hdr", sep="")
    output_file <- paste(getwd(), "/data/CLIMATE/", country,"/", tiles[a],"/Tmax",y,"_",tiles[a],".hdr", sep="")
    try(GET(url, authenticate('guest', ""), write_disk(output_file, overwrite = TRUE)))
    
    url <- paste("ftp://palantir.boku.ac.at/Public/ClimateData/TiledClimateData/", gsub("_", "/", tiles[a]),"/Tmin",y,"_",tiles[a],".hdr", sep="")
    output_file <- paste(getwd(), "/data/CLIMATE/", country,"/", tiles[a],"/Tmin",y,"_",tiles[a],".hdr", sep="")
    try(GET(url, authenticate('guest', ""), write_disk(output_file, overwrite = TRUE)))
    
    url <- paste("ftp://palantir.boku.ac.at/Public/ClimateData/TiledClimateData/", gsub("_", "/", tiles[a]),"/Prcp",y,"_",tiles[a],".hdr", sep="")
    output_file <- paste(getwd(), "/data/CLIMATE/", country,"/", tiles[a],"/Prcp",y,"_",tiles[a],".hdr", sep="")
    try(GET(url, authenticate('guest', ""), write_disk(output_file, overwrite = FALSE)))
    
    url <- paste("ftp://palantir.boku.ac.at/Public/ClimateData/TiledClimateData/", gsub("_", "/", tiles[a]),"/Prcp",y,"_",tiles[a],".tif", sep="")
    output_file <- paste(getwd(), "/data/CLIMATE/", country,"/", tiles[a],"/Prcp",y,"_",tiles[a],".tif", sep="")
    try(GET(url, authenticate('guest', ""), write_disk(output_file, overwrite = FALSE)))
  }
}

# calculate mean temperature for each tile and year (days are in the band of the tif)
print("calculate mean temperature for each tile and year")
for(a in 1:length(tiles)){
  print(tiles[a])
  for(y in survey_dates){
    print(y)
    # calculate mean temperature for each tile and year (days are in the band of the tif)
    yr.max <- stack(paste(getwd(), "/data/CLIMATE/", country,"/",tiles[a],"/Tmax",y,"_",tiles[a],".tif",sep=""))
    yr.min <- stack(paste(getwd(), "/data/CLIMATE/", country,"/",tiles[a],"/Tmin",y,"_",tiles[a],".tif",sep=""))
    yr.daily.mean <- ((yr.max + yr.min)/2)/100
    if(!dir.exists(paste0(dirclimate, country, "/Tmean/"))) dir.create(paste0(dirclimate, country, "/Tmean/"))
    if(!dir.exists(paste0(dirclimate, country, "/Tmean/daily/"))) dir.create(paste0(dirclimate, country, "/Tmean/daily/"))
    save(yr.daily.mean, file=paste(dirclimate, country,"/Tmean/daily/",y,"_",tiles[a],".RData",sep=""))
    
    # calculate mean daily temperature
    yr.mean <- mean(yr.daily.mean)
    if(!dir.exists(paste0(dirclimate, country, "/Tmean/yearly/"))) dir.create(paste0(dirclimate, country, "/Tmean/yearly/"))
    save(yr.mean, file=paste(dirclimate, country,"/Tmean/yearly/",y,"_",tiles[a],".RData",sep=""))
  }
}

for(a in 1:length(tiles)){
  mat_stack <- stack()
  for(y in survey_dates){
    load(file=paste(dirclimate, country,"/Tmean/yearly/",y,"_",tiles[a],".RData",sep=""))
    yr.mean[yr.mean < -80] <- NA
    mat_stack <- stack(mat_stack, yr.mean)
  }
  save(mat_stack, file=paste(dirclimate, country,"/Tmean/yearly/mat_", tiles[a],".RData",sep=""))
}

# Mean annual precipitation
print("Mean annual precipitations")
for(a in 1:length(tiles)){
  map_stack <- stack()
  print(tiles[a])
  for(y in survey_dates){
    print(y)
    # calculate sum of precipitation (days are in the band of the tif)
    prec <- stack(paste(dirclimate, country,"/",tiles[a],"/Prcp",y,"_",tiles[a],".tif",sep=""))
    yr.sum <- sum(prec)
    yr.sum[yr.sum < -80] <- NA
    yr.sum <- yr.sum/100
    map_stack <- stack(map_stack, yr.sum)
  }
  if(!dir.exists(paste0(dirclimate, country, "/Prec/"))) dir.create(paste0(dirclimate, country, "/Prec/"))
  if(!dir.exists(paste0(dirclimate, country, "/Prec/map/"))) dir.create(paste0(dirclimate, country, "/Prec/map/"))
  save(map_stack, file=paste(dirclimate, country,"/Prec/map/yearly_",tiles[a],".RData",sep=""))
}

# sum the degrees over 5.5 in a year
print("sum the degrees over 5.5 in a year")
for(a in 1:length(tiles)){
  print(tiles[a])
  sgdd_stack <- stack()
  for(y in survey_dates){
    print(y)
    load(file=paste(dirclimate, country,"/Tmean/daily/",y,"_",tiles[a],".RData",sep=""))
    sgdd <- calc(yr.daily.mean, sgdd_fun)
    sgdd <- sum(sgdd)
    sgdd_stack <- stack(sgdd_stack, sgdd)
  }
  if(!dir.exists(paste0(dirclimate, country, "/Tmean/sgdd/"))) dir.create(paste0(dirclimate, country, "/Tmean/sgdd/"))
  save(sgdd_stack, file=paste(dirclimate, country,"/Tmean/sgdd/yearly_",tiles[a],".RData",sep=""))
}


# spring frosts
print("spring frosts")
for(a in 1:length(tiles)){
  print(tiles[a])
  sf_stack <- stack()
  for(y in survey_dates){
    print(y)
    yr.min <- stack(paste(dirclimate, country,"/",tiles[a],"/Tmin",y,"_",tiles[a],".tif",sep=""))
    spring.frost <- calc(yr.min[[61:152]], spring_frost_fun)
    spring.frost <- sum(spring.frost)
    sf_stack <- stack(sf_stack, spring.frost)
  }
  if(!dir.exists(paste0(dirclimate, country, "/Frost/"))) dir.create(paste0(dirclimate, country, "/Frost/"))
  save(sf_stack, file=paste(dirclimate, country,"/Frost/yearly_",tiles[a],".RData",sep=""))
}


# calculate yearly precipitation for the summmer months (mm/year) - for SWS
print("calculate yearly precipitation for the summmer months (mm/year) - for SWS")
for(a in 1:length(tiles)){
  summer_p_stack <- stack()
  print(tiles[a])
  for(y in survey_dates){
    print(y)
    #get the date from the index of the layer
    prec <- stack(paste(dirclimate, country,"/",tiles[a],"/Prcp",y,"_",tiles[a],".tif",sep=""))
    m <- vector()
    for(i in 1:nlayers(prec)){
      m[i] <- format(strptime(paste(y, i), format="%Y %j"), format = "%m")
    }
    indices <- as.numeric(m)
    
    prec_summer <- prec[[which(indices %in% c(6,7,8))]]
    prec_summer[prec_summer < -800] <- NA
    prec_summer <- prec_summer/100
    prec_summer <- sum(prec_summer)
    
    summer_p_stack <- stack(summer_p_stack, prec_summer)
  }
  save(summer_p_stack, file=paste(dirclimate, country,"/Prec/map/yearly_summer_",tiles[a],".RData",sep=""))
}
