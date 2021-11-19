
extract_varX <- function(in_df, varX_raster, varX="mat", buff_cor=NULL){
  df <- cbind(in_df[is.na(in_df[[varX]]),c("plotcode",'longitude','latitude')], 
              extract(varX_raster,
                      in_df[is.na(in_df[[varX]]),c("longitude","latitude")],
                      buffer=buff_cor, fun=mean))
  df <- as.data.frame(df)
  colnames(df) <- c("plotcode",'longitude','latitude', varX)
  return(df)
}

extract_climate_variable_from_coords <- function(coords, climate_yearly_stack, 
                                                 stack_start_Year, stack_end_Year,
                                                 start_Year, end_Year){
  climate_extract <- extract(climate_yearly_stack, coords)
  if(is.na(climate_extract[1])){
    for(b in seq(from=5000, to=25000, by=5000)){
      print(b)
      climate_extract <- extract(climate_yearly_stack, coords, buffer=b, fun=mean)
      if(!is.na(climate_extract[1])) break
    }
  }
  if(!is.na(climate_extract[1])){
    # take the mean of the survey period
    yearly_climate <- as.data.frame(cbind(as.vector(climate_extract), 
                                          as.vector(stack_start_Year:stack_end_Year)))
    colnames(yearly_climate) <- c('climate_extract','year')
    start_yr <- start_Year-2
    survey_period <- seq(start_yr, end_Year)
    
    return(mean(yearly_climate$climate_extract[yearly_climate$year %in% survey_period]))
  }
  return(NA)
}

extract_water_availablity_from_coords <- function(coords, pet_stack, prec_stack, start_Year, end_Year){
  
  pet <- extract(pet_stack, coords)
  prec <- extract(prec_stack, coords)
  if(is.na(prec[1])){
    for(b in seq(from=5000, to=25000, by=5000)){
      print(b)
      prec <- extract(prec_stack, coords, buffer=b, fun=mean)
      if(!is.na(prec[1])) break
    }
  }
  monthly_wai <- as.data.frame(cbind(as.vector(pet), as.vector(prec)))
  colnames(monthly_wai) <- c("monthly_pet","monthly_prec")
  monthly_wai$wai <- (monthly_wai$monthly_prec-monthly_wai$monthly_pet)/monthly_wai$monthly_pet
  monthly_wai$ws <- monthly_wai$monthly_prec/monthly_wai$monthly_pet
  monthly_wai$year <- as.numeric(substr(names(s), 2, 5))
  monthly_wai$month <- as.numeric(substr(names(s), 7, length(names(s))))
  
  monthly_wai <- monthly_wai[monthly_wai$monthly_pet>0,]  
  survey_period <- seq(start_Year, end_Year)
  monthly_wai <- monthly_wai[monthly_wai$year %in% survey_period,]
  return(cbind(mean(monthly_wai$wai, na.rm=T), 
               mean(monthly_wai$ws[monthly_wai$month %in% c(6,7,8)], na.rm=T)))
}


extract_summer_water_stress_from_coords <- 
  function(coords, pet_stack, prec_stack, stack_start_Year, stack_end_Year, 
           start_Year, end_Year){
    
    pet <- extract(pet_stack, coords)
    prec <- extract(prec_stack, coords)
    if(is.na(prec[1])){
      for(b in seq(from=5000, to=25000, by=5000)){
        print(b)
        prec <- extract(prec_stack, coords, buffer=b, fun=mean)
        if(!is.na(prec[1])) break
      }
    }
    summer_ws <- as.data.frame(cbind(as.vector(pet), as.vector(prec),
                                      as.vector(stack_start_Year:stack_end_Year)))
    colnames(summer_ws) <- c("pet","prec", "year")
    summer_ws$sws <- summer_ws$prec/summer_ws$pet
    summer_ws$p_pet <- summer_ws$prec-summer_ws$pet
    
    survey_period <- seq(start_Year, end_Year)
    summer_ws <- summer_ws[summer_ws$year %in% survey_period,]
    return(cbind(mean(summer_ws$sws, na.rm=T), mean(summer_ws$p_pet, na.rm=T)))
  }


extract_yearly_water_availablity_from_coords <- 
  function(coords, pet_stack, prec_stack, stack_start_Year, stack_end_Year, 
           start_Year, end_Year){
    
    pet <- extract(pet_stack, coords)
    prec <- extract(prec_stack, coords)
    if(is.na(prec[1])){
      for(b in seq(from=5000, to=25000, by=5000)){
        print(b)
        prec <- extract(prec_stack, coords, buffer=b, fun=mean)
        if(!is.na(prec[1])) break
      }
    }
    yearly_wai <- as.data.frame(cbind(as.vector(pet), as.vector(prec),
                                      as.vector(stack_start_Year:stack_end_Year)))
    colnames(yearly_wai) <- c("pet","prec", "year")
    yearly_wai$wai <- (yearly_wai$prec-yearly_wai$pet)/yearly_wai$pet
    yearly_wai$p_pet <- yearly_wai$prec-yearly_wai$pet
    
    survey_period <- seq(start_Year, end_Year)
    yearly_wai <- yearly_wai[yearly_wai$year %in% survey_period,]
    return(cbind(mean(yearly_wai$wai, na.rm=T), mean(yearly_wai$p_pet, na.rm=T)))
  }


# sum the degrees over 5.5 in a year
sgdd_fun <- function(x){ 
  x[x < -80] <- NA; 
  x[!is.na(x)] <- x[!is.na(x)]-5.5;
  x[x < 0] <- 0; 
  return(x)
}

spring_frost_fun <- function(x){ 
  x[x < -8000] <- NA; 
  x[!is.na(x)] <- ifelse(x[!is.na(x)] < 0, 1, 0) 
  return(x)
}
