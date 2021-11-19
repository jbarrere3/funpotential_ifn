###################################################################################
# France: Moreno climate data download, preparation and extraction for each plot
###################################################################################
#setwd("~/Documents/Projects/sAPROPOS/data/climate/moreno")
dirclimate <- paste0(getwd(), "/data/CLIMATE/")

library(raster)
rasterOptions(maxmemory = 1e+09)
library(RCurl)
library(httr)

source("R/CLIMATE/climate_data_prep_functions.R")

# French tiles
#tiles <- c('B_7','C_6','C_7','C_8','D_6','D_7','D_8','E_7','E_8','E_9')
tiles <- c('B_7','C_6')
survey_dates <- 2000:2006
country <- 'france'

#############################
# 1. Download and preparation
#############################
source("R/CLIMATE/moreno_climate_prep.R")

#######################################
# 2. Extract climate data for each plot
#######################################
load(file=paste0(dirclimate, "climate_moreno_tile_extents.RData"))

save_file_name <- paste(dirclimate, country, "/climate_data.RData", sep="")
surveys_start <- 2000
surveys_end <- 2006

pet.cru <- stack(paste0(dirclimate, "CRU/PET/cru_ts4.05.1981.1990.pet.dat.nc"))
pet.country <- pet.cru[[20:32]]

pet.summer.cru <- stack("/Volumes/FREESIAS/Data/Spatial/Climate/CRU/v3.24.01/PET/1981.2015.pet_summer.yearly.dat_eur.nc")
pet.summer.country <- pet.summer.cru[[20:32]]

### French plots
plots <- read.csv(file="../../FunDiv_plots_Nadja.csv", header=TRUE, stringsAsFactors=FALSE)
plots_country <- plots[plots$country=='FG',c("plotcode",'longitude','latitude', 'surveydate1','surveydate2')]
plots_country$end_year <- as.numeric(format(as.Date(plots_country$surveydate2),"%Y"))
plots_country$start_year <- plots_country$end_year -7
plots_country$tile <- NA
plots_country$mat <- NA; plots_country$sgdd <- NA; plots_country$spring_frosts <- NA; plots_country$map <- NA 
plots_country$wai_yearly <- NA
plots_country$p_pet_yearly <- NA
plots_country$sws_yearly <- NA
plots_country$p_pet_summer <- NA

source("../../../analysis/moreno_climate_extract.R")

