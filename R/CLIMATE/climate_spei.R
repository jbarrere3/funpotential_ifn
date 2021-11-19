setwd("~/Documents/Projects/sAPROPOS/data")

library("sp")
library("ncdf4")
library("chron")
library("RColorBrewer")
library("lattice")
library("raster")
library("rgdal")
library("maptools")
library("matrixStats")

plots <- read.csv(file="FunDiv_plots_Nadja.csv", header=TRUE, stringsAsFactors=FALSE)
# remove the plots in the Canary Islands
plots <- plots[plots$latitude>30,]

plots_all <- plots[,c("plotcode",'country','longitude','latitude', 'surveydate1','surveydate2')]
plots_all$start_year <- as.numeric(format(as.Date(plots_all$surveydate1),"%Y"))
plots_all$end_year <- as.numeric(format(as.Date(plots_all$surveydate2),"%Y"))

ncpath <- "/Volumes/FREESIAS/Data/Spatial/Climate/CSIC/SPEI/"
ncname <- "spei12"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "spei"  # note: tmp means temperature (not temporary)
# open a netCDF file
ncin <- nc_open(ncfname)
print(ncin)
names(ncin)

# get longitude and latitude
lon <- ncvar_get(ncin, "lon", start=c(319), count=c(111))
dim(lon)
lat <- ncvar_get(ncin, "lat", start=c(227), count=c(98))
dim(lat)
# get time
time <- ncvar_get(ncin, "time", start=c(960), count=c(420))
dim(time)

##################################################
## Get the input variable (SPEI) and its attributes, 
## and verify the size of the array
spei.array <- ncvar_get(ncin, dname, start=c(319, 227, 960),
                        count=c(111, 98, 420), verbose = TRUE)
dlname <- ncatt_get(ncin, dname,"long_name")
dunits <- ncatt_get(ncin, dname, "units")
fillvalue <- ncatt_get(ncin, dname, "_FillValue")
dim(spei.array) 

#################################################
###CONVERSION OF THE NETCDF FILES TO DATA.FRAMES
#################################################
# Replace netCDF fillvalues with R NAs
#values of a variable that are either missing or simply not available 
#(i.e. ocean grid points in a terrestrial data set) are flagged 
#using specific fill values (_FillValue) or missing values (missing_value)
spei.array[spei.array == fillvalue$value] <- NA
length(na.omit(as.vector(spei.array[, , 1])))


#################################################################################
##Convert the whole array to a data frame, and calculate the annual mean
#create a long vector tmp.vec.long using the as.vector() reshaping function, 
#verify its length, which should be 4568760
spei.vec.long <- as.vector(spei.array)
length(spei.vec.long)

#Then reshape that vector into a 4568760 by 12 matrix using the matrix() function, 
##verify its dimensions, which should be 4568760 by 12.

nlon <- dim(lon)
nlat <- dim(lat)
nt <- dim(time)

#spei ann has dimension of 10878 (no. pixels) by 420
str(spei.vec.long)

spei.ann <- matrix(spei.vec.long, nrow = nlon * nlat, ncol = nt)
dim(spei.ann)
str(spei.ann)

#there is the montly information for each pixel by 12 months by 35 years
head(na.omit(spei.ann))
12*35
12*30
#Create data frame from the spei.ann matrix.
lonlat <- as.matrix(expand.grid(lon, lat))
spei.df02 <- data.frame(cbind(lonlat, spei.ann))
mm<-c("speiJan", "speiFeb", "speiMar", "speiApr", "speiMay", 
      "speiJun", "speiJul", "speiAug", "speiSep", "speiOct", "speiNov", "speiDec")
mmm<-rep(mm,35)
names(spei.df02) <- c("lon", "lat", mmm)
options(width = 96)
head(na.omit(spei.df02, 20))

#Get annual mean
str(spei.df02)
dim(spei.df02)
ss <- spei.df02[,3:422]

b <- matrix(nrow=35, ncol=10878 ) 
b <-apply(b, 1,as.numeric)
year_data <- data.frame(b)
colnames(year_data) <- paste("yr", 1:35, sep = "")
str(year_data)

# loop through the 35 years
for(i in 1:35){
  z <- i*12
  j <- z-11
  aa <- ss[,j:z]
  str(aa)
  year_data[i] <- rowMeans(aa)
  summary(year_data[,i])
}

long <- spei.df02[1]
lat <- spei.df02[2]

final_data <- data.frame(long, lat, year_data)

save(final_data, file="SPEI12_year81_15.RData")

# match the spei to the plots
coords <- data.frame(lon=plots_all$longitude,lat=plots_all$latitude)

dfr <- rasterFromXYZ(final_data) 
SPEIplot <- extract(dfr, coords)
SPEIplot <- as.data.frame(SPEIplot)

# select the years and bind to the plot data
SPEIplot <- SPEIplot[,c(1:35)]
spei_plots_all <- cbind(plots_all, SPEIplot)
colnames(spei_plots_all)[9:43] <- 1981:2015

# calculate the min and mean spei of the years between the two surveys
spei_plots_all$min_spei_survey_years <- NA
spei_plots_all$mean_spei_survey_years <- NA
# French plots missing start year
spei_plots_all$start_year[spei_plots_all$country=='FG'] <- spei_plots_all$end_year[spei_plots_all$country=='FG'] - 7

for(i in 1:nrow(spei_plots_all)){
  if(is.na(spei_plots_all$mean_spei_survey_years[i])){
    print(spei_plots_all$plotcode[i])
    yrs <- as.character(spei_plots_all$start_year[i]:spei_plots_all$end_year[i])
    
    spei_plots_all$mean_spei_survey_years[i] <- rowMeans(spei_plots_all[i, yrs], na.rm=TRUE)
    spei_plots_all$min_spei_survey_years[i] <- rowMins(as.matrix(spei_plots_all[i, yrs]), na.rm=TRUE)
  }
}

save(spei_plots_all, file="climate_data_spei_plots.RData")
