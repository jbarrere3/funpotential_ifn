setwd("~/Documents/Projects/sAPROPOS/data")

########################################################################
# Based on inventory plot survey dates (Moreno data)
########################################################################

load("climate/moreno/spain/climate_data.RData")
plots_es <- plots_country
load("climate/moreno/france/climate_data.RData")
plots_fg <- plots_country
# check that start and end dates aren't reversed
plots_fg <- plots_fg[,c(1:5,7,6,8:17)]
load("climate/moreno/germany/climate_data.RData")
plots_de <- plots_country
load("climate/moreno/finland/climate_data.RData")
plots_fi <- plots_country
load("climate/moreno/sweden/climate_data.RData")
load("climate/moreno/wallonia/climate_data.RData")
plots_wa <- plots_country

plots_fi$country <- 'FI'
plots_sw$country <- 'SW'
plots_de$country <- 'DE'
plots_fg$country <- 'FR'
plots_es$country <- 'ES'
plots_wa$country <- 'WA'

colnames(plots_fi);colnames(plots_sw);colnames(plots_de);colnames(plots_fg);colnames(plots_es);colnames(plots_wa)
moreno_climate_data <- rbind(plots_de, plots_es, plots_fg, plots_fi, plots_sw, plots_wa)

###########################
# Original world clim data
###########################

climate.wc <- read.csv(file="FunDiv_plots_climate_and_dem_Nadja.csv", header=TRUE)
head(climate.wc[climate.wc$country=='SW',])
head(climate.wc[climate.wc$country=='FG',])
head(climate.wc[climate.wc$country=='WA',])
climate.wc <- climate.wc[,c("plotcode","bio1","bio12","wai")]
climate.wc$mat.wc <- climate.wc$bio1/10
colnames(climate.wc)[3:4] <- c("map.wc","wai.wc")

climate_plots <- merge(moreno_climate_data, climate.wc, by=c("plotcode"), all.x=TRUE)

save(climate_plots, file="sapropos_climate_data_dec17.RData")
write.csv(climate_plots, file="FunDiv_plots_climate_sapropos_dec17.csv", row.names = F)

colnames(climate_plots)
columns <- c(1:8,17:21,9:16)
climate_plots <- climate_plots[,columns]
colnames(climate_plots)
colnames(climate_plots)[18:21] <- c('wai','p_pet_yr','sws','p_pet_summer')
save(climate_plots, file="sapropos_climate_data_dec17.RData")

write.csv(climate_plots, file="FunDiv_plots_climate_sapropos_dec17.csv", row.names = F)


# merge the SPEI 0.5 degree data
load(file="sapropos_climate_data_dec17.RData")
load(file="climate_data_spei_plots.RData")

climate_plots <- merge(climate_plots, spei_plots_all[,c('plotcode','min_spei_survey_years','mean_spei_survey_years')], by=c("plotcode"), all.x=TRUE)
colnames(climate_plots)[22:23] <- c('spei_min','spei_mean')

save(climate_plots, file="sapropos_climate_data_dec17.RData")

write.csv(climate_plots, file="FunDiv_plots_climate_sapropos_dec17.csv", row.names = F)

###########################
# Plotting
###########################
library(raster)

# points from scratch
climate_plots <- climate_plots[!is.nan(climate_plots$spei_mean),]
coords = cbind(climate_plots$longitude, climate_plots$latitude)
# make spatial data frame
spdf = SpatialPointsDataFrame(coords, climate_plots)

require(maptools)
require(classInt)
require(RColorBrewer)
library(rgeos)
library(rgdal)

load(file="raster/europe_all_small.RData")
europe.all <- europe.all.small
cp <- as(extent(-10.47757, 32, 35.26787, 70.07531),"SpatialPolygons")
proj4string(cp) <- CRS(proj4string(europe.all))
eur <- gIntersection(europe.all, cp, byid=TRUE)
europe.all <- eur

plotclr <- rev(brewer.pal(5, "Spectral"))

png("../results/figures/climate_survey_periods.png", width=15, height=20, units="cm", res=900)
par(mfrow = c(3, 2), oma=c(1, 1, 1, 1))
par(mar=c(0, 0, 0, 0) + 0.1)

# MAT 
breaks=c(round(quantile(spdf$mat, probs = seq(0, 1, by = 0.1),na.rm=T),2))
cuts <- classIntervals(spdf$mat, style="fixed",
                       fixedBreaks=breaks)
colcode <- findColours(cuts, plotclr)

plot(spdf, col=colcode, pch=19, cex=0.2)
legend("topleft", legend=names(attr(colcode, "table")), 
       fill=attr(colcode, "palette"), bty="n", cex=0.8, title="MAT") 
plot(europe.all, border="dark grey", add=TRUE, lwd=0.4)

# MAP 
breaks=c(round(quantile(spdf$map, probs = seq(0, 1, by = 0.1),na.rm=T),2))
cuts <- classIntervals(spdf$map, style="fixed",
                       fixedBreaks=breaks)
colcode <- findColours(cuts, plotclr)

plot(spdf, col=colcode, pch=19, cex=0.2)
legend("topleft", legend=names(attr(colcode, "table")), 
       fill=attr(colcode, "palette"), bty="n", cex=0.8, title="MAP") 
plot(europe.all, border="dark grey", add=TRUE, lwd=0.4)

# sum growing degrees days
breaks=c(round(quantile(spdf$sgdd, probs = seq(0, 1, by = 0.1),na.rm=T),2))
cuts <- classIntervals(spdf$sgdd, style="fixed",
                       fixedBreaks=breaks)
colcode <- findColours(cuts, plotclr)

plot(spdf, col=colcode, pch=19, cex=0.2)
legend("topleft", legend=names(attr(colcode, "table")), 
       fill=attr(colcode, "palette"), bty="n", cex=0.8, title="SGDD") 
plot(europe.all, border="dark grey", add=TRUE, lwd=0.4)

# Spring frosts
breaks=c(round(quantile(spdf$spring_frosts, probs = seq(0, 1, by = 0.1),na.rm=T),2))
cuts <- classIntervals(spdf$spring_frosts, style="fixed",
                       fixedBreaks=breaks)
colcode <- findColours(cuts, plotclr)

plot(spdf, col=colcode, pch=19, cex=0.2)
legend("topleft", legend=names(attr(colcode, "table")), 
       fill=attr(colcode, "palette"), bty="n", cex=0.8, title="Spring frosts (days)") 
plot(europe.all, border="dark grey", add=TRUE, lwd=0.4)


# WAI
breaks=c(round(quantile(spdf$wai, probs = seq(0, 1, by = 0.1),na.rm=T),2))
cuts <- classIntervals(spdf$wai, style="fixed",
                       fixedBreaks=breaks)
colcode <- findColours(cuts, plotclr)

plot(spdf, col=colcode, pch=19, cex=0.2)
legend("topleft", legend=names(attr(colcode, "table")), 
       fill=attr(colcode, "palette"), bty="n", cex=0.8, title="WAI") 
plot(europe.all, border="dark grey", add=TRUE, lwd=0.4)


# SWS
breaks=c(round(quantile(spdf$sws, probs = seq(0, 1, by = 0.1),na.rm=T),2))
cuts <- classIntervals(spdf$sws, style="fixed",
                       fixedBreaks=breaks)
colcode <- findColours(cuts, plotclr)

plot(spdf, col=colcode, pch=19, cex=0.2)
legend("topleft", legend=names(attr(colcode, "table")), 
       fill=attr(colcode, "palette"), bty="n", cex=0.8, title="SWS") 
plot(europe.all, border="dark grey", add=TRUE, lwd=0.4)

# mean SPEI 12
breaks=c(round(quantile(spdf$spei_mean, probs = seq(0, 1, by = 0.1),na.rm=T),2))
cuts <- classIntervals(spdf$spei_mean, style="fixed",
                       fixedBreaks=breaks)
colcode <- findColours(cuts, plotclr)

plot(spdf, col=colcode, pch=19, cex=0.2)
legend("topleft", legend=names(attr(colcode, "table")), 
       fill=attr(colcode, "palette"), bty="n", cex=0.8, title="Mean SPEI") 
plot(europe.all, border="dark grey", add=TRUE, lwd=0.4)



# min SPEI 12
breaks=c(round(quantile(spdf$spei_min, probs = seq(0, 1, by = 0.1), na.rm=T),2))
cuts <- classIntervals(spdf$spei_min, style="fixed",
                       fixedBreaks=breaks)
colcode <- findColours(cuts, plotclr)

plot(spdf, col=colcode, pch=19, cex=0.2)
legend("topleft", legend=names(attr(colcode, "table")), 
       fill=attr(colcode, "palette"), bty="n", cex=0.8, title="Minimum SPEI") 
plot(europe.all, border="dark grey", add=TRUE, lwd=0.4)

dev.off()


countries <- levels(as.factor(climate_plots$country))
png("../results/figures/climate_compare_wc_moreno.png", width=15, height=20, units="cm", res=900)
par(mfrow = c(3, 1))

plot(climate_plots$mat, climate_plots$mat.wc, pch=16, cex=0.6, col=as.factor(climate_plots$country),
     main="MAT", ylab="WC MAT", xlab="Survey period MAT")
legend("topright", legend=countries, col=as.factor(countries), pch=16)
plot(climate_plots$map, climate_plots$map.wc, pch=16, cex=0.6, col=as.factor(climate_plots$country),
     main="MAP (mm/yr)", ylab="WC MAP", xlab="Survey period MAP")
legend("topright", legend=countries, col=as.factor(countries), pch=16)
plot(climate_plots$wai, climate_plots$wai.wc, pch=16, cex=0.6, col=as.factor(climate_plots$country),
     main="WAI", ylab="WC WAI (mm/yr)", xlab="Survey period WAI  (monthly mm/day)")
legend("topright", legend=countries, col=as.factor(countries), pch=16)

dev.off()

######################################
### WAI and SWS based on CRU
######################################
load(file="climate_wai_max_stress_survey_dates.RData")
#climate.water$wai <- as.numeric(as.character(climate.water$wai))
#climate.water$sws <- as.numeric(as.character(climate.water$sws))
#colnames(climate.water) <- c('plotcode','wai_cru','sws_cru')

load(file="sapropos_climate_data_oct17.RData")

#clim_dat <- merge(climate_plots, climate.water, by=c("plotcode"))

countries <- levels(as.factor(climate_plots$country))
png("../results/figures/climate_compare_CRU_moreno.png", width=15, height=20, units="cm", res=900)
par(mfrow = c(2, 1))
plot(climate_plots$wai, climate_plots$wai_cru, pch=16, cex=0.6, col=as.factor(climate_plots$country),
     main="WAI", ylab="CRU", xlab="CRU and Moreno")
legend("topleft", legend=countries, col=as.factor(countries), pch=16)
plot(climate_plots$sws, climate_plots$sws_cru, pch=16, cex=0.6, col=as.factor(climate_plots$country),
     main="Summer water stress", ylab="CRU", xlab="CRU and Moreno")
legend("topleft", legend=countries, col=as.factor(countries), pch=16)

dev.off()


plot(climate_plots$wai.wc, climate_plots$wai_cru, pch=16, cex=0.6, col=as.factor(climate_plots$country),
     main="WAI", ylab="CRU", xlab="WorldClim")
legend("topleft", legend=countries, col=as.factor(countries), pch=16)

png("../results/figures/climate_wai_wc_moreno_yearly.png", width=7, height=10, units="cm", res=900)

countries <- levels(as.factor(climate_plots$country))
plot(climate_plots$wai.wc, climate_plots$wai_yearly, pch=16, cex=0.6, col=as.factor(climate_plots$country),
     main="WAI (mm/yr)", ylab="CRU and Moreno", xlab="WorldClim")
legend("topleft", legend=countries, col=as.factor(countries), pch=16)
abline(0,1)

dev.off()

png("../results/figures/climate_wai_and_sws.png", 15, height=20, units="cm", res=900)
par(mfrow = c(2, 1))
countries <- levels(as.factor(climate_plots$country))
plot(climate_plots$wai, climate_plots$sws, pch=16, cex=0.6, col=as.factor(climate_plots$country),
     main="WAI and SWS", ylab="SWS", xlab="WAI")
legend("topleft", legend=countries, col=as.factor(countries), pch=16)

countries <- levels(as.factor(climate_plots$country))
plot(climate_plots$p_pet_yr, climate_plots$p_pet_summer, pch=16, cex=0.6, col=as.factor(climate_plots$country),
     main="MAP - PET", ylab="Summer month average", xlab="Yearly average")
legend("topleft", legend=countries, col=as.factor(countries), pch=16)

dev.off()

