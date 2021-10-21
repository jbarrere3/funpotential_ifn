read_FUNDIV_plot_data <- function(data_path = "data/FUNDIV"){
  plots <- read.csv(file.path(data_path, "FunDiv_plots_Nadja.csv"),
                    stringsAsFactors=FALSE)
  clims <- read.csv(file.path(data_path,
                              "FunDiv_plots_climate_sapropos.csv"),
                    stringsAsFactors=FALSE)
  clims <- clims %>% select(-longitude, -latitude, -country )

  # merging data
  data <- left_join( plots, clims, by = "plotcode")
  return(data)
}

read_FUNDIV_tree_data <- function(data_path = "data/FUNDIV", remove_harv = TRUE) {
  #function to read and select data
  require(dplyr)
  trees <- read.csv(file.path(data_path, "FunDiv_trees_Nadja.csv"),
                    stringsAsFactors=FALSE)
  # combine Betula pendula and pubescens with the Betula genus (ids 46 and 47 with 48)
  print("read tree")
  trees$speciesid[trees$speciesid %in% c(46,47)] <- 48
  species <- read.csv(file.path(data_path, "FunDiv_species_Nadja.csv"),
                      stringsAsFactors=FALSE, fileEncoding = "cp1252")
  # remove characters at the end of the species name
  species$species[species$id ==277] <- "pubescens"
  species <- species %>% mutate(sp = paste(genus, species)) %>%
      dplyr::select(c(id, sp))
  species[species$id == 277, "sp"] <- "Quercus pubescens"
  species[species$id == 48, "sp"] <- "Betula"
  plots <- read.csv(file.path(data_path, "FunDiv_plots_Nadja.csv"),
                    stringsAsFactors=FALSE)
   # not sure if this is only me but I have encoding problem wityh the species csv file
  plots <- plots %>%
      dplyr::select(-country, -ba_ha1, -ba_ha2,-surveydate1, -surveydate2)
  clims <- read.csv(file.path(data_path,
                              "FunDiv_plots_climate_sapropos.csv"),
                    stringsAsFactors=FALSE)
    print("read climatic data")
  clims <- clims %>% select(-longitude, -latitude, -country )
  # merging data
  data <- left_join(trees, species, by = c("speciesid" = "id"))
  data <- left_join(data, plots, by = "plotcode")
  data <- left_join(data, clims, by = "plotcode")

  data <-  data %>% group_by(plotcode) %>%
      mutate(BATOT_ha1=sum(ba_ha1))

  ## funBASUP <- function(df) {
  ##     data.frame(treecode = df$treecode,
  ##                BASUP = sapply(df$dbh1,
  ##                               function(x, dd) sum(dd$ba_ha1[dd$dbh1>x]),
  ##                               dd = df))
  ## }
  ## dfBASUP <- data %>% do(funBASUP(.))
  ## data <- left_join(data, dfBASUP, by = "treecode")


  if (remove_harv){
  # remove plot with harvesting
  data <- group_by(data, plotcode) %>%
           mutate(N_harv = sum(treestatus_th == 3)) %>%
           filter(N_harv <1)
  # remove French plots where management has been recorded, we don't have details
  # on individually harvested trees, the management column contains a 1 for those French
  # plots in which management has been recorded
  data$management[is.na(data$management)] <- 0
  data <- filter(data, management ==0)
  # remove plots with harvesting based on plot code
  harvs <- read.csv(file.path(data_path, "FunDivEUROPE_plot_management.csv"),
                    stringsAsFactors=FALSE)
  plots_with_harv <- harvs$plotcode[harvs$management2 >0 & !is.na(harvs$management2)]
  data <- filter(data, ! plotcode %in% plots_with_harv)
  }
  # remove plots in the Canary Islands
  data <- filter(data, latitude >30)
  # calculate the number per hectare from the weight
  data$n_ha1 <-  NA
  data$n_ha1[data$country %in% c('ES','FI','SW','WA') &
             data$treestatus_th>1 &
             !is.na(data$weight1)] <- 1/((data$weight1[data$country %in%
                                                       c('ES','FI','SW','WA') &
                                                       data$treestatus_th>1 &
                                                       !is.na(data$weight1)]*
                                          data$weight1[data$country %in%
                                                       c('ES','FI','SW','WA') &
                                                       data$treestatus_th>1 &
                                                       !is.na(data$weight1)]*
                                          3.14159265)/10000)
  data$n_ha2 <-  NA
  data$n_ha2[data$country %in% c('ES','FI','SW','WA') &
             data$treestatus_th<3 &
             !is.na(data$weight2)] <- 1/((data$weight2[data$country %in%
                                                       c('ES','FI','SW','WA') &
                                                       data$treestatus_th<3 &
                                                       !is.na(data$weight2)]*
                                          data$weight2[data$country %in%
                                                       c('ES','FI','SW','WA') &
                                                       data$treestatus_th<3 &
                                                       !is.na(data$weight2)]*
                                          3.14159265)/10000)
  data$n_ha1[data$country =='DE'] <- data$weight1[data$country =='DE']
  data$n_ha2[data$country =='DE'] <- data$weight2[data$country =='DE']
  # weight is missing in France recompute it base on ba and ba_ha

  data$n_ha1[data$country =='FG' &
               data$treestatus_th>1]  <- data$ba_ha1[data$country =='FG' &
                                              data$treestatus_th>1]/data$ba1[data$country =='FG' &
                                                                                data$treestatus_th>1]
  data$n_ha2[data$country =='FG' &
               data$treestatus_th<3]  <- data$ba_ha2[data$country =='FG' &
                                                        data$treestatus_th<3]/data$ba2[data$country =='FG' &
                                                        data$treestatus_th<3]
  #compute competitors basal area
  data <- data %>% mutate(BATOTcomp = BATOT_ha1 - ba_ha1)
  # CHANGE FG to FR
  data <- data %>% mutate(country=replace(country, country=='FG', 'FR')) %>%
     as.data.frame()

  data <- filter(data, !is.na(BATOT_ha1) & !is.na(surveydate2)& !is.na(surveydate1))
  data$treecode2 <- paste0("code_", seq_len(length.out = nrow(data)))
  ## remove unused variables
  if (remove_harv){
  data <- select(data, -c(treecode, speciesid, height1, height2, ba1,
                          ba_ha1, ba2, bachange_ha_yr, dbh1_mod,
                          ba1_mod, bachange_ha_yr_mod, weight1, weight2,
                          biome, surveydate1, end_year, tile,
                          bio1, BATOT_ha1, N_harv, n_ha1))
  }else{
  data <- select(data, -c(treecode, speciesid, height1, height2, ba1,
                          ba_ha1, ba2, bachange_ha_yr, dbh1_mod,
                          ba1_mod, bachange_ha_yr_mod, weight1, weight2,
                          biome, surveydate1, end_year, tile,
                          bio1, BATOT_ha1, n_ha1))
  }
    return(data)
}

LatLongSp <- function(spsel, df){
df %>% dplyr::filter(sp == spsel) %>%  select(plotcode, longitude, latitude)
}

FUNDIV_data_for_sdm <- function(df,sps){
ll <- lapply(sps$sp, LatLongSp, df = df)
names(ll) <- sps$sp
saveRDS(ll, file.path("output", "sdm_data.rds"))
}


read_FUNDIV_rege_data <- function(format = "wide", data_path = "data") {
  if(!format %in% c("wide", "long")) stop("format must be 'wide' or 'long'")
  #function to read and select data
  require(dplyr)
  rege <- read.csv(file.path(data_path, "FunDiv_regeneration_4sizeclasses.csv"),
                    stringsAsFactors=FALSE)
  # replace NA per 0. This is right Sophia ? Sophie: yes I assume so.
  rege$count_1[is.na(rege$count_1)] <- 0
  rege$count_2[is.na(rege$count_2)] <- 0
  rege$n_ha_1[is.na(rege$n_ha_1)] <- 0
  rege$n_ha_2[is.na(rege$n_ha_2)] <- 0
  # combine Betula pendula and pubescens with the Betula genus (ids 46 and 47 with 48)
  rege$speciesid[rege$speciesid %in% c(46,47)] <- 48
  species <- read.csv(file.path(data_path, "FunDiv_species_Nadja.csv"),
                      stringsAsFactors=FALSE)
  # remove characters at the end of the species name
  species$species[species$id ==277] <- "pubescens"
  species <- species %>% mutate(sp = paste(genus, species)) %>%
      select(c(id, sp))
  species[species$id == 48, "sp"] <- "Betula"
  species[species$id == 277, "sp"] <- "Quercus pubescens"
  rege <- left_join(rege, species, by = c("speciesid" = "id"))
  rege <-  filter(rege, sp != " ")
  rege <- as.data.frame(group_by(rege, plotcode, sp, sizecategory) %>%
                              summarize(count_1 = sum(count_1, na.rm = TRUE),
                                        n_ha_1 = sum(n_ha_1, na.rm = TRUE),
                                        count_2 = sum(count_2, na.rm = TRUE),
                                        n_ha_2 = sum(n_ha_2, na.rm = TRUE)) %>%
                        ungroup())

  # spread data by size class
  require(tidyr)
  rege_c1<- rege %>% mutate(obs_id = paste(plotcode, sp),
                          count_1_cat = sizecategory) %>%
            select(c(plotcode, sp, obs_id, count_1_cat, count_1)) %>%
      spread(count_1_cat, count_1, fill = 0, sep = "_")
  rege_c2<- rege %>% mutate(obs_id = paste(plotcode, sp),
                          count_2_cat = sizecategory) %>%
                          ungroup() %>%
            select(c(obs_id, count_2_cat, count_2)) %>%
      spread(count_2_cat, count_2, fill = 0, sep = "_")
  rege_nh1<- rege %>% mutate(obs_id = paste(plotcode, sp),
                          nha_1_cat = sizecategory) %>%
                          ungroup() %>%
            select(c(obs_id, nha_1_cat, n_ha_1)) %>%
      spread(nha_1_cat, n_ha_1, fill = 0, sep = "_")
  rege_nh2<- rege %>% mutate(obs_id = paste(plotcode, sp),
                          nha_2_cat = sizecategory) %>%
                          ungroup() %>%
            select(c(obs_id, nha_2_cat, n_ha_2)) %>%
      spread(nha_2_cat, n_ha_2, fill = 0, sep = "_")

  rege <- left_join(rege_c1, rege_c2, by = "obs_id")
  rege <- left_join(rege, rege_nh1, by = "obs_id")
  rege <- left_join(rege, rege_nh2, by = "obs_id")

  plots <- read.csv(file.path(data_path, "FunDiv_plots_Nadja.csv"),
                    stringsAsFactors=FALSE)
  clims <- read.csv(file.path(data_path,
                              "FunDiv_plots_climate_sapropos.csv"),
                    stringsAsFactors=FALSE)
  clims <- clims %>% select(-longitude, -latitude, -country )

  # merging data
  data <- left_join(rege, plots, by = "plotcode")
  data <- left_join(data, clims, by = "plotcode")

  # remove plots in the Canary Islands
  data <- filter(data, latitude >30)

  # remove plots with evidence of harvesting
  adult_trees <- read.csv(file.path(data_path, "FunDiv_trees_Nadja.csv"),
                    stringsAsFactors=FALSE)
  adult_trees <- group_by(adult_trees, plotcode) %>%
                      mutate(N_harv = sum(treestatus_th == 3)) %>%
                      filter(N_harv <1) %>%
                      select(plotcode)
  plot_sel <- unique(adult_trees$plotcode)
  data <- dplyr::filter(data, plotcode %in% plot_sel)

  if(format == "long"){
  ## get data in long
  data <-  gather(data, vart_t_sizecat, var, count_1_cat_1:nha_2_cat_4) %>%
      separate(vart_t_sizecat, c("vart", "census","cat", "sizecat"), "_") %>%
      spread(vart, var) %>% select(-cat)
  }
return(data)
}

combine_adult_and_sapling_data <- function(adult_data, sapling_data, data_path="data"){
  # exclude Wallonia and France as we don't have regeneration data from them
  adult_data <- filter(adult_data, country %in% c('DE','ES','SW','FI'))

  adult_plots <- group_by(adult_data, plotcode, sp) %>%
                    summarize(adult_ba_ha1 = sum(ba_ha1, na.rm = TRUE),
                              adult_ba_ha2 = sum(ba_ha2, na.rm = TRUE))

  ingrowth_plots <- group_by(adult_data, plotcode, sp) %>%
                        filter(treestatus_th == 1) %>%
                        summarize(ingrowth_n_ha2 = sum(n_ha2, na.rm = TRUE))

  adult_plots <- full_join(adult_plots,
                           ingrowth_plots,
                              by = c("plotcode","sp"))

  sapling_data <- sapling_data %>% select(plotcode:nha_2_cat_4)
  adult_saplings <- full_join(adult_plots,
                              sapling_data,
                              by = c("plotcode","sp"))
  # replace the NAs with 0, I don't like using the indexes it was the simplest solution
  adult_saplings[, 3:4][is.na(adult_saplings[ , 3:4])] <- 0
  adult_saplings[, 6:21][is.na(adult_saplings[ , 6:21])] <- 0
  plots <- read.csv(file.path(data_path, "FunDiv_plots_Nadja.csv"),
                    stringsAsFactors=FALSE)
  clims <- read.csv(file.path(data_path,
                              "FunDiv_plots_climate_sapropos.csv"),
                    stringsAsFactors=FALSE)
  # merging data
  data <- left_join(adult_saplings, plots, by = "plotcode")
  data <- left_join(data, clims, by = "plotcode")

  # CHANGE FG to FR
  data <- data %>% mutate(country=replace(country, country=='FG', 'FR')) %>%
     as.data.frame()

  return(data)
}

read_COMPADRE_sp <- function(data_path = "data"){
  COMPADRE_sp <- read.csv(file.path(data_path, "COMPADRE_species_list.csv"),
                          stringsAsFactors=FALSE)
  return(COMPADRE_sp)
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

com_sp <- function(data, COMPADRE_sp){
  sp_FD <- unique(data$sp)
  sp_CP <- trim(unique(COMPADRE_sp$SpeciesAccepted))
  return(sp_FD[sp_FD  %in% sp_CP])
}



Abundance_Bin_Breaks <- function(data, var = "bio1", breaks=clim_breaks, sp="Fagus sylvatica"){
  all <- data[data$treestatus_th %in% c(2, 4), ]
# abundance
  allplots <- tapply(all$treecode, list(all$plotcode,all$sp), length)
  allabun <- apply(allplots, 1, sum, na.rm=T)
  speciesabun <- allplots[,colnames(allplots) == sp] / allabun
  speciesabun[is.na(speciesabun)] <- 0
  plotclim <- tapply(all[[var]], all$plotcode, mean)
  plot.size_cut <- cut(plotclim, breaks=clim_breaks, include.lowest = TRUE, right = FALSE)
  abunperc <- aggregate(speciesabun,by=list(plot.size_cut), FUN = quantile, probs = c(0.5,0.6,0.7,0.8,0.9,0.95), simplify=T)
# basal area
  allplotsba <- tapply(all$ba_ha1, list(all$plotcode,all$sp), sum)
  allba <- apply(allplotsba, 1, sum, na.rm=T)
  speciesba <- allplotsba[,colnames(allplotsba) == sp] / allba
  speciesba[is.na(speciesba)] <- 0
  baperc <- aggregate(speciesba,by=list(plot.size_cut), FUN = quantile, probs = c(0.5,0.6,0.7,0.8,0.9,0.95), simplify=T)

  out <- data.frame(abunperc$x, baperc$x)
  colnames(out) <- c("abun50","abun60","abun70","abun80","abun90","abun95","ba50",
                     "ba60","ba70","ba80","ba90","ba95")
  return(out)
}


select_species <- function(data, nlim_tree  = 2000, nlim_plot = 500){
  # select the species with sufficient number of individuals
  data <- filter(data, sp != " " & treestatus_th != 1)
  nplot <- group_by(data, sp) %>% summarise(n_plot= n_distinct(plotcode))
  abund <- group_by(data, sp) %>% summarise(n_tree= n())
  res <- filter(arrange(left_join(nplot, abund, "sp"), desc(n_tree)),  n_tree> nlim_tree & n_plot > nlim_plot)
  #exclude exotic and castanea
  print("species selected")
  res <- filter(res, !sp %in% c("Pinus radiata", "Eucalyptus globulus",
                                "Castanea sativa",
                                "Pseudotsuga menziesii",
                                "Robinia pseudacacia",
                                "Quercus "))
  return(res)
}

species_climate_range<- function(data, var = "sgdd", sp="Fagus sylvatica", by.break = 50){
  # select climate range for the species
  # df <- data[data$sp == sp, ]
  # return(unique(quantile(df[[var]], probs = seq(0, 1, by = by.break), na.rm=TRUE)))
  df <- data[data$sp == sp & !is.na(data$sgdd), ]
  return(seq(min(df$sgdd),max(df$sgdd), by=by.break))
}

species_climate_range_sgdd <- function(data, sp="Fagus sylvatica", by.break = 50){
  df <- data[data$sp == sp & !is.na(data$sgdd), ]
  return(seq(min(df$sgdd),max(df$sgdd), by=by.break))
}

# function to compute the mid bin
compute_mid_bins <- function(x, br){
  ints  <- findInterval(x, br, rightmost.closed  = TRUE)
  (br[ints] + br[ints + 1]) / 2
}
compute_mid_bins <- function(x){
  str <- matrix(unlist(strsplit(x,",")), nrow=length(x), byrow=T)
  start <- as.numeric(substr(str[,1],2,nchar(str[,1])))
  end <- as.numeric(substr(str[,2],1,nchar(str[,2])-1))
  return(start+((end-start)/2))
}

growth_climate_bins_country <- function(data, var = "sgdd", br, sp="Fagus sylvatica"){
  # returns the growth by country in each climate bin for the given species
  require(tidyr)
  df <- data[data$treestatus_th %in% 2 & data$sp == sp, ]

  growth <- as.data.frame(tapply(((df$dbh2 - df$dbh1)/df$yearsbetweensurveys),
                                 list(df$country,cut(df[[var]], breaks=br,
                                include.lowest = TRUE, right = FALSE)),
                                median,
                                na.rm = TRUE))
  ncol <- ncol(growth)
  growth$country <- as.factor(rownames(growth))
  growth.long <- gather(growth, climate.bin, growth, 1:ncol, factor_key=FALSE)
  nobs <- as.data.frame(tapply(df$dbh2,
                               list(df$country,cut(df[[var]], breaks=br, include.lowest = TRUE, right = FALSE)),
                               length))
  nobs$country <- as.factor(rownames(nobs))
  nobs.long <- gather(nobs, climate.bin, growth.obs, 1:ncol, factor_key=FALSE)

  # growth_bin <- left_join(pcl.long, growth.long,  by = c("country", "climate.bin"))
  growth_bin <- left_join(growth.long, nobs.long,  by = c("country", "climate.bin"))
  growth_bin$country <- as.character(growth_bin$country)
  # remove bins with no growth observations (due to the gather function)
  growth_bin <- filter(growth_bin, !is.na(growth.obs))
  return(growth_bin)
}

survival_climate_bins_country <- function(data, var = "bio1", br, sp="Fagus sylvatica"){
  # returns the survival rate by country in each climate bin for the given species
  df <- data[data$treestatus_th %in% c(2, 4) & data$sp == sp, ]
  df$climate_cut <- cut(df[[var]], breaks=br, include.lowest = TRUE, right = FALSE)

  times <- sort(unique(df$yearsbetweensurveys))
  df$yearsbetweensurveys <- factor(df$yearsbetweensurveys)
  df_surv <- df[df$treestatus_th %in% 2 , ]
  df_dead <- df[df$treestatus_th %in% 4 , ]
  df_surv_dead <- df[df$treestatus_th %in% c(2,4), ]
  surv <- table(df_surv$climate_cut, df_surv$yearsbetweensurveys, df_surv$country)
  surv_dead <- table(df_surv_dead$climate_cut, df_surv_dead$yearsbetweensurveys, df_surv_dead$country)
  surv_rate <- surv/(surv_dead)
  surv_rate[surv_dead ==0] <-  NA
  for (i in seq_len(length(times))){
    surv_rate[ , i, ] <- surv_rate[ , i, ]^(1/times[i])
  }
  surv_bin <- as.data.frame(apply(surv_rate, MARGIN = c(1,3), FUN = median, na.rm = TRUE))
  surv_bin$climate.bin <- as.factor(rownames(surv_bin))
  surv.long <- gather(surv_bin, country, survival, 1:dim(surv_rate)[3], factor_key=FALSE)
  surv.long$climate.bin <- as.character(surv.long$climate.bin)

  nobs <- as.data.frame(apply(surv_dead, MARGIN = c(1,3), sum, na.rm = TRUE))
  nobs$climate.bin <- as.factor(rownames(nobs))
  nobs.long <- gather(nobs, country, surv.obs, 1:dim(surv_rate)[3], factor_key=FALSE)
  nobs.long$climate.bin <- as.character(nobs.long$climate.bin)

  survival_bin <- left_join(surv.long, nobs.long,  by = c("country", "climate.bin"))
  # remove bins with no survival observations (due to the gather function)
  survival_bin <- filter(survival_bin, surv.obs >0)
  return(survival_bin)
}

abundance_climate_bins_country <- function(data, var = "bio1", br, sp="Fagus sylvatica"){
  # returns the abunance and basal area (total and relative) of the species in each climate bin in each country
  # abundance
  df <- data[data$treestatus_th %in% c(2, 4), ]
  plot.abund <- aggregate(n_ha1 ~ country + plotcode + sp, data=df, FUN=sum)
  plot.clim <- as.data.frame(tapply(df[[var]], list(df$plotcode), mean))
  plot.clim$plotcode <- as.factor(rownames(plot.clim))
  colnames(plot.clim)[1] <- "mean.climate"

  plot.abund <- merge(plot.abund, plot.clim, by=c("plotcode"))
  plot.abund$climate.bin <- cut(plot.abund$mean.climate, breaks=br, include.lowest = TRUE, right = FALSE)

  abund.perc <- aggregate(n_ha1 ~ country + climate.bin, data=plot.abund[plot.abund$sp==sp,], FUN = quantile,
                          probs = c(0.5,0.6,0.7,0.8,0.9,0.95))

  abund.perc <- cbind(abund.perc$country, abund.perc$climate.bin, as.data.frame(abund.perc$n_ha1))
  colnames(abund.perc)[1:2] <- c('country','climate.bin')

  # basal area
  plot.ba <- aggregate(ba_ha1 ~ country + plotcode + sp, data=df, FUN=sum)
  plot.ba <- merge(plot.ba, plot.clim, by=c("plotcode"))
  plot.ba$climate.bin <- cut(plot.ba$mean.climate, breaks=br, include.lowest = TRUE, right = FALSE)

  ba.perc <- aggregate(ba_ha1 ~ country + climate.bin, data=plot.ba[plot.ba$sp==sp,], FUN = quantile,
                       probs = c(0.5,0.6,0.7,0.8,0.9,0.95))

  ba.perc <- cbind(ba.perc$country, ba.perc$climate.bin, as.data.frame(ba.perc$ba_ha1))
  colnames(ba.perc)[1:2] <- c('country','climate.bin')

  # relative abundance - based on number per hectare
  plot.abund.sp <- aggregate(n_ha1 ~ country + plotcode, data=df[df$sp==sp,], FUN=sum)
  plot.abund.all <- aggregate(n_ha1 ~ country + plotcode, data=df, FUN=sum)

  rel.abund <- merge(plot.abund.sp, plot.abund.all, by=c('plotcode','country'), all.y=T)
  colnames(rel.abund)[3:4] <- c('sp.abund','total.abund')
  rel.abund$sp.abund[is.na(rel.abund$sp.abund)] <- 0
  rel.abund$rel.abund <- rel.abund$sp.abund/rel.abund$total.abund

  plot.clim <- as.data.frame(tapply(df[[var]], list(df$plotcode), mean))
  plot.clim$plotcode <- as.factor(rownames(plot.clim))
  colnames(plot.clim)[1] <- "mean.climate"

  rel.abund <- merge(rel.abund, plot.clim, by=c("plotcode"))
  rel.abund$climate.bin <- cut(rel.abund$mean.climate, breaks=br, include.lowest = TRUE, right = FALSE)

  rel.abund.perc <- aggregate(rel.abund ~ country + climate.bin, data=rel.abund, FUN = quantile,
                              probs = c(0.5,0.6,0.7,0.8,0.9,0.95))

  rel.abund.perc <- cbind(rel.abund.perc$country,
                          rel.abund.perc$climate.bin,
                          as.data.frame(rel.abund.perc$rel.abund))
  colnames(rel.abund.perc)[1:2] <- c('country','climate.bin')

  # relative abundance - based on basal area
  plot.ba.sp <- aggregate(ba_ha1 ~ country + plotcode, data=df[df$sp==sp,], FUN=sum)
  plot.ba.all <- aggregate(ba_ha1 ~ country + plotcode, data=df, FUN=sum)

  rel.ba <- merge(plot.ba.sp, plot.ba.all, by=c('plotcode','country'), all.y=T)
  colnames(rel.ba)[3:4] <- c('sp.ba','total.ba')
  rel.ba$sp.ba[is.na(rel.ba$sp.ba)] <- 0
  rel.ba$rel.ba <- rel.ba$sp.ba/rel.ba$total.ba

  plot.clim <- as.data.frame(tapply(df[[var]], list(df$plotcode), mean))
  plot.clim$plotcode <- as.factor(rownames(plot.clim))
  colnames(plot.clim)[1] <- "mean.climate"

  rel.ba <- merge(rel.ba, plot.clim, by=c("plotcode"))
  rel.ba$climate.bin <- cut(rel.ba$mean.climate, breaks=br, include.lowest = TRUE, right = FALSE)

  rel.ba.perc <- aggregate(cbind(rel.ba, total.ba) ~ country + climate.bin, data=rel.ba,
                           FUN = quantile,
                           probs = c(0.5,0.6,0.7,0.8,0.9,0.95))

  rel.ba.perc <- cbind(rel.ba.perc$country, rel.ba.perc$climate.bin,
                       as.data.frame(rel.ba.perc$rel.ba),
                       as.data.frame(rel.ba.perc$total.ba))
  colnames(rel.ba.perc) <- c('country','climate.bin',
                             "relba50","relba60","relba70","relba80","relba90","relba95",
                             "totba50","totba60","totba70","totba80","totba90","totba95")
  rel.abund.perc$country <- as.character(rel.abund.perc$country);
  rel.abund.perc$climate.bin <- as.character(rel.abund.perc$climate.bin)
  abund.perc$country <- as.character(abund.perc$country);
  abund.perc$climate.bin <- as.character(abund.perc$climate.bin)
  ba.perc$country <- as.character(ba.perc$country);
  ba.perc$climate.bin <- as.character(ba.perc$climate.bin)
  rel.ba.perc$country <- as.character(rel.ba.perc$country);
  rel.ba.perc$climate.bin <- as.character(rel.ba.perc$climate.bin)

  out <- left_join(abund.perc, ba.perc,  by = c("country", "climate.bin"))
  out <- left_join(out, rel.abund.perc,  by = c("country", "climate.bin"))
  out <- left_join(out, rel.ba.perc,  by = c("country", "climate.bin"))

  colnames(out) <- c("country","climate.bin","abun50","abun60","abun70","abun80","abun90","abun95",
                     "ba50","ba60","ba70","ba80","ba90","ba95",
                     "relabund50","relabund60","relabund70","relabund80","relabund90","relabund95",
                     "relba50","relba60","relba70","relba80","relba90","relba95",
                     "totba50","totba60","totba70","totba80","totba90","totba95")
  return(out)
}

probabilitypresence_climate_bins_country <- function(data, var = "bio1", br,
                                                     sp="Fagus sylvatica"){
  data <- data[data$treestatus_th %in% c(2, 4), ]
  # probability of presence - I changed that to be based on number of plot per climate bin
  Abund_all<- as.data.frame(tapply(data$n_ha1,
                                list(data$plotcode, data$sp),
                                sum, na.rm = TRUE))
  PresAbs <- data.frame(plotcode = rownames(Abund_all),
                        PA = Abund_all[[sp]] >0)
  PresAbs$plotcode <- as.character(PresAbs$plotcode)
  PresAbs <- left_join(PresAbs, data[!duplicated(data$plotcode),
                                     c("plotcode", "country", var)],
                       by = "plotcode")
  Abund <- data.frame(plotcode = rownames(Abund_all),
                      Abundance = Abund_all[[sp]])
  Abund$plotcode <- as.character(Abund$plotcode)
  Abund <- left_join(Abund,
                     data[!duplicated(data$plotcode),
                          c("plotcode", "country", var)],
                     by = "plotcode")
  Abund$Abundance[is.na(Abund$Abundance)] <-  0
  PresAbs$PA[is.na(PresAbs$PA)] <-  FALSE
  PA2 <- rep(0, length.out = length(PresAbs$PA))
  PA2[PresAbs$PA] <-  1
  PresAbs$PA <- PA2
  PresAbs$count <-  1
  Pres <- tapply(PresAbs$PA,
                 list(PresAbs$country,cut(PresAbs[[var]],
                                          breaks=br, include.lowest = TRUE, right = FALSE)),
                 sum)
  Tot <- tapply(PresAbs$count,
                 list(PresAbs$country,cut(PresAbs[[var]],
                                          breaks=br, include.lowest = TRUE, right = FALSE)),
                 sum)
  ProbP <- as.data.frame(Pres/Tot)
  ProbP$country <- as.factor(rownames(ProbP))
  nc <- ncol(ProbP)-1
  ProbP.long <- gather(ProbP, climate.bin, ProbPres, 1:nc, factor_key=FALSE)

  abund <- as.data.frame(tapply(Abund$Abundance,
                                list(Abund$country,cut(Abund[[var]], breaks=br,
                                include.lowest = TRUE, right = FALSE)),
                                mean, na.rm = TRUE))
  abund$country <- as.factor(rownames(abund))
  nc <- ncol(abund)-1
  abund.long <- gather(abund, climate.bin, abundance, 1:nc, factor_key=FALSE)

  # remove bins with no observations (due to the gather function)
  abund.long$ProbPres <- ProbP.long$ProbPres
  abund.long <- filter(abund.long, !is.na(abundance))
  return(abund.long)
}

climate_bin_population_metrics <- function(data, climate.var = "sgdd",
                                           sp="Pinus sylvestris")
{
  # calculate the population metrics for the given species and climate variable
  # select climate range for the species
  br <- species_climate_range(data, var = climate.var, sp=sp)
  growth <- growth_climate_bins_country(data, var = climate.var, br=br, sp=sp)
  surv <- survival_climate_bins_country(data, var = climate.var, br=br, sp=sp)

  t1 <- merge(growth, surv, by=c('country','climate.bin'))

  t1$passing.time <- 200/t1$growth
  t1$life.expectancy <- 1/(1-t1$survival)
  t1$life.expectancy[is.infinite(t1$life.expectancy)] <- NA
  t1$lifetime.growth <- t1$growth*t1$life.expectancy

  abund <- abundance_climate_bins_country(data, var = climate.var, br=br, sp=sp)
  pp <- probabilitypresence_climate_bins_country(data, var = climate.var,
                                                 br=br, sp=sp)
  t1 <- merge(t1, abund, by=c('country','climate.bin'), all.x = T)
  t1 <- merge(t1, pp, by=c('country','climate.bin'), all.x = T)

  t1$climate.bin.mid <- compute_mid_bins(t1$climate.bin)
  t1$mean.climate <- t1$climate.bin.mid
  t1$sp <- sp
  t1$climate.var <- climate.var

  return(t1)
}

climate_bin_population_metrics_all_species_df <- function(data, sps,
                                                       climate.var = "sgdd"){
  # calculate the population metrics for the selected species and climate variable
  out.df <- NULL
  for(i in 1:nrow(sps)){
    print(sps$sp[i])
    df <- climate_bin_population_metrics(data, climate.var, sp=sps$sp[i])
    if(is.null(out.df)){
      out.df <- df
    }else{
      out.df <- rbind(out.df, df)
    }
  }
  return(out.df)
}

sapling_metrics <- function(data, selsp="Pinus sylvestris"){
  # returns for each plot and species: ratio of size category 2 saplings to adults
  # ratio of n_ha in the second size class to the first
  # nha of ingrowth trees
  # ratio of nha of ingrowth to adult trees
  require(tidyr)
  df <- data[data$sp==selsp,]

  if(nrow(df[df$nha_2_cat_2>0,]) >0) {

    # add a minimum adult baha to all plots to account for plots without any adults
    min_adult_baha2 <- aggregate(adult_ba_ha2 ~ country, data=df[df$adult_ba_ha2>0,], FUN=min)
    colnames(min_adult_baha2)[2] <- "min_adult_baha2"
    df <- merge(df, min_adult_baha2, by=c("country"))
    df$adult_ba_ha2 <- df$adult_ba_ha2 + df$min_adult_baha2

    # Metric 1: ratio of saplings (size cat 2) to adults in the plots
    df$saplings2_adults <- df$nha_2_cat_2/df$adult_ba_ha2

    # add a minimum nha to all plots to account for plots without any category 1 saplings
    min_naha2 <- aggregate(nha_2_cat_1 ~ country, data=df[df$nha_2_cat_1>0,], FUN=min)
    colnames(min_naha2)[2] <- "min_sapling1_nha2"
    df <- merge(df, min_naha2, by=c("country"))
    df$nha_2_cat_1 <- df$nha_2_cat_1 + df$min_sapling1_nha2

    # Metric 2: ratio of saplings (size cat 2) to saplings (size cat 1)
    df$saplings2_saplings1 <- df$nha_2_cat_2/df$nha_2_cat_1

    # Metric 3: baha of ingrowth trees
    df$ingrowth_n_ha2[is.na(df$ingrowth_n_ha2)] <- 0

    # Metric 4: baha of ingrowth trees / baha in first survey of adult trees
    df$ingrowth_adult <- df$ingrowth_n_ha2/df$adult_ba_ha2

    df <- df %>% dplyr::select(-count_1_cat_1, -count_1_cat_2, -count_1_cat_3, -count_1_cat_4,
                               -count_2_cat_1, -count_2_cat_2, -count_2_cat_3, -count_2_cat_4,
                               -nha_1_cat_1, -nha_1_cat_2, -nha_1_cat_3, -nha_1_cat_4,
                               -nha_2_cat_1, -nha_2_cat_2, -nha_2_cat_3, -nha_2_cat_4,
                               -min_adult_baha2, -min_sapling1_nha2)
    return(df)
  }else{return(NULL)}
}

sapling_metrics_all_species_df <- function(adult_data, sapling_data, sps){
  data <- combine_adult_and_sapling_data(adult_data=adult_data, sapling_data=sapling_data)
  out.df <- NULL
  for(i in 1:nrow(sps)){
    print(sps$sp[i])
    df <- sapling_metrics(data, selsp=sps$sp[i])
    if(!is.null(df)){
      if(is.null(out.df)){
        out.df <- df
      }else{
        out.df <- rbind(out.df, df)
      }
    }
  }
  return(out.df)
}


compute_N_ingrowth <- function(data){
  data <- filter(data, country %in% c('FR','ES','SW','FI'))
  vars_plot <- c("plotcode", "speciesid", "country", "sp", "cluster", "longitude",
                 "latitude", "yearsbetweensurveys", "biome", "management",
                 "surveydate1", "surveydate2", "start_year", "end_year",
                 "tile", "mat", "sgdd", "spring_frosts", "map", "bio1", "map.wc",
                 "wai.wc", "mat.wc", "wai", "p_pet_yr", "sws", "p_pet_summer",
                 "spei_min", "spei_mean", "N_harv", "n_ha1", "n_ha2", "BATOTcomp")
  colNums <- match(vars_plot,names(data))
  data_plots <- data %>% filter(!duplicated(paste(plotcode, sp))) %>% select(colNums)
  table(data$treestatus_th[data$country == "FR" & data$dbh1 <100 &
                     data$treestatus_th %in% c(1,2) & !is.na(data$dbh1)])
  ingrowth_plots <- group_by(data, plotcode, sp) %>%
                        filter(treestatus_th %in% c(1,2)) %>%
                        summarize(tot_n_ha2 = sum(n_ha2, na.rm = TRUE),
                                  ingrowth_n_ha2 = sum(n_ha2[treestatus_th == 1 & dbh2<130], na.rm = TRUE))

    ## TODO improve selection of recruited tree so far focus only on small tree but this is very weak (FOR France is their a problem??)
 df <- left_join(data_plots, ingrowth_plots, by = c("plotcode", "sp"))
    return(df)
}
