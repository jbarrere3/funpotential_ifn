#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_analysis.R  
#' @description R script containing all functions relative to data
#               analysis
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Section 1 - Bayesian mortality model ####
#' @description Group of functions used to fit the mortality model in jags
#' @authors Julien BARRERE
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Format the French data as input for the mortality model for a specific species
#' @param FUNDIV_tree FUNDIV tree table
#' @param FUNDIV_climate climate index for each FUNDIV plots
#' @param disturbance_per_plot dataframe containging the type, area and year of each 
#'                                       disturbance intercepting a FUNDIV plot buffer
#' @param NFI_plot_remeasure Table containing NFI plot data for remeasured plots
#' @param species.in character specifying the species for which to run the model (ex: "Abies alba")
format_FrenchNFI_mortality <- function(FUNDIV_tree, FUNDIV_climate, disturbance_per_plot, 
                                       NFI_plot_remeasure, species.in){
  
  print(species.in)
  
  #%%%%%%%%%%%%%%%%
  # Step 1: Format the data for the model
  #%%%%%%%%%%%%%%%%
  
  print("---Format dataset")
  # Format data
  data <- FUNDIV_tree %>%
    filter(country == "FR" & treestatus_th %in% c(2, 3, 4, 5)) %>%
    mutate(h = case_when(treestatus_th == 3 ~ 1, T ~ 0), 
           d = case_when(treestatus_th %in% c(4, 5) ~ 1, T ~ 0), 
           a = case_when(treestatus_th == 2 ~ 1, T ~ 0)) %>%
    rename(comp = BATOTcomp, dbh = dbh1, years = yearsbetweensurveys, id = treecode2) %>%
    dplyr::select(id, plotcode, sp, h, d, a, years, dbh, comp) %>%
    merge(subset(FUNDIV_climate, select = c("plotcode", "wai", "sgdd")), 
          by = "plotcode", all.x = T, all.y = F) %>%
    merge((disturbance_per_plot %>% 
             group_by(plotcode) %>% 
             summarise(DA = sum(coverage, na.rm = T), 
                       DS = sum(severity*coverage, na.rm = T)/sum(coverage, na.rm = T))), 
          by = "plotcode", all.x = T, all.y = F) %>%
    filter(!is.na(sgdd)) %>%
    filter(!is.na(DA))  %>%
    dplyr::select(-id)
  
  # Normalize variables between 0 and 1
  var.normalize <- c("dbh", "comp", "wai", "sgdd", "DA", "DS")
  data <- cbind(data[, !(colnames(data) %in% var.normalize)], 
                sapply(data[, var.normalize], 
                       function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm=T)))) 
  
  
  
  #%%%%%%%%%%%%%%%%
  # Step 2: Compute prior values for harvest conditional probabilities
  #%%%%%%%%%%%%%%%%
  
  print("---Compute prior values for harvest conditional probabilities")
  
  # Initialize list of parameters
  param <- list()
  
  # Format dataset
  data_priorharvest <- data %>%
    dplyr::select(plotcode, dbh, a, h) %>%
    merge((NFI_plot_remeasure %>%
             mutate(plotcode = paste0(idp, "_FG"), 
                    D = case_when(nincid5 > 0 ~ 1, nincid5 == 0 ~ 0), 
                    DS = case_when(incid5 == 0 ~ 0, 
                                   incid5 == 1 ~ 0.125, 
                                   incid5 == 2 ~ 0.375, 
                                   incid5 == 3 ~ 0.625, 
                                   incid5 == 4 ~ 0.875)) %>%
             dplyr::select(plotcode, D, DS)), 
          by = "plotcode", all.x = T, all.y = F) %>%
    filter(!is.na(D)) %>%
    filter(!(D == 1 & is.na(DS))) %>%
    filter(!(a == 1 & D == 1)) %>% # Remove alive trees in disturbed plots
    dplyr::select(-a)
  
  # Fit model for undisturbed plots
  mod.undist <- glm(h ~ dbh, 
                    data = subset(data_priorharvest, D == 0), 
                    family = binomial(link = "logit"))
  # Get parameters e0 and e1 (cf. doc Multinomial harvest model)
  param$e0 <- as.numeric(mod.undist$coefficients[1])
  param$e1 <- as.numeric(mod.undist$coefficients[2])
  
  # Fit model for disturbed plots
  mod.dist <- glm(h ~ dbh*DS, 
                  data = subset(data_priorharvest, D == 1), 
                  family = binomial(link = "logit"))
  
  # Get parameters d0, d1, d2 and d3 (cf. doc Multinomial harvest model)
  param$d0 <- as.numeric(mod.dist$coefficients[1])
  param$d1 <- as.numeric(mod.dist$coefficients[2])
  param$d2 <- as.numeric(mod.dist$coefficients[3])
  param$d3 <- as.numeric(mod.dist$coefficients[4])
  
  # Restrict dataset to studied species
  data <- data %>%
    filter(sp == species.in) %>%
    dplyr::select(-sp) %>%
    arrange(plotcode)
  
  # separate data at plot and tree level
  data_plots <- data %>%
    dplyr::select(plotcode, years, wai, sgdd, DA, DS) %>%
    distinct() %>%
    mutate(beg = NA_real_, end = NA_real_)
  for(i in 1:dim(data_plots)[1]){
    trees_i <- which(data$plotcode == data_plots$plotcode[i])
    data_plots$beg[i] <- trees_i[1]
    data_plots$end[i] <- trees_i[length(trees_i)]
  }
  data_trees <- subset(data, select = c("h", "d", "a", "dbh", "comp"))
  
  
  
  #%%%%%%%%%%%%%%%%
  # Step 3 Write and fit bayesian model
  #%%%%%%%%%%%%%%%%
  
  # Inputs
  data_jags <- list(
    # At the plot level
    Nplots = dim(data_plots)[1],
    wai = data_plots$wai,
    sgdd = data_plots$sgdd,
    DA = data_plots$DA,
    DS = data_plots$DS,
    beg = data_plots$beg,
    end = data_plots$end,
    # At the tree level
    h = data_trees$h,
    d = data_trees$d,
    a = data_trees$a,
    dbh = data_trees$dbh,
    comp = data_trees$comp,
    # priors
    d0 = param$d0,
    d1 = param$d1,
    d2 = param$d2,
    d3 = param$d3,
    e0 = param$e0,
    e1 = param$e1)
  
  return(data_jags)
  
}

#' Generate simulated data for mortality model
#' @param nDataPerPlots numeric: number of trees per plots
#' @param nPlots Number of plots
#' @param ref.data.in Reference dataset to generate realistic data
generate_data_jags <- function(nDataPerPlots, nPlots, ref.data.in){
  # Initialize output
  out <- list()
  # Generate prior values
  out$parameters <- list()
  out$parameters$a0 <- rnorm(nPlots, mean = runif(1, -5, 5), sd = 0.1)
  out$parameters$a1 <- rnorm(nPlots, mean = runif(1, -5, 5), sd = 0.1)
  out$parameters$b0 <- rnorm(nPlots*nDataPerPlots, mean = runif(1, -5, 5), sd = 0.1)
  out$parameters$b1 <- rnorm(nPlots*nDataPerPlots, mean = runif(1, -5, 5), sd = 0.1)
  out$parameters$b2 <- rnorm(nPlots*nDataPerPlots, mean = runif(1, -5, 5), sd = 0.1)
  out$parameters$b3 <- rnorm(nPlots*nDataPerPlots, mean = runif(1, -5, 5), sd = 0.1)
  out$parameters$b4 <- rnorm(nPlots*nDataPerPlots, mean = runif(1, -5, 5), sd = 0.1)
  out$parameters$c0 <- rnorm(nPlots*nDataPerPlots, mean = runif(1, -5, 5), sd = 0.1)
  out$parameters$c1 <- rnorm(nPlots*nDataPerPlots, mean = runif(1, -5, 5), sd = 0.1)
  out$parameters$c2 <- rnorm(nPlots*nDataPerPlots, mean = runif(1, -5, 5), sd = 0.1)
  
  fitbeta <-  function(x, n){
    x.in <- x; x.in[which(x == 1)] = 0.99; x.in[which(x == 0)] = 0.01
    fit.in <- fitdistr(x.in, "beta", list(shape1 = 100*mean(x.in), shape2 = 100*(1 - mean(x.in))))
    rbeta(n, as.numeric(fit.in$estimate[1]), as.numeric(fit.in$estimate[2]))
  } 
  
  # Generate data at the plot level
  data.plots.in <- data.frame(plot.id = paste0("plot_", c(1:nPlots)), 
                              sgdd = fitbeta(ref.data.in$sgdd, nPlots), 
                              wai = fitbeta(ref.data.in$wai, nPlots), 
                              DA = fitbeta(ref.data.in$DA, nPlots), 
                              DS = fitbeta(ref.data.in$DS, nPlots), 
                              D = NA_integer_, 
                              beg = c(0:(nPlots - 1))*nDataPerPlots + 1, 
                              end = c(1:nPlots)*nDataPerPlots)
  for(j in 1:nPlots){
    data.plots.in$D[j] <- rbinom(1, 1, exp(out$parameters$a0 + out$parameters$a1*data.plots.in$DA[j])/
                                   (1 + exp(out$parameters$a0 + out$parameters$a1*data.plots.in$DA[j])))
  }
  
  # Generate data at tree level
  data.tree.in <- data.frame(tree.id = paste0("tree_", c(1:(nPlots*nDataPerPlots))), 
                             plot.id = rep(paste0("plot_", c(1:nPlots)), array(nDataPerPlots, dim = nPlots)), 
                             dbh = fitbeta(ref.data.in$dbh, nPlots*nDataPerPlots), 
                             comp = fitbeta(ref.data.in$comp, nPlots*nDataPerPlots), 
                             h = NA_integer_, 
                             d = NA_integer_, 
                             a = NA_integer_, 
                             b0 = out$parameters$b0, 
                             b1 = out$parameters$b1, 
                             b2 = out$parameters$b2, 
                             b3 = out$parameters$b3, 
                             b4 = out$parameters$b4, 
                             c0 = out$parameters$c0, 
                             c1 = out$parameters$c1, 
                             c2 = out$parameters$c2) %>%
    merge(data.plots.in, by = "plot.id", all.x = T, all.y = F) %>%
    # Compute intermediate probabilities based on parameters and data
    mutate(pdD = c0 + c1*DS + c2*DS*dbh, 
           pdBM = b0 + b1*dbh + b2*comp + b3*sgdd + b4*wai, 
           phdD = ref.data.in$d0 + ref.data.in$d1*dbh + ref.data.in$d2*DS + ref.data.in$d3*dbh*DS, 
           phadBM = ref.data.in$e0 + ref.data.in$e1*dbh) %>%
    # Apply inverse logit function to constrain probabilities between 0 and 1
    mutate(pdD = exp(pdD)/(1 + exp(pdD)), 
           pdBM = exp(pdBM)/(1 + exp(pdBM)), 
           phdD = exp(phdD)/(1 + exp(phdD)), 
           phadBM = exp(phadBM)/(1 + exp(phadBM))) %>%
    mutate(pdDj = 1 - (1 - pdD)*(1 - pdBM)) %>%
    # Compute probability to be alive, dead or harvested
    mutate(ph = (1 - D)*phadBM + D*pdDj*phdD, 
           pd = (1 - D)*pdBM*(1 - phadBM) + D*pdDj*(1 - phdD)) %>%
    mutate(pa = 1 - (ph + pd))
  
  for(i in 1:(nPlots*nDataPerPlots)){
    data.tree.in[i, c("h", "d", "a")] <-
      as.integer(rmultinom(n = 1, size = 1, 
                           prob = as.numeric(data.tree.in[i, c("ph", "pd", "pa")])))}
  
  # Data output
  out$data <- list(
    # At the plot level
    Nplots = nPlots,
    wai = data.plots.in$wai,
    sgdd = data.plots.in$sgdd,
    DA = data.plots.in$DA,
    DS = data.plots.in$DS,
    beg = data.plots.in$beg,
    end = data.plots.in$end,
    # At the tree level
    h = data.tree.in$h,
    d = data.tree.in$d,
    a = data.tree.in$a,
    dbh = data.tree.in$dbh,
    comp = data.tree.in$comp,
    # priors
    d0 = ref.data.in$d0,
    d1 = ref.data.in$d1,
    d2 = ref.data.in$d2,
    d3 = ref.data.in$d3,
    e0 = ref.data.in$e0,
    e1 = ref.data.in$e1)
  
  return(out)
}

#' Fit the mortality model
#' @param data_jags.in List contianing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @return A rjags object
fit_mortality <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin){
  # Model
  model <- 
    "model{
  
  # Loop on all plots
  for(j in 1:Nplots){
    D[j] ~ dbern(pD[j])
    logit(pD[j]) = a0 + a1*DA[j]
    
    # Loop on trees
    for(i in beg[j]:end[j]){
      
      # Probability to observe harvested tree
      h[i] ~ dbern(ph[i])
      ph[i] = (1 - D[j])*phadBM[i] + D[j]*pdDj[i]*phdD[i]
      
      # Probability to observe dead tree
      d[i] ~ dbern(pd[i])
      pd[i] = (1 - D[j])*pdBM[i]*(1 - phadBM[i]) + D[j]*pdDj[i]*(1 - phdD[i])
      
      # Probability to observe alive tree
      a[i] ~ dbern(pa[i])
      pa[i] = (1 - D[j])*(1 - pdBM[i])*(1 - phadBM[i]) + D[j]*(1 - pdDj[i])
      
      # Probability to die in disturbed plots
      pdDj[i] = 1 - (1 - pdBM[i])*(1 - pdD[i])
      
      # Probability to die from background mortality
      logit(pdBM[i]) = b0 + b1*dbh[i] + b2*comp[i] + b3*sgdd[j] + b4*wai[j]
      
      # Probability to die from a disturbance
      logit(pdD[i]) = c0 + c1*DS[j] + c2*dbh[i]*DS[j]
      
      # Probability to be harvested knowing that the tree is dead from a disturbance
      logit(phdD[i]) = d0 + d1*dbh[i] + d2*DS[j] + d3*dbh[i]*DS[j]
      
      # Probability to be harvested in undisturbed plots
      logit(phadBM[i]) = e0 + e1*dbh[i]
    }
  }
  
  # Priors
  a0 ~ dnorm(0, 0.01)
  a1 ~ dnorm(0, 0.01)
  b0 ~ dnorm(0, 0.01)
  b1 ~ dnorm(0, 0.01)
  b2 ~ dnorm(0, 0.01)
  b3 ~ dnorm(0, 0.01)
  b4 ~ dnorm(0, 0.01)
  c0 ~ dnorm(0, 0.01)
  c1 ~ dnorm(0, 0.01)
  c2 ~ dnorm(0, 0.01)
  

  }"
  
  
  
  # Function to initialize priors
  initjags <- function(){
    return(list(a0 = runif(1, -10, 10),
                a1 = runif(1, -10, 10),
                b0 = runif(1, -10, 10),
                b1 = runif(1, -10, 10),
                b2 = runif(1, -10, 10),
                b3 = runif(1, -10, 10),
                b4 = runif(1, -10, 10),
                c0 = runif(1, -10, 10),
                c1 = runif(1, -10, 10),
                c2 = runif(1, -10, 10)))
  }
  
  # Run model
  out <- R2jags::jags(data = data_jags.in,
                      param = c("a0", "a1", "b0", "b1", "b2", "b3", 
                                "b4", "c0", "c1", "c2"),
                      inits = list(initjags(),
                                   initjags(),
                                   initjags()),
                      model.file = textConnection(model),
                      n.chains = n.chains,
                      n.iter = n.iter,
                      n.burnin = n.burn,
                      n.thin = n.thin,
                      DIC = TRUE, 
                      progress.bar = "text")
  return(out)
}



#' Format the French data as input for the mortality model for a specific species with true information on Dj
#' @param FUNDIV_tree FUNDIV tree table
#' @param FUNDIV_climate climate index for each FUNDIV plots
#' @param disturbance_per_plot dataframe containging the type, area and year of each 
#'                                       disturbance intercepting a FUNDIV plot buffer
#' @param NFI_plot_remeasure Table containing NFI plot data for remeasured plots
#' @param species.in character specifying the species for which to run the model (ex: "Abies alba")
format_FrenchNFI_mortality_2 <- function(FUNDIV_tree, FUNDIV_climate, disturbance_per_plot, 
                                         NFI_plot_remeasure, species.in){
  
  print(species.in)
  
  #%%%%%%%%%%%%%%%%
  # Step 1: Format the data for the model
  #%%%%%%%%%%%%%%%%
  
  print("---Format dataset")
  # Format data
  data <- FUNDIV_tree %>%
    filter(country == "FR" & treestatus_th %in% c(2, 3, 4, 5)) %>%
    mutate(h = case_when(treestatus_th == 3 ~ 1, T ~ 0), 
           d = case_when(treestatus_th %in% c(4, 5) ~ 1, T ~ 0), 
           a = case_when(treestatus_th == 2 ~ 1, T ~ 0)) %>%
    rename(comp = BATOTcomp, dbh = dbh1, years = yearsbetweensurveys, id = treecode2) %>%
    dplyr::select(id, plotcode, sp, h, d, a, years, dbh, comp) %>%
    merge(subset(FUNDIV_climate, select = c("plotcode", "wai", "sgdd")), 
          by = "plotcode", all.x = T, all.y = F) %>%
    merge((disturbance_per_plot %>% 
             group_by(plotcode) %>% 
             summarise(DA = sum(coverage, na.rm = T), 
                       DS = sum(severity*coverage, na.rm = T)/sum(coverage, na.rm = T))), 
          by = "plotcode", all.x = T, all.y = F) %>%
    filter(!is.na(sgdd)) %>%
    filter(!is.na(DA))  %>%
    dplyr::select(-id)
  
  # Normalize variables between 0 and 1
  var.normalize <- c("dbh", "comp", "wai", "sgdd", "DA", "DS")
  data <- cbind(data[, !(colnames(data) %in% var.normalize)], 
                sapply(data[, var.normalize], 
                       function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm=T)))) 
  
  
  
  #%%%%%%%%%%%%%%%%
  # Step 2: Compute prior values for harvest conditional probabilities
  #%%%%%%%%%%%%%%%%
  
  print("---Compute prior values for harvest conditional probabilities")
  
  # Initialize list of parameters
  param <- list()
  
  # Format dataset
  data_priorharvest <- data %>%
    dplyr::select(plotcode, dbh, a, h) %>%
    merge((NFI_plot_remeasure %>%
             mutate(plotcode = paste0(idp, "_FG"), 
                    D = case_when(nincid5 > 0 ~ 1, nincid5 == 0 ~ 0), 
                    DS = case_when(incid5 == 0 ~ 0, 
                                   incid5 == 1 ~ 0.125, 
                                   incid5 == 2 ~ 0.375, 
                                   incid5 == 3 ~ 0.625, 
                                   incid5 == 4 ~ 0.875)) %>%
             dplyr::select(plotcode, D, DS)), 
          by = "plotcode", all.x = T, all.y = F) %>%
    filter(!is.na(D)) %>%
    filter(!(D == 1 & is.na(DS))) %>%
    filter(!(a == 1 & D == 1)) %>% # Remove alive trees in disturbed plots
    dplyr::select(-a)
  
  # Fit model for undisturbed plots
  mod.undist <- glm(h ~ dbh, 
                    data = subset(data_priorharvest, D == 0), 
                    family = binomial(link = "logit"))
  # Get parameters e0 and e1 (cf. doc Multinomial harvest model)
  param$e0 <- as.numeric(mod.undist$coefficients[1])
  param$e1 <- as.numeric(mod.undist$coefficients[2])
  
  # Fit model for disturbed plots
  mod.dist <- glm(h ~ dbh*DS, 
                  data = subset(data_priorharvest, D == 1), 
                  family = binomial(link = "logit"))
  
  # Get parameters d0, d1, d2 and d3 (cf. doc Multinomial harvest model)
  param$d0 <- as.numeric(mod.dist$coefficients[1])
  param$d1 <- as.numeric(mod.dist$coefficients[2])
  param$d2 <- as.numeric(mod.dist$coefficients[3])
  param$d3 <- as.numeric(mod.dist$coefficients[4])
  
  # Restrict dataset to studied species
  data <- data %>%
    filter(sp == species.in) %>%
    dplyr::select(-sp) %>%
    arrange(plotcode)
  
  ## separate data at plot and tree level
  # Format plot level data
  data_plots <- data %>%
    dplyr::select(plotcode, years, wai, sgdd, DA, DS) %>%
    merge((data_priorharvest %>% dplyr::select(plotcode, D)), 
          by = "plotcode", all.x = T, all.y = F) %>%
    filter(!is.na(D)) %>%
    distinct() %>%
    mutate(beg = NA_real_, end = NA_real_)
  # Format tree level data
  data_trees <- data %>%
    filter(plotcode %in% data_plots$plotcode) %>%
    dplyr::select(h, d, a, dbh, comp)
  # Add at which lines of data_tree the plot starts and ends 
  for(i in 1:dim(data_plots)[1]){
    trees_i <- which(data$plotcode == data_plots$plotcode[i])
    data_plots$beg[i] <- trees_i[1]
    data_plots$end[i] <- trees_i[length(trees_i)]
  }
  
  
  #%%%%%%%%%%%%%%%%
  # Step 3 Write and fit bayesian model
  #%%%%%%%%%%%%%%%%
  
  # Inputs
  data_jags <- list(
    # At the plot level
    Nplots = dim(data_plots)[1],
    wai = data_plots$wai,
    sgdd = data_plots$sgdd,
    D = data_plots$D,
    DS = data_plots$DS,
    beg = data_plots$beg,
    end = data_plots$end,
    # At the tree level
    h = data_trees$h,
    d = data_trees$d,
    a = data_trees$a,
    dbh = data_trees$dbh,
    comp = data_trees$comp,
    # priors
    d0 = param$d0,
    d1 = param$d1,
    d2 = param$d2,
    d3 = param$d3,
    e0 = param$e0,
    e1 = param$e1)
  
  return(data_jags)
  
}

#' Generate simulated data for mortality model with true data of Dj
#' @param nDataPerPlots numeric: number of trees per plots
#' @param nPlots Number of plots
#' @param ref.data.in Reference dataset to generate realistic data
generate_data_jags_2 <- function(nDataPerPlots, nPlots, ref.data.in){
  # Initialize output
  out <- list()
  # Generate prior values
  out$parameters <- list()
  out$parameters$b0 <- rnorm(nPlots*nDataPerPlots, mean = runif(1, -5, 5), sd = 0.1)
  out$parameters$b1 <- rnorm(nPlots*nDataPerPlots, mean = runif(1, -5, 5), sd = 0.1)
  out$parameters$b2 <- rnorm(nPlots*nDataPerPlots, mean = runif(1, -5, 5), sd = 0.1)
  out$parameters$b3 <- rnorm(nPlots*nDataPerPlots, mean = runif(1, -5, 5), sd = 0.1)
  out$parameters$b4 <- rnorm(nPlots*nDataPerPlots, mean = runif(1, -5, 5), sd = 0.1)
  out$parameters$c0 <- rnorm(nPlots*nDataPerPlots, mean = runif(1, -5, 5), sd = 0.1)
  out$parameters$c1 <- rnorm(nPlots*nDataPerPlots, mean = runif(1, -5, 5), sd = 0.1)
  out$parameters$c2 <- rnorm(nPlots*nDataPerPlots, mean = runif(1, -5, 5), sd = 0.1)
  
  fitbeta <-  function(x, n){
    x.in <- x; x.in[which(x == 1)] = 0.99; x.in[which(x == 0)] = 0.01
    fit.in <- fitdistr(x.in, "beta", list(shape1 = 100*mean(x.in), shape2 = 100*(1 - mean(x.in))))
    rbeta(n, as.numeric(fit.in$estimate[1]), as.numeric(fit.in$estimate[2]))
  } 
  
  # Generate data at the plot level
  data.plots.in <- data.frame(plot.id = paste0("plot_", c(1:nPlots)), 
                              sgdd = fitbeta(ref.data.in$sgdd, nPlots), 
                              wai = fitbeta(ref.data.in$wai, nPlots), 
                              DS = fitbeta(ref.data.in$DS, nPlots), 
                              D = rbinom(nPlots, 1, mean(ref.data.in$D)), 
                              beg = c(0:(nPlots - 1))*nDataPerPlots + 1, 
                              end = c(1:nPlots)*nDataPerPlots)
  
  # Generate data at tree level
  data.tree.in <- data.frame(tree.id = paste0("tree_", c(1:(nPlots*nDataPerPlots))), 
                             plot.id = rep(paste0("plot_", c(1:nPlots)), array(nDataPerPlots, dim = nPlots)), 
                             dbh = fitbeta(ref.data.in$dbh, nPlots*nDataPerPlots), 
                             comp = fitbeta(ref.data.in$comp, nPlots*nDataPerPlots), 
                             h = NA_integer_, 
                             d = NA_integer_, 
                             a = NA_integer_, 
                             b0 = out$parameters$b0, 
                             b1 = out$parameters$b1, 
                             b2 = out$parameters$b2, 
                             b3 = out$parameters$b3, 
                             b4 = out$parameters$b4, 
                             c0 = out$parameters$c0, 
                             c1 = out$parameters$c1, 
                             c2 = out$parameters$c2) %>%
    merge(data.plots.in, by = "plot.id", all.x = T, all.y = F) %>%
    # Compute intermediate probabilities based on parameters and data
    mutate(pdD = c0 + c1*DS + c2*DS*dbh, 
           pdBM = b0 + b1*dbh + b2*comp + b3*sgdd + b4*wai, 
           phdD = ref.data.in$d0 + ref.data.in$d1*dbh + ref.data.in$d2*DS + ref.data.in$d3*dbh*DS, 
           phadBM = ref.data.in$e0 + ref.data.in$e1*dbh) %>%
    # Apply inverse logit function to constrain probabilities between 0 and 1
    mutate(pdD = exp(pdD)/(1 + exp(pdD)), 
           pdBM = exp(pdBM)/(1 + exp(pdBM)), 
           phdD = exp(phdD)/(1 + exp(phdD)), 
           phadBM = exp(phadBM)/(1 + exp(phadBM))) %>%
    mutate(pdDj = 1 - (1 - pdD)*(1 - pdBM)) %>%
    # Compute probability to be alive, dead or harvested
    mutate(ph = (1 - D)*phadBM + D*pdDj*phdD, 
           pd = (1 - D)*pdBM*(1 - phadBM) + D*pdDj*(1 - phdD)) %>%
    mutate(pa = 1 - (ph + pd))
  
  for(i in 1:(nPlots*nDataPerPlots)){
    data.tree.in[i, c("h", "d", "a")] <-
      as.integer(rmultinom(n = 1, size = 1, 
                           prob = as.numeric(data.tree.in[i, c("ph", "pd", "pa")])))}
  
  # Data output
  out$data <- list(
    # At the plot level
    Nplots = nPlots,
    wai = data.plots.in$wai,
    sgdd = data.plots.in$sgdd,
    D = data.plots.in$D,
    DS = data.plots.in$DS,
    beg = data.plots.in$beg,
    end = data.plots.in$end,
    # At the tree level
    h = data.tree.in$h,
    d = data.tree.in$d,
    a = data.tree.in$a,
    dbh = data.tree.in$dbh,
    comp = data.tree.in$comp,
    # priors
    d0 = ref.data.in$d0,
    d1 = ref.data.in$d1,
    d2 = ref.data.in$d2,
    d3 = ref.data.in$d3,
    e0 = ref.data.in$e0,
    e1 = ref.data.in$e1)
  
  return(out)
}



#' Fit the mortality model with true data on Dj
#' @param data_jags.in List contianing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @return A rjags object
fit_mortality_2 <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin){
  # Model
  model <- 
    "model{
  
  # Loop on all plots
  for(j in 1:Nplots){
    # Loop on trees
    for(i in beg[j]:end[j]){
      
      # Probability to observe harvested tree
      h[i] ~ dbern(ph[i])
      ph[i] = (1 - D[j])*phadBM[i] + D[j]*pdDj[i]*phdD[i]
      
      # Probability to observe dead tree
      d[i] ~ dbern(pd[i])
      pd[i] = (1 - D[j])*pdBM[i]*(1 - phadBM[i]) + D[j]*pdDj[i]*(1 - phdD[i])
      
      # Probability to observe alive tree
      a[i] ~ dbern(pa[i])
      pa[i] = (1 - D[j])*(1 - pdBM[i])*(1 - phadBM[i]) + D[j]*(1 - pdDj[i])
      
      # Probability to die in disturbed plots
      pdDj[i] = 1 - (1 - pdBM[i])*(1 - pdD[i])
      
      # Probability to die from background mortality
      logit(pdBM[i]) = b0 + b1*dbh[i] + b2*comp[i] + b3*sgdd[j] + b4*wai[j]
      
      # Probability to die from a disturbance
      logit(pdD[i]) = c0 + c1*DS[j] + c2*dbh[i]*DS[j]
      
      # Probability to be harvested knowing that the tree is dead from a disturbance
      logit(phdD[i]) = d0 + d1*dbh[i] + d2*DS[j] + d3*dbh[i]*DS[j]
      
      # Probability to be harvested in undisturbed plots
      logit(phadBM[i]) = e0 + e1*dbh[i]
    }
  }
  
  # Priors
  b0 ~ dnorm(0, 0.01)
  b1 ~ dnorm(0, 0.01)
  b2 ~ dnorm(0, 0.01)
  b3 ~ dnorm(0, 0.01)
  b4 ~ dnorm(0, 0.01)
  c0 ~ dnorm(0, 0.01)
  c1 ~ dnorm(0, 0.01)
  c2 ~ dnorm(0, 0.01)
  

  }"
  
  
  
  # Function to initialize priors
  initjags <- function(){
    return(list(b0 = runif(1, -10, 10),
                b1 = runif(1, -10, 10),
                b2 = runif(1, -10, 10),
                b3 = runif(1, -10, 10),
                b4 = runif(1, -10, 10),
                c0 = runif(1, -10, 10),
                c1 = runif(1, -10, 10),
                c2 = runif(1, -10, 10)))
  }
  
  # Run model
  out <- R2jags::jags(data = data_jags.in,
                      param = c("b0", "b1", "b2", "b3", 
                                "b4", "c0", "c1", "c2"),
                      inits = list(initjags(),
                                   initjags(),
                                   initjags()),
                      model.file = textConnection(model),
                      n.chains = n.chains,
                      n.iter = n.iter,
                      n.burnin = n.burn,
                      n.thin = n.thin,
                      DIC = TRUE, 
                      progress.bar = "text")
  return(out)
}





