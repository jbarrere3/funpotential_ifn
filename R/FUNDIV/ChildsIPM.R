###########################################################
#### FUNCTION TO BUILD IPM AND COMPUTE DEMO METRICS

## LOAD DATA WITH ENSEMBLE MODEL PROBABILITY OF PRESENCE FOR EACH SPECIES AND EACH PLOTS
load_data_pred <- function(path_data_pred = "data/data_proba_pred.rds"){
df <- readRDS(path_data_pred)
return(df)
}


################################################################################
## calculate the GL nodes and weights
################################################################################

# Gauss-Legendre quadrature nodes and weights on interval (L,U)
gaussQuadInt <- function(L, U, order=7) {
require(statmod)
  out <- gauss.quad(order) # GL is the default
  w <- out$weights; x <- out$nodes
  weights <- 0.5 * (U - L) * w
  nodes <- 0.5 * (U + L) + 0.5 * (U - L) * x
  return(list(weights = weights, nodes = nodes))
}


#################
### FUNCTIONS TO BUILD THE IPM

# Survival function:
s_x<- function(size, res, list_covs) {
  require(data.table)
  list_covs[["id"]] <- seq.int(length(size))
  country <- list_covs$country
  list_covs$country <- country[1]
  df_pred_p <- do.call(data.frame, list_covs)
  df_pred_p$size <- size[df_pred_p$id]
  df_pred_p$logsize <- log(df_pred_p$size)
  df_pred_p$yearsbetweensurveys <- 1
  df_pred_p$dead <- 1
  df_pred_p$country <- factor(df_pred_p$country, levels = country)
  # to try to speed up in matrix
  MM <- model.matrix(as.formula(res$formula), df_pred_p)
  if(length(country)>1){
  MM <- cbind("(Intercept)" = rep(1, length(size)),
              MM[, !grepl("country", colnames(MM)),drop=FALSE])
  }
  pred_c<- res$family$linkinv(MM %*% res$params_m)
  rm(MM)
  return(as.vector(1- pred_c)*(1-0.00508))
}
## (1-0.00508) is the proportion of tree not harvested per year in average
## across all species and countries this is coming from
## harvesting_rate line "mean"

# Conditional growth density function
g_x1x <- function(x1, x, res, list_covs) {
  require(lme4)
  require(data.table)
  pred.y <- rep(0, length.out = length(x1))
  index.y <- seq_len(length.out = length(x1))
  pos_sel <- x1>x
  if(sum(pos_sel) >0){
   country <- list_covs$country
   list_covs$country <- country[1]
   list_covs[["id"]] <- index.y[pos_sel]
   df_pred_p<- do.call(data.frame, list_covs)
   if (length(x)>1){
   df_pred_p$size <- x[df_pred_p$id]
   }else{
   df_pred_p$size <- x
   }
   df_pred_p$logsize <- log(df_pred_p$size)
   df_pred_p$logincr <- 1
   df_pred_p$country <- factor(df_pred_p$country, levels = country)
   MM <- model.matrix(as.formula(res$formula), df_pred_p)
   if(length(country)>1){
   MM <- cbind("(Intercept)" = rep(1, nrow(MM)),
               MM[, !grepl("country", colnames(MM)), drop = FALSE])
   }
   pred_c<- MM %*% res$params_m
   incr <- log((x1-x)[pos_sel])
   pred.y[pos_sel] <- dnorm(incr, mean = pred_c,
                            sd = res$sigma)*1/((x1 -x)[pos_sel])
 }
  return(pred.y)
}

# Growth kernel with mean growth computed outside the function
g_x1x_m<- function(x1, x, g_m, g_sigma) {
  pred.y <- rep(0, length.out = length(x1))
  index.y <- seq_len(length.out = length(x1))
  pos_sel <- x1>x

   d_x1_x <- x1-x
   d_x1_x[d_x1_x<0] <-  0
   incr <- log((x1-x)[pos_sel])
   pred.y[pos_sel] <- dnorm(incr, mean = g_m[pos_sel] ,
                            sd = g_sigma)*1/((x1 -x)[pos_sel])
  return(pred.y)
}

# mean growth
g_mean <- function(x, res, list_covs) {
  require(lme4)
  require(data.table)
   country <- list_covs$country
   list_covs$country <- country[1]
   list_covs[["size"]] <- x
   df_pred_p<- do.call(data.frame, list_covs)
   df_pred_p$logsize <- log(df_pred_p$size)
   df_pred_p$logincr <- 1
   df_pred_p$country <- factor(df_pred_p$country, levels = country)
   MM <- model.matrix(as.formula(res$formula), df_pred_p)
   if(length(country)>1){
   MM <- cbind("(Intercept)" = rep(1, nrow(MM)),
               MM[, !grepl("country", colnames(MM)), drop = FALSE])
   }
   pred_c<- MM %*% res$params_m
  return(pred_c)
}


################################################################################
################################################################################
## Kernel computation
################################################################################
################################################################################

################################################################################
## usual mid-point kernel function
################################################################################
P_x1_x <- function(x1, x, g_res, s_res, list_covs) {
  g_x1x(x1, x, g_res, list_covs) * s_x(x, s_res, list_covs)
}



# Gauss Legendre intergration on both dimension with eviction correction
mk_P_e_intMulti_Diag<- function(m, L, U, g_res, s_res, list_covs, level = 420,
                              diag_tresh= 50, correction = "ceiling") {
  h <- (U - L) / m
  meshpts <- L + ((1:m) - 1/2) * h
  df <- expand.grid(idx1=1:m, idx=1:m)
  if(!correction %in% c("ceiling", "sizeExtremes", "constant", "none"))
      stop("correction should be in ceiling, sizeExtremes, constant, or none")
  require(dplyr)
  df <- df %>% mutate(id = 1:n(),
                      G = 0,
                      x1 = meshpts[idx1],
                      x = meshpts[idx],
                      pos_sel = x1>=x,
                      diag_dist = sqrt(((x+x1)/2-x)^2+(x1-(x+x1)/2)^2),
                      diag_sel = diag_dist < diag_tresh & pos_sel,
                      up_diag_sel_mb = pos_sel & diag_dist >= diag_tresh)
  require(data.table)
  out1 <- gaussQuadInt(-h/2, h/2, 3)
  out2 <- gaussQuadInt(-h/2, h/2, floor(level/3))
  # create data for mean growth along x for the integration on dim 3
  quadx <- expand.grid(idx=1:m,
                       map1=seq.int(length(out1$weights)))
  quadx <- quadx %>% mutate(x = meshpts[idx] + out1$nodes[map1],
                            idx_3 = paste0(idx, c("A", "B", "C")[map1]))
  quadx <- as.data.table(quadx)
  quadx <- quadx[ , g_m:= g_mean(x, res = g_res,
                                 list_covs = list_covs)]
  g_mt <- quadx$g_m
  names(g_mt) <- quadx$idx_3
  # create data for integration on 3 x floor(level)
  quad <- expand.grid(id=df$id[df$diag_sel],
                      map1=seq.int(length(out1$weights)),
                      map2=seq.int(length(out2$weights)))
  # geat g_mean for all x values
  quad <- quad %>% mutate(x = df$x[id] + out1$nodes[map1],
                          x1 = df$x1[id] + out2$nodes[map2],
                          idx = df$idx[id],
                          idx1 = df$idx1[id],
                          idx_3 = paste0(idx, c("A", "B", "C")[map1]),
                          weights1 = out1$weights[map1],
                          weights2 = out2$weights[map2],
                          g_m= g_mt[idx_3])
  quad <- as.data.table(quad)
  # compute growth kernel
  quad <- quad[ , fvals := g_x1x_m(x1, x, g_m = g_m, g_sigma = g_res$sigma)]
  quad <- quad[ , G:= fvals *weights1*weights2]
  quad <- quad[, Gs:=sum(G), by = id]
  uniquad <- subset(unique(quad,by = "id"))
  rm(quad)
  res_GLe<- uniquad$Gs/h
  rm(uniquad)
  # mid-bin intergation for the rest
  df2 <- df %>% filter(up_diag_sel_mb)
  df2 <- as.data.table(df2)
  df2 <- df2[ , G := g_x1x(x1, x, res = g_res,
                                 list_covs = list_covs)*h]
  df$G[df$diag_sel] <- res_GLe
  df$G[df$up_diag_sel_mb] <- df2$G
  rm(df2, res_GLe)
  G <- df$G
  dim(G) <- c(m,m)
  rm(df)
  # eviction correction as in IPMpack
  if (correction == "constant"){# Based on IPMpack
   nvals <- colSums(G)
   P <- t((t(G)/nvals) * s_x(size = meshpts,
                             res = s_res,
                             list_covs = list_covs))
  }
  if (correction == "sizeExtremes"){# Based on IPMpack
   select_size_t<- meshpts > (U- 2*diag_tresh)
   DiffNvals <- pmax(1- colSums(G), 0)
   G[m, select_size_t] <- G[m, select_size_t] + DiffNvals[select_size_t]
   P <- G * s_x(size = meshpts, res = s_res, list_covs = list_covs)
  }
  if (correction == "ceiling"){# Based on Williams et al. 2012 Ecology
   DiffNvals <- pmax(1- colSums(G), 0)
   sx <- s_x(c(meshpts,U), res = s_res, list_covs = list_covs)
   G <- cbind(G, rep(0, length.out = length(meshpts)))
   G <- rbind(G, c(DiffNvals,1))
   P<- t(t(G) * sx)
  }
  if (correction == "none"){
   P <- t(t(G) * s_x(size = meshpts, res = s_res, list_covs = list_covs))
  }
  return(P)
}



#########################################################
# Arnaud version to speed up code with a band matrix and a matrix application of the weight


  fun_growth_mean <- function(g_res, list_covs, mesh_x){
  params_gr <- g_res$params_m
  params_i <- params_gr[!grepl("ntercept",names(params_gr)) &
                        !grepl("size",names(params_gr))] # Parameters without size
  params_i_inter <- params_i[grepl(":",names(params_i))] # Parameters with interaction
  L1 <- sub('.*:','',names(params_i_inter))
  L2 <- sub(':.*','',names(params_i_inter))
  K_i_inter <- sum(unlist(list_covs[L1])*unlist(list_covs[L2])*params_i_inter)
  params_i_no <- params_i[!(grepl(":",names(params_i)))] # Parameters without interaction
  K_i <- params_gr[grepl("ntercept",names(params_gr))] +
      sum(unlist(list_covs[names(params_i_no)]) * params_i_no) +
      K_i_inter # K for plot i
  ### Interaction between size and covariates has to be written as size:cov
  params_size_inter <- params_gr[grepl('size',names(params_gr)) &
                                 !grepl('logsize',names(params_gr)) &
                                 grepl(':',names(params_gr))]
  L <- sub('.*:','',names(params_size_inter))
  K_size <- params_gr['size'] + sum(unlist(list_covs[L]) * params_size_inter)
  params_logsize_inter <- params_gr[grepl('logsize',names(params_gr)) &
                                    grepl(':',names(params_gr))]
  L <- sub('.*:','',names(params_logsize_inter))
  K_logsize <- params_gr['logsize'] +
      sum(unlist(list_covs[L]) * params_logsize_inter)
  K_logsize[is.na(K_logsize)] <- 0
  mu_gr <- K_i + K_size * mesh_x + K_logsize * log(mesh_x)
  return(mu_gr)
  }

fun_surv_mean <- function(s_res, list_covs, mesh_x, weights1){
  params_sv <- s_res$params_m
  sig_sv <- s_res$sigma
  params_i <- params_sv[!grepl("ntercept",names(params_sv)) &
                        !grepl("size",names(params_sv))] # Parameters without size
  params_i_inter <- params_i[grepl(":",names(params_i))] # Parameters with interaction
  L1 <- sub('.*:','',names(params_i_inter))
  L2 <- sub(':.*','',names(params_i_inter))
  K_i_inter <- sum(unlist(list_covs[L1])*unlist(list_covs[L2])*params_i_inter)
  params_i_no <- params_i[!(grepl(":",names(params_i)))] # Parameters without interaction
  K_i <- params_sv[grepl("ntercept",names(params_sv))] +
      sum(unlist(list_covs[names(params_i_no)]) * params_i_no) +
      K_i_inter # K for plot i
  ### Interaction between size and covariates has to be written as size:cov
  params_size_inter <- params_sv[grepl('size',names(params_sv)) &
                                 !grepl('logsize',names(params_sv)) &
                                 grepl(':',names(params_sv))]
  L <- sub('.*:','',names(params_size_inter))
  K_size <- params_sv['size'] + sum(unlist(list_covs[L]) * params_size_inter)
  params_logsize_inter <- params_sv[grepl('logsize',names(params_sv)) &
                                    grepl(':',names(params_sv))]
  L <- sub('.*:','',names(params_logsize_inter))
  K_logsize <- params_sv['logsize'] +
      sum(unlist(list_covs[L]) * params_logsize_inter)
  K_logsize[is.na(K_logsize)] <- 0 # usefull?
  K_j <- sum(c(1,as.numeric(as.character(list_covs[names(params_sv)[2:length(params_sv)]]))) *
             params_sv,na.rm=TRUE)
  P_sv <- s_res$family$linkinv(K_i  + K_size * mesh_x + K_logsize * log(mesh_x))
  m <- length(mesh_x)/3
  P_sv <- P_sv[1:m] * weights1[1] +
      P_sv[(m+1):(2*m)] * weights1[2] +
      P_sv[(2*m+1):(3*m)] * weights1[3] # Avoid if no [1]
  return((1-P_sv)*(1 - 0.00508))
}
## (1-0.00508) is the proportion of tree not harvested per year in average
## across all species and countries this is coming from
## harvesting_rate line "mean"

fun_mid_int <- function(m,mesh_x_t, N_int, g_res, h, list_covs){
  require(dplyr)
  df <- expand.grid(dim1=  1:m, dim2=1:m) #mesh_x_t)
  df$x1 <- mesh_x_t[df$dim1]
  df$x <- mesh_x_t[df$dim2]
  df <- df %>% mutate(id = 1:n(),
                      P = 0,
                      triang_sel = dim1>dim2,
                      diag_dist = dim1 - dim2,
                      up_diag_sel_mb = diag_dist>=N_int)
  require(data.table)
  df2 <- df %>% filter(up_diag_sel_mb)
  df2 <- as.data.table(df2)
  df2 <- df2[ , P := g_x1x(x1, x, res = g_res,
                           list_covs = list_covs)*h]
  df$P[df$up_diag_sel_mb] <- df2$P
  P2<- df$P
  dim(P2) <- c(m,m)
  rm(df, df2)
return(P2)
}


# Gauss-Legendre 2 dimensions integration
# P_int percentage of the distance to the diagonal in cells where the GL integration is applied
mk_P_GL_2 <- function(m, L, U, g_res, s_res, list_covs,
                       diag_tresh= 50,level=420,correction="none",WMat){
  if(! correction %in% c("constant", "ceiling", "sizeExtremes", "none"))
      stop("correction must be in constant, ceiling, sizeExtremes, or none")
  level <- floor(level/3)
  # we inegrate on 3 diemnsion along the size at t and level/3 along size at t+1
  h <- (U - L) / m
  #build weight for GL integration on the two dim
  out1 <- gaussQuadInt(-h/2, h/2, 3) # For x integration
  weights1 <- out1$weights / sum(out1$weights) #equivalent to devided by h
  out2 <- gaussQuadInt(-h/2, h/2, level) # For x1 integration
  mesh_x <- seq(L+h/2,U-h/2,length.out=m)
  N_int <- sum((mesh_x-min(mesh_x))<diag_tresh)
  mesh_x_t <- mesh_x
  # vector for integraion on dim 3
  mesh_x <- as.vector(outer(mesh_x,out1$nodes,'+'))

### Growth
  mu_gr <- fun_growth_mean(g_res, list_covs, mesh_x)
  sig_gr <- g_res$sigma

### Survival
  P_sv <-   fun_surv_mean(s_res, list_covs, mesh_x, weights1)
### Create matrix for GL integration on dimension level
  mesh_x1B <- as.vector(t(outer(seq(0,(N_int-1)*h,length.out=N_int),out2$nodes,'+')))
  mesh_x1A <- mesh_x1B - out1$nodes[1] # to resacle the position on x
  mesh_x1B <- mesh_x1B - out1$nodes[2]
  mesh_x1C <- mesh_x1B - out1$nodes[3]

  temp <- function(d_x1_x,mu,sig=sig_gr){
    out <- rep(0, length.out = length(d_x1_x))
    out[d_x1_x>0] <- dnorm(log(d_x1_x[d_x1_x>0]),mu[d_x1_x>0],sig)*1/d_x1_x[d_x1_x>0]
    return(out)
  }

  P_incr <- outer(mesh_x1A,mu_gr[1:m],'temp') * weights1[1] +
      outer(mesh_x1B,mu_gr[(m+1):(2*m)],'temp') * weights1[2] +
      outer(mesh_x1C,mu_gr[(2*m+1):(3*m)],'temp') * weights1[3]

  if(missing(WMat)){
     WMat <- build_weight_matrix(out2$weights,N_int)
  }
  P_incr <- t(WMat) %*% P_incr
  P <- matrix(0,ncol=m,nrow=m)
  for (k in (1:m)){
    ind_k <- k:min(k+N_int-1,m)
    P[ind_k,k] <- P_incr[1:min(m-k+1,N_int),k]
  }

## ADD mid point integration for the rest of the triangular matrix
  P2 <- fun_mid_int(m ,mesh_x_t ,N_int ,g_res ,h, list_covs)
  P <- P+P2
  gc()
  if(correction == "constant"){
# Based on IPMpack
    nvals <- colSums(P)
    P <- t((t(P)/nvals))
    P <- P * (matrix(rep(t(P_sv),m),ncol=m,nrow=m))
  }

  if(correction ==  "sizeExtremes"){# Based on IPMpack
   selectsize_t <- (N_int + 0:(m-1)) > m
   DiffNvals <- pmax(1- colSums(P), 0)
   P[m,selectsize_t] <- P[m, selectsize_t] + DiffNvals[selectsize_t]
   P <- t(t(P) *P_sv)
  }
  if(correction ==  "none"){# Based on IPMpack
   P <- t(t(P) *P_sv)
  }
  if (correction == "ceiling"){# Based on Williams et al. 2012 Ecology integral towards infinity is not explicitely calculated
   P_sv_U<-   fun_surv_mean(s_res, list_covs, U +out1$nodes , weights1)
   nvals <- colSums(P)
   P <- rbind(P,pmax((1-nvals),0))
   P <- cbind(P,c(rep(0,m),1))
   P <- t(t(P) * c(P_sv, P_sv_U))
  }
  return(P)
}

build_weight_matrix <- function(weight,N_int){
  N_sub <- length(weight)
  ct <- ceiling(matrix(rep(1:(N_sub*N_int),N_sub),ncol=N_int,nrow=N_int*N_sub)/N_sub)
  ct2 <- t(matrix(rep(1:N_int,N_int*N_sub),ncol=N_int*N_sub,nrow=N_int))
  ind <- 0 * ct
  ind[ct==ct2] <- 1
  for (k in 1:N_int){
    ind[which(ind[,k]==1),k] <- weight
  }
  return(ind)
}


#############################################################################
### GENERAL FUNCTION TO BUILD IPM FROM FITTED OBJECT
get_mean_country_params <- function(params){
  intercept <- mean(params[ grep("country", names(params))])
  other_params <- params[!grepl("country", names(params))]
  params_m <- c("intercept" = intercept, other_params)

}
## FIT S AND G MODELS
fit_SG_models <- function(spsel, data, vars_s, vars_g){
  dir.create("output", showWarnings = FALSE)
  output_file <- file.path("output",paste(spsel,
                                          paste(unique(c(vars_s, vars_g)),
                                                collapse = "."),
                                          "rds", sep = "."))
  df <- format_data_survival(data, spsel)
  sv <- survival_size_climate(df, vars_s, strip = TRUE)
  if(length(unique(df$country))>1){
  params_s_m <- get_mean_country_params(sv$coefficients)
  }else{
  params_s_m <- sv$coefficients
  names(params_s_m)[1] <- "intercept"
  }
  s_res <- list(params_m = params_s_m,
                formula = sv$formula, family = sv$family )
  df <- format_data_growth(data[data$country %in% unique(df$country), ], spsel)
  # keep only country for wich mortality is fitted
  gr <- growth_size_climate(df, vars_g)
  if(length(unique(df$country))>1){
  params_g_m <- get_mean_country_params(fixef(gr))
  }else{
  params_g_m <- fixef(gr)
  names(params_g_m)[1] <- "intercept"
  }
  g_res <- list(params_m = params_g_m,
                formula = lme4:::nobars(formula(gr)), sigma = sigma(gr))
  list_res <- list(df = df , sv = s_res, gr = g_res)
  saveRDS(list_res, file = output_file)
  return(output_file)
}


## FIT S AND G MODELS selectiing best model based on AIC
fit_SG_models_best2<- function(spsel, data, saveRDS_TF = TRUE){
  dir.create("output", showWarnings = FALSE)
  output_file <- file.path("output",paste(spsel,"best2",
                                          "rds", sep = "."))
  df <- format_data_survival(data, spsel)
  sv1 <- survival_size_climate(df, c("sgddb", "waib"), strip = TRUE)
  sv2 <- survival_size_climate(df, c("sgddb", "waib", "size:waib",
                                       "logsize:waib",
                                        "size:sgddb", "logsize:sgddb"),
                               strip = TRUE)
  sv3 <- survival_size_climate(df, c("sgddb", "waib", "size:waib",
                                       "logsize:waib",
                                        "size:sgddb", "logsize:sgddb",
                                       "BATOTcomp:waib", "BATOTcomp:sgddb"),
                               strip = TRUE)
  sv4 <- survival_size_climate(df, c("sgdd", "wai", "sgdd2","wai2"),
                               strip = TRUE)
  sv5 <- survival_size_climate(df, c("sgdd", "wai", "sgdd2","wai2", "size:wai",
                                       "logsize:wai",
                                        "size:sgdd", "logsize:sgdd"),
                               strip = TRUE)
  sv6 <- survival_size_climate(df, c("sgdd", "wai", "sgdd2","wai2", "size:wai",
                                       "logsize:wai",
                                        "size:sgdd", "logsize:sgdd",
                                       "BATOTcomp:wai", "BATOTcomp:sgdd"),
                               strip = TRUE)
  sv <- list(sv1, sv2, sv3, sv4, sv5, sv6)[[which.min(c(sv1$aic, sv2$aic, sv3$aic,
                                                        sv4$aic, sv5$aic, sv6$aic))]]
  if(length(unique(df$country))>1){
  params_s_m <- get_mean_country_params(sv$coefficients)
  }else{
  params_s_m <- sv$coefficients
  names(params_s_m)[1] <- "intercept"
  }
  s_res <- list(params_m = params_s_m,
                formula = sv$formula, family = sv$family )
  df <- format_data_growth(data[data$country %in% unique(df$country), ], spsel)
  # keep only country for wich mortality is fitted
  gr1 <- growth_size_climate(df, c("sgddb", "waib"))
  gr2 <- growth_size_climate(df, c("sgddb", "waib", "size:waib",
                                       "logsize:waib",
                                        "size:sgddb", "logsize:sgddb"))
  gr3 <- growth_size_climate(df, c("sgddb", "waib", "size:waib",
                                       "logsize:waib",
                                        "size:sgddb", "logsize:sgddb",
                                       "BATOTcomp:waib", "BATOTcomp:sgddb"))
  gr4 <- growth_size_climate(df, c("sgdd", "wai", "sgdd2", "wai2"))
  gr5 <- growth_size_climate(df, c("sgdd", "wai", "sgdd2", "wai2",
                                   "size:wai", "logsize:wai",
                                   "size:sgdd", "logsize:sgdd"))
  gr6 <- growth_size_climate(df, c("sgdd", "wai", "sgdd2", "wai2",
                                   "size:wai", "logsize:wai",
                                   "size:sgdd", "logsize:sgdd",
                                   "BATOTcomp:wai", "BATOTcomp:sgdd"))
  gr <- list(gr1, gr2, gr3, gr4, gr5, gr6)[[which.min(c(AIC(gr1), AIC(gr2),
                                                        AIC(gr3),
                                                        AIC(gr4), AIC(gr5),
                                                        AIC(gr6)))]]
  grREML  <- update(gr, REML = TRUE)
  if(length(unique(df$country))>1){
  params_g_m <- get_mean_country_params(fixef(grREML))
  }else{
  params_g_m <- fixef(grREML)
  names(params_g_m)[1] <- "intercept"
  }
  g_res <- list(params_m = params_g_m,
                formula = lme4:::nobars(formula(grREML)), sigma = sigma(grREML))
  list_res <- list(df = df , sv = s_res, gr = g_res)
  if(saveRDS_TF){
  saveRDS(list_res, file = output_file)
  return(output_file)
  }else{
  return(list_res)
  }
}



## FIT S AND G MODELS selectiing best model based on AIC keeping only param and sd
fit_SG_models_best2_keep <- function(spsel, data){

  df <- format_data_survival(data, spsel)
  sv1 <- survival_size_climate(df, c("sgddb", "waib"), strip = F)
  sv2 <- survival_size_climate(df, c("sgddb", "waib", "size:waib",
                                       "logsize:waib",
                                        "size:sgddb", "logsize:sgddb"),
                               strip = F)
  sv3 <- survival_size_climate(df, c("sgddb", "waib", "size:waib",
                                       "logsize:waib",
                                        "size:sgddb", "logsize:sgddb",
                                       "BATOTcomp:waib", "BATOTcomp:sgddb"),
                               strip = F)
  sv4 <- survival_size_climate(df, c("sgdd", "wai", "sgdd2","wai2"),
                               strip = F)
  sv5 <- survival_size_climate(df, c("sgdd", "wai", "sgdd2","wai2", "size:wai",
                                       "logsize:wai",
                                        "size:sgdd", "logsize:sgdd"),
                               strip = F)
  sv6 <- survival_size_climate(df, c("sgdd", "wai", "sgdd2","wai2", "size:wai",
                                       "logsize:wai",
                                        "size:sgdd", "logsize:sgdd",
                                       "BATOTcomp:wai", "BATOTcomp:sgdd"),
                               strip = F)
    sv <- list(sv1, sv2, sv3, sv4, sv5, sv6)[[which.min(c(sv1$aic, sv2$aic,
                                                          sv3$aic,
                                                          sv4$aic, sv5$aic,
                                                          sv6$aic))]]
  sv_summary  <- summary(sv)$coefficients[, 1:2]
  df <- format_data_growth(data[data$country %in% unique(df$country), ], spsel)
  # keep only country for wich mortality is fitted
  gr1 <- growth_size_climate(df, c("sgddb", "waib"))
  gr2 <- growth_size_climate(df, c("sgddb", "waib", "size:waib",
                                       "logsize:waib",
                                        "size:sgddb", "logsize:sgddb"))
  gr3 <- growth_size_climate(df, c("sgddb", "waib", "size:waib",
                                       "logsize:waib",
                                        "size:sgddb", "logsize:sgddb",
                                       "BATOTcomp:waib", "BATOTcomp:sgddb"))
  gr4 <- growth_size_climate(df, c("sgdd", "wai", "sgdd2", "wai2"))
  gr5 <- growth_size_climate(df, c("sgdd", "wai", "sgdd2", "wai2",
                                   "size:wai", "logsize:wai",
                                   "size:sgdd", "logsize:sgdd"))
  gr6 <- growth_size_climate(df, c("sgdd", "wai", "sgdd2", "wai2",
                                   "size:wai", "logsize:wai",
                                   "size:sgdd", "logsize:sgdd",
                                   "BATOTcomp:wai", "BATOTcomp:sgdd"))
    gr <- list(gr1, gr2, gr3, gr4, gr5, gr6)[[which.min(c(AIC(gr1), AIC(gr2),
                                                          AIC(gr3),
                                                          AIC(gr4), AIC(gr5),
                                                          AIC(gr6)))]]
  grREML  <- update(gr, REML = TRUE)
  gr_summary  <- summary(grREML)$coefficients[, 1:2]
  return(list(gr = gr_summary, sv = sv_summary))
}

### COMPUTE VARIABLE IMPORTANCE  with package ingredients
variables_importance_sg <- function(spsel, data){
    require(lme4)
  df <- format_data_survival(data, spsel)
  sv1 <- survival_size_climate(df, c("sgddb", "waib"), strip = F)
  sv2 <- survival_size_climate(df, c("sgddb", "waib", "size:waib",
                                       "logsize:waib",
                                        "size:sgddb", "logsize:sgddb"),
                               strip = F)
  sv3 <- survival_size_climate(df, c("sgddb", "waib", "size:waib",
                                       "logsize:waib",
                                        "size:sgddb", "logsize:sgddb",
                                       "BATOTcomp:waib", "BATOTcomp:sgddb"),
                               strip = F)
  sv4 <- survival_size_climate(df, c("sgdd", "wai", "sgdd2","wai2"),
                               strip = F)
  sv5 <- survival_size_climate(df, c("sgdd", "wai", "sgdd2","wai2", "size:wai",
                                       "logsize:wai",
                                        "size:sgdd", "logsize:sgdd"),
                               strip = F)
  sv6 <- survival_size_climate(df, c("sgdd", "wai", "sgdd2","wai2", "size:wai",
                                       "logsize:wai",
                                        "size:sgdd", "logsize:sgdd",
                                       "BATOTcomp:wai", "BATOTcomp:sgdd"),
                               strip = F)
  sv <- list(sv1, sv2, sv3, sv4, sv5, sv6)[[which.min(c(sv1$aic, sv2$aic,
                                                        sv3$aic,
                                                        sv4$aic, sv5$aic,
                                                        sv6$aic))]]


  explain_sv <- DALEX::explain(sv,
                               data = df[,c("plotcode", "size", "logsize",
                                            "country", "BATOTcomp", "sgdd",
                                            "wai", "sgdd2","wai2",
                                            "sgddb", "waib",
                                            "yearsbetweensurveys")],
                         y = df$dead)

  var_imp_sv <- ingredients::feature_importance(explain_sv,
                               variable_groups =
                                    list("climate" = c("sgddb", "sgdd", "sgdd2",
                                                       "waib", "wai", "wai2"),
                                          "size" = c("size", "logsize"),
                                         "BATOTcomp" = c("BATOTcomp")),
                               B = 20, type = "difference")
  svB  <- var_imp_sv$dropout_loss[2:4]/sum(var_imp_sv$dropout_loss[2:4])
  var_imp_sv_auc <- ingredients::feature_importance(explain_sv,
                               variable_groups =
                                    list("climate" = c("sgddb", "sgdd", "sgdd2",
                                                       "waib", "wai", "wai2"),
                                          "size" = c("size", "logsize"),
                                         "BATOTcomp" = c("BATOTcomp")),
                               loss_function = DALEX::loss_one_minus_auc,
                               B = 20, type = "difference")
  sv_auc  <- var_imp_sv_auc$dropout_loss[2:4]/sum(var_imp_sv_auc$dropout_loss[2:4])

  df <- format_data_growth(data[data$country %in% unique(df$country), ], spsel)
  # keep only country for wich mortality is fitted
  gr1 <- growth_size_climate(df, c("sgddb", "waib"))
  gr2 <- growth_size_climate(df, c("sgddb", "waib", "size:waib",
                                       "logsize:waib",
                                        "size:sgddb", "logsize:sgddb"))
  gr3 <- growth_size_climate(df, c("sgddb", "waib", "size:waib",
                                       "logsize:waib",
                                        "size:sgddb", "logsize:sgddb",
                                       "BATOTcomp:waib", "BATOTcomp:sgddb"))
  gr4 <- growth_size_climate(df, c("sgdd", "wai", "sgdd2", "wai2"))
  gr5 <- growth_size_climate(df, c("sgdd", "wai", "sgdd2", "wai2",
                                   "size:wai", "logsize:wai",
                                   "size:sgdd", "logsize:sgdd"))
  gr6 <- growth_size_climate(df, c("sgdd", "wai", "sgdd2", "wai2",
                                   "size:wai", "logsize:wai",
                                   "size:sgdd", "logsize:sgdd",
                                   "BATOTcomp:wai", "BATOTcomp:sgdd"))
  gr <- list(gr1, gr2, gr3, gr4, gr5, gr6)[[which.min(c(AIC(gr1), AIC(gr2),
                                                        AIC(gr3),
                                                        AIC(gr4), AIC(gr5),
                                                        AIC(gr6)))]]
  grREML  <- update(gr, REML = TRUE)
  explain_gr <- DALEX::explain(grREML,
                         data = df[,c("plotcode", "size", "logsize",
                                      "country", "BATOTcomp", "sgdd",
                                      "wai", "sgdd2","wai2",
                                      "sgddb", "waib",
                                      "yearsbetweensurveys")],
                         y = df$logincr)
  var_imp_gr <- ingredients::feature_importance(explain_gr,
                               variable_groups =
                                    list("climate" = c("sgddb", "sgdd", "sgdd2",
                                                        "waib", "wai", "wai2"),
                                          "size" = c("size", "logsize"),
                                         "BATOTcomp" = c("BATOTcomp")),
                               B = 20, type = "difference")
  grB  <- var_imp_gr$dropout_loss[2:4]/sum(var_imp_gr$dropout_loss[2:4])
  names(grB)  <- names(svB)  <- names(sv_auc)  <- c("climate", "BATOTcomp",
                                                  "size")
  return(list(gr = grB, sv = svB, sv_auc = sv_auc))
}


# Format Variable Importancey

format_var_imp <- function(list_var_imp){
    tab_gr <- data.frame(species = names(list_var_imp),
                         t(sapply(list_var_imp, function(ll) ll$gr)))
    tab_sv <- data.frame(species = names(list_var_imp),
                         t(sapply(list_var_imp, function(ll) ll$sv)))
    tab_sv_auc <- data.frame(species = names(list_var_imp),
                             t(sapply(list_var_imp, function(ll) ll$sv_auc)))
    saveRDS(list(gr = tab_gr, sv = tab_sv, sv_auc = tab_sv_auc),
            file = "output/varImportanceBest2A.rds")
}

# To compute cross validation of growth and survival models
fit_SG_models_best2_cross<- function(spsel, data, saveRDS_TF = FALSE){
  dir.create("output", showWarnings = FALSE)
  output_file <- file.path("output",paste(spsel,"best2",
                                          "rds", sep = "."))
  df <- format_data_survival(data, spsel)
  sv1 <- survival_size_climate(df, c("sgddb", "waib"), strip = TRUE)
  sv2 <- survival_size_climate(df, c("sgddb", "waib", "size:waib",
                                       "logsize:waib",
                                        "size:sgddb", "logsize:sgddb"),
                               strip = TRUE)
  sv3 <- survival_size_climate(df, c("sgddb", "waib", "size:waib",
                                       "logsize:waib",
                                        "size:sgddb", "logsize:sgddb",
                                       "BATOTcomp:waib", "BATOTcomp:sgddb"),
                               strip = TRUE)
  sv4 <- survival_size_climate(df, c("sgdd", "wai", "sgdd2","wai2"),
                               strip = TRUE)
  sv5 <- survival_size_climate(df, c("sgdd", "wai", "sgdd2","wai2", "size:wai",
                                       "logsize:wai",
                                        "size:sgdd", "logsize:sgdd"),
                               strip = TRUE)
  sv6 <- survival_size_climate(df, c("sgdd", "wai", "sgdd2","wai2", "size:wai",
                                       "logsize:wai",
                                        "size:sgdd", "logsize:sgdd",
                                       "BATOTcomp:wai", "BATOTcomp:sgdd"),
                               strip = TRUE)
  sv_sel <- which.min(c(sv1$aic, sv2$aic, sv3$aic,
                        sv4$aic, sv5$aic, sv6$aic))
  sv <- list(sv1, sv2, sv3, sv4, sv5, sv6)[[sv_sel]]
  df <- format_data_growth(data[data$country %in% unique(df$country), ], spsel)
  # keep only country for wich mortality is fitted
  gr1 <- growth_size_climate(df, c("sgddb", "waib"))
  gr2 <- growth_size_climate(df, c("sgddb", "waib", "size:waib",
                                       "logsize:waib",
                                        "size:sgddb", "logsize:sgddb"))
  gr3 <- growth_size_climate(df, c("sgddb", "waib", "size:waib",
                                       "logsize:waib",
                                        "size:sgddb", "logsize:sgddb",
                                       "BATOTcomp:waib", "BATOTcomp:sgddb"))
  gr4 <- growth_size_climate(df, c("sgdd", "wai", "sgdd2", "wai2"))
  gr5 <- growth_size_climate(df, c("sgdd", "wai", "sgdd2", "wai2",
                                   "size:wai", "logsize:wai",
                                   "size:sgdd", "logsize:sgdd"))
  gr6 <- growth_size_climate(df, c("sgdd", "wai", "sgdd2", "wai2",
                                   "size:wai", "logsize:wai",
                                   "size:sgdd", "logsize:sgdd",
                                   "BATOTcomp:wai", "BATOTcomp:sgdd"))
  gr_sel <- which.min(c(AIC(gr1), AIC(gr2), AIC(gr3),
                      AIC(gr4), AIC(gr5), AIC(gr6)))
  gr <- list(gr1, gr2, gr3, gr4, gr5, gr6)[[gr_sel]]
  grREML  <- update(gr, REML = TRUE)
  list_res <- list(df = df , sv = sv, gr = grREML,
                   gr_sel = gr_sel, sv_sel = sv_sel)
  return(list_res)
}

# compute weight to resample data more equally along the climatic gradient
resample_data_pca_weight <- function(spsel, data, clim_pca, n_breaks = 40, prop = 0.7){
  data <- left_join(data, clim_pca$df_PC, by = "plotcode")

  data <- data %>% filter(sp == spsel)
  data_fit <- data %>% mutate(PC1_cut = cut(PC1, breaks = n_breaks, labels = FALSE)) %>%
                group_by(PC1_cut) %>% mutate(p_PC1_cut = 1/n()) %>% ungroup() %>%
                sample_frac(size = prop, weight = p_PC1_cut)
  data_eval <- data[!data$treecode2 %in% data_fit$treecode2, ]
  data_eval <- data_eval[data_eval$country %in% unique(data_fit$country), ]
  return(list(fit = data_fit, eval = data_eval))
}

## FIT S AND G MODELS with a resampling of the data
fit_SG_models_resample<- function(spsel, data, clim_pca, vars_s, vars_g){
  list_df <- resample_data_pca_weight(spsel, data, clim_pca)
  df <- format_data_survival(list_df$fit, spsel)
  sv <- survival_size_climate(df, vars_s, strip = TRUE)
  if(length(unique(df$country))>1){
  params_s_m <- get_mean_country_params(sv$coefficients)
  }else{
  params_s_m <- sv$coefficients
  }
  s_res <- list(params_m = params_s_m,
                formula = sv$formula, family = sv$family )
  df <- format_data_growth(list_df$fit[list_df$fit$country %in% unique(df$country), ], spsel)
  # keep only country for wich mortality is fitted
  gr <- growth_size_climate(df, vars_g)
  if(length(unique(df$country))>1){
  params_g_m <- get_mean_country_params(fixef(gr))
  }else{
  params_g_m <- fixef(gr)
  }
  g_res <- list(params_m = params_g_m,
                formula = lme4:::nobars(formula(gr)), sigma = sigma(gr))
  list_res <- list(df = df , sv = s_res, gr = g_res)
  return(list_res)
}

# cross validate the growth and survival models
crossvalidate_SG_models_resample<- function(spsel, data, clim_pca, vars_s, vars_g){
  list_df <- resample_data_pca_weight(spsel, data, clim_pca)
  df <- format_data_survival(list_df$fit, spsel)
  sv <- survival_size_climate(df, vars_s, strip = TRUE)
  df_eval<- format_data_survival(list_df$eval, spsel)
  pred <- predict(sv, newdata = df_eval, type = "response")
  df <- format_data_growth(list_df$fit[list_df$fit$country %in% unique(df$country), ], spsel)
  df_eval<- format_data_growth(list_df$eval[list_df$eval$country %in% unique(df$country), ], spsel)
  auc_sv <- as.numeric(pROC::auc(df_eval$dead, pred))
  # keep only country for wich mortality is fitted
  df <- format_data_growth(list_df$fit[list_df$fit$country %in% unique(df$country), ], spsel)
  df_eval<- format_data_growth(list_df$eval[list_df$eval$country %in% unique(df$country), ], spsel)
  gr <- growth_size_climate(df, vars_g)
  pred <- predict(gr, newdata = df_eval, re.form = NA)
  rmsd_gr <- sqrt(mean((pred - df_eval$logincr)^2))
  nrmsd_gr <- rmsd_gr/(max(df_eval$logincr) - min(df_eval$logincr))
  return(list(auc_sv, rmsd_gr, nrmsd_gr))
}


# FIT S AND G MODELS with a resampling of the data with AIC based best model selection
fit_SG_models_resample_best2<- function(spsel, data, clim_pca, vars_s, vars_g){
  list_df <- resample_data_pca_weight(spsel, data, clim_pca)
  df <- format_data_survival(list_df$fit, spsel)
  list_res <- fit_SG_models_best2(spsel, df, saveRDS_TF = FALSE)
  return(list_res)
}


# cross validate the growth and survival models with AIC based best model selection
crossvalidate_SG_models_resample_best2<- function(spsel, data, clim_pca, vars_s, vars_g){
  list_df <- resample_data_pca_weight(spsel, data, clim_pca)
  df <- format_data_survival(list_df$fit, spsel)
  df_eval_sv<- format_data_survival(list_df$eval[list_df$eval$country %in% unique(df$country), ], spsel)
  df_eval_gr<- format_data_growth(list_df$eval[list_df$eval$country %in% unique(df$country), ], spsel)
  list_res <- fit_SG_models_best2_cross(spsel, df, saveRDS_TF = FALSE)
  pred <- predict(list_res$sv, newdata = df_eval_sv, type = "response")
  auc_sv <- as.numeric(pROC::auc(df_eval_sv$dead, pred))
  pred <- predict(list_res$gr, newdata = df_eval_gr, re.form = NA)
  rmsd_gr <- sqrt(mean((pred - df_eval_gr$logincr)^2))
  nrmsd_gr <- rmsd_gr/(max(df_eval_gr$logincr) - min(df_eval_gr$logincr))
  rmsd_gre <- sqrt(mean((exp(pred) - exp(df_eval_gr$logincr))^2))
  nrmsd_gre <- rmsd_gr/(max(exp(df_eval_gr$logincr)) - min(exp(df_eval_gr$logincr)))
  aic_gr <- AIC(list_res$gr)
  aic_sv <- AIC(list_res$sv)
  sv_sel <- list_res$sv_sel
  gr_sel <- list_res$gr_sel
 return(c(auc_sv, rmsd_gr, nrmsd_gr, rmsd_gre, nrmsd_gre, aic_sv, aic_gr, sv_sel, gr_sel))
}


# cross validate the growth and survival models with AIC based best model selection
crossvalidate_SG_models_resample_best2_residual<- function(spsel, data, clim_pca, vars_s, vars_g){
  list_df <- resample_data_pca_weight(spsel, data, clim_pca)
  df <- format_data_survival(list_df$fit, spsel)
  df_resid_gr<- format_data_growth(list_df$fit[list_df$fit$country %in% unique(df$country), ], spsel)
  list_res <- fit_SG_models_best2_cross(spsel, df, saveRDS_TF = FALSE)
  pred <- predict(list_res$gr, newdata = df_resid_gr, re.form = NA)
  df_resid_gr$logincr_pred <- pred
 return(df_resid_gr)
}

## fit model with resample for plot 
pred_SG_models_resample_best2<- function(spsel, data, clim_pca){
  list_df <- resample_data_pca_weight(spsel, data, clim_pca)
  df <- format_data_survival(list_df$fit, spsel)
  list_res <- fit_SG_models_best2_cross(spsel, df, saveRDS_TF = FALSE)
  df <- format_data_survival(data, spsel)
  vars <- c("size", "sgdd", "wai", "BATOTcomp")
  pred_list <- vector("list")
  for (var in vars){
   res_gr <- predict_for_var_mean(var, vars, list_res$df, list_res$gr, df)
   res_sv <- predict_for_var_surv_mean(var, vars, list_res$df, list_res$sv, df)
   df_pred <- data.frame(seq_var = res_gr$seq_var,
                        gr = exp(res_gr$pred_m),
                        sv = res_sv$pred)
   df_pred$var <- var
   df_pred$index <- seq_len(length(res_gr$seq_var))
   pred_list[[var]] <- df_pred
  }
  df_res <- bind_rows(pred_list)
  return(df_res)
}


## cross validate growth and survival model N times
crossvalide_SG_resample_best2_N<- function(spsel= "Fagus sylvatica",
                                           data,
                                           clim_pca,
                                           vars_g = c("sgddb", "waib"),
                                           vars_s = c("sgddb", "waib"),
                                           N = 100){
 mat_cross <- matrix(NA, ncol = 9, nrow = N)
 colnames(mat_cross) <- c("auc_sv", "rmsd_gr", "nrmsd_gr", "rmsd_gre", "nrmsd_gre","AIC_sv", "AIC_gr", "sv_sel", "gr_sel")
 for (i in 1:N){
  mat_cross[i,] <- crossvalidate_SG_models_resample_best2(spsel, data, clim_pca, vars_s, vars_g)
 }
 print(spsel)
 return(mat_cross)
}

crossvalide_SG_resample_best2_residual_N<- function(spsel= "Fagus sylvatica",
                                           data,
                                           clim_pca,
                                           vars_g = c("sgddb", "waib"),
                                           vars_s = c("sgddb", "waib"),
                                           N = 100){
     dd <- crossvalidate_SG_models_resample_best2_residual(spsel, data, clim_pca, vars_s, vars_g)
     dd$repet <- 1
 print(spsel)
 return(dd)
}

## cross validate growth and survival model N times
pred_SG_resample_best2_N<- function(spsel= "Fagus sylvatica",
                                      data,
                                      clim_pca,
                                      N = 100){
 list_pred_resample<- vector("list")
 for (i in 1:N){
  temp <- pred_SG_models_resample_best2(spsel, data, clim_pca)
  temp$resample <- i
  list_pred_resample[[i]] <- temp
 }
 print(spsel)
 df_res <- bind_rows(list_pred_resample)
 return(df_res)
}


# FUNCTIONS TO FIT ALL RESPECTIVE VERSION OF THE MODEL FOR ALL SPECIES
#base model
make_fit_all_sp <- function(sps,data,
                             vars_g = c("sgddb", "waib"),
                             vars_s = c("sgddb", "waib")){
 ll <- lapply(sps$sp, fit_SG_models, data = data, vars_g = vars_g, vars_s = vars_s)
 names(ll) <- sps$sp
 return(ll)
}
# interactions size : clim
make_fit_inter_all_sp <- function(sps,data,
                                  vars_g = c("sgddb", "waib", "size:waib",
                                             "logsize:waib",
                                        "size:sgddb", "logsize:sgddb"),
                                  vars_s = c("sgddb", "waib", "size:waib",
                                             "logsize:waib",
                                        "size:sgddb", "logsize:sgddb")){
 ll <- lapply(sps$sp, fit_SG_models, data = data, vars_g = vars_g, vars_s = vars_s)
 names(ll) <- sps$sp
 return(ll)
}
# interactions size : clim and BATOTcomp : clim
make_fit_inter2_all_sp <- function(sps,data,
                                   vars_g = c("sgddb", "waib", "size:waib",
                                              "logsize:waib",
                                        "size:sgddb", "logsize:sgddb",
                                        "BATOTcomp:waib", "BATOTcomp:sgddb"),
                                   vars_s = c("sgddb", "waib", "size:waib",
                                              "logsize:waib",
                                        "size:sgddb", "logsize:sgddb",
                                        "BATOTcomp:waib", "BATOTcomp:sgddb")){
 ll <- lapply(sps$sp, fit_SG_models, data = data, vars_g = vars_g, vars_s = vars_s)
 names(ll) <- sps$sp
 return(ll)
}

# best model with interaction and poly
make_fit_best2_all_sp <- function(sps,data){
 ll <- lapply(sps$sp, fit_SG_models_best2, data = data)
 names(ll) <- sps$sp
 return(ll)
}

make_fit_best2_all_sp_keep <- function(sps,data){
 ll <- lapply(sps$sp, fit_SG_models_best2_keep, data = data)
 names(ll) <- sps$sp
 return(ll)
}

make_var_imp_best2_all_sp <- function(sps,data){
 ll <- lapply(sps$sp, variables_importance_sg, data = data)
 names(ll) <- sps$sp
 return(ll)
}


# GENERATE DATA FOR IPM FOR ONE PREDICTION
generate_data_pred <- function(vars_g, vars_s, df, vars_value = NA){
 vars_base <- c("BATOTcomp")
 vars_pred <- c(vars_base,unique(c(vars_g, vars_s)))
 #REMOVE INTERACTIONS
 vars_pred <- vars_pred[!grepl(":", vars_pred)]
 list_m <- lapply(df[ , vars_pred],  mean, na.rm = TRUE)
 list_m[["country"]] <- unique(df$country)
 if(any(!is.na(vars_value))){
   list_m[names(vars_value)] <- vars_value
 }
list_m$sgdd <- 1/list_m$sgddb
list_m$wai <- 1/list_m$waib - 1
list_m$sgdd2 <- list_m$sgdd^2
list_m$wai2 <- list_m$wai^2
 return(list_m)
}

generate_data_pred2 <- function(i, df, df2){
 list_m <-as.list(df[i,])
 list_m[["country"]] <- unique(df2$country)
list_m$sgdd <- 1/list_m$sgddb
list_m$wai <- 1/list_m$waib - 1
list_m$sgdd2 <- list_m$sgdd^2
list_m$wai2 <- list_m$wai^2
 return(list_m)
}



###################################
###################################
# Build IPM for plot i

# First version
make_IPM_i_e <- function(i, data_plots_pred, fit_sg,
                             spsel= "Fagus sylvatica",
                             vars_g = c("sgddb", "waib"),
                             vars_s = c("sgddb", "waib"),
                             m_size = 700,
                             level = 420,
                             correction = "ceiling"){
  minSize <- 100*0.9
  maxSize <- 1.1* max(fit_sg$df$dbh1, na.rm = TRUE)
  h <- (maxSize - minSize) / m_size
  meshpts <- minSize + ((1:m_size) - 1/2) * h
  list_m <- generate_data_pred(vars_g, vars_s, fit_sg$df,
                               vars_value = data_plots_pred[i, ])
  P <- mk_P_e_intMulti_Diag(m = m_size, L = minSize, U = maxSize,
                          g_res = fit_sg$gr, s_res = fit_sg$sv,
                          list_covs = list_m, diag_tresh= 50,
                          level = level, correction = correction)
  if(correction == "ceiling") {meshpts <- c(meshpts, maxSize)}
 return(list(P = P, meshpts = meshpts, list_m = list_m))
}

# Arnaud version
make_IPM_GL_2_i <- function(i, data_plots_pred, fit_sg,
                             spsel= "Fagus sylvatica",
                             vars_g = c("sgddb", "waib"),
                             vars_s = c("sgddb", "waib"),
                             m_size = 700,
                             level = 420,
			     correction = "ceiling",
			     diag_tresh = 50,
			     WMat){
  minSize <- 100*0.9
  maxSize <- 1.1* max(fit_sg$df$dbh1, na.rm = TRUE)
  h <- (maxSize - minSize) / m_size
  meshpts <- seq(minSize+h/2, maxSize-h/2, length.out = m_size)
  list_m <- generate_data_pred2(i, data_plots_pred, fit_sg$df)
  if (!missing(WMat)){
    P <- mk_P_GL_2(m = m_size, L = minSize, U = maxSize,
                          g_res = fit_sg$gr, s_res = fit_sg$sv,
			  list_covs = list_m, diag_tresh = diag_tresh,
			  level = level, correction = correction,WMat=WMat)
  }else{
    P <- mk_P_GL_2(m = m_size, L = minSize, U = maxSize,
                          g_res = fit_sg$gr, s_res = fit_sg$sv,
			  list_covs = list_m, diag_tresh = diag_tresh,
			  level = level, correction = correction)
  }
  if(correction == "ceiling") {meshpts <- c(meshpts, maxSize)}
 return(list(P = P, meshpts = meshpts, list_m = list_m))
}



###########################################################################
### FUNCTIONS TO EXPLORE IPM SIZE EFFECT ON METRICS ESTIMATION
###########################################################################

make_IPM_mSize_compare<- function(spsel= "Fagus sylvatica",
                           list_fit_sg,
                           vars_g = c("sgddb", "waib"),
                           vars_s = c("sgddb", "waib"),
                           m_size = 700,
                           level = 420){
print(paste("start ", spsel))
fit_sg <- readRDS(list_fit_sg[[spsel]])
data_plots_pred <- expand.grid(sgddb = quantile(fit_sg$df$sgddb,
                                                probs = c(0.5),
                                                na.rm = TRUE),
                               waib = quantile(fit_sg$df$waib,
                                                probs = c(0.5),
                                                na.rm = TRUE),
                               BATOTcomp = quantile(fit_sg$df$BATOTcomp,
                                                probs = c(0.5),
                                                na.rm = TRUE))
list_IPM_50<- make_IPM_i_e(1, data_plots_pred, fit_sg,
                   spsel= spsel,
                   vars_g = vars_g,
                   vars_s = vars_s,
                   m_size = 50,
                   level = level,
                   correction = "ceiling")
res_demo50 <- derive_demo_metrics_i(list_IPM_50,
                                   fit_sg, chosenSize = 600)
res_demo50 <- as.data.frame(res_demo50)
res_demo50$m_size <- 50
rm(list_IPM_50)
gc()

list_IPM_100<- make_IPM_i_e(1, data_plots_pred, fit_sg,
                   spsel= spsel,
                   vars_g = vars_g,
                   vars_s = vars_s,
                   m_size = 100,
                   level = level,
                   correction = "ceiling")
res_demo100 <- derive_demo_metrics_i(list_IPM_100,
                                   fit_sg, chosenSize = 600)
res_demo100 <- as.data.frame(res_demo100)
res_demo100$m_size <- 100
rm(list_IPM_100)
gc()
list_IPM_200<- make_IPM_i_e(1, data_plots_pred, fit_sg,
                   spsel= spsel,
                   vars_g = vars_g,
                   vars_s = vars_s,
                   m_size = 200,
                   level = level,
                   correction = "ceiling")
res_demo200 <- derive_demo_metrics_i(list_IPM_200,
                                   fit_sg, chosenSize = 600)
res_demo200 <- as.data.frame(res_demo200)
res_demo200$m_size <- 200
rm(list_IPM_200)
gc()
list_IPM_500<- make_IPM_i_e(1, data_plots_pred, fit_sg,
                   spsel= spsel,
                   vars_g = vars_g,
                   vars_s = vars_s,
                   m_size = 500,
                   level = level,
                   correction = "ceiling")
res_demo500 <- derive_demo_metrics_i(list_IPM_500, fit_sg,
                                   chosenSize = 600)
res_demo500 <- as.data.frame(res_demo500)
res_demo500$m_size <- 500
rm(list_IPM_500)
gc()

list_IPM_800<- make_IPM_i_e(1, data_plots_pred, fit_sg,
                   spsel= spsel,
                   vars_g = vars_g,
                   vars_s = vars_s,
                   m_size = 800,
                   level = level,
                   correction = "ceiling")
res_demo800 <- derive_demo_metrics_i(list_IPM_800, fit_sg,
                                   chosenSize = 600)
res_demo800 <- as.data.frame(res_demo800)
res_demo800$m_size <- 800
rm(list_IPM_800)
gc()

list_IPM_1100<- make_IPM_i_e(1, data_plots_pred, fit_sg,
                   spsel= spsel,
                   vars_g = vars_g,
                   vars_s = vars_s,
                   m_size = 1100,
                   level = level,
                   correction = "ceiling")
res_demo1100 <- derive_demo_metrics_i(list_IPM_1100, fit_sg,
                                   chosenSize = 600)
res_demo1100 <- as.data.frame(res_demo1100)
res_demo1100$m_size <- 1100
rm(list_IPM_1100)
gc()

list_IPM_2000<- make_IPM_i_e(1, data_plots_pred, fit_sg,
                   spsel= spsel,
                   vars_g = vars_g,
                   vars_s = vars_s,
                   m_size = 2000,
                   level = level,
                   correction = "ceiling")
res_demo2000 <- derive_demo_metrics_i(list_IPM_2000, fit_sg,
                                   chosenSize = 600)
res_demo2000 <- as.data.frame(res_demo2000)
res_demo2000$m_size <- 2000
rm(list_IPM_2000)
gc()

list_IPM_3000<- make_IPM_i_e(1, data_plots_pred, fit_sg,
                   spsel= spsel,
                   vars_g = vars_g,
                   vars_s = vars_s,
                   m_size = 3000,
                   level = level,
                   correction = "ceiling")
res_demo3000 <- derive_demo_metrics_i(list_IPM_3000, fit_sg,
                                   chosenSize = 600)
res_demo3000 <- as.data.frame(res_demo3000)
res_demo3000$m_size <- 3000
rm(list_IPM_3000)
gc()

## list_IPM_4000<- make_IPM_i_e(1, data_plots_pred, fit_sg,
##                    spsel= spsel,
##                    vars_g = vars_g,
##                    vars_s = vars_s,
##                    m_size = 4000,
##                    level = level,
##                    correction = "ceiling")
## res_demo4000 <- derive_demo_metrics_i(list_IPM_4000, fit_sg,
##                                    chosenSize = 600)
## res_demo4000 <- as.data.frame(res_demo4000)
## res_demo4000$m_size <- 4000
## rm(list_IPM_4000)
## gc()

res <- rbind(res_demo50, res_demo100,
             res_demo200, res_demo500,
             res_demo800, res_demo1100,
             res_demo2000, res_demo3000)#,
#             res_demo4000)
return(res)
gc()
}

format_compare_mSize <- function(i, spvec, ll){
sp <- spvec[i]
df <- ll[[i]]
df$sp <- sp
return(df)
}

## EVALUATE SENSITIVITY TO IPM MATRIX SIZE
make_IPM_metrics_compare_m_all_sp <- function(sps, list_fit_sg, m_size = 700,
                                                 level = 420){
 l <- lapply(sps$sp, make_IPM_mSize_compare,
             list_fit_sg = list_fit_sg, m_size = m_size ,
             level = level)
 return(l)
}

format_all_IPM_metrics_compare_m_all_sp <- function(sps, ll_compare, name = "base"){

df <- bind_rows(lapply(seq_len(length(sps$sp)), format_compare_mSize, sps$sp, ll_compare))
saveRDS(df, file = file.path("output", paste0("IPM_compare_mSize_", name, ".rds")))
}



###########################################################################
## FUNCTIONS to DERIVE Demographic metrics from IPM
###########################################################################

# Mean lifespan
survivorshipcurve <- function(P,meshpts, maxAge = 700){
h <- diff(meshpts)[1]
nBigMatrix <- nrow(P)
e <- matrix(1,nrow=1,ncol=dim(P)[1])
la <- rep(NA,maxAge)
offspring_prob <- c(0, 0, 0, 1, rep(0, nBigMatrix -4))
# TODO DO SOMETHING BETTER FOR INITIAL SIZE ?
# We can calculate la[1], survival to age 1, as
la[1] <- sum((e %*% P)*offspring_prob)
#Later survivorships require P^a so let's do the calculation recursively
Pa <- P
for(a in 2:maxAge){
	Pa=Pa %*% P
	la[a]= sum((e %*% Pa)*offspring_prob)
}
la <- c(1,la)
pa<- la[2:(maxAge+1)]/la[1:maxAge]
return(list(la=la, pa=pa))
}

survivorshipcurve_evic<- function(P,meshpts, maxAge = 1000){
h <- diff(meshpts)[1]
nBigMatrix <- nrow(P)
e <- matrix(1,nrow=1,ncol=dim(P)[1])
la <- matrix(NA, length(meshpts),maxAge)
offspring_prob <- c(0, 0, 0, 1, rep(0, nBigMatrix -4))
# TODO DO SOMETHING BETTER FOR INITIAL SIZE ?
# We can calculate la[1], survival to age 1, as
la[, 1] <- as.numeric(P%*%offspring_prob)
#Later survivorships require P^a so let's do the calculation recursively
for(a in 2:maxAge){
    la[, a]= as.numeric(P%*%la[, a-1])
}
return(la)
}


# from https://github.com/ipmbook/first-edition
lifespan <- function(P,meshpts){
nBigMatrix <- nrow(P)
N <- solve(diag(nBigMatrix)-P, sparse = TRUE, tol = .Machine$double.eps/10000)
e <- matrix(1,nrow=1,ncol=dim(P)[1])
offspring_prob <- c(1, rep(0, nBigMatrix -1))
mean.age.death <- round(sum((e %*% N)*offspring_prob) -1,4)
Var.nu <- round(sum((e %*% (2 * N %*% N - N))* offspring_prob)  -
                   sum((e %*% N * offspring_prob ))^2,4)
return(list(mean_lifespan = mean.age.death,
            variance_lifespan = Var.nu,
            CV_lifespan = sqrt(Var.nu)/mean.age.death))
}


# size at death kernel from https://github.com/ipmbook/first-editions
size_at_death_kernel <- function(P, meshpts, list_covs, s_res){
h <- diff(meshpts)[1]
nBigMatrix <- nrow(P)
N <- solve(diag(nBigMatrix)-P, sparse = TRUE, tol = .Machine$double.eps/10000)
e <- matrix(1,nrow=1,ncol=dim(P)[1])
offspring_prob <- c(1, rep(0, nBigMatrix -1))
Omega <- (1-s_x(meshpts,s_res, list_covs)) * N

# each of the columns defines a probability distribution and so should sum to 1, let's check
## round(e %*% Omega,4)
#all 1's as it should be.

# then the distribution of sizes at death for the offspring distribution is
dist.size.death <- (Omega %*% offspring_prob)

# 1st the mean is
mean.size.death     <-round(sum(dist.size.death*meshpts),4)
# 2nd the variance is
var.size.death <- round(sum(dist.size.death*meshpts*meshpts) -
          sum(dist.size.death*meshpts)*sum(dist.size.death*meshpts),4)

return(list(mean.size.death = mean.size.death,
            var.size.death = var.size.death,
            Omega = Omega))
}

# From Cochran and Ellner 1992
passageTime <- function(chosenSize,P, meshpts){
 P <- as.matrix(P)
 require(MASS)
 loc <- which(abs(chosenSize-meshpts) ==
 	     min(abs(chosenSize - meshpts)),arr.ind=TRUE)[1]
 matrix.dim <- length(P[1,])
 Tprime <- P
 Tprime[,loc] <- 0
 N <- solve(diag(matrix.dim)-Tprime)
 N2 <- N %*% N
 time.to.absorb <- (N2/N)[loc, ]

 return(time.to.absorb)
}

#' Calculates Keyfitz' entropy from https://github.com/jonesor/Rage/blob/master/R/kEntropy.R Need to check with Rob
kEntropyIPM <- function(P,meshpts, maxAge = 1000, trapeze = FALSE){
h <- diff(meshpts)[1]
nBigMatrix <- nrow(P)
e <- matrix(1,nrow=1,ncol=dim(P)[1])
la <- rep(NA,maxAge)
offspring_prob <- c(1, rep(0, nBigMatrix -1))
# We can calculate la[1], survival to age 1, as
la[1] <- sum((e %*% P)*offspring_prob);
#Later survivorships require P^a so let's do the calculation recursively
Pa <- P
for(a in 2:maxAge){
	Pa=Pa %*% P
	la[a]= sum((e %*% Pa)*offspring_prob)
}
la <- c(1,la)
la <- la[!is.na(la) & la !=0]
# if zero lead to problem in log(la). Compute kE on the non zero part (check with Rob)
if(trapeze == TRUE){
  ma <- function(x,n=2){stats::filter(x,rep(1/n,n), sides=2)}
  la2 <- na.omit(as.vector(ma(la)))
  return(-sum(la2*log(la2))/sum(la2))
}else{
    return(-sum(la*log(la))/sum(la))
}
}

##########################################
## function to compute all demo metrics with error handling
derive_demo_metrics_i <- function(ll, fit_sg, chosenSize = 600, kE_TF = FALSE){
res_lifespan <- tryCatch(lifespan(ll$P, ll$meshpts),
                         error = function(e) return(list(mean_lifespan = NA,
                                                         variance_lifespan = NA,
                                                         CV_lifespan = NA)))
SizePassageTime <- min(max(ll$meshpts), chosenSize)
res_passageTime <- tryCatch(passageTime(chosenSize = SizePassageTime,P = ll$P,
                                        meshpts = ll$meshpts),
                            error = function(e) return(NA))
if (all(is.na(res_passageTime))){
res_passageTime_100 <- NA
}else{
res_passageTime_100 <- res_passageTime[which.min(abs(ll$meshpts - 100))]
}

#print("passage time")
res_sizeDeath <- tryCatch(size_at_death_kernel(P = ll$P, meshpts = ll$meshpts,
                                               list_covs = ll$list_m,
                                               s_res = fit_sg$sv)[1:2],
                          error = function(e) return(list(mean.size.death = NA,
                                                          var.size.death = NA)))
#print("size death made")
if(kE_TF){
res_kE <- kEntropyIPM(ll$P,ll$meshpts, maxAge = 1000)
#print("kE made")
return(list(lifespan = res_lifespan,
            passageTime_100 = res_passageTime_100,
            SizePassageTime = SizePassageTime,
            sizeDeath = res_sizeDeath,
            kE = res_kE))
}else{
return(list(lifespan = res_lifespan,
            passageTime_100 = res_passageTime_100,
            SizePassageTime = SizePassageTime,
            sizeDeath = res_sizeDeath))
}
}

#########################
#########################
## FUNCTIONS TO DERIVE DEMO METRICS FOR DIFFERENT DATA

# COMPUTE FOR ONE RAW OF THE DATA FOR PREDICTION
# Build IPM and derive demo metrics
make_IPM_e_metrics_i<- function(i, data_plots_pred, fit_sg,
                             spsel= "Fagus sylvatica",
                             vars_g = c("sgddb", "waib"),
                             vars_s = c("sgddb", "waib"),
                             m_size = 700,
                             level = 25,
                             chosenSize = 600,
                             fun_IPM_i = make_IPM_i_e,
                             kE_TF = FALSE){
 res_IPM <- fun_IPM_i(i, data_plots_pred, fit_sg,
                               spsel, vars_g, vars_s,
                               m_size, level)
 # print(paste("IPMe made for i = " , i, spsel))
 res_demo <- derive_demo_metrics_i(res_IPM, fit_sg, chosenSize = chosenSize,
                                   kE_TF = kE_TF)
 # add mean growth and survival
 list_m <- generate_data_pred(vars_g, vars_s, fit_sg$df,
                              vars_value = data_plots_pred[i, ])
 g_m <- g_mean(c(150, 500),fit_sg$gr, list_m)
 res_demo$g_m_150 <- g_m[1]
 res_demo$g_m_500 <- g_m[2]
 s_m <- s_x(c(150, 500), fit_sg$sv, list_m)
 res_demo$s_m_150 <- s_m[1]
 res_demo$s_m_500 <- s_m[2]
 res_demo <- as.data.frame(res_demo)
 # TEST IF PB IN IPM a value on the diag greater than 1 will lead to trouble and make no sense
 if (any(diag(res_IPM$P)>1)){
   res_demo[1,!grepl("_m_", names(res_demo))] <- NA
 }
# print(paste("all metrics made for i = " , i, spsel))
 return(res_demo)
 gc()
}

########################################################
####### BUILD IPM DEMO PRED FOR DIFFERENT DATA, ALL RAWS
########################################################

# for the mean sgdd and wai of the species
make_IPM_e_metrics_mean<- function(spsel= "Fagus sylvatica",
                           list_fit_sg,
                           vars_g = c("sgddb", "waib"),
                           vars_s = c("sgddb", "waib"),
                           m_size = 700,
                           level = 420,
                           fun_IPM_i = make_IPM_i_e){
fit_sg <- readRDS(list_fit_sg[[spsel]])
data_plots_pred <- expand.grid(sgddb = 1/quantile(1/fit_sg$df$sgddb,
                                    probs = 0.5, na.rm = TRUE),
                                waib = 1/(quantile(1/fit_sg$df$waib -1,
                                    probs = 0.5, na.rm = TRUE)+1),
                                BATOTcomp = c(0,30))
data_plots_pred$sgdd <- 1/data_plots_pred$sgddb
data_plots_pred$wai <- 1/data_plots_pred$waib - 1
data_plots_pred$sgdd2 <- data_plots_pred$sgdd^2
data_plots_pred$wai2 <- data_plots_pred$wai^2
list_demo1 <- lapply(seq_len(length.out = nrow(data_plots_pred)),
                    make_IPM_e_metrics_i, data_plots_pred, fit_sg,
                    spsel= spsel,
                    vars_g = vars_g,
                    vars_s = vars_s,
                    m_size = m_size, level = level, fun_IPM_i = fun_IPM_i, kE_TF = TRUE)
# add mean surv and mean growth
list_m <- generate_data_pred(vars_g, vars_s, fit_sg$df,
                               vars_value = data_plots_pred[1, ])
s_0 <- s_x(200, fit_sg$sv, list_m)
g_0 <- exp(g_mean(200, fit_sg$gr, list_m))
list_m <- generate_data_pred(vars_g, vars_s, fit_sg$df,
                               vars_value = data_plots_pred[2, ])
s_30 <- s_x(200, fit_sg$sv, list_m)
g_30 <- exp(g_mean(200, fit_sg$gr, list_m))

df1 <- do.call("rbind", lapply(list_demo1, function(x) as.data.frame(x)))
df1 <- cbind(data_plots_pred, df1)
df1$g_mean <- c(g_0, g_30)
df1$s_mean <- c(s_0, s_30)
print(paste(spsel, m_size, " done mean"))
return(df1)
}

# For a sequence of sgdd and wai
make_IPM_e_metrics_2Vars<- function(spsel= "Fagus sylvatica",
                           list_fit_sg,
                           vars_g = c("sgddb", "waib"),
                           vars_s = c("sgddb", "waib"),
                           m_size = 700,
                           level = 420,
                           TwoClimVar = FALSE,
                           fun_IPM_i = make_IPM_i_e, kE_TF = FALSE){
fit_sg <- readRDS(list_fit_sg[[spsel]])
#sgdd
data_plots_pred1 <- expand.grid(sgddb = 1/seq(from = quantile(1/fit_sg$df$sgddb,
                                                  probs = 0.03, na.rm = TRUE),
                                              to = quantile(1/fit_sg$df$sgddb,
                                                  probs = 0.97, na.rm = TRUE),
                                              length.out = 30),
                               waib = mean( fit_sg$df$waib, na.rm = TRUE),
                               BATOTcomp = c(0,30))

data_plots_pred1$sgdd <- 1/data_plots_pred1$sgddb
data_plots_pred1$wai <- 1/data_plots_pred1$waib - 1
data_plots_pred1$sgdd2 <- data_plots_pred1$sgdd^2
data_plots_pred1$wai2 <- data_plots_pred1$wai^2

list_demo1 <- lapply(seq_len(length.out = nrow(data_plots_pred1)),
                    make_IPM_e_metrics_i, data_plots_pred1, fit_sg,
                    spsel= spsel,
                    vars_g = vars_g,
                    vars_s = vars_s,
                    m_size = m_size, level = level, fun_IPM_i = fun_IPM_i, kE_TF = kE_TF)
df1 <- do.call("rbind", lapply(list_demo1, function(x) as.data.frame(x)))
df1 <- cbind(data_plots_pred1, df1)

print(paste(spsel, m_size, " done sgddb"))
if(TwoClimVar){
#wai
data_plots_pred2 <- expand.grid(waib = 1/(seq(from = quantile(1/fit_sg$df$waib - 1, probs = 0.03,na.rm = TRUE),
                                              to = quantile(1/fit_sg$df$waib - 1, probs = 0.97, na.rm = TRUE),
                                              length.out = 30) +1),
                               sgddb = mean( fit_sg$df$sgddb, na.rm = TRUE),
                                BATOTcomp = c(0,30))
data_plots_pred2$sgdd <- 1/data_plots_pred2$sgddb
data_plots_pred2$wai <- 1/data_plots_pred2$waib - 1
data_plots_pred2$sgdd2 <- data_plots_pred2$sgdd^2
data_plots_pred2$wai2 <- data_plots_pred2$wai^2
list_demo2 <- lapply(seq_len(length.out = nrow(data_plots_pred2)),
                    make_IPM_e_metrics_i, data_plots_pred2, fit_sg,
                    spsel= spsel,
                    vars_g = vars_g,
                    vars_s = vars_s,
                    m_size = m_size,
                    level = level, fun_IPM_i = fun_IPM_i, kE_TF = kE_TF)
df2 <- do.call("rbind", lapply(list_demo2, function(x) as.data.frame(x)))
df2 <- cbind(data_plots_pred2, df2)

print(paste(spsel, m_size, " done waib"))
return(list(sgdd = df1, wai = df2))
}else{
return(df1)
}
}

# For a sequence of axis 1 of pca of climate
make_IPM_e_metrics_pca <- function(spsel= "Fagus sylvatica",
                           clim_pca,
                           list_fit_sg,
                           vars_g = c("sgddb", "waib"),
                           vars_s = c("sgddb", "waib"),
                           m_size = 700,
                           level = 420,
                           fun_IPM_i = make_IPM_i_e, kE_TF = FALSE){
fit_sg <- readRDS(list_fit_sg[[spsel]])
#pca
df <- left_join(fit_sg$df, clim_pca$df_PC, by = "plotcode")
seq_PC1 <- seq(from = quantile(df$PC1, probs = 0.03, na.rm = TRUE),
               to = quantile(df$PC1, probs = 0.97, na.rm = TRUE),
               length.out = 30)
mat <- as.matrix(data.frame(PC1 = seq_PC1,
                            PC2 = rep(mean(clim_pca$df_PC$PC2, na.rm = TRUE), 30)))
rot <- clim_pca$pca$rotation
scal <- clim_pca$pca$scale
cent <- clim_pca$pca$center
pred_back<- mat %*% t(rot)
sgdd_seq <- pred_back[, 1]*scal[1]+cent[1]
wai_seq <- pred_back[, 2]*scal[2]+cent[2]
sgddb_seq <- 1/sgdd_seq
waib_seq <- 1/(wai_seq+1)

data_plots_pred <- rbind(data.frame(sgddb = sgddb_seq, waib = waib_seq,
                                    sgdd = sgdd_seq, wai = wai_seq,
                                    sgdd2 = sgdd_seq^2, wai2 = wai_seq^2,
                                    BATOTcomp = 0),
                         data.frame(sgddb = sgddb_seq, waib = waib_seq,
                                    sgdd = sgdd_seq, wai = wai_seq,
                                    sgdd2 = sgdd_seq^2, wai2 = wai_seq^2,
                                    BATOTcomp = 30))
list_demo <- lapply(seq_len(length.out = nrow(data_plots_pred)),
                    make_IPM_e_metrics_i, data_plots_pred, fit_sg,
                    spsel= spsel,
                    vars_g = vars_g,
                    vars_s = vars_s,
                    m_size = m_size, level = level,
                    fun_IPM_i = fun_IPM_i, kE_TF = kE_TF)
df <- do.call("rbind", lapply(list_demo, function(x) as.data.frame(x)))
df <- cbind(data_plots_pred, df)
df$PC1 <- seq_PC1
return(df)
}



#######################################################
## WITH a resampling of the data repeated 100 times
make_IPM_e_metrics_pca_resample<- function(spsel= "Fagus sylvatica",
                                           data,
                                           clim_pca,
                                           name,
                                           vars_g = c("sgddb", "waib"),
                                           vars_s = c("sgddb", "waib"),
                                           m_size = 700,
                                           level = 420,
                                           fun_IPM_i = make_IPM_i_e,
                                           kE_TF = FALSE){
res_list <- vector("list")
df <- format_data_survival(data, spsel)
df <- format_data_growth(data[data$country %in% unique(df$country), ], spsel)
df <- left_join(df, clim_pca$df_PC, by = "plotcode")
seq_PC1 <- seq(from = quantile(df$PC1, probs = 0.03, na.rm = TRUE),
             to = quantile(df$PC1, probs = 0.97, na.rm = TRUE),
             length.out = 30)
mat <- as.matrix(data.frame(PC1 = seq_PC1,
                          PC2 = rep(mean(clim_pca$df_PC$PC2, na.rm = TRUE),
                              30)))
rot <- clim_pca$pca$rotation
scal <- clim_pca$pca$scale
cent <- clim_pca$pca$center
pred_back<- mat %*% t(rot)
sgdd_seq <- pred_back[, 1]*scal[1]+cent[1]
wai_seq <- pred_back[, 2]*scal[2]+cent[2]
sgddb_seq <- 1/sgdd_seq
waib_seq <- 1/(wai_seq+1)

data_plots_pred <- rbind(data.frame(sgddb = sgddb_seq, waib = waib_seq,
                                  BATOTcomp = 0),
                       data.frame(sgddb = sgddb_seq, waib = waib_seq,
                                  BATOTcomp = 30))
data_plots_pred$sgdd <- 1/data_plots_pred$sgddb
data_plots_pred$wai <- 1/data_plots_pred$waib - 1
data_plots_pred$sgdd2 <- data_plots_pred$sgdd^2
data_plots_pred$wai2 <- data_plots_pred$wai^2

for (i in 1:100){
fit_sg <- fit_SG_models_resample(spsel, data, clim_pca, vars_s, vars_g)

list_demo <- lapply(seq_len(length.out = nrow(data_plots_pred)),
                    make_IPM_e_metrics_i, data_plots_pred, fit_sg,
                    spsel= spsel,
                    vars_g = vars_g,
                    vars_s = vars_s,
                    m_size = m_size, level = level,
                    fun_IPM_i = fun_IPM_i, kE_TF = kE_TF)

df <- do.call("rbind", lapply(list_demo, function(x) as.data.frame(x)))
df <- cbind(data_plots_pred, df)
df$PC1 <- seq_PC1
res_list[[i]] <- df
}
file_name <- paste0("output/resample_pca", name, spsel, ".rds")
saveRDS(res_list, file = file_name)
print(paste0("resample done for ", spsel))
return(file_name)
}


# resampling and selection of best model for each resample based on AIC

make_IPM_e_metrics_pca_resample_best2_i<-  function(spsel, data, clim_pca,
                                                   vars_s, vars_g, N_seq,
                                                   m_size, level,
                                                   fun_IPM_i, kE_TF,
                                                    data_plots_pred){
fit_sg <- fit_SG_models_resample_best2(spsel, data, clim_pca, vars_s, vars_g)
list_demo <- lapply(seq_len(length.out = nrow(data_plots_pred)),
                    make_IPM_e_metrics_i, data_plots_pred, fit_sg,
                    spsel= spsel,
                    vars_g = vars_g,
                    vars_s = vars_s,
                    m_size = m_size, level = level,
                    fun_IPM_i = fun_IPM_i, kE_TF = kE_TF)

gc()
df <- do.call("rbind", lapply(list_demo, function(x) as.data.frame(x)))
df <- cbind(data_plots_pred, df)
return(df)
gc()
}

make_IPM_e_metrics_pca_resample_best2<- function(spsel= "Fagus sylvatica",
                                           data,
                                           clim_pca,
                                           name,
                                           vars_g = c("sgddb", "waib"),
                                           vars_s = c("sgddb", "waib"),
                                           m_size = 700,
                                           level = 420,
                                           fun_IPM_i = make_IPM_i_e,
                                           kE_TF = FALSE){
print(paste0("start build IPM resample best2 species : ", spsel))
## df <- format_data_survival(data, spsel)
## df <- format_data_growth(data[data$country %in% unique(df$country), ], spsel)
## df <- left_join(df, clim_pca$df_PC, by = "plotcode")
## seq_PC1 <- seq(from = quantile(df$PC1, probs = 0.03, na.rm = TRUE),
##              to = quantile(df$PC1, probs = 0.97, na.rm = TRUE),
##              length.out = 30)
## mat <- as.matrix(data.frame(PC1 = seq_PC1,
##                           PC2 = rep(mean(clim_pca$df_PC$PC2, na.rm = TRUE),
##                               30)))
## rot <- clim_pca$pca$rotation
## scal <- clim_pca$pca$scale
## cent <- clim_pca$pca$center
## pred_back<- mat %*% t(rot)
## sgdd_seq <- pred_back[, 1]*scal[1]+cent[1]
## wai_seq <- pred_back[, 2]*scal[2]+cent[2]
## sgddb_seq <- 1/sgdd_seq
## waib_seq <- 1/(wai_seq+1)

## data_plots_pred <- rbind(data.frame(sgddb = sgddb_seq, waib = waib_seq,
##                                   BATOTcomp = 0, PC1 = seq_PC1),
##                        data.frame(sgddb = sgddb_seq, waib = waib_seq,
##                                   BATOTcomp = 30, PC1 = seq_PC1))
## data_plots_pred$sgdd <- 1/data_plots_pred$sgddb
## data_plots_pred$wai <- 1/data_plots_pred$waib - 1
## data_plots_pred$sgdd2 <- data_plots_pred$sgdd^2
## data_plots_pred$wai2 <- data_plots_pred$wai^2

## res_list <- vector("list")
##     N_seq <-  30
## for (i in 1:100){
##    df <- make_IPM_e_metrics_pca_resample_best2_i(spsel, data, clim_pca,
##                                                          vars_s, vars_g,
##                                                          N_seq = N_seq,
##                                                          m_size = m_size,
##                                                          level = level,
##                                                          fun_IPM_i = fun_IPM_i,
##                                                          kE_TF = kE_TF,
##                                                          data_plots_pred = data_plots_pred)
##     print(paste0("done for ", spsel, i))
##     res_list[[i]] <- df
##     gc()
## }
## print(paste0("resample done for ", spsel))
file_name <- paste0("output/resample_pca", name, spsel, ".rds")
## saveRDS(res_list, file = file_name)
return(file_name)
## rm(res_list)
gc()
}



# along sgdd and wai sequence
make_IPM_e_metrics_2Vars_resample_best2_i<- function(spsel, data, clim_pca,
                                                     vars_s, vars_g,
                                                     N_seq, m_size, level,
                                                     fun_IPM_i, kE_TF,
                                                     data_plots_pred1, data_plots_pred2){
fit_sg <- fit_SG_models_resample_best2(spsel, data, clim_pca, vars_s, vars_g)
list_demo1 <- lapply(seq_len(length.out = nrow(data_plots_pred1)),
                    make_IPM_e_metrics_i, data_plots_pred1, fit_sg,
                    spsel= spsel,
                    vars_g = vars_g,
                    vars_s = vars_s,
                    m_size = m_size, level = level, fun_IPM_i = fun_IPM_i,
                    kE_TF = kE_TF)
df1 <- do.call("rbind", lapply(list_demo1, function(x) as.data.frame(x)))
df1 <- cbind(data_plots_pred1, df1)
#wai
list_demo2 <- lapply(seq_len(length.out = nrow(data_plots_pred2)),
                    make_IPM_e_metrics_i, data_plots_pred2, fit_sg,
                    spsel= spsel,
                    vars_g = vars_g,
                    vars_s = vars_s,
                    m_size = m_size,
                    level = level, fun_IPM_i = fun_IPM_i, kE_TF = kE_TF)
gc()
df2 <- do.call("rbind", lapply(list_demo2, function(x) as.data.frame(x)))
df2 <- cbind(data_plots_pred2, df2)
return(list(df1, df2))
gc()
}

make_IPM_e_metrics_2Vars_resample_best2<- function(spsel= "Fagus sylvatica",
                                           data,
                                           clim_pca,
                                           name,
                                           vars_g = c("sgddb", "waib"),
                                           vars_s = c("sgddb", "waib"),
                                           m_size = 700,
                                           level = 420,
                                           fun_IPM_i = make_IPM_i_e,
                                           kE_TF = FALSE){
res_list_sgdd<- vector("list")
res_list_wai<- vector("list")
N_seq <- 30
df <- format_data_survival(data, spsel)
df <- format_data_growth(data[data$country %in% unique(df$country), ], spsel)
df <- left_join(df, clim_pca$df_PC, by = "plotcode")
data_plots_pred1 <- expand.grid(sgddb = 1/seq(from = quantile(1/df$sgddb,
                                                  probs = 0.03, na.rm = TRUE),
                                              to = quantile(1/df$sgddb,
                                                  probs = 0.97, na.rm = TRUE),
                                              length.out = N_seq),
                               waib = mean( df$waib, na.rm = TRUE),
                               BATOTcomp = c(0,30))
data_plots_pred1$sgdd <- 1/data_plots_pred1$sgddb
data_plots_pred1$wai <- 1/data_plots_pred1$waib - 1
data_plots_pred1$sgdd2 <- data_plots_pred1$sgdd^2
data_plots_pred1$wai2 <- data_plots_pred1$wai^2
data_plots_pred2 <- expand.grid(waib = 1/(seq(from = quantile(1/df$waib - 1, probs = 0.03,na.rm = TRUE),
                                              to = quantile(1/df$waib - 1, probs = 0.97, na.rm = TRUE),
                                              length.out = N_seq) +1),
                               sgddb = mean( df$sgddb, na.rm = TRUE),
                                BATOTcomp = c(0,30))
data_plots_pred2$sgdd <- 1/data_plots_pred2$sgddb
data_plots_pred2$wai <- 1/data_plots_pred2$waib - 1
data_plots_pred2$sgdd2 <- data_plots_pred2$sgdd^2
data_plots_pred2$wai2 <- data_plots_pred2$wai^2

for (i in 1:100){
   df <- tryCatch(make_IPM_e_metrics_2Vars_resample_best2_i(spsel, data,
                                                         clim_pca,
                                                         vars_s, vars_g,
                                                         N_seq = 30,
                                                         m_size = m_size,
                                                         level = level,
                                                         fun_IPM_i = fun_IPM_i,
                                                         kE_TF = kE_TF,
                                                         data_plots_pred1,
                                                         data_plots_pred2),
                          error = function(e) return(list(NA, NA)))
   res_list_sgdd[[i]] <- df[[1]]
   res_list_wai[[i]] <- df[[2]]
   gc()
}
print(paste(spsel, m_size, " done resample 100"))
file_name <- paste0("output/resample_2Vars", name, spsel, ".rds")
saveRDS(list(sgdd = res_list_sgdd, wai = res_list_wai), file = file_name)
return(file_name)
rm(res_list_sgdd, res_list_wai)
gc()
}

#########################
#### FUNCTIONS to EXTRACT SPECIES DISTRIBUTION CARACTERISTIC

#' Weighted quantile
#'
#' Function copied from **spatstat** package.
#'
#' @param x Vector of values
#' @param w Vector of weights
#' @param probs Vector of probabilities
#' @param na.rm Ignore missing data?
#' @export
weighted.quantile <- function(x, w, probs=seq(0,1,0.25), na.rm=TRUE) {
  x <- as.numeric(as.vector(x))
  w <- as.numeric(as.vector(w))
  if(anyNA(x) || anyNA(w)) {
    ok <- !(is.na(x) | is.na(w))
    x <- x[ok]
    w <- w[ok]
  }
  stopifnot(all(w >= 0))
  if(all(w == 0)) stop("All weights are zero", call.=FALSE)
  #'
  oo <- order(x)
  x <- x[oo]
  w <- w[oo]
  Fx <- cumsum(w)/sum(w)
  #'
  result <- numeric(length(probs))
  for(i in seq_along(result)) {
    p <- probs[i]
    lefties <- which(Fx <= p)
    if(length(lefties) == 0) {
      result[i] <- x[1]
    } else {
      left <- max(lefties)
      result[i] <- x[left]
      if(Fx[left] < p && left < length(x)) {
        right <- left+1
        y <- x[left] + (x[right]-x[left]) * (p-Fx[left])/(Fx[right]-Fx[left])
        if(is.finite(y)) result[i] <- y
      }
    }
  }
  names(result) <- paste0(format(100 * probs, trim = TRUE), "%")
  return(result)
}

# mean position on axis 1 of pca
get_mean_sp_pca_i<- function(spsel, data, data_pred,
                                  clim_pca){
spsel
if (spsel == "Betula") {
plotsel <- data_pred$plotcode[data_pred[, gsub(" ", ".", spsel)] ==1]
Betulasp <- c("Betula.pendula_Proba", "Betula.pubescens_Proba")
proba <- apply(data_pred[ , Betulasp], 1, mean)
names(proba) <-  data_pred$plotcode
}else{
if(spsel == "Pinus uncinata"){
spsel2 <- "Pinus mugo"
}else{
spsel2 <- spsel
}
plotsel <- data_pred$plotcode[data_pred[, gsub(" ", ".", spsel)] ==1]
proba <- data_pred[,
                   paste0(gsub(" ", ".", spsel2), "_Proba")]
names(proba) <- data_pred$plotcode
}

df_t <- data[!duplicated(data$plotcode), ]
df <- merge(df_t, data.frame(plotcode = names(proba),
                             proba = proba,
                             PresAbsObs = data_pred[, gsub(" ", ".", spsel)]),
            by = "plotcode")
#pca
df<- left_join(df, clim_pca$df_PC, by = "plotcode")
mean_PC1 <- mean(df$PC1[df$PresAbsObs == 1], na.rm = TRUE)
median_PC1 <- median(df$PC1[df$PresAbsObs == 1], na.rm = TRUE)
ql_PC1 <- quantile(df$PC1[df$PresAbsObs == 1], probs = 0.05, na.rm = TRUE)
qh_PC1 <- quantile(df$PC1[df$PresAbsObs == 1], probs = 0.95, na.rm = TRUE)
mean_PC1_proba<- weighted.mean(df$PC1, df$proba, na.rm = TRUE)
median_PC1_proba <- weighted.quantile(df$PC1, df$proba, probs = 0.5, na.rm = TRUE)
ql_PC1_proba <- weighted.quantile(df$PC1, df$proba, probs = 0.05, na.rm = TRUE)
qh_PC1_proba <- weighted.quantile(df$PC1, df$proba, probs = 0.95, na.rm = TRUE)
res <- data.frame(sp = spsel,
                  mean_PC1 = mean_PC1, median_PC1 = median_PC1,
                  ql_PC1 = ql_PC1, qh_PC1 = qh_PC1,
                  mean_PC1_proba = mean_PC1_proba, median_PC1_proba= median_PC1_proba,
                  ql_PC1_proba = ql_PC1_proba, qh_PC1_proba = qh_PC1_proba)
return(res)
}

#######################
## PREDICT IPM METRICS FOR ALL SPECIES

# Predic IPM metrics along sgdd and wai
make_IPM_e_metrics_all_sp_TwoClimVar_BaseA<- function(sps, list_fit_sg, name = "BaseA",
                                                    m_size = 700, level = 420,
                                                    fun_IPM_i = make_IPM_GL_2_i){
 library(future.apply)
 plan(multiprocess, workers = 10)
 set.seed(123)
 l <- future_lapply(sps$sp, make_IPM_e_metrics_2Vars, list_fit_sg = list_fit_sg, m_size = m_size ,
             level = level, TwoClimVar = TRUE, fun_IPM_i = fun_IPM_i,
             future.seed = TRUE)
 saveRDS(l, file = file.path("output", paste0("MetricsTwoClimVar", name, ".rds")))
}


make_IPM_e_metrics_all_sp_TwoClimVar_Inter2A<- function(sps, list_fit_sg, name = "Inter2A",
                                                    m_size = 700, level = 420,
                                                    fun_IPM_i = make_IPM_GL_2_i){
 library(future.apply)
 set.seed(123)
 plan(multiprocess, workers = 10)
 l <- future_lapply(sps$sp, make_IPM_e_metrics_2Vars, list_fit_sg = list_fit_sg, m_size = m_size ,
             level = level, TwoClimVar = TRUE, fun_IPM_i = fun_IPM_i,
             future.seed = TRUE)
 saveRDS(l, file = file.path("output", paste0("MetricsTwoClimVar", name, ".rds")))
}


make_IPM_e_metrics_all_sp_TwoClimVar_Best2A<- function(sps, list_fit_sg, name = "Best2A",
                                                    m_size = 700, level = 420,
                                                    fun_IPM_i = make_IPM_GL_2_i){
 library(future.apply)
 set.seed(123)
 plan(multiprocess, workers = 10)
 l <- future_lapply(sps$sp, make_IPM_e_metrics_2Vars, list_fit_sg = list_fit_sg, m_size = m_size ,
             level = level, TwoClimVar = TRUE, fun_IPM_i = fun_IPM_i,
             future.seed = TRUE)
 saveRDS(l, file = file.path("output", paste0("MetricsTwoClimVar", name, ".rds")))
}


make_IPM_e_metrics_all_sp_2Vars_Best2A_Resample<- function(sps, clim_pca, data,
                                               name = "Best2A",
                                               m_size = 700, level = 420,
                                               fun_IPM_i = make_IPM_GL_2_i){
 library(future.apply)
 set.seed(123)
 plan(multiprocess, workers = 7)
 l <- future_lapply(sps$sp, make_IPM_e_metrics_2Vars_resample_best2,
                    data = data,
             clim_pca = clim_pca, name = name,
             vars_g = c("sgddb", "waib", "size:waib", "logsize:waib",
                                        "size:sgddb", "logsize:sgddb",
                 "BATOTcomp:waib", "BATOTcomp:sgddb"),
             vars_s = c("sgddb", "waib", "size:waib", "logsize:waib",
                                        "size:sgddb", "logsize:sgddb",
                 "BATOTcomp:waib", "BATOTcomp:sgddb"),
             m_size = m_size ,
             level = level, fun_IPM_i = fun_IPM_i,
             future.seed = TRUE)
 saveRDS(l, file = file.path("output", paste0("Metrics2VarsResample",
                name, ".rds")))
}

##  PREDIC IPM METRICS ALONG AXIS ONE OF CLIMATIC PCA
make_IPM_e_metrics_all_sp_pca_BaseA<- function(sps, clim_pca, list_fit_sg, name = "BaseA",
                                                    m_size = 700, level = 420,
                                                    fun_IPM_i = make_IPM_GL_2_i){
 l <- lapply(sps$sp, make_IPM_e_metrics_pca, clim_pca = clim_pca,
             list_fit_sg = list_fit_sg, m_size = m_size ,
             level = level, fun_IPM_i = fun_IPM_i)
 saveRDS(l, file = file.path("output", paste0("MetricsPCAVar", name, ".rds")))
}

make_IPM_e_metrics_all_sp_pca_BaseA_Resample<- function(sps, clim_pca, data,
                                               name = "BaseA",
                                               m_size = 700, level = 420,
                                               fun_IPM_i = make_IPM_GL_2_i){
 library(future.apply)
 set.seed(123)
 plan(multiprocess, workers = 10)
 l <- future_lapply(sps$sp, make_IPM_e_metrics_pca_resample, data = data,
             clim_pca = clim_pca, name = name,
              m_size = m_size ,
             level = level, fun_IPM_i = fun_IPM_i,
             future.seed = TRUE)
 saveRDS(l, file = file.path("output", paste0("MetricsPCAVarResample",
                name, ".rds")))
}


make_IPM_e_metrics_all_sp_pca_Inter2A<- function(sps, clim_pca, list_fit_sg, name = "Inter2A",
                                                    m_size = 700, level = 420,
                                                    fun_IPM_i = make_IPM_GL_2_i){
 library(future.apply)
 set.seed(123)
 plan(multiprocess, workers = 8)
 l <- future_lapply(sps$sp, make_IPM_e_metrics_pca,
             clim_pca = clim_pca,
             list_fit_sg = list_fit_sg, m_size = m_size ,
             level = level, fun_IPM_i = fun_IPM_i,
             future.seed = TRUE)
 saveRDS(l, file = file.path("output", paste0("MetricsPCAVar", name, ".rds")))
}

make_IPM_e_metrics_all_sp_pca_Inter2A_Resample<- function(sps, clim_pca, data,
                                               name = "Inter2A",
                                               m_size = 700, level = 420,
                                               fun_IPM_i = make_IPM_GL_2_i){
 library(future.apply)
 set.seed(123)
 plan(multiprocess, workers = 10)
 l <- future_lapply(sps$sp, make_IPM_e_metrics_pca_resample, data = data,
             clim_pca = clim_pca, name = name,
             vars_g = c("sgddb", "waib", "size:waib", "logsize:waib",
                                        "size:sgddb", "logsize:sgddb",
                 "BATOTcomp:waib", "BATOTcomp:sgddb"),
             vars_s = c("sgddb", "waib", "size:waib", "logsize:waib",
                                        "size:sgddb", "logsize:sgddb",
                 "BATOTcomp:waib", "BATOTcomp:sgddb"),
             m_size = m_size ,
             level = level, fun_IPM_i = fun_IPM_i,
             future.seed = TRUE)
 saveRDS(l, file = file.path("output", paste0("MetricsPCAVarResample",
                name, ".rds")))
}


make_IPM_e_metrics_all_sp_pca_Best2A<- function(sps, clim_pca, list_fit_sg, name = "Best2A",
                                                    m_size = 700, level = 420,
                                                    fun_IPM_i = make_IPM_GL_2_i){
 library(future.apply)
 plan(multiprocess, workers = 10)
 l <- future_lapply(sps$sp, make_IPM_e_metrics_pca, clim_pca = clim_pca,
             list_fit_sg = list_fit_sg, m_size = m_size ,
             level = level, fun_IPM_i = fun_IPM_i)
 saveRDS(l, file = file.path("output", paste0("MetricsPCAVar", name, ".rds")))
}

make_IPM_e_metrics_all_sp_pca_Best2A_Resample<- function(sps, clim_pca, data,
                                               name = "Best2A",
                                               m_size = 700, level = 420,
                                               fun_IPM_i = make_IPM_GL_2_i){
 library(future.apply)
 set.seed(123)
 print("start all species")
 plan(multiprocess, workers = 27)
 spvec <- sps$sp
 names(spvec) <- spvec
 l <- future_lapply(spvec, make_IPM_e_metrics_pca_resample_best2, data = data,
             clim_pca = clim_pca, name = name,
             vars_g = c("sgddb", "waib", "size:waib", "logsize:waib",
                                        "size:sgddb", "logsize:sgddb",
                 "BATOTcomp:waib", "BATOTcomp:sgddb"),
             vars_s = c("sgddb", "waib", "size:waib", "logsize:waib",
                                        "size:sgddb", "logsize:sgddb",
                 "BATOTcomp:waib", "BATOTcomp:sgddb"),
             m_size = m_size ,
             level = level, fun_IPM_i = fun_IPM_i,
             future.seed = TRUE)
 saveRDS(l, file = file.path("output", paste0("MetricsPCAVarResample",
                name, ".rds")))
}


###############################################################
## PREDICT IPM METRICS AT MEDIAN CLIMATE OF SPECIES
make_IPM_e_metrics_all_sp_mean<- function(sps, list_fit_sg, m_size = 700, level = 420,
                                                      fun_IPM_i = make_IPM_i_e){
 ## library(future.apply)
 ## plan(multiprocess, workers = 10)
 l <- lapply(sps$sp, make_IPM_e_metrics_mean, list_fit_sg = list_fit_sg, m_size = m_size ,
             level = level, fun_IPM_i = fun_IPM_i)
 return(l)
}

make_IPM_e_metrics_all_sp_mean_A<- function(sps, list_fit_sg, m_size = 700, level = 420,
                                                fun_IPM_i = make_IPM_GL_2_i){
 ## library(future.apply)
 ## plan(multiprocess, workers = 10)
 l <- lapply(sps$sp,make_IPM_e_metrics_mean, list_fit_sg= list_fit_sg, m_size = m_size ,
             level= level, fun_IPM_i= fun_IPM_i)
 print("done")
 return(l)
}



## COMPUTE CROSSVALIDATION ALL SPECIES
crossvalidate_all_sp_pca_Best2A_Resample<- function(sps, clim_pca, data,
                                               name = "Best2A"){
 library(future.apply)
 set.seed(123)
 plan(multiprocess, workers = 20)

 l <- future_lapply(sps$sp, crossvalide_SG_resample_best2_N, data = data,
             clim_pca = clim_pca,
             vars_g = c("sgddb", "waib", "size:waib", "logsize:waib",
                        "size:sgddb", "logsize:sgddb",
                        "BATOTcomp:waib", "BATOTcomp:sgddb"),
             vars_s = c("sgddb", "waib", "size:waib", "logsize:waib",
                        "size:sgddb", "logsize:sgddb",
                        "BATOTcomp:waib", "BATOTcomp:sgddb"),
             future.seed = TRUE)
 return(l)
}


crossvalidate_all_sp_pca_Best2A_residual_Resample<- function(sps, clim_pca, data,
                                               name = "Best2A"){
 library(future.apply)
 set.seed(123)
 plan(multiprocess, workers = 20)

 l <- future_lapply(sps$sp, crossvalide_SG_resample_best2_residual_N, data = data,
             clim_pca = clim_pca,
             vars_g = c("sgddb", "waib", "size:waib", "logsize:waib",
                        "size:sgddb", "logsize:sgddb",
                        "BATOTcomp:waib", "BATOTcomp:sgddb"),
             vars_s = c("sgddb", "waib", "size:waib", "logsize:waib",
                        "size:sgddb", "logsize:sgddb",
                        "BATOTcomp:waib", "BATOTcomp:sgddb"),
             future.seed = TRUE)
 return(l)
}


format_residual_resample_sp <-  function(df){
   df$resid <- df$logincr - df$logincr_pred
   return(df)
}

format_residual_resample <-  function(l){
  ll <- lapply(l, format_residual_resample_sp)
  return(ll)
}

pred_all_sp_pca_Best2A_Resample<- function(sps, clim_pca, data,
                                               name = "Best2A"){
 library(future.apply)
 set.seed(123)
 plan(multiprocess, workers = 15)

 l <- future_lapply(sps$sp, pred_SG_resample_best2_N, data = data,
             clim_pca = clim_pca,
            future.seed = TRUE)
 return(l)
}


format_crossvalidate_all<- function(sps, ll,
                                    name = "Best2A"){
 f <- function(i,sps,ll) {
     df_t <- as.data.frame(ll[[i]])
     df_t$sp <- sps$sp[i]
     return(df_t)
 }

 df<-bind_rows(lapply(seq_len(length.out = length(sps$sp)), f, sps, ll))
 res <- df %>% group_by(sp) %>% summarise(auc_sv = mean(auc_sv, na.rm = TRUE),
                                          rmsd_gr = mean(rmsd_gr, na.rm = TRUE),
                                          nrmsd_gr = mean(nrmsd_gr, na.rm = TRUE),
                                          rmsd_gre = mean(rmsd_gre, na.rm = TRUE),
                                          nrmsd_gre = mean(nrmsd_gre, na.rm = TRUE))
 saveRDS(res, file = file.path("output", paste0("crossvalidate",name, ".rds")))
}

format_crossmodelselect_all<- function(sps, ll,
                                    name = "Best2A"){
 f <- function(i,sps,ll) {
     df_t <- as.data.frame(ll[[i]])
     df_t$sp <- sps$sp[i]
     return(df_t)
 }

 df<-bind_rows(lapply(seq_len(length.out = length(sps$sp)), f, sps, ll))
 tab_sv <- table(df$sp, factor(df$sv_sel, levels = 1:6))
 colnames(tab_sv) <- paste0("model_", colnames(tab_sv))
 df_sv <-  data.frame(sp = rownames(tab_sv), as.data.frame.matrix(tab_sv))
 tab_gr <- table(df$sp, factor(df$gr_sel, levels = 1:6))
 colnames(tab_gr) <- paste0("model_", colnames(tab_gr))
 df_gr <-  data.frame(sp = rownames(tab_gr), as.data.frame.matrix(tab_gr))
 saveRDS(list(sv = df_sv, gr = df_gr), file = file.path("output", paste0("crossmodelselect",name, ".rds")))
}


format_varImportance_all<- function(sps, ll,
                                    name = "Best2A"){
 f <- function(i,sps,ll) {
     df_t <- as.data.frame(ll[[i]])
     df_t$sp <- sps$sp[i]
     return(df_t)
 }

 df<-bind_rows(lapply(seq_len(length.out = length(sps$sp)), f, sps, ll))

 res <- df %>% group_by(sp, var) %>% summarise(varImpGr = range(gr)[2] - range(gr)[1],
                                               varImpSv = range(sv)[2] - range(sv)[1])
 tab_sv <- table(df$sp, factor(df$sv_sel, levels = 1:6))
 colnames(tab_sv) <- paste0("model_", colnames(tab_sv))
 df_sv <-  data.frame(sp = rownames(tab_sv), as.data.frame.matrix(tab_sv))
 tab_gr <- table(df$sp, factor(df$gr_sel, levels = 1:6))
 colnames(tab_gr) <- paste0("model_", colnames(tab_gr))
 df_gr <-  data.frame(sp = rownames(tab_gr), as.data.frame.matrix(tab_gr))
 saveRDS(list(sv = df_sv, gr = df_gr), file = file.path("output", paste0("crossmodelselect",name, ".rds")))
}


##############################################################
##############################################################
########## FORMAT IPM PREDICTONS
##############################################################
##############################################################


##############################################
########## FORMAT OUT PER SPECIES WITH PROBA PRESENCE
#### MERGE PREDICTED IPM WITH PROBA OF PRESENCE

format_IPM_metrics_TwoClimVar<- function(i, ll, spvec, data, data_pred){
spsel <- spvec[i]
if (spsel == "Betula") {
plotsel <- data_pred$plotcode[data_pred[, gsub(" ", ".", spsel)] ==1]
Betulasp <- c("Betula.pendula_Proba", "Betula.pubescens_Proba")
proba <- apply(data_pred[, Betulasp],
               1, mean)
names(proba) <-  data_pred$plotcode
}else{
 if(spsel == "Pinus uncinata"){
  spsel2 <- "Pinus mugo"
 }else{
  spsel2 <- spsel
 }
 plotsel <- data_pred$plotcode[data_pred[, gsub(" ", ".", spsel)] ==1]
 proba <- data_pred[,
                    paste0(gsub(" ", ".", spsel2), "_Proba")]
 names(proba) <-  data_pred$plotcode
}
df_t<- data[ , ]
df_t <- df_t[!duplicated(df_t$plotcode), ]
df <- merge(df_t, data.frame(plotcode = names(proba), proba = proba),
            by = "plotcode")
df_sgdd<- ll[[i]]$sgdd
df_wai<- ll[[i]]$wai
# get good cut point to match IPM predict
seq_sgdd<- unique(1/df_sgdd$sgddb)
seq_wai<- unique(1/df_wai$waib - 1)
mid_bin_sgdd<- mean(seq_sgdd[-1] - seq_sgdd[-length(seq_sgdd)])/2
mid_bin_wai<- mean(seq_wai[-1] - seq_wai[-length(seq_wai)])/2
breaks_sgdd <- seq(from = min(seq_sgdd) -mid_bin_sgdd,
                to = max(seq_sgdd)+mid_bin_sgdd,
                length.out = length(seq_sgdd)+1)
breaks_wai <- seq(from = min(seq_wai) -mid_bin_wai,
                to = max(seq_wai)+mid_bin_wai,
                length.out = length(seq_wai)+1)
df <- df %>% mutate(sgdd_cut = cut(sgdd, breaks = breaks_sgdd, label = FALSE),
                    wai_cut = cut(wai, breaks = breaks_wai, label = FALSE))
data_sgdd <- df %>% filter(sgdd > min(breaks_sgdd) &
                           sgdd < max(breaks_sgdd) ) %>%
                    group_by(sgdd_cut) %>% summarise(sgdd = mean(sgdd,
                                                         na.rm = TRUE),
                                                     proba_m = mean(proba,
                                                         na.rm = TRUE),
                                                     proba_ql = quantile(proba,
                                                         probs = 0.05,
                                                         na.rm = TRUE),
                                                     proba_qh = quantile(proba,
                                                         probs = 0.95,
                                                         na.rm = TRUE))
data_wai <- df %>% filter(wai > min(breaks_wai) & wai < max(breaks_wai) ) %>%
                   group_by(wai_cut) %>% summarise(wai = mean(wai,
                                                       na.rm = TRUE),
                                                   proba_m = mean(proba,
                                                       na.rm = TRUE),
                                                   proba_ql = quantile(proba,
                                                       probs = 0.05,
                                                       na.rm = TRUE),
                                                   proba_qh = quantile(proba,
                                                       probs = 0.95,
                                                       na.rm = TRUE))
# Merge with proba
df_sgdd$sp <- spsel
df_sgdd_t <- df_sgdd %>%
    gather(variable, value, -c(sp, BATOTcomp, sgddb, waib, sgdd, wai, sgdd2, wai2)) %>%
    mutate(sgdd = 1/sgddb)
data_sgdd$sgddB <- data_sgdd$sgdd
data_sgdd$sgdd <- sort(round(unique(df_sgdd_t$sgdd)))
df_sgdd_t$sgdd <- round(df_sgdd_t$sgdd)
df_sgdd_tot <- left_join(df_sgdd_t, data_sgdd, by = "sgdd")
df_sgdd_tot$clim_val<- df_sgdd_tot$sgdd
df_sgdd_tot$clim_var<- "sgdd"
df_sgdd_tot <- df_sgdd_tot %>% dplyr::select(-sgddB, -sgdd, -sgdd_cut, -sgddb, -wai, -waib, -sgdd2, -wai2)

df_wai$sp <- spsel
df_wai_t <- df_wai %>%
    gather(variable, value, -c(sp, BATOTcomp, sgddb, waib, sgdd, wai, sgdd2, wai2)) %>%
    mutate(wai = 1/waib-1)
data_wai$waiB <- data_wai$wai
data_wai$wai <- sort(round(unique(df_wai_t$wai), 3))
df_wai_t$wai <- round(df_wai_t$wai, 3)
df_wai_tot <- left_join(df_wai_t, data_wai, by = "wai")
df_wai_tot$clim_val<- df_wai_tot$wai
df_wai_tot$clim_var<- "wai"
df_wai_tot <- df_wai_tot %>% dplyr::select(-waiB, -wai, -wai_cut, -sgddb, -sgdd, -wai, -waib, -sgdd2, -wai2)
df_tot <- rbind(df_sgdd_tot, df_wai_tot)
return(df_tot)
}


extract_proba_per_pca_seq <- function(spsel, data_pred, data, df_pca, clim_pca){
 if (spsel == "Betula") {
   Betulasp <- c("Betula.pendula_Proba", "Betula.pubescens_Proba")
   proba <- apply(data_pred[ , Betulasp], 1, mean)
   names(proba) <-  data_pred$plotcode
 }else{
   if(spsel == "Pinus uncinata"){
    spsel2 <- "Pinus mugo"
   }else{
    spsel2 <- spsel
   }
    proba <- data_pred[,
                       paste0(gsub(" ", ".", spsel2), "_Proba")]
    names(proba) <- data_pred$plotcode
 }
 df_t <- data[!duplicated(data$plotcode), ]
 df <- left_join(df_t,
                 data.frame(plotcode = names(proba),
                            proba = proba,
                            PresAbsObs = data_pred[, gsub(" ", ".", spsel)]),
                 by = "plotcode")
 #pca
 df<- left_join(df, clim_pca$df_PC, by = "plotcode")
 seq_PC1 <- unique(df_pca$PC1)
 step_seq <- seq_PC1[2] - seq_PC1[1]
 seq_PC1_old<- seq_PC1
 seq_PC1 <- c(seq_PC1[1] - (5:1*step_seq), seq_PC1, seq_PC1[30] + (1:5*step_seq))
 mid_bin_PC1<- mean(seq_PC1[-1] - seq_PC1[-length(seq_PC1)])/2
 breaks_PC1 <- seq(from = min(seq_PC1) -mid_bin_PC1,
                 to = max(seq_PC1)+mid_bin_PC1,
                 length.out = length(seq_PC1)+1)
 # get good cut point to match IPM predict
 df <- df %>% mutate(PC1_cut = cut(PC1, breaks = breaks_PC1, label = FALSE))
 data_PC1 <- df %>% filter(PC1 > min(breaks_PC1) &
                           PC1 < max(breaks_PC1) ) %>%
                     group_by(PC1_cut) %>% summarise( proba_m = mean(proba,
                                                          na.rm = TRUE),
                                                      proba_ql = quantile(proba,
                                                          probs = 0.05,
                                                          na.rm = TRUE),
                                                      proba_qh = quantile(proba,
                                                          probs = 0.95,
                                                          na.rm = TRUE),
                                                     proba_obs = mean(PresAbsObs,
                                                          na.rm = TRUE)) %>%
     ungroup()
 data_PC1$PC1 <- seq_PC1[data_PC1$PC1_cut]
 return(data_PC1)
}



format_IPM_metrics_pca<- function(i, ll, spvec, data, data_pred,
                                  clim_pca, data_PC1_mean){
spsel <- spvec[i]
df_pca<- ll[[i]]
df_pca$sgdd2 <- NULL
df_pca$wai2 <- NULL
data_PC1 <- extract_proba_per_pca_seq(spsel, data_pred, data, df_pca, clim_pca)
# merge with proba
df_pca$sp <- spsel
df_pca_t <- left_join(df_pca, data_PC1,  by = "PC1")
df_pca_tot <- df_pca_t %>%
    gather(variable, value, -c(sp, BATOTcomp, PC1, sgddb, waib, sgdd, wai))
df_pca_tot$PC1_mean <- data_PC1_mean$mean_PC1[data_PC1_mean$sp == spsel]
df_pca_tot$PC1_median <- data_PC1_mean$median_PC1[data_PC1_mean$sp == spsel]
df_pca_tot$PC1_ql <- data_PC1_mean$ql_PC1[data_PC1_mean$sp == spsel]
df_pca_tot$PC1_qh <- data_PC1_mean$qh_PC1[data_PC1_mean$sp == spsel]
df_pca_tot <- df_pca_tot %>%
    mutate(PC1_std = (PC1- PC1_mean)/(max(PC1) - min(PC1)),
           sp = spsel)
print(spsel)
return(df_pca_tot)
}

# Resample pca
format_IPM_metrics_pca_Resample <- function(i, ll, spvec, data, data_pred,
                                  clim_pca, data_PC1_mean){
spsel <- spvec[i]
df_pca<- ll[[i]]
data_PC1 <- extract_proba_per_pca_seq(spsel, data_pred, data, df_pca, clim_pca)
# merge with proba
df_pca$sp <- spsel
#df_pca$PC1 <- c(seq_PC1_old, seq_PC1_old)
df_pca_t <- left_join(data_PC1, df_pca, by = "PC1")
df_pca_tot <- df_pca_t %>% dplyr::select( -PC1_cut, -clim_code) %>%
                   gather(variable, value, -c(sp, BATOTcomp, PC1,
                                              sgddb, waib))
df_pca_tot$PC1_mean <- data_PC1_mean$mean_PC1[data_PC1_mean$sp == spsel]
df_pca_tot$PC1_median <- data_PC1_mean$median_PC1[data_PC1_mean$sp == spsel]
df_pca_tot$PC1_ql <- data_PC1_mean$ql_PC1[data_PC1_mean$sp == spsel]
df_pca_tot$PC1_qh <- data_PC1_mean$qh_PC1[data_PC1_mean$sp == spsel]
df_pca_tot <- df_pca_tot %>%
    mutate(PC1_std = (PC1- PC1_mean)/(max(PC1) - min(PC1)),
           sp = spsel)
print(spsel)
return(df_pca_tot)
}


###############################################
##### EXTRACT RESAMPLE SUMMARY

extract_resample_summary_2Vars<- function(res_t){
 res <- res_t %>% dplyr::group_by(clim_code) %>%
                dplyr::summarise(sgddb = mean(sgddb, na.rm = TRUE),
                                 waib = mean(waib, na.rm = TRUE),
                                 BATOTcomp = mean(BATOTcomp, na.rm = TRUE),
                                 lifespan.mean_lifespan_m= mean(lifespan.mean_lifespan, na.rm = TRUE),
                                 passageTime_100_m= mean(passageTime_100, na.rm = TRUE),
                                 SizePassageTime_m= mean(SizePassageTime, na.rm = TRUE),
                                 sizeDeath.mean.size.death_m = mean(sizeDeath.mean.size.death, na.rm = TRUE),
                                 g_m_150_m = mean(g_m_150, na.rm = TRUE),
                                 g_m_500_m = mean(g_m_500, na.rm = TRUE),
                                 s_m_150_m = mean(s_m_150, na.rm = TRUE),
                                 s_m_500_m = mean(s_m_500, na.rm = TRUE),
                                 lifespan.mean_lifespan_ql= quantile(lifespan.mean_lifespan, probs = 0.05,
                                                                     na.rm = TRUE),
                                 passageTime_100_ql = quantile(passageTime_100, probs = 0.05, na.rm = TRUE),
                                 SizePassageTime_ql = quantile(SizePassageTime, probs = 0.05, na.rm = TRUE),
                                 sizeDeath.mean.size.death_ql = quantile(sizeDeath.mean.size.death,
                                                                         probs = 0.05, na.rm = TRUE),
                                 g_m_150_ql = quantile(g_m_150, probs = 0.05, na.rm = TRUE),
                                 g_m_500_ql = quantile(g_m_500, probs = 0.05, na.rm = TRUE),
                                 s_m_150_ql = quantile(s_m_150, probs = 0.05, na.rm = TRUE),
                                 s_m_500_ql = quantile(s_m_500, probs = 0.05, na.rm = TRUE),
                                 lifespan.mean_lifespan_qh= quantile(lifespan.mean_lifespan,
                                                                     probs = 0.95, na.rm = TRUE),
                                 passageTime_100_qh = quantile(passageTime_100, probs = 0.95, na.rm = TRUE),
                                 SizePassageTime_qh = quantile(SizePassageTime, probs = 0.95, na.rm = TRUE),
                                 sizeDeath.mean.size.death_qh = quantile(sizeDeath.mean.size.death,
                                                                         probs = 0.95, na.rm = TRUE),
                                 g_m_150_qh = quantile(g_m_150, probs = 0.95, na.rm = TRUE),
                                 g_m_500_qh = quantile(g_m_500, probs = 0.95, na.rm = TRUE),
                                 s_m_150_qh = quantile(s_m_150, probs = 0.95, na.rm = TRUE),
                                 s_m_500_qh = quantile(s_m_500, probs = 0.95, na.rm = TRUE))

res
}

extract_resample_summary_pca<- function(res_t){
 res <- res_t %>% dplyr::group_by(clim_code) %>%
                dplyr::summarise(sgddb = mean(sgddb, na.rm = TRUE),
                                 waib = mean(waib, na.rm = TRUE),
                                 BATOTcomp = mean(BATOTcomp, na.rm = TRUE),
                                 lifespan.mean_lifespan_m= mean(lifespan.mean_lifespan, na.rm = TRUE),
                                 passageTime_100_m= mean(passageTime_100, na.rm = TRUE),
                                 SizePassageTime_m= mean(SizePassageTime, na.rm = TRUE),
                                 sizeDeath.mean.size.death_m = mean(sizeDeath.mean.size.death, na.rm = TRUE),
                                 g_m_150_m = mean(g_m_150, na.rm = TRUE),
                                 g_m_500_m = mean(g_m_500, na.rm = TRUE),
                                 s_m_150_m = mean(s_m_150, na.rm = TRUE),
                                 s_m_500_m = mean(s_m_500, na.rm = TRUE),
                                 lifespan.mean_lifespan_ql= quantile(lifespan.mean_lifespan, probs = 0.05,
                                                                     na.rm = TRUE),
                                 passageTime_100_ql = quantile(passageTime_100, probs = 0.05, na.rm = TRUE),
                                 SizePassageTime_ql = quantile(SizePassageTime, probs = 0.05, na.rm = TRUE),
                                 sizeDeath.mean.size.death_ql = quantile(sizeDeath.mean.size.death,
                                                                         probs = 0.05, na.rm = TRUE),
                                 g_m_150_ql = quantile(g_m_150, probs = 0.05, na.rm = TRUE),
                                 g_m_500_ql = quantile(g_m_500, probs = 0.05, na.rm = TRUE),
                                 s_m_150_ql = quantile(s_m_150, probs = 0.05, na.rm = TRUE),
                                 s_m_500_ql = quantile(s_m_500, probs = 0.05, na.rm = TRUE),
                                 lifespan.mean_lifespan_qh= quantile(lifespan.mean_lifespan,
                                                                     probs = 0.95, na.rm = TRUE),
                                 passageTime_100_qh = quantile(passageTime_100, probs = 0.95, na.rm = TRUE),
                                 SizePassageTime_qh = quantile(SizePassageTime, probs = 0.95, na.rm = TRUE),
                                 sizeDeath.mean.size.death_qh = quantile(sizeDeath.mean.size.death,
                                                                         probs = 0.95, na.rm = TRUE),
                                 g_m_150_qh = quantile(g_m_150, probs = 0.95, na.rm = TRUE),
                                 g_m_500_qh = quantile(g_m_500, probs = 0.95, na.rm = TRUE),
                                 s_m_150_qh = quantile(s_m_150, probs = 0.95, na.rm = TRUE),
                                 s_m_500_qh = quantile(s_m_500, probs = 0.95, na.rm = TRUE),
                                 PC1 = mean(PC1))
res
}

process_resample_2Vars<-  function(name){
 list_resample <- readRDS(name)
 seq_i <- seq_len(length.out = length(list_resample$sgdd))
 f <- function(i, ll){
  ll[[i]]$resample <- i
  ll[[i]]$clim_code<- seq_len(length.out = nrow(ll[[i]]))
  ll[[i]]
 }
 res_sgdd_l <- lapply(seq_i, f, list_resample$sgdd)
 res_wai_l <- lapply(seq_i, f, list_resample$wai)
 res_sgdd_t <- dplyr::bind_rows(res_sgdd_l)
 res_wai_t <- dplyr::bind_rows(res_wai_l)
 res_sgdd <- extract_resample_summary_2Vars(res_sgdd_t)
 res_wai <- extract_resample_summary_2Vars(res_wai_t)
 return(list(res_sgdd, res_wai))
}

process_resample_pca<-  function(name){
 list_resample <- readRDS(name)
 seq_i <- seq_len(length.out = length(list_resample))
 f <- function(i, ll){
  ll[[i]]$resample <- i
  ll[[i]]$clim_code<- seq_len(length.out = nrow(ll[[i]]))
  ll[[i]]
 }
 res_l <- lapply(seq_i, f, list_resample)
 res_t <- dplyr::bind_rows(res_l)
 res <- extract_resample_summary_pca(res_t)

 res
}


process_IPM_e_metrics_all_sp_pca_Resample<- function(ll_name,
                                               name = "BaseA"){
 ll <- readRDS(ll_name)
 l <- lapply(ll, process_resample_pca)
 saveRDS(l, file = file.path("output", paste0("MetricsPCAVarResampleProcess",
                name, ".rds")))
}


process_IPM_e_metrics_all_sp_2Vars_Resample<- function(ll_name,
                                               name = "BaseA"){
 ll <- readRDS(ll_name)
 l <- lapply(ll, process_resample_2Vars)
 saveRDS(l, file = file.path("output", paste0("Metrics2VarsResampleProcess",
                name, ".rds")))
}


## MERGE to keep all 100 resample and do stat

merge_resample <-  function(name){
 list_resample <- readRDS(name)
 seq_i <- seq_len(length.out = length(list_resample))
 f <- function(i, ll){
  ll[[i]]$resample <- i
  ll[[i]]$clim_code<- seq_len(length.out = nrow(ll[[i]]))
  ll[[i]]
 }
 res_l <- lapply(seq_i, f, list_resample)
 res_t <- dplyr::bind_rows(res_l)
return(res_t)
}


format_IPM_metrics_pca_each_resample<- function(i, ll, spvec, data, data_pred,
                                  clim_pca, data_PC1_mean){
spsel <- spvec[i]
df_pca<- ll[[i]]
df_pca$sgdd2 <- NULL
df_pca$wai2 <- NULL
df_pca$clim_code <- c(1:30, 1:30)[df_pca$clim_code]
data_PC1 <- extract_proba_per_pca_seq(spsel, data_pred, data, df_pca, clim_pca)
# merge with proba
df_pca$sp <- spsel
df_pca_t <- left_join(df_pca, data_PC1,  by = "PC1")
df_pca_t$PC1_mean <- data_PC1_mean$mean_PC1[data_PC1_mean$sp == spsel]
df_pca_t$PC1_median <- data_PC1_mean$median_PC1[data_PC1_mean$sp == spsel]
df_pca_t$PC1_ql <- data_PC1_mean$ql_PC1[data_PC1_mean$sp == spsel]
df_pca_t$PC1_qh <- data_PC1_mean$qh_PC1[data_PC1_mean$sp == spsel]
# reshape format
df_pca_tot <- df_pca_t %>%
    gather(variable, value, -c(sp, BATOTcomp, PC1, sgddb, waib, sgdd, wai,
                               resample, clim_code, PC1_mean,
                               PC1_median, PC1_ql, PC1_qh))

print(spsel)
return(df_pca_tot)
}


merge_IPM_e_metrics_all_sp_pca_Resample <- function(ll_name, name = "Best2", sps, data,
                                                    data_pred, clim_pca, data_PC1_mean){
 ll <- readRDS(ll_name)
 l <- lapply(ll, merge_resample)
 seq_i <- seq_len(length.out = length(l))
 l2 <- lapply(seq_i, format_IPM_metrics_pca_each_resample, ll = l,
                  spvec = sps$sp, data = data,
                  data_pred = data_pred, clim_pca = clim_pca,
                  data_PC1_mean = data_PC1_mean)

 res <- dplyr::bind_rows(l2)
 saveRDS(res, file = file.path("output", paste0("MetricsPCAVarResampleMerge",
                name, ".rds")))
}

#### MERGE EDGES

format_IPM_metrics_edge_pca_each_resample<- function(i, ll, spvec, data, data_pred,
                                  clim_pca, data_PC1_mean){
spsel <- spvec[i]
df_pca<- ll[[i]]
df_pca$sgdd2 <- NULL
df_pca$wai2 <- NULL
df_pca$clim_code <- c(1:30, 1:30)[df_pca$clim_code]
data_PC1 <- extract_proba_per_pca_seq(spsel, data_pred, data, df_pca, clim_pca)
# merge with proba
df_pca$sp <- spsel
df_pca_t <- left_join(df_pca, data_PC1,  by = "PC1")
df_pca_t$PC1_mean <- data_PC1_mean$mean_PC1[data_PC1_mean$sp == spsel]
df_pca_t$PC1_median <- data_PC1_mean$median_PC1[data_PC1_mean$sp == spsel]
df_pca_t$PC1_ql <- data_PC1_mean$ql_PC1[data_PC1_mean$sp == spsel]
df_pca_t$PC1_qh <- data_PC1_mean$qh_PC1[data_PC1_mean$sp == spsel]
# get 2 value closest to median
dd <- df_pca_t %>% mutate(dist = abs(PC1 - PC1_median) ,
                          dist %in% sort(unique(dist))[1:2])
median_index <- unique(dd$clim_code[dd$dist %in% sort(unique(dd$dist))[1:2]])
low_index <- 1:2
high_index <- 29:30
# reshape format
df_pca_tot <- df_pca_t %>%
    gather(variable, value, -c(sp, BATOTcomp, PC1, sgddb, waib, sgdd, wai,
                               resample, clim_code, PC1_mean,
                               PC1_median, PC1_ql, PC1_qh))

# keep only low median and high
df_pca_tot <- df_pca_tot %>% filter(clim_code %in% c(low_index, median_index, high_index)) %>%
                                    mutate(clim_code2= case_when(clim_code %in% low_index~ "low",
                                                                 clim_code %in% median_index~ "median",
                                                                 clim_code %in% high_index~ "high"))
# Compute mean per clim_code2 and BATOTcomp
df_pca_extract <- df_pca_tot %>% group_by(clim_code2, BATOTcomp, resample, variable) %>%
      summarise(sp = unique(sp),
                PC1 = mean(PC1),
                sgddb = mean(sgddb),
                sgdd = mean(sgdd),
                waib = mean(waib),
                wai = mean(wai),
                value = mean(value)) %>% ungroup()
print(spsel)
return(df_pca_extract)
}


merge_IPM_e_metrics_edge_all_sp_pca_Resample <- function(ll_name, name = "Best2", sps, data,
                                                    data_pred, clim_pca, data_PC1_mean){
 ll <- readRDS(ll_name)
 l <- lapply(ll, merge_resample)
 seq_i <- seq_len(length.out = length(l))
 l2 <- lapply(seq_i, format_IPM_metrics_edge_pca_each_resample, ll = l,
                  spvec = sps$sp, data = data,
                  data_pred = data_pred, clim_pca = clim_pca,
                  data_PC1_mean = data_PC1_mean)

 res <- dplyr::bind_rows(l2)
 df2 <- res %>% filter(clim_code2 == "median") %>% arrange(PC1)
 res$sp <- as.factor(res$sp)
 res$sp <- factor(res$sp, levels = unique(df2$sp))

 saveRDS(res, file = file.path("output", paste0("MetricsPCAVarResampleMergeEdge",
                name, ".rds")))
}




#######################
### EXTRACT PARAMS

extract_params <- function(spsel= "Fagus sylvatica",
                            list_fit_sg){
fit_sg <- readRDS(list_fit_sg[[spsel]])
params_sv <- fit_sg$sv$params
names(params_sv) <- paste0("sv_", names(params_sv))
params_gr <- fit_sg$gr$params
names(params_gr) <- paste0("gr_", names(params_gr))
return(data.frame(sp = spsel, as.list(params_sv), as.list(params_gr)))
}

extract_params_all_sp<- function(sps, list_fit_sg){
 l <- lapply(sps$sp,extract_params , list_fit_sg= list_fit_sg)
 df <- rbind_list(l)
 print(names(df))
 return(df)
}


#################################
### SPECIES DISTRIBUTION ALONG PC1
## get species dist along CLIMATE PCA axis one
format_data_mean_metrics <- function(list_mean, sps){
 names(list_mean) <- sps$sp
 f <- function(spsel, list_m) {
  df <- list_m[[spsel]]
  df$sp <- spsel
  return(df)
 }
 ll <- lapply(sps$sp,f, list_m = list_mean)
 df <- do.call("rbind", ll)
 return(df)
}

## PC1 Dist
make_PC1_mean_sp<- function(sps, data, data_pred, clim_pca){
 l <- lapply(sps$sp, get_mean_sp_pca_i, data = data, data_pred = data_pred,
             clim_pca = clim_pca)
 res <- do.call("rbind", l)
 return(res)
}

#####################################################
#####################################################
## EXTRACT METRICS A THE TWO EDGES AND CENTRE OF DISTRIBUTION

extract_summary_IPM_metrics_all_sp_pca <- function(df){
    df <- df[df$sp != "Juniperus thurifera", ]
    df$value[!is.na(df$value) & df$value < 1.1 &
             df$variable == "passageTime_100"] <-  NA
    df$value[!is.na(df$value) & df$value == 0] <-  NA
    df$value[df$variable == "g_m_150"] <-
        exp(df$value[df$variable == "g_m_150"])
    df$value[df$variable == "g_m_500"] <-
        exp(df$value[df$variable == "g_m_500"])
    df <- df %>% filter(!variable %in% c("lifespan.variance_lifespan",
                                         "lifespan.CV_lifespan",
                                         "SizePassageTime",
                                         "sizeDeath.var.size.death",
                                         "proba_qh", "proba_ql", "PC1_cut",
                                         "proba_obs", "s_m_500", "g_m_500"))
    df$PC1_median_cut <- cut(df$PC1_median, 4)
    df$variable <- factor(df$variable, levels = c("g_m_150", "s_m_150",
                                                  "lifespan.mean_lifespan",
                                                  "passageTime_100",
                                                  "sizeDeath.mean.size.death",
                                                  "proba_m"))
    df$value[df$value <0.02 &
             df$variable %in% c("passageTime_100", "lifespan.mean_lifespan",
                                "s_m_150")] <- NA
    df <- df[df$sp != "Pinus pinea", ]
    df <- df[!is.na(df$value), ]
    df1 <- df[df$BATOTcomp == 0 & !is.na(df$BATOTcomp),] # extract demo at clim with max prob pres
    df_extract1 <- df1 %>% #filter(variable == "proba_m") %>%
        group_by(sp) %>%
        summarise(PC1_max = PC1[variable == "proba_m"][value[variable == "proba_m"] ==
                                                       max(value[variable == "proba_m"])],
                  PC1_median = mean(PC1_median),
                  PC1_l = mean(sort(PC1[variable == "proba_m"])[1:2]),
                  PC1_h = mean(sort(PC1[variable == "proba_m"],
                                    decreasing = TRUE)[1:2]),
                  # extract PC1 position with max proba, median, and quantile low and high
                  g_m_150_max_NoComp = value[variable == "g_m_150"][PC1[variable == "g_m_150"] == PC1_max],
                  g_m_150_med_NoComp = value[variable == "g_m_150"][
                      which.min(abs(PC1[variable == "g_m_150"] - PC1_median)) ],
                  PC1_med = PC1[variable == "proba_m"][which.min(abs(PC1[variable == "proba_m"] - PC1_median)) ],
                  g_m_150_l_NoComp = mean(value[variable == "g_m_150"][PC1[variable == "g_m_150"] %in%
                                                               sort(PC1[variable == "proba_m"])[1:2]],
                                  na.rm = TRUE),
                  g_m_150_h_NoComp = mean(value[variable == "g_m_150"][PC1[variable == "g_m_150"] %in%
                                                               sort(PC1[variable == "proba_m"],
                                                                    decreasing = TRUE)[1:2]],
                                  na.rm = TRUE),
                  s_m_150_max_NoComp = value[variable == "s_m_150"][PC1[variable == "s_m_150"] == PC1_max],
                  s_m_150_med_NoComp = value[variable == "s_m_150"][
                      which.min(abs(PC1[variable == "s_m_150"] - PC1_median)) ],
                  s_m_150_l_NoComp = mean(value[variable == "s_m_150"][PC1[variable == "s_m_150"] %in%
                                                               sort(PC1[variable == "proba_m"])[1:2]],
                                  na.rm = TRUE),
                  s_m_150_h_NoComp = mean(value[variable == "s_m_150"][PC1[variable == "s_m_150"] %in%
                                                               sort(PC1[variable == "proba_m"],
                                                                    decreasing = TRUE)[1:2]],
                                  na.rm = TRUE),
                  lifespan_max_NoComp = value[variable == "lifespan.mean_lifespan"][PC1[variable ==
                                                                                 "lifespan.mean_lifespan"] ==
                                                                             PC1_max],
                  lifespan_med_NoComp = value[variable == "lifespan.mean_lifespan"][
                      which.min(abs(PC1[variable == "lifespan.mean_lifespan"] - PC1_median)) ],
                  lifespan_l_NoComp = mean(value[variable == "lifespan.mean_lifespan"][PC1[variable ==
                                                                                 "lifespan.mean_lifespan"] %in%
                                                                         sort(PC1[variable == "proba_m"])[1:2]],
                                    na.rm = TRUE),
                  lifespan_h_NoComp = mean(value[variable == "lifespan.mean_lifespan"][PC1[variable ==
                                                                                "lifespan.mean_lifespan"] %in%
                                                                               sort(PC1[variable == "proba_m"],
                                                                                     decreasing = TRUE)[1:2]],
                                    na.rm = TRUE),
                  passageTime_max_NoComp = value[variable == "passageTime_100"]
                        [which.min(abs(PC1[variable == "passageTime_100"] - PC1_max))],
                  passageTime_med_NoComp = value[variable == "passageTime_100"][
                      which.min(abs(PC1[variable == "passageTime_100"] - PC1_median)) ],
                  passageTime_l_NoComp = mean(value[variable == "passageTime_100"][PC1[variable ==
                                                                                "passageTime_100"] %in%
                                                                         sort(PC1[variable == "proba_m"])[1:2]],
                                       na.rm = TRUE),
                  passageTime_h_NoComp = mean(value[variable == "passageTime_100"][PC1[variable ==
                                                                                "passageTime_100"] %in%
                                                                            sort(PC1[variable == "proba_m"],
                                                                                 decreasing = TRUE)[1:2]],
                                       na.rm = TRUE),
                  sizeDeath_max_NoComp = value[variable == "sizeDeath.mean.size.death"][PC1[variable ==
                                                                            "sizeDeath.mean.size.death"] ==
                                                                                 PC1_max],
                  sizeDeath_med_NoComp = value[variable == "sizeDeath.mean.size.death"][
                      which.min(abs(PC1[variable == "sizeDeath.mean.size.death"] - PC1_median)) ],
                  sizeDeath_l_NoComp = mean(value[variable == "sizeDeath.mean.size.death"][PC1[variable ==
                                                                              "sizeDeath.mean.size.death"] %in%
                                                                        sort(PC1[variable == "proba_m"])[1:2]],
                                     na.rm = TRUE),
                  sizeDeath_h_NoComp = mean(value[variable == "sizeDeath.mean.size.death"][PC1[variable ==
                                                                            "sizeDeath.mean.size.death"] %in%
                                                                            sort(PC1[variable == "proba_m"],
                                                                                      decreasing = TRUE)[1:2]],
                                     na.rm = TRUE),
                  PC1_median = mean(PC1_median, na.rm = TRUE)
                  ) %>% ungroup()
    df2 <- df[df$BATOTcomp > 0 & !is.na(df$BATOTcomp),] # extract demo at clim with max prob pres
    df_extract2 <- df2 %>% #filter(variable == "proba_m") %>%
        group_by(sp) %>%
        summarise(PC1_max = PC1[variable == "proba_m"][value[variable == "proba_m"] ==
                                                       max(value[variable == "proba_m"])],
                  PC1_median = mean(PC1_median),
                  PC1_l = mean(sort(PC1[variable == "proba_m"])[1:2]),
                  PC1_h = mean(sort(PC1[variable == "proba_m"],
                                    decreasing = TRUE)[1:2]),
                  g_m_150_max_Comp = value[variable == "g_m_150"][PC1[variable == "g_m_150"] == PC1_max],
                  g_m_150_med_Comp = value[variable == "g_m_150"][
                      which.min(abs(PC1[variable == "g_m_150"] - PC1_median)) ],
                  PC1_med = PC1[variable == "proba_m"][which.min(abs(PC1[variable == "proba_m"] - PC1_median)) ],
                  g_m_150_l_Comp = mean(value[variable == "g_m_150"][PC1[variable == "g_m_150"] %in%
                                                               sort(PC1[variable == "proba_m"])[1:2]],
                                  na.rm = TRUE),
                  g_m_150_h_Comp = mean(value[variable == "g_m_150"][PC1[variable == "g_m_150"] %in%
                                                               sort(PC1[variable == "proba_m"],
                                                                    decreasing = TRUE)[1:2]],
                                  na.rm = TRUE),
                  s_m_150_max_Comp = value[variable == "s_m_150"][PC1[variable == "s_m_150"] == PC1_max],
                  s_m_150_med_Comp = value[variable == "s_m_150"][
                      which.min(abs(PC1[variable == "s_m_150"] - PC1_median)) ],
                  s_m_150_l_Comp = mean(value[variable == "s_m_150"][PC1[variable == "s_m_150"] %in%
                                                               sort(PC1[variable == "proba_m"])[1:2]],
                                  na.rm = TRUE),
                  s_m_150_h_Comp = mean(value[variable == "s_m_150"][PC1[variable == "s_m_150"] %in%
                                                               sort(PC1[variable == "proba_m"],
                                                                    decreasing = TRUE)[1:2]],
                                  na.rm = TRUE),
                  lifespan_max_Comp = value[variable == "lifespan.mean_lifespan"][PC1[variable ==
                                                                                 "lifespan.mean_lifespan"] ==
                                                                             PC1_max],
                  lifespan_med_Comp = value[variable == "lifespan.mean_lifespan"][
                      which.min(abs(PC1[variable == "lifespan.mean_lifespan"] - PC1_median)) ],
                  lifespan_l_Comp = mean(value[variable == "lifespan.mean_lifespan"][PC1[variable ==
                                                                                 "lifespan.mean_lifespan"] %in%
                                                                         sort(PC1[variable == "proba_m"])[1:2]],
                                    na.rm = TRUE),
                  lifespan_h_Comp = mean(value[variable == "lifespan.mean_lifespan"][PC1[variable ==
                                                                                "lifespan.mean_lifespan"] %in%
                                                                               sort(PC1[variable == "proba_m"],
                                                                                     decreasing = TRUE)[1:2]],
                                    na.rm = TRUE),
                  passageTime_max_Comp = value[variable == "passageTime_100"][
                                                   which.min(abs(PC1[variable == "passageTime_100"] -
                                                                 PC1_max))],
                  passageTime_med_Comp = value[variable == "passageTime_100"][
                      which.min(abs(PC1[variable == "passageTime_100"] - PC1_median)) ],
                  passageTime_l_Comp = mean(value[variable == "passageTime_100"][PC1[variable ==
                                                                                "passageTime_100"] %in%
                                                                         sort(PC1[variable == "proba_m"])[1:2]],
                                       na.rm = TRUE),
                  passageTime_h_Comp = mean(value[variable == "passageTime_100"][PC1[variable ==
                                                                                "passageTime_100"] %in%
                                                                            sort(PC1[variable == "proba_m"],
                                                                                 decreasing = TRUE)[1:2]],
                                       na.rm = TRUE),
                  sizeDeath_max_Comp = value[variable == "sizeDeath.mean.size.death"][PC1[variable ==
                                                                            "sizeDeath.mean.size.death"] ==
                                                                                 PC1_max],
                  sizeDeath_med_Comp = value[variable == "sizeDeath.mean.size.death"][
                      which.min(abs(PC1[variable == "sizeDeath.mean.size.death"] - PC1_median)) ],
                  sizeDeath_l_Comp = mean(value[variable == "sizeDeath.mean.size.death"][PC1[variable ==
                                                                              "sizeDeath.mean.size.death"] %in%
                                                                        sort(PC1[variable == "proba_m"])[1:2]],
                                     na.rm = TRUE),
                  sizeDeath_h_Comp = mean(value[variable == "sizeDeath.mean.size.death"][PC1[variable ==
                                                                            "sizeDeath.mean.size.death"] %in%
                                                                            sort(PC1[variable == "proba_m"],
                                                                                      decreasing = TRUE)[1:2]],
                                     na.rm = TRUE),
                  PC1_median = mean(PC1_median, na.rm = TRUE)
                  ) %>% ungroup()%>% select( -PC1_max, -PC1_l, -PC1_h, -PC1_median)
    df_extract <- left_join(df_extract1, df_extract2, by = "sp")

    df_ratio<- df_extract %>% mutate(g_m_150_lp_NoComp= -( g_m_150_l_NoComp - g_m_150_max_NoComp)/g_m_150_max_NoComp,
                                 g_m_150_hp_NoComp = -(g_m_150_h_NoComp-g_m_150_max_NoComp)/g_m_150_max_NoComp,
                                 s_m_150_lp_NoComp = -(s_m_150_l_NoComp-s_m_150_max_NoComp)/s_m_150_max_NoComp,
                                 s_m_150_hp_NoComp = -(s_m_150_h_NoComp-s_m_150_max_NoComp)/s_m_150_max_NoComp,
                                 lifespan_lp_NoComp = -(lifespan_l_NoComp-lifespan_max_NoComp)/lifespan_max_NoComp,
                                 lifespan_hp_NoComp = -(lifespan_h_NoComp-lifespan_max_NoComp)/lifespan_max_NoComp,
                                 passageTime_lp_NoComp = -(passageTime_l_NoComp-passageTime_max_NoComp)/passageTime_max_NoComp,
                                 passageTime_hp_NoComp = -(passageTime_h_NoComp-passageTime_max_NoComp)/passageTime_max_NoComp,
                                 sizeDeath_lp_NoComp = -(sizeDeath_l_NoComp-sizeDeath_max_NoComp)/sizeDeath_max_NoComp,
                                 sizeDeath_hp_NoComp = -(sizeDeath_h_NoComp-sizeDeath_max_NoComp)/sizeDeath_max_NoComp
                                 ) %>%
        dplyr::select(-PC1_max, -PC1_l, -PC1_h,
               -sizeDeath_l_NoComp, -sizeDeath_h_NoComp, -sizeDeath_max_NoComp,
               -passageTime_l_NoComp, -passageTime_h_NoComp,-passageTime_max_NoComp,
               -lifespan_l_NoComp, -lifespan_h_NoComp, -lifespan_max_NoComp,
               -s_m_150_l_NoComp, -s_m_150_h_NoComp, -s_m_150_max_NoComp,
               -g_m_150_l_NoComp, -g_m_150_h_NoComp, -g_m_150_max_NoComp)

return(list(df_ratio, df_extract))
}



#####################################################################
### FUNCTION TO SELECT EDGE WITH A DECREASE OF PRESENCE ABSENCE
#####################################################################

select_edges <- function(name){
df <- readRDS(name)
dd <- df %>% filter(variable == "proba_obs") %>%
        group_by(sp, clim_code2) %>%
        summarise(value = mean(value),
                  PC1 = mean(PC1)) %>%
        gather(variable, value2,  -(sp:clim_code2)) %>%
        unite(variable2, c("clim_code2", "variable")) %>%
        spread(variable2, value2)
df2 <- df %>% filter(clim_code2 == "median") %>% arrange(PC1)
dd$sp <- as.factor(dd$sp)
dd$sp <- factor(dd$sp, levels = unique(df2$sp))
res <- data.frame(sp = dd$sp,
                  select_low_TF = dd$low_value/dd$median_value < 0.9,
                  select_high_TF = dd$high_value/dd$median_value < 0.9)#,
#                  ratio_low = dd$low_value/dd$median_value,
#                  ratio_high = dd$high_value/dd$median_value)
df_m <- left_join(df, res, by = "sp")
return(df_m)
}


##############################
#############################
## AOV of edges demo metrics
fun_fit_aov_sp <- function(spsel, df){
 res1 <- aov(value ~ edge_type, df[df$sp == spsel, ])
 p_value <- anova(res1)$"Pr(>F)"[1]
 p_value
}

aov_demo_edge_type <- function(df, var = "s_m_150", edge = "low"){
if(!edge %in% c("low", "high")) stop("error edge should be either low or high")
df <- df %>% filter(BATOTcomp == 0 ) %>% rename(edge_type = clim_code2)
if(var %in% c("g_m_150", "g_m_500")) df$value <- exp(df$value)
if (edge == "low"){
df_t <- df %>% filter(variable == var & select_low_TF == TRUE & edge_type %in% c("low", "median"))
}
if (edge == "high"){
df_t <- df %>% filter(variable == var & select_high_TF == TRUE & edge_type %in% c("high", "median"))
}
#aov per species
l_res <- sapply(unique(df_t$sp), FUN =fun_fit_aov_sp, df_t)
names(l_res) <- unique(df_t$sp)
#aov overall species with random species effect
res2 <- lmerTest::lmer(value ~ edge_type +(1+edge_type|sp), df_t)
random_res <- anova(res2)$"Pr(>F)"
pvalues <- list(species = l_res, random = random_res)
return(pvalues)
}

# loop over all variables and edges
aov_demo_edge_type_all <- function(df, vars = c("g_m_150",
                                                  "lifespan.mean_lifespan",
                                                  "passageTime_100",
                                                  "s_m_150")){
list_sp <- vector("list")
list_r <- vector("list")
 for (var in vars){
   list_res_l<- aov_demo_edge_type(df, var, edge= "low")
   list_res_h<- aov_demo_edge_type(df, var, edge= "high")
   list_sp[[var]] <- list(low = list_res_l$species, high = list_res_h$species)
   list_r[[var]] <- list(low = list_res_l$random, high = list_res_h$random)
 }
return(list(species = list_sp, random = list_r))
}

## INCLUDE THE TWO LEVELS OF COMPETITIONS

aov_demo_edge_type_BATOTcomp<- function(df, var = "s_m_150", edge = "low"){
if(!edge %in% c("low", "high")) stop("error edge should be either low or high")
df <- df %>% rename(edge_type = clim_code2)
if(var %in% c("g_m_150", "g_m_500")) df$value <- exp(df$value)
if (edge == "low"){
df_t <- df %>% filter(variable == var & select_low_TF == TRUE & edge_type %in% c("low", "median"))
}
if (edge == "high"){
df_t <- df %>% filter(variable == var & select_high_TF == TRUE & edge_type %in% c("high", "median"))
}
#aov per species
l_res_l<- sapply(unique(df_t$sp), FUN =fun_fit_aov_sp, df_t[df_t$BATOTcomp == 0 , ])
names(l_res_l) <- unique(df_t$sp)
l_res_h<- sapply(unique(df_t$sp), FUN =fun_fit_aov_sp, df_t[df_t$BATOTcomp == 30 , ])
names(l_res_h) <- unique(df_t$sp)
#aov overall species with random species effect
res2_l<- lmerTest::lmer(value ~ edge_type +(1+edge_type|sp), df_t[df_t$BATOTcomp == 0 , ])
random_res_l<- anova(res2_l)$"Pr(>F)"
res2_h<- lmerTest::lmer(value ~ edge_type +(1+edge_type|sp), df_t[df_t$BATOTcomp == 30 , ])
random_res_h<- anova(res2_h)$"Pr(>F)"
pvalues <- list(species = list(Comp0 = l_res_l, Comp30 = l_res_h),
                random = list(Comp0 = random_res_l, Comp30 = random_res_h))
return(pvalues)
}

# loop over all variables and edges
aov_demo_edge_type_BATOTcomp_all <- function(df, vars = c("g_m_150",
                                                  "lifespan.mean_lifespan",
                                                  "passageTime_100",
                                                  "s_m_150")){
list_sp <- vector("list")
list_r <- vector("list")
 for (var in vars){
   list_res_l<- aov_demo_edge_type_BATOTcomp(df, var, edge= "low")
   list_res_h<- aov_demo_edge_type_BATOTcomp(df, var, edge= "high")
   list_sp[[var]] <- list(low = list_res_l$species, high = list_res_h$species)
   list_r[[var]] <- list(low = list_res_l$random, high = list_res_h$random)
 }
return(list(species = list_sp, random = list_r))
}



# format species pvalues table
format_table_pvalues_sp <- function(ll_res){
f <- function(var, ll){
 llt <- ll[[var]]
 df_l <- data.frame(sp = names(llt[["low"]]), edge_type = "low", p_val= llt[["low"]])
 df_h <- data.frame(sp = names(llt[["high"]]), edge_type = "high", p_val= llt[["high"]])
 df <- bind_rows(df_l, df_h)
 df$variable <- var
 df
}
df_pvals <- bind_rows(lapply(names(ll_res$species), f, ll_res$species))
count_sp_signif <- df_pvals %>% group_by(variable,edge_type) %>%
    summarise(n_signif_sp = sum(p_val<0.05, na.rm = TRUE)/ sum(!is.na(p_val)))
return(list(df_pvals, count_sp_signif))
}


#### TEST LINK BETWEEN EDGE DIFF AND PC1
pval_demo_diff_per_PC1<- function(df,
                                    vars = c("g_m_150", "s_m_150",
                                             "passageTime_100",
                                             "lifespan.mean_lifespan")){
df <- df %>% filter(variable %in% vars) %>% rename(edge_type = clim_code2)
df$value[df$variable %in% c("g_m_150", "g_m_500")] <- exp(df$value[df$variable %in% c("g_m_150", "g_m_500")])
df_low <- df %>% filter(select_low_TF == TRUE & edge_type%in% c("low", "median") & BATOTcomp == 0) %>%
    group_by(sp, resample, variable) %>%
    summarise(rel_demo_change = (value[edge_type == "low"] -
                                  value[edge_type == "median"])/value[edge_type == "median"],
              ratio_demo_change = value[edge_type == "low"]/value[edge_type == "median"],
              PC1 = mean(PC1[edge_type == "median"])) %>%
    mutate(edge_type = "low") %>% ungroup()
df_high <- df %>% filter(select_high_TF == TRUE & edge_type%in% c("high", "median") & BATOTcomp == 0) %>%
    group_by(sp, resample,variable) %>%
    summarise(rel_demo_change = (value[edge_type == "high"] -
                                  value[edge_type == "median"])/value[edge_type == "median"],
              ratio_demo_change = value[edge_type == "high"]/value[edge_type == "median"],
              PC1 = mean(PC1[edge_type == "median"])) %>%
    mutate(edge_type = "high") %>% ungroup()
df_p <- rbind(df_low, df_high)
df_p$edge_type <- factor(df_p$edge_type)
df_p$edge_type <- factor(df_p$edge_type,levels(df_p$edge_type)[c(2,1)])
df_p$variable <- factor(df_p$variable)
df_p$variable <- factor(df_p$variable,levels(df_p$variable)[c(1,4,3,2)])
##################################
### compute regression lines

# format data
df_var_mean<- df_p %>% group_by(sp, edge_type, variable) %>%
    summarise(var = var(rel_demo_change, na.rm = TRUE),
              mean = mean(rel_demo_change, na.rm = TRUE),
              PC1 = mean(PC1, na.rm = TRUE))
df_facet <- expand.grid(variable = unique(df_var_mean$variable),
                        edge_type = unique(df_var_mean$edge_type))
# lm with weight based on inverse of variance
p_vals = sapply(seq_len(nrow(df_facet)),
                function(i,df_f, df) {
                 df_t <- df[df$variable == df_f$variable[i] &
                            df$edge_type== df_f$edge_type[i], ]
                 anova(lm(mean ~ PC1, data = df_t,
                          weights = 1/df_t$var))$'Pr(>F)'[1]},
                df_facet, df_var_mean)
df_facet$p_vals <- p_vals
df_facet
}



###############################################
#### TEST LINK BETWEEN EDGE DIFF AND Traits
pval_demo_diff_per_traits <- function(df, df_traits,
                                    traits_sel = c("Nmass",
                                                    "Wood_density",
                                                   "PI50", "Leaf_size_cm2"),
                                    vars = c("g_m_150", "s_m_150",
                                             "passageTime_100",
                                             "lifespan.mean_lifespan"),
                                    name = "open"){
df <- df %>% filter(variable %in% vars) %>% rename(edge_type = clim_code2)
df$value[df$variable %in%
         c("g_m_150", "g_m_500")] <- exp(df$value[df$variable %in%
                                                  c("g_m_150", "g_m_500")])
df_low <- df %>% filter(select_low_TF == TRUE &
                        edge_type%in% c("low", "median") & BATOTcomp == 0) %>%
    group_by(sp, resample, variable) %>%
    summarise(rel_demo_change = (value[edge_type == "low"] -
                                  value[edge_type == "median"])/
              value[edge_type == "median"],
              ratio_demo_change = value[edge_type == "low"]/
              value[edge_type == "median"],
              PC1 = mean(PC1[edge_type == "median"])) %>%
    mutate(edge_type = "low") %>% ungroup()
df_high <- df %>% filter(select_high_TF == TRUE &
                         edge_type%in% c("high", "median") &
                         BATOTcomp == 0) %>%
    group_by(sp, resample,variable) %>%
    summarise(rel_demo_change = (value[edge_type == "high"] -
                                  value[edge_type == "median"])/
              value[edge_type == "median"],
              ratio_demo_change = value[edge_type == "high"]/
              value[edge_type == "median"],
              PC1 = mean(PC1[edge_type == "median"])) %>%
    mutate(edge_type = "high") %>% ungroup()
df_p <- rbind(df_low, df_high)
df_p$edge_type <- factor(df_p$edge_type)
df_p$edge_type <- factor(df_p$edge_type,levels(df_p$edge_type)[c(2,1)])
df_p$variable <- factor(df_p$variable)
df_p$variable <- factor(df_p$variable,levels(df_p$variable)[c(1,4,3,2)])
##################################
### test link of demo diff with traits

# format data
df_var_mean<- df_p %>% group_by(sp, edge_type, variable) %>%
    summarise(var = var(rel_demo_change, na.rm = TRUE),
              mean = mean(rel_demo_change, na.rm = TRUE),
              PC1 = mean(PC1, na.rm = TRUE))
# exclude genus mean
df_merge <- left_join(df_var_mean, df_traits, by = c("sp" = "sp"))
df_facet <- expand.grid(variable = unique(df_var_mean$variable),
                        edge_type = unique(df_var_mean$edge_type))
# lm with weight based on inverse of variance
for (tt in traits_sel){
# exclude genus mean
df_t <- df_merge
df_t[df_t[[paste0(tt, "_GenusM")]] == "YES" &
     !is.na(df_t[[paste0(tt, "_GenusM")]]), tt] <- NA
 p_vals = sapply(seq_len(nrow(df_facet)),
                 function(i,df_f, df) {
                  df_t <- df[df$variable == df_f$variable[i] &
                             df$edge_type== df_f$edge_type[i], ]
                  df_t[[tt]] <- as.vector(scale(df_t[[tt]]))
                  df_t$mean <- as.vector(scale(df_t$mean))
                  ff <- as.formula(paste0("mean ~ ", tt, ""))
                  res <- lm(ff, data = df_t,
                           weights = 1/df_t$var)
                  c(anova(res)$'Pr(>F)'[1], summary(res)$r.squared,
                    coef(res)[2],
                    sjstats::std_beta(res)[1, c("std.estimate",
                                                "conf.low", "conf.high")])},
                 df_facet, df_t)
 res_t<- as.data.frame(t(p_vals))
 names(res_t) <- paste0(tt, c("_pval", "_adjR2", "_coef", "_std.estimate", "_conf.low", "_conf.high"))
 df_facet <- cbind(df_facet, res_t)
}
saveRDS(df_facet, file = file.path("output",
                                   paste0("reg_demo_edge_vs_traits_resample_",
                                          name,".rds")))
}


###############################################
#### TEST LINK BETWEEN EDGE DIFF AND Traits



pval_demo_diff_per_traits_phylo <- function(df, df_traits, tree,
                                    traits_sel = c("Nmass",
                                                    "Wood_density",
                                                   "PI50", "Leaf_size_cm2"),
                                    vars = c("g_m_150", "s_m_150",
                                             "passageTime_100",
                                             "lifespan.mean_lifespan"),
                                    name = "open"){
df <- df %>% filter(variable %in% vars) %>% rename(edge_type = clim_code2)
df$value[df$variable %in%
         c("g_m_150", "g_m_500")] <- exp(df$value[df$variable %in%
                                                  c("g_m_150", "g_m_500")])
df_low <- df %>% filter(select_low_TF == TRUE &
                        edge_type%in% c("low", "median") & BATOTcomp == 0) %>%
    group_by(sp, resample, variable) %>%
    summarise(rel_demo_change = (value[edge_type == "low"] -
                                  value[edge_type == "median"])/
              value[edge_type == "median"],
              ratio_demo_change = value[edge_type == "low"]/
              value[edge_type == "median"],
              PC1 = mean(PC1[edge_type == "median"])) %>%
    mutate(edge_type = "low") %>% ungroup()
df_high <- df %>% filter(select_high_TF == TRUE &
                         edge_type%in% c("high", "median") &
                         BATOTcomp == 0) %>%
    group_by(sp, resample,variable) %>%
    summarise(rel_demo_change = (value[edge_type == "high"] -
                                  value[edge_type == "median"])/
              value[edge_type == "median"],
              ratio_demo_change = value[edge_type == "high"]/
              value[edge_type == "median"],
              PC1 = mean(PC1[edge_type == "median"])) %>%
    mutate(edge_type = "high") %>% ungroup()
df_p <- rbind(df_low, df_high)
df_p$edge_type <- factor(df_p$edge_type)
df_p$edge_type <- factor(df_p$edge_type,levels(df_p$edge_type)[c(2,1)])
df_p$variable <- factor(df_p$variable)
df_p$variable <- factor(df_p$variable,levels(df_p$variable)[c(1,3,4,2)])
##################################
### test link of demo diff with traits

# format data
df_var_mean<- df_p %>% group_by(sp, edge_type, variable) %>%
    summarise(var = var(rel_demo_change, na.rm = TRUE),
              mean = mean(rel_demo_change, na.rm = TRUE),
              PC1 = mean(PC1, na.rm = TRUE))
# exclude genus mean
df_merge <- left_join(df_var_mean, df_traits, by = c("sp" = "sp"))
df_facet <- expand.grid(variable = unique(df_var_mean$variable),
                        edge_type = unique(df_var_mean$edge_type))
# lm with weight based on inverse of variance
for (tt in traits_sel){
# exclude genus mean
df_t <- df_merge
df_t <- df_t[!is.na(df_t[[tt]]), ]
df_t <- df_t[order(df_t[[tt]]), ]
df_t[df_t[[paste0(tt, "_GenusM")]] == "YES" &
     !is.na(df_t[[paste0(tt, "_GenusM")]]), tt] <- NA
 p_vals  <-  sapply(seq_len(nrow(df_facet)),
                 function(i,df_f, df) {
                  require(nlme)
                  df_t <- df[df$variable == df_f$variable[i] &
                             df$edge_type== df_f$edge_type[i], ]
                  df_t[[tt]] <- as.vector(scale(df_t[[tt]]))
                  df_t$mean <- as.vector(scale(df_t$mean))
                  ff <- as.formula(paste0("mean ~ ", tt, ""))
                  Drop <- tree$tip.label[-c(which(tree$tip.label %in% df_t$sp))]
                  tree_t<-drop.tip(tree,Drop)
                  row.names(df_t) <- df_t$sp
                  res <- tryCatch(gls(ff, data = df_t,
                                      weights = ~var,
                                      correlation = corPagel(1, tree_t),
                                      method = "ML"),
                                  error = function(e) "NO FIT")
                  type_fit <- "PGLS Lambda"
                  if(res == "NO FIT"){
                  res <- gls(ff, data = df_t,
                             weights = ~var,
                             correlation = corBrownian(phy = tree_t),
                             method = "ML")
                  type_fit <- "PGLS Brownian"
                  }
                  interv <- intervals(res)
                  c(anova(res)$'p-value'[2],
                    interv$coef[2, 1:3], type_fit)
 },
 df_facet, df_t)
 res_t<- as.data.frame(t(p_vals))
 names(res_t) <- paste0(tt, c("_pval", "_conf.low",
                              "_std.estimate", "_conf.high", "_type_fit"))
 indx <- !grepl("_type_fit", names(res_t))
 res_t[indx] <- lapply(res_t[indx], function(x) as.numeric(as.character(x)))
 df_facet <- cbind(df_facet, res_t)
}
saveRDS(df_facet, file = file.path("output",
                                   paste0("reg_demo_edge_vs_traits_resample_phylo_",
                                          name,".rds")))
}


##############################
##############################
## EXPLORE LINKS WITH TRAITS
cor_traits_params<- function(df2, df3){
require(vegan)

df2 <- df2[df2$Latin_name!= "Pinus pinea", ]
df3 <- df3[df3$sp != "Pinus pinea", ]
df2 <-  df2[match(df3$sp, df2$Latin_name), ]
rownames(df2) <- df2$sp
df3 <- as.data.frame(df3)
rownames(df3) <- df3$Latin_name


mat_pvalue_sv<- matrix(NA, nrow = 5, ncol = 12)
rownames(mat_pvalue_sv) <- names(df2)[2:6]
colnames(mat_pvalue_sv) <- names(df3)[2:13]
mat_cor_sv<- matrix(NA, nrow = 5, ncol = 12)
rownames(mat_cor_sv) <- names(df2)[2:6]
colnames(mat_cor_sv) <- names(df3)[2:13]


for (i in 1:5){
 for(j in 1:12){
  rest <- cor.test(df2[, i+1], df3[, j+1])
  mat_pvalue_sv[i, j] <- rest$p.value
  mat_cor_sv[i, j] <- rest$estimate
 }
}


mat_pvalue_gr<- matrix(NA, nrow = 5, ncol = 12)
rownames(mat_pvalue_gr) <- names(df2)[2:6]
colnames(mat_pvalue_gr) <- names(df3)[14:25]
mat_cor_gr<- matrix(NA, nrow = 5, ncol = 12)
rownames(mat_cor_gr) <- names(df2)[2:6]
colnames(mat_cor_gr) <- names(df3)[14:25]


for (i in 1:5){
 for(j in 1:12){
  rest <- cor.test(df2[, i+1], df3[, j+13])
  mat_pvalue_gr[i, j] <- rest$p.value
  mat_cor_gr[i, j] <- rest$estimate
 }
}

return(list(mat_pvalue_gr = mat_pvalue_gr,
            mat_pvalue_sv = mat_pvalue_sv,
            mat_cor_gr= mat_cor_gr,
            mat_cor_sv = mat_cor_sv))
}



procrustes_traits_dist_summary<- function(list_1, df2){
  require(vegan)
  df1 <- list_1[[1]]
  df1b <- list_1[[2]]
  df1 <- df1[df1$sp != "Pinus pinea", ]
  df1b <- df1b[df1b$sp != "Pinus pinea", ]
  df2 <- df2[df2$Latin_name!= "Pinus pinea", ]
  df2 <-  df2[match(df1$sp, df2$Latin_name), ]
  rownames(df1) <- df1$sp
  df <- left_join(df1, df2, by = c("sp" = "Latin_name"))
  ## test correlation of ratio with traits
  vars <- c("g_m_150_hp_NoComp", "g_m_150_lp_NoComp", "s_m_150_hp_NoComp",
            "s_m_150_lp_NoComp", "lifespan_hp_NoComp", "lifespan_lp_NoComp",
            "passageTime_hp_NoComp", "passageTime_lp_NoComp",
            "sizeDeath_hp_NoComp", "sizeDeath_lp_NoComp")
  mat_cor <-  matrix(NA, nrow = length(vars), ncol = 4)
  mat_pvalue <-  matrix(NA, nrow = length(vars), ncol = 4)
  row.names(mat_cor) <- row.names(mat_pvalue) <- vars
  colnames(mat_cor) <- colnames(mat_pvalue) <- c("Wood.density.mean", "SLA.mean", "Leaf.N.mean","Seed.mass")
  traits <- c("Wood.density.mean", "SLA.mean", "Leaf.N.mean","Seed.mass")
  for (i in traits){
    for (j in vars){
   res <- cor.test(df[[j]], df[[i]], , method = "spearman")
   mat_pvalue[j,i] <- res$p.value
   mat_cor[j,i] <- res$estimate
  }
  }
  #test if ratio less than 1
  remove_outliers <- function(x, na.rm = TRUE, ...) {
    qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
    H <- 1.5 * IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    y
  }
  mat_wt <-  matrix(NA, nrow = length(vars), ncol = 3)
  row.names(mat_wt) <- vars
  colnames(mat_wt) <- c("wilcox.test.less", "wilcox.test.less.no.out", "wilcox.test.two.sided")
  for (i in vars){
   y <- remove_outliers(df[[i]])
   mat_wt[i, 1] <- wilcox.test(df[[i]], mu = 0, alternative = "less")$p.value
   mat_wt[i, 3] <- wilcox.test(df[[i]], mu = 0)$p.value
   mat_wt[i, 2] <- wilcox.test(y, mu = 0, alternative = "less")$p.value
  }

  ## pdf("test.pdf", width = 20, height = 10)
  ## par(mfrow = c(2,5), mar = c(1,1,1,1))
  ## for (i in vars){
  ## boxplot(df[[i]], outline = FALSE, notch = TRUE)
  ## abline(h = 0, col = "red")
  ## }
  ## dev.off()

  # ANALYSIS ON PCA OF DIST AND TRAITS OR TRAIT PER TRAIT
  print(paste0("Number of species with NA ", sum(!complete.cases(df1))))
  df1 <- df1[complete.cases(df1), ]
  df2 <-  df2[df2$Latin_name %in% df1$sp, ]
  ## TODO SOLVE PROBLEM WITH SPECIES WITH NA IN PASSAGE TIME
  protest_IPM_dist_traits <- vegan::protest(prcomp(dplyr::select(df1, -sp, -PC1_median), scale. = TRUE),
                                            prcomp(dplyr::select(df2, -Latin_name), scale. = TRUE))

  pca_IPM <- prcomp(dplyr::select(df1, -sp, -PC1_median), scale. = TRUE)
  envfit_pca_IPM_traits<- envfit(pca_IPM, dplyr::select(df2, -Latin_name), perm = 999)

  return(list(mat_pvalue = mat_pvalue,
              mat_cor = mat_cor,
              mat_wt = mat_wt,
              protest_IPM_dist_traits = protest_IPM_dist_traits ,
              envfit_pca_IPM_traits = envfit_pca_IPM_traits))
}

procrustes_traits_params<- function(df3, df2){
require(vegan)

df2 <- df2[df2$Latin_name!= "Pinus pinea", ]
df2 <-  df2[match(df1$sp, df2$Latin_name), ]
df3 <- df3[df3$sp != "Pinus pinea", ]
df3 <-  df3[match(df1$sp, df3$sp), ]
rownames(df3) <- df3$sp
df3 <- as.data.frame(df3)

protest_params_traits <- vegan::protest(prcomp(dplyr::select(df3, -sp), scale. = TRUE),
                                        prcomp(dplyr::select(df2, -Latin_name), scale. = TRUE))

pca_params <- prcomp(dplyr::select(df3, -sp), scale. = TRUE)
envfit_pca_params_traits <- envfit(pca_params, dplyr::select(df2, -Latin_name), perm = 999)

df3b <- df3 %>%  select_if(grepl("wai", names(.)) |
                             grepl("sgdd", names(.)) |
                           names(.) == "sp")
pca_paramsC <- prcomp(dplyr::select(df3b, -sp), scale. = TRUE)
envfit_pca_paramsClim_traits <- envfit(pca_paramsC, dplyr::select(df2, -Latin_name), perm = 999)

df3b <- df3 %>%  select_if(grepl("sv_", names(.)) | names(.) == "sp")
pca_params <- prcomp(dplyr::select(df3b, -sp), scale. = TRUE)
envfit_pca_param_sv_traits <- envfit(pca_params, dplyr::select(df2, -Latin_name), perm = 999)

df3b <- df3 %>%  select_if(grepl("gr_", names(.)) | names(.) == "sp")
pca_params <- prcomp(dplyr::select(df3b, -sp), scale. = TRUE)
envfit_pca_param_gr_traits<- envfit(pca_params, dplyr::select(df2, -Latin_name), perm = 999)

return(list(protest_params_traits = protest_params_traits,
            envfit_pca_params_traits = envfit_pca_params_traits,
            envfit_pca_params_traits = envfit_pca_params_traits,
            envfit_pca_params_sv_traits = envfit_pca_params_sv_traits,
            envfit_pca_params_gr_traits = envfit_pca_params_gr_traits,
            envfit_pca_paramsClim_traits = envfit_pca_paramsClim_traits))
}




format_IPM_metrics_all_sp_TwoClimVar<- function(file, sps, data, data_pred){
ll <- readRDS(file.path("output", file))
seq_i <- seq_len(length.out = length(sps$sp))
list_df <- lapply(seq_i, format_IPM_metrics_TwoClimVar, ll = ll,
                  spvec = sps$sp, data = data, data_pred = data_pred)
 df <- bind_rows(list_df)
 return(df)
}

format_IPM_metrics_all_sp_pca<- function(file, sps, data, data_pred,
                                         clim_pca, data_PC1_mean){
ll <- readRDS(file.path("output", file))
seq_i <- seq_len(length.out = length(sps$sp))
list_df <- lapply(seq_i, format_IPM_metrics_pca, ll = ll,
                  spvec = sps$sp, data = data,
                  data_pred = data_pred, clim_pca = clim_pca,
                  data_PC1_mean = data_PC1_mean)
 df2 <- bind_rows(list_df)
 return(df2)
}

format_IPM_metrics_all_sp_pca_Resample<- function(file, sps, data, data_pred,
                                         clim_pca, data_PC1_mean){
ll <- readRDS(file.path("output", file))
seq_i <- seq_len(length.out = length(sps$sp))
list_df <- lapply(seq_i, format_IPM_metrics_pca_Resample, ll = ll,
                  spvec = sps$sp, data = data,
                  data_pred = data_pred, clim_pca = clim_pca,
                  data_PC1_mean = data_PC1_mean)
 df2 <- bind_rows(list_df)
 return(df2)
}



