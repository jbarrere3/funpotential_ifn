###############################
###############################
### R functions to fit growth and survival models

# SMALL GLM OBJECT

stripGlmLR <- function(cm) {
  cm$y = c()
  cm$model = c()

  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr$qr = c()
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()
  cm$family$variance = c()
  cm$family$dev.resids = c()
  cm$family$validmu = c()
  cm$family$simulate = c()
  cm$offset <-  c()
  #attr(cm$terms,".Environment") = c()
  attr(cm$formula,".Environment") = c()

  cm
}


## COMPUTE HARVESTING RATE

format_data_harvest <- function(data, spsel){
  df <- data %>% filter(sp == spsel) %>% filter(!is.na(yearsbetweensurveys) &
                                                treestatus_th!= 1 &
                                                !country %in% c("ES","FR"))
  if(nrow(df)>10){
    # keep alive dead and harvested tree
    df$harv <- 1
    df$harv[df$treestatus_th %in% c(2, 4,5)] <- 0
    if(sum(df$harv ==1)>1){
      # remove country with zero mortality as this seems to lead
      # to a problem in the estimation
      DA_country <- table(df$harv, df$country)
      country_to_remove <- colnames(DA_country)[DA_country[2, ] ==0]
      if(length(country_to_remove)>0){
        df <- df[df$country != country_to_remove, ]
      }
      df <- df %>% filter(!is.na(wai) & !is.na(sgdd)
                          & !is.na(BATOTcomp))
      return(df)
    }else{
      return(NA)
    }
  }else{
     return(NA)
  }
}



compute_sp_harvesting_rate <- function(spsel,  data){
    require(dplyr)
    print(spsel)
      df <- format_data_harvest(data, spsel)
    if(!is.na(df)){
      res <- glm(harv ~ offset(log(yearsbetweensurveys)),
                 family = binomial(link = "cloglog"),
                 data = df, control = glm.control(maxit = 50))
      preds <- predict(res, newdata = data.frame(harv = 0,
                                                 yearsbetweensurveys = 1),
                       se.fit = TRUE, type = "link")
      critval <- 1.96 ## approx 95% CI
      upr <- preds$fit + (critval * preds$se.fit)
      lwr <- preds$fit - (critval * preds$se.fit)
      fit <- preds$fit
      fit2 <- res$family$linkinv(fit)
      upr2 <- res$family$linkinv(upr)
      lwr2 <- res$family$linkinv(lwr)
    return(data.frame(sp = rep(spsel, length(fit)),
                      mean = fit2, lwr = lwr2, upr = upr2,
                      nobs = nrow(df)))
    }else{
    return(data.frame(sp = spsel,
                      mean = NA, lwr = NA, upr = NA,
                      nobs = NA))
    }
}


compute_sp_harvesting_rate_all_sp<- function(sps,  data){
  df <- data %>% filter(treestatus_th %in% c(2, 3, 4, 5) & !is.na(wai) & !is.na(sgdd)) %>% filter(!country %in% c("WA", "FI")) %>%
      filter(!(country == "FR" & surveydate2< 2009))
  df <- df %>% group_by(plotcode) %>% mutate(ncut = sum(treestatus_th %in% c(3)),
                                             ntot = n(),
                                             pcut = ncut/ntot) %>%
                                     filter(pcut <= 0.8) # remove clear cut
  ll <- lapply(sps$sp, compute_sp_harvesting_rate, data = df)
  df_harv <- dplyr::bind_rows(ll)
# all sp
  df <- df %>% filter(!is.na(yearsbetweensurveys))
  # keep alive dead and harvested tree
  df$harv <- 1
  df$harv[df$treestatus_th %in% c(2, 4,5)] <- 0
  df <- df %>% filter(!is.na(wai) & !is.na(sgdd)
                      & !is.na(BATOTcomp))
      res <- glm(harv ~ offset(log(yearsbetweensurveys)),
                 family = binomial(link = "cloglog"),
                 data = df, control = glm.control(maxit = 50))
      preds <- predict(res, newdata = data.frame(harv = 0,
                                                 yearsbetweensurveys = 1),
                       se.fit = TRUE, type = "link")
      critval <- 1.96 ## approx 95% CI
      upr <- preds$fit + (critval * preds$se.fit)
      lwr <- preds$fit - (critval * preds$se.fit)
      fit <- preds$fit
      fit2 <- res$family$linkinv(fit)
      upr2 <- res$family$linkinv(upr)
      lwr2 <- res$family$linkinv(lwr)
 res_allsp <- data.frame(sp = rep("mean", length(fit)),
                      mean = fit2, lwr = lwr2, upr = upr2,
                      nobs = nrow(df))
return(rbind(df_harv, res_allsp))
}


################################
################################
### FORMAT DATA

# function to remove outlier based on http://ctfs.si.edu/Public/CTFSRPackage/index.php/web/topics/growth~slash~growth.r/trim.growth
fun.remove.outlier <-  function(G, D, year, err.limit = 1.5 ){
TF <- (G <= 25) & (D+year*G) > (D -err.limit*(0.006214*D + 0.9036)) &
      !is.na(G) & !is.na(D) & !is.na(year)
return(TF)
}

format_data_growth <- function(data, spsel){
 df <- data %>% dplyr::filter(treestatus_th == 2) %>% # keep trees alive only
              mutate(size = dbh1,
                     incr = (dbh2-dbh1)/yearsbetweensurveys,
                     sizeNext = dbh1 +incr) %>%
              dplyr::filter(!is.na(incr) & !is.na(size) & !is.na(BATOTcomp) & !is.na(sgdd) & !is.na(wai))
  df <- df[df$country != "WA", ]
  df <- df[fun.remove.outlier(df$incr, df$dbh1, df$yearsbetweensurveys), ]
  df <- df %>% dplyr::filter(sp == spsel) %>%
              dplyr::filter(incr > 0  &
                     incr <  quantile(incr, probs = 0.98))
  # remove outlier and negative growth
  df$logincr <- log(df$incr)
  df$logsize <- log(df$size)
  df <- df %>% arrange(sgdd) %>% mutate(sgddb= 1/(sgdd),
                                        waib = 1/(wai + 1),
                                        sgdd2 = sgdd^2,
                                        wai2 = wai^2,
                                        spei_meanb = 1/(spei_mean + 1.3),
                                        logsgdd= log(sgdd),
                                        logwai = log(wai + 1),
                                        logspei_mean = log(spei_mean + 1.3))
return(df)
}

format_data_survival <- function(data, spsel){
df_m <- read.csv("data/FunDivEUROPE_plot_management.csv")
  data <-  data[data$plotcode %in%  df_m$plotcode[df_m$management2==0 |
                                                  is.na(df_m$management2)], ]

  df <- data[data$treestatus_th %in% c(2, 4) & data$sp == spsel, ]
  df <- df[!is.na(df$yearsbetweensurveys), ]
  df <- df[df$country != "WA", ]
  df <- df[!(df$country == "FR" & df$surveydate2< 2009), ]
  df <- df[!is.na(df$wai) & !is.na(df$sgdd), ]
  # keep alive and naturaly dead tree
  df$size <- df$dbh1
  df$dead <- df$treestatus_th
  df$dead[df$treestatus_th == 4] <- 1
  df$dead[df$treestatus_th == 2] <- 0
  df$logsize <- log(df$size)
  df$surv <- 1 - df$dead
  # remove country with zero mortality as this seems to lead to a problem in the estimation
  DA_country <- table(df$dead, df$country)
  country_to_remove <- colnames(DA_country)[DA_country[2, ] ==0]
  if(length(country_to_remove)>0){
    df <- df[df$country != country_to_remove, ]
  }
  df <- df %>% arrange(sgdd) %>% filter(!is.na(wai) & !is.na(sgdd)
                                        & !is.na(BATOTcomp))  %>%
             mutate(sgddb= 1/sgdd,
                    waib = 1/(wai + 1),
                    sgdd2 = sgdd^2,
                    wai2 = wai^2,
                    spei_meanb = 1/(spei_mean + 1.3),
                    logsgdd= log(sgdd),
                    logwai = log(wai + 1),
                    logspei_mean = log(spei_mean + 1.3))
  return(df)
}



# BUILD AN AIC TABLE TO CHOSE CLIMATIC VARIABLE
fit_one_var <- function(var, data){
require(lme4)
require(MuMIn)
if(length(table(data$country))>1){
form <- as.formula(paste("logincr~ country + 0 + size + logsize + BATOTcomp +", var, "+(1|plotcode)"))
}else{
form <- as.formula(paste("logincr~ size + logsize + BATOTcomp +", var, "+(1|plotcode)"))
}
res <-  lmer(formula = form, data)
return(c(AIC(res), MuMIn::r.squaredGLMM(res)[1]))
}

fit_one_var_error<- function(var, data){
tryCatch(fit_one_var(var, data),
      error=function(e) e)
}

fit_clim_growth <- function(spsel, data){
df <-format_data_growth(data, spsel)
clim_vars <- c("sgdd", "wai",
               "sgddb", "waib",
               "logsgdd", "logwai",
               "poly(sgdd, 2)", "poly(wai,2)",
               "poly(sgdd, 3)", "poly(wai,3)")
df <- df[complete.cases(df[, names(df) %in% clim_vars]), ]
#exclude alll obs with at least one missing var
l_select <- lapply(clim_vars, fit_one_var_error, data = df)
return(l_select)
}

fit_all_clim_growth <- function(sps, data){
ll_select <- lapply(sps$sp, fit_clim_growth, data)
ll_select
}


print_aic<- function(sps, ll){
clim_vars <- c("sgdd", "wai",
               "sgddb", "waib",
               "logsgdd", "logwai",
               "poly(sgdd, 2)", "poly(wai,2)",
               "poly(sgdd, 3)", "poly(wai,3)")
resAIC <- t(sapply(ll, function(x) sapply(x, function(z) z[1])))
resR2 <-  t(sapply(ll, function(x) sapply(x, function(z) z[2])))
colnames(resAIC) <- colnames(resR2) <- clim_vars
rownames(resAIC) <- rownames(resR2) <- sps$sp
resAIC_sgdd <- resAIC[, c(1,3,5,7)]
resAIC_wai <- resAIC[, c(1,3,5,7)+1]
DeltaAIC_sgdd <- resAIC_sgdd - apply(resAIC_sgdd, 1, min)
DeltaAIC_wai <- resAIC_wai - apply(resAIC_wai, 1, min)
WeightAIC_sgdd<- exp(-0.5*DeltaAIC_sgdd)/apply(exp(-0.5*DeltaAIC_sgdd), 1, sum)
WeightAIC_wai<- exp(-0.5*DeltaAIC_wai)/apply(exp(-0.5*DeltaAIC_wai), 1, sum)
resR2_sgdd <- resR2[, c(1,3,5,7)]
resR2_wai <- resR2[, c(1,3,5,7)+1]
return(list(sgdd = list(resAIC_sgdd, DeltaAIC_sgdd, WeightAIC_sgdd, resR2_sgdd),
            wai = list(resAIC_wai, DeltaAIC_wai, WeightAIC_wai, resR2_wai)))
}

aic_per_var <- function(gg){
select_wai <- grepl("wai", colnames(gg[[1]]))
Daic_wai <- rbind(gg[[1]][ , select_wai]- apply(gg[[1]][ , select_wai], 1, min),
                 apply((gg[[1]][ , select_wai]- apply(gg[[1]][ , select_wai],
                                                      1, min)),
                       2, sum)/nrow(gg[[1]]))
R2_wai <- rbind(-gg[[2]][ , select_wai]+ apply(gg[[2]][ , select_wai], 1, max),
                 apply((-gg[[2]][ , select_wai]+apply(gg[[2]][ , select_wai],
                                                      1, max)),
                       2, mean))
rownames(R2_wai) <- rownames(Daic_wai) <- c(rownames(gg[[1]]), "tot")
select_sgdd <- grepl("sgdd", colnames(gg[[1]]))
Daic_sgdd <- rbind(gg[[1]][ , select_sgdd]- apply(gg[[1]][ , select_sgdd],
                                                1, min),
                   apply((gg[[1]][, select_sgdd]- apply(gg[[1]][ , select_sgdd],
                                                        1, min)),
                         2, sum)/nrow(gg[[1]]))
R2_sgdd <- rbind(-gg[[2]][, select_sgdd]+apply(gg[[2]][, select_sgdd], 1, max),
                 apply((-gg[[2]][ , select_sgdd]+apply(gg[[2]][ , select_sgdd],
                                                      1, max)),
                       2, mean))
rownames(R2_sgdd) <- rownames(Daic_sgdd) <- c(rownames(gg[[1]]), "tot")
par(mfrow = c(2,2), mar = c(2,2,2,2))
barplot(R2_wai["tot", ])
barplot(R2_sgdd["tot", ])
barplot(Daic_wai["tot", ])
barplot(Daic_sgdd["tot", ])
}

# select climatic variables for survival
fit_one_var_s<- function(var, data){
if(length(table(data$country))>1){
form <- as.formula(paste("dead~ country + 0 + size + logsize + BATOTcomp +",
                         var,"+ offset(log(yearsbetweensurveys))"))
}else{
form <- as.formula(paste("dead~ size + logsize + BATOTcomp +",
                         var, "+offset(log(yearsbetweensurveys))"))
}
res <- glm(form,
           family = binomial(link = "cloglog"),
           data = data, control = glm.control(maxit = 50))
return(AIC(res))
}

fit_one_var_error_s<- function(var, data){
tryCatch(fit_one_var_s(var, data),
      error=function(e) e)
}

fit_clim_survival <- function(spsel, data){
df <-format_data_survival(data, spsel)
clim_vars <- c("sgdd", "wai",
               "sgddb", "waib",
               "logsgdd", "logwai",
               "poly(sgdd, 2)", "poly(wai,2)",
               "poly(sgdd, 3)", "poly(wai,3)")
df <- df[complete.cases(df[, names(df) %in% clim_vars]), ]
#exclude alll obs with at least one missing var
l_select <- lapply(clim_vars, fit_one_var_error_s, data = df)
return(l_select)
}

fit_all_clim_survival <- function(sps, data){
ll_select <- lapply(sps$sp, fit_clim_survival, data)
ll_select
}

print_aic_s<- function(sps, ll){
clim_vars <- c("sgdd", "wai",
               "sgddb", "waib",
               "logsgdd", "logwai",
               "poly(sgdd, 2)", "poly(wai,2)",
               "poly(sgdd, 3)", "poly(wai,3)")
resAIC <- t(sapply(ll, function(x) sapply(x, function(z) z[1])))
colnames(resAIC) <- clim_vars
return(resAIC)
}

aic_per_var_s<- function(ss, spsel){
select_wai <- grepl("wai", colnames(ss))
Daic_wai <- rbind(ss[ , select_wai]- apply(ss[ , select_wai], 1, min),
                 apply((ss[ , select_wai]- apply(ss[ , select_wai],
                                                      1, min))==0,
                       2, sum)/nrow(ss))
select_sgdd <- grepl("sgdd", colnames(ss))
Daic_sgdd <- rbind(ss[ , select_sgdd]- apply(ss[ , select_sgdd],
                                                1, min),
                   apply((ss[, select_sgdd]- apply(ss[ , select_sgdd],
                                                        1, min))==0,
                         2, sum)/nrow(ss))
par(mfrow = c(1,2), mar = c(2,2,2,2))
barplot(Daic_wai[28, ])
barplot(Daic_sgdd[28, ])
}

######################################
######################################
#### TEST AND SELECT INTERACTION

growth_size_climate_inter_test<- function(spsel = "Pinus sylvestris", data) {
      require(MuMIn)
   options(na.action = "na.fail")
print(spsel)
    require(lmerTest)
#spsel <- "Quercus suber" #"Pinus sylvestris" "Fagus sylvatica"
df <- format_data_growth(data, spsel)
var_to_scale <- c("size", "logsize", "BATOTcomp", "sgddb", "waib")
df[, var_to_scale] <- scale(df[, var_to_scale])
if(length(unique(df$country))>1){
    res <-  lmer(logincr~ country + size+ logsize+sgddb  +
                  waib +BATOTcomp +
                  BATOTcomp:sgddb + size:sgddb+logsize:sgddb +
                  BATOTcomp:waib + size:waib+logsize:waib+
                 (1|plotcode)
             , df)
   anova_res<- anova(res)
   anova_pvalt <-  anova_res[["Pr(>F)"]]
   names(anova_pvalt) <- rownames( anova_res)
   anova_pval <- anova_pvalt[-c(1)]
   dd <- dredge(res, fixed = c("country", "size", "logsize", "BATOTcomp",
                               "sgddb", "waib"))
  }else{
    res <-  lmer(logincr~ size+ logsize+sgddb  +
                  waib +BATOTcomp +
                  BATOTcomp:sgddb + size:sgddb+logsize:sgddb +
                  BATOTcomp:waib + size:waib+logsize:waib+
                 (1|plotcode)
               , df)
   anova_res<- anova(res)
   anova_pvalt <-  anova_res[["Pr(>F)"]]
   names(anova_pvalt) <- rownames( anova_res)
   anova_pval <- anova_pvalt
   dd <- dredge(res, fixed = c("size", "logsize", "BATOTcomp", "sgddb", "waib"))
  }
    print(paste(spsel, "done"))
  return(list(anova_pval ,
              as.data.frame(dd[, ])[1, !names(dd) %in% c("df", "logLik",
                                                         "AICc", "delta",
                                                         "weight")],
              as.data.frame(dd[, ])[1, names(dd) %in% c("weight")]))
}

growth_size_climate_inter_test_clim<- function(spsel = "Pinus sylvestris", data) {
      require(MuMIn)
   options(na.action = "na.fail")
print(spsel)
    require(lmerTest)
#spsel <- "Quercus suber" #"Pinus sylvestris" "Fagus sylvatica"
df <- format_data_growth(data, spsel)
var_to_scale <- c("size", "logsize", "BATOTcomp", "sgddb", "waib")
df[, var_to_scale] <- scale(df[, var_to_scale])
if(length(unique(df$country))>1){
    res <-  lmer(logincr~ country + size+ logsize+sgddb  +
                  waib +BATOTcomp + waib:sgddb +
                 (1|plotcode)
             , df)
   anova_res<- anova(res)
   anova_pvalt <-  anova_res[["Pr(>F)"]]
   names(anova_pvalt) <- rownames( anova_res)
   anova_pval <- anova_pvalt[-c(1)]
   dd <- dredge(res, fixed = c("country", "size", "logsize", "BATOTcomp",
                               "sgddb", "waib"))
  }else{
    res <-  lmer(logincr~ size+ logsize+sgddb  +
                  waib +BATOTcomp + waib:sgddb +
                 (1|plotcode)
               , df)
   anova_res<- anova(res)
   anova_pvalt <-  anova_res[["Pr(>F)"]]
   names(anova_pvalt) <- rownames( anova_res)
   anova_pval <- anova_pvalt
   dd <- dredge(res, fixed = c("size", "logsize", "BATOTcomp", "sgddb", "waib"))
  }
    print(paste(spsel, "done"))
  return(list(anova_pval ,
              as.data.frame(dd[, ])[1, !names(dd) %in% c("df", "logLik",
                                                         "AICc", "delta",
                                                         "weight")],
              as.data.frame(dd[, ])[1, names(dd) %in% c("weight")]))
}


growth_size_climate_inter_test_poly3<- function(spsel = "Pinus sylvestris",
                                                data) {
      require(MuMIn)
   options(na.action = "na.fail")
print(spsel)
    require(lmerTest)
df <- format_data_growth(data, spsel)
var_to_scale <- c("size", "logsize", "BATOTcomp", "sgdd", "wai")
df[, var_to_scale] <- scale(df[, var_to_scale])
if(length(unique(df$country))>1){
    res <-  lmer(logincr~ country + size+ logsize+poly(sgdd,3)  +
                  poly(wai,3) +BATOTcomp +
                  BATOTcomp:poly(sgdd,3) + size:poly(sgdd,3)+logsize:poly(sgdd,3) +
                  BATOTcomp:poly(wai,3) + size:poly(wai,3)+logsize:poly(wai,3)+
                 (1|plotcode)
             , df)
   anova_res<- anova(res)
   anova_pvalt <-  anova_res[["Pr(>F)"]]
   names(anova_pvalt) <- rownames( anova_res)
   anova_pval <- anova_pvalt[-c(1)]
   dd <- dredge(res, fixed = c("country", "size", "logsize", "BATOTcomp",
                               "poly(sgdd,3)", "poly(wai,3)"))
  }else{
    res <-  lmer(logincr~ size+ logsize+poly(sgdd,3)  +
                  poly(wai,3) +BATOTcomp +
                  BATOTcomp:poly(sgdd,3) + size:poly(sgdd,3)+
                 logsize:poly(sgdd,3) +
                  BATOTcomp:poly(wai,3) + size:poly(wai,3)+logsize:poly(wai,3)+
                 (1|plotcode)
               , df)
   anova_res<- anova(res)
   anova_pvalt <-  anova_res[["Pr(>F)"]]
   names(anova_pvalt) <- rownames( anova_res)
   anova_pval <- anova_pvalt
   dd <- dredge(res, fixed = c("size", "logsize", "BATOTcomp", "poly(sgdd,3)", "poly(wai,3)"))
  }
    print(paste(spsel, "done"))
  return(list(anova_pval ,
              as.data.frame(dd[, ])[1, !names(dd) %in% c("df", "logLik",
                                                         "AICc", "delta",
                                                         "weight")],
              as.data.frame(dd[, ])[1, names(dd) %in% c("weight")]))
}

table_growth_size_climate_inter_test_all<- function(sps, data) {
    tab <- lapply(sps$sp, growth_size_climate_inter_test, data)
 return(tab)
}

table_growth_size_climate_inter_test_all_poly3<- function(sps, data) {
    tab <- lapply(sps$sp, growth_size_climate_inter_test_poly3, data)
 return(tab)
}

table_growth_size_climate_inter_test_all_clim<- function(sps, data) {
    tab <- lapply(sps$sp, growth_size_climate_inter_test_clim, data)
 return(tab)
}

table_growth_size_climate_inter_test_all2<- function(sps, tab) {
 mat_pval <- t(sapply(tab, function(x) x[[1]]))
 TOT_pval <-  apply(mat_pval<0.05, MARGIN = 2, sum)/nrow(mat_pval)
 mat_pval2 <- rbind(mat_pval<0.05, TOT_pval)
 row.names(mat_pval2) <-  c(sps$sp, "Tot percentage")
 f <- function(x){
  v <- unlist(x[[2]])
  return(v[!names(v) %in% c("country", "(Intercept)")])
 }
 mat_select <- t(sapply(tab,f ))
 TOT_select <-  apply(!is.na(mat_select), MARGIN = 2, sum)/nrow(mat_pval)
 mat_select2 <- rbind(!is.na(mat_select), TOT_select)
 row.names(mat_select2) <-  c(sps$sp, "Tot percentage")
 row.names(mat_pval) <-  sps$sp
require(knitr)
print( kable(mat_select2[, !colnames(mat_select2) %in% c("size", "logsize",
                                                         "BATOTcomp", "sgddb",
                                                         "waib")],
       caption = "Selected variables in best models", digits = 2))
print( kable(mat_pval2[, !colnames(mat_pval2) %in% c("size", "logsize",
                                                     "BATOTcomp", "sgddb",
                                                     "waib")],
       caption = "Significative variables in full models", digits = 2))

return(list(pval = mat_pval2, var_select = mat_select2))
}


###############
## Survival interaction test

survival_size_climate_inter_test<- function(spsel = "Pinus sylvestris", data) {
   require(MuMIn)
   options(na.action = "na.fail")
df <- format_data_survival(data, spsel)
var_to_scale <- c("size", "logsize", "BATOTcomp", "sgdd", "waib")
df[, var_to_scale] <- scale(df[, var_to_scale])
if(length(unique(df$country))>1){
    res <- glm(dead~ country + size + logsize +
               sgddb+ waib +BATOTcomp+
               sgddb:size+sgddb:logsize+ sgddb:BATOTcomp+
               waib:size+waib:logsize+ waib:BATOTcomp+
               offset(log(yearsbetweensurveys)),
               family = binomial(link = "cloglog"),
               data = df, control = glm.control(maxit = 50))
   anova_res<- anova(res, test = "Chisq")
   anova_pvalt <-  anova_res[["Pr(>Chi)"]]
   names(anova_pvalt) <- rownames( anova_res)
   anova_pval <- anova_pvalt[-(1:2)]
   dd <- dredge(res, fixed = c("country", "size", "logsize", "BATOTcomp", "sgddb",
                               "waib", "offset(log(yearsbetweensurveys))"))
  }else{
    res <- glm(dead~ size + logsize +
               sgddb+ waib +BATOTcomp+
               sgddb:size+sgddb:logsize+ sgddb:BATOTcomp+
               waib:size+waib:logsize+ waib:BATOTcomp+
               offset(log(yearsbetweensurveys)),
               family = binomial(link = "cloglog"),
               data = df, control = glm.control(maxit = 50))
   anova_res<- anova(res, test = "Chisq")
   anova_pvalt <-  anova_res[["Pr(>Chi)"]]
   names(anova_pvalt) <- rownames( anova_res)
   anova_pval <- anova_pvalt[-(1)]
   dd <- dredge(res, fixed = c("size", "logsize", "BATOTcomp", "sgddb",
                               "waib", "offset(log(yearsbetweensurveys))"))
  }
    print(paste(spsel, "done"))
  return(list(anova_pval,
              as.data.frame(dd[, ])[1, !names(dd) %in% c("df", "logLik", "AICc", "delta", "weight")],
              as.data.frame(dd[, ])[1, names(dd) %in% c("weight")]))
}

survival_size_climate_inter_test_clim<- function(spsel = "Pinus sylvestris", data) {
   require(MuMIn)
   options(na.action = "na.fail")
df <- format_data_survival(data, spsel)
var_to_scale <- c("size", "logsize", "BATOTcomp", "sgdd", "wai")
df[, var_to_scale] <- scale(df[, var_to_scale])
if(length(unique(df$country))>1){
    res <- glm(dead~ country + size + logsize +
               sgddb+ waib +BATOTcomp+ sgddb:waib+
               offset(log(yearsbetweensurveys)),
               family = binomial(link = "cloglog"),
               data = df, control = glm.control(maxit = 50))
   anova_res<- anova(res, test = "Chisq")
   anova_pvalt <-  anova_res[["Pr(>Chi)"]]
   names(anova_pvalt) <- rownames( anova_res)
   anova_pval <- anova_pvalt[-(1:2)]
   dd <- dredge(res, fixed = c("country", "size", "logsize", "BATOTcomp", "sgddb",
                               "waib", "offset(log(yearsbetweensurveys))"))
  }else{
    res <- glm(dead~ size + logsize +
               sgddb+ waib +BATOTcomp+ sgddb:waib+
               offset(log(yearsbetweensurveys)),
               family = binomial(link = "cloglog"),
               data = df, control = glm.control(maxit = 50))
   anova_res<- anova(res, test = "Chisq")
   anova_pvalt <-  anova_res[["Pr(>Chi)"]]
   names(anova_pvalt) <- rownames( anova_res)
   anova_pval <- anova_pvalt[-(1)]
   dd <- dredge(res, fixed = c("size", "logsize", "BATOTcomp", "sgddb",
                               "waib", "offset(log(yearsbetweensurveys))"))
  }
    print(paste(spsel, "done"))
  return(list(anova_pval,
              as.data.frame(dd[, ])[1, !names(dd) %in% c("df", "logLik", "AICc", "delta", "weight")],
              as.data.frame(dd[, ])[1, names(dd) %in% c("weight")]))
}

survival_size_climate_inter_test_poly3<- function(spsel = "Pinus sylvestris", data) {
   require(MuMIn)
   options(na.action = "na.fail")
df <- format_data_survival(data, spsel)
var_to_scale <- c("size", "logsize", "BATOTcomp", "sgdd", "wai")
df[, var_to_scale] <- scale(df[, var_to_scale])
if(length(unique(df$country))>1){
    res <- glm(dead~ country + size + logsize +
               poly(sgdd,3)+ poly(wai,3) +BATOTcomp+
               poly(sgdd,3):size+poly(sgdd,3):logsize+ poly(sgdd,3):BATOTcomp+
               poly(wai,3):size+poly(wai,3):logsize+ poly(wai,3):BATOTcomp+
               offset(log(yearsbetweensurveys)),
               family = binomial(link = "cloglog"),
               data = df, control = glm.control(maxit = 50))
   anova_res<- anova(res, test = "Chisq")
   anova_pvalt <-  anova_res[["Pr(>Chi)"]]
   names(anova_pvalt) <- rownames( anova_res)
   anova_pval <- anova_pvalt[-(1:2)]
   dd <- dredge(res, fixed = c("country", "size", "logsize", "BATOTcomp", "poly(sgdd,3)",
                               "poly(wai,3)", "offset(log(yearsbetweensurveys))"))
  }else{
    res <- glm(dead~ size + logsize +
               poly(sgdd,3)+ poly(wai,3) +BATOTcomp+
               poly(sgdd,3):size+poly(sgdd,3):logsize+ poly(sgdd,3):BATOTcomp+
               poly(wai,3):size+poly(wai,3):logsize+ poly(wai,3):BATOTcomp+
               offset(log(yearsbetweensurveys)),
               family = binomial(link = "cloglog"),
               data = df, control = glm.control(maxit = 50))
   anova_res<- anova(res, test = "Chisq")
   anova_pvalt <-  anova_res[["Pr(>Chi)"]]
   names(anova_pvalt) <- rownames( anova_res)
   anova_pval <- anova_pvalt[-(1)]
   dd <- dredge(res, fixed = c("size", "logsize", "BATOTcomp", "poly(sgdd,3)",
                               "poly(wai,3)", "offset(log(yearsbetweensurveys))"))
  }
    print(paste(spsel, "done"))
  return(list(anova_pval,
              as.data.frame(dd[, ])[1, !names(dd) %in% c("df", "logLik", "AICc", "delta", "weight")],
              as.data.frame(dd[, ])[1, names(dd) %in% c("weight")]))
}

table_survival_size_climate_inter_test_all<- function(sps, data) {
 tab <- lapply(sps$sp, survival_size_climate_inter_test, data)
 return(tab)
}

table_survival_size_climate_inter_test_all_poly3<- function(sps, data) {
 tab <- lapply(sps$sp, survival_size_climate_inter_test_poly3, data)
 return(tab)
}

table_survival_size_climate_inter_test_all_clim<- function(sps, data) {
 tab <- lapply(sps$sp, survival_size_climate_inter_test_clim, data)
 return(tab)
}

table_survival_size_climate_inter_test_all2<- function(sps, tab) {
 mat_pval <- t(sapply(tab, function(x) x[[1]]))
 TOT_pval <-  apply(mat_pval<0.05, MARGIN = 2, sum)/nrow(mat_pval)
 mat_pval2 <- rbind(mat_pval<0.05, TOT_pval)
 row.names(mat_pval2) <-  c(sps$sp, "Tot percentage")
 f <- function(x){
  v <- unlist(x[[2]])
  return(v[!names(v) %in% c("country", "(Intercept)")])
 }
 mat_select <- t(sapply(tab,f ))
 TOT_select <-  apply(!is.na(mat_select), MARGIN = 2, sum)/nrow(mat_pval)
 mat_select2 <- rbind(!is.na(mat_select), TOT_select)
 row.names(mat_select2) <-  c(sps$sp, "Tot percentage")
 row.names(mat_pval) <-  sps$sp
require(knitr)
mat_sel <- mat_select2[, !colnames(mat_select2) %in% c("size", "logsize", "BATOTcomp",
                                                       "sgddb", "waib", "offset(log(yearsbetweensurveys))")]
print( kable(mat_sel,
       caption = "Selected variables in best models", digits = 2))
mat_p <- mat_pval2[, !colnames(mat_pval2) %in% c("size", "logsize", "BATOTcomp", "sgddb",
                                                 "waib", "offset(log(yearsbetweensurveys))")]
print( kable(mat_p,
       caption = "All Significative variables in full models", digits = 2))

return(list(pvalue = mat_p, var_select = mat_sel))
}




##################################
##################################
###### FIT GROWTH MODELS
##

growth_size_climate <- function(df, vars) {
    require(lme4)
f1 <- paste0("logincr~ country + 0 + size+ logsize + BATOTcomp +",
             paste0(vars, collapse = " + "),
             " + (1|plotcode)")
f2 <- paste0("logincr~ size+ logsize + BATOTcomp +",
             paste(vars, collapse = " + "),
             " + (1|plotcode)")
if(length(unique(df$country))>1){
    res <-  lmer(f1, df, REML = FALSE)
  }else{
    res <-  lmer(f2, df, REML = FALSE)
  }
  return(res)
}

predict_for_var <- function(var, vars, df, res){
  require(lme4)
  seq_var <- seq(from = quantile(df[[var]], probs = 0.00001),
                  to = quantile(df[[var]], probs = 0.99999),
                  length.out = 100)
  #pred for partial residual
  names_0 <- vars[vars != var]
  vec_0 <- rep(0, length.out = length(names_0))
  names(vec_0) <- names_0
  list_temp <- as.list(vec_0)
  list_temp[[var]] <- seq_var
  list_temp[["country"]] <- unique(df$country)
  df_pred_p<- do.call(expand.grid, list_temp)

  if(var == "size"){
      df_pred_p$logsize <- log(df_pred_p$size)
  }else{
      df_pred_p$logsize <- 0
  }
  df_pred_i <- df_pred_p
  df_pred_i[, vars] <- 0
  df_pred_i$logsize <- 0
  df_pred_p$id <-  as.integer(factor(df_pred_p[[var]]))
  df_pred_i$id <-  as.integer(factor(df_pred_i[[var]]))

  # residual
  df_pred_pr <-  df
  df_pred_pr[, var] <-  0
  if(var == "size"){
      df_pred_pr$logsize <- 0
  }
  #
  #pred for mean
  vec_m <- apply(df[ ,names_0 ], MARGIN = 2, mean, na.rm = TRUE)
  list_temp <- as.list(vec_m)
  list_temp[[var]] <- as.vector(seq_var)
  list_temp[["country"]] <- unique(df$country)
  df_pred_pm<- do.call(expand.grid, list_temp)
  df_pred_pm$logsize <- log(df_pred_pm$size)

  df_pred_pm$id <-  as.integer(factor(df_pred_pm[[var]]))
  pred_c <- predict(res, newdata = df_pred_p,
                    re.form=NA) - predict(res, newdata = df_pred_i,
                                          re.form=NA)
  lin_pred_res<- tapply(pred_c, INDEX = df_pred_p$id, mean)
  pred_cm <- predict(res, newdata = df_pred_pm,
                             re.form=NA)
  lin_pred_m<- tapply(pred_cm, INDEX = df_pred_pm$id, mean)
  pr <- df$logincr - (predict(res, newdata = df_pred_pr, re.form=NA))
return(list(pred=lin_pred_res, pr = pr,
            pred_m = lin_pred_m,
            seq_var = seq_var))
}


predict_for_var_mean<- function(var, vars, df, res, df2){
  require(lme4)
  seq_var <- seq(from = quantile(df2[[var]], probs = 0.00001),
                  to = quantile(df2[[var]], probs = 0.99999),
                  length.out = 100)
  #pred for mean
  names_0 <- vars[vars != var]
  vec_m <- apply(df[ ,names_0 ], MARGIN = 2, mean, na.rm = TRUE)
  list_temp <- as.list(vec_m)
  list_temp[[var]] <- as.vector(seq_var)
  list_temp[["country"]] <- unique(df$country)
  df_pred_pm<- do.call(expand.grid, list_temp)
  df_pred_pm$logsize <- log(df_pred_pm$size)

  df_pred_pm<- df_pred_pm %>%
                mutate(sgddb= 1/(sgdd),
                       waib = 1/(wai + 1),
                       sgdd2 = sgdd^2,
                       wai2 = wai^2)

  df_pred_pm$id <-  as.integer(factor(df_pred_pm[[var]]))
  pred_cm <- predict(res, newdata = df_pred_pm,
                             re.form=NA)
  lin_pred_m<- tapply(pred_cm, INDEX = df_pred_pm$id, mean)
return(list(pred_m = lin_pred_m,
            seq_var = seq_var))
}


predict_for_inter_xy<- function(x, y, vars, df, res){
  require(lme4)
  seq_var <- seq(from = quantile(df[[x]], probs = 0.001),
                  to = quantile(df[[x]], probs = 0.999),
                  length.out = 100)
  level_y <- quantile(df[[y]], probs = c(0.2, 0.4, 0.6, 0.8))
  #pred for partial residual
  names_m <- vars[!vars %in% c(x,y)]
  vec_m <- lapply(names_m, function(x, df) mean(df[[x]]) , df = df)
  names(vec_m) <- names_m
  list_temp <- as.list(vec_m)
  list_temp[[x]] <- seq_var
  list_temp[[y]] <- level_y
  list_temp[["country"]] <- unique(df$country)
  df_pred_p<- do.call(expand.grid, list_temp)
  df_pred_p$logsize <- log(df_pred_p$size)
  df_pred_p$id <-  as.integer(factor(paste(df_pred_p[[x]],df_pred_p[[y]]), ordered = FALSE))
  ## add category to obs
  breakss <- c(min(df[[y]])*0.99,
               (level_y[-1] + level_y[-length(level_y)])/2,
              max(df[[y]])*1.01)
  df$levels_y<- level_y[cut(df[[y]], breakss, labels = FALSE)]

  # PRED
  pred_c <- predict(res, newdata = df_pred_p,
                    re.form=NA)
  lin_pred_res<- tapply(pred_c, INDEX = df_pred_p$id, mean)

  seq_var <- tapply(df_pred_p[[x]], INDEX = df_pred_p$id, mean)
  levels_y_p <- tapply(df_pred_p[[y]], INDEX = df_pred_p$id, mean)
  df_pred <- data.frame(id = names(lin_pred_res), var = seq_var, levels_y = levels_y_p, pred = lin_pred_res)
  df_pred <- df_pred[order(df_pred$var), ]
  pr <- df$logincr
return(list(pred=df_pred$pred,
            pr = pr,
            seq_var = df_pred$var,
            levels_y_p = df_pred$levels_y,
            levels_y_pr= df$levels_y))
}


plot_one_vars <- function(var, vars, df, res, spsel, title = FALSE,
                          type = "partial_resdiual"){
list_pred <- predict_for_var(var, vars, df, res)
xlabb <- var
x <- list_pred$seq_var
# for transformed variables
if (var == "waib"){
xlabb <- "wai"
x <- 1/x - 1
var <- "wai"
}
if (var == "sgddb"){
xlabb <- "sgdd"
x <- 1/x
var <- "sgdd"
}
if (var == "spei_meanb"){
xlabb <- "spei_mean"
x <- 1/x -1.3
var <- "spei_mean"
}
if (var == "logwai"){
xlabb <- "wai"
x <- exp(x) - 1
var <- "wai"
}
if (var == "logsgdd"){
xlabb <- "sgdd"
x <- exp(x)
var <- "sgdd"
}
if (var == "logspei_mean"){
xlabb <- "spei_mean"
x <- exp(x) -1.3
var <- "spei_mean"
}
sig <- rep(sigma(res),length(x))
if(type == "partial_residual"){
p <- plot.bin2d(df[[var]], list_pred$pr,
                 xlab = xlabb, ylab = "Partial residual of log(incr)")+
         geom_line(aes(x, y),
                   data.frame(x = x,
                              y = list_pred$pred),
                   size = 1.5, colour = "red") +
         geom_line(aes(x, y),
                   data.frame(x = x,
                              y = list_pred$pred +2*sig),
                   size = 1.5, colour = "red",
                   linetype = "dashed") +
         geom_line(aes(x, y),
                   data.frame(x = x,
                              y = list_pred$pred -2*sig),
                   size = 1.5, colour = "red",
                   linetype = "dashed")

}
if(type == "mean"){
p <- plot.bin2d(df[[var]], df$logincr,
                 xlab = xlabb, ylab = "log(incr)")+
         geom_line(aes(x, y),
                   data.frame(x = x,
                              y = list_pred$pred_m),
                   size = 1.5, colour = "red")
}

if(title == TRUE) p <- p +  ggtitle(spsel)

return(p)
}


plot_inter_vars <- function(var, var_y, vars, df, res, spsel){
list_pred <- predict_for_inter_xy(var, var_y, vars, df, res)
xlabb <- var
x <- list_pred$seq_var
# for transformed variables
if (var == "waib"){
xlabb <- "wai"
x <- 1/x - 1
var <- "wai"
}
if (var == "sgddb"){
xlabb <- "sgdd"
x <- 1/x
var <- "sgdd"
}
if (var == "spei_meanb"){
xlabb <- "spei_mean"
x <- 1/x -1.3
var <- "spei_mean"
}
if (var == "logwai"){
xlabb <- "wai"
x <- exp(x) - 1
var <- "wai"
}
if (var == "logsgdd"){
xlabb <- "sgdd"
x <- exp(x)
var <- "sgdd"
}
if (var == "logspei_mean"){
xlabb <- "spei_mean"
x <- exp(x) -1.3
var <- "spei_mean"
}

df$pr <-  list_pred$pr
df$levels_y<- list_pred$levels_y_pr
plot_list <-  vector("list")
for (j in seq_len(length.out = length(unique(df$levels_y)))){
i <- unique(df$levels_y)[j]
plot_list[[j]]<- plot.bin2d(df[[var]][df$levels_y == i], list_pred$pr[df$levels_y == i],
                            xlab = xlabb, ylab = "log(incr)")+
              geom_line(aes(x, y),
                   data.frame(x = x[list_pred$levels_y_p== i],
                              y = list_pred$pred[list_pred$levels_y_p== i]),
                   size = 1.5, colour = "red") +
              ggtitle(paste0(var_y, " = ", round(i,5)))
}

return(plot_list)
}




growth_size_climate_plot_pred<- function(spsel = "Pinus sylvestris",
                                         vars = c("sgddb", "waib"),
                                         data, type = "partial_residual"){
  df <-format_data_growth(data, spsel)
  res <- growth_size_climate(df, vars)
  vars <- gsub(",[0-9]\\)", "", gsub("poly\\(", "", vars))
  fixed_vars <- c("size", "logsize", "BATOTcomp", vars)
  p1 <- plot_one_vars(var = "size", vars = fixed_vars,
                      df = df, res = res,
                      spsel =spsel, type = type)
  p2 <- plot_one_vars(var = "BATOTcomp", vars = fixed_vars,
                      df = df, res = res,
                      spsel =spsel, type = type)
  p3 <- plot_one_vars(var = vars[1], vars = fixed_vars,
                      df = df, res = res,
                      spsel =spsel, type = type)
  p4 <- plot_one_vars(var = vars[2], vars = fixed_vars,
                      df = df, res = res,
                      spsel =spsel, type = type, title = TRUE)

  print(multiplot(p1, p2, p3, p4,
                  cols = 2))
}

growth_size_climate_plot_pred_inter<- function(spsel = "Pinus sylvestris",
                                                   vars = c("waib","sgddb"),
                                                   inters = c("size:waib", "logsize:waib",
                                                              "size:sgddb", "logsize:sgddb"),
                                                   data){
  df <-format_data_growth(data, spsel)
  res <- growth_size_climate(df, c(vars, inters))
  vars <- gsub(",[0-9]\\)", "", gsub("poly\\(", "", vars))
  inters1 <- gsub(":.*","", inters)
  inters2 <- gsub(".*:","", inters)
  inters1b <- inters1[inters1 != "logsize"]
  inters2b <- inters2[inters1 != "logsize"]
  fixed_vars <- c("size", "logsize", "BATOTcomp", vars)
  p1 <- plot_one_vars(var = "BATOTcomp", vars = fixed_vars,
                      df = df, res = res,
                      spsel =spsel, type = "mean")
  df0 <- data.frame()
  p2 <- ggplot(df0) + geom_point() + theme_bw() +
      geom_text(data = data.frame(text = spsel),
                aes(label = text), size = 5, x = 0.5, y = 0.5)
  p3 <- ggplot(df0) + geom_point() + theme_bw()
  list_inter_plots <- vector("list")

  for (i in seq_len(length.out = length(inters1b))){
  list_inter_plots[[i]] <- plot_inter_vars(inters1b[i], inters2b[i], vars = fixed_vars, df, res, spsel)
  }
  list_inter <- do.call(c,list_inter_plots)
  list_plots <- do.call(c,list(list(p1, p2, p3, p3), list_inter))
  print(multiplot(plotlist = list_plots,
                  layout = matrix(seq_len(length.out = length(list_plots)), ncol = 4, byrow = TRUE)))
}


growth_size_climate_plot_fit<- function(spsel = "Pinus sylvestris",
                                        data,
                                        vars = c("sgddb", "waib")){
  df <-format_data_growth(data, spsel)
  res <- growth_size_climate(df, vars)

  #test dist of residual with DHARMa from Florian
  require(DHARMa)
  simulationOutput <- simulateResiduals(fittedModel = res, n = 250)
  n <- nrow(df)
  m <- (1:n)/(n + 1)
  df_qq <-  data.frame(m = m , stdres = sort(simulationOutput$scaledResiduals))
  pqq <- ggplot(df_qq, aes(x = m, y = stdres)) + geom_point() +
          geom_abline(intercept = 0, slope = 1, colour = "red")+
          xlab("expected") + ylab("Observed")+
          theme_bw() +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.position = "none")+
          ggtitle(spsel)
  p0 <- plot.bin2d( log(df$incr), predict(res, type = "response"),
             xlab = "observed", ylab = "predicted")+
         geom_abline(intercept = 0, slope = 1, colour = "red")
  p1 <- plot.bin2d( df$size, resid(res),
           xlab = "size", ylab = "residual")

  png(file.path("figures", paste0(gsub(" ", "_", spsel),"fit_ggplot.png")),
       width = 960, height = 480)
  multiplot(pqq, p0, p1,cols = 2)
  dev.off()
}




growth_size_climate_all_pred_ggplot_waib<- function(sps, data) {
  pdf(file.path("figures", "all_growth_pred_ggplot_waib.pdf"))
    l <- lapply(sps$sp, growth_size_climate_plot_pred, data = data)
  dev.off()
}


growth_size_climate_all_pred_ggplot_waib_inter<- function(sps, data) {
  pdf(file.path("figures", "all_growth_pred_ggplot_waib_inter.pdf"))
    l <- lapply(sps$sp, growth_size_climate_plot_pred_inter, data = data)
  dev.off()
}

growth_size_climate_all_pred_ggplot_waib_inter_clim<- function(sps, data) {
  pdf(file.path("figures", "all_growth_pred_ggplot_waib_inter_clim.pdf"))
    l <- lapply(sps$sp, growth_size_climate_plot_pred_inter, data = data,
                inters = c("size:waib", "logsize:waib","size:sgddb", "logsize:sgddb", "waib:sgddb"))
  dev.off()
}


growth_size_climate_all_pred_ggplot_m_waib<- function(sps, data) {
  pdf(file.path("figures", "all_growth_pred_ggplot_m_waib.pdf"))
    l <- lapply(sps$sp, growth_size_climate_plot_pred, data = data, type = "mean")
  dev.off()
}

growth_size_climate_all_pred_ggplot_poly3wai<- function(sps, data) {
  pdf(file.path("figures", "all_growth_pred_ggplot_poly3wai.pdf"))
    l <- lapply(sps$sp, growth_size_climate_plot_pred,
                vars = c("poly(sgdd,3)", "poly(wai,3)"),
                data = data)
  dev.off()
}

growth_size_climate_all_pred_ggplot_m_poly3wai<- function(sps, data) {
  pdf(file.path("figures", "all_growth_pred_ggplot_m_poly3wai.pdf"))
    l <- lapply(sps$sp, growth_size_climate_plot_pred,
                vars = c("poly(sgdd,3)", "poly(wai,3)"),
                data = data, type = "mean")
  dev.off()
}

growth_size_climate_all_pred_ggplot_poly2wai<- function(sps, data) {
  pdf(file.path("figures", "all_growth_pred_ggplot_poly2wai.pdf"))
    l <- lapply(sps$sp, growth_size_climate_plot_pred,
                vars = c("poly(sgdd,2)", "poly(wai,2)"),
                data = data)
  dev.off()
}

growth_size_climate_all_pred_ggplot_m_poly2wai<- function(sps, data) {
  pdf(file.path("figures", "all_growth_pred_ggplot_m_poly2wai.pdf"))
    l <- lapply(sps$sp, growth_size_climate_plot_pred,
                vars = c("poly(sgdd,2)", "poly(wai,2)"),
                data = data, type = "mean")
  dev.off()
}

growth_size_climate_all_pred_ggplot_wai<- function(sps, data) {
  pdf(file.path("figures", "all_growth_pred_ggplot_wai.pdf"))
    l <- lapply(sps$sp, growth_size_climate_plot_pred, vars = c("sgdd", "wai"), data = data)
  dev.off()
}

growth_size_climate_all_pred_ggplot_m_wai<- function(sps, data) {
  pdf(file.path("figures", "all_growth_pred_ggplot_m_wai.pdf"))
    l <- lapply(sps$sp, growth_size_climate_plot_pred, vars = c("sgdd", "wai"),
                data = data, type = "mean")
  dev.off()
}

growth_size_climate_all_pred_ggplot_logwai<- function(sps, data) {
  pdf(file.path("figures", "all_growth_pred_ggplot_logwai.pdf"))
    l <- lapply(sps$sp, growth_size_climate_plot_pred,
                vars = c("logsgdd", "logwai"), data = data)
  dev.off()
}

growth_size_climate_all_pred_ggplot_m_logwai<- function(sps, data) {
  pdf(file.path("figures", "all_growth_pred_ggplot_m_logwai.pdf"))
    l <- lapply(sps$sp, growth_size_climate_plot_pred, vars = c("logsgdd", "logwai"),
                data = data, type = "mean")
  dev.off()
}


growth_size_climate_all_fit_ggplot<- function(sps, data) {
    l <- lapply(sps$sp, growth_size_climate_plot_fit, data)
    system("convert  figures/*fit_ggplot.png figures/all_growth_fit_ggplot.pdf")
    system("rm  figures/*fit_ggplot.png")
}

###########################
###########################
### FIT SURVIVAL MODELS

predict_for_var_surv<- function(var, vars, df, res){
  seq_var <- seq(from = quantile(df[[var]], probs = 0.00001),
                  to = quantile(df[[var]], probs = 0.99999),
                  length.out = 100)
  names_m <- vars[vars != var]
  list_m <- lapply(df[, names_m], median)
  list_m[[var]] <- seq_var
  list_m[["country"]] <- unique(df$country)
  df_pred_p<- do.call(expand.grid, list_m)
  df_pred_p$logsize <- log(df_pred_p$size)
  df_pred_p$id <-  as.integer(factor(df_pred_p[[var]]))
  df_pred_p$yearsbetweensurveys <- 1
  pred_c <- predict(res, newdata = df_pred_p, type = "response")
  pred_m<- tapply(pred_c, INDEX = df_pred_p$id, mean)

return(list(pred = pred_m, var_seq = seq_var))
}

predict_for_var_surv_mean<- function(var, vars, df, res, df2){
  seq_var <- seq(from = quantile(df2[[var]], probs = 0.00001),
                  to = quantile(df2[[var]], probs = 0.99999),
                  length.out = 100)
  names_m <- vars[vars != var]
  list_m <- lapply(df[, names_m], mean)
  list_m[[var]] <- seq_var
  list_m[["country"]] <- unique(df$country)
  df_pred_p<- do.call(expand.grid, list_m)
  df_pred_p$logsize <- log(df_pred_p$size)
  df_pred_p$id <-  as.integer(factor(df_pred_p[[var]]))
  df_pred_p$yearsbetweensurveys <- 1
  df_pred_p<- df_pred_p %>%
                mutate(sgddb= 1/(sgdd),
                       waib = 1/(wai + 1),
                       sgdd2 = sgdd^2,
                       wai2 = wai^2)
  pred_c <- predict(res, newdata = df_pred_p, type = "response")
  pred_m<- tapply(pred_c, INDEX = df_pred_p$id, mean)

return(list(pred = pred_m, var_seq = seq_var))
}



predict_for_inter_xy_surv<- function(x,y, vars, df, res){
  seq_var <- seq(from = quantile(df[[x]], probs = 0.001),
                  to = quantile(df[[x]], probs = 0.999),
                  length.out = 100)
  level_y <- quantile(df[[y]], probs = c(0.2, 0.4, 0.6, 0.8))
  #pred for partial residual
  names_m <- vars[!vars %in% c(x,y)]
  list_temp <- lapply(df[, names_m], mean)
  list_temp[[x]] <- seq_var
  list_temp[[y]] <- level_y
  list_temp[["country"]] <- unique(df$country)
  df_pred_p<- do.call(expand.grid, list_temp)
  df_pred_p$logsize <- log(df_pred_p$size)
  df_pred_p$id <-  as.integer(factor(paste(df_pred_p[[x]],df_pred_p[[y]]), ordered = FALSE))
  df_pred_p$yearsbetweensurveys <- 1
  pred_c <- predict(res, newdata = df_pred_p, type = "response")
  pred_m<- tapply(pred_c, INDEX = df_pred_p$id, mean)
  seq_var <- tapply(df_pred_p[[x]], INDEX = df_pred_p$id, mean)
  levels_y_p <- tapply(df_pred_p[[y]], INDEX = df_pred_p$id, mean)
  df_pred <- data.frame(id = names(pred_m), var = seq_var, levels_y = levels_y_p, pred = pred_m)
  df_pred <- df_pred[order(df_pred$var), ]

return(list(pred = df_pred$pred, var_seq = df_pred$var, levels_y = df_pred$levels_y))
}


survival_size_climate <- function(df, vars = c("waib", "sgddb"), strip = FALSE){
f1 <- paste0("dead~ country + 0 + size+ logsize + BATOTcomp +",
             paste0(vars, collapse = " + "),
             " + offset(log(yearsbetweensurveys))")
f2 <- paste0("dead~ size+ logsize + BATOTcomp +",
             paste0(vars, collapse = " + "),
             " + offset(log(yearsbetweensurveys))")

  if(length(unique(df$country))>1){
#    print("include country")
    res <- glm(f1,
               family = binomial(link = "cloglog"),
               data = df, control = glm.control(maxit = 50))
  }else{
    res <- glm(f2,
               family = binomial(link = "cloglog"),
               data = df, control = glm.control(maxit = 50))
  }
if (strip) res <- stripGlmLR(res)
return(res)
}

survival_size_stan<- function(df, country_TF){
# based on code from Needham et al. 2018 PRSBS
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
n.iter <- 700
max.size <- max(df$size)
lowr2 <- ifelse(max.size < 200, -0.5, -0.1)
upr2 <- ifelse(max.size < 200, -0.01, -0.00001)
thresh <- ifelse(max.size < 50, median(size),max.size*0.2)
if(country_TF){
surv.data <- list(N = length(df$size),
                  C = length(unique(df$country)),
                  Ccode = as.vector(unclass(factor(df$country))),
                  dbh = df$size,
                  surv = 1 - df$dead,
                  time = df$yearsbetweensurveys,
                  thresh = thresh,
                  lowr2 = lowr2,
                  upr2 = upr2)

surv.model <- stan(file = 'R/survCountry.stan',
                   data = surv.data, chains = 0)
surv.fit <- stan(fit = surv.model, data = surv.data, iter = n.iter,
                 chains = 3, control = list(adapt_delta = 0.98), cores = 3)
# extract the parameters
surv.mcmc <- lapply(extract(surv.fit), median)
surv.chains <- do.call(cbind, extract(surv.fit))
}else{
surv.data <- list(N = length(df$size),
                  dbh = df$size,
                  surv = 1 - df$dead,
                  time = df$yearsbetweensurveys,
                  thresh = thresh,
                  lowr2 = lowr2,
                  upr2 = upr2)

surv.model <- stan(file = 'R/surv.stan',
                   data = surv.data, chains = 0)
surv.fit <- stan(fit = surv.model, data = surv.data, iter = n.iter,
                 chains = 3, control = list(adapt_delta = 0.98), cores = 3)
# extract the parameters
surv.mcmc <- lapply(extract(surv.fit), median)
surv.chains <- do.call(cbind, extract(surv.fit))
}
return(list(fit = surv.mcmc, thresh = thresh))
}

predict_surv_stan <- function(res, df){
pred1 <- 1- mean(res$fit$K)/(1+exp(-res$fit$r1*(seq(100, res$thresh, length.out = 100) - res$fit$p1)))
pred2 <- 1- mean(res$fit$K)/(1+exp(-res$fit$r2*(seq(res$thresh, max(df$size), length.out = 100) - res$fit$p2)))
var_seq <- c(seq(100, res$thresh, length.out = 100),
             seq(res$thresh, max(df$size), length.out = 100))
return(list(pred = c(pred1, pred2), var_seq = var_seq))
}

survival_size <- function(df){

f1 <- "dead~ country + 0 + size+ logsize + offset(log(yearsbetweensurveys))"
f2 <- "dead~ size+ logsize + offset(log(yearsbetweensurveys))"
f1p2 <- "dead~ country + 0 + poly(size,2) + offset(log(yearsbetweensurveys))"
f2p2 <- "dead~ poly(size,2) + offset(log(yearsbetweensurveys))"
f1p2l <- "dead~ country + 0 + poly(size,2)+ logsize + offset(log(yearsbetweensurveys))"
f2p2l <- "dead~ poly(size,2)+ logsize + offset(log(yearsbetweensurveys))"
  require(gam)
  if(length(unique(df$country))>1){
    print("include country")
    res <- glm(f1,
               family = binomial(link = "cloglog"),
               data = df, control = glm.control(maxit = 50))

    resT <- glm(dead ~country + 0 + logsize + offset(log(yearsbetweensurveys)),
               family = binomial(link = "cloglog"),
               data = df, control = glm.control(maxit = 50))
    o<-segmented(resT,seg.Z=~logsize, psi = log(max(df$size)*0.2))
    resgam <- gam(dead~ country + 0 + s(size, 4) + offset(log(yearsbetweensurveys)),
                family = binomial(link = "cloglog"),
               data = df)
    resp2 <- glm(f1p2,
               family = binomial(link = "cloglog"),
               data = df, control = glm.control(maxit = 50))
    resp2l <- glm(f1p2l,
               family = binomial(link = "cloglog"),
               data = df, control = glm.control(maxit = 50))
    resW <- glm(dead~ country + 0 + size+ logsize + offset(log(yearsbetweensurveys)),
                family = binomial(link = "cloglog"),
               data = df,
               weights = df$size/100)
    resp2W <- glm(dead~ country + 0 + poly(size,2) + offset(log(yearsbetweensurveys)),
                family = binomial(link = "cloglog"),
               data = df,
               weights = df$size/100)
    resp2lW <- glm(dead~ country + 0 + poly(size,2)+ logsize + offset(log(yearsbetweensurveys)),
                family = binomial(link = "cloglog"),
               data = df,
               weights = df$size/100)
    resStan <- survival_size_stan(df, country_TF = TRUE)
  }else{
    res <- glm(f2,
               family = binomial(link = "cloglog"),
               data = df, control = glm.control(maxit = 50))
    resT <- glm(dead~logsize + offset(log(yearsbetweensurveys)),
               family = binomial(link = "cloglog"),
               data = df, control = glm.control(maxit = 50))
    o<-segmented(resT,seg.Z=~logsize, psi = log(max(df$size)*0.2))
    resgam <- gam(dead~  s(size, 4) + offset(log(yearsbetweensurveys)),
                family = binomial(link = "cloglog"),
               data = df)
    resp2 <- glm(f2p2,
               family = binomial(link = "cloglog"),
               data = df, control = glm.control(maxit = 50))
    resp2l <- glm(f2p2l,
               family = binomial(link = "cloglog"),
               data = df, control = glm.control(maxit = 50))
    resW <- glm(dead~  size+ logsize + offset(log(yearsbetweensurveys)),
                family = binomial(link = "cloglog"),
               data = df,
               weights = df$size/100)
    resp2W <- glm(dead~ poly(size,2) + offset(log(yearsbetweensurveys)),
                family = binomial(link = "cloglog"),
               data = df,
               weights = df$size/100)
    resp2lW <- glm(dead~ poly(size,2)+ logsize + offset(log(yearsbetweensurveys)),
                family = binomial(link = "cloglog"),
               data = df,
               weights = df$size/100)
    resStan <- survival_size_stan(df, country_TF = FALSE)
  }

print(coef(res))
print(coef(resW))
return(list(res, o, resgam, resp2, resp2l, resW, resp2W, resp2lW, resStan))
}


plot_surv_bin_eval <- function(df, res, ncuts = 30, legend_TF = FALSE, spsel = NA){
 cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02")
 names(cols) <- c("DE", "ES", "FI", "FR", "SW", "WA")
 df_bin <- Surv_Var_Bin(df, var = "proba", ncuts = ncuts, type_breaks = "quant",
                              KeepCountry = TRUE)
 par(mfrow = c(1,2))
 plot(df_bin$var, df_bin$mean, pch = 16,
      cex = df_bin$nobs/max(df_bin$nobs)+0.5,
      col = cols[as.character(df_bin$country)],
      xlab = "Predicted Proba of mortality", ylab = "Proba of mortality",
      ylim = range(0, df_bin$surv, df_bin$lwr, df_bin$upr, na.rm = TRUE))
  segments(df_bin$var, df_bin$lwr ,df_bin$var, df_bin$upr,
           col = cols[as.character(df_bin$country)])
 abline(a = 0, b = 1, lwd = 2)
 if(legend_TF){
    legend(quantile(df_bin$var, probs = 0.8),
           quantile(df_bin$upr, probs = 0.99),
           legend = names(cols),
           col =cols,
           pch = 16, bty = "n")
    title(main = spsel)
 }
  df_mean <- Compute_Surv_CI_mean(df, KeepC = TRUE)
  vec <- tapply(df$proba, INDEX = df$country, mean, na.rm = TRUE)
  seq_c <- seq_len(length.out = nrow(df_mean))
 plot(seq_c, df_mean$mean, col = cols[as.character(df_mean$country)],
       xaxt = "n", ylab = "mean proba of mortality", pch = 16,
       ylim = range(0, df_bin$surv, df_bin$lwr, df_bin$upr, na.rm = TRUE))
  segments(seq_c, df_mean$lwr ,seq_c, df_mean$upr,
           col = cols[as.character(df_mean$country)], lwd =2)
  axis(side = 1, at = seq_c, labels = FALSE, tck = -0.01, xlab = "country")
  points(seq_c, vec[df_mean$country], pch = 4, cex = 2, col = cols[as.character(df_mean$country)])

}



plot_surv_bin_pred <- function(var, vars, df, res, ncuts = 30, legend_TF = FALSE, spsel = NA){
cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02")
names(cols) <- c("DE", "ES", "FI", "FR", "SW", "WA")
pp <- predict_for_var_surv(var = var, vars,  df, res)
varb <-  var
xp <- pp$var_seq
if(var ==  "sgddb"){
 varb <- "sgdd"
 xp <- 1/xp
 }
if(var ==  "waib"){
 varb <- "wai"
 xp <- 1/xp - 1
}
if (var == "spei_meanb"){
 varb <- "spei_mean"
 xp <- 1/xp -1.3
}
if (var == "logwai"){
 varb <- "wai"
 xp <- exp(xp) - 1
}
if (var == "logsgdd"){
 varb <- "sgdd"
 xp <- exp(xp)
}
if (var == "logspei_mean"){
 varb <- "spei_mean"
 xp <- exp(xp) -1.3
}
df_bin <- Surv_Var_Bin(df, var = varb, ncuts = ncuts, type_breaks = "equal",
                             KeepCountry = TRUE)
plot(df_bin$var, df_bin$mean, pch = 16,
     cex = log(df_bin$nobs/max(df_bin$nobs, na.rm = TRUE)+0.1)+2.8,
     col = cols[as.character(df_bin$country)],
     xlab = varb, ylab = "Proba of mortality",
     ylim = range(0, df_bin$surv,  df_bin$lwr,pp$pred, na.rm = TRUE)) # df_bin$upr,
segments(df_bin$var, df_bin$lwr ,df_bin$var, df_bin$upr,
         col = cols[as.character(df_bin$country)])
lines(xp, pp$pred, col = "red", lwd = 3)
if(legend_TF){
    legend(quantile(df_bin$var, probs = 0.8, na.rm = TRUE),
           quantile(df_bin$upr, probs = 0.99, na.rm = TRUE),
           legend = names(cols),
           col =cols,
           pch = 16, bty = "n")
    title(main = spsel)
 }
}




plot_surv_bin_pred_inter<- function(var, vary, vars, df, res, ncuts = 30, legend_TF = FALSE, spsel = NA){
cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02")
names(cols) <- c("DE", "ES", "FI", "FR", "SW", "WA")
pp <- predict_for_inter_xy_surv(x = var, y = vary, vars,  df, res)
varb <-  var
xp <- pp$var_seq
yp <- pp$levels_y

if(var ==  "sgddb"){
 varb <- "sgdd"
 xp <- 1/xp
 }
if(var ==  "waib"){
 varb <- "wai"
 xp <- 1/xp - 1
}
if (var == "spei_meanb"){
 varb <- "spei_mean"
 xp <- 1/xp -1.3
}
if (var == "logwai"){
 varb <- "wai"
 xp <- exp(xp) - 1
}
if (var == "logsgdd"){
 varb <- "sgdd"
 xp <- exp(xp)
}
if (var == "logspei_mean"){
 varb <- "spei_mean"
 xp <- exp(xp) -1.3
}

levels_y <- unique(yp)
df_bin <- Surv_Inter_Bin(df, var_x= varb, var_y = vary, var_y_levels = levels_y,
                         ncuts = ncuts, type_breaks = "quant",
                         KeepCountry = TRUE)

for (i in seq_len(length.out = length(levels_y)-1)){
plot(df_bin$var_x[df_bin$var_y == levels_y[i]], df_bin$mean[df_bin$var_y == levels_y[i]], pch = 16,
     cex = df_bin$nobs[df_bin$var_y == levels_y[i]]/max(df_bin$nobs)+1,
     col = cols[as.character(df_bin$country[df_bin$var_y == levels_y[i]])],
     xlab = varb, ylab = "Proba of mortality",
     ylim = range(0, df_bin$surv, df_bin$lwr, df_bin$upr, pp$pred, na.rm = TRUE))
segments(df_bin$var_x[df_bin$var_y == levels_y[i]], df_bin$lwr[df_bin$var_y == levels_y[i]] ,
         df_bin$var_x[df_bin$var_y == levels_y[i]], df_bin$upr[df_bin$var_y == levels_y[i]],
         col = cols[as.character(df_bin$country[df_bin$var_y == levels_y[i]])])
lines(xp[yp == levels_y[i]], pp$pred[yp == levels_y[i]], col = "red", lwd = 3)
abline(v = mean(df[[varb]], na.rm = TRUE),lwd= 2)
title(main = paste(vary, " = ", round(levels_y[i],5)))
}
i <- length(unique(yp))
plot(df_bin$var_x[df_bin$var_y == levels_y[i]], df_bin$mean[df_bin$var_y == levels_y[i]], pch = 16,
     cex = df_bin$nobs[df_bin$var_y == levels_y[i]]/max(df_bin$nobs)+1,
     col = cols[as.character(df_bin$country[df_bin$var_y == levels_y[i]])],
     xlab = varb, ylab = "Proba of mortality",
     ylim = range(0, df_bin$surv, df_bin$lwr, df_bin$upr, pp$pred, na.rm = TRUE))
segments(df_bin$var_x[df_bin$var_y == levels_y[i]], df_bin$lwr[df_bin$var_y == levels_y[i]] ,
         df_bin$var_x[df_bin$var_y == levels_y[i]], df_bin$upr[df_bin$var_y == levels_y[i]],
         col = cols[as.character(df_bin$country[df_bin$var_y == levels_y[i]])])
lines(xp[yp == levels_y[i]], pp$pred[yp == levels_y[i]], col = "red", lwd = 3)
if(legend_TF){
    legend(quantile(df_bin$var_x, probs = 0.8),
           max(df_bin$upr),
           legend = names(cols),
           col =cols,
           pch = 16, bty = "n")
    title(main = paste(vary, " = ", round(levels_y[i],5), " ", spsel))
}
}





plot_survival_eval<- function(spsel = "Pinus sylvestris",
                                data,
                                vars = c("waib", "sgddb"),
                                plot_type = "pred_surv") {
  df <- format_data_survival(data, spsel)
  df2 <- df
  res <- survival_size_climate(df, vars)
  df2$yearsbetweensurveys <- 1
  df$proba <- res$family$linkinv(predict(res, newdata = df2, type = "link"))
  plot_surv_bin_eval(df, res, ncuts = 20, legend_TF = TRUE, spsel = spsel)
}


plot_survival_size<- function(spsel = "Pinus sylvestris",
                                         data) {
    df <- format_data_survival(data, spsel)
    print(spsel)
    res <- survival_size(df)
    plot_surv_bin_pred(var= "size", vars = c("BATOTcomp"),  df, res[[1]],
                       ncuts = 50, legend_TF = TRUE, spsel = spsel)
    pp_o<- predict_for_var_surv(var = "size", vars = c("BATOTcomp"),  df, res[[2]])
    pp_gam<- predict_for_var_surv(var = "size", vars = c("BATOTcomp"),  df, res[[3]])
    ## pp_p2<- predict_for_var_surv(var = "size", vars = c("BATOTcomp"),  df, res[[4]])
    pp_p2l<- predict_for_var_surv(var = "size", vars = c("BATOTcomp"),  df, res[[5]])
    pp_W<- predict_for_var_surv(var = "size", vars = c("BATOTcomp"),  df, res[[6]])
    ## pp_p2W<- predict_for_var_surv(var = "size", vars = c("BATOTcomp"),  df, res[[7]])
    pp_p2lW<- predict_for_var_surv(var = "size", vars = c("BATOTcomp"),  df, res[[8]])
    pp_stan <- predict_surv_stan(res[[9]], df)
    ## lines(pp_p2$var_seq, pp_p2$pred, col = "blue", lwd = 3)
    lines(pp_o$var_seq, pp_o$pred, col = "blue", lty = 1, lwd = 6)
    lines(pp_gam$var_seq, pp_gam$pred, col = "orange", lty = 1, lwd = 3)
    lines(pp_p2l$var_seq, pp_p2l$pred, col = "green", lwd = 3)
    lines(pp_W$var_seq, pp_W$pred, col = "red", lty = 2, lwd = 3)
    ## lines(pp_p2W$var_seq, pp_p2W$pred, col = "blue", lty = 2, lwd = 3)
    lines(pp_p2lW$var_seq, pp_p2lW$pred, col = "green", lty = 2, lwd = 3)
    lines(pp_stan$var_seq, pp_stan$pred, col = "purple", lty = 2, lwd = 3)

}



plot_survival_size_climate<- function(spsel = "Pinus sylvestris",
                                         data,
                                         vars = c("waib", "sgddb"),
                                         plot_type = "pred_surv") {
  df <- format_data_survival(data, spsel)
  res <- survival_size_climate(df, vars)
  #remove interactions
  vars <- vars[!grepl(":", vars)]
  if(grepl("poly", vars[1])) vars <- c("wai", "sgdd")
  vars2 <- c("size", "BATOTcomp", vars)
  if(plot_type == "pred_surv"){
    par(mfrow = c(2, 2))
    plot_surv_bin_pred(var= "size", vars = vars2,  df, res)
    plot_surv_bin_pred(var= "BATOTcomp", vars = vars2,  df, res)
    plot_surv_bin_pred(var= vars[1], vars = vars2,  df, res)
    plot_surv_bin_pred(var= vars[2], vars = vars2,  df, res, legend_TF = TRUE, spsel = spsel)
  }
  if(plot_type == "fit_surv"){
      require(DHARMa)
    simulationOutput <- simulateResiduals(res)
    n <- nrow(df)
    m <- (1:n)/(n + 1)
    df_qq <-  data.frame(m = m , stdres = sort(simulationOutput$scaledResiduals))
    p1 <- plot.bin2d( simulationOutput$fittedPredictedResponse,
               simulationOutput$scaledResidualsNormal,
               xlab = "Predicted value", ylab = "Standardized residual")+
           geom_abline(intercept = 0, slope = 0, colour = "red") +ggtitle(spsel)
    p2 <- plot.bin2d( df$size,
               simulationOutput$scaledResidualsNormal,
               xlab = "size", ylab = "Standardized residual")+
           geom_abline(intercept = 0, slope = 0, colour = "red")
    p3 <- plot.bin2d( df$BATOTcomp,
               simulationOutput$scaledResidualsNormal,
               xlab = "BATOTcomp", ylab = "Standardized residual")+
           geom_abline(intercept = 0, slope = 0, colour = "red")
    p4 <- plot.bin2d( df[[vars[1]]],
               simulationOutput$scaledResidualsNormal,
               xlab = vars[1], ylab = "Standardized residual")+
           geom_abline(intercept = 0, slope = 0, colour = "red")
    p5 <- plot.bin2d( df[[vars[2]]],
               simulationOutput$scaledResidualsNormal,
               xlab = vars[2], ylab = "Standardized residual")+
           geom_abline(intercept = 0, slope = 0, colour = "red")
    png(file.path("figures", paste0(gsub(" ", "_", spsel),
                                    "fit_survival.png")),
         width = 960, height = 480)
    multiplot(p1, p2, p3, p4, p5 ,cols = 2)
    dev.off()
  }
}


plot_survival_size_climate_inter<- function(spsel = "Pinus sylvestris",
                                         data,
                                         vars = c("waib", "sgddb"),
                                         inters = c("size:waib", "logsize:waib",
                                                    "size:sgddb", "logsize:sgddb",
                                                    "BATOTcomp:waib", "BATOTcomp:sgddb")) {
  df <- format_data_survival(data, spsel)
  res <- survival_size_climate(df, c(vars, inters))
  vars <- gsub(",[0-9]\\)", "", gsub("poly\\(", "", vars))
  inters1 <- gsub(":.*","", inters)
  inters2 <- gsub(".*:","", inters)
  inters1b <- inters1[inters1 != "logsize"]
  inters2b <- inters2[inters1 != "logsize"]
  fixed_vars <- c("size", "logsize", "BATOTcomp", vars)
  list_inter_plots <- vector("list")
  par(mfrow = c(4, 4))
  for (i in seq_len(length.out = length(inters1b))){
  tryCatch({
      plot_surv_bin_pred_inter(var = inters1b[i], vary = inters2b[i], vars = fixed_vars,
                             df, res, legend_TF = TRUE, spsel = spsel)
  }, error=function(e){})
  }

}


survival_size_pred_all<- function(sps, data) {
 pdf(file.path("figures", "all_pred_survival_size.pdf"))
 l <- lapply(sps$sp, plot_survival_size, data)
 dev.off()
}


survival_size_climate_pred_all_waib<- function(sps, data) {
 pdf(file.path("figures", "all_pred_survival_waib.pdf"))
 l <- lapply(sps$sp, plot_survival_size_climate, data)
 dev.off()
}

survival_size_climate_pred_all_waib_inter<- function(sps, data) {
 pdf(file.path("figures", "all_pred_survival_waib_inter.pdf"))
 l <- lapply(sps$sp, plot_survival_size_climate, data, vars = c("sgddb", "waib",
                                                               "size:waib", "logsize:waib",
                                                               "size:sgddb", "logsize:sgddb"))
 dev.off()
}

survival_size_climate_pred_all_waib_inter2<- function(sps, data) {
 pdf(file.path("figures", "all_pred_survival_waib_inter2.pdf"))
 l <- lapply(sps$sp, plot_survival_size_climate, data,vars = c("sgddb", "waib",
                                                               "size:waib", "logsize:waib",
                                                               "size:sgddb", "logsize:sgddb",
                                                               "BATOTcomp:waib", "BATOTcomp:sgddb"))
 dev.off()
}


survival_size_climate_eval_all_waib<- function(sps, data) {
 pdf(file.path("figures", "all_eval_survival_waib.pdf"))
 l <- lapply(sps$sp, plot_survival_eval, data)
 dev.off()
}

survival_size_climate_eval_inter_all_waib<- function(sps, data) {
 pdf(file.path("figures", "all_eval_inter_survival_waib.pdf"))
 l <- lapply(sps$sp, plot_survival_eval, data,
             vars = c("waib", "sgddb","size:waib", "logsize:waib",
                      "size:sgddb", "logsize:sgddb"))
 dev.off()
}

survival_size_climate_eval_inter2_all_waib<- function(sps, data) {
 pdf(file.path("figures", "all_eval_inter2_survival_waib.pdf"))
 l <- lapply(sps$sp, plot_survival_eval, data,
             vars = c("waib", "sgddb","size:waib", "logsize:waib",
                      "size:sgddb", "logsize:sgddb",
                      "BATOTcomp:waib", "BATOTcomp:sgddb"))
 dev.off()
}



survival_size_climate_pred_all_waib_inter<- function(sps, data) {
 pdf(file.path("figures", "all_pred_survival_waib_inter.pdf"))
 l <- lapply(sps$sp, plot_survival_size_climate_inter, data)
 dev.off()
}

survival_size_climate_pred_all_waib_inter_clim<- function(sps, data) {
 pdf(file.path("figures", "all_pred_survival_waib_inter_clim.pdf"))
 l <- lapply(sps$sp, plot_survival_size_climate_inter, data,
             inters = c("size:waib", "logsize:waib",
                        "size:sgddb", "logsize:sgddb",
                        "BATOTcomp:waib", "BATOTcomp:sgddb", "waib:sgddb"))
 dev.off()
}

survival_size_climate_pred_all_logwai<- function(sps, data) {
 pdf(file.path("figures", "all_pred_survival_logwai.pdf"))
 l <- lapply(sps$sp, plot_survival_size_climate, data, vars = c("logwai", "logsgdd"))
 dev.off()
}

survival_size_climate_pred_all_poly3<- function(sps, data) {
 pdf(file.path("figures", "all_pred_survival_poly3.pdf"))
 l <- lapply(sps$sp, plot_survival_size_climate, data, vars = c("poly(wai,3)", "poly(sgdd,3)"))
 dev.off()
}


survival_size_climate_fit_all<- function(sps, data) {
 l <- lapply(sps$sp, plot_survival_size_climate, data, plot_type = "fit_surv")
    system("convert  figures/*fit_survival.png figures/all_fit_survival.pdf")
    system("rm  figures/*fit_survival.png")
}


