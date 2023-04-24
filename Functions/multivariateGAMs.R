mvGams <- function(dataList, primaryvarList, secondaryvarList, 
                   beta=FALSE, modelPath, modelIter=1000, modelWarmup=500){
  
  lapply(1:length(dataList), function(n){
    
    print(n)
    model.data <- dataList[[n]]
    model.data$primary <- model.data[,primaryvarList[n]]
    model.data$secondary <- model.data[,secondaryvarList[n]]
    
    print(paste0("Running model ", n, ": ", primaryvarList[n], " & ", secondaryvarList[n]))
    
    # # dummy variable to switch off random effect smoothers for prediction
    # model.data$dummy = 1
    # 
    # we need to include difference in dates as covariate for beta diversity models
    if(!beta){
      
      # region level model
      region.m <- brm(mvbind(primary, secondary) ~ s(pred.date) + (1|site/locality),
                           data=model.data, cores = 4, chains = 4,
                      iter = modelIter, warmup=modelWarmup, control = list(adapt_delta = 0.99),
                      prior = c(#set_prior("cauchy(0,2)", class = "sd", group = "site", coef = "Intercept"),
                                #set_prior("cauchy(0,2)", class = "sd", group = "site:locality", coef = "Intercept"),
                                set_prior("normal(0,3)", class = "Intercept", resp="primary"),
                                set_prior("normal(0,3)", class = "Intercept", resp="secondary"),
                                set_prior("normal(0,3)", class = "b", resp="primary"),
                                set_prior("normal(0,3)", class = "b", resp="secondary")))
      
      # region.m <- gam(list(primary ~ s(pred.date, k=9) + s(site, locality, bs="re", by=dummy),
      #                      secondary ~ s(pred.date, k=9) + s(site, locality, bs="re", by=dummy)), 
      #                 family=mvn(d=2), data=model.data, method="REML")
      
      # site level model (called 'locality' in plot)
      site.m  <- brm(mvbind(primary, secondary) ~ s(pred.date, by = site) + site + (1|locality),
                       data=model.data, cores = 4, chains = 4,
                     iter = modelIter, warmup=modelWarmup, control = list(adapt_delta = 0.99),
                     prior = c(#set_prior("cauchy(0,2)", class = "sd", group = "site", coef = "Intercept"),
                         #set_prior("cauchy(0,2)", class = "sd", group = "site:locality", coef = "Intercept"),
                         set_prior("normal(0,3)", class = "Intercept", resp="primary"),
                         set_prior("normal(0,3)", class = "Intercept", resp="secondary"),
                         set_prior("normal(0,3)", class = "b", resp="primary"),
                         set_prior("normal(0,3)", class = "b", resp="secondary")))
        
      # site.m <- gam(list(primary ~ s(pred.date, by=site, k=9) + site + s(locality, bs="re", by=dummy),
      #                    secondary ~ s(pred.date, by=site, k=9) + site + s(locality, bs="re", by=dummy)), 
      #               family=mvn(d=2), data=model.data, method="REML")
      
      # locality level model (called 'site' in plot)
      locality.m  <- brm(mvbind(primary, secondary) ~ s(pred.date, by = locality) + locality,
                     data=model.data, cores = 4, chains = 4,
                     iter = modelIter, warmup=modelWarmup, control = list(adapt_delta = 0.99),
                     prior = c(#set_prior("cauchy(0,2)", class = "sd", group = "site", coef = "Intercept"),
                       #set_prior("cauchy(0,2)", class = "sd", group = "site:locality", coef = "Intercept"),
                       set_prior("normal(0,3)", class = "Intercept", resp="primary"),
                       set_prior("normal(0,3)", class = "Intercept", resp="secondary"),
                       set_prior("normal(0,3)", class = "b", resp="primary"),
                       set_prior("normal(0,3)", class = "b", resp="secondary")))
      # 
      # locality.m <- gam(list(primary ~ s(pred.date, by=locality, k=4) + locality,
      #                        secondary ~ s(pred.date, by=locality, k=4) + locality), 
      #                   family=mvn(d=2), data=model.data, method="REML")
    } else {
      
      # region level model
      region.m <- brm(mvbind(primary, secondary) ~ s(pred.date) + date.diff + (1|site/locality),
                      data=model.data, cores = 4, chains = 4,
                      iter = modelIter, warmup=modelWarmup, control = list(adapt_delta = 0.99),
                      prior = c(#set_prior("cauchy(0,2)", class = "sd", group = "site", coef = "Intercept"),
                        #set_prior("cauchy(0,2)", class = "sd", group = "site:locality", coef = "Intercept"),
                        set_prior("normal(0,3)", class = "Intercept", resp="primary"),
                        set_prior("normal(0,3)", class = "Intercept", resp="secondary"),
                        set_prior("normal(0,3)", class = "b", resp="primary"),
                        set_prior("normal(0,3)", class = "b", resp="secondary")))
      # 
      # region.m <- gam(list(primary ~ s(pred.date, k=9) + date.diff + s(site, locality, bs="re", by=dummy),
      #                      secondary ~ s(pred.date, k=9) + date.diff + s(site, locality, bs="re", by=dummy)), 
      #                 family=mvn(d=2), data=model.data, method="REML")
      
      # site level model (called 'locality' in plot)
      site.m  <- brm(mvbind(primary, secondary) ~ s(pred.date, by = site) + date.diff + site + (1|locality),
                     data=model.data, cores = 4, chains = 4,
                     iter = modelIter, warmup=modelWarmup, control = list(adapt_delta = 0.99),
                     prior = c(#set_prior("cauchy(0,2)", class = "sd", group = "site", coef = "Intercept"),
                       #set_prior("cauchy(0,2)", class = "sd", group = "site:locality", coef = "Intercept"),
                       set_prior("normal(0,3)", class = "Intercept", resp="primary"),
                       set_prior("normal(0,3)", class = "Intercept", resp="secondary"),
                       set_prior("normal(0,3)", class = "b", resp="primary"),
                       set_prior("normal(0,3)", class = "b", resp="secondary")))
      # site.m <- gam(list(primary ~ s(pred.date, by=site, k=9) + site + date.diff +  s(locality, bs="re", by=dummy),
      #                    secondary ~ s(pred.date, by=site, k=9) + site + date.diff + s(locality, bs="re", by=dummy)), 
      #               family=mvn(d=2), data=model.data, method="REML")
      # 
      # locality level model (called 'site' in plot)
      locality.m  <- brm(mvbind(primary, secondary) ~ s(pred.date, by = locality) + date.diff + locality,
                         data=model.data, cores = 4, chains = 4,
                         iter = modelIter, warmup=modelWarmup, control = list(adapt_delta = 0.99),
                         prior = c(#set_prior("cauchy(0,2)", class = "sd", group = "site", coef = "Intercept"),
                           #set_prior("cauchy(0,2)", class = "sd", group = "site:locality", coef = "Intercept"),
                           set_prior("normal(0,3)", class = "Intercept", resp="primary"),
                           set_prior("normal(0,3)", class = "Intercept", resp="secondary"),
                           set_prior("normal(0,3)", class = "b", resp="primary"),
                           set_prior("normal(0,3)", class = "b", resp="secondary")))
      
      # locality.m <- gam(list(primary ~ s(pred.date, by=locality, k=4) + locality + date.diff,
      #                        secondary ~ s(pred.date, by=locality, k=4)  + locality + date.diff), 
      #                   family=mvn(d=2), data=model.data, method="REML")
      # 
    }
    
    # predictions
    pred.df <- data.frame(pred.date = seq(round(min(model.data$pred.date),-1), 
                                          max(model.data$pred.date), 
                                          10),
                          date.diff = mean(dataList[[2]]$date.diff, na.rm=TRUE),
                          site=NA, locality = NA)
    
    region.pred.df <- pred.df
    
    region.preds <- predict(region.m, newdata=region.pred.df, re.formula=NA, allow_new_levels=TRUE)
    
    site.pred.df <- data.frame(pred.date = rep(pred.df$pred.date, 3),
                               date.diff = mean(dataList[[2]]$date.diff, na.rm=TRUE),
                               site=rep(levels(model.data$site), each=nrow(pred.df)), 
                               locality = NA)
    
    site.preds <- predict(site.m, newdata=site.pred.df, re.formula=NA, allow_new_levels=TRUE)
    
    locality.pred.df <- data.frame(pred.date = rep(pred.df$pred.date, length(levels(dataList[[1]]$locality))),
                                   date.diff = mean(dataList[[2]]$date.diff, na.rm=TRUE),
                                   locality=rep(levels(model.data$locality), each=nrow(pred.df)))
    
    locality.preds <- predict(locality.m, newdata=locality.pred.df, re.formula=NA)
    
    # temporal autocorrelation test
    # residR <- resid(region.m)
    # siteR <- resid(site.m)
    # localityR <- resid(locality.m)
    # taP <- lmtest::dwtest(residR[,1] ~ 1, order.by = region.m$model$pred.date)
    # taS <- lmtest::dwtest(residR[,2] ~ 1, order.by = region.m$model$pred.date)
    # print(taP$p.value)
    # print(taS$p.value)
    # 
    # taP <- lmtest::dwtest(siteR[,1] ~ 1, order.by = site.m$model$pred.date)
    # taS <- lmtest::dwtest(siteR[,2] ~ 1, order.by = site.m$model$pred.date)
    # print(taP$p.value)
    # print(taS$p.value)
    # 
    # 
    # taP <- lmtest::dwtest(localityR[,1] ~ 1, order.by = locality.m$model$pred.date)
    # taS <- lmtest::dwtest(localityR[,2] ~ 1, order.by = locality.m$model$pred.date)
    # print(taP$p.value)
    # print(taS$p.value)
    
    tempList = list(models = list(region.m, site.m, locality.m),
         preds = list(list(region.pred.df, region.preds),
                      list(site.pred.df, site.preds),
                      list(locality.pred.df, locality.preds)))
    
    saveRDS(tempList, paste0(modelPath, "/", primaryvarList[n],"-", secondaryvarList[n],".rds"))
    
    return(tempList)
  })
  
}