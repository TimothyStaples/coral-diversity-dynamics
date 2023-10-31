mvGamsFreq <- function(dataList, primaryvarList, secondaryvarList, dateDiff = FALSE){
  
  lapply(1:length(dataList), function(n){
    
    print(n)
    model.data <- dataList[[n]]
    model.data$primary <- model.data[,primaryvarList[n]]
    model.data$secondary <- model.data[,secondaryvarList[n]]
    model.data$date.diff <- log(model.data$date.diff)
    
    # dummy variable to switch off random effect smoothers for prediction
    model.data$dummy = 1
    
    model.data = model.data[complete.cases(model.data[,c("primary", "secondary")]),]
    
    # we need to include difference in dates as covariate for beta diversity models
    if(!dateDiff){
      # region level model
      region.m <- gam(list(primary ~ s(pred.date, k=9, bs="ps") + s(site, locality, bs="re", by=dummy),
                           secondary ~ s(pred.date, k=9, bs="ps") + s(site, locality, bs="re", by=dummy)), 
                      family=mvn(d=2), data=model.data, method="REML")
      
      # site level model (called 'locality' in plot)
      site.m <- gam(list(primary ~ s(pred.date, by=site,  bs="ps", k=9) + site + s(locality, bs="re", by=dummy),
                         secondary ~ s(pred.date, by=site,  bs="ps", k=9) + site + s(locality, bs="re", by=dummy)), 
                    family=mvn(d=2), data=model.data, method="REML")
      
      # locality level model (called 'site' in plot)
      locality.m <- gam(list(primary ~ s(pred.date, by=locality, bs="ps", k=4) + locality,
                             secondary ~ s(pred.date, by=locality, bs="ps", k=4) + locality), 
                        family=mvn(d=2), data=model.data, method="REML")
      
      
        } else {
      
      # region level model
      region.m <- gam(list(primary ~ s(pred.date, k=9, bs="ps") + date.diff + s(site, locality, bs="re", by=dummy),
                           secondary ~ s(pred.date, k=9, bs="ps") + date.diff + s(site, locality, bs="re", by=dummy)), 
                      family=mvn(d=2), data=model.data, method="REML")
      
      # site level model (called 'locality' in plot)
      site.m <- gam(list(primary ~ s(pred.date, by=site, k=9, bs="ps") + site + date.diff +  s(locality, bs="re", by=dummy),
                         secondary ~ s(pred.date, by=site, k=9, bs="ps") + site + date.diff + s(locality, bs="re", by=dummy)), 
                    family=mvn(d=2), data=model.data, method="REML")
      
      # locality level model (called 'site' in plot)
      locality.m <- gam(list(primary ~ s(pred.date, by=locality, k=4, bs="ps") + locality + date.diff,
                             secondary ~ s(pred.date, by=locality, k=4, bs="ps")  + locality + date.diff), 
                        family=mvn(d=2), data=model.data, method="REML")
      
    }
    
    # predictions
    pred.df <- data.frame(pred.date = seq(round(min(model.data$pred.date),-1), 
                                          max(model.data$pred.date), 
                                          10),
                          date.diff = mean(model.data$date.diff, na.rm=TRUE),
                          top.date.diff = mean(model.data$top.date.diff, na.rm=TRUE),
                          site="Kanomi", locality = "Bonah River", dummy=0)
    
    region.pred.df <- cbind(pred.df, predict(region.m, newdata=pred.df, se.fit=TRUE))
    
    site.pred.df <- data.frame(pred.date = rep(pred.df$pred.date, 3),
                               date.diff = mean(model.data$date.diff, na.rm=TRUE),
                               site=rep(levels(model.data$site), each=nrow(pred.df)), 
                               locality = "Bonah River", dummy=0)
    
    site.pred.df <- cbind(site.pred.df,
                          predict(site.m, newdata=site.pred.df, se.fit=TRUE))
    
    locality.pred.df <- data.frame(pred.date = rep(pred.df$pred.date, length(levels(dataList[[1]]$locality))),
                                   date.diff = mean(model.data$date.diff, na.rm=TRUE),
                                   locality=rep(levels(model.data$locality), each=nrow(pred.df)))
    locality.pred.df <- cbind(locality.pred.df,
                              predict(locality.m, newdata=locality.pred.df, se.fit=TRUE))
    
    # temporal autocorrelation test
    # residR <- resid(region.m)
    # siteR <- resid(site.m)
    # localityR <- resid(locality.m)
    # taP <- lmtest::dwtest(residR[,1] ~ 1, order.by = region.m$model$pred.date)
    # taS <- lmtest::dwtest(residR[,2] ~ 1, order.by = region.m$model$pred.date)
    # 
    # taP <- lmtest::dwtest(siteR[,1] ~ 1, order.by = site.m$model$pred.date)
    # taS <- lmtest::dwtest(siteR[,2] ~ 1, order.by = site.m$model$pred.date)
    # 
    # 
    # taP <- lmtest::dwtest(localityR[,1] ~ 1, order.by = locality.m$model$pred.date)
    # taS <- lmtest::dwtest(localityR[,2] ~ 1, order.by = locality.m$model$pred.date)
 

    list(models = list(region.m, site.m, locality.m),
         preds = list(region.pred.df, site.pred.df, locality.pred.df))
    
  })
  
}
