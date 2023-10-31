mvGamsFreqRegion <- function(dataList, primaryvarList, secondaryvarList, dateDiff = FALSE){
  
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
      region.m <- gam(list(primary ~ s(pred.date, k=9, bs="cs") + s(site, locality, bs="re", by=dummy),
                           secondary ~ s(pred.date, k=9, bs="cs") + s(site, locality, bs="re", by=dummy)), 
                      family=mvn(d=2), data=model.data, method="REML")
      
      
    } else {
      
      # region level model
      region.m <- gam(list(primary ~ s(pred.date, k=9, bs="cs") + date.diff + s(site, locality, bs="re", by=dummy),
                           secondary ~ s(pred.date, k=9, bs="cs") + date.diff + s(site, locality, bs="re", by=dummy)), 
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
    
    
    list(models = list(region.m),
         preds = list(region.pred.df))
    
  })
  
}
