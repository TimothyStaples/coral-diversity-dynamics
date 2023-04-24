model.sites <- function(data,
                       response,
                       knots,
                       gam.family="gaussian",
                       date.diff=FALSE){

require(mgcv)
require(gamm4)
  
data$response <- data[,response]

if(gam.family=="betar"){
  data$raw.response <- data$response
  data$response <- beta.tr(data$response)
}

data$scaled.date <- scale(data$pred.date)

if(date.diff){
  
  if(response == "top.diss"){
    
    date.diff.scaled <- scale(data$top.date.diff)
    pred.date.diff = 0 
  } else {
    
    date.diff.scaled <- scale(data$date.diff)  
    pred.date.diff = (100 - attr(date.diff.scaled, "scaled:center")) /
      attr(date.diff.scaled, "scaled:scale")
  }
  
  data$difference <- as.vector(date.diff.scaled)
  attr(data$difference, "scaled:center") <- NULL
  attr(data$difference, "scaled:scale") <- NULL
  
  model <-   gamm(response ~ s(scaled.date, k=knots, by=site) + site + difference,
                  random=list(locality = ~1), data=data, 
                  family=gam.family,
                  method="REML")
} else {
  pred.date.diff = NA
  model <-   gamm(response ~ s(scaled.date, k=knots, by=site) + site,
                  random=list(locality = ~1), data=data, 
                  family=gam.family,
                  method="REML")
}

# Predict trends over time. Predictions have to use the same time-slices as the
# overall model so comparisons make sense
pred.date.seq <- seq(min(data$scaled.date),
                     max(data$scaled.date),
                     length.out=200)
raw.date.seq <- seq(min(data$pred.date),
                    max(data$pred.date),
                    length.out=200)

# Predict trends over time
pred.df <- do.call("rbind", lapply(unique(data$site), 
                                   function(site){
                                     print(site)
                                     temp.data <- data[as.character(data$site) == site, ]
                                     
                                     date.logic <- pred.date.seq > min(temp.data$scaled.date) &
                                       pred.date.seq < max(temp.data$scaled.date)
                                     
                                     data.frame(scaled.date = pred.date.seq[date.logic],
                                                raw.date = raw.date.seq[date.logic],
                                                site=site,
                                                locality = NA,
                                                difference=pred.date.diff)
                                   }))

trend.preds <- cbind(pred.df,
                     predict(model$gam, 
                             newdata=pred.df, se.fit=TRUE, type="link"))
trend.preds$upper <- trend.preds$fit + 1.96*trend.preds$se.fit
trend.preds$lower <- trend.preds$fit - 1.96*trend.preds$se.fit

preds <- trend.preds



return(list(model, preds, data))

}