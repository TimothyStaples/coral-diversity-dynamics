model.localities.relative <- function(data,
                             response,
                             knots,
                             gam.family,
                             date.diff=FALSE){
  
require(mgcv)

data$response <- data[,response]
data$dummy <- 1
data$scaled.date <- scale(data$top.date.diff)

if(date.diff){
  
data$difference <- data$date.diff
  
model <- gam(response ~ s(scaled.date, k=knots, by=locality) + locality + 
               difference +
               s(site, bs="re", by=dummy), data=data, family = gam.family,
             method="REML")
} else {
  model <- gam(response ~ s(scaled.date, k=knots, by=locality) + locality + 
                 s(site, bs="re", by=dummy), data=data, family = gam.family,
               method="REML")
}


# Predict trends over time. Predictions have to use the same time-slices as the
# overall model so comparisons make sense
pred.date.seq <- seq(min(data$scaled.date),
                     max(data$scaled.date),
                     length.out=200)
raw.date.seq <- seq(min(data$top.date.diff),
                     max(data$top.date.diff),
                     length.out=200)

pred.df <- do.call("rbind", lapply(unique(data$locality), 
                                   function(loc){

  temp.data <- data[as.character(data$locality) == loc, ]
  
  date.logic <- pred.date.seq > min(temp.data$scaled.date) &
                pred.date.seq < max(temp.data$scaled.date)
  
  data.frame(scaled.date = pred.date.seq[date.logic],
             raw.date = raw.date.seq[date.logic],
             locality=loc,
             site=temp.data$site[1],
             dummy=0,
             difference=0)
}))

trend.preds <- cbind(pred.df,
                     predict(model, 
                             newdata=pred.df, se.fit=TRUE, type="link"))

trend.preds$upper <- trend.preds$fit + 1.96*trend.preds$se.fit
trend.preds$lower <- trend.preds$fit - 1.96*trend.preds$se.fit

preds <- trend.preds

return(list(model, preds, data))

}