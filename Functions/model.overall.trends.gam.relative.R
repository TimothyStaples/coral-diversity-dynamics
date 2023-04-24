model.overall.trend.relative <- function(data,
                             response,
                             knots,
                             gam.family="gaussian",
                             date.diff=FALSE){

require(mgcv)

data$response <- data[,response]
data$dummy <- 1
data$scaled.date <- scale(data$top.date.diff)

if(gam.family=="betar"){
  data$raw.response <- data$response
  data$response <- beta.tr(data$response)
}

if(date.diff){
  
  if(response == "top.diss"){
    
    date.diff.scaled <- scale(data$top.date.diff)
    
  } else {
  
    date.diff.scaled <- scale(data$date.diff)  
    
  }
  
  data$difference <- as.vector(date.diff.scaled)
  attr(data$difference, "scaled:center") <- NULL
  attr(data$difference, "scaled:scale") <- NULL
  
  model <-   gamm(response ~ s(scaled.date, k=knots) + difference,
                  random=list(site = ~1, locality = ~1), data=data, 
                  family=gam.family,
                  method="REML")
} else {
model <-   gamm(response ~ s(scaled.date, k=knots),
                random=list(site = ~1, locality = ~1), data=data, 
                family=gam.family,
                method="REML")
}

# Predict trends over time
pred.df <-  data.frame(scaled.date = seq(min(data$scaled.date),
                                         max(data$scaled.date),
                                         length.out=200),
                       raw.date = seq(min(data$top.date.diff),
                                     max(data$top.date.diff),
                                     length.out=200),
                       locality=NA, 
                       site=NA,
                       difference=0)

trend.preds <- cbind(pred.df,
                     predict(model$gam, 
                             newdata=pred.df, se.fit=TRUE, type="link"))
trend.preds$upper <- trend.preds$fit + 1.96*trend.preds$se.fit
trend.preds$lower <- trend.preds$fit - 1.96*trend.preds$se.fit

preds <- cbind(trend.preds)

return(list(model, preds, data))

}