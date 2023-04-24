model.overall.trend1 <- function(data,
                             response,
                             knots,
                             gam.family="gaussian",
                             date.diff=FALSE){

require(mgcv)

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
  
  model <-   gamm(response ~ s(scaled.date, k=knots) + difference,
                  random=list(site = ~1, locality = ~1), data=data, 
                  family=gam.family,
                  method="REML")
} else {
pred.date.diff = NA
model <-   gamm(response ~ s(scaled.date, k=knots),
                random=list(site = ~1, locality = ~1), data=data, 
                family=gam.family,
                method="REML")
}

# Predict trends over time
pred.df <-  data.frame(scaled.date = seq(min(data$scaled.date),
                                         max(data$scaled.date),
                                         length.out=200),
                       raw.date = seq(min(data$pred.date),
                                     max(data$pred.date),
                                     length.out=200),
                       locality=NA, 
                       site=NA,
                       difference=pred.date.diff)

trend.preds <- cbind(pred.df,
                     predict(model$gam, 
                             newdata=pred.df, se.fit=TRUE, type="link"))
trend.preds$upper <- trend.preds$fit + 1.96*trend.preds$se.fit
trend.preds$lower <- trend.preds$fit - 1.96*trend.preds$se.fit

# Predict slope of trends over time
# 
# # create prediction matrix offset from initial matrix by small amount
# 
# pred.df.eps <- pred.df
# eps <- (max(pred.df$scaled.date)-min(pred.df$scaled.date))*1e-5 # finite difference interval
# pred.df.eps$scaled.date <- pred.df.eps$scaled.date + eps
# 
# # create linear predictor matrices
# X0<-predict(model$gam, newdata=pred.df, type="lpmatrix")
# X1<-predict(model$gam, newdata=pred.df.eps, type="lpmatrix")
# 
# Xp <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives
# 
# # get slope estimates and CIs across spline, and determine if 95% CIs cross 0.
# A <- model
# 
# Xi <- Xp*0 
# Xi[,-1] <- Xp[,-1] ## Xi%*%coef(b) = smooth deriv i
# df <- Xi%*%coef(A$gam)              ## ith smooth derivative 
# df.sd <- rowSums(Xi%*%A$gam$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
#   
# slopes <- data.frame(slope = df,
#                slope.se = df.sd,
#                slope.upper = df + 1.96*df.sd,
#                slope.lower = df - 1.96*df.sd)

preds <- trend.preds

return(list(model, preds, data))

}