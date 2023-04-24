model.overall.trend <- function(data,
                             response,
                             knots,
                             family="gaussian"){

require(mgcv)
require(gamm4)
  
data$response <- data[,response]

if(family=="betar"){
  data$raw.response <- data$response
  data$response <- beta.tr(data$response)
}

data$scaled.date <- scale(data$pred.date)

error.trap<-tryCatch(
  gamm(response ~ s(scaled.date, k=knots),
       random=list(site = ~1, locality = ~1), data=data, family=family),
  error=function(e) e)

if(!inherits(error.trap, "error")){
  
  model <- gamm(response ~ s(scaled.date, k=knots),
                random=list(site = ~1, locality = ~1), data=data)
} else {
  
  model <- list(gam = gam(response ~ s(scaled.date, k=knots), 
                          data=data))
  
}

# Predict trends over time
pred.df <-  data.frame(scaled.date = seq(min(data$scaled.date),
                                         max(data$scaled.date),
                                         length.out=200),
                       raw.date = seq(min(data$pred.date),
                                     max(data$pred.date),
                                     length.out=200),
                       locality=NA, 
                       site=NA)

trend.preds <- cbind(pred.df,
                     predict(model$gam, 
                             newdata=pred.df, se.fit=TRUE, type="link"))
trend.preds$upper <- trend.preds$fit + 1.96*trend.preds$se.fit
trend.preds$lower <- trend.preds$fit - 1.96*trend.preds$se.fit

# Predict slope of trends over time

# create prediction matrix offset from initial matrix by small amount

pred.df.eps <- pred.df
eps <- (max(pred.df$scaled.date)-min(pred.df$scaled.date))*1e-5 # finite difference interval
pred.df.eps$scaled.date <- pred.df.eps$scaled.date + eps

# create linear predictor matrices
X0<-predict(model$gam, newdata=pred.df, type="lpmatrix")
X1<-predict(model$gam, newdata=pred.df.eps, type="lpmatrix")

Xp <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives

# get slope estimates and CIs across spline, and determine if 95% CIs cross 0.
A <- model

Xi <- Xp*0 
Xi[,-1] <- Xp[,-1] ## Xi%*%coef(b) = smooth deriv i
df <- Xi%*%coef(A$gam)              ## ith smooth derivative 
df.sd <- rowSums(Xi%*%A$gam$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  
slopes <- data.frame(slope = df,
               slope.se = df.sd,
               slope.upper = df + 1.96*df.sd,
               slope.lower = df - 1.96*df.sd)

preds <- cbind(trend.preds, slopes)



return(list(model, preds, data))

}