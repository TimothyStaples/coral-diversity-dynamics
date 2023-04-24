model.localities <- function(data,
                             response,
                             knots,
                             family){
  
require(mgcv)
require(gamm4)
  
data$response <- data[,response]
data$scaled.date <- scale(data$pred.date)

error.trap<-tryCatch(
  gamm(response ~ s(scaled.date, k=knots, by=locality) + locality + site,
       random=list(site = ~1), data=data, family=family),
  error=function(e) e)

if(!inherits(error.trap, "error")){

model <- gamm(response ~ s(scaled.date, k=knots, by=locality) + locality,
              random=list(site = ~1), data=data, family=family)
} else {

model <- list(gam = gam(response ~ s(scaled.date, k=knots, by=locality) + locality + site, 
                        data=data, family=family))
  
}

# Predict trends over time. Predictions have to use the same time-slices as the
# overall model so comparisons make sense
pred.date.seq <- seq(min(data$scaled.date),
                     max(data$scaled.date),
                     length.out=200)
raw.date.seq <- seq(min(data$pred.date),
                     max(data$pred.date),
                     length.out=200)

pred.df <- do.call("rbind", lapply(unique(data$locality), 
                                   function(loc){
  print(loc)
  temp.data <- data[as.character(data$locality) == loc, ]
  
  date.logic <- pred.date.seq > min(temp.data$scaled.date) &
                pred.date.seq < max(temp.data$scaled.date)
  
  data.frame(scaled.date = pred.date.seq[date.logic],
             raw.date = raw.date.seq[date.logic],
             locality=loc,
             site = temp.data$site[1])
}))

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
slopes<-do.call("rbind", lapply(unique(pred.df$locality), function(x){
  
  group.cols<-which(grepl(paste0(x), colnames(Xp)))
  
  Xi <- Xp*0 
  Xi[,group.cols] <- Xp[,group.cols] ## Xi%*%coef(b) = smooth deriv i
  df <- Xi%*%coef(A$gam)              ## ith smooth derivative 
  df.sd <- rowSums(Xi%*%A$gam$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  
  return(data.frame(slope = df,
               slope.se = df.sd,
               slope.upper = df + 1.96*df.sd,
               slope.lower = df - 1.96*df.sd)[which(pred.df$locality == x), ])
  
}))

preds <- cbind(trend.preds, slopes)

return(list(model, preds, data))

}