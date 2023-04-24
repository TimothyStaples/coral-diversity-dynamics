model.localities1 <- function(data,
                             response,
                             knots,
                             gam.family,
                             date.diff=FALSE){
  
  require(mgcv)
  
  data$response <- data[,response]
  data$scaled.date <- scale(data$pred.date)
  data$dummy <- 1
  
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
    
    model <- gam(response ~ s(scaled.date, k=knots, by=locality) + locality + difference,
                 data=data, family = gam.family, method="REML")
  } else {
    pred.date.diff = NA
    model <- gam(response ~ s(scaled.date, k=knots, by=locality) + locality,
                 data=data, family = gam.family,
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
  
  pred.df <- do.call("rbind", lapply(unique(data$locality), 
                                     function(loc){
                                       print(loc)
                                       temp.data <- data[as.character(data$locality) == loc, ]
                                       
                                       date.logic <- pred.date.seq > min(temp.data$scaled.date) &
                                         pred.date.seq < max(temp.data$scaled.date)
                                       
                                       data.frame(scaled.date = pred.date.seq[date.logic],
                                                  raw.date = raw.date.seq[date.logic],
                                                  locality=loc,
                                                  site=temp.data$site[1],
                                                  difference=pred.date.diff)
                                     }))
  
  trend.preds <- cbind(pred.df,
                       predict(model, 
                               newdata=pred.df, se.fit=TRUE, type="link"))
  
  trend.preds$upper <- trend.preds$fit + 1.96*trend.preds$se.fit
  trend.preds$lower <- trend.preds$fit - 1.96*trend.preds$se.fit
  
  # Predict slope of trends over time
  
  # create prediction matrix offset from initial matrix by small amount
  
  # pred.df.eps <- pred.df
  # eps <- (max(pred.df$scaled.date)-min(pred.df$scaled.date))*1e-5 # finite difference interval
  # pred.df.eps$scaled.date <- pred.df.eps$scaled.date + eps
  # 
  # # create linear predictor matrices
  # X0<-predict(model, newdata=pred.df, type="lpmatrix")
  # X1<-predict(model, newdata=pred.df.eps, type="lpmatrix")
  # 
  # Xp <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives
  # 
  # # get slope estimates and CIs across spline, and determine if 95% CIs cross 0.
  # A <- model
  # slopes<-do.call("rbind", lapply(unique(pred.df$locality), function(x){
  #   
  #   group.cols<-which(grepl(paste0(x), colnames(Xp)))
  #   
  #   Xi <- Xp*0 
  #   Xi[,group.cols] <- Xp[,group.cols] ## Xi%*%coef(b) = smooth deriv i
  #   df <- Xi%*%coef(A)              ## ith smooth derivative 
  #   df.sd <- rowSums(Xi%*%A$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  #   
  #   return(data.frame(slope = df,
  #                     slope.se = df.sd,
  #                     slope.upper = df + 1.96*df.sd,
  #                     slope.lower = df - 1.96*df.sd)[which(pred.df$locality == x), ])
  #   
  # }))
  
  preds <- trend.preds
  
  return(list(model, preds, data))
}