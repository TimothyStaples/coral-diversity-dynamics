beta.model <- function(data, 
                       test.model = FALSE, 
                       knots = -1,
                       name,
                       ts.col){
  
  require(mgcv)
  
  data$dummy = 1 # to switch off random effects for prediction
  
  # six models
  # - top dissim absolute time
  # - top dissim relative time
  # - bottom dissim absolute time
  # - bottom dissim relative time
  # - seq dissim absolute time
  # - seq dissim relative time
  model.list <- lapply(1:6, function(n){
    
    print(n)
    
    response <- c("t.beta", "b.beta", "s.beta")[(n+1) %/% 2]
    time.var <- rep(c("age", "age.from.top"), 3)[n]
    
    sub.data <- data
    sub.data$response <- beta.tr(sub.data[,response])
    sub.data$time.var <- sub.data[,time.var]
    sub.data <- sub.data[!is.na(sub.data$response),]
    
    overall.model <- gam(response ~ s(time.var, bs="cr", k=knots) + s(ts, bs="re", by=dummy),
                     family=betar(),
                     data=sub.data)
  
    ts.model <- gam(response ~ s(time.var, bs="cr", k=knots, by=ts) + ts,
                    family=betar(),
                    data=sub.data)
    
    return(list(df = sub.data,
                overall.model = overall.model,
                ts.model = ts.model))
    
    })

  
  # plot models ####
  
  pdf(date.wrap(paste0("./plots/beta models (", name, ")"), ".pdf"),
      height=4.25, width=3.5)
  par(mfrow=c(3,2), mar=c(0,0,0,0), oma=c(2.5,4,1,1), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=0)
  
  pred.list <- lapply(1:length(model.list), function(n){
    sub.data <- model.list[[n]]$df
    sub.omodel <- model.list[[n]]$overall.model
    sub.tmodel <- model.list[[n]]$ts.model
    
    if(n %% 2 == 0){
      xlims <- range(sub.data$time.var)
    } else { 
      xlims <- rev(range(sub.data$time.var))
      }
  
    plot(x = NULL, y = NULL, 
         xlim = xlims,
         ylim = c(0, 0.99), axes = FALSE,
         xlab = "", ylab = "", yaxs="i")
    
    if(n == 5){
      axis(side=1, at=seq(6000,9000,500), mgp=c(3,0,0))
      mtext(side=1, line=1, text="Years BP", cex=0.8)
      } else {
      axis(side=1, at=seq(6000,9000,500), labels=NA)}
    if(n == 6){
      axis(side=1, at=seq(0,2000,500), mgp=c(3,0,0))
      mtext(side=1, line=1, text = "Years from first time-point", cex=0.8)
      } else {
      axis(side=1, at=seq(0,2000,500), labels=NA)}
    axis(side=1, at=seq(0,10000,100), labels=NA, tcl=-0.125)
    
    if(n %in% c(1,3,5)){
      axis(side=2, at=seq(0,1,0.2), las=1)
      mtext(side=2, line=1.75, 
            text=paste0("Dissimilarity from\n", 
                        c("first", "last", "previous")[(n+1)/2], " time-point"),
            cex=0.8)
      
      }  else {
    axis(side=2, at=seq(0,1,0.2), labels=NA)}
    axis(side=2, at=seq(0,1,0.1), labels=NA, tcl=-0.125)
    
    pred.df <- data.frame(time.var = seq(min(sub.data$time.var),
                                         max(sub.data$time.var),
                                         len=200),
                          ts = sub.data$ts[1],
                          dummy=0)
    
    tpreds <- do.call("rbind", lapply(unique(data$ts), function(ts){
      print(ts)
      
      ts.df <- sub.data[sub.data$ts == ts,]
      if(nrow(ts.df) < 2){return(NULL)}
      
      temp.col <- ts.col[unique(data$ts) == ts]
      
      ts.pred <- pred.df[pred.df$time.var >= min(ts.df$time.var) &
                           pred.df$time.var <= max(ts.df$time.var), ]
      if(nrow(ts.pred) < 2){return(NULL)}
      
      ts.pred$ts = ts

            ts.pred <- cbind(ts.pred,
                       as.data.frame(predict(sub.tmodel, newdata=ts.pred, se.fit=TRUE)))
      
      ts.pred$raw.fit <- plogis(ts.pred$fit)
      ts.pred$upper <- plogis(ts.pred$fit + 1.96*ts.pred$se.fit)
      ts.pred$lower <- plogis(ts.pred$fit - 1.96*ts.pred$se.fit)
      
      lines(x = ts.pred$time.var,
            y = ts.pred$raw.fit, lwd=0.75, col=temp.col)  
      
      return(ts.pred)
      
    }))
      
    opreds <- cbind(pred.df,
                    as.data.frame(predict(sub.omodel, newdata=pred.df, se.fit=TRUE)))
    
    opreds$raw.fit <- plogis(opreds$fit)
    opreds$upper <- plogis(opreds$fit + 1.96*opreds$se.fit)
    opreds$lower <- plogis(opreds$fit - 1.96*opreds$se.fit)
    
  polygon(x=c(opreds$time.var, rev(opreds$time.var)),
          y=c(opreds$upper, rev(opreds$lower)),
          border=NA, col=rgb(0.75,0.75,0.75,0.5))
  
  lines(x = opreds$time.var,
        y = opreds$raw.fit, lwd=1.5)  

  box()
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.935, "y"),
       labels=paste0("(", letters[n],")"),
       font=2, adj=0)
  
  return(list(overall.preds = opreds,
              ts.preds = tpreds))
  
  })  
  
  dev.off()
  
  
  
  return(list(data = data,
              model.list = model.list,
              pred.list = pred.list))
}