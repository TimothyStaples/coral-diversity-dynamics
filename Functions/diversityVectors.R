diversityVectors <- function(rich.df, turnover.df){
  
  rich.vect <- do.call("rbind", lapply(split(rich.df, rich.df$locality),
                                       function(x){
                                         
                                         # sort x into order
                                         x.sort <- x[order(x$pred.date, decreasing=TRUE),]
                                         
                                         x.sort$d.coordH0 = c(NA, diff(x.sort$coordH0))
                                         x.sort$d.ratioH0 = c(NA, diff(x.sort$ratioH0))
                                         x.sort$d.coordH1 = c(NA, diff(x.sort$coordH1))
                                         x.sort$d.ratioH1 = c(NA, diff(x.sort$ratioH1))
                                         x.sort$d.coordH2 = c(NA, diff(x.sort$coordH2))
                                         x.sort$d.ratioH2 = c(NA, diff(x.sort$ratioH2))
                                         
                                         x.sort$pre.coordH0 <- c(NA, x.sort$coordH0[-nrow(x.sort)])
                                         x.sort$pre.ratioH0 <- c(NA, x.sort$ratioH0[-nrow(x.sort)])
                                         x.sort$pre.coordH1 <- c(NA, x.sort$coordH1[-nrow(x.sort)])
                                         x.sort$pre.ratioH1 <- c(NA, x.sort$ratioH1[-nrow(x.sort)])
                                         x.sort$pre.coordH2 <- c(NA, x.sort$coordH2[-nrow(x.sort)])
                                         x.sort$pre.ratioH2 <- c(NA, x.sort$ratioH2[-nrow(x.sort)])

                                         x.sort$age.diff <- c(NA, abs(diff(x.sort$pred.date)))
                                         
                                         return(x.sort)
                                       }))
  
  turnover.vect <- do.call("rbind", lapply(split(turnover.df, turnover.df$locality),
                                           function(x){
                                             
                                             # sort x into order
                                             x.sort <- x[order(x$pred.date, decreasing=TRUE),]
                                             
                                             # now calculate vectors
                                             x.sort$d.diss <- c(NA, diff(x.sort$dissimilarity))
                                             x.sort$d.cons <- c(NA, diff(x.sort$conservation))
                                             x.sort$pre.diss <- c(NA, x.sort$dissimilarity[-nrow(x.sort)])
                                             x.sort$pre.cons <- c(NA, x.sort$conservation[-nrow(x.sort)])
                                             
                                             x.sort$d.turn <- c(NA, diff(x.sort$turnover))
                                             x.sort$d.pers <- c(NA, diff(x.sort$persistance))
                                             x.sort$pre.turn <- c(NA, x.sort$turnover[-nrow(x.sort)])
                                             x.sort$pre.pers <- c(NA, x.sort$persistance[-nrow(x.sort)])
                                             
                                             x.sort$age.diff <- c(NA, abs(diff(x.sort$pred.date)))
                                             
                                             return(x.sort)
                                           }))
  
  data.list <- list(rich.vect, rich.vect, rich.vect, 
                    turnover.vect, turnover.vect)
  
  vector.list <- lapply(1:length(data.list), function(n){
    
    print(n)
    
    rich.vect <- data.list[[n]]
    primary <- c("d.coordH0", "d.coordH1", "d.coordH2", "d.diss", "d.turn")[n]
    secondary <- c("d.ratioH2", "d.ratioH2", "d.ratioH2", "d.cons", "d.pers")[n]
    
    # turn rich vectors into angles
    rich.vect$div.angle <- cont.angles(rich.vect[,primary], rich.vect[,secondary])
    rich.vect$div.mag <- sqrt(rich.vect[,primary]^2 + rich.vect[,secondary]^2)
    rich.vect$primary <- rich.vect[,primary]
    rich.vect$secondary <- rich.vect[,secondary]
    rich.vect$dummy=1
    
    aC <- circular(rich.vect$div.angle[!is.na(rich.vect$div.angle)], 
                   units="degrees")
    angle.dens <- density(aC, bw=10, kernel="vonmises")
    angle.dens <- data.frame(x = as.numeric(angle.dens$x), 
                             y = as.numeric(angle.dens$y))
    
    rich.vect$log.mag <- log(rich.vect$div.mag+1)
    
    angle.gam <- gam(log.mag ~ s(div.angle, bs="cc", k=-1) + age.diff + s(site, locality, bs="re", by=dummy), 
                     data=rich.vect[!is.na(log(rich.vect$div.mag)),])
    # plot(angle.gam)
    # simRes <- simulateResiduals(angle.gam)
    # plot(simRes)
    
    pred.df <- data.frame(div.angle = seq(0,359,1), 
                          age.diff=mean(rich.vect$age.diff, na.rm=TRUE),
                          site="a", locality="a", dummy=0)
    pred.df <- cbind(pred.df,
                     as.data.frame(predict(angle.gam, newdata=pred.df, se.fit=TRUE)))
    pred.df$upper <- exp(pred.df$fit + 1.96 * pred.df$se.fit)-1
    pred.df$lower <- exp(pred.df$fit - 1.96 * pred.df$se.fit)-1
    pred.df$fit <- exp(pred.df$fit) - 1
    
    return(list(data=rich.vect,
                model=angle.gam,
                model.pred=pred.df,
                dens=angle.dens))
  })
  
  return(vector.list)
  
}