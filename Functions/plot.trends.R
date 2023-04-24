plot.trends <- function(file.name,
                        overall.model,
                        locality.model,
                        axis.label,
                        overall.ylims,
                        locality.ylims,
                        int.ylims,
                        slope.ylims,
                        transform){

# decompose differences in slope and predicted value
  
if(transform=="logit"){
  
  locality.model[[2]][,c("fit", "upper", "lower")] = 
    sapply(locality.model[[2]][,c("fit", "upper", "lower")], plogis)
    
}
  
pred.var <- do.call("rbind", lapply(overall.model[[2]]$raw.date, function(date){
    
    sub.preds <- locality.model[[2]][locality.model[[2]]$raw.date == date,]
    
    if(dim(sub.preds)[1] ==0){return(0)}
    
    int.ss <- (sub.preds$fit - overall.model[[2]]$fit[overall.model[[2]]$raw.date == date])^2
    slope.ss <- (sub.preds$slope - overall.model[[2]]$slope[overall.model[[2]]$raw.date == date])^2
    
    return(data.frame(int.ss = mean(int.ss),
                      slope.ss = mean(slope.ss)))
    
  }))
  
pdf(date.wrap(file.name, ".pdf"), height=4.5, width = 4.25)

split.screen(rbind(c(0.125,0.95,0.8,0.95),
                   c(0.125,0.95,0.25,0.8),
                   c(0.125,0.95,0.1,0.25)))

xlimits <- summary(overall.model[[2]]$raw.date)[c(6,1)] * c(1, 0.975)

screen(1)
par(mar=c(0,0,0,0), ps=8, tcl = -0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, 
     ylim=int.ylims,
     xlim=xlimits,
     axes=FALSE, xlab="", ylab="", yaxs="i")

axis(side=1, labels=NA)
axis(side=2, at=pretty(seq(int.ylims[1], int.ylims[2], length.out=10), 2))
axis(side=2, at=pretty(seq(int.ylims[1], int.ylims[2], length.out=10), 5), labels=NA, tcl = -0.125)
mtext(side=2, text="SD", line=1.25, las=0)

polygon(x=c(par("usr")[1],
            overall.model[[2]]$raw.date,
            par("usr")[2], par("usr")[2]),
        y=c(0, pred.var$int.ss, 0, 0),
        col="goldenrod")
box()

close.screen(1)

screen(2)
par(mar=c(0,0,0,0), ps=8, tcl = -0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, 
     xlim=xlimits,
     ylim=overall.ylims,
     axes=FALSE, xlab="", ylab="")

axis(side=1, labels=NA)
axis(side=2, at=pretty(locality.model[[2]]$fit, 4))
axis(side=2, at=pretty(locality.model[[2]]$fit, 8), labels=NA, tcl = -0.125)
mtext(side=2, text=axis.label, 
      line=1.5, las=0)

with(overall.model[[2]],
     polygon(x=c(raw.date, rev(raw.date)),
             y=c(upper, rev(lower)), 
             border=NA, col="grey80"))

lines(overall.model[[2]]$fit ~ overall.model[[2]]$raw.date, 
      type="l", lwd=2, col="grey50")

sapply(locs[locs.order], function(x){
  
  preds <- locality.model[[2]][locality.model[[2]]$locality == x, ]
  temp.site <- as.numeric(huon_site$site[as.character(huon_site$locality) == x])[1]
  temp.col <- c("blue","red","darkgreen")[temp.site]
  
  lines(preds$fit ~ preds$raw.date, type="l", col=temp.col)
  
})

sapply(locs[locs.order], function(x){
  
  preds <- locality.model[[2]][locality.model[[2]]$locality == x, ]
  temp.site <- as.numeric(huon_site$site[as.character(huon_site$locality) == x])[1]
  temp.col <- c("blue","red","darkgreen")[temp.site]
  
  text(y=preds$fit[1], x=preds$raw.date[1], pos=4,
       labels = (1:length(locs))[locs[locs.order] == x], 
       offset=0.2, font=2, col=temp.col)
  
})

box()
close.screen(2)

screen(3)
par(mar=c(0,0,0,0), ps=8, tcl = -0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, 
     ylim=slope.ylims,
     xlim=xlimits, axes=FALSE, xlab="", ylab="", yaxs="i")

axis(side=1, mgp=c(3,0.1,0))
mtext(side=1, line=1.15, text="Years before present")

axis(side=2, at=pretty(pred.var$slope.ss, 2))
axis(side=2, at=pretty(pred.var$slope.ss, 5), labels=NA, tcl = -0.125)
mtext(side=2, text="Slope SD", line=1.25, las=0)

polygon(x=c(par("usr")[1],
            overall.model[[2]]$raw.date,
            par("usr")[2], par("usr")[2]),
        y=c(0, pred.var$slope.ss, 0, 0),
        col="darkolivegreen1")
box()
close.screen(3)
close.screen(all.screens=TRUE)

# Locality-specific plots
par(mar=c(0,0,0,0), oma=c(3,3,4,2), ps=8, tcl = -0.25, mgp=c(3,0.5,0), las=1, mfrow=c(3,3))

ylimits<-summary(locality.model[[3]]$response)[c(1,6)] * c(1, 1.1)
xlimits<-summary(locality.model[[3]]$pred.date)[c(6,1)]

sapply(1:length(locs), function(n){
  
  x <- locs[locs.order][n]
  
  temp.data <- locality.model[[3]][locality.model[[3]]$locality == x, ]
  temp.pred <- locality.model[[2]][locality.model[[2]]$locality == x, ]
  
  plot(x=NULL, y=NULL, xlim=xlimits, ylim=locality.ylims, xlab="", ylab="",
       axes=FALSE)
  
  temp.site <- as.numeric(huon_site$site[as.character(huon_site$locality) == x])[1]
  temp.site.name <- c("Hubegong", "Kanomi", "Sialum")[temp.site]
  temp.col <- c("blue","red","darkgreen")[temp.site]
  temp.rgb <- col2rgb(temp.col)/255
  
  rect(xleft=par("usr")[1], 
       xright=par("usr")[2], 
       ybottom=par("usr")[3], 
       ytop=par("usr")[4],
       border=NA, col=rgb(temp.rgb[1], temp.rgb[2], temp.rgb[3], 0.2))
  
  if(n %in% c(3,6,9)){
  mtext(side=4, line=0.25, las=0, text=paste0("Site ", temp.site, ": ", temp.site.name),
        col=temp.col)
  }
  
  if(n %in% c(7:9)){
    axis(side=1, mgp=c(3,0.2,0))
  } else { axis(side=1, labels=NA)}
  
  if(n %in% c(1,4,7)){
    axis(side=2, at=pretty(locality.model[[3]]$response, 4))
  } else { axis(side=2,  at=pretty(locality.model[[3]]$response, 4), labels=NA)}

  if(n == 4){mtext(side=2, line=1.5, text=axis.label, las=0, cex=0.85)}  
  if(n == 8){mtext(side=1, line=1.25, text="Years before present", cex=0.85)}  
  
  polygon(x=c(temp.pred$raw.date, rev(temp.pred$raw.date)),
          y=c(temp.pred$upper, rev(temp.pred$lower)),
          border=NA, col="grey70")
  
  points(temp.data$response ~ temp.data$pred.date, pch=16)
  lines(y=temp.pred$fit, x=temp.pred$raw.date)
  text(x=relative.axis.point(0.035, "x"),
       y=relative.axis.point(0.91, "y"),
       labels=paste0(n, ": ", x), font=2, adj=0, cex=1.15)  

  box()
})

dev.off()

return(pred.var)

}