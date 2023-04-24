betaExamplePlot <- function(path){
  
  circum.point <- function(cx, cy, r, angle){
    data.frame(x = cx + r * cos(angle * (pi/180)),
               y = cy + r * sin(angle * (pi/180)))
  }
  
  pdf(path, height=6, width=4, useDingbats = FALSE)
  
  split.screen(rbind(c(0.175,0.575,0.74,0.98),
                     c(0.575,0.975,0.74,0.98),
                     c(0.175,0.575,0.5,0.74),
                     c(0.575,0.975,0.5,0.74),
                     c(0.175,0.975,0.075,0.435)))
  
  rad.mod <- list(45,
                  45+cumsum(rep(45, 6)),
                  c(45+cumsum(rep(60, 4)), 310, 320),
                  rep(c(-40,40), 3) + c(5, -5, -15, 15, -5, 5))
  
  core.rads = seq(0.05,0.45,len=6)
  cols <- c("black", "grey60", "orange", "purple")
  
  ord.points <- lapply(1:4, function(n){
    
    screen(n)
    par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0,0), las=1)
    
    rads <- rad.mod[[n]]
    temp.col <- cols[n]
    
    plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", 
         asp=1, xlab="", ylab="", axes=FALSE)
    
    if(n %in% c(1,3)){axis(side=2, mgp=c(3,0.5,0), at=seq(0,0.8,0.2))
      mtext(side=2, line=1.5, text="nMDS2", las=0)
    } else {axis(side=2,labels=NA)}
    if(n %in% 3:4){axis(side=1, at=seq(0,0.8,0.2))
      mtext(side=1, line=0.75, text="nMDS1")} else {
        axis(side=1,labels=NA)}
    
    a <- lapply(core.rads, function(x){draw.circle(x=0.5,y=0.5, r=x, border="grey90", lwd=0.5)})
    points(x=0.5, y=0.5, cex=1.25, pch=16, col=temp.col)
    
    # baseline lines
    straight.point <- rbind(cbind(x=0.5, y=0.5), circum.point(cx=0.5,cy=0.5, r=core.rads, angle=rads))
    
    segments(x0 = straight.point$x[1], y0 = straight.point$y[1],
             x1 = straight.point$x[-1], y1 = straight.point$y[-1],
             col=temp.col, lty="31")
    
    segments(x0=straight.point$x[-length(straight.point$x)],
             x1=straight.point$x[-1],
             y0=straight.point$y[-length(straight.point$x)],
             y1=straight.point$y[-1], col=temp.col)
    points(x=straight.point$x, y=straight.point$y, 
           pch=c(16, rep(21, nrow(straight.point)-1)), bg="white", col=temp.col)
    text(x=straight.point$x, y=straight.point$y,
         labels=1:nrow(straight.point), 
         col=c("white", rep(temp.col, nrow(straight.point)-1)), cex=0.6)
    
    text(x=relative.axis.point(0.025, "x"),
         y=relative.axis.point(ifelse(n <3, 0.935, 0.915), "y"),
         labels=paste0("(", LETTERS[n], ")"), font=2, adj=0)
    box()
    close.screen(n)
    return(straight.point)
    
  })
  
  screen(5)
  par(mar=c(0,0,0,0), ps=8, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  plot(x=NULL, y=NULL, xlim=c(0,0.5), ylim=c(0,0.65), xlab="", ylab="", xaxs="i",
       yaxs="i", axes=FALSE)
  
  axis(side=1,mgp=c(3,0,0))
  mtext(side=1, line=0.85, text=expression("Directional "*beta*" diversity"))
  mtext(side=1, line=1.3, text="(Dissimilarity from time series baseline)")
  axis(side=2)
  mtext(side=2, line=2.15, text=expression("Sequential "*beta*" diversity"), las=0)
  mtext(side=2, line=1.65, text="(Dissimilarity from previous time point)", las=0)
  
  agg.dists <- lapply(1:length(ord.points), function(n){
    temp.points <- ord.points[[n]]
    temp.dist <- as.matrix(dist(temp.points))
    
    dist.base <- temp.dist[-1,1]
    dist.along <- diag(temp.dist[-1,-ncol(temp.dist)])
    return(cbind(dist.base, dist.along))
  })
  
  sapply(1:nrow(agg.dists[[1]]), function(n){
    
    point.dists <- t(sapply(agg.dists, function(x){x[n,]}))
    segments(x0=min(point.dists[,1]), x1 = max(point.dists[,1]),
             y0=min(point.dists[,2]), y1=max(point.dists[,2]),
             lty="31", col="grey80")
  })
  
  sapply(1:length(ord.points), function(n){
    
    dist.base <- agg.dists[[n]][,1]
    dist.along <- agg.dists[[n]][,2]
    
    segments(x0=dist.base[-length(dist.base)],
             x1=dist.base[-1],
             y0=dist.along[-length(dist.along)],
             y1=dist.along[-1], col=cols[n])
    
    points(x=dist.base, y=dist.along, pch=21, bg="white", cex=1.5, col=cols[n])
    text(x=dist.base, y=dist.along, labels=2:(length(dist.along)+1), cex=0.75, col=cols[n])  
    
    text(x=rev(dist.base)[1], y=rev(dist.along)[1],
         labels=paste0("(",LETTERS[n],")"), font=2, col=cols[n], pos=4)
  })
  
  text(x=relative.axis.point(0.01, "x"),
       y=relative.axis.point(0.95,"y"),
       labels="(E)", font=2, adj=0)
  box()
  close.screen(5)
  
  close.screen(all.screens=TRUE)
  dev.off()
  
}