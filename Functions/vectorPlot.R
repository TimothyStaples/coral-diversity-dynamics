vectorPlot <- function(path,
                       vector.list, xlabs, ylabs, 
                       xlims, ylims, maglim, denslim,
                       vec.axis){
  
  pdf(path,
      height=7.25, width=5.8, useDingbats = FALSE)
  
  colMat <- col2rgb(c("#D4C36A","#116416","#D46D6A",
                      "#413075","#D4C36A")) /255
  colRGB <- rgb(colMat[1,], colMat[2,], colMat[3,], 1)
  
  col.grad <- colorRampPalette(rev(colRGB), alpha=TRUE)(360)
  
  split.screen(rbind(c(0.53,0.96,0.86,0.99), # alpha density
                     c(0.53,0.96,0.73,0.86),
                     c(0.125,0.45,0.73,0.99),
                     
                     c(0.53,0.96,0.53,0.66),
                     c(0.53,0.96,0.40,0.53), 
                     c(0.125,0.45,0.40,0.66),
                     
                     c(0.53,0.96,0.2,0.33),
                     c(0.53,0.96,0.07,0.2),
                     c(0.125,0.45,0.07,0.33)))
  
  sapply(1:length(vector.list), function(n){
    
    rich.vect <- vector.list[[n]]$data
    angle.gam <- vector.list[[n]]$model
    pred.df <- vector.list[[n]]$model.pred
    angle.dens <- vector.list[[n]]$dens
    radii <- list(c(2.9,2.4,2.7,2.4,2.9),
                  c(0.9,0.775,0.85,0.775,0.9),
                  c(0.875,0.725,0.825,0.725,0.875))[[n]]

    vec.axis <- vec.axis[[n]]
    
    screen.mod <- (n-1)*3
    screen(screen.mod + 1)
    par(mar=c(0,0,0,0), oma=c(0,0,0,0), 
        tcl=-0.25, ps=8, mgp=c(3,0.5,0), las=1)
    plot(x=NULL,y=NULL, xlim=c(0,359),ylim=denslim[[n]],
         ylab="", xlab="", xaxs="i", yaxs="i", axes=FALSE)
    axis(side=2, at=seq(0,0.4,0.1))
    mtext(side=2, text="Density", las=0, line=1.5)
    image(x=seq(0,359,1), y=c(par("usr")[3],par("usr")[4]),
          z=matrix(seq(0,358,1), ncol=1), add=TRUE,
          col=col.grad, useRaster=TRUE)
    polygon(x=c(par("usr")[1], angle.dens$x, par("usr")[2]),
            y=c(par("usr")[4], angle.dens$y, par("usr")[4]),
            col="white")
    # dividers <- seq(45,360,90)
    # centers <- c(22.5, 90, 180, 270, 360-22.5)
    # abline(v=seq(45,360,90), col=rgb(0.5,0.5,0.5,1), lwd=2, lty="31")
    # text(x = centers, y=angle.dens$y[match(round(centers,0), round(angle.dens$x,0))],
    #      pos=3,
    #      labels=paste0("(", c("1a",2,3,4,"1b"), ")"),
    #      col=col.grad[round(centers)], font=2)
    text(x=relative.axis.point(0.01, "x"),
         y=relative.axis.point(0.90, "y"),
         labels=paste0("(",LETTERS[c(2,5,8)][n],")"), font=2, adj=0)
    box()
    close.screen(screen.mod+1)
    
    screen(screen.mod+2)
    par(mar=c(0,0,0,0), oma=c(0,0,0,0), 
        tcl=-0.25, ps=8, mgp=c(3,0.5,0), las=1)
    plot(x=NULL,y=NULL, xlim=c(0,359),ylim=maglim[[n]],
         ylab="", xlab="", xaxs="i", yaxs="i", axes=FALSE)
    axis(side=2, at=pretty(seq(par("usr")[3], par("usr")[4], len=50), 4))
    axis(side=1, mgp=c(3,0,0), at=seq(0,360,60))
    axis(side=1, mgp=c(3,0,0), at=seq(0,360,30), tcl=-0.125, labels=NA)
    
    mtext(side=1, line=0.95, text=expression("Angle ("*degree*")"))
    mtext(side=2, text="Vector length", las=0, line=1.5)
    
    # make colored confidence intervals
    temp.col <- col2rgb(col.grad)/255
    image(x=seq(0,359,1), y=c(par("usr")[3],par("usr")[4]),
          z=matrix(seq(0,358,1), ncol=1), add=TRUE,
          col=rgb(temp.col[1,], temp.col[2,], temp.col[3,], 0.3), useRaster=TRUE)
    polygon(x=c(par("usr")[1], par("usr")[1], pred.df$div.angle, par("usr")[2], par("usr")[2]),
            y=c(par("usr")[4], pred.df$upper[1], pred.df$upper, rev(pred.df$upper)[1], par("usr")[4]),
            col="white", border=NA)
    polygon(x=c(par("usr")[1], par("usr")[1], pred.df$div.angle, par("usr")[2], par("usr")[2]),
            y=c(par("usr")[3], pred.df$lower[1], pred.df$lower, rev(pred.df$lower)[1], par("usr")[3]),
            col="white", border=NA)
    
    a<-sapply(2:nrow(pred.df), function(n){
      lines(x=pred.df$div.angle[(n-1):n], y=pred.df$fit[(n-1):n], 
            lwd=2, col=col.grad[n])
    })
    
    text(x=relative.axis.point(0.01, "x"),
         y=relative.axis.point(0.9, "y"),
         labels=paste0("(",LETTERS[c(3,6,9)][n],")"), font=2, adj=0)
    axis(side=3, labels=NA, at=seq(0,360,60), tcl=0.25)
    axis(side=3, mgp=c(3,0,0), at=seq(0,360,30), tcl=0.125, labels=NA)
    box()
    close.screen(screen.mod+2)
    
    screen(screen.mod+3)
    par(mar=c(0,0,0,0), oma=c(0,0,0,0), 
        tcl=-0.25, ps=8, mgp=c(3,0.5,0), las=1)
    plot(x=NULL, y=NULL, xlim=xlims[[n]], 
         ylim=ylims[[n]], asp=1,
         axes=FALSE, xlab="", ylab="")
    
    axis(side=1, at=vec.axis, mgp=c(3,0,0))
    axis(side=2, at=vec.axis)
    
    mtext(side=1, line=0.95, text=xlabs[[n]])
    mtext(side=2, line=1.5,
          text=ylabs[[n]], las=0)
    
    # dividers
    segments(x0=0, y0=0,
             x1=c(relative.axis.point(0.85,"x"), 0, relative.axis.point(0.15, "x"),0),
             y1=c(0, relative.axis.point(0.85, "y"), 0, relative.axis.point(0.15, "y")),
             col="grey60")
    text(x=c(relative.axis.point(0.85,"x"), 0, relative.axis.point(0.15, "x"),0),
         y=c(0, relative.axis.point(0.85, "y"), 0, relative.axis.point(0.15, "y")),
         labels=c(expression("0"*degree),
                  expression("90"*degree),
                  expression("180"*degree),
                  expression("270"*degree)),
         pos=c(4,3,2,1), offset=0.25, col="grey60")
    
    rich.vect$angle.round <- round(rich.vect$div.angle, 0)
    rich.vect.nona <- rich.vect[!is.na(rich.vect$primary),]
    Arrows(x0 = 0, x1 = rich.vect.nona$primary, 
           y0 = 0, y1 = rich.vect.nona$secondary,
           arr.type="triangle", arr.width=0.1, arr.length=0.1,
           col = col.grad[match(rich.vect.nona$angle.round, seq(0,369,1))])
    
    text(x=relative.axis.point(0.01, "x"),
         y=relative.axis.point(0.96, "y"),
         labels=paste0("(",LETTERS[c(1,4,7)][n],")"), font=2, adj=0)
    # mtext(side=3, line=0, 
    #       text=c(expression(bold(alpha*" Diversity")),
    #              expression(bold("Directional "*beta*" diversity")),
    #              expression(bold("Sequential "*beta*" diversity")))[n])
    box()
    close.screen(screen.mod+3)
    
  })
  close.screen(all.screens=TRUE)
  dev.off() 
}