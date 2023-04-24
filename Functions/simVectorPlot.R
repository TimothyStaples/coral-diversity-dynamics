simVectorPlot <- function(path,
                       vector.list){
  
  pdf(path,
      height=8.4, width=6.75, useDingbats = FALSE)
  
  colMat <- col2rgb(c("#D4C36A","#116416","#D46D6A",
                      "#413075","#D4C36A")) /255
  colRGB <- rgb(colMat[1,], colMat[2,], colMat[3,], 1)
  
  col.grad <- colorRampPalette(rev(colRGB), alpha=TRUE)(360)
  colramp <- colorRamp(rev(colRGB), alpha=TRUE)
  colsim <- colorRamp(rgb(colMat[1,], colMat[2,], colMat[3,], 0.1), alpha=TRUE)
  
  split.screen(rbind(c(0.53,0.96,0.86,0.99), # alpha density
                     c(0.53,0.96,0.73,0.86),
                     c(0.125,0.45,0.73,0.99),
                     
                     c(0.53,0.96,0.53,0.66),
                     c(0.53,0.96,0.40,0.53), 
                     c(0.125,0.45,0.40,0.66),
                     
                     c(0.53,0.96,0.2,0.33),
                     c(0.53,0.96,0.07,0.2),
                     c(0.125,0.45,0.07,0.33)))
  sapply(1:length(vector.list[[3]]), function(n){
    print(n)
    rich.vect.list <- lapply(vector.list, function(x){x[[n]]$data})
    angle.gams <- lapply(vector.list, function(x){x[[n]]$model})
    pred.list <- lapply(vector.list, function(x){x[[n]]$model.pred})
    angle.dens <- lapply(vector.list, function(x){x[[n]]$dens})
    
    dens.lims <- c(0.37, 0.37, 0.37)[n]
    mag.lims <- c(1.9, 0.38, 0.48)[n]
    vec.lims <- c(2.5, 0.8, 0.75)[n]
    radii <- list(c(2.9,2.4,2.7,2.4,2.9),
                  c(0.9,0.775,0.85,0.775,0.9),
                  c(0.875,0.725,0.825,0.725,0.875))[[n]]
    vec.axis <- list(seq(-3,3,1),
                     seq(-1,1,0.25),
                     seq(-1,1,0.25))[[n]]
    
    screen.mod <- (n-1)*3
    screen(screen.mod + 1)
    par(mar=c(0,0,0,0), oma=c(0,0,0,0), 
        tcl=-0.25, ps=8, mgp=c(3,0.5,0), las=1)
    
    plot(x=NULL,y=NULL, xlim=c(0,359),ylim=c(0,dens.lims),
         ylab="", xlab="", xaxs="i", yaxs="i", axes=FALSE)
    
    axis(side=2, at=seq(0,0.4,0.1))
    mtext(side=2, text="Density", las=0, line=1.5)
    
    a <- sapply(angle.dens, function(x){
      
      segments(x0 = x$x[-nrow(x)],
               x1 = x$x[-1],
               y0 = x$y[-nrow(x)],
               y1 = x$y[-1],
               col=rgb(colramp(x$x/max(x$x))[,1]/255,
                       colramp(x$x/max(x$x))[,2]/255,
                       colramp(x$x/max(x$x))[,3]/255,
                       0.05))
      
    })
    
    meanDens <- rowMeans(sapply(angle.dens, function(x){x$y}))
    x <- angle.dens[[1]]
    segments(x0 = x$x[-nrow(x)],
             x1 = x$x[-1],
             y0 = meanDens[-nrow(x)],
             y1 = meanDens[-1],
             col=rgb(colramp(x$x/max(x$x))[,1]/255,
                     colramp(x$x/max(x$x))[,2]/255,
                     colramp(x$x/max(x$x))[,3]/255,
                     1), lwd=2)
    
    text(x=relative.axis.point(0.01, "x"),
         y=relative.axis.point(0.90, "y"),
         labels=paste0("(",LETTERS[c(2,5,8)][n],")"), font=2, adj=0)
    box()
    close.screen(screen.mod+1)
    
    screen(screen.mod+2)
    par(mar=c(0,0,0,0), oma=c(0,0,0,0), 
        tcl=-0.25, ps=8, mgp=c(3,0.5,0), las=1)
    plot(x=NULL,y=NULL, xlim=c(0,359),ylim=c(0,mag.lims),
         ylab="", xlab="", xaxs="i", yaxs="i", axes=FALSE)
    axis(side=2, at=pretty(seq(par("usr")[3], par("usr")[4], len=50), 4))
    axis(side=1, mgp=c(3,0,0), at=seq(0,360,60))
    axis(side=1, mgp=c(3,0,0), at=seq(0,360,30), tcl=-0.125, labels=NA)
    
    mtext(side=1, line=0.75, text=expression("Angle ("*degree*")"))
    mtext(side=2, text="Vector length", las=0, line=1.5)
    
    sapply(pred.list, function(x){
      
      segments(x0 = x$div.angle[-nrow(x)],
               x1 = x$div.angle[-1],
               y0 = x$fit[-nrow(x)],
               y1 = x$fit[-1],
               col=rgb(colramp(x$div.angle/max(x$div.angle))[,1]/255,
                       colramp(x$div.angle/max(x$div.angle))[,2]/255,
                       colramp(x$div.angle/max(x$div.angle))[,3]/255,
                       0.05))
          })
    meanMag <- rowMeans(sapply(pred.list, function(x){x$fit}))
    x <- pred.list[[1]]
    
    segments(x0 = x$div.angle[-nrow(x)],
             x1 = x$div.angle[-1],
             y0 = meanMag[-nrow(x)],
             y1 = meanMag[-1],
             col=rgb(colramp(x$div.angle/max(x$div.angle))[,1]/255,
                     colramp(x$div.angle/max(x$div.angle))[,2]/255,
                     colramp(x$div.angle/max(x$div.angle))[,3]/255,
                     1), lwd=2)
    
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
    plot(x=NULL, y=NULL, xlim=c(-vec.lims,vec.lims), 
         ylim=c(-vec.lims,vec.lims), asp=1,
         axes=FALSE, xlab="", ylab="")
    
    axis(side=1, at=vec.axis, mgp=c(3,0,0))
    axis(side=2, at=vec.axis)
    
    primary.lab <- list(list(bquote(Delta*" Primary "*alpha),bquote('("complexity")')),
                        list(bquote(Delta*" Primary directional "*beta), bquote('("dissimilarity")')),
                        list(bquote(Delta*" Primary sequential "*beta), bquote('("turnover")')))[[n]]
    
    secondary.lab <- list(list(bquote(Delta*" Secondary "*alpha),bquote('("redundancy")')),
                          list(bquote(Delta*" Secondary directional "*beta), bquote('("conservation")')),
                          list(bquote(Delta*" Secondary sequential "*beta), bquote('("persistence")')))[[n]]
    
    mtext(side=1, line=c(0.75, 1.35), text=do.call(expression, primary.lab))
    mtext(side=2, line=c(ifelse(n==1, 1.6, 2.6),
                         ifelse(n==1, 1, 2)),
          text=do.call(expression, secondary.lab), las=0)
    
    # dividers
    segments(x0=0,
             x1=c(relative.axis.point(0.85,"x"), 0, relative.axis.point(0.15, "x"),0),
             y0=0,
             y1=c(0, relative.axis.point(0.85, "y"), 0, relative.axis.point(0.15, "y")),
             col="grey60")
    text(x=c(relative.axis.point(0.85,"x"), 0, relative.axis.point(0.15, "x"),0),
         y=c(0, relative.axis.point(0.85, "y"), 0, relative.axis.point(0.15, "y")),
         labels=c(expression("0"*degree),
                  expression("90"*degree),
                  expression("180"*degree),
                  expression("270"*degree)),
         pos=c(4,3,2,1), offset=0.25, col="grey60")
    
    rich.vect <- do.call("rbind", lapply(rich.vect.list, function(x){
    
      x$angle.round <- round(x$div.angle, 0)
      x[!is.na(x$primary),]
      
    }))
      
    with(rich.vect[sample(1:nrow(rich.vect), round(nrow(rich.vect)*1)),],
    Arrows(x0 = 0, x1 = primary, 
           y0 = 0, y1 = secondary,
           arr.type="triangle", arr.width=0.1, arr.length=0.1,
           col=rgb(colramp(angle.round/360)[,1]/255,
                   colramp(angle.round/360)[,2]/255,
                   colramp(angle.round/360)[,3]/255,
                   0.25)))
    
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