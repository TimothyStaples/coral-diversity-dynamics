diversityTrendsPlotsF <- function(loc.color,
                                 fullPlotPath,
                                 model.list,
                                 rich.df,
                                 turnover.df,
                                 xlims, ylims,
                                 xlabs, ylabs,
                                 
                                 vector.list,
                                 v.xlabs, v.ylabs,
                                 v.xlims, v.ylims,
                                 v.maglim, v.denslim,
                                 v.vec.axis){
  
  library(shape)
  
  preds <- model.list[[1]]$preds[[1]]
  
  # funtion to plot trend arrow in diversity space as well as segments for age
  main.trend <- function(preds, col){
    
    fit1 <- preds$fit.1
    fit2 <- preds$fit.2
    pred.date <- preds$pred.date
    
    lines(fit2 ~ fit1, col=col)  
    Arrows(x0=fit1[3], x1=fit1[1],
           y0=fit2[3], y1=fit2[1],
           arr.type="triangle", arr.length=0.1,
           arr.width=0.1, col=col)
    
    sapply(max(pred.date) - seq(0,5000,250), function(seg.n){
      
      temp.match <- which(pred.date == seg.n)
      base.age <- max(pred.date)
      diff.age <- base.age - seg.n
      
      if(length(temp.match)>0){
        
        target.fit1 <- fit1[temp.match - c(1,0)]
        target.fit2 <- fit2[temp.match - c(1,0)]
        
        if(length(target.fit1)<2 | length(target.fit2) <2){return(NULL)}
        arrows(x0=target.fit1[1], x1=target.fit1[2],
               y0=target.fit2[1], y1=target.fit2[2],
               angle=145, col=col,
               length=ifelse(diff.age %% 1000 == 0, 0.025,0.025),
               lwd=ifelse(diff.age %% 1000 == 0, 0.5,0.5))
        
      }
      
    })
    
  }
  
  pdf(fullPlotPath, height=7*1.45, width=9*1.65, useDingbats = FALSE)

  
  split.screen(rbind(
    
    # div trends
    c(0.06,0.25,0.7,0.97),
    c(0.25,0.44,0.7,0.97),
    c(0.06,0.25,0.375,0.645),
    c(0.25,0.44,0.375,0.645),
    c(0.06,0.25,0.05,0.32),
    c(0.25,0.44,0.05,0.32),
    
    # vector trends
    
    c(0.8,0.99,0.83,0.95),
    c(0.8,0.99,0.725,0.83),
    c(0.55,0.74,0.725,0.95),
    
    c(0.8,0.99,0.53,0.65),
    c(0.8,0.99,0.4,0.53), 
    c(0.55,0.74,0.4,0.65),
    
    c(0.8,0.99,0.2,0.325),
    c(0.8,0.99,0.075,0.2),
    c(0.55,0.74,0.075,0.325)
  ))
  
  # DIV TRENDS ####
  a <- sapply(1:length(model.list), function(n1){
    
    temp.preds <- model.list[[n1]]$preds[[1]]
    x.lab <- xlabs[[n1]]
    y.lab <- ylabs[[n1]]
    
    x.lims<-xlims[[n1]]
    y.lims<-ylims[[n1]]
    
    # Region-level trends ####
    screen(1 + (2*(n1-1)))
    par(mar=c(0,0,0,0), ps=12, las=1, tcl=-0.25, mgp=c(3,0.5,0))
    
    
    print(n1)
    plot(x=NULL, y=NULL, xlim=x.lims, ylim=y.lims, asp=1,
         xlab="",ylab="", yaxs="i", xaxs="i", axes=FALSE)
    axis(side=1, mgp=c(3,0.2,0))
    axis(side=2, las=1)
    mtext(side=2, line=1.5, text=y.lab, las=0)
    mtext(side=1, line=1.5, text=x.lab, at=par("usr")[2])
    
    if(n1==1){mtext(side=3, line=0.2, text="Regional trends", font=2)}
    abline(a=0, b=1, lty="31", col="grey80")
    
    main.trend(temp.preds, "black")
    text(x=relative.axis.point(0.03, "x"),
         y=relative.axis.point(0.935,"y"),
         labels=paste0("(",LETTERS[c(1,3,5)[n1]],")"), adj=0, font=2)
    
    box()
    close.screen(1 + (2*(n1-1)))
    # Locality-level trends ####
    
    screen(2 + (2*(n1-1)))
    par(mar=c(0,0,0,0), ps=12, las=1, tcl=-0.25, mgp=c(3,0.5,0))
    
    plot(x=NULL, y=NULL, xlim=x.lims, ylim=y.lims, asp=1,
         xlab="",ylab="", yaxs="i", xaxs="i", axes=FALSE)
    axis(side=1, mgp=c(3,0.2,0))
    axis(side=2, las=1, labels=NA)
    abline(a=0, b=1, lty="31", col="grey80")
    
    temp.preds <- model.list[[n1]]$preds[[3]]
    
    a <- sapply(levels(loc.color$locality), function(x){
      
      print(x)
      site.range <- range(rich.df$pred.date[rich.df$locality==x])
      
      localBin <- temp.preds$locality == x &
                  temp.preds$pred.date >= site.range[1] &
                  temp.preds$pred.date <= site.range[2]
      
      temp.pred <- temp.preds[localBin,]
    
      site.col = loc.color$col[loc.color$locality == x]
      
      main.trend(temp.pred, site.col)
      
      text(x=temp.pred$fit.1[1],
           y=temp.pred$fit.2[1],
           pos=4, offset=0.25,
           labels=locality.number$number[locality.number$locality == x],
           col=site.col)
      
    })
    
    if(n1==1){mtext(side=3, line=0.2, text="Site trends", font=2)}
    
    text(x=relative.axis.point(0.03, "x"),
         y=relative.axis.point(0.935,"y"),
         labels=paste0("(",LETTERS[c(2,4,6)[n1]],")"), adj=0, font=2)
    
    box()
    close.screen(2 + (2*(n1-1)))
  })
  
  # VECTOR TRENDS ####
  
  colMat <- col2rgb(c("#D4C36A","#116416","#D46D6A",
                      "#413075","#D4C36A")) /255
  colRGB <- rgb(colMat[1,], colMat[2,], colMat[3,], 1)
  
  col.grad <- colorRampPalette(rev(colRGB), alpha=TRUE)(360)
  
  a <- sapply(1:length(vector.list), function(n){
    
    rich.vect <- vector.list[[n]]$data
    angle.gam <- vector.list[[n]]$model
    pred.df <- vector.list[[n]]$model.pred
    angle.dens <- vector.list[[n]]$dens
    radii <- list(c(2.9,2.4,2.7,2.4,2.9),
                  c(0.9,0.775,0.85,0.775,0.9),
                  c(0.875,0.725,0.825,0.725,0.875))[[n]]
    
    vec.axis <- v.vec.axis[[n]]
    
    screen.mod <- 6 + (n-1)*3
    
    screen(screen.mod + 1)
    par(mar=c(0,0,0,0), ps=12, las=1, tcl=-0.25, mgp=c(3,0.5,0))
    
    plot(x=NULL,y=NULL, xlim=c(0,359),ylim=v.denslim[[n]],
         ylab="", xlab="", xaxs="i", yaxs="i", axes=FALSE)
    axis(side=2, at=seq(0,0.4,0.1))
    mtext(side=2, text="Density", las=0, line=1.7)
    image(x=seq(0,359,1), y=c(par("usr")[3],par("usr")[4]),
          z=matrix(seq(0,358,1), ncol=1), add=TRUE,
          col=col.grad, useRaster=TRUE)
    polygon(x=c(par("usr")[1], angle.dens$x, par("usr")[2]),
            y=c(par("usr")[4], angle.dens$y, par("usr")[4]),
            col="white")

    text(x=relative.axis.point(0.03, "x"),
         y=relative.axis.point(0.86, "y"),
         labels=paste0("(",LETTERS[6+c(2,5,8)][n],")"), font=2, adj=0)
    box()
    close.screen(screen.mod+1)
    
    screen(screen.mod+2)
    par(mar=c(0,0,0,0), ps=12, las=1, tcl=-0.25, mgp=c(3,0.5,0))
    
    plot(x=NULL,y=NULL, xlim=c(0,359),ylim=v.maglim[[n]],
         ylab="", xlab="", xaxs="i", yaxs="i", axes=FALSE)
    axis(side=2, at=pretty(seq(par("usr")[3], par("usr")[4], len=50), 4))
    axis(side=1, mgp=c(3,0.2,0), at=seq(0,360,60))
    axis(side=1, mgp=c(3,0.2,0), at=seq(0,360,30), tcl=-0.125, labels=NA)
    
    mtext(side=1, line=1.2, text=expression("Angle ("*degree*")"))
    mtext(side=2, text="Vector length", las=0, line=1.7)
    
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
    
    text(x=relative.axis.point(0.03, "x"),
         y=relative.axis.point(0.86, "y"),
         labels=paste0("(",LETTERS[6+c(3,6,9)][n],")"), font=2, adj=0)
    axis(side=3, labels=NA, at=seq(0,360,60), tcl=0.25)
    axis(side=3, mgp=c(3,0.2,0), at=seq(0,360,30), tcl=0.125, labels=NA)
    box()
    close.screen(screen.mod+2)
    
    screen(screen.mod+3)
    par(mar=c(0,0,0,0), ps=12, las=1, tcl=-0.25, mgp=c(3,0.5,0))
    
    plot(x=NULL, y=NULL, xlim=v.xlims[[n]], 
         ylim=v.ylims[[n]], asp=1,
         axes=FALSE, xlab="", ylab="")
    
    axis(side=1, at=vec.axis, mgp=c(3,0.2,0))
    axis(side=2, at=vec.axis)
    
    mtext(side=1, line=1.2, text=v.xlabs[[n]])
    mtext(side=2, line=2,
          text=v.ylabs[[n]], las=0)
    
    # dividers
    segments(x0=0, y0=0,
             x1=c(relative.axis.point(0.8,"x"), 0, relative.axis.point(0.2, "x"),0),
             y1=c(0, relative.axis.point(0.85, "y"), 0, relative.axis.point(0.15, "y")),
             col="grey60")
    text(x=c(relative.axis.point(0.8,"x"), 0, relative.axis.point(0.2, "x"),0),
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
    
    text(x=relative.axis.point(0.03, "x"),
         y=relative.axis.point(0.935, "y"),
         labels=paste0("(",LETTERS[6+c(1,4,7)][n],")"), font=2, adj=0)
    # mtext(side=3, line=0, 
    #       text=c(expression(bold(alpha*" Diversity")),
    #              expression(bold("Directional "*beta*" diversity")),
    #              expression(bold("Sequential "*beta*" diversity")))[n])
    box()
    close.screen(screen.mod+3)
    
  })
  close.screen(all.screens=TRUE)
  dev.off() 
  
  
  close.screen(all.screens=TRUE)
  dev.off()
  
}