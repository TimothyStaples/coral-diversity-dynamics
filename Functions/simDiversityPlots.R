simDiversityPlots <- function(loc.color,
                              simModelList,
                              regionPlotPath,
                              sitePlotPath,
                              rich.df){
  
  library(shape)
  
  # funtion to plot trend arrow in diversity space as well as segments for age
  main.trend <- function(preds, col){
    lines(preds$fit.2 ~ preds$fit.1, col=col)  
    Arrows(x0=preds$fit.1[3], x1=preds$fit.1[1],
           y0=preds$fit.2[3], y1=preds$fit.2[1],
           arr.type="triangle", arr.length=0.1,
           arr.width=0.1, col=col)
    
    sapply(max(preds$pred.date) - seq(0,5000,250), function(seg.n){
      
      temp.match <- which(preds$pred.date == seg.n)
      base.age <- max(preds$pred.date)
      diff.age <- base.age - seg.n
      
      if(length(temp.match)>0){
        
        target <- preds[temp.match - c(1,0),]
        if(nrow(target)<2){return(NULL)}
        arrows(x0=target$fit.1[1], x1=target$fit.1[2],
               y0=target$fit.2[1], y1=target$fit.2[2],
               angle=145, col=col,
               length=ifelse(diff.age %% 1000 == 0, 0.025,0.025),
               lwd=ifelse(diff.age %% 1000 == 0, 0.5,0.5))
        
      }
      
    })
    
  }
  
  # Regional plot only ####
  
  pdf(regionPlotPath, height=6.4, width=2.425, useDingbats = FALSE)
  
  par(mfrow=c(3,1), mar=c(3.5,5,0.25,1), oma=c(0,0,0,0), 
      ps=10, tcl=-0.25, mgp=c(3,0.5,0))
  
  # set up each plot
  a<-sapply(1:length(simModelList[[1]]), function(n1){
    
    print(n1)
    
  x.lab <- list(list(bquote("Primary "*alpha),bquote('("complexity")')),
                list(bquote("Primary directional "*beta), bquote('("dissimilarity")')),
                list(bquote("Primary sequential "*beta), bquote('("turnover")')))[[n1]]
  
  y.lab <- list(list(bquote("Secondary "*alpha),bquote('("redundancy")')),
                list(bquote("Secondary directional "*beta), bquote('("conservation")')),
                list(bquote("Secondary sequential "*beta), bquote('("persistence")')))[[n1]]
  
  x.lims<-list(c(4.4,7.8), c(0.35,1.125), c(0.25,0.87))[[n1]]
  y.lims<-list(c(-1,1), c(-0.3,0.3), c(-0.2,0.2))[[n1]]
  
  # Region-level trends ####
  plot(x=NULL, y=NULL, xlim=x.lims, ylim=y.lims, asp=1,
       xlab="",ylab="", yaxs="i", xaxs="i", axes=FALSE)
  axis(side=1, mgp=c(3,0.1,0))
  axis(side=2, las=1)
  mtext(side=2, line=c(3,2), text=do.call(expression, y.lab), cex=0.9, las=0)
  mtext(side=1, line=c(1.5,2.5), text=do.call(expression, x.lab), cex=0.9)
  abline(h=0, lty="31", col="grey80")

  # each sim plot
  simPreds <- lapply(1:length(simModelList), function(n){
  
    model.list <- simModelList[[n]][[3]]
    temp.preds <- model.list[[n1]]$preds[[1]]
    
    main.trend(temp.preds, rgb(0.5,0.5,0.5,0.1))
    
    return(temp.preds)
  })
  
  # now 'average' trends across the sims
  simMeans <- cbind(rowMeans(sapply(simPreds, function(x){x$fit.1})),
                    rowMeans(sapply(simPreds, function(x){x$fit.2})))
  simMeanPred <- simPreds[[1]]          
  simMeanPred$fit.1 = simMeans[,1]
  simMeanPred$fit.2 = simMeans[,2]
  
  main.trend(simMeanPred, "black")
  
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.95,"y"),
       labels=paste0("(",LETTERS[n1],")"), adj=0, font=2)
  
  box()
  
  })
  dev.off()
  
  # Site-level trends ####
  
  pdf(sitePlotPath, height=6.4, width=2.425, useDingbats = FALSE)
  
  par(mfrow=c(3,1), mar=c(3.5,5,0.25,1), oma=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0))
  
  a<-sapply(1:length(simModelList[[1]]), function(n1){
    
    x.lab <- list(list(bquote("Primary "*alpha),bquote('("complexity")')),
                  list(bquote("Primary directional "*beta), bquote('("dissimilarity")')),
                  list(bquote("Primary sequential "*beta), bquote('("turnover")')))[[n1]]
    
    y.lab <- list(list(bquote("Secondary "*alpha),bquote('("redundancy")')),
                  list(bquote("Secondary directional "*beta), bquote('("conservation")')),
                  list(bquote("Secondary sequential "*beta), bquote('("persistence")')))[[n1]]
    
    x.lims<-list(c(4.4,7.8), c(0.35,1.125), c(0.25,0.87))[[n1]]
    y.lims<-list(c(-1,1), c(-0.3,0.3), c(-0.2,0.2))[[n1]]
    
    plot(x=NULL, y=NULL, xlim=x.lims, ylim=y.lims, asp=1,
         xlab="",ylab="", yaxs="i", xaxs="i", axes=FALSE)
    axis(side=1, mgp=c(3,0.1,0))
    axis(side=2, las=1)
    mtext(side=2, line=c(3,2), text=do.call(expression, y.lab), cex=0.9, las=0)
    mtext(side=1, line=c(1.5,2.5), text=do.call(expression, x.lab), cex=0.9)
    abline(h=0, lty="31", col="grey80")
    
    a <- sapply(levels(loc.color$locality), function(x){
    
    simPreds <- lapply(1:length(simModelList), function(n){
      
      model.list <- simModelList[[n]][[3]]
      temp.preds <- model.list[[n1]]$preds[[3]]
      temp.pred <- temp.preds[temp.preds$locality == x, ]
      site.range <- range(rich.df$pred.date[rich.df$locality==x])
      temp.pred <- temp.pred[temp.pred$pred.date >= site.range[1] &
                               temp.pred$pred.dat <= site.range[2],]
      return(temp.pred)
    })
      
    site.col = loc.color$col[loc.color$locality == x]
    simCol <- c(col2rgb(site.col)/255, 0.1)
    simCol <- rgb(simCol[1], simCol[2], simCol[3], 0.025)
      
    sapply(simPreds, function(x){main.trend(x, simCol)})
      
    simMeans <-  cbind(rowMeans(sapply(simPreds, function(x){x$fit.1})),
                        rowMeans(sapply(simPreds, function(x){x$fit.2})))
    simMeanPred <- simPreds[[1]]          
    simMeanPred$fit.1 = simMeans[,1]
    simMeanPred$fit.2 = simMeans[,2]
    
    main.trend(simMeanPred, site.col)
    
    text(x=simMeanPred$fit.1[1],
           y=simMeanPred$fit.2[1],
           pos=4, offset=0.25,
           labels=locality.number$number[locality.number$locality == x],
           col=site.col)
      
    })

    text(x=relative.axis.point(0.02, "x"),
         y=relative.axis.point(0.95,"y"),
         labels=paste0("(",LETTERS[n1],")"), adj=0, font=2)
    box()
    
  })
  

  
  dev.off()
  
  
}