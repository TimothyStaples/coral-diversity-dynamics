divRotationPlot <- function(path){
  
  
  richness.col <- "grey70"
  top.beta.col <- "grey70"
  seq.beta.col <- "grey70"
  sign.locs = c(0.17,0.83)
  lower.pos = 0.2
  upper.pos = 0.8
  lower.inner = 0.415
  upper.inner = 0.585
  lower.outer = 0.25
  upper.outer = 0.8
  arrow.length = 0.425
  fade.rad = 0.19
  
  pdf(path,height=2.7, width=5.5, useDingbats = FALSE)
  
  par(mfcol=c(1,2), mar=c(2,2.5,0.5,1), oma=c(1,0.5,0,0),
      ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
  
  #                   Alpha diagonal ####
  plot(x=NULL, y=NULL, 
       xlim=c(1,6), ylim=c(1,6), xlab="", ylab="", axes=FALSE,
       xaxs="i", yaxs="i")
  axis(side=1, mgp=c(3,0.2,0))
  axis(side=2)
  mtext(side=2, line=1.5, text=expression("Taxonomic diversity"), las=0)
  mtext(side=1, line=1.35, text=expression("Functional diversity"))
  
  # artifiial axes back-transformed onto raw scales
  a1 <- seq(relative.axis.point(0.2, "x"),
            relative.axis.point(0.8, "x"), len=8)
  a0 = a1 - (relative.axis.point(0.01, "x") - par("usr")[1])
  a2 = a1 + (relative.axis.point(0.01, "x") - par("usr")[1])
  
  segments(x0 = a1[1], x1 = rev(a1)[1],
           y0 = a1[1], y1 = rev(a1)[1], col=richness.col)
  
  segments(x0 = a1, x1 = a0,
           y0 = a1, y1 = a2, col=richness.col)
  
  segments(x0 = a1[1], x1 = rev(a1)[1],
           y0 = rev(a1)[1], y1 = a1[1], col=richness.col)
  
  segments(x0 = a1, x1 = a2,
           y0 = rev(a1), y1 = rev(a2), col=richness.col)
  
  text(x=relative.axis.point(sign.locs[1],"x"), y=relative.axis.point(sign.locs[1],"y"), 
       labels="-", font=2, col=richness.col, cex=1.2)
  text(x=relative.axis.point(sign.locs[2],"x"), y=relative.axis.point(sign.locs[2],"y"), 
       labels="+", font=2, col=richness.col, cex=1.2)
  text(x=relative.axis.point(sign.locs[1],"x"), y=relative.axis.point(sign.locs[2],"y"), 
       labels="+", font=2, col=richness.col, cex=1.2)
  text(x=relative.axis.point(sign.locs[2],"x"), y=relative.axis.point(sign.locs[1],"y"), 
       labels="-", font=2, col=richness.col, cex=1.2)
  
  text(x=relative.axis.point(0.38, "x"), 
       y=relative.axis.point(0.3, "y"), 
       srt=45, labels ="Coordinated", col=richness.col)
  
  text(x=relative.axis.point(0.3, "x"), 
       y=relative.axis.point(0.62, "y"), 
       srt=-45, labels ="Ratio", col=richness.col)
  
  x.scores <- seq(par("usr")[1], par("usr")[2], len=50)
  y.scores <- seq(par("usr")[3], par("usr")[4], len=50)
  test.scores <- expand.grid(x.scores, y.scores)
  alpha.diag <- rotate.data(test.scores[,1], test.scores[,2], 45)
  
  text(x=relative.axis.point(0.065, "x"),
       y = relative.axis.point(0.95, "y"),
       labels="(A)", font=2)
  
  box()
  
  #                    Alpha rotated ####
  
  plot(x=NULL, y=NULL, xlim=range(alpha.diag[,"x.pos"]), ylim=range(alpha.diag[,"y.pos"]),
       xlab="", ylab="", axes=FALSE)
  
  axis(side=1, mgp=c(3,0.2,0), at=seq(2,8,1))
  mtext(side=1, line=1.35, text="Coordinated rotated axis")
  
  axis(side=2)
  mtext(side=2, line=1.5, text="Ratio rotated axis", las=0)
  
  origin <- c(relative.axis.point(0.5, "x"), 
              relative.axis.point(0.5, "y"))
  arrow.x <- relative.axis.point(arrow.length, "x") - par("usr")[1]
  arrow.y <- relative.axis.point(arrow.length, "y") - par("usr")[3]
  arrow.x <- origin[1] + arrow.x * cos(seq(0,359,45) * (pi/180))
  arrow.y <- origin[2] + arrow.y * sin(seq(0,359,45) * (pi/180))
  
  points(x=arrow.x, y=arrow.y, pch=16, col="grey70")
  
  segments(x0=origin[1], y0=origin[2], x1=arrow.x, y1 = arrow.y,
           col="grey70", lty="31")
  
  text(x=arrow.x, y= arrow.y, labels=1:8, pos=c(1,1,2,1,1,1,2,1), 
       col="grey60", cex=0.8, adj=0.5, offset=0.4)
  
  fade.a <- relative.axis.point(fade.rad, "x") - par("usr")[1]
  draw.ellipse(x=relative.axis.point(0.5,"x"),
               y=relative.axis.point(0.5,"y"),
               a=fade.a*1.9, b=fade.a, 
               border=NA, col=rgb(1,1,1,0.85))
  # dominant arrows
  Arrows(x0=relative.axis.point(0.5, "x"), 
         y0=relative.axis.point(0.5, "y"),
         x1=relative.axis.point(c(0.5,0.5,lower.inner,upper.inner), "x"),
         y1=relative.axis.point(c(lower.inner,upper.inner,0.5,0.5), "y"),
         arr.type="triangle", arr.length=0.1, arr.width=0.1)
  
  text(x=relative.axis.point(c(0.5,0.5,lower.inner,upper.inner), "x"),
       y=relative.axis.point(c(lower.inner,upper.inner,0.5,0.5), "y"),
       labels=c("Func > Tax", "Tax > Func", 
                "Coord. -", "Coord. +"),
       pos=c(1,3,2,4), cex=0.8, offset=0.4)
  
  box()
  text(x=relative.axis.point(0.065, "x"),
       y = relative.axis.point(0.95, "y"),
       labels="(B)", font=2)
  
  #                       Close plot ####
  
  dev.off()
  
  
  
}