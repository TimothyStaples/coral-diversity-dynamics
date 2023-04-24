siteMapPlot <- function(path){
  
  # Rotated site map with age trends ####
  #                 import date data ####
  
  Tpng<-read.csv("/home/timothy/Dropbox/Tim/Post-doc/Research projects/huon_cleaning/raw.datafiles/huon.transect.csv",
                 stringsAsFactors = FALSE)
  Tloc<-Tpng[!duplicated(Tpng$locality) & Tpng$locality %in% huon_coral$locality,]
  Tloc$latitude<-as.numeric(gsub("\\+AC0","",Tloc$latitude))
  Tloc$longitude<-as.numeric(Tloc$longitude)
  Tloc$color <- c("blue", "red", "darkgreen")[as.factor(Tloc$site)]
  Tloc <- Tloc[order(Tloc$latitude),]
  Tloc$site.num <- 1:nrow(Tloc)
  
  #            import and rotate map ####
  
  rotateCoords <- function(crds, angle=0, center= c(min(crds[,1]),min(crds[,2]))) {
    co <- cos(-angle*pi/180)
    si <- sin(-angle*pi/180)
    adj <- matrix(rep(center,nrow(crds)),ncol=2,byrow=TRUE)
    crds <- crds-adj
    cbind(co * crds[,1] - si * crds[,2],
          si * crds[,1] + co * crds[,2]) + adj
  }
  
  png_outline <- readShapeSpatial("/home/timothy/Dropbox/Tim/Post-doc/Research projects/huon_cleaning/raw.datafiles/shape_files/PNG_adm/PNG_adm0.shp")
  png_elev2 <- raster("/home/timothy/Dropbox/Tim/Post-doc/Research projects/huon_cleaning/raw.datafiles/shape_files/PNG_elev/ASTGTM2_S07E147/ASTGTM2_S07E147_dem.tif")
  
  proj4string(png_outline) <- crs(png_elev2)
  
  rotate.angle<-35
  
  locs <- Tloc[,c("latitude", "longitude", "locality", "site.num", "site", "color")]
  locs <- locs[complete.cases(locs),]
  coordinates(locs)<-locs[,c("longitude", "latitude")]
  proj4string(locs) <- proj4string(png_elev2)
  
  png_outline_rotate<-elide(png_outline, rotate=rotate.angle)
  
  png_outline_simp <- crop(png_outline, extent(140,154,-15,-3))
  
  # png_river_rotate<-elide(png_river_smooth, 
  #                         rotate=rotate.angle, 
  #                         bb=bbox(png_outline))
  
  locs_rotate<-rotateCoords(cbind(locs$longitude, locs$latitude),
                            angle=rotate.angle,
                            center=bbox(png_outline)[,1])
  
  northarrow<-readPicture("Plots/North_Pointer_rotateCairo.svg")
  
  #                         Plot map ####
  
  pdf(path, height=3.5, width=5.5)
  par(mar=c(0.5,2,0.5,0.5), ps=6, tck=-0.025, mgp=c(3,0.5,0), las=1)
  
  split.screen(rbind(c(0.055,0.5,0.025,0.975),
                     c(0.255,0.405,0.8,0.975),
                     
                     c(0.54,0.69,0.70,0.90),
                     c(0.69,0.84,0.70,0.90),
                     c(0.84,0.99,0.70,0.90),
                     
                     c(0.54,0.69,0.40,0.60),
                     c(0.69,0.84,0.40,0.60),
                     c(0.84,0.99,0.40,0.60),
                     
                     c(0.54,0.69,0.1,0.30),
                     c(0.69,0.84,0.1,0.30),
                     c(0.84,0.99,0.1,0.30)))
  
  screen(1)
  par(mar=c(0,0,0,0))
  plot(x=NULL, y=NULL, xlim=c(149.45,149.75), ylim=c(-11.3,-10.9), axes=FALSE,
       xaxs="i", yaxs="i", asp=1)
  
  plot.extent<-par("usr")
  
  # outline
  plot(crop(png_outline_rotate, extent(plot.extent)),
       xlim=c(147.7,148.1), ylim=c(-6.7,-6.3),
       border=NA, col="grey80", bg="white", add=TRUE, asp=1)
  
  # plot(crop(png_river_rotate, extent(plot.extent)),
  #      col="grey30", lwd=1, add=TRUE, asp=1)
  
  # sites
  
  # offset several points
  locs_mod <- locs_rotate
  locs_mod[locs$locality == "Kwambu Platforms",] = 
    locs_mod[locs$locality == "Kwambu Platforms",] + c(0.0075, -0.01)
  locs_mod[locs$locality == "Loto Beach",] = 
    locs_mod[locs$locality == "Loto Beach",] + c(0.01, 0.01)
  locs_mod[locs$locality == "NW Point",] = 
    locs_mod[locs$locality == "NW Point",] + c(0.015, -0.01)
  locs_mod[locs$locality == "Midway Cove",] = 
    locs_mod[locs$locality == "Midway Cove",] + c(0.01, -0.0175)
  
  sapply(1:dim(locs_mod)[1], function(x){
    text(y=locs_mod[x,2], 
         x=locs_mod[x,1],
         labels=paste0(locs$site.num[x], ": ", locs$locality[x]), 
         pos=4, cex=0.85, col=locs$color[x],
         offset=ifelse(locs_mod[x,1] != locs_rotate[x,1], 0.1, 0.25))
  })
  
  segments(x0=locs_rotate[,1],
           x1=locs_mod[,1],
           y0=locs_rotate[,2],
           y1=locs_mod[,2], col=locs$color)
  
  box(lwd=2, col="white")
  
  points(locs_rotate[,2] ~ locs_rotate[,1], pch=21, cex=0.7, bg=locs$color)
  
  # site indicators
  
  hub.axis <- locs_rotate[locs$locality %in% 
                            c("Loto Beach", "Midway Cove"),2] + c(-0.025, 0.025)
  axis(side=4, at= hub.axis, labels=NA, tck=0.025, line=-3.15, col="blue")
  mtext(side=4, line=-3.4, at= mean(hub.axis), las=0,
        text="Hubegong", font=2, col="blue")
  
  kanomi.axis <- locs_rotate[locs$locality %in% 
                               c("Pukau NW", "Bonah River"),2] + c(-0.01, 0.01)
  axis(side = 4, at = kanomi.axis, labels=NA, tck=0.025, line=-3.75, col="red")
  mtext(side=4, line=-4, at = mean(kanomi.axis), text="Kanomi", font=2, col="red", las=0)
  
  sialum.axis <- sort(locs_rotate[locs$locality %in% 
                                    c("Kilisairo River NW", "Kwambu Platforms"),2]) + c(-0.015, 0.025)
  axis(side=4, at = sialum.axis, labels=NA, tck=0.025, line=-3, col = "darkgreen")
  mtext(side = 4, line = -3.25, at = mean(sialum.axis), text="Sialum", font=2, col = "darkgreen", las=0)
  
  # AXIS 1
  
  # rotate points along the bottom axis of the original plot
  axis1.seq<-seq(147.5, 147.85, 0.05)
  axis1.thres<- -6.335
  
  bottom.ends<-rotateCoords(cbind(axis1.seq, axis1.thres),
                            angle=rotate.angle, center=bbox(png_outline)[,1])
  
  bottom.ticks<-rotateCoords(cbind(axis1.seq, -6.34),
                             angle=rotate.angle, center=bbox(png_outline)[,1])
  
  # create white space in front of plots
  polygon(x=c(par("usr")[1], bottom.ends[,1], par("usr")[1]),
          y=c(par("usr")[3], bottom.ends[,2], par("usr")[3]),
          border=NA, col="white")
  
  # add line on
  segments(x0=bottom.ends[1,1], x1=bottom.ends[dim(bottom.ends)[1],1],
           y0=bottom.ends[1,2], y1=bottom.ends[dim(bottom.ends)[1],2])
  
  # add axis tick marks
  segments(x0=bottom.ends[,1], x1=bottom.ticks[,1],
           y0=bottom.ends[,2], y1=bottom.ticks[,2])
  
  # add axis text
  par(xpd=NA)
  text(x=bottom.ticks[4:7,1], y=bottom.ticks[4:7,2]-0.005, pos=2,
       labels=paste0(sprintf("%#.2f", axis1.seq[4:7]), "°E"),
       offset=0.2)
  par(xpd=FALSE)
  
  # AXIS 2
  axis2.seq<-seq(-6.3, -6, 0.05)
  
  side.ends<-rotateCoords(cbind(147.55, axis2.seq),
                          angle=rotate.angle, center=bbox(png_outline)[,1])
  side.ticks<-rotateCoords(cbind(147.545, axis2.seq),
                           angle=rotate.angle, center=bbox(png_outline)[,1])
  
  # create white space in front of plots
  polygon(x=c(par("usr")[1], side.ends[,1], par("usr")[1]),
          y=c(par("usr")[4], side.ends[,2], par("usr")[4]),
          border=NA, col="white")
  
  segments(x0=side.ends[1,1], x1=side.ends[dim(side.ends)[1],1],
           y0=side.ends[1,2], y1=side.ends[dim(side.ends)[1],2])
  
  segments(x0=side.ends[,1], x1=side.ticks[,1],
           y0=side.ends[,2], y1=side.ticks[,2])
  
  # add axis text
  par(xpd=NA)
  text(x=side.ticks[3:6,1], y=side.ticks[3:6,2]+0.005, pos=2,
       labels=paste0(sprintf(axis2.seq[3:6], fmt="%#.2f"), "°S"),
       offset=0.2)
  
  
  segments(x0=relative.axis.point(0.002, "x"), 
           x1=relative.axis.point(0.002, "x"),
           y0=relative.axis.point(0.601, "y"), 
           y1=relative.axis.point(0.271, "y"))
  
  par(xpd=FALSE)
  
  # SCALE BAR & NORTH ARROW
  
  grid.picture(northarrow, x=0.145,y=0.7425, distort=FALSE, width=0.035)
  
  par(xpd=NA)
  scalebar(d=5, xy=c(relative.axis.point(0.05, "x"),
                     relative.axis.point(0.35, "y")), 
           lonlat=TRUE, label=c("","5 km", ""), lwd=0.5)
  par(xpd=FALSE)
  
  text(x=relative.axis.point(0, "x"),
       y = relative.axis.point(0.975, "y"),
       labels = paste0("(A)"),
       adj=0, font =2)
  
  close.screen(1)
  
  screen(2)
  par(mar=c(0,0,0,0), ps=6, tcl=-0.125, mgp=c(3,0.5,0), las=1)
  
  png_extent<-as.vector(extent(png_outline_simp))
  
  plot(png_outline_simp, lwd=0.5, col="grey80", border=NA,
       xlim=c(142.5,152), ylim=c(-11,-3))  
  
  text(x=143.5, y=-6.25, label="PNG", font=2, col="grey45")
  
  box()
  
  
  axis(side=1, at=seq(140,155,5), 
       labels=paste0(seq(140,155,5), "°E"), mgp=c(3,-0.2,0))
  axis(side=2, at=seq(-10,0,5), labels=paste0(seq(-10,0,5), "°S"),
       mgp=c(3,0.2,0))
  
  arrows(x1=mean(locs$longitude)+0.25,
         y1=mean(locs$latitude)-0.25,
         x0=150, y0=-9, lwd=2,
         length=0.05)
  
  close.screen(2)
  
  date.order <- c("Kwambu Platforms", "Kilasairo River SE", "Kilasairo River NW",
                  "Pukau NW", "Paradise Springs", "Bonah River",
                  "Midway Cove", "NW Point", "Loto Beach")
  
  raw.dates <- read.csv("./raw.datafiles/raw_huon_dates.csv")
  
  sapply(1:length(date.order), function(n){
    
    test.ylims <- rbind(c(8750,5750),
                        c(8800,6000),
                        c(9550,7250))[(n-1) %/% 3 + 1,]
    
    screen(n+2)
    par(mar=c(0,0,0,0), las=0, ps=6, tcl=-0.125, mgp=c(3,0.5,0), las=1)
    
    temp.loc <- date.order[n]
    date.locs <- raw.dates[raw.dates$locality == temp.loc,]
    date.locs <- date.locs[date.locs$sigma < 1000, ]
    site.locs <- raw.dates[raw.dates$site == date.locs$site[1],]
    
    temp.col <- c("darkgreen", "red", "blue")[((n-1) %/% 3) + 1]
    
    poly.col <- col2rgb(temp.col)/255
    
    plot(x=NULL, y=NULL, xlim=test.ylims, ylim=c(-5,14), axes=FALSE, xlab="", ylab="")
    
    if(n %in% c(1,4,7)){
      axis(side=2, at=seq(-5,15,5), mgp=c(3,0.2,0))
      mtext(side=2, line=1, text="Transect height (m)", las=0)
      } else {
      axis(side=2, at=seq(-5,15,5), labels=NA)
    }
    
    if(n %in% c(2,5,8)){
      mtext(side=1, line=0.325, text="Thousands of years before present")
      
    }
    axis(side=1, at = seq(6000,9000,1000), labels=6:9, mgp=c(3,-0.2,0))
    axis(side=1, at = seq(6000,9000,500), labels=NA, tcl=-0.125)
    
    
    segments(y0 = date.locs$height.from.dist, y1 = date.locs$height.from.dist,
             x0 = date.locs$date + date.locs$sigma, 
             x1 = date.locs$date - date.locs$sigma, col="grey70")
    
    points(date.locs$height.from.dist ~ date.locs$date, pch=16, cex=0.5, col="grey70")
    
    bacon.dirs <- paste0("/home/timothy/Dropbox/Tim/Post-doc/Research projects/huon_cleaning/Bacon_runs/",
                         date.order[n])
    bacon.files <- list.files(bacon.dirs, pattern="_ages")
    bacon.ages <- read.table(paste0(bacon.dirs, "/", bacon.files), header=TRUE)
    bacon.ages$height <- -((bacon.ages$depth / 100) - max(date.locs$height.from.dist))
    
    polygon(y=c(bacon.ages$height, rev(bacon.ages$height)),
            x = c(bacon.ages$max, rev(bacon.ages$min)),
            col=rgb(poly.col[1], poly.col[2], poly.col[3], 0.3),
            border=NA)
    
    lines(y=bacon.ages$height, x=bacon.ages$mean, lwd=1, col=temp.col)
    
    sub.tran <- huon_site[huon_site$locality == temp.loc, ]
    
    # tran.width <- abs(relative.axis.point(0.03, "y") - relative.axis.point(0, "y"))
    # segments(x0=sub.tran$height.from.dist,
    #          x1=sub.tran$height.from.dist,
    #          y0=sub.tran$transect.age + tran.width,
    #          y1=sub.tran$transect.age - tran.width,
    #          col=temp.col)
    
    text(x=relative.axis.point(0.02, "x"),
         y = relative.axis.point(0.915, "y"),
         labels = paste0("(", LETTERS[n+1], "): Site ", 
                         locality.number$number[locality.number$locality == temp.loc]),
         adj=0, font =2)
    box()
    
    close.screen(n+2)
    
    
  })
  close.screen(all.screens=TRUE)
  
  dev.off()
  
  
}