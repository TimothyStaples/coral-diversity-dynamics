point.arc<-function(x0, x1, y0, y1, offset, position="above", flip=FALSE, ...){
  
  # FUNCTION PURPOSE #
  # Create a circle segment (arc) between two points. Function returns the x and
  # y points along the arc, with 1 point for every degree. These data can be fed
  # straight into points() to plot the arc.
  
  # The function works by calculating the origin of an imaginary circle between
  # the two points, with a radius given by the distance and the offset arguments.
  # It then calculates points along the circumference of that circle, using the
  # two points provided as start and stop points.
  
  # This function creates a symmetrical arc even if the plot window is not square,
  # or the axes are on different scales. To do this, it needs to scale the circle
  # to the dimensions of the plot device, so it should be called after opening a 
  # plot device, as it requires some information from the plot device to make sure
  # the arc is symmetrical.
  
  # REQUIRED ARGUMENTS #
  # x0, x1, y0, y1: Plot points to start and end arc, as per segments().
  
  # offset:   Indicates the magnitude of the radius of the circle the the arc is
  #           plotted along. Bigger numbers mean will be offset. Bigger numbers 
  #           mean flatter arc. A single number greater or equal to 0. An offset
  #           of 0 will plot an arc where the circle diameter is the distance
  #           between the two points.
  
  # OPTIONAL ARGUMENTS #
  
  # position: Tells the function whether the arc should run above or below the 
  #           given points. A character string of "above" or "below",
  
  # flip:     When the angle between the two points is ~180 degrees, the arc 
  #           sometimes goes the wrong way. Set this to TRUE to flip
  #           the arc (mirror image it).
  
  # ERROR CHECKING #
  
  # are points numbers?
  if(!is.numeric(c(x0,x1,y0,y1))){
    stop("Start and end points must be numbers.")
  }
  
    # Offset must be zero or greater.
  if(sum(offset >=0) < length(offset)){
    stop("offset must be positive.")
  }
  
  # FUNCTION CODE #
  
  # because these arcs might be plot on axes with different units, the first thing
  # we need to do is convert our points to a proportional score across the axis
  # which we can then back-transform at the end
  
  # get size of plot window in cm
  plot.window <- dev.size(units="cm")
  
  # get axis limits in plot units
  x.axis<-par("usr")[1:2]
  y.axis<-par("usr")[3:4]
  
  # transform each point into cm scores, using the origin of the plot as our
  # 0cm point.
  transformed.points <- cbind( (x0-x.axis[1]) / (x.axis[2]-x.axis[1]) * plot.window[1],
                               (x1-x.axis[1]) / (x.axis[2]-x.axis[1]) * plot.window[1],
                               (y0-y.axis[1]) / (y.axis[2]-y.axis[1]) * plot.window[2],
                               (y1-y.axis[1]) / (y.axis[2]-y.axis[1]) * plot.window[2])
  
  # now we have the points in the correct units, we can calculate our arc 
  
  # get mid point between two points
  if(dim(transformed.points)[1] > 1){
  mid.point<-cbind(rowMeans(transformed.points[,c(1,2)]),
                   rowMeans(transformed.points[,c(3,4)]))
  } else {
    mid.point<-cbind(mean(transformed.points[,c(1,2)]),
                     mean(transformed.points[,c(3,4)]))
  }
               
  # now use some right-angled triangle trig to get coordinates of
  # offset mid point
  
  # 1. find slope perpendicular to slope from point1 (or 2) to mid.point
  slope.b <- -1 / ((transformed.points[,4] - mid.point[,2]) /
                  (transformed.points[,2] - mid.point[,1]))
  
  # 2. get position of circle origin to calculate arc, based on whether you
  # want the arc to be above or below the 2 connecting points. If arc should be
  # above, offset point needs to be below, and vice versa.
  if(position=="above"){
    offset.point <- mid.point - matrix(offset*c(1/sqrt(1+slope.b^2),
                                        slope.b / sqrt(1+slope.b^2)),
                                      ncol=2)
  }
  
  if(position=="below"){
    offset.point<- mid.point + matrix(offset*c(1/sqrt(1+slope.b^2),
                                               slope.b / sqrt(1+slope.b^2)),
                                      ncol=2)
  }
  
  # we need radius from circle origin to the points
  radius<-sqrt((transformed.points[,1]-offset.point[,1])^2 + 
               (transformed.points[,2]-offset.point[,2])^2)
  
  # and angles from circle origin to each original point
  angles<-cbind(atan2(transformed.points[,3] - offset.point[,2],
                      transformed.points[,1] - offset.point[,1]),
                atan2(transformed.points[,4] - offset.point[,2],
                      transformed.points[,2] - offset.point[,1]))
  
  angles.deg<-ifelse(angles<0,
                     angles*180/pi + 360,
                     angles*180/pi)
  
  # calculate points along arc using atan2 (default is 1 degree increments)
  # calculate in both directions, and use the smallest vector (the shortest
  # distance between points)
  
  # we need to rename our start and end degrees into max and min so our sequence
  # always runs from small degrees to large degrees
  min.deg<-apply(angles.deg, 1, min)
  max.deg<-apply(angles.deg, 1, max)
  
  if(position=="above" & !flip |
     position=="below" & flip){
    deg.vect<-mapply(mind=min.deg, maxd=max.deg, 
                       function(mind, maxd){seq(mind, maxd, length.out=200)})
  }
  
  if(position=="below" & !flip |
     position=="above" & flip){
    
    deg.vect<-mapply(mind=min.deg, maxd=max.deg, 
                       function(mind, maxd){c(seq(mind, 360, length.out=100),
                                              seq(1, maxd, length.out=100))})
  }

  # this converts points back to device units (cms). arc.points is a list of two
  # components: x scores for each arc, and y scores for each arc.
  vect=as.data.frame(deg.vect)[,1]
  radius=radius[1]
  offset=as.data.frame(t(offset.point))[,1]
  
  arc.points<-mapply(vect=as.data.frame(deg.vect),
                     radius=radius,
                     offset=as.data.frame(t(offset.point)),
         function(vect, radius, offset){
        
           temp.points<-cbind(
             offset[1] + radius * cos((vect * pi) / (180)),
             offset[2] + radius * sin((vect * pi) / (180))
             )
         
           # now we need to convert our points along the arc back into plot units
           arc.x<-(temp.points[,1] / plot.window[1]) * 
                  (x.axis[2]-x.axis[1]) + x.axis[1]
           arc.y<-(temp.points[,2] / plot.window[2]) * 
                  (y.axis[2]-y.axis[1]) + y.axis[1]
           
           points(x=arc.x, y=arc.y, type="l")
           
           return(cbind(arc.x, arc.y))
           }, SIMPLIFY=FALSE)

  invisible(arc.points)
}