alphaDiversityCalc <- function(siteData, 
                               genusMatProp,
                               grMatProp){
  
  tempGen <- as.list(as.data.frame(t(genusMatProp)))
  tempGen <- lapply(tempGen, function(x){x[x>0]})
  tempGr <- as.list(as.data.frame(t(grMatProp)))
  tempGr <- lapply(tempGr, function(x){x[x>0]})
  
  tempDf <- data.frame(transect = rownames(genusMatProp),
                       genusH0 = sapply(tempGen, hillCalc, l=1),
                       genusH1 = sapply(tempGen, hillCalc, l=1e-6),
                       genusH2 = sapply(tempGen, hillCalc, l=-1),
                       grH0 = sapply(tempGr, hillCalc, l=1),
                       grH1 = sapply(tempGr, hillCalc, l=1e-6),
                       grH2 = sapply(tempGr, hillCalc, l=-1))

  h0 <- rotate.data(x=tempDf$grH0,
                    y=tempDf$genusH0, angle=45)[,3:4]
  colnames(h0) = c("coordH0", "ratioH0")
  
  h1 <- rotate.data(x=tempDf$grH1,
                    y=tempDf$genusH1, angle=45)[,3:4]
  colnames(h1) = c("coordH1", "ratioH1")
  
  h2 <- rotate.data(x=tempDf$grH2,
                    y=tempDf$genusH2, angle=45)[,3:4]
  colnames(h2) = c("coordH2", "ratioH2")
  
  extraDf <- do.call("rbind",
                         lapply(split(siteData,
                                      f=siteData$locality),
                                function(x){
                                  
                                  x <- x[order(x$transect.age, decreasing=TRUE),]
                                  x$date.diff <- c(NA, x$transect.age[-nrow(x)] - x$transect.age[-1])
                                  x$top.date.diff <- x$transect.age[1] - x$transect.age
                                  
                                  return(x)
                                }))
  
  divDf <- cbind(tempDf, h0, h1, h2)
  
  divDf <- merge(divDf, extraDf,
                 by.x="transect", by.y="transect",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)
  
  colnames(divDf)[colnames(divDf)=="transect.age"] = "pred.date"
  
  return(divDf)
  
}
