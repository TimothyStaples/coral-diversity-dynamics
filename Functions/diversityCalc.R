alphaDiversityCalc <- function(siteData, 
                               genusMatProp,
                               grMatProp){
  
  tempGen <- as.list(as.data.frame(t(genusMatProp)))
  tempGen <- lapply(tempGen, function(x){x[x>0]})
  tempGr <- as.list(as.data.frame(t(grMatProp)))
  tempGr <- lapply(tempGr, function(x){x[x>0]})
  
  tempDf <- data.frame(transect = rownames(genusMatProp),
                       genusH0 = sapply(tempGen, renyi, hill=TRUE, scales=0),
                       genusH1 = sapply(tempGen, renyi, hill=TRUE, scales=1),
                       genusH2 = sapply(tempGen, renyi, hill=TRUE, scales=2),
                       grH0 = sapply(tempGr, renyi, hill=TRUE, scales=0),
                       grH1 = sapply(tempGr, renyi, hill=TRUE, scales=1),
                       grH2 = sapply(tempGr, renyi, hill=TRUE, scales=2))

  extraDf <- do.call("rbind",
                         lapply(split(siteData,
                                      f=siteData$locality),
                                function(x){
                                  
                                  x <- x[order(x$transect.age, decreasing=TRUE),]
                                  x$date.diff <- c(NA, x$transect.age[-nrow(x)] - x$transect.age[-1])
                                  x$top.date.diff <- x$transect.age[1] - x$transect.age
                                  
                                  return(x)
                                }))
  
  divDf <- merge(tempDf, extraDf,
                 by.x="transect", by.y="transect",
                 all.x=TRUE, all.y=FALSE, sort=FALSE)
  
  colnames(divDf)[colnames(divDf)=="transect.age"] = "pred.date"
  
  return(divDf)
  
}
