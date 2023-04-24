betaDiversityCalc <- function(siteData, 
                              genusMatProp,
                              grMatProp, method){
  
  rownames(genusMatProp) = siteData$transect
  rownames(grMatProp) = siteData$transect
  
  turnover.df <- do.call("rbind",
                         lapply(split(siteData,
                                      f=siteData$locality),
                                function(x){

                                  x <- x[order(x$transect.age, decreasing=TRUE),]
                                  
                                  gen.loc.mat <- genusMatProp[match(as.character(x$transect),
                                                                    rownames(genusMatProp)),]
                                  gen.loc.mat <- gen.loc.mat[, colSums(gen.loc.mat) >0]
                                  gen.dist.mat <- as.matrix(vegdist(gen.loc.mat, method=method))
                                  
                                  gr.loc.mat <- grMatProp[match(as.character(x$transect),
                                                                  rownames(grMatProp)),]
                                  gr.loc.mat <- gr.loc.mat[, colSums(gr.loc.mat) >0]
                                  gr.dist.mat <- as.matrix(vegdist(gr.loc.mat, method=method))
                                  
                                  x$top.gen.diss <- c(NA, gen.dist.mat[1,-1])
                                  x$seq.gen.diss <- c(NA, diag(gen.dist.mat[-1,-dim(gen.dist.mat)[1]]))
                                  x$top.gr.diss <- c(NA, gr.dist.mat[1,-1])
                                  x$seq.gr.diss <- c(NA, diag(gr.dist.mat[-1,-dim(gr.dist.mat)[1]]))
                                  
                                  x$top.gen.diss.diff <- x$top.gen.diss - x$top.gen.diss[2]
                                  x$top.gr.diss.diff <- x$top.gr.diss - x$top.gr.diss[2]
                                  x$seq.gen.diss.diff <- x$seq.gen.diss - x$seq.gen.diss[2]
                                  x$seq.gr.diss.diff <- x$seq.gr.diss - x$seq.gr.diss[2]
                                  
                                  x$date.diff <- c(NA, x$transect.age[-nrow(x)] - x$transect.age[-1])
                                  x$top.date.diff <- x$transect.age[1] - x$transect.age
                                  
                                  return(x)
                                }))
  
  top.beta.rotate <- rotate.data(x=turnover.df$top.gr.diss,
                                 y=turnover.df$top.gen.diss,
                                 angle=45)[,3:4]
  colnames(top.beta.rotate) <- c("dissimilarity", "conservation")
  
  turnover.df <- cbind(turnover.df, top.beta.rotate)
  
  seq.beta.rotate <- rotate.data(x=turnover.df$seq.gr.diss,
                                 y=turnover.df$seq.gen.diss,
                                 angle=45)[,3:4]
  colnames(seq.beta.rotate) <- c("turnover", "persistance")
  
  colnames(turnover.df)[colnames(turnover.df)=="transect.age"] = "pred.date"
  
  return(cbind(turnover.df, seq.beta.rotate))
  
}