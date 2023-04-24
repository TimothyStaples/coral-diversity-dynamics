ts.diversity <- function(data, rare.sample){
  
  require(vegan)
  
  # data needs to have four columns, in the following order:
  # (1) ts identifier 
  # (3) ages
  # (4) taxon
  # (5) taxa abundance
  
  colnames(data) = c("ts", "age", "taxa", "abund")
  
  ts <- unique(data[,1])
  
  ts.list <- lapply(ts, function(x){
    
    # subset just data from one ts
    sub.df <- data[data$ts == x, ]
    
    # remove any NA or 0 abundance data
    sub.df <- sub.df[!(is.na(sub.df$taxa) | 
                     is.na(sub.df$abund) |
                     sub.df$abund == 0), ]
    
    # aggregate counts for each taxa in each age bracket
    sub.df$cell.ident <- paste0(sub.df$age, ".", sub.df$taxa)
    agg.df <- do.call("rbind", lapply(split(sub.df, f=sub.df$cell.ident),
                                      function(y){
                                      
        return(data.frame(ts = y$ts[1],
                          age = y$age[1],
                          taxa = y$taxa[1],
                          cell.ident = y$cell.ident[1],
                          abund = sum(y$abund)))  
                                        
                                      }))
    
    # create empty site-species matrix to fill with abundances
    #     rows are 'ages', reprenting segments of time-series
    #     columns are 'taxa'
    sub.ages <- sort(unique(agg.df$age), decreasing=TRUE)
    sub.taxa <- sort(unique(agg.df$taxa))
    
    exp <- expand.grid(ages = sub.ages, taxa = sub.taxa)
    
    exp$abund <- agg.df$abund[match(paste0(exp$ages, ".", exp$taxa),
                                    agg.df$cell.ident)]
    exp$abund[is.na(exp$abund)] = 0 
    
    ssmat <- matrix(exp$abund, 
                    nrow=length(sub.ages),
                    ncol = length(sub.taxa),
                    dimnames = list(sub.ages, sub.taxa))
    
    prop.ssmat <- prop.table(ssmat, 1)

    # now calculate diversity at each time.step
    ts.alpha <- NA #rarefy(ssmat, rare.sample)
    
    # beta diversity comes from dissim matrix
    ts.dissim <- as.matrix(vegdist(prop.ssmat, method="jaccard"))

    ts.top.beta <- c(NA, ts.dissim[-1,1])
    ts.bottom.beta <- c(ts.dissim[-nrow(ts.dissim),ncol(ts.dissim)], NA)
    ts.seq.beta <- diag(ts.dissim[-1,-ncol(ts.dissim)])
    
    return.data <- agg.df
    return.ssmat <- prop.ssmat
    return.div <- data.frame(ts = agg.df$ts[1],
                             age = as.numeric(rownames(ssmat)),
                             alpha = ts.alpha,
                             t.beta = ts.top.beta,
                             b.beta = ts.bottom.beta,
                             s.beta = c(NA, ts.seq.beta))
    return.div$age.lag = c(NA, abs(diff(return.div$age)))
    return.div$age.from.top <- return.div$age[1] - return.div$age
    
    return(list(data=return.data,
                raw.ssmat = ssmat,
                prop.ssmat=prop.ssmat,
                div.df = return.div))  
    
    
  })
  
  return.list <- list(agg.data = do.call("rbind", lapply(ts.list, function(x){x$data})),
                      div.data = do.call("rbind", lapply(ts.list, function(x){x$div.df})),
                      raw.ssmat = lapply(ts.list, function(x){x$raw.ssmat}),
                      prop.ssmat = lapply(ts.list, function(x){x$prop.ssmat}))

  return(return.list)  
}