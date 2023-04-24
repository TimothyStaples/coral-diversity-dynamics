make.boundary.mat<-function(huon, species){
  
  # 1. Calculate observed interactions ####
  
  # add species to data-frame and get unique levels
  sp.levels<-levels(as.factor(species))
  huon$int.sp<-as.factor(species)
  
  # get genus presence at each intercept boundary in each transect
  obs_spint_mat_list<-lapply(split(huon, 
                             f=paste0(huon$transect, ":", huon$sub.tr)),
                       function(x){
      # order measurements
      x<-x[order(x$tape),]
      # create position vector
      x$pos<-1:dim(x)[1]
        
      # overwrite intercepts for non-whole coral
      x$int.sp[x$whole !="W"]=NA
                     
      # two columns with each row = two adjacent intercepts
      temp.df<-data.frame(x$int.sp[match(1:(max(x$pos)-1), x$pos)],
                          x$int.sp[match(2:max(x$pos), x$pos)])
      colnames(temp.df)=NULL
      
      # set up empty matrix
      
      empty.mat <- matrix(0, nrow=length(sp.levels), ncol=nrow(temp.df),
                          dimnames=list(as.character(sp.levels), 
                                        1:nrow(temp.df)))
      for(i in 1:dim(temp.df)[1]){
      
        if(!is.na(temp.df[i,1])){
        empty.mat[unlist(temp.df[i,1]),i] = 1
        
        }
        
        if(!is.na(temp.df[i,2])){
          empty.mat[unlist(temp.df[i,2]),i] = 1
        }
        
      }
    
      return(empty.mat)
      
      })
  
  # remove empty intercept boundaries (where there are 0 whole coral touching)
  obs_spint_mat_list.z <- lapply(obs_spint_mat_list, function(mat){
    mat[, colSums(mat)==2]
    })

  #obs_spint_mat <- do.call("cbind", obs_spint_mat_list)
  obs_spint_mat.z <- do.call("cbind", obs_spint_mat_list.z)
  
  obs_spint_mat.z <- obs_spint_mat.z[rowSums(obs_spint_mat.z)>0, ]
  
  return(obs_spint_mat.z)
}