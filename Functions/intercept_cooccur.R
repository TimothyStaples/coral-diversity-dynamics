intercept_cooccur<-function(huon, species){
  
  ## AUTHOR: Timothy Staples
  ## DATE: 19/03/2018
  
  ## FUNCTION PURPOSE:
  ## Takes a dataset of fossil coral transects from Huon Peninsula
  ## - calculates observed interactions between "species"
  ## - randomly shuffles intercepts along transects to record null interactions
  ## - calculates "species" with significantly higher or lower interaction rates
  
  ## ARGUMENTS:
  ## huon:      dataframe of fossil coral transect intercepts
  ## species:   vector created from dataframe of species to use to calculate
  ##            interactions
  
  # 1. Calculate observed interactions ####
  
  # add species to data-frame and get unique levels
  sp.levels<-levels(as.factor(species))
  huon$int.sp<-as.factor(species)
  
    obs_spsp_mat_list<-lapply(split(huon, 
                             f=paste0(huon$transect, ":", huon$sub.tr)),
                       function(x){
                         
      # order measurements
      x<-x[order(x$tape),]
      # create position vector
      x$pos<-1:dim(x)[1]
        
      # overwrite intercepts for non-whole coral
      x$int.sp[x$whole.bin==0]=NA
                     
      # two columns with each row = two adjacent intercepts
      temp.df<-data.frame(x$int.sp[match(1:(max(x$pos)-1), x$pos)],
                          x$int.sp[match(2:max(x$pos), x$pos)])
      colnames(temp.df)=NULL
      
      # table interactions and return species x species matrix
      as.matrix(table(x$int.sp[match(1:(max(x$pos)-1), x$pos)],
                                x$int.sp[match(2:max(x$pos), x$pos)]))
      
      })
    
    obs_spsp_mat<-Reduce("+", obs_spsp_mat_list) 
    
    obs_pcorr<-partial_correlation(as.data.frame.matrix(obs_spsp_mat), method="correlation")
  
    # 2. Calculate expected interactions ####
    
    # 2.1 create intercept x species matrix for each sub-transect ####
    intercept_sp_mat<-lapply(split(huon, f=paste0(huon$transect,
                                               ":",
                                               huon$sub.tr)),
                       function(subtr){
                         
                         # make position column
                         subtr<-subtr[order(subtr$tape),]
                         subtr$pos<-1:dim(subtr)[1]
                         
                         # make empty intercept x species matrix
                         temp.mat<-matrix(0, 
                                          nrow=length(sp.levels),
                                          ncol=dim(subtr)[1],
                                          dimnames=list(sp.levels,
                                                        subtr$pos))
                         
                         # ignore non-whole coral
                         subtr$int.sp[subtr$whole.bin==0]=NA
                         
                         temp.mat[cbind(subtr$int.sp,
                                        subtr$pos)]=1
                         
                         # remove non-present species
                         temp.mat<-temp.mat[rowSums(temp.mat)>0,]
                         
                         return(temp.mat)
                         
                       }) 
    
      # 2.2 create array of random communities  ####
    
        # start by randomising intercepts along transect
    int_sp_permlist<-lapply(1:length(intercept_sp_mat), function(mat){
      print(mat)
      mat<-intercept_sp_mat[[mat]]
      
      # empty species x intercept matrix (with all species)
      storage_array<-array(0,
                           dim=c(length(sp.levels),
                                 dim(mat)[2],
                                 999),
                           dimnames=list(sp.levels,
                                         1:dim(mat)[2],
                                         1:999))
      
      perms<-permatswap(mat, times=999)$perm
      
      for(i in 1:999){
        
        storage_array[match(rownames(perms[[i]]),
                            rownames(storage_array[,,i])),
                      ,
                      i] = perms[[i]]
        
        rownames(perms[[i]])[1]
        perms[[i]][1,]
        
      }
        
      return(storage_array)
       
    })

    # re-add missing rows to each intercept x species matrix
    
          # 2.3 sum each transect's interactions for a given permutation ####
    perm.mats<-lapply(int_sp_permlist, function(permlist){permlist[,,1]})

    spsp_perm_array<-array(unlist(lapply(1:999, function(perm){
      print(perm)
      
      # extract the matrix slice of a permutation from each array
      perm.mats<-lapply(int_sp_permlist, function(permlist){permlist[,,perm]})
    
      # convert intercept x species array into a species x species array. 

     Reduce("+", lapply(perm.mats, function(slice){
        
        index<-which(slice==1, arr.ind=TRUE)
        
        temp.df<-as.data.frame(cbind(rep(NA, max(index[,2])),
                                     rep(NA, max(index[,2]))))
        
        temp.df[index[,2]+1,1] <- rownames(slice)[index[,1]]
        
        temp.df[index[,2],2] <- rownames(slice)[index[,1]]
        
        temp.df[,1]<-factor(temp.df[,1], levels=sp.levels)
        temp.df[,2]<-factor(temp.df[,2], levels=sp.levels)
        
        as.matrix(table(temp.df))
      }))
      
    })),
    
    dim=c(length(sp.levels), length(sp.levels), 999),
    dimnames=list(sp.levels, sp.levels, 1:999))
    # each matrix slice in this object represents a single permutation


          # 2.4 Calculate partial correlations for each permutation.
    
    # We need to have a look here at errors, potentially from trying to 
    # estimate the partial correlation of species that DO NOT interact in the
    # permutation (or too many perhaps). We could try trimming this down a
    # spot.
    spsp_perm_pcorr<-array(
      
      #apply(spsp_perm_array, 3, function(perm){
      sapply(1:dim(spsp_perm_array)[3], function(perm){
        
        print(perm)
        perm<-spsp_perm_array[,,perm]
      
      partial_correlation(perm, method="shrinkage")
        
        }),
      
      dim=c(length(sp.levels), length(sp.levels), 999),
      dimnames=list(sp.levels, sp.levels, 1:999))
    
    # Then we need four things: 
    
    # absolute differences between expected and observed pcorrs (with negatives
    # converted to be positives)
    abs_pcorr_diffs<- -1 * sweep(abs(spsp_perm_pcorr),  # absolute null pcorrs
                                 c(1,2),  # for each matrix slice
                                 abs(obs_pcorr), # absolute observed pcorrs
                                 "-") # difference
    
    # mean pcorr of null sims
    mean_perm_pcorr<-apply(spsp_perm_pcorr, c(1,2), mean)
    
    # sd of pcorr of null sims
    sd_perm_pcorr<-apply(spsp_perm_pcorr, c(1,2), sd)
    
    # mean effect size, standardised by sd of permutations
    ses_pcorr <- (obs_pcorr - mean_perm_pcorr) / sd_perm_pcorr
    
    # get proportion of obs - exp differences.
    pcor_pvalues <- apply(abs_pcorr_diffs, c(1, 2), function(x) {
      mean(x < 0)
    })
    
    # diagonals are intraspecies, which are converted to NAs
    #diag(pcor_pvalues) <- NA
    
    # adjust p-values for multiple comparison
    pcor_pvalues_adjusted <- matrix(p.adjust(pcor_pvalues, 
                                             method = "fdr"), 
                                    nrow = length(sp.levels), 
                                    ncol = length(sp.levels),
                                    dimnames=list(sp.levels, sp.levels))
    
    return(list(adjusted_pvalues = pcor_pvalues_adjusted,
                raw_pvalues = pcor_pvalues,
                ses = ses_pcorr,
                obs_pcorr = obs_pcorr,
                perm_pcorr = spsp_perm_pcorr))
    
}