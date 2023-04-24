# in this revised version, I realised that I was running
# my partial correlations on a n x n species-species matrix,
# not an n x m species-site matrix. This is a problem because
# it's not how Ben Bolker does it in netassoc. So my solution (so far)
# is to count a neighbor interaction as a 'site', so our columns will
# be as many as we have potential interactions. Shorter transects will
# have 0s for the final columns. Hopefully this will help the partial
# correlation.

# Our 'sites' here will be the boundary between two intercepts.

intercept_cooccur1<-function(huon, species, perms){
  
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
    mat[, colSums(mat)>0]
  })
  
  #obs_spint_mat <- do.call("cbind", obs_spint_mat_list)
  obs_spint_mat.z <- do.call("cbind", obs_spint_mat_list.z)
  
  #obs_pcorr<-partial_correlation(obs_spint_mat, method="shrinkage")
  obs_pcorr<-partial_correlation(obs_spint_mat.z, method="shrinkage")
  
  # 2. Calculate expected interactions ####
  
  # 2.2 create array of random communities  ####

  int_sp_permlist <- lapply(obs_spint_mat_list, function(mat){
  
    permatswap(mat, times=perms)  
    
    })
  
  # 2.3 create a single species-boundary matrix for each permutation ####
  perm.mats <- lapply(1:perms, function(n){
    
    tempmat <- do.call("cbind", lapply(int_sp_permlist, function(permlist){
      permlist$perm[[n]]
    }))
    
    tempmat[, colSums(tempmat)>0]
    
  })
  
    # 2.4 Calculate partial correlations for each permutation. ####
    
    # We need to have a look here at errors, potentially from trying to 
    # estimate the partial correlation of species that DO NOT interact in the
    # permutation (or too many perhaps). We could try trimming this down a
    # spot.
    spsp_perm_pcorr<-array(
      
      sapply(1:length(perm.mats), function(perm){
        
        print(perm)
        perm<-perm.mats[[perm]]
      
      partial_correlation(perm, method="shrinkage")
        
        }),
      
      dim=c(length(sp.levels), length(sp.levels), perms),
      dimnames=list(sp.levels, sp.levels, 1:perms))
    
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
