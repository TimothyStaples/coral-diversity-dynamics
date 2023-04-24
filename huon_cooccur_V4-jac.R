## ##################################################### ####
# Coexistence of species in Huon Holocence transect data ####
# Author:    Timothy L Staples                           ####
# ###################################################### ####
# WORKING DIRECTORY ####

rm(list = ls())
opSys <- Sys.info()["sysname"]
setwd(paste0(ifelse(opSys == "Linux",
                    "/home/timothy/",
                    "/Users/uqtstapl/"),
             "Dropbox/Tim/Post-doc/Research projects/huon_cooccur"))


# FUNCTIONS ####

sapply(list.files(path = "./Functions", pattern = ".R", full.names = TRUE),
       source)

date.wrap <- function(string, ext){
  paste0(string, " ", Sys.Date(), ext)
}

rotateCoords <- function(crds, angle=0, center= c(min(crds[,1]),min(crds[,2]))) {
  co <- cos(-angle*pi/180)
  si <- sin(-angle*pi/180)
  adj <- matrix(rep(center,nrow(crds)),ncol=2,byrow=TRUE)
  crds <- crds-adj
  cbind(co * crds[,1] - si * crds[,2],
        si * crds[,1] + co * crds[,2]) + adj
}

`%near%` <- function(a,b){
  tolerance = .Machine$double.eps^0.5
  abs(a - b) < tolerance
}

test.angles <- data.frame(x=c(1,0,-1,0),
                          y=c(0,1,0,-1))
test.angles$angle <- cont.angles(test.angles$x, test.angles$y)

cont.angles <- function(x,y){
  # converts differenced vector into angles on a continuous 360 degree scale
  angles <- -atan2(x, y) * (180/pi)
  angles <- ifelse(angles > 0, angles, 360 - abs(angles))
  rot.angles <- angles + 90
  rot.angles <- ifelse(rot.angles >= 360, rot.angles - 360, rot.angles)
  return(rot.angles)
}

cut.centers <- function(x, type){
  
  if(type=="center"){
    return(rowMeans(cbind(as.numeric( sub("\\((.+),.*", "\\1", x) ),
                          as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", x) ))))
  }
  
  if(type=="lower"){
    return(as.numeric( sub("\\((.+),.*", "\\1", x) ))
  }
  
  if(type=="upper"){
    return(as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", x) ))
  }
}


# PACKAGES ####

package.loader(c("vegan", "lme4", "MASS", "DHARMa", "shape", "betapart",
                 "circular", "maptools", "sp", "rgdal", "rgeos", "fields",
                 "grImport2", "mgcv", "plotrix", "iNEXT", "performance",
                 "TeachingDemos"))

# ####
# DATA IMPORT ####
#        Huon peninsula fossil transects ####

huon_tr <- read.csv("./raw.datafiles/huon_intercept_data.csv", stringsAsFactors = TRUE)
huon_site <- read.csv("./raw.datafiles/huon_transect_data.csv", stringsAsFactors = TRUE)

colnames(huon_site)[match(c("age.mean", "age.max", "age.min"),
                          colnames(huon_site))] = c("transect.age", "date.udi", "date.ldi")

huon_tr <- merge(huon_tr, huon_site[, c("transect", "site", "locality",
                                     "transect.age", "date.ldi", "date.udi")],
               by.x = "transect", by.y = "transect")

# Just want hermatypic coral
huon_tr$species.fact[huon_tr$genus == "Sinularia" & !is.na(huon_tr$genus)] = "aherm.coral"

huon_coral<-droplevels(huon_tr[huon_tr$species.fact=="herm.coral" &
                                 !is.na(huon_tr$species.fact),])

# get growth form counts for each tape in each transect
huon_coral$tr.tape<-paste0(huon_coral$transect, "_", huon_coral$tape)

# subset only sites with reasonable temporal resolution
loc.rep<-table(huon_site$locality[!duplicated(huon_site$transect)])
locs<-names(loc.rep[loc.rep>=5])
locs<-locs[locs != "Sang River NW"]
locs.order <- c(2,3,4,1,8,9,5,7,6)

huon_coral<-droplevels(huon_coral[huon_coral$locality %in% locs,])

huon_coral$genus[huon_coral$genus == "Acropora/Isopora"] = "Acropora"
huon_coral$genus[huon_coral$genus == "Gonipora"] = "Goniopora"
huon_coral <- huon_coral[!is.na(huon_coral$locality),]
huon_coral <- droplevels(huon_coral)

locality.number <- huon_coral[!duplicated(huon_coral$locality),c("locality","site")]
locality.number <- locality.number[order(locality.number$locality),]
locality.number$number <- c(6,9,8,7,3,1,2,5,4)

# transect-level dataframe
huon_site <- huon_coral[!duplicated(huon_coral$transect),
                        c("transect", "site", "locality", "transect.age")]
  
#                 Growth form assessment ####

# growth form tables ####

# there are a number of field/taxonomic growth form assessments that
# resulted in multiple Veron growth forms. We need to decide how to
# assess these.

# Step 1 - get the synonymized growth form assessments.
huon_coral$genSp <- paste0(huon_coral$genus, " ", huon_coral$species)
gensptable <- as.matrix(table(huon_coral$growth.comb, huon_coral$genSp))
gensptable <- gensptable[order(rowSums(gensptable),decreasing=TRUE),]

genspDf <- do.call("rbind", lapply(1:nrow(gensptable), function(n1){
  print(n1)
  x <- gensptable[n1,]
  x <- sort(x[x>0], decreasing=TRUE)
  x <- x[x/sum(x) > 0.01]
  xCat <- paste0(sapply(1:length(x), function(n){paste0(names(x)[n], " (", round((x[n] / sum(x))*100, 1),"%)")}), collapse=": ")
  
  return(data.frame(grCat = rownames(gensptable)[n1], grSum = sum(x), spCat = xCat))
}))

gentable <- as.matrix(table(huon_coral$growth.comb, huon_coral$genus))
gentable <- gentable[order(rowSums(gentable),decreasing=TRUE),]
gentable <- gentable[match(rownames(gensptable), rownames(gentable)),]

rownames(gensptable) == rownames(gentable)

genDf <- do.call("rbind", lapply(1:nrow(gentable), function(n1){
  print(n1)
  x <- gentable[n1,]
  x <- sort(x[x>0], decreasing=TRUE)
  x <- x[x/sum(x) > 0.01]
  xCat <- paste0(sapply(1:length(x), function(n){paste0(names(x)[n], " (", round((x[n] / sum(x))*100, 1),"%)")}), collapse=": ")
  
  return(data.frame(grCat = rownames(gentable)[n1], grSum = sum(x), spCat = xCat))
}))

combDf <- cbind(genDf, genspDf$spCat)
write.csv(combDf, "./outputs/combinedGrowthSummary.csv")


fullTable <- as.data.frame(table(huon_coral$growth.comb, huon_coral$genSp))
fullTable <- fullTable[fullTable$Freq > 0,]
fullTable <- fullTable[order(fullTable$Freq, decreasing=TRUE),]
colnames(fullTable) <- c("growth form", "species", "count")
write.csv(fullTable, "./outputs/fullGrowthSummary.csv")

# assign growth forms ####

# start by getting the most common terms by genus

# growth names by genus
huon_coral$growth.comb <- gsub("\\?", "unknown", as.character(huon_coral$growth.comb))
gr.simp <- gsub('[[:punct:] ]+', " ", as.character(huon_coral$growth.comb))
gr.words <- strsplit(gr.simp, " ")
gr.simp <- sapply(strsplit(gr.simp, " "),  function(x){paste(unique(x), collapse=" ")})

gr.words <- do.call('rbind', lapply(1:length(gr.words), function(n){
  
  if(length(gr.words[[n]]) == 0){return(NULL)}
  cbind(as.character(huon_coral$genus[n]),
        unique(as.vector(gr.words[[n]])))
}))

# table genus and words
gr.word.tab <- table(gr.words[,1], gr.words[,2])
gr.word.tab[1:10,1:10]

apply(gr.word.tab, 1, function(x){
  sort(x[x>0], decreasing=TRUE)
}, simplify=FALSE)

apply(gr.word.tab, 2, function(x){
  sort(x[x>0], decreasing=TRUE)
}, simplify=FALSE)


# create assignment for each growth form
gr.cats <- read.csv("./raw.datafiles/growth.grep.new.csv")

growth.mat <- sapply(1:nrow(gr.cats), function(n){
  
  grPos <- grepl(gr.cats$pos.arg[n], tolower(gr.simp))
  
  if(is.na(gr.cats$neg.arg[n])){
  return(grPos)
  } else {
  return((grPos + !grepl(gr.cats$neg.arg[n], tolower(gr.simp))) == 2)
  }
})
colnames(growth.mat) = gr.cats$growth.form

# Step 2 - if multiple, collapse to most likely of the two for the intercept genus
#          (E.g., if Acropora intercept listed as "massive" and "branching_open", record
#           as "branching_open")

# if the other potential growth from is branching, favour the other, indicating
# a difference that contributes to functional diversity
branchingCols <- grepl("branching", colnames(growth.mat))

# remove unknowns, solitary if there's a non-unknown ID (probably unknown fallen into the
# "acropora" branching id column)
growth.mat[rowSums(growth.mat[,colnames(growth.mat) != "unknown"]) > 0,"unknown"] = 0

singles <- which((growth.mat>0)==1 & rowSums(growth.mat>0)==1, arr.ind = TRUE)

length(unique(singles[,1]))/nrow(growth.mat)

growth.form <- rep("unknown", nrow(growth.mat))
growth.form[singles[,1]] = colnames(growth.mat)[singles[,2]]
growth.form[huon_coral$growth.comb %in% c("", "base")] = "unknown"
growth.form[huon_coral$growth.comb==""] = "unknown"

growth.form[growth.form == "unknown" &
            huon_coral$genSp %in% c("Acropora hyacinthus", "Acropora solitaryensis",
                                    "Acropora cytherea",
                                    "Acropora divaricata")] = "table_plate_laminar"

growth.form[growth.form == "unknown" &
              huon_coral$genSp %in% c("Acropora digitifera", "Acropora monticulosa",
                                      "Acropora cerealis", "Acropora nana",
                                      "Acropora echinata", "Acropora nasuta",
                                      "Acropora florida", "Acropora pulchra",
                                      "Acropora gemmifera", "Acropora robusta",
                                      "Acropora intermedia", "Acropora selago",
                                      "Acropora lutkeni", "Acropora valida",
                                      "Acropora milleopora", "Isopora cuneata",
                                      "Pocillopora damicornis",
                                      "Stylophora pistillata", "Acropora abrotanoides",
                                      "Acropora acuminata", "Acropora humilis",
                                      "Acropora papillare", "Acropora valenciennesi",
                                      "Pocillopora grandis", "Pocillopora meandrina",
                                      "Pocillopora verrucosa", "Porites annae",
                                      "Porites cylindrica", "Porites nigrescens")] = "branching"


growth.form[growth.form == "unknown" &
              huon_coral$genSp %in% c("Isopora palifera")] = "branching_stout" 

growth.form[growth.form == "unknown" & 
              huon_coral$genSp %in% c("Acropora minuta",
                                      "Galaxea astreata",
                                      "Hydnophora exesa",
                                      "Pavona explanulata", "Porites lichen",
                                      "Porites vaughani", "Stylocoeniella armata")] = "encrusting"


growth.form[growth.form == "unknown" & 
              huon_coral$genSp %in% c("Galaxea fascicularis", "Montipora informis",
                                      "Fungia fungites", "Cyphastrea microphthalma",
                                      "Cyphastrea serailia", "Favites flexuosa",
                                      "Favites halicora", "Goniastrea edwardsi",
                                      "Goniastrea stelligera", "Hydnophora microconos",
                                      "Lobophyllia corymbosa", "Lobophyllia hemprichii",
                                      "Montipora corbettensis", "Montipora grisea",
                                      "Montipora turgescens", "Pleuractis paumotensis")] = "massive_columnar_submassive"

## override other assessments of Galaxea
growth.form[huon_coral$genSp %in% c("Galaxea fascicularis", "Galaxea NA")] = "massive_columnar_submassive"


# add growth form column and remove unknown column
huon_coral$growth.form <- as.factor(growth.form)
huon_coral_all <- huon_coral
huon_coral <- droplevels(huon_coral[huon_coral$growth.form != "unknown",])

# TIME-SPACE CORRELATION ####

# Let's see how much of our time-space we actually have filled, and whether
# there is any confounding of these effects (if some sections of the study region
# were all from the earlier part of the overall time-series, for example)

huon_loc <- huon_site[!duplicated(huon_site$locality),]

rownames(huon_loc) = huon_loc$locality

huon_loc$locality <- as.character(huon_loc$locality)
huon_loc$distance <- sapply(as.character(huon_loc$locality), function(x){

  gc.dist(lat1= huon_loc["Midway Cove", "latitude"],
          lon1 = huon_loc["Midway Cove", "longitude"],
          lat2= huon_loc[x, "latitude"],
          lon2= huon_loc[x, "longitude"])  
  
})

huon_subsite <- droplevels(merge(huon_site, huon_loc[,c("distance","locality")],
                   all.x=TRUE, all.y=TRUE))
huon_subsite <- droplevels(huon_subsite[huon_subsite$locality %in% huon_coral$locality,])

pdf(date.wrap("./plots/sample scale", ".pdf"), height=3.5, width=5)
par(mar=c(3,3,1,1), ps=8, tcl = -0.25, mgp=c(3,0.5,0), las=1)

plot(huon_subsite$distance ~ huon_subsite$transect.age, pch=21,
     bg=c("blue","red","darkgreen")[huon_subsite$site], xaxt="n",
     xlim=c(9500,6000), ylim=c(0,30), type="n", xlab="", ylab="")

sapply(split(huon_subsite, f=huon_subsite$locality), function(x){
 
  segments(x0=min(x$transect.age),
           x1=max(x$transect.age),
           y0=x$distance[1],
           y1=x$distance[1]) 
  
})
points(huon_subsite$distance ~ huon_subsite$transect.age, pch=21,
       bg=c("blue","red","darkgreen")[huon_subsite$site])

axis(side=1, mgp=c(3,0.2,0))
mtext(side=1, line=1.25, text="Years before present")
mtext(side=2, line=1.5, text="Distance from southern-most locality (km)", las=0)

dev.off()

# DATA PREP ####
#                 Compositional Matrices ####

genusMatCount <- table(huon_coral$transect, droplevels(huon_coral$genus))
grMatCount <- table(huon_coral$transect, droplevels(huon_coral$growth.form))

genusMatProp <- do.call("rbind", 
                          lapply(split(huon_coral, f=huon_coral$transect), 
                                 function(x){
                                   
                                   print(x$transect[1])
                                   x <- x[x$intercept != "unknown", ]
                                   
                                   x$sqrt.int <- sqrt(x$int.delta)
                                   temp <- tapply(x$sqrt.int, x$genus, sum, na.rm=TRUE) / sum(x$sqrt.int[!is.na(x$genus)], na.rm=TRUE)
                                   ifelse(is.na(temp), 0, temp)
                                 }))

grMatProp <- do.call("rbind", 
                       lapply(split(huon_coral, f=huon_coral$transect), 
                              function(x){
                                
                                x <- x[x$intercept != "unknown", ]
                                
                                x$sqrt.int <- sqrt(x$int.delta)
                                temp <- tapply(x$sqrt.int, x$growth.form, sum, na.rm=TRUE) / 
                                  sum(x$sqrt.int[!is.na(x$growth.form)], na.rm=TRUE)
                                ifelse(is.na(temp), 0, temp)
                              }))

#                        Alpha diversity ####

richDf <- alphaDiversityCalc(siteData = huon_site,
                             genusMatProp = genusMatProp,
                             grMatProp = grMatProp)

#                         Beta diversity ####

turnoverDf <- betaDiversityCalc(siteData = huon_site,
                                genusMatProp = genusMatProp,
                                grMatProp = grMatProp,
                                method="jaccard")

# MODELS ####
#                Multivariate GAM models ####

# alpha.model.list <- mvGams(dataList = list(richDf, richDf, richDf),
#                            primaryvarList = c("coordH0", "coordH1", "coordH2"),
#                            secondaryvarList = c("ratioH0", "ratioH1", "ratioH2"),
#                            beta=FALSE, modelIter = 10000, modelWarmup=1000,
#                            modelPath = "./Outputs")
# 
# beta.model.list <- mvGams(dataList = list(turnoverDf, turnoverDf),
#                            primaryvarList = c("dissimilarity", "turnover"),
#                            secondaryvarList = c("conservation", "persistance"),
#                           beta=TRUE, modelIter = 10000, modelWarmup=1000,
#                           modelPath = "./Outputs")

alpha.model.list <- mvGamsFreq(dataList = list(richDf, richDf, richDf),
                           primaryvarList = c("genusH0", "genusH1", "genusH2"),
                           secondaryvarList = c("grH0", "grH1", "grH2"),
                           dateDiff=FALSE)

beta.model.list <- mvGamsFreq(dataList = list(turnoverDf, turnoverDf),
                          primaryvarList = c("top.gen.diss", "seq.gen.diss"),
                          secondaryvarList = c("top.gr.diss", "seq.gr.diss"),
                          dateDiff=TRUE)

# list structure is: (1) diversity type, (2) models & predictions, (3) spatial scale
sapply(1:3, function(n){
  
  p.table <- summary(alpha.model.list[[n]]$models[[1]])$p.table
  s.table <- summary(alpha.model.list[[n]]$models[[1]])$s.table
  
  write.csv(p.table, date.wrap(paste0("./outputs/Region-Ptable-H",n-1), ".csv"))
  write.csv(s.table, date.wrap(paste0("./outputs/Region-Stable-H",n-1), ".csv"))
  
  p.table <- summary(alpha.model.list[[n]]$models[[3]])$p.table
  s.table <- summary(alpha.model.list[[n]]$models[[3]])$s.table
  
  write.csv(p.table, date.wrap(paste0("./outputs/Site-Ptable-H",n-1), ".csv"))
  write.csv(s.table, date.wrap(paste0("./outputs/Site-Stable-H",n-1), ".csv"))

})
sapply(1:2, function(n){
  
  p.table <- summary(beta.model.list[[n]]$models[[1]])$p.table
  s.table <- summary(beta.model.list[[n]]$models[[1]])$s.table
  
  write.csv(p.table, date.wrap(paste0("./outputs/Region-Ptable-",c("Dir","Seq")[n]), ".csv"))
  write.csv(s.table, date.wrap(paste0("./outputs/Region-Stable-",c("Dir","Seq")[n]), ".csv"))
  
  p.table <- summary(beta.model.list[[n]]$models[[3]])$p.table
  s.table <- summary(beta.model.list[[n]]$models[[3]])$s.table
  
  write.csv(p.table, date.wrap(paste0("./outputs/Site-Ptable-",c("Dir","Seq")[n]), ".csv"))
  write.csv(s.table, date.wrap(paste0("./outputs/Site-Stable-",c("Dir","Seq")[n]), ".csv"))
  
})

compare_performance(alpha.model.list[[1]]$models[[1]],
                    alpha.model.list[[2]]$models[[1]],
                    alpha.model.list[[3]]$models[[1]])
compare_performance(alpha.model.list[[1]]$models[[3]],
                    alpha.model.list[[2]]$models[[3]],
                    alpha.model.list[[3]]$models[[3]])

compare_performance(beta.model.list[[1]]$models[[1]],
                    beta.model.list[[2]]$models[[1]])
compare_performance(beta.model.list[[1]]$models[[3]],
                    beta.model.list[[2]]$models[[3]])


#                      Diversity vectors ####

vectorList <- diversityVectors(rich.df = richDf,
                               turnover.df = turnoverDf)

# NULL MODELS ####
#                      Markov Simulation ####    

simIter = 99

# generate simulated turnover data for each site
markovSimList <- lapply(unique(huon_site$locality), function(loc){
  
  print(as.character(loc))
  
  loc_sub <- droplevels(huon_site[huon_site$locality == loc,])
  
  genData <- genusMatProp[as.character(loc_sub$transect),]
  grData <- grMatProp[as.character(loc_sub$transect),]
  
  tDates <- sort(loc_sub$transect.age, decreasing=TRUE)
  
  # Simulations
  genSims <- markovSimm(nrow(genData), genData, tDates, simIter)
  grSims <- markovSimm(nrow(grData), grData, tDates, simIter)
  
  # Pair and calculate beta diversities
  simAlpha <- lapply(1:simIter, function(n){
    
    tempGen <- genSims[[n]]
    rownames(tempGen) = loc_sub$transect
    tempGr <- grSims[[n]]
    rownames(tempGr) = loc_sub$transect
    
    alphaDiversityCalc(loc_sub, tempGen, tempGr)
  })
  simBeta <- lapply(1:simIter, function(n){
    print(n)
    if(sum(is.na(genSims[[n]])) > 0 | sum(is.na(grSims[[n]]))){
      return(NULL)
    }
    betaDiversityCalc(loc_sub, genSims[[n]], grSims[[n]])
  })
  simBeta <- simBeta[!sapply(simBeta, is.null)]
  
  return(list(simAlpha, simBeta, genSims, grSims))
  
})

markovAlphaModels <- lapply(1:simIter, function(n){
  print(n)
  tempData <- do.call("rbind", lapply(markovSimList, 
                                      function(x){x[[1]][[n]]}))
  
  # make sure the simulation didn't break (NA values)
  if(sum(is.na(tempData$genusH0)) + sum(is.na(tempData$grH0)) > 0){
    print("Simulation fail, no model")
    return(NULL)
  }
  
  a <- mvGamsFreq(dataList = list(tempData, tempData, tempData),
             primaryvarList = c("genusH0", "genusH1", "genusH2"),
             secondaryvarList = c("grH0", "grH1", "grH2"),
             dateDiff=FALSE)
  
return(list(a[[1]]$preds, a[[2]]$preds, a[[3]]$preds))
  
  })
markovAlphaModels <- markovAlphaModels[!sapply(markovAlphaModels, is.null)]

markovBetaModels <- lapply(1:simIter, function(n){
  print(n)
  
  simLengths <- sapply(markovSimList, function(x){length(x[[2]])})
  
  if(sum(n > simLengths) > 0){return(NULL)}
  tempData <- do.call("rbind", lapply(markovSimList, 
                                      function(x){x[[2]][[n]]}))
  
  a <- mvGamsFreq(dataList = list(tempData, tempData),
             primaryvarList = c("top.gen.diss", "seq.gen.diss"),
             secondaryvarList = c("top.gr.diss", "seq.gr.diss"),
             dateDiff = TRUE)
  
  return(list(a[[1]]$preds, a[[2]]$preds))
})

#                     Matrix permutation ####

# generate simulated turnover data for each site
permatSimList <- lapply(unique(huon_site$locality), function(loc){
  
  print(as.character(loc))
  
  loc_sub <- droplevels(huon_site[huon_site$locality == loc,])
  
  genData <- genusMatProp[as.character(loc_sub$transect),]
  grData <- grMatProp[as.character(loc_sub$transect),]
  
  #swaps can only occur on individuals, so we treat each 0.01% of sampling
  #as an individual in a count simulation
  genDat <- round(genData*1000)
  grDat <- round(grData*1000)
  
  # Simulations
  genSims <- permatswap(genDat, fixedmar="both", times=simIter, 
                        method="quasiswap")$perm
  grSims <- permatswap(grDat, fixedmar="both", times=simIter, 
                       method="quasiswap")$perm
  
  # Pair and calculate beta diversities
  simAlpha <- lapply(1:simIter, function(n){
    
    tempGen <- prop.table(genSims[[n]], 1)
    rownames(tempGen) = loc_sub$transect
    tempGr <- prop.table(grSims[[n]], 1)
    rownames(tempGr) = loc_sub$transect
    
    alphaDiversityCalc(loc_sub, tempGen, tempGr)
  })
  
  # Pair and calculate beta diversities
  simBeta <- lapply(1:simIter, function(n){
    betaDiversityCalc(loc_sub, 
                      prop.table(genSims[[n]], 1), 
                      prop.table(grSims[[n]],1))
  })
  
  return(list(simAlpha, simBeta))
  
})

# now use each set of data to run turnover models

permatAlphaModels <-  lapply(1:simIter, function(n){
  print(n)
  tempData <- do.call("rbind", lapply(permatSimList, 
                                      function(x){x[[1]][[n]]}))
  
  a <- mvGamsFreq(dataList = list(tempData, tempData, tempData),
             primaryvarList = c("genusH0", "genusH1", "genusH2"),
             secondaryvarList = c("grH0", "grH1", "grH2"))
  
  return(list(a[[1]]$preds, a[[2]]$preds, a[[3]]$preds))
  
})
permatBetaModels <- lapply(1:simIter, function(n){
  print(n)
  tempData <- do.call("rbind", lapply(permatSimList, function(x){x[[2]][[n]]}))
  
  a <- mvGamsFreq(dataList = list(tempData, tempData),
             primaryvarList = c("top.gen.diss", "seq.gen.diss"),
             secondaryvarList = c("top.gr.diss", "seq.gr.diss"))
  
  return(list(a[[1]]$preds, a[[2]]$preds))
  
})

# save models
saveRDS(markovSimList, "./outputs/markovSimList.RDS")
saveRDS(permatSimList, "./outputs/permatSimList.RDS")

saveRDS(markovAlphaModels, "./outputs/markovAlphaModels.RDS")
saveRDS(markovBetaModels, "./outputs/markovBetaModels.RDS")

saveRDS(permatAlphaModels, "./outputs/permatAlphaModels.RDS")
saveRDS(permatBetaModels, "./outputs/permatBetaModels.RDS")

markovSimList <- readRDS("./outputs/markovSimList.RDS")
permatSimList <- readRDS("./outputs/permatSimList.RDS")
markovAlphaModels <- readRDS("./outputs/markovAlphaModels.RDS")
markovBetaModels <- readRDS("./outputs/markovBetaModels.RDS")
permatAlphaModels <- readRDS("./outputs/permatAlphaModels.RDS")
permatBetaModels <- readRDS("./outputs/permatBetaModels.RDS")

# null region plots ####
plot(x=NULL, y=NULL, xlim=c(0,20), ylim=c(0,20))

markovTaxfit <- do.call("rbind", lapply(markovAlphaModels, function(x){
  x[[1]][[1]]$fit.1
  }))
markovTaxMean <- colMeans(markovTaxfit)
markovFunfit <- do.call("rbind", lapply(markovAlphaModels, function(x){
  x[[1]][[1]]$fit.2
}))
markovFunMean <- colMeans(markovFunfit)

lapply(markovAlphaModels, function(x){
  lines(x[[1]][[1]]$fit.2 ~   x[[1]][[1]]$fit.1, col=rgb(0.8,0.8,0.8,0.3))
})
lines(markovFunMean ~ markovTaxMean)

permTaxfit <- do.call("rbind", lapply(permatAlphaModels, function(x){
  x[[1]][[1]]$fit.1
}))
permTaxMean <- colMeans(permTaxfit)
permFunfit <- do.call("rbind", lapply(permatAlphaModels, function(x){
  x[[1]][[1]]$fit.2
}))
permFunMean <- colMeans(permFunfit)

lapply(permatAlphaModels, function(x){
  lines(x[[1]][[1]]$fit.2 ~   x[[1]][[1]]$fit.1, col=rgb(1,0.1,0.1,0.3))
})
lines(permFunMean ~ permTaxMean)


# PLOTS ####

#   Dir vs seq beta diversity comparison ####

betaExamplePlot("./Plots/turnover ords.pdf")

#                               Site map ####

library(raster)
siteMapPlot("./Plots/huon map rotated.pdf")

# region level trends ####

#                                  alpha ####
# funtion to plot trend arrow in diversity space as well as segments for age
main.trend <- function(preds, col){
  
  fit1 <- preds$fit.1
  fit2 <- preds$fit.2
  pred.date <- preds$pred.date
  
  lines(fit2 ~ fit1, col=col)  
  Arrows(x0=fit1[3], x1=fit1[1],
         y0=fit2[3], y1=fit2[1],
         arr.type="triangle", arr.length=0.1,
         arr.width=0.1, col=col)
  
  sapply(max(pred.date) - seq(0,5000,250), function(seg.n){
    
    temp.match <- which(pred.date == seg.n)
    base.age <- max(pred.date)
    diff.age <- base.age - seg.n
    
    if(length(temp.match)>0){
      
      target.fit1 <- fit1[temp.match - c(1,0)]
      target.fit2 <- fit2[temp.match - c(1,0)]
      
      if(length(target.fit1)<2 | length(target.fit2) <2){return(NULL)}
      arrows(x0=target.fit1[1], x1=target.fit1[2],
             y0=target.fit2[1], y1=target.fit2[2],
             angle=145, col=col,
             length=ifelse(diff.age %% 1000 == 0, 0.025,0.025),
             lwd=ifelse(diff.age %% 1000 == 0, 0.5,0.5))
      
    }
    
  })
  
}

raw.plot <- function(ylims){
  par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0,0), las=1)
  plot(x=NULL, y=NULL, xlim=c(9400,5900), ylim=ylims,
       yaxs="i", xlab="", ylab="", axes=FALSE)
  axis(side=1, at=seq(0,10000,500), labels=NA, tcl=-0.125)
}

raw.trend <- function(preds, fit.var, line.col, poly.col){
  
  preds$trend <- preds[,fit.var]
  preds$se <- preds[, paste0("se.fit.",substr(fit.var, 
                                              nchar(fit.var), 
                                              nchar(fit.var)))]
  preds$upper <- preds$trend + 1.96 * preds$se
  preds$lower <- preds$trend - 1.96 * preds$se
  
  polygon(x=c(preds$pred.date, rev(preds$pred.date)),
          y=c(preds$upper, rev(preds$lower)),
          border=NA, col=poly.col)
  lines(preds$trend ~ preds$pred.date, lwd=1.5, col=line.col)
  Arrows(x0=preds$pred.date[2], x1=preds$pred.date[1],
         y0=preds$trend[2], y1=preds$trend[1], col=line.col,
         arr.type="triangle", arr.length=0.1, arr.width=0.1)
  
  sapply(max(preds$pred.date) - seq(0,5000,250), function(seg.n){
    
    temp.match <- which(preds$pred.date == seg.n)
    base.age <- max(preds$pred.date)
    diff.age <- base.age - seg.n
    
    if(length(temp.match)>0){
      
      target <- preds[temp.match - c(1,0),]
      if(nrow(target)<2){return(NULL)}
      arrows(x0=target$pred.date[1], x1=target$pred.date[2],
             y0=target$trend[2], y1=target$trend[2],
             angle=90, col=line.col,
             length=ifelse(diff.age %% 1000 == 0, 0.075,0.025),
             lwd=ifelse(diff.age %% 1000 == 0, 1,0.5))
      
    }
    
  })
  
}

pdf(date.wrap("./Plots/alpha region trends",".pdf"), 
    height=8.5, width=5.75, useDingbats = FALSE)

height <- ((0.97-0.06)-0.1) / 3
bots <- c(0.06) + height * c(0:2) + c(0, 0.05, 0.1)
tops <- c(0.06) + height * c(1:3)

split.screen(rbind(c(0.08,0.48,0.72,0.99),
                   c(0.08,0.48,0.39,0.66),
                   c(0.08,0.48,0.06,0.33),
                   
                   c(0.59,0.99,0.855,0.99),
                   c(0.59,0.99,0.72,0.855),
                   c(0.59,0.99,0.525,0.66),
                   c(0.59,0.99,0.39,0.525),
                   c(0.59,0.99,0.195,0.33),
                   c(0.59,0.99,0.06,0.195)))

loc.color = loc.color
model.list = c(alpha.model.list, alpha.model.list, alpha.model.list)
xlabs=list(expression("Tax. "*alpha*" diversity"),
           expression("Tax. "*alpha*" diversity"),
           expression("Tax. "*alpha*" diversity"))
ylabs=list(expression("Fun. "*alpha*" diversity"),
           expression("Fun. "*alpha*" diversity"),
           expression("Fun. "*alpha*" diversity"))
xlims=list(c(0,13),c(0,6),c(0,5))
ylims=list(c(0,13),c(0,6),c(0,5))

sapply(1:3, function(n){
                        
screen(n)
print(n)
x.lims<-xlims[[n]]
y.lims<-ylims[[n]]
x.lab<-xlabs[[n]]
y.lab<-ylabs[[n]]

if(n <= 3){
temp.preds <- alpha.model.list[[n]]$preds[[1]]
} else {
temp.preds <- beta.model.list[[n-3]]$preds[[1]]
}  

par(mar=c(0,0,0,0), ps=10, las=1, tcl=-0.25, mgp=c(3,0.5,0))

plot(x=NULL, y=NULL, xlim=x.lims, ylim=y.lims, asp=1,
     xlab="",ylab="", yaxs="i", xaxs="i", axes=FALSE)
axis(side=1, mgp=c(3,0.2,0))
axis(side=2, las=1)
mtext(side=2, line=1.25, text=y.lab, las=0)
mtext(side=1, line=1.25, text=x.lab)

abline(a=0, b=1, lty="31", col="grey80")

main.trend(temp.preds, "black")
text(x=relative.axis.point(0.03, "x"),
     y=relative.axis.point(0.935,"y"),
     labels=paste0("(",LETTERS[c(1,4,7)][n],")"), adj=0, font=2)
text(x=relative.axis.point(0.125, "x"),
     y=relative.axis.point(0.935,"y"),
     labels=paste0("Hill q = ", n-1), adj=0)

box()
close.screen(n)

screen(n+(n-1)+3)

fit1.ylims <- c(min(temp.preds$fit.1 - 1.96 * temp.preds$se.fit.1),
                max(temp.preds$fit.1 + 1.96 * temp.preds$se.fit.1)) * c(0.9,1.1)

raw.plot(fit1.ylims)
axis(side=1, at=seq(6000,9000,1000), labels=NA)
axis(side=2, mgp=c(3,0.5,0))
mtext(side=2, line=1.25, text=x.lab, las=0)
raw.trend(preds=temp.preds,
          fit.var="fit.1",line.col="black", poly.col="grey80")
text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.9, "y"),
     labels=paste0("(",LETTERS[c(2,5,8)][n],")"), font=2, adj=0)
box()
close.screen(n+(n-1)+3)

screen(n+(n-1)+4)

fit1.ylims <- c(min(temp.preds$fit.2 - 1.96 * temp.preds$se.fit.2),
                max(temp.preds$fit.2 + 1.96 * temp.preds$se.fit.2)) * c(0.9,1.1)

raw.plot(fit1.ylims)
axis(side=1, at=seq(6000,9000,1000), labels=format(seq(6000,9000,1000), big.mark=","), mgp=c(3,0.2,0))
mtext(side=1, line=1, text="Time (years before present)")
axis(side=2, mgp=c(3,0.5,0))
mtext(side=2, line=1.75, text=y.lab, las=0)
raw.trend(preds=temp.preds,
          fit.var="fit.2",line.col="black", poly.col="grey80")
text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.9, "y"),
     labels=paste0("(",LETTERS[c(3,6,9)][n],")"), font=2, adj=0)
box()
close.screen(n+(n-1)+4)

})
close.screen(all.screens=TRUE)
dev.off()

#                                   beta ####

pdf(date.wrap("./Plots/beta region trends jaccard",".pdf"), 
    height=8.5, width=5.75, useDingbats = FALSE, bg="white")

split.screen(rbind(c(0.08,0.48,0.72,0.99),
                   c(0.08,0.48,0.39,0.66),
                   c(0.08,0.48,0.06,0.33),
                   
                   c(0.59,0.99,0.855,0.99),
                   c(0.59,0.99,0.72,0.855),
                   c(0.59,0.99,0.525,0.66),
                   c(0.59,0.99,0.39,0.525),
                   c(0.59,0.99,0.195,0.33),
                   c(0.59,0.99,0.06,0.195)))

loc.color = loc.color
model.list = c(beta.model.list, beta.model.list)
xlabs=list(expression("Tax. "*beta*" diversity"),
           expression("Tax. "*beta*" diversity"))
ylabs=list(expression("Fun. "*beta*" diversity"),
           expression("Fun. "*beta*" diversity"))
xlims=list(c(0.2,0.7), c(0.2,0.7))
ylims=list(c(0.2,0.7), c(0.2,0.7))

sapply(1:2, function(n){
  
  screen(n)
  print(n)
  x.lims<-xlims[[n]]
  y.lims<-ylims[[n]]
  x.lab<-xlabs[[n]]
  y.lab<-ylabs[[n]]
  
  temp.preds <- beta.model.list[[n]]$preds[[1]]

  
  par(mar=c(0,0,0,0), ps=10, las=1, tcl=-0.25, mgp=c(3,0.5,0))
  
  plot(x=NULL, y=NULL, xlim=x.lims, ylim=y.lims, asp=1,
       xlab="",ylab="", yaxs="i", xaxs="i", axes=FALSE)
  axis(side=1, mgp=c(3,0.2,0))
  axis(side=2, las=1)
  mtext(side=2, line=1.5, text=y.lab, las=0)
  mtext(side=1, line=1.25, text=x.lab)
  
  abline(a=0, b=1, lty="31", col="grey80")
  
  main.trend(temp.preds, "black")
  text(x=relative.axis.point(0.03, "x"),
       y=relative.axis.point(0.935,"y"),
       labels=paste0("(",LETTERS[c(1,4,7)][n],")"), adj=0, font=2)
  text(x=relative.axis.point(0.125, "x"),
       y=relative.axis.point(0.935,"y"),
       labels=c("Directional (baseline)",
                "Sequential (pairwise)")[n], adj=0)
  
  box()
  close.screen(n)
  
  screen(n+(n-1)+3)
  
  fit1.ylims <- c(min(temp.preds$fit.1 - 1.96 * temp.preds$se.fit.1),
                  max(temp.preds$fit.1 + 1.96 * temp.preds$se.fit.1)) * c(0.9,1.1)
  
  raw.plot(fit1.ylims)
  axis(side=1, at=seq(6000,9000,1000), labels=NA)
  axis(side=2, mgp=c(3,0.5,0))
  mtext(side=2, line=2, text=x.lab, las=0)
  raw.trend(preds=temp.preds,
            fit.var="fit.1",line.col="black", poly.col="grey80")
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.9, "y"),
       labels=paste0("(",LETTERS[c(2,5,8)][n],")"), font=2, adj=0)
  box()
  close.screen(n+(n-1)+3)
  
  screen(n+(n-1)+4)
  
  fit1.ylims <- c(min(temp.preds$fit.2 - 1.96 * temp.preds$se.fit.2),
                  max(temp.preds$fit.2 + 1.96 * temp.preds$se.fit.2)) * c(0.9,1.1)
  
  raw.plot(fit1.ylims)
  axis(side=1, at=seq(6000,9000,1000), labels=format(seq(6000,9000,1000), big.mark=","), mgp=c(3,0.2,0))
  mtext(side=1, line=1, text="Time (years before present)")
  axis(side=2, mgp=c(3,0.5,0))
  mtext(side=2, line=2, text=y.lab, las=0)
  raw.trend(preds=temp.preds,
            fit.var="fit.2",line.col="black", poly.col="grey80")
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.9, "y"),
       labels=paste0("(",LETTERS[c(3,6,9)][n],")"), font=2, adj=0)
  box()
  close.screen(n+(n-1)+4)
  
})
close.screen(all.screens=TRUE)
dev.off()

# site trends ####

loc.color <- data.frame(locality=unique(huon_coral$locality))
loc.color$col <- c("blue","red","darkgreen")[huon_tr$site[match(loc.color$locality,
                                                                huon_tr$locality)]]

xlabs=list(expression("Tax. "*alpha*" diversity"),
           expression("Tax. "*alpha*" diversity"),
           expression("Tax. "*alpha*" diversity"),
           expression("Tax. "*beta*" diversity"),
           expression("Tax. "*beta*" diversity"))
ylabs=list(expression("Fun. "*alpha*" diversity"),
           expression("Fun. "*alpha*" diversity"),
           expression("Fun. "*alpha*" diversity"),
           expression("Fun. "*beta*" diversity"),
           expression("Fun. "*beta*" diversity"))
xlims=list(c(3,16),c(1.5,7),c(1,5.5), c(0,1), c(-0.2,1))
ylims=list(c(3.5,5.5),c(2,4.5),c(1.5,4), c(-0.2,1), c(-0.2,1))

pdf(date.wrap("./Plots/site trends jaccard",".pdf"), height=8.5, width=5.75, bg="white")

split.screen(rbind(c(0.09,0.49,0.72,0.99),
                   c(0.09,0.49,0.39,0.66),
                   c(0.09,0.49,0.06,0.33),
                   
                   c(0.59,0.99,0.72,0.99),
                   c(0.59,0.99,0.39,0.66)))

sapply(1:5, function(n){
  
  print(n)
  screen(n)
  x.lims<-xlims[[n]]
  y.lims<-ylims[[n]]
  x.lab<-xlabs[[n]]
  y.lab<-ylabs[[n]]
  
  if(n <= 3){
    temp.preds <- alpha.model.list[[n]]$preds[[3]]
  } else {
    temp.preds <- beta.model.list[[n-3]]$preds[[3]]
  }  
  
  par(mar=c(0,0,0,0), ps=10, las=1, tcl=-0.25, mgp=c(3,0.5,0))
  
  plot(x=NULL, y=NULL, xlim=x.lims, ylim=y.lims, #asp=1,
       xlab="",ylab="", yaxs="i", xaxs="i", axes=FALSE)
  axis(side=1, mgp=c(3,0.2,0))
  axis(side=2, las=1)
  mtext(side=2, line=1.5, text=y.lab, las=0)
  mtext(side=1, line=1.25, text=x.lab)
  
  abline(a=0, b=1, lty="31", col="grey80")
  
  a <- sapply(levels(loc.color$locality), function(x){
    
    print(x)
    site.range <- range(richDf$pred.date[richDf$locality==x])
    
    localBin <- temp.preds$locality == x &
      temp.preds$pred.date >= site.range[1] &
      temp.preds$pred.date <= site.range[2]
    
    temp.pred <- temp.preds[localBin,]
    
    site.col = loc.color$col[loc.color$locality == x]
    
    main.trend(temp.pred, site.col)
    
    text(x=temp.pred$fit.1[1],
         y=temp.pred$fit.2[1],
         pos=4, offset=0.25,
         labels=locality.number$number[locality.number$locality == x],
         col=site.col)
    
  })
  
  text(x=relative.axis.point(0.03, "x"),
       y=relative.axis.point(0.935,"y"),
       labels=paste0("(",LETTERS[c(1,3,5,2,4)][n],")"), adj=0, font=2)
  
  if(n <= 3){
  text(x=relative.axis.point(0.125, "x"),
       y=relative.axis.point(0.935,"y"),
       labels=paste0("Hill q = ", n-1), adj=0)
  } else {
    text(x=relative.axis.point(0.125, "x"),
         y=relative.axis.point(0.935,"y"),
         labels=c("Directional (baseline)",
                "Sequential (pairwise)")[n-3], adj=0)
  }
  box()
  close.screen(n)
  
})
close.screen(all.screens=TRUE)
dev.off()
                      

#                 Diversity trends plots ####

loc.color <- data.frame(locality=unique(huon_coral$locality))
loc.color$col <- c("blue","red","darkgreen")[huon_tr$site[match(loc.color$locality,
                                                                huon_tr$locality)]]
dev.off()
diversityTrendsPlotsF(loc.color = loc.color,
                     model.list = alpha.model.list,
                     rich.df = richDf,
                     turnover.df = turnoverDf,
                     fullPlotPath = date.wrap("./Plots/multi gam results alpha",".pdf"),
                     xlabs=list(expression("Tax "*alpha*" diversity (Hill q=0)"),
                                expression("Tax "*alpha*" diversity (Hill q=1)"),
                                expression("Tax "*alpha*" diversity (Hill q=2)")),
                     ylabs=list(expression("Fun "*alpha*" diversity (Hill q=0)"),
                                expression("Fun "*alpha*" diversity (Hill q=1)"),
                                expression("Fun "*alpha*" diversity (Hill q=2)")),
                     xlims=list(c(6,14),c(2.2,8),c(1.5,6.7)),
                     ylims=list(c(-1,8),c(-2.2,2.2),c(-1.5,1.5)),
                     
                     vector.list = vectorList[1:3],
                     v.xlabs=list(expression(Delta*" Coordinated "*alpha*" diversity (Hill q=0)"),
                                expression(Delta*" Coordinated "*alpha*" diversity (Hill q=1)"),
                                expression(Delta*" Coordinated "*alpha*" diversity (Hill q=2)")),
                     v.ylabs=list(expression(Delta*" Ratio "*alpha*" diversity (Hill q=0)"),
                                expression(Delta*" Ratio "*alpha*" diversity (Hill q=1)"),
                                expression(Delta*" Ratio "*alpha*" diversity (Hill q=2)")),
                     v.xlims=list(c(-7,7),c(-5,5),c(-4,4)),
                     v.ylims=list(c(-7,7),c(-5,5),c(-4,4)),
                     v.maglim = list(c(0,5.9), c(0,2.9), c(0,2.9)),
                     v.denslim = list(c(0,0.5),c(0,0.4),c(0,0.35)),
                     v.vec.axis <- list(seq(-7,7,2),
                                      seq(-4.5,4.5,2),
                                      seq(-4.5,4.5,2)))


dev.off()
diversityTrendsPlotsF(loc.color = loc.color,
                      model.list = beta.model.list,
                      rich.df = richDf,
                      turnover.df = turnoverDf,
                      fullPlotPath = date.wrap("./Plots/multi gam results beta",".pdf"),
                      xlabs=list(expression("Coordinated directional "*beta*" diversity"),
                                 expression("Coordinated sequential "*beta*" diversity")),
                      ylabs=list(expression("Ratio directional "*beta*" diversity"),
                                 expression("Ratio sequential "*beta*" diversity")),
                      xlims=list(c(0.1,0.9),c(0.1,0.7)),
                      ylims=list(c(-0.35,0.35),c(-0.25,0.3)),
                      
                      vector.list = vectorList[4:5],
                      v.xlabs=list(expression(Delta*" Coordinated directional "*beta*" diversity"),
                                 expression(Delta*" Coordinated sequential "*beta*" diversity")),
                      v.ylabs=list(expression(Delta*" Ratio directional "*beta*" diversity"),
                                 expression(Delta*" Ratio sequential "*beta*" diversity")),
                      v.xlims=list(c(-0.85,0.85),c(-0.85,0.85)),
                      v.ylims=list(c(-0.85,0.85),c(-0.85,0.85)),
                      v.maglim = list(c(0,0.36), c(0,0.36)),
                      v.denslim = list(c(0,0.4),c(0,0.32)),
                      v.vec.axis <- list(seq(-1,1,0.25),
                                       seq(-1,1,0.25)))


# diversityTrendsPlots(loc.color = loc.color,
#                      model.list = alpha.model.list,
#                      rich.df = richDf,
#                      turnover.df = turnoverDf,
#                      fullPlotPath = date.wrap("./Plots/multi gam results_new1",".pdf"),
#                      regionPlotPath = date.wrap("./Plots/regional gam results1",".pdf"),
#                      sitePlotPath = date.wrap("./Plots/site gam results1",".pdf"),
#                      xlabs=list(expression("Coordinated "*alpha*" diversity (Hill q = 0)"),
#                                 expression("Coordinated "*alpha*" diversity (Hill q = 1)"),
#                                 expression("Coordinated "*alpha*" diversity (Hill q = 2)")),
#                      ylabs=list(expression("Ratio "*alpha*" diversity (Hill q = 0)"),
#                                 expression("Ratio "*alpha*" diversity (Hill q = 1)"),
#                                 expression("Ratio "*alpha*" diversity (Hill q = 2)")),
#                      xlims=list(c(4,20),c(3,7.8),c(2,7)),
#                      ylims=list(c(-5,5),c(-2.2,2.2),c(-1.5,1.5)))

#                            Vector plot ####

vectorPlot(path = date.wrap("./Plots/combined vector alpha",".pdf"),
           vector.list = vectorList[1:3],
           xlabs=list(expression("Coordinated "*alpha*" diversity (Hill q = 0)"),
                      expression("Coordinated "*alpha*" diversity (Hill q = 1)"),
                      expression("Coordinated "*alpha*" diversity (Hill q = 2)")),
           ylabs=list(expression("Ratio "*alpha*" diversity (Hill q = 0)"),
                      expression("Ratio "*alpha*" diversity (Hill q = 1)"),
                      expression("Ratio "*alpha*" diversity (Hill q = 2)")),
           xlims=list(c(-6.5,6.5),c(-5,5),c(-4,4)),
           ylims=list(c(-6.5,6.5),c(-5,5),c(-4,4)),
           maglim = list(c(0,5), c(0,3), c(0,2.5)),
           denslim = list(c(0,0.5),c(0,0.4),c(0,0.35)),
           vec.axis <- list(seq(-6,6,2),
                            seq(-4,4,2),
                            seq(-4,4,2)))

vectorPlot(path = date.wrap("./Plots/combined vector beta",".pdf"),
          

#                        Trends with CIs ####

loc.color <- data.frame(locality=unique(huon_coral$locality))
loc.color$col <- c("blue","red","darkgreen")[huon_tr$site[match(loc.color$locality,
                                                                huon_tr$locality)]]

#              raw plots with error bars ####

# Primary axis plots
sapply(1:length(alpha.model.list), function(n1){
  
  print(n1)
  
  region.preds <- alpha.model.list[[n1]]$preds[[1]]
  site.preds <- alpha.model.list[[n1]]$preds[[2]]
  locality.preds <- alpha.model.list[[n1]]$preds[[3]]
  raw.data <- dataList[[n1]]
  
  fit1.ylims<-rbind(c(3,24), c(0.5,10), c(0,8))[n1,]
  fit1.ylab <- c(expression("Coordinated "*alpha*" diversity (Hill q = 0)"),
                 expression("Coordinated "*alpha*" diversity (Hill q = 1)"),
                 expression("Coordinated "*alpha*" diversity (Hill q = 2)"))[n1]
  fit1.var <- c("coordH0", "coordH1", "coordH2")[n1]
  
  x.marks <- seq(6000,9000,1000)
  
  # PRIMARY AXIS PLOTS
  
  # set up split.screen. Locality is 3x3 panel, with Site and region as extra side plots
  pdf(date.wrap(paste0("./Plots/raw div ",fit1.var),".pdf"), 
      height=4.5, width=7, useDingbats = FALSE)
  
  panel3.x <- seq(0.425,0.99,len=4)
  panel3.y <- seq(0.08,0.97,len=4)
  panels <- expand.grid(1:3, 1:3)
  panel.mat <- cbind(panel3.x[-4][panels[,1]],
                     panel3.x[-1][panels[,1]],
                     panel3.y[-4][panels[,2]],
                     panel3.y[-1][panels[,2]])
  
  split.screen(rbind(c(0.065,0.35,0.525,0.97),
                     c(0.065,0.35,0.08,0.525),
                     panel.mat))
  
  # region plot
  screen(1)
  raw.plot(fit1.ylims)
  axis(side=1, at=x.marks, labels=NA)
  axis(side=2, mgp=c(3,0.5,0))
  mtext(side=2, line=1.25, text=fit1.ylab, las=0)
  raw.trend(preds=region.preds,
            fit.var="fit.1",line.col="black", poly.col="grey80")
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.95, "y"),
       labels=paste0("(",LETTERS[1],") Region trend"), font=2, adj=0)
  box()
  close.screen(1)
  
  # site plot 
  screen(2)
  raw.plot(fit1.ylims)
  axis(side=1, labels=NA)
  axis(side=2, mgp=c(3,0.5,0))
  mtext(side=2, line=1.25, text=fit1.ylab, las=0)
  axis(side=1, at=x.marks, labels=format(x.marks, big.mark=","))
  mtext(side=1, line=0.75, text="Years before present")
  a <- sapply(levels(huon_coral$site), function(x){
    
    print(x)
    temp.pred <- site.preds[site.preds$site == x, ]
    site.range <- range(raw.data$pred.date[raw.data$site==x &
                                             !is.na(raw.data[,fit1.var])])
    temp.pred <- temp.pred[temp.pred$pred.date >= site.range[1] &
                             temp.pred$pred.dat <= site.range[2],]
    
    site.col = c("blue", "red", "darkgreen")[levels(raw.data$site) == x]
    poly.col <- col2rgb(site.col)/255
    poly.col <- rgb(poly.col[1], poly.col[2], poly.col[3], 0.2)
    
    raw.trend(temp.pred, "fit.1", site.col, poly.col)
    
    text(x=temp.pred$pred.date[1],
         y=temp.pred$fit.1[1],
         pos=4, offset=0.25,
         labels=substr(temp.pred$site[1], 1, 1),
         col=site.col)
    
  })
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.935, "y"),
       labels=paste0("(",LETTERS[2],") Locality trends"), font=2, adj=0)
  box()
  close.screen(2)
  
  # locality subplots
  sapply(1:9, function(n2){
    screen(2+n2)
    
    raw.plot(fit1.ylims)
    if(n2%in% c(1,4,7)){axis(side=2, mgp=c(3,0.5,0))
    } else {axis(side=2, labels=NA)
    }
    if(n2==4){mtext(side=2, line=1.25, text=fit1.ylab, las=0)}
    
    if(n2 <= 3){
      axis(side=1, at=x.marks, labels=format(x.marks, big.mark=","))
    } else {
      axis(side=1, at=x.marks, labels=NA)
    }
    if(n2==2){mtext(side=1, line=0.75, text="Years before present")}
    if(n2==8){mtext(side=3, line=-0.1, text="Site trends", font=2)}
    
    temp.loc <- locality.number$locality[locality.number$number==n2]
    
    temp.pred <- locality.preds[locality.preds$locality == temp.loc, ]
    site.range <- range(raw.data$pred.date[raw.data$locality==temp.loc &
                                             !is.na(raw.data[,fit1.var])])
    temp.pred <- temp.pred[temp.pred$pred.date >= site.range[1] &
                             temp.pred$pred.dat <= site.range[2],]
    
    site.col = loc.color$col[loc.color$locality == temp.loc]
    poly.col <- col2rgb(site.col)/255
    point.col <- rgb(poly.col[1], poly.col[2], poly.col[3], 0.4) 
    poly.col <- rgb(poly.col[1], poly.col[2], poly.col[3], 0.2)
    
    temp.raw <- raw.data[raw.data$locality==temp.loc,]
    temp.raw$fit.var <- temp.raw[,fit1.var]
    points(temp.raw$fit.var ~ temp.raw$pred.date, pch=16, col=point.col, cex=0.8)
    raw.trend(temp.pred, "fit.1", site.col, poly.col)
    
    text(y=temp.pred$fit.1[1],
         x=temp.pred$pred.date[1],
         pos=ifelse(n2==8, 3, 4), offset=0.25,
         labels=locality.number$number[temp.loc],
         col=site.col)
    
    box()
    text(x=relative.axis.point(0.02, "x"),
         y=relative.axis.point(ifelse(n2>6, 0.94, 0.925), "y"),
         labels=paste0("(",LETTERS[c(9,10,11,6,7,8,3,4,5)[n2]],")"), font=2, adj=0)
    close.screen(2+n2)
  })
  
  close.screen(all.screens=TRUE)
  dev.off()
  
})

# Secondary axis plots
sapply(1:length(alpha.model.list), function(n1){
  
  print(n1)
  
  region.preds <- alpha.model.list[[n1]]$preds[[1]]
  site.preds <- alpha.model.list[[n1]]$preds[[2]]
  locality.preds <- alpha.model.list[[n1]]$preds[[3]]
  raw.data <- dataList[[n1]]
  
  fit1.ylims<-rbind(c(-2,7.5), c(-3.5,3.5), c(-2,2))[n1,]
  fit1.ylab <- c(expression("Ratio "*alpha*" diversity (Hill q = 0)"),
                 expression("Ratio "*alpha*" diversity (Hill q = 1)"),
                 expression("Ratio "*alpha*" diversity (Hill q = 2)"))[n1]
  fit1.var <- c("ratioH0", "ratioH1", "ratioH2")[n1]
  
  x.marks <- seq(6000,9000,1000)
  
  # set up split.screen. Locality is 3x3 panel, with Site and region as extra side plots
  pdf(date.wrap(paste0("./Plots/raw div ",fit1.var),".pdf"), 
      height=4.5, width=7, useDingbats = FALSE)
  
  panel3.x <- seq(0.425,0.99,len=4)
  panel3.y <- seq(0.08,0.97,len=4)
  panels <- expand.grid(1:3, 1:3)
  panel.mat <- cbind(panel3.x[-4][panels[,1]],
                     panel3.x[-1][panels[,1]],
                     panel3.y[-4][panels[,2]],
                     panel3.y[-1][panels[,2]])
  
  split.screen(rbind(c(0.065,0.35,0.525,0.97),
                     c(0.065,0.35,0.08,0.525),
                     panel.mat))
  
  # region plot
  screen(1)
  raw.plot(fit1.ylims)
  abline(h=0, lty="31", col="grey80")
  axis(side=1, at=x.marks, labels=NA)
  axis(side=2, mgp=c(3,0.5,0))
  mtext(side=2, line=1.5, text=fit1.ylab, las=0)
  raw.trend(preds=region.preds,
            fit.var="fit.2",line.col="black", poly.col=rgb(0.3,0.3,0.3,0.2))
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.95, "y"),
       labels=paste0("(",LETTERS[1],") Region trend"), font=2, adj=0)
  box()
  close.screen(1)
  
  # site plot 
  screen(2)
  raw.plot(fit1.ylims)
  abline(h=0, lty="31", col="grey80")
  axis(side=1, labels=NA)
  axis(side=2, mgp=c(3,0.5,0))
  mtext(side=2, line=1.5, text=fit1.ylab, las=0)
  axis(side=1, at=x.marks, labels=format(x.marks, big.mark=","))
  mtext(side=1, line=0.75, text="Years before present")
  a <- sapply(levels(huon_coral$site), function(x){
    
    print(x)
    temp.pred <- site.preds[site.preds$site == x, ]
    site.range <- range(raw.data$pred.date[raw.data$site==x &
                                             !is.na(raw.data[,fit1.var])])
    temp.pred <- temp.pred[temp.pred$pred.date >= site.range[1] &
                             temp.pred$pred.dat <= site.range[2],]
    
    site.col = c("blue", "red", "darkgreen")[levels(raw.data$site) == x]
    poly.col <- col2rgb(site.col)/255
    poly.col <- rgb(poly.col[1], poly.col[2], poly.col[3], 0.2)
    
    raw.trend(temp.pred, "fit.2", site.col, poly.col)
    
    text(x=temp.pred$pred.date[1],
         y=temp.pred$fit.2[1],
         pos=4, offset=0.25,
         labels=substr(temp.pred$site[1], 1, 1),
         col=site.col)
    
  })
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.935, "y"),
       labels=paste0("(",LETTERS[2],") Locality trends"), font=2, adj=0)
  box()
  close.screen(2)
  
  # locality subplots
  sapply(1:9, function(n2){
    screen(2+n2)
    
    raw.plot(fit1.ylims)
    abline(h=0, lty="31", col="grey80")
    if(n2%in% c(1,4,7)){axis(side=2, mgp=c(3,0.5,0))
    } else {axis(side=2, labels=NA)
    }
    if(n2==4){mtext(side=2, line=1.5, text=fit1.ylab, las=0)}
    
    if(n2 <= 3){
      axis(side=1, at=x.marks, labels=format(x.marks, big.mark=","))
    } else {
      axis(side=1, at=x.marks, labels=NA)
    }
    if(n2==2){mtext(side=1, line=0.75, text="Years before present")}
    if(n2==8){mtext(side=3, line=-0.1, text="Site trends", font=2)}
    
    temp.loc <- locality.number$locality[locality.number$number==n2]
    
    temp.pred <- locality.preds[locality.preds$locality == temp.loc, ]
    site.range <- range(raw.data$pred.date[raw.data$locality==temp.loc &
                                             !is.na(raw.data[,fit1.var])])
    temp.pred <- temp.pred[temp.pred$pred.date >= site.range[1] &
                             temp.pred$pred.dat <= site.range[2],]
    
    site.col = loc.color$col[loc.color$locality == temp.loc]
    poly.col <- col2rgb(site.col)/255
    point.col <- rgb(poly.col[1], poly.col[2], poly.col[3], 0.4) 
    poly.col <- rgb(poly.col[1], poly.col[2], poly.col[3], 0.2)
    
    temp.raw <- raw.data[raw.data$locality==temp.loc,]
    temp.raw$fit.var <- temp.raw[,fit1.var]
    points(temp.raw$fit.var ~ temp.raw$pred.date, pch=16, col=point.col, cex=0.8)
    raw.trend(temp.pred, "fit.2", site.col, poly.col)
    
    text(y=temp.pred$fit.2[1],
         x=temp.pred$pred.date[1],
         pos=ifelse(n2==8, 3, 4), offset=0.25,
         labels=locality.number$number[temp.loc],
         col=site.col)
    
    box()
    text(x=relative.axis.point(0.02, "x"),
         y=relative.axis.point(ifelse(n2>6, 0.94, 0.925), "y"),
         labels=paste0("(",LETTERS[c(9,10,11,6,7,8,3,4,5)[n2]],")"), font=2, adj=0)
    close.screen(2+n2)
  })
  
  close.screen(all.screens=TRUE)
  dev.off()
  
})

# Primary axis plots
dataList = list(turnoverDf, turnoverDf)
sapply(1:length(beta.model.list), function(n1){
  
  print(n1)
  
  region.preds <- beta.model.list[[n1]]$preds[[1]]
  site.preds <- beta.model.list[[n1]]$preds[[2]]
  locality.preds <- beta.model.list[[n1]]$preds[[3]]
  raw.data <- dataList[[n1]]
  
  fit1.ylims<-rbind(c(0,1), c(0,0.8))[n1,]
  fit1.ylab <- c(expression("Coordinated directional "*beta*" diversity"),
                 expression("Coordinated sequential "*alpha*" diversity"))[n1]
  fit1.var <- c("dissimilarity", "turnover")[n1]
  
  x.marks <- seq(6000,9000,1000)
  
  # PRIMARY AXIS PLOTS
  
  # set up split.screen. Locality is 3x3 panel, with Site and region as extra side plots
  pdf(date.wrap(paste0("./Plots/raw div ",fit1.var),".pdf"), 
      height=4.5, width=7, useDingbats = FALSE)
  
  panel3.x <- seq(0.425,0.99,len=4)
  panel3.y <- seq(0.08,0.97,len=4)
  panels <- expand.grid(1:3, 1:3)
  panel.mat <- cbind(panel3.x[-4][panels[,1]],
                     panel3.x[-1][panels[,1]],
                     panel3.y[-4][panels[,2]],
                     panel3.y[-1][panels[,2]])
  
  split.screen(rbind(c(0.065,0.35,0.525,0.97),
                     c(0.065,0.35,0.08,0.525),
                     panel.mat))
  
  # region plot
  screen(1)
  raw.plot(fit1.ylims)
  axis(side=1, at=x.marks, labels=NA)
  axis(side=2, mgp=c(3,0.5,0))
  mtext(side=2, line=1.25, text=fit1.ylab, las=0)
  raw.trend(preds=region.preds,
            fit.var="fit.1",line.col="black", poly.col="grey80")
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.95, "y"),
       labels=paste0("(",LETTERS[1],") Region trend"), font=2, adj=0)
  box()
  close.screen(1)
  
  # site plot 
  screen(2)
  raw.plot(fit1.ylims)
  axis(side=1, labels=NA)
  axis(side=2, mgp=c(3,0.5,0))
  mtext(side=2, line=1.25, text=fit1.ylab, las=0)
  axis(side=1, at=x.marks, labels=format(x.marks, big.mark=","))
  mtext(side=1, line=0.75, text="Years before present")
  a <- sapply(levels(huon_coral$site), function(x){
    
    print(x)
    temp.pred <- site.preds[site.preds$site == x, ]
    site.range <- range(raw.data$pred.date[raw.data$site==x &
                                             !is.na(raw.data[,fit1.var])])
    temp.pred <- temp.pred[temp.pred$pred.date >= site.range[1] &
                             temp.pred$pred.dat <= site.range[2],]
    
    site.col = c("blue", "red", "darkgreen")[levels(raw.data$site) == x]
    poly.col <- col2rgb(site.col)/255
    poly.col <- rgb(poly.col[1], poly.col[2], poly.col[3], 0.2)
    
    raw.trend(temp.pred, "fit.1", site.col, poly.col)
    
    text(x=temp.pred$pred.date[1],
         y=temp.pred$fit.1[1],
         pos=4, offset=0.25,
         labels=substr(temp.pred$site[1], 1, 1),
         col=site.col)
    
  })
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.935, "y"),
       labels=paste0("(",LETTERS[2],") Locality trends"), font=2, adj=0)
  box()
  close.screen(2)
  
  # locality subplots
  sapply(1:9, function(n2){
    screen(2+n2)
    
    raw.plot(fit1.ylims)
    if(n2%in% c(1,4,7)){axis(side=2, mgp=c(3,0.5,0))
    } else {axis(side=2, labels=NA)
    }
    if(n2==4){mtext(side=2, line=1.25, text=fit1.ylab, las=0)}
    
    if(n2 <= 3){
      axis(side=1, at=x.marks, labels=format(x.marks, big.mark=","))
    } else {
      axis(side=1, at=x.marks, labels=NA)
    }
    if(n2==2){mtext(side=1, line=0.75, text="Years before present")}
    if(n2==8){mtext(side=3, line=-0.1, text="Site trends", font=2)}
    
    temp.loc <- locality.number$locality[locality.number$number==n2]
    
    temp.pred <- locality.preds[locality.preds$locality == temp.loc, ]
    site.range <- range(raw.data$pred.date[raw.data$locality==temp.loc &
                                             !is.na(raw.data[,fit1.var])])
    temp.pred <- temp.pred[temp.pred$pred.date >= site.range[1] &
                             temp.pred$pred.dat <= site.range[2],]
    
    site.col = loc.color$col[loc.color$locality == temp.loc]
    poly.col <- col2rgb(site.col)/255
    point.col <- rgb(poly.col[1], poly.col[2], poly.col[3], 0.4) 
    poly.col <- rgb(poly.col[1], poly.col[2], poly.col[3], 0.2)
    
    temp.raw <- raw.data[raw.data$locality==temp.loc,]
    temp.raw$fit.var <- temp.raw[,fit1.var]
    points(temp.raw$fit.var ~ temp.raw$pred.date, pch=16, col=point.col, cex=0.8)
    raw.trend(temp.pred, "fit.1", site.col, poly.col)
    
    text(y=temp.pred$fit.1[1],
         x=temp.pred$pred.date[1],
         pos=ifelse(n2==8, 3, 4), offset=0.25,
         labels=locality.number$number[temp.loc],
         col=site.col)
    
    box()
    text(x=relative.axis.point(0.02, "x"),
         y=relative.axis.point(ifelse(n2>6, 0.94, 0.925), "y"),
         labels=paste0("(",LETTERS[c(9,10,11,6,7,8,3,4,5)[n2]],")"), font=2, adj=0)
    close.screen(2+n2)
  })
  
  close.screen(all.screens=TRUE)
  dev.off()
  
})

# Secondary axis plots
sapply(1:length(beta.model.list), function(n1){
  
  print(n1)
  
  region.preds <- beta.model.list[[n1]]$preds[[1]]
  site.preds <- beta.model.list[[n1]]$preds[[2]]
  locality.preds <- beta.model.list[[n1]]$preds[[3]]
  raw.data <- dataList[[n1]]
  
  fit1.ylims<-rbind(c(-0.3,0.3), c(-0.3,0.3))[n1,]
  fit1.ylab <- c(expression("Ratio directional "*beta*" diversity"),
                 expression("Ratio sequential "*beta*" diversity"))[n1]
  fit1.var <- c("conservation", "persistance")[n1]
  
  x.marks <- seq(6000,9000,1000)
  
  # set up split.screen. Locality is 3x3 panel, with Site and region as extra side plots
  pdf(date.wrap(paste0("./Plots/raw div ",fit1.var),".pdf"), 
      height=4.5, width=7, useDingbats = FALSE)
  
  panel3.x <- seq(0.425,0.99,len=4)
  panel3.y <- seq(0.08,0.97,len=4)
  panels <- expand.grid(1:3, 1:3)
  panel.mat <- cbind(panel3.x[-4][panels[,1]],
                     panel3.x[-1][panels[,1]],
                     panel3.y[-4][panels[,2]],
                     panel3.y[-1][panels[,2]])
  
  split.screen(rbind(c(0.065,0.35,0.525,0.97),
                     c(0.065,0.35,0.08,0.525),
                     panel.mat))
  
  # region plot
  screen(1)
  raw.plot(fit1.ylims)
  abline(h=0, lty="31", col="grey80")
  axis(side=1, at=x.marks, labels=NA)
  axis(side=2, mgp=c(3,0.5,0))
  mtext(side=2, line=1.5, text=fit1.ylab, las=0)
  raw.trend(preds=region.preds,
            fit.var="fit.2",line.col="black", poly.col=rgb(0.3,0.3,0.3,0.2))
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.95, "y"),
       labels=paste0("(",LETTERS[1],") Region trend"), font=2, adj=0)
  box()
  close.screen(1)
  
  # site plot 
  screen(2)
  raw.plot(fit1.ylims)
  abline(h=0, lty="31", col="grey80")
  axis(side=1, labels=NA)
  axis(side=2, mgp=c(3,0.5,0))
  mtext(side=2, line=1.5, text=fit1.ylab, las=0)
  axis(side=1, at=x.marks, labels=format(x.marks, big.mark=","))
  mtext(side=1, line=0.75, text="Years before present")
  a <- sapply(levels(huon_coral$site), function(x){
    
    print(x)
    temp.pred <- site.preds[site.preds$site == x, ]
    site.range <- range(raw.data$pred.date[raw.data$site==x &
                                             !is.na(raw.data[,fit1.var])])
    temp.pred <- temp.pred[temp.pred$pred.date >= site.range[1] &
                             temp.pred$pred.dat <= site.range[2],]
    
    site.col = c("blue", "red", "darkgreen")[levels(raw.data$site) == x]
    poly.col <- col2rgb(site.col)/255
    poly.col <- rgb(poly.col[1], poly.col[2], poly.col[3], 0.2)
    
    raw.trend(temp.pred, "fit.2", site.col, poly.col)
    
    text(x=temp.pred$pred.date[1],
         y=temp.pred$fit.2[1],
         pos=4, offset=0.25,
         labels=substr(temp.pred$site[1], 1, 1),
         col=site.col)
    
  })
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.935, "y"),
       labels=paste0("(",LETTERS[2],") Locality trends"), font=2, adj=0)
  box()
  close.screen(2)
  
  # locality subplots
  sapply(1:9, function(n2){
    screen(2+n2)
    
    raw.plot(fit1.ylims)
    abline(h=0, lty="31", col="grey80")
    if(n2%in% c(1,4,7)){axis(side=2, mgp=c(3,0.5,0))
    } else {axis(side=2, labels=NA)
    }
    if(n2==4){mtext(side=2, line=1.5, text=fit1.ylab, las=0)}
    
    if(n2 <= 3){
      axis(side=1, at=x.marks, labels=format(x.marks, big.mark=","))
    } else {
      axis(side=1, at=x.marks, labels=NA)
    }
    if(n2==2){mtext(side=1, line=0.75, text="Years before present")}
    if(n2==8){mtext(side=3, line=-0.1, text="Site trends", font=2)}
    
    temp.loc <- locality.number$locality[locality.number$number==n2]
    
    temp.pred <- locality.preds[locality.preds$locality == temp.loc, ]
    site.range <- range(raw.data$pred.date[raw.data$locality==temp.loc &
                                             !is.na(raw.data[,fit1.var])])
    temp.pred <- temp.pred[temp.pred$pred.date >= site.range[1] &
                             temp.pred$pred.dat <= site.range[2],]
    
    site.col = loc.color$col[loc.color$locality == temp.loc]
    poly.col <- col2rgb(site.col)/255
    point.col <- rgb(poly.col[1], poly.col[2], poly.col[3], 0.4) 
    poly.col <- rgb(poly.col[1], poly.col[2], poly.col[3], 0.2)
    
    temp.raw <- raw.data[raw.data$locality==temp.loc,]
    temp.raw$fit.var <- temp.raw[,fit1.var]
    points(temp.raw$fit.var ~ temp.raw$pred.date, pch=16, col=point.col, cex=0.8)
    raw.trend(temp.pred, "fit.2", site.col, poly.col)
    
    text(y=temp.pred$fit.2[1],
         x=temp.pred$pred.date[1],
         pos=ifelse(n2==8, 3, 4), offset=0.25,
         labels=locality.number$number[temp.loc],
         col=site.col)
    
    box()
    text(x=relative.axis.point(0.02, "x"),
         y=relative.axis.point(ifelse(n2>6, 0.94, 0.925), "y"),
         labels=paste0("(",LETTERS[c(9,10,11,6,7,8,3,4,5)[n2]],")"), font=2, adj=0)
    close.screen(2+n2)
  })
  
  close.screen(all.screens=TRUE)
  dev.off()
  
})

#                  Alpha simulation comparison plot ####

pdf(date.wrap("./plots/Alpha simulation", ".pdf"), height=6, width=6)
par(mfrow=c(3,3), mar=c(0,0,0,0), oma=c(4,4,1,1), ps=10, tcl=-0.25, las=1, mgp=c(3,0.5,0))

sapply(9:1, function(n){

  plot(x=NULL, y=NULL, 
       xlim=c(9500, 5000), 
       ylim=c(0,16), xaxt="n", yaxt="n", yaxs="i", xlab="", ylab="")
  
  # locality obs data
  tempLoc <- locs[locs.order[n]]
  locRange <- range(huon_site$transect.age[huon_site$locality == tempLoc])
  
  tempObs <- alpha.model.list[[2]][[2]][[3]]
  tempObs <- tempObs[tempObs$locality == tempLoc,]
  tempObs <- tempObs[tempObs$pred.date > locRange[1] &
                     tempObs$pred.date < locRange[2],]
    
  
  # PLOT MARKOV TRENDS
  
  markCol <- col2rgb("purple")/255
  markovTrends <- sapply(markovAlphaModels, function(x){
    
    if(is.null(x)){return(NULL)}
    xSub <- x[[2]][[3]]
    xSub <- xSub[xSub$locality == tempLoc,]
    xSub <- xSub[xSub$pred.date > locRange[1] &
                   xSub$pred.date < locRange[2],]
    
      lines(xSub$fit.1 ~ xSub$pred.date, lwd=1, col=rgb(markCol[1],markCol[2],markCol[3],0.02))
      return(xSub$fit.1)
    })
  
  polygon(x=c(tempObs$pred.date, rev(tempObs$pred.date)),
          y=c(apply(markovTrends, 1, quantile, prob = 0.025),
              rev(apply(markovTrends, 1, quantile, prob = 0.975))),
          border=NA, col=rgb(markCol[1],markCol[2],markCol[3],0.25))
  lines(rowMeans(markovTrends) ~ 
          tempObs$pred.date, col=rgb(markCol[1],markCol[2],markCol[3],1), lwd=1.5)
  text(y=rowMeans(markovTrends)[1],
       x=tempObs$pred.date[1],
       col=rgb(markCol[1],markCol[2],markCol[3],1), labels="Markov", font=2, pos=4)
  
  # PLOT PERM TRENDS
  permCol <- col2rgb("sienna4")/255
  permatTrends <- sapply(permatAlphaModels, function(x){
    
    xSub <- x[[2]][[3]]
    xSub <- xSub[xSub$locality == tempLoc,]
    xSub <- xSub[xSub$pred.date > locRange[1] &
                   xSub$pred.date < locRange[2],]
    
    lines(xSub$fit.1 ~ xSub$pred.date, lwd=1, col=rgb(permCol[1],permCol[2],permCol[3],0.02))
    return(xSub$fit.1)
  })
  
  polygon(x=c(tempObs$pred.date, rev(tempObs$pred.date)),
          y=c(apply(permatTrends, 1, quantile, prob = 0.025),
              rev(apply(permatTrends, 1, quantile, prob = 0.975))),
          border=NA, col=rgb(permCol[1],permCol[2],permCol[3],0.25))
  lines(rowMeans(permatTrends) ~ 
          tempObs$pred.date, col=rgb(permCol[1],permCol[2],permCol[3],1), lwd=1.5)
  text(y=rowMeans(permatTrends)[1],
       x=tempObs$pred.date[1],
       col=rgb(permCol[1],permCol[2],permCol[3],1), labels="Perm", font=2, pos=4)
  
  # PLOT OBS TRENDS
  polygon(x=c(tempObs$pred.date, rev(tempObs$pred.date)),
         y=c(tempObs$fit.1 +
               1.96 * tempObs$se.fit.1, rev(tempObs$fit.1 -
                   1.96 * tempObs$se.fit.1)),
         border=NA, col=rgb(0,0,0,0.25))
  
  lines(tempObs$fit.1 ~ tempObs$pred.date, lwd=2)
    
  text(x=relative.axis.point(0.02, "x"),
         y=relative.axis.point(0.95, "y"),
       labels=paste0("(", LETTERS[n], ") Site ", 10-n, ": ", tempLoc), font=2, adj=0)
  
  if(n %in% 1:3){axis(side=1, at=seq(6000,9500,500),
                   labels=format(seq(6000,9500,500), big.mark=","),
                   mgp=c(3,0.2,0))
    
    if(n==2){
    mtext(side=1, line=1.5, text="Years before present", cex=0.8)
    }
  } else{axis(side=1, labels=NA)}
  
  if(n %in% c(3,6,9)){axis(side=2)
    
    if(n==6){
    mtext(side=2, line=1.5, text = expression("Taxonomic "*alpha*" diversity (Hill q=1)"), cex=0.8, las=0)
    }
    
  } else{axis(side=2, labels=NA)}

})
dev.off()

pdf(date.wrap("./plots/Alpha simulation ratio", ".pdf"), height=6, width=6)
par(mfrow=c(3,3), mar=c(0,0,0,0), oma=c(4,4,1,1), ps=10, tcl=-0.25, las=1, mgp=c(3,0.5,0))

sapply(9:1, function(n){
  
  plot(x=NULL, y=NULL, 
       xlim=c(9500, 5000), 
       ylim=c(0,6), xaxt="n", yaxt="n", yaxs="i", xlab="", ylab="")
  
  # locality obs data
  tempLoc <- locs[locs.order[n]]
  locRange <- range(huon_site$transect.age[huon_site$locality == tempLoc])
  
  tempObs <- alpha.model.list[[2]][[2]][[3]]
  tempObs <- tempObs[tempObs$locality == tempLoc,]
  tempObs <- tempObs[tempObs$pred.date > locRange[1] &
                       tempObs$pred.date < locRange[2],]
  
  
  # PLOT MARKOV TRENDS
  
  markCol <- col2rgb("purple")/255
  markovTrends <- sapply(markovAlphaModels, function(x){
    
    xSub <- x[[2]][[3]]
    xSub <- xSub[xSub$locality == tempLoc,]
    xSub <- xSub[xSub$pred.date > locRange[1] &
                   xSub$pred.date < locRange[2],]
    
    lines(xSub$fit.2 ~ xSub$pred.date, lwd=1, col=rgb(markCol[1],markCol[2],markCol[3],0.02))
    return(xSub$fit.2)
  })
  
  polygon(x=c(tempObs$pred.date, rev(tempObs$pred.date)),
          y=c(apply(markovTrends, 1, quantile, prob = 0.025),
              rev(apply(markovTrends, 1, quantile, prob = 0.975))),
          border=NA, col=rgb(markCol[1],markCol[2],markCol[3],0.25))
  lines(rowMeans(markovTrends) ~ 
          tempObs$pred.date, col=rgb(markCol[1],markCol[2],markCol[3],1), lwd=1.5)
  text(y=rowMeans(markovTrends)[1],
       x=tempObs$pred.date[1],
       col=rgb(markCol[1],markCol[2],markCol[3],1), labels="Markov", font=2, pos=4)
  
  # PLOT PERM TRENDS
  permCol <- col2rgb("sienna4")/255
  permatTrends <- sapply(permatAlphaModels, function(x){
    
    xSub <- x[[2]][[3]]
    xSub <- xSub[xSub$locality == tempLoc,]
    xSub <- xSub[xSub$pred.date > locRange[1] &
                   xSub$pred.date < locRange[2],]
    
    lines(xSub$fit.2 ~ xSub$pred.date, lwd=1, col=rgb(permCol[1],permCol[2],permCol[3],0.02))
    return(xSub$fit.2)
  })
  
  polygon(x=c(tempObs$pred.date, rev(tempObs$pred.date)),
          y=c(apply(permatTrends, 1, quantile, prob = 0.025),
              rev(apply(permatTrends, 1, quantile, prob = 0.975))),
          border=NA, col=rgb(permCol[1],permCol[2],permCol[3],0.25))
  lines(rowMeans(permatTrends) ~ 
          tempObs$pred.date, col=rgb(permCol[1],permCol[2],permCol[3],1), lwd=1.5)
  text(y=rowMeans(permatTrends)[1],
       x=tempObs$pred.date[1],
       col=rgb(permCol[1],permCol[2],permCol[3],1), labels="Perm", font=2, pos=4)
  
  # PLOT OBS TRENDS
  polygon(x=c(tempObs$pred.date, rev(tempObs$pred.date)),
          y=c(tempObs$fit.2 +
                1.96 * tempObs$se.fit.2, rev(tempObs$fit.2 -
                                               1.96 * tempObs$se.fit.2)),
          border=NA, col=rgb(0,0,0,0.25))
  
  lines(tempObs$fit.2 ~ tempObs$pred.date, lwd=2)
  
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.95, "y"),
       labels=paste0("(", LETTERS[n], ") Site ", 10-n, ": ", tempLoc), font=2, adj=0)
  
  if(n %in% 1:3){axis(side=1, at=seq(6000,9500,500),
                      labels=format(seq(6000,9500,500), big.mark=","),
                      mgp=c(3,0.2,0))
    
    if(n==2){
      mtext(side=1, line=1.5, text="Years before present", cex=0.8)
    }
  } else{axis(side=1, labels=NA)}
  
  if(n %in% c(3,6,9)){axis(side=2)
    
    if(n==6){
      mtext(side=2, line=1.5, text = expression("Functional "*alpha*" diversity (Hill q=1)"), cex=0.8, las=0)
    }
    
  } else{axis(side=2, labels=NA)}
  
})
dev.off()

#                   Beta simulation comparison plot ####

pdf(date.wrap("./plots/Beta simulation coord dir", ".pdf"), height=6, width=6)
par(mfrow=c(3,3), mar=c(0,0,0,0), oma=c(4,4,1,1), ps=10, tcl=-0.25, las=1, mgp=c(3,0.5,0))

sapply(9:1, function(n){
  
  plot(x=NULL, y=NULL, 
       xlim=c(9500, 5000), 
       ylim=c(0,1), xaxt="n", yaxt="n", yaxs="i", xlab="", ylab="")
  
  # locality obs data
  tempLoc <- locs[locs.order[n]]
  locRange <- range(huon_site$transect.age[huon_site$locality == tempLoc])
  
  tempObs <- beta.model.list[[1]][[2]][[3]]
  tempObs <- tempObs[tempObs$locality == tempLoc,]
  tempObs <- tempObs[tempObs$pred.date > locRange[1] &
                       tempObs$pred.date < locRange[2],]
  
  # PLOT MARKOV TRENDS
  
  markCol <- col2rgb("purple")/255
  markovTrends <- sapply(markovBetaModels, function(x){
    
    xSub <- x[[1]][[3]]
    xSub <- xSub[xSub$locality == tempLoc,]
    xSub <- xSub[xSub$pred.date > locRange[1] &
                   xSub$pred.date < locRange[2],]
    if(is.null(xSub$fit.1)){return(NULL)}
    
    lines(xSub$fit.1 ~ xSub$pred.date, lwd=1, col=rgb(markCol[1],markCol[2],markCol[3],0.02))
    return(xSub$fit.1)
  })
  markovTrends <- t(do.call("rbind", markovTrends[!sapply(markovTrends, is.null)]))
  
  polygon(x=c(tempObs$pred.date, rev(tempObs$pred.date)),
          y=c(apply(markovTrends, 1, quantile, prob = 0.025),
              rev(apply(markovTrends, 1, quantile, prob = 0.975))),
          border=NA, col=rgb(markCol[1],markCol[2],markCol[3],0.25))
  lines(rowMeans(markovTrends) ~ 
          tempObs$pred.date, col=rgb(markCol[1],markCol[2],markCol[3],1), lwd=1.5)
  text(y=rowMeans(markovTrends)[1],
       x=tempObs$pred.date[1],
       col=rgb(markCol[1],markCol[2],markCol[3],1), labels="Markov", font=2, pos=4)
  
  # PLOT PERM TRENDS
  permCol <- col2rgb("chocolate1")/255
  permatTrends <- sapply(permatBetaModels, function(x){
    
    xSub <- x[[1]][[3]]
    xSub <- xSub[xSub$locality == tempLoc,]
    xSub <- xSub[xSub$pred.date > locRange[1] &
                   xSub$pred.date < locRange[2],]
    
    lines(xSub$fit.1 ~ xSub$pred.date, lwd=1, col=rgb(permCol[1],permCol[2],permCol[3],0.02))
    return(xSub$fit.1)
  })
  
  polygon(x=c(tempObs$pred.date, rev(tempObs$pred.date)),
          y=c(apply(permatTrends, 1, quantile, prob = 0.025),
              rev(apply(permatTrends, 1, quantile, prob = 0.975))),
          border=NA, col=rgb(permCol[1],permCol[2],permCol[3],0.25))
  lines(rowMeans(permatTrends) ~ 
          tempObs$pred.date, col=rgb(permCol[1],permCol[2],permCol[3],1), lwd=1.5)
  text(y=rowMeans(permatTrends)[1],
       x=tempObs$pred.date[1],
       col=rgb(permCol[1],permCol[2],permCol[3],1), labels="Perm", font=2, pos=4)
  
  # PLOT OBS TRENDS
  polygon(x=c(tempObs$pred.date, rev(tempObs$pred.date)),
          y=c(tempObs$fit.1 +
                1.96 * tempObs$se.fit.1, rev(tempObs$fit.1 -
                                               1.96 * tempObs$se.fit.1)),
          border=NA, col=rgb(0,0,0,0.25))
  
  lines(tempObs$fit.1 ~ tempObs$pred.date, lwd=2)
  
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.95, "y"),
       labels=paste0("(", LETTERS[n], ") Site ", 10-n, ": ", tempLoc), font=2, adj=0)
  
  if(n %in% 1:3){axis(side=1, at=seq(6000,9500,500),
                      labels=format(seq(6000,9500,500), big.mark=","),
                      mgp=c(3,0.2,0))
    
    if(n==2){
      mtext(side=1, line=1.5, text="Years before present", cex=0.8)
    }
  } else{axis(side=1, labels=NA)}
  
  if(n %in% c(3,6,9)){axis(side=2)
    
    if(n==6){
      mtext(side=2, line=1.5, text = expression("Taxonomic directional "*beta*" diversity"), cex=0.8, las=0)
    }
    
  } else{axis(side=2, labels=NA)}
  
})
dev.off()


pdf(date.wrap("./plots/Beta simulation ratio dir", ".pdf"), height=6, width=6)
par(mfrow=c(3,3), mar=c(0,0,0,0), oma=c(4,4,1,1), ps=10, tcl=-0.25, las=1, mgp=c(3,0.5,0))

sapply(9:1, function(n){
  
  plot(x=NULL, y=NULL, 
       xlim=c(9500, 5000), 
       ylim=c(-0.6,0.6), xaxt="n", yaxt="n", yaxs="i", xlab="", ylab="")
  
  # locality obs data
  tempLoc <- locs[locs.order[n]]
  locRange <- range(huon_site$transect.age[huon_site$locality == tempLoc])
  
  tempObs <- beta.model.list[[1]][[2]][[3]]
  tempObs <- tempObs[tempObs$locality == tempLoc,]
  tempObs <- tempObs[tempObs$pred.date > locRange[1] &
                       tempObs$pred.date < locRange[2],]
  
  # PLOT MARKOV TRENDS
  
  markCol <- col2rgb("purple")/255
  markovTrends <- sapply(markovBetaModels, function(x){
    
    xSub <- x[[1]][[3]]
    xSub <- xSub[xSub$locality == tempLoc,]
    xSub <- xSub[xSub$pred.date > locRange[1] &
                   xSub$pred.date < locRange[2],]
    
    if(is.null(xSub$fit.2)){return(NULL)}
    
    lines(xSub$fit.2 ~ xSub$pred.date, lwd=1, col=rgb(markCol[1],markCol[2],markCol[3],0.02))
    return(xSub$fit.2)
  })
  markovTrends <- t(do.call("rbind", markovTrends[!sapply(markovTrends, is.null)]))
  
  polygon(x=c(tempObs$pred.date, rev(tempObs$pred.date)),
          y=c(apply(markovTrends, 1, quantile, prob = 0.025),
              rev(apply(markovTrends, 1, quantile, prob = 0.975))),
          border=NA, col=rgb(markCol[1],markCol[2],markCol[3],0.25))
  lines(rowMeans(markovTrends) ~ 
          tempObs$pred.date, col=rgb(markCol[1],markCol[2],markCol[3],1), lwd=1.5)
  text(y=rowMeans(markovTrends)[1],
       x=tempObs$pred.date[1],
       col=rgb(markCol[1],markCol[2],markCol[3],1), labels="Markov", font=2, pos=4)
  
  # PLOT PERM TRENDS
  permCol <- col2rgb("chocolate1")/255
  permatTrends <- sapply(permatBetaModels, function(x){
    
    xSub <- x[[1]][[3]]
    xSub <- xSub[xSub$locality == tempLoc,]
    xSub <- xSub[xSub$pred.date > locRange[1] &
                   xSub$pred.date < locRange[2],]
    
    lines(xSub$fit.2 ~ xSub$pred.date, lwd=1, col=rgb(permCol[1],permCol[2],permCol[3],0.02))
    return(xSub$fit.2)
  })
  
  polygon(x=c(tempObs$pred.date, rev(tempObs$pred.date)),
          y=c(apply(permatTrends, 1, quantile, prob = 0.025),
              rev(apply(permatTrends, 1, quantile, prob = 0.975))),
          border=NA, col=rgb(permCol[1],permCol[2],permCol[3],0.25))
  lines(rowMeans(permatTrends) ~ 
          tempObs$pred.date, col=rgb(permCol[1],permCol[2],permCol[3],1), lwd=1.5)
  text(y=rowMeans(permatTrends)[1],
       x=tempObs$pred.date[1],
       col=rgb(permCol[1],permCol[2],permCol[3],1), labels="Perm", font=2, pos=4)
  
  # PLOT OBS TRENDS
  polygon(x=c(tempObs$pred.date, rev(tempObs$pred.date)),
          y=c(tempObs$fit.2 +
                1.96 * tempObs$se.fit.2, rev(tempObs$fit.2 -
                                               1.96 * tempObs$se.fit.2)),
          border=NA, col=rgb(0,0,0,0.25))
  
  lines(tempObs$fit.2 ~ tempObs$pred.date, lwd=2)
  
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.95, "y"),
       labels=paste0("(", LETTERS[n], ") Site ", 10-n, ": ", tempLoc), font=2, adj=0)
  
  if(n %in% 1:3){axis(side=1, at=seq(6000,9500,500),
                      labels=format(seq(6000,9500,500), big.mark=","),
                      mgp=c(3,0.2,0))
    
    if(n==2){
      mtext(side=1, line=1.5, text="Years before present", cex=0.8)
    }
  } else{axis(side=1, labels=NA)}
  
  if(n %in% c(3,6,9)){axis(side=2)
    
    if(n==6){
      mtext(side=2, line=1.75, text = expression("Functional directional "*beta*" diversity"), cex=0.8, las=0)
    }
    
  } else{axis(side=2, labels=NA)}
  
})
dev.off()


pdf(date.wrap("./plots/Beta simulation coord seq", ".pdf"), height=6, width=6)
par(mfrow=c(3,3), mar=c(0,0,0,0), oma=c(4,4,1,1), ps=10, tcl=-0.25, las=1, mgp=c(3,0.5,0))

sapply(1:9, function(n){
  
  plot(x=NULL, y=NULL, 
       xlim=c(9500, 5000), 
       ylim=c(0,1), xaxt="n", yaxt="n", yaxs="i", xlab="", ylab="")
  
  # locality obs data
  tempLoc <- locs[locs.order[n]]
  locRange <- range(huon_site$transect.age[huon_site$locality == tempLoc])
  
  tempObs <- beta.model.list[[2]][[2]][[3]]
  tempObs <- tempObs[tempObs$locality == tempLoc,]
  tempObs <- tempObs[tempObs$pred.date > locRange[1] &
                       tempObs$pred.date < locRange[2],]
  
  # PLOT MARKOV TRENDS
  
  markCol <- col2rgb("purple")/255
  markovTrends <- sapply(markovBetaModels, function(x){
    
    xSub <- x[[2]][[3]]
    xSub <- xSub[xSub$locality == tempLoc,]
    xSub <- xSub[xSub$pred.date > locRange[1] &
                   xSub$pred.date < locRange[2],]
    
    lines(xSub$fit.1 ~ xSub$pred.date, lwd=1, col=rgb(markCol[1],markCol[2],markCol[3],0.02))
    return(xSub$fit.1)
  })
  
  polygon(x=c(tempObs$pred.date, rev(tempObs$pred.date)),
          y=c(apply(markovTrends, 1, quantile, prob = 0.025),
              rev(apply(markovTrends, 1, quantile, prob = 0.975))),
          border=NA, col=rgb(markCol[1],markCol[2],markCol[3],0.25))
  lines(rowMeans(markovTrends) ~ 
          tempObs$pred.date, col=rgb(markCol[1],markCol[2],markCol[3],1), lwd=1.5)
  text(y=rowMeans(markovTrends)[1],
       x=tempObs$pred.date[1],
       col=rgb(markCol[1],markCol[2],markCol[3],1), labels="Markov", font=2, pos=4)
  
  # PLOT PERM TRENDS
  permCol <- col2rgb("chocolate1")/255
  permatTrends <- sapply(permatBetaModels, function(x){
    
    xSub <- x[[2]][[3]]
    xSub <- xSub[xSub$locality == tempLoc,]
    xSub <- xSub[xSub$pred.date > locRange[1] &
                   xSub$pred.date < locRange[2],]
    
    lines(xSub$fit.1 ~ xSub$pred.date, lwd=1, col=rgb(permCol[1],permCol[2],permCol[3],0.02))
    return(xSub$fit.1)
  })
  
  polygon(x=c(tempObs$pred.date, rev(tempObs$pred.date)),
          y=c(apply(permatTrends, 1, quantile, prob = 0.025),
              rev(apply(permatTrends, 1, quantile, prob = 0.975))),
          border=NA, col=rgb(permCol[1],permCol[2],permCol[3],0.25))
  lines(rowMeans(permatTrends) ~ 
          tempObs$pred.date, col=rgb(permCol[1],permCol[2],permCol[3],1), lwd=1.5)
  text(y=rowMeans(permatTrends)[1],
       x=tempObs$pred.date[1],
       col=rgb(permCol[1],permCol[2],permCol[3],1), labels="Perm", font=2, pos=4)
  
  # PLOT OBS TRENDS
  polygon(x=c(tempObs$pred.date, rev(tempObs$pred.date)),
          y=c(tempObs$fit.1 +
                1.96 * tempObs$se.fit.1, rev(tempObs$fit.1 -
                                               1.96 * tempObs$se.fit.1)),
          border=NA, col=rgb(0,0,0,0.25))
  
  lines(tempObs$fit.1 ~ tempObs$pred.date, lwd=2)
  
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.95, "y"),
       labels=paste0("(", LETTERS[n], ") ", tempLoc), font=2, adj=0)
  
  if(n %in% 7:9){axis(side=1, at=seq(6000,9500,500),
                      labels=format(seq(6000,9500,500), big.mark=","),
                      mgp=c(3,0.2,0))
    
    if(n==8){
      mtext(side=1, line=1.5, text="Years before present", cex=0.8)
    }
  } else{axis(side=1, labels=NA)}
  
  if(n %in% c(1,4,7)){axis(side=2)
    
    if(n==4){
      mtext(side=2, line=1.5, text = expression("Taxonomic sequential "*beta*" diversity"), cex=0.8, las=0)
    }
    
  } else{axis(side=2, labels=NA)}
  
})
dev.off()


pdf(date.wrap("./plots/Beta simulation ratio seq", ".pdf"), height=6, width=6)
par(mfrow=c(3,3), mar=c(0,0,0,0), oma=c(4,4,1,1), ps=10, tcl=-0.25, las=1, mgp=c(3,0.5,0))

sapply(1:9, function(n){
  
  plot(x=NULL, y=NULL, 
       xlim=c(9500, 5000), 
       ylim=c(-0.6,0.6), xaxt="n", yaxt="n", yaxs="i", xlab="", ylab="")
  
  # locality obs data
  tempLoc <- locs[locs.order[n]]
  locRange <- range(huon_site$transect.age[huon_site$locality == tempLoc])
  
  tempObs <- beta.model.list[[2]][[2]][[3]]
  tempObs <- tempObs[tempObs$locality == tempLoc,]
  tempObs <- tempObs[tempObs$pred.date > locRange[1] &
                       tempObs$pred.date < locRange[2],]
  
  # PLOT MARKOV TRENDS
  
  markCol <- col2rgb("purple")/255
  markovTrends <- sapply(markovBetaModels, function(x){
    
    xSub <- x[[2]][[3]]
    xSub <- xSub[xSub$locality == tempLoc,]
    xSub <- xSub[xSub$pred.date > locRange[1] &
                   xSub$pred.date < locRange[2],]
    
    lines(xSub$fit.2 ~ xSub$pred.date, lwd=1, col=rgb(markCol[1],markCol[2],markCol[3],0.02))
    return(xSub$fit.2)
  })
  
  polygon(x=c(tempObs$pred.date, rev(tempObs$pred.date)),
          y=c(apply(markovTrends, 1, quantile, prob = 0.025),
              rev(apply(markovTrends, 1, quantile, prob = 0.975))),
          border=NA, col=rgb(markCol[1],markCol[2],markCol[3],0.25))
  lines(rowMeans(markovTrends) ~ 
          tempObs$pred.date, col=rgb(markCol[1],markCol[2],markCol[3],1), lwd=1.5)
  text(y=rowMeans(markovTrends)[1],
       x=tempObs$pred.date[1],
       col=rgb(markCol[1],markCol[2],markCol[3],1), labels="Markov", font=2, pos=4)
  
  # PLOT PERM TRENDS
  permCol <- col2rgb("chocolate1")/255
  permatTrends <- sapply(permatBetaModels, function(x){
    
    xSub <- x[[2]][[3]]
    xSub <- xSub[xSub$locality == tempLoc,]
    xSub <- xSub[xSub$pred.date > locRange[1] &
                   xSub$pred.date < locRange[2],]
    
    lines(xSub$fit.2 ~ xSub$pred.date, lwd=1, col=rgb(permCol[1],permCol[2],permCol[3],0.02))
    return(xSub$fit.2)
  })
  
  polygon(x=c(tempObs$pred.date, rev(tempObs$pred.date)),
          y=c(apply(permatTrends, 1, quantile, prob = 0.025),
              rev(apply(permatTrends, 1, quantile, prob = 0.975))),
          border=NA, col=rgb(permCol[1],permCol[2],permCol[3],0.25))
  lines(rowMeans(permatTrends) ~ 
          tempObs$pred.date, col=rgb(permCol[1],permCol[2],permCol[3],1), lwd=1.5)
  text(y=rowMeans(permatTrends)[1],
       x=tempObs$pred.date[1],
       col=rgb(permCol[1],permCol[2],permCol[3],1), labels="Perm", font=2, pos=4)
  
  # PLOT OBS TRENDS
  polygon(x=c(tempObs$pred.date, rev(tempObs$pred.date)),
          y=c(tempObs$fit.2 +
                1.96 * tempObs$se.fit.2, rev(tempObs$fit.2 -
                                               1.96 * tempObs$se.fit.2)),
          border=NA, col=rgb(0,0,0,0.25))
  
  lines(tempObs$fit.2 ~ tempObs$pred.date, lwd=2)
  
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.95, "y"),
       labels=paste0("(", LETTERS[n], ") ", tempLoc), font=2, adj=0)
  
  if(n %in% 7:9){axis(side=1, at=seq(6000,9500,500),
                      labels=format(seq(6000,9500,500), big.mark=","),
                      mgp=c(3,0.2,0))
    
    if(n==8){
      mtext(side=1, line=1.5, text="Years before present", cex=0.8)
    }
  } else{axis(side=1, labels=NA)}
  
  if(n %in% c(1,4,7)){axis(side=2)
    
    if(n==4){
      mtext(side=2, line=1.75, text = expression("Functional sequential "*beta*" diversity"), cex=0.8, las=0)
    }
    
  } else{axis(side=2, labels=NA)}
  
})
dev.off()

#            Matrix of genus/growth form ####

gr.gen.mat <- table(huon_coral_all$genus, huon_coral_all$growth.form)
gr.gen.mat <- gr.gen.mat[rowSums(gr.gen.mat) > 0, ]
gen.count <- rowSums(gr.gen.mat)

#Other <- colSums(gr.gen.mat[gen.count < 10, ])
#commons <- gr.gen.mat[gen.count >= 10, ]
#final.mat <- rbind(commons, Other)
final.mat <- gr.gen.mat
final.mat <- final.mat[, c(2,1,3:ncol(final.mat))]

rownames(final.mat)[rownames(final.mat) == "Other"] = paste0("Other (n=", length(Other), ")")
gen.count <- rowSums(final.mat) / sum(final.mat)
gr.count <- colSums(final.mat) / sum(final.mat)

gr.label.mat <- data.frame(col = colnames(final.mat),
                           label = c("Stout\nbranching",
                                     "Branching",
                                     "Encrusting", 
                                     "Massive,\nSubmassive,\nColumnar", 
                                     "Tabular,\nPlate,\nLaminar", 
                                     "Unknown"),
                           pos = 1:ncol(final.mat), stringsAsFactors = FALSE)

final.mat <- final.mat[nrow(final.mat):1, order(gr.label.mat$pos)]
colnames(final.mat) <- gr.label.mat$label[order(gr.label.mat$pos)]
gr.count <- gr.count[order(gr.label.mat$pos)]
#final.mat <- final.mat / rowSums(final.mat)

pdf(date.wrap("./Plots/genus growth form matrix", ".pdf"), height=7.5, width=6)
par(mar=c(0,0,0,0), ps=10)
plot.new()

gen.y <- 0.1
gr.x <- 0.36
end.x <- 0.95
end.y <- 1

coord.mat <- expand.grid(x.cent = seq(gr.x, end.x, len=ncol(final.mat)),
                         y.cent = seq(gen.y, end.y, len=nrow(final.mat)))
x.width <- 0.5 * diff(coord.mat[1:2,1])
y.width <- 0.5 * max(diff(coord.mat[,2]))

coord.mat$xleft <- coord.mat$x.cent - x.width
coord.mat$xright <- coord.mat$x.cent + x.width
coord.mat$ybottom <- coord.mat$y.cent - y.width
coord.mat$ytop <- coord.mat$y.cent + y.width

# color grads
main.grad <- c("white",
               colorRampPalette(c(rgb(0.9,0.9,1),"blue","darkblue","black"), alpha=TRUE)(max(final.mat)))

genus.perc <- rev(round(gen.count, 2))*100
gr.perc <- rev(round(gr.count, 2))*100
genus.grad <- rev(colorRampPalette(c("red4", "white"))(max(genus.perc)+1))
gr.grad <- colorRampPalette(c("white", "darkgreen"))(max(gr.perc)+1)

# grid lines
segments(y0 = unique(coord.mat$y.cent),
         y1 = unique(coord.mat$y.cent),
         x0 =  gr.x - 3.5 * x.width - 0.0075,
         x1 = gr.x - 3.5 * x.width)

segments(y0 = gen.y -  0.045,
         y1 = gen.y - 0.055,
         x0 =  unique(coord.mat$x.cent),
         x1 = unique(coord.mat$x.cent))

# vertical genus rectangles
rect(xleft = gr.x - 3.5 * x.width,
     xright = gr.x - 1.5 * x.width,
     ybottom = seq(gen.y, end.y, len=nrow(final.mat)) - y.width,
     ytop = seq(gen.y, end.y, len=nrow(final.mat)) + y.width,
     col= genus.grad[genus.perc+1],)

rect(xleft = seq(gr.x, end.x, len=ncol(final.mat)) - x.width,
     xright = seq(gr.x, end.x, len=ncol(final.mat)) + x.width,
     ybottom = gen.y - 4 * y.width,
     ytop = gen.y - 2 * y.width,
     col= gr.grad[rev(gr.perc)+1])

# big matrix
rect(xleft = coord.mat$xleft,
     xright = coord.mat$xright,
     ybottom = coord.mat$ybottom,
     ytop = coord.mat$ytop,
     col= main.grad[t(final.mat+1)])

# Props
final.cols <- t(final.mat)
text(x = coord.mat$x.cent,
     y = coord.mat$y.cent,
     labels = final.cols, cex=0.75,
     col = ifelse(final.cols > 300 | final.cols == 0, "white", "black"), font=2)

# gen labels
sapply(1:nrow(final.mat), function(n){
  
  text(x = gr.x - 3.5 * x.width,
       y = seq(gen.y, end.y, len=nrow(final.mat))[n],
       labels = rownames(final.mat)[n], pos = 2, offset=0.5,
       font = 3)
})

# gen text
text(x = gr.x - 2.5 * x.width,
     y = seq(gen.y, end.y, len=nrow(final.mat)),
     labels = sprintf("%.2f", round(rev(gen.count), 3)*100), cex=0.75,
     col = ifelse(rev(gen.count) > 0.5, "white", "black"), font=2)

# gr labels
text(x = seq(gr.x, end.x, len=ncol(final.mat)),
     y = gen.y -  0.06,
     labels = colnames(final.mat), srt=0, adj=c(0.5,1))

# growth prop
text(x = seq(gr.x, end.x, len=ncol(final.mat)),
     y =  gen.y - 3 * y.width,
     labels = sprintf("%.2f", round(gr.count, 3)*100), cex=0.75,
     col = ifelse(gr.count > 0.25, "white", "black"), font=2)

# Totals text

text(x=gr.x - 2.5 * x.width, 
     y = gen.y - 3 * y.width, 
     labels="Totals\n(%)", font=2)

dev.off()


fullTable <- as.data.frame(table(huon_coral_all$growth.comb, huon_coral_all$genSp))
fullTable <- merge(fullTable, huon_coral_all[!duplicated(huon_coral_all$growth.comb),c("growth.form", "growth.comb")],
                   by.x="Var1", by.y="growth.comb", all.x=TRUE, all.y=FALSE, sort=FALSE)
fullTable <- fullTable[fullTable$Freq > 0,]
fullTable <- fullTable[order(fullTable$Freq, decreasing=TRUE),]
colnames(fullTable) <- c("growth assess", "species", "count", "growth.form")
write.csv(fullTable, "./outputs/fullGrowthSummaryFormed.csv")



fullTable <- as.data.frame(table(huon_coral_all$growth.form, huon_coral_all$genSp))
fullTable <- fullTable[fullTable$Freq > 0,]
fullTable <- fullTable[order(fullTable$Freq, decreasing=TRUE),]
write.csv(fullTable, "./outputs/speciesGrowth.csv")

#   Summary stats for manuscript ####

# transect lengths
huon_tr_sub <- huon_tr[huon_tr$locality %in% locs,]
tran.length <- tapply(huon_tr_sub$int.delta, huon_tr_sub$transect, sum)
summary(tran.length[tran.length > 7])


# number of transects per locality
table(droplevels(huon_site$locality[huon_site$locality %in% locs]))

test <- genusMatProp[rownames(genusMatProp) %in% huon_site$transect[huon_site$locality == "Loto Beach"],]
test[,colSums(test)>0]

test <- grMatProp[rownames(grMatProp) %in% huon_site$transect[huon_site$locality == "Loto Beach"],]
test[,colSums(test)>0]

test <- richDf[richDf$locality == "Loto Beach",]
test[order(test$pred.date),]

huon_site[huon_site$locality == "Bonah River",]

# time between transects

sd(turnoverDf$date.diff, na.rm=TRUE)

# inflation of top dissims
summary(turnoverDf$seq.gen.diss / turnoverDf$top.gen.diss)
summary(turnoverDf$seq.gr.diss / turnoverDf$top.gr.diss)


# ####
# Compare dir and seq beta vectors ####

turnover.vect$beta.diff1 <- turnover.vect$d.diss - turnover.vect$d.turn
turnover.vect$beta.diff2 <- turnover.vect$d.cons - turnover.vect$d.pers

cor(turnover.vect[,c("d.diss", "d.turn")], use="complete.obs")
cor.test(x=turnover.vect$d.diss, y=turnover.vect$d.turn)
plot(turnover.vect$d.diss ~ turnover.vect$d.turn)

cor(turnover.vect[,c("d.cons", "d.pers")], use="complete.obs")
cor.test(x=turnover.vect$d.cons, y=turnover.vect$d.pers)
plot(turnover.vect$d.cons ~ turnover.vect$d.pers)

# Axis correlations ####
#       Compare tax and fun alpha richness ####
H0cor <- lmer(grH0 ~ genusH0 + (1|site/locality), data=richDf)
r2(H0cor)
summary(H0cor)
H0corC <- summary(H0cor)$coefficients
cbind(H0corC[,1] - 1.96 * H0corC[,2],
      H0corC[,1] + 1.96 * H0corC[,2])

H1cor <- lmer(grH1 ~ genusH1 + (1|site/locality), data=richDf)
r2(H1cor)
summary(H1cor)
H1corC <- summary(H1cor)$coefficients
cbind(H1corC[,1] - 1.96 * H1corC[,2],
      H1corC[,1] + 1.96 * H1corC[,2])

H2cor <- lmer(grH2 ~ genusH2 + (1|site/locality), data=richDf)
r2(H2cor)
summary(H2cor)
H2corC <- summary(H2cor)$coefficients
cbind(H2corC[,1] - 1.96 * H2corC[,2],
      H2corC[,1] + 1.96 * H2corC[,2])

dirBcor <- lmer(top.gr.diss ~ top.gen.diss + (1|site/locality), data=turnoverDf)
r2(dirBcor)
summary(dirBcor)
dirBcorC <- summary(dirBcor)$coefficients
cbind(dirBcorC[,1] - 1.96 * dirBcorC[,2],
      dirBcorC[,1] + 1.96 * dirBcorC[,2])

seqBcor <- lmer(seq.gr.diss ~ seq.gen.diss + (1|site/locality), data=turnoverDf)
r2(seqBcor)
summary(seqBcor)
seqBcorC <- summary(seqBcor)$coefficients
cbind(seqBcorC[,1] - 1.96 * seqBcorC[,2],
      seqBcorC[,1] + 1.96 * seqBcorC[,2])

#       directional versus sequential ####
taxSeqDir <- lmer(top.gen.diss ~ seq.gen.diss + (1|site/locality),
               data=turnoverDf)
r2(taxSeqDir)
taxSeqDirC <- summary(taxSeqDir)$coefficients
cbind(taxSeqDirC[,1] - 1.96 * taxSeqDirC[,2],
      taxSeqDirC[,1] + 1.96 * taxSeqDirC[,2])

funSeqDir <- lmer(top.gr.diss ~ seq.gr.diss + (1|site/locality),
                    data=turnoverDf)
r2(funSeqDir)
funSeqDirC <- summary(funSeqDir)$coefficients
cbind(funSeqDirC[,1] - 1.96 * funSeqDirC[,2],
      funSeqDirC[,1] + 1.96 * funSeqDirC[,2])

# ####
# Genus addition over time ####
uptime <- model.list[[1]][[2]][[1]]

pdf("./Plots/genus occupancy.pdf", height=7, width=4, useDingbats = FALSE)

split.screen(rbind(c(0.23,0.98,0.065,0.85),
                   c(0.23,0.98,0.85,0.99)))

screen(1)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, las=1, mgp=c(3,0.5,0))
plot(x=NULL, y=NULL, 
     xlim=rev(range(richDf$pred.date)), 
     ylim=rev(c(1.25, ncol(genusMatProp)+0.75)),
     axes=FALSE, xlab="", ylab="")

uptime1 <- which(diff(uptime$fit.1) < 0)[1]
downtime1 <- which(diff(uptime$fit.1) > 0)
downtime1 <- downtime1[downtime1 > uptime1][1]
uptime2 <- which(diff(uptime$fit.1) < 0)
uptime2 <- uptime2[uptime2 > downtime1][1]

rect(xleft=uptime$pred.date[c(1, uptime1, downtime1, uptime2)],
     xright = uptime$pred.date[c(uptime1, downtime1, uptime2, nrow(uptime))],
     ybottom = par("usr")[3],
     ytop = par("usr")[4], border="grey80",
     col=c(rgb(1,0,0,0.1), rgb(0,0,1,0.1))[c(1,2,1,2)])

sapply(1:ncol(genusMatProp), function(n){
  
  print(n)
  gen.sub <- genusMatProp[,n]
  
  if(sum(gen.sub > 0) > 1){
  segments(x0=min(richDf$pred.dat[gen.sub > 0]),
           x1=max(richDf$pred.dat[gen.sub > 0]),
           y0=n, y1=n, lwd=1.5)
  }
  
  points(y=rep(n, sum(gen.sub > 0)), x=richDf$pred.dat[gen.sub > 0],
         lwd=0.5, pch=21, bg=rgb(0,0,0,1), cex=0.5)
         #cex=0.5+ sqrt(gen.sub[gen.sub>0])
  })
axis(side=2, at=1:ncol(genusMatProp), 
     labels=colnames(genusMatProp), font=3, las=1)
axis(side=1, mgp=c(3,0,0), at = seq(6000,9500,500),
     labels=format(seq(6000,9500,500), big.mark=","))
box()
mtext(side=1, line=1, text="Years before present")
close.screen(1)

screen(2)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, las=1, mgp=c(3,0.5,0))
plot(x=NULL, y=NULL, 
     xlim=rev(range(richDf$pred.date)), 
     ylim=c(4.9,8),
     axes=FALSE, xlab="", ylab="")

rect(xleft=uptime$pred.date[c(1, uptime1, downtime1, uptime2)],
     xright = uptime$pred.date[c(uptime1, downtime1, uptime2, nrow(uptime))],
     ybottom = par("usr")[3],
     ytop = par("usr")[4], border="grey80",
     col=c(rgb(1,0,0,0.1), rgb(0,0,1,0.1))[c(1,2,1,2)])

polygon(x=c(uptime$pred.date, rev(uptime$pred.date)),
        y=c(uptime$fit.1 + 1.96*uptime$se.fit.1,
            rev(uptime$fit.1 - 1.96*uptime$se.fit.1)),
        border=NA, col=rgb(0.5,0.5,0.5,0.5))
lines(x=uptime$pred.date, y=uptime$fit.1, lwd=2)

box()

axis(side=1, labels=NA)
axis(side=2)
mtext(side=2, line=1.5, text="Coordinated\nalpha diversity", las=0)
close.screen(2)
dev.off()

# Growth from occupancy ####


uptime <- model.list[[1]][[2]][[1]]

pdf("./Plots/growth form occupancy.pdf", height=7, width=4, useDingbats = FALSE)

split.screen(rbind(c(0.23,0.98,0.065,0.85),
                   c(0.23,0.98,0.85,0.99)))

screen(1)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, las=1, mgp=c(3,0.5,0))
plot(x=NULL, y=NULL, 
     xlim=rev(range(richDf$pred.date)), 
     ylim=rev(c(1.25, ncol(grMatProp)+0.75)),
     axes=FALSE, xlab="", ylab="")

uptime1 <- which(diff(uptime$fit.1) < 0)[1]
downtime1 <- which(diff(uptime$fit.1) > 0)
downtime1 <- downtime1[downtime1 > uptime1][1]
uptime2 <- which(diff(uptime$fit.1) < 0)
uptime2 <- uptime2[uptime2 > downtime1][1]

rect(xleft=uptime$pred.date[c(1, uptime1, downtime1, uptime2)],
     xright = uptime$pred.date[c(uptime1, downtime1, uptime2, nrow(uptime))],
     ybottom = par("usr")[3],
     ytop = par("usr")[4], border="grey80",
     col=c(rgb(1,0,0,0.1), rgb(0,0,1,0.1))[c(1,2,1,2)])

sapply(1:ncol(grMatProp), function(n){
  print(n)
  gen.sub <- grMatProp[,n]
  
  if(sum(gen.sub > 0) > 1){
    segments(x0=min(richDf$pred.dat[gen.sub > 0]),
             x1=max(richDf$pred.dat[gen.sub > 0]),
             y0=n, y1=n, lwd=1.5)
  }
  
  points(y=rep(n, sum(gen.sub > 0)), x=richDf$pred.dat[gen.sub > 0],
         lwd=0.5, pch=21, bg=rgb(0,0,0,1), cex=0.5)
  #cex=0.5+ sqrt(gen.sub[gen.sub>0])
})
axis(side=2, at=1:ncol(grMatProp), 
     labels=colnames(grMatProp), font=3, las=1)
axis(side=1, mgp=c(3,0,0), at = seq(6000,9500,500),
     labels=format(seq(6000,9500,500), big.mark=","))
box()
mtext(side=1, line=1, text="Years before present")
close.screen(1)

screen(2)
par(mar=c(0,0,0,0), ps=8, tcl=-0.25, las=1, mgp=c(3,0.5,0))
plot(x=NULL, y=NULL, 
     xlim=rev(range(richDf$pred.date)), 
     ylim=c(4.9,8),
     axes=FALSE, xlab="", ylab="")

rect(xleft=uptime$pred.date[c(1, uptime1, downtime1, uptime2)],
     xright = uptime$pred.date[c(uptime1, downtime1, uptime2, nrow(uptime))],
     ybottom = par("usr")[3],
     ytop = par("usr")[4], border="grey80",
     col=c(rgb(1,0,0,0.1), rgb(0,0,1,0.1))[c(1,2,1,2)])

polygon(x=c(uptime$pred.date, rev(uptime$pred.date)),
        y=c(uptime$fit.1 + 1.96*uptime$se.fit.1,
            rev(uptime$fit.1 - 1.96*uptime$se.fit.2)),
        border=NA, col=rgb(0.5,0.5,0.5,0.5))
lines(x=uptime$pred.date, y=uptime$fit.1, lwd=2)

box()

axis(side=1, labels=NA)
axis(side=2)
mtext(side=2, line=1.5, text="Coordinated\nalpha diversity", las=0)
close.screen(2)
dev.off()

# PLOTS ####
# ####
# GAM PLOTS ####
# ####
# Spine Plot ####
     
# ####
# Example ordination ####
     
# ####
# ####
# ####
# Div rotation graph ####
# ####
     
pdf(date.wrap("./plots/diversity axis rotation plot", ".pdf"), height=7.15, width=5,
    useDingbats = FALSE)
     
     par(mfrow=c(3,2), mar=c(2,3,1,1), oma=c(2,1,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
     
     plot(x=NULL, y=NULL, xlim=c(1,10), ylim=c(1,10), xlab="", ylab="", asp=1,
          xaxt="n")
     axis(side=1, mgp=c(3,0.2,0))
     
     mtext(side=1, line=1.25, text="Growth form diversity", cex=0.8)
     mtext(side=2, line=1.5, text="Genus diversity", las=0, cex=0.8)
     
     text(x = relative.axis.point(c(0.025), "x"), y = relative.axis.point(c(0.025), "y"),
          labels = c("Simple\nsystem"), font=2,  col="grey80", adj=c(0,0))
     
     text(x = relative.axis.point(c(0.975), "x"), y = relative.axis.point(c(0.975), "y"),
          labels = c("Complex\nsystem"), font=2,  
          col="grey30", adj=c(1,1))
     
     Arrows(x0=relative.axis.point(c(0.2), "x"),
            y0=relative.axis.point(c(0.2), "y"),
            x1=relative.axis.point(c(0.8), "x"),
            y1=relative.axis.point(c(0.8), "y"),
            arr.type="triangle", arr.width=0.25, arr.length=0.25,
            col="black", lwd=2)
     
     text(x=relative.axis.point(0.19, "x"),
          y=relative.axis.point(0.22, "y"), srt=45,
          labels="Axis of complexity", adj=c(0,0), col="black", cex=0.9)
     
     gradient.rect(xleft=relative.axis.point(0.15, "x"),
                   xright=relative.axis.point(0.85, "x"),
                   ybottom=relative.axis.point(0.15, "y"),
                   ytop = relative.axis.point(0.85, "y"),
                   col = colorRampPalette(c(rgb(1,1,1,0.6), rgb(1,1,1,0.1)), alpha=TRUE)(5000),
                   border=NA)
     
     Arrows(x0=relative.axis.point(c(0.8), "x"),
            y0=relative.axis.point(c(0.2), "y"),
            x1=relative.axis.point(c(0.2), "x"),
            y1=relative.axis.point(c(0.8), "y"),
            arr.type="triangle", arr.width=0.25, arr.length=0.25,
            col="red", lwd=2)
     
     text(x = relative.axis.point(c(0.975), "x"), y = relative.axis.point(c(0.025), "y"),
          labels = c("Vulnerable\nsystem"), font=2,  
          col=rgb(255/255, 150/255, 150/255), adj=c(1,0))
     
     text(x = relative.axis.point(c(0.025), "x"), y = relative.axis.point(c(0.975), "y"),
          labels = c("Redundant\nsystem"), font=2,  col="red", adj=c(0,1))
     
     text(x=relative.axis.point(0.81, "x"),
          y=relative.axis.point(0.22, "y"), srt=-45,
          labels="Axis of redundancy", adj=c(1,0), col="red", cex=0.9)
     
     
     gradient.rect(xleft=relative.axis.point(0.15, "x"),
                   xright=relative.axis.point(0.85, "x"),
                   ybottom=relative.axis.point(0.15, "y"),
                   ytop = relative.axis.point(0.85, "y"),
                   col = colorRampPalette(c(rgb(1,1,1,0.85), rgb(1,1,1,0.1)), alpha=TRUE)(5000),
                   border=NA, gradient = "y")
     
     plot(x=NULL, y=NULL, xlim=c(4.25,8.5), ylim=c(-1,1.5), xlab="", ylab="", asp=1,
          xaxt="n")
     
     axis(side=1, mgp=c(3,0.25,0))
     mtext(side=1, line=1.25, text="Complexity", cex=0.8)
     mtext(side=2, line=1.5, text="Redundancy", las=0, cex=0.8)
     
     library(shape)
     
     loc.color <- data.frame(locality=unique(huon_coral$locality))
     loc.color$col <- c("blue","red","darkgreen")[huon_tr$site[match(loc.color$locality,
                                                                     huon_tr$locality)]]
     rect(xleft=par("usr")[1], xright=par("usr")[2],
          ybottom=par("usr")[3], ytop=par("usr")[4],
          border=NA, col=rgb(1,1,1,0.7))
     
     gr.pred <- grrich.all.m[[2]]
     gen.pred <- genrich.all.m[[2]]
     
     complexity <- (gr.pred$fit * cos(45 * (pi/180))) + (gen.pred$fit * sin(45 * (pi/180)))
     redundancy <- (-gr.pred$fit * sin(45 * (pi/180))) + (gen.pred$fit * cos(45 * (pi/180)))
     
     Arrows(x0=complexity[2],
            y0=redundancy[2],
            x1=complexity[1],
            y1=redundancy[1],
            arr.type="triangle", arr.width=0.25, arr.length=0.25)
     
     lines(redundancy ~ complexity, lwd=1.5)
     
     sapply(unique(huon_coral$locality), function(x){
       
       gr.pred <- grrich.loc.m[[2]]
       gr.pred <- gr.pred[gr.pred$locality == x, ]
       
       gen.pred <- genrich.loc.m[[2]]
       gen.pred <- gen.pred[gen.pred$locality == x, ]
       
       complexity <- (gr.pred$fit * cos(45 * (pi/180))) + (gen.pred$fit * sin(45 * (pi/180)))
       redundancy <- (-gr.pred$fit * sin(45 * (pi/180))) + (gen.pred$fit * cos(45 * (pi/180)))
       
       lines(redundancy ~ complexity,
             col=loc.color$col[loc.color$locality==x])  
       
       Arrows(x0=complexity[2],
              y0=redundancy[2],
              x1=complexity[1],
              y1=redundancy[1],
              arr.type="triangle", arr.width=0.125, arr.length=0.125,
              col=loc.color$col[loc.color$locality==x])
       
       text(x=complexity[1],
            y=redundancy[1],
            pos=4, offset=0.25,
            labels=locality.number$number[locality.number$locality == x],
            col=loc.color$col[loc.color$locality==x])
       
     })
     
     # convert top similarity onto overlap and replacement axes
     
     plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", asp=1,
          xaxt="n", yaxs="i", xaxs="i")
     axis(side=1, mgp=c(3,0.2,0))
     
     mtext(side=1, line=1.25, text="Growth form dissimilarity to first composition", cex=0.8)
     mtext(side=2, line=2, text="Genus disimilarity to first composition", las=0, cex=0.8)
     
     text(x = relative.axis.point(c(0.025), "x"), y = relative.axis.point(c(0.025), "y"),
          labels = c("Identical\ncomposition"), font=2,  col="grey80", adj=c(0,0))
     
     text(x = relative.axis.point(c(0.975), "x"), y = relative.axis.point(c(0.975), "y"),
          labels = c("Different\ncomposition"), font=2,  
          col="grey30", adj=c(1,1))
     
     Arrows(x0=relative.axis.point(c(0.2), "x"),
            y0=relative.axis.point(c(0.2), "y"),
            x1=relative.axis.point(c(0.8), "x"),
            y1=relative.axis.point(c(0.8), "y"),
            arr.type="triangle", arr.width=0.25, arr.length=0.25,
            col="black", lwd=2)
     
     text(x=relative.axis.point(0.19, "x"),
          y=relative.axis.point(0.22, "y"), srt=45,
          labels="Axis of dissimilarity", adj=c(0,0), col="black", cex=0.9)
     
     gradient.rect(xleft=relative.axis.point(0.15, "x"),
                   xright=relative.axis.point(0.85, "x"),
                   ybottom=relative.axis.point(0.15, "y"),
                   ytop = relative.axis.point(0.85, "y"),
                   col = colorRampPalette(c(rgb(1,1,1,0.6), rgb(1,1,1,0.1)), alpha=TRUE)(5000),
                   border=NA)
     
     Arrows(x0=relative.axis.point(c(0.8), "x"),
            y0=relative.axis.point(c(0.2), "y"),
            x1=relative.axis.point(c(0.2), "x"),
            y1=relative.axis.point(c(0.8), "y"),
            arr.type="triangle", arr.width=0.25, arr.length=0.25,
            col="red", lwd=2)
     
     text(x = relative.axis.point(c(0.975), "x"), y = relative.axis.point(c(0.025), "y"),
          labels = c("Taxonomic\nconservation"), font=2,  
          col=rgb(255/255, 150/255, 150/255), adj=c(1,0))
     
     text(x = relative.axis.point(c(0.025), "x"), y = relative.axis.point(c(0.975), "y"),
          labels = c("Functional\nconservation"), font=2,  col="red", adj=c(0,1))
     
     text(x=relative.axis.point(0.81, "x"),
          y=relative.axis.point(0.22, "y"), srt=-45,
          labels="Axis of conservation", adj=c(1,0), col="red", cex=0.9)
     
     gradient.rect(xleft=relative.axis.point(0.15, "x"),
                   xright=relative.axis.point(0.85, "x"),
                   ybottom=relative.axis.point(0.15, "y"),
                   ytop = relative.axis.point(0.85, "y"),
                   col = colorRampPalette(c(rgb(1,1,1,0.85), rgb(1,1,1,0.1)), alpha=TRUE)(5000),
                   border=NA, gradient = "y")
     
     plot(x=NULL, y=NULL, xlim=c(0.35,1.1), ylim=c(-0.375,0.35), xlab="", ylab="", asp=1,
          xaxt="n")
     
     axis(side=1, mgp=c(3,0.25,0))
     mtext(side=1, line=1.25, text="Dissimilarity", cex=0.8)
     mtext(side=2, line=2, text="Conservation", las=0, cex=0.8)
     
     library(shape)
     
     loc.color <- data.frame(locality=unique(huon_coral$locality))
     loc.color$col <- c("blue","red","darkgreen")[huon_tr$site[match(loc.color$locality,
                                                                     huon_tr$locality)]]
     rect(xleft=par("usr")[1], xright=par("usr")[2],
          ybottom=par("usr")[3], ytop=par("usr")[4],
          border=NA, col=rgb(1,1,1,0.7))
     
     gr.pred <- gr.top.beta.all.m[[2]]
     gen.pred <- genus.top.beta.all.m[[2]]
     
     overlap <- (gr.pred$fit * cos(45 * (pi/180))) + (gen.pred$fit * sin(45 * (pi/180)))
     replacement <- (-gr.pred$fit * sin(45 * (pi/180))) + (gen.pred$fit * cos(45 * (pi/180)))
     
     Arrows(x0=overlap[2],
            y0=replacement[2],
            x1=overlap[1],
            y1=replacement[1],
            arr.type="triangle", arr.width=0.25, arr.length=0.25)
     
     lines(replacement ~ overlap, lwd=1.5)
     
     sapply(unique(huon_coral$locality), function(x){
       
       gr.pred <- gr.top.beta.loc.m[[2]]
       gr.pred <- gr.pred[gr.pred$locality == x, ]
       
       gen.pred <- genus.top.beta.loc.m[[2]]
       gen.pred <- gen.pred[gen.pred$locality == x, ]
       
       complexity <- (gr.pred$fit * cos(45 * (pi/180))) + (gen.pred$fit * sin(45 * (pi/180)))
       redundancy <- (-gr.pred$fit * sin(45 * (pi/180))) + (gen.pred$fit * cos(45 * (pi/180)))
       
       lines(redundancy ~ complexity,
             col=loc.color$col[loc.color$locality==x])  
       
       Arrows(x0=complexity[2],
              y0=redundancy[2],
              x1=complexity[1],
              y1=redundancy[1],
              arr.type="triangle", arr.width=0.125, arr.length=0.125,
              col=loc.color$col[loc.color$locality==x])
       
       text(x=complexity[1],
            y=redundancy[1],
            pos=4, offset=0.25,
            labels=locality.number$number[locality.number$locality == x],
            col=loc.color$col[loc.color$locality==x])
       
     })
     
     
     plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", asp=1,
          xaxt="n", yaxs="i", xaxs="i")
     axis(side=1, mgp=c(3,0.2,0))
     
     mtext(side=1, line=1.5, text="Growth form dissimilarity (per 100 years)", cex=0.8)
     mtext(side=2, line=2, text="Genus dissimilarity (per 100 years)", las=0, cex=0.8)
     
     text(x = relative.axis.point(c(0.025), "x"), y = relative.axis.point(c(0.025), "y"),
          labels = c("Slow, coordinated\nturnover"), font=2,  col="grey80", adj=c(0,0))
     
     text(x = relative.axis.point(c(0.975), "x"), y = relative.axis.point(c(0.975), "y"),
          labels = c("Rapid, coordinated\nturnover"), font=2,  
          col="grey30", adj=c(1,1))
     
     Arrows(x0=relative.axis.point(c(0.2), "x"),
            y0=relative.axis.point(c(0.2), "y"),
            x1=relative.axis.point(c(0.8), "x"),
            y1=relative.axis.point(c(0.8), "y"),
            arr.type="triangle", arr.width=0.25, arr.length=0.25,
            col="black", lwd=2)
     
     text(x=relative.axis.point(0.19, "x"),
          y=relative.axis.point(0.22, "y"), srt=45,
          labels="Axis of turnover", adj=c(0,0), col="black", cex=0.9)
     
     gradient.rect(xleft=relative.axis.point(0.15, "x"),
                   xright=relative.axis.point(0.85, "x"),
                   ybottom=relative.axis.point(0.15, "y"),
                   ytop = relative.axis.point(0.85, "y"),
                   col = colorRampPalette(c(rgb(1,1,1,0.6), rgb(1,1,1,0.1)), alpha=TRUE)(5000),
                   border=NA)
     
     Arrows(x0=relative.axis.point(c(0.8), "x"),
            y0=relative.axis.point(c(0.2), "y"),
            x1=relative.axis.point(c(0.2), "x"),
            y1=relative.axis.point(c(0.8), "y"),
            arr.type="triangle", arr.width=0.25, arr.length=0.25,
            col="red", lwd=2)
     
     text(x = relative.axis.point(c(0.975), "x"), y = relative.axis.point(c(0.025), "y"),
          labels = c("Taxonomic\npersistence"), font=2,
          col=rgb(255/255, 150/255, 150/255), adj=c(1,0))
     
     text(x = relative.axis.point(c(0.025), "x"), y = relative.axis.point(c(0.975), "y"),
          labels = c("Functional\npersistence"), font=2,
          col="red", adj=c(0,1))
     
     text(x=relative.axis.point(0.81, "x"),
          y=relative.axis.point(0.22, "y"), srt=-45,
          labels="Axis of persistence", adj=c(1,0), col="red", cex=0.9)
     
     gradient.rect(xleft=relative.axis.point(0.15, "x"),
                   xright=relative.axis.point(0.85, "x"),
                   ybottom=relative.axis.point(0.15, "y"),
                   ytop = relative.axis.point(0.85, "y"),
                   col = colorRampPalette(c(rgb(1,1,1,0.85), rgb(1,1,1,0.1)), alpha=TRUE)(5000),
                   border=NA, gradient = "y")
     
     plot(x=NULL, y=NULL, xlim=c(0.3,0.85), ylim=c(-0.15,0.25), xlab="", ylab="", asp=1,
          xaxt="n")
     
     axis(side=1, mgp=c(3,0.25,0))
     mtext(side=1, line=1.5, text="Turnover", cex=0.8)
     mtext(side=2, line=2, text="Persistence", las=0, cex=0.8)
     
     library(shape)
     
     loc.color <- data.frame(locality=unique(huon_coral$locality))
     loc.color$col <- c("blue","red","darkgreen")[huon_tr$site[match(loc.color$locality,
                                                                     huon_tr$locality)]]
     rect(xleft=par("usr")[1], xright=par("usr")[2],
          ybottom=par("usr")[3], ytop=par("usr")[4],
          border=NA, col=rgb(1,1,1,0.7))
     
     gr.pred <- gr.region.beta.all.m[[2]]
     gen.pred <- gr.region.beta.all.m[[2]]
     
     overlap <- (gr.pred$fit * cos(45 * (pi/180))) + (gen.pred$fit * sin(45 * (pi/180)))
     replacement <- (-gr.pred$fit * sin(45 * (pi/180))) + (gen.pred$fit * cos(45 * (pi/180)))
     
     Arrows(x0=overlap[2],
            y0=replacement[2],
            x1=overlap[1],
            y1=replacement[1],
            arr.type="triangle", arr.width=0.25, arr.length=0.25)
     
     lines(replacement ~ overlap, lwd=1.5)
     
     sapply(unique(huon_coral$locality), function(x){
       
       gr.pred <- gr.region.beta.loc.m[[2]]
       gr.pred <- gr.pred[gr.pred$locality == x, ]
       
       gen.pred <- genus.region.beta.loc.m[[2]]
       gen.pred <- gen.pred[gen.pred$locality == x, ]
       
       complexity <- (gr.pred$fit * cos(45 * (pi/180))) + (gen.pred$fit * sin(45 * (pi/180)))
       redundancy <- (-gr.pred$fit * sin(45 * (pi/180))) + (gen.pred$fit * cos(45 * (pi/180)))
       
       lines(redundancy ~ complexity,
             col=loc.color$col[loc.color$locality==x])  
       
       Arrows(x0=complexity[2],
              y0=redundancy[2],
              x1=complexity[1],
              y1=redundancy[1],
              arr.type="triangle", arr.width=0.125, arr.length=0.125,
              col=loc.color$col[loc.color$locality==x])
       
       text(x=complexity[1],
            y=redundancy[1],
            pos=4, offset=0.25,
            labels=locality.number$number[locality.number$locality == x],
            col=loc.color$col[loc.color$locality==x])
       
     })
     
     dev.off()
     

#                         Overall diversity trends ####
# this is a 2x2 plot of alpha and beta diversity for growth form and genus richness.

pdf(date.wrap("./plots/gen diversity models", ".pdf"), height=5.25, width=6.25, useDingbats = FALSE)

par(mfrow=c(3,3), mar=c(0,0,0,0), oma=c(3.5,5,2,1), 
    ps=10, mgp=c(3,0.5,0), tcl = -0.25, las=1)

loc.models <- list(genrich.loc.m,
                   genus.top.beta.loc.m,
                   genus.region.beta.loc.m)

site.models <- list(genrich.site.m,
                   genus.top.beta.site.m,
                   genus.region.beta.site.m)

all.models <- list(genrich.all.m,
                   genus.top.beta.all.m,
                   genus.region.beta.all.m)

sapply(1:3, function(n){

  # set up plot
  if(n == 1){ylims=c(1.75,7)}
  if(n == 2){ylims=c(0,1.05)}
  if(n == 3){ylims=c(0,0.8)}
  
  # LOCALITY TRENDS 
  
  plot(x=NULL, y=NULL, xlim=c(9500,5900), ylim=ylims, 
       axes=FALSE, yaxs="i", xlab="", ylab="")
  
  axis(side=2)
  mtext(side=2, line=2, 
          text=c("Rarefied richness (n = 14)",
                 "Dissimilarity to\nfirst time-point",
                 "Dissimilarity (per 100 years)")[n],
          cex=0.8, las=0)
  
    axis(side=1, at=seq(6000,9000,500), labels=NA, tcl= -0.125)
    axis(side=1, at=seq(6000,9000,1000), labels=NA)  
  
  if(n == 3){
    par(xpd=NA)
    text(x=seq(6000,9000,1000),
         y=relative.axis.point(-0.05, "y"),
         labels=paste0(seq(6,9,1),",000"), 
         srt=30, adj=1)
    par(xpd=FALSE)
    }
  
  if(n == 1){mtext(side=3, line=0.2, text="Locality trends", font=2, cex=0.8)}
  
  # plot locality trajectories
  temp.loc <- loc.models[[n]]
  loc.preds <- temp.loc[[2]]
  temp.site <- site.models[[n]]
  site.preds <- temp.site[[2]]
  
  loc.preds$temp.col <- c("blue","red","darkgreen")[loc.preds$site]
  site.preds$temp.col <- c("blue","red","darkgreen")[site.preds$site] 
  ci.cols <- col2rgb(c("blue","red","darkgreen"))/255
  
  sapply(split(loc.preds, f=loc.preds$locality), function(x){
    
    lines(x$fit ~ x$raw.date, col=x$temp.col[1])
    
  })
  
  sapply(split(loc.preds, f=loc.preds$locality), function(x){
    
    text(x=x$raw.date[1],
         y=x$fit[1],
         pos=4, offset=0.25,
         labels=locality.number$number[locality.number$locality == x$locality[1]],
         col=x$temp.col[1])
  })
  
  text(x=relative.axis.point(0.025, "x"),
       y=relative.axis.point(0.93, "y"),
       labels=paste0("(", letters[c(1,4,7)[n]], ")"), font=2, adj=0)
  
  box()
  
  # SITE TRENDS ###
  plot(x=NULL, y=NULL, xlim=c(9500,5900), ylim=ylims, 
       axes=FALSE, yaxs="i", xlab="", ylab="")
  
  axis(side=1, at=seq(6000,9000,500), labels=NA, tcl= -0.125)
  axis(side=1, at=seq(6000,9000,1000), labels=NA)
  if(n == 3){
    par(xpd=NA)
    text(x=seq(6000,9000,1000),
         y=relative.axis.point(-0.05, "y"),
         labels=paste0(seq(6,9,1),",000"), 
         srt=30, adj=1)
    par(xpd=FALSE)
    
    mtext(side=1, line=2.15, text="Years before present", cex=0.8)
    
    
  }
  
  if(n == 1){mtext(side=3, line=0.2, text="Site trends", font=2, cex=0.8)}
  
  sapply(1:3, function(n){
   
    x <- split(site.preds, f=site.preds$site)[[n]]
    
    polygon(y=c(x$upper, rev(x$lower)),
            x=c(x$raw.date, rev(x$raw.date)),
            border=NA, 
            col=rgb(ci.cols[1,n], ci.cols[2,n], ci.cols[3,n], 0.2))
    
    lines(x$fit ~ x$raw.date, col=x$temp.col[1], lwd=1.5)
    
  })
  
  sapply(split(site.preds, f=site.preds$site), function(x){
    
    text(x=x$raw.date[1],
         y=x$fit[1],
         pos=4, offset=0.25,
         labels=substr(as.character(x$site[1]), 1, 1),
         col=x$temp.col[1])
  })
  box()
  text(x=relative.axis.point(0.025, "x"),
       y=relative.axis.point(0.93, "y"),
       labels=paste0("(", letters[c(2,5,8)[n]], ")"), font=2, adj=0)
  
  plot(x=NULL, y=NULL, xlim=c(9500,5900), ylim=ylims, 
       axes=FALSE, yaxs="i", xlab="", ylab="")
  box()
  
  axis(side=1, at=seq(6000,9000,500), labels=NA, tcl= -0.125)
  axis(side=1, at=seq(6000,9000,1000), labels=NA)
  if(n == 3){
    par(xpd=NA)
    text(x=seq(6000,9000,1000),
         y=relative.axis.point(-0.05, "y"),
         labels=paste0(seq(6,9,1),",000"), 
         srt=30, adj=1)
    par(xpd=FALSE)
    
  }
  
  
  # plot overall trajectory
  temp.all <- all.models[[n]]
  all.preds <- temp.all[[2]]
  
  if(n == 1){mtext(side=3, line=0.2, text="Region trends", font=2, cex=0.8)}
  
  polygon(x=c(all.preds$raw.date, rev(all.preds$raw.date)),
          y=c(all.preds$upper, rev(all.preds$lower)),
          col=rgb(0.5,0.5,0.5,0.2), border=NA)
  
  lines(all.preds$fit ~ all.preds$raw.date, lwd=2)
  
  text(x=relative.axis.point(0.025, "x"),
       y=relative.axis.point(0.93, "y"),
       labels=paste0("(", letters[c(3,6,9)[n]], ")"), font=2, adj=0)
  
  box()
  
})

dev.off()

#                         Relative diversity trends ####

pdf(date.wrap("./plots/relative diversity models", ".pdf"), 
    height=5.25, width=4.25, useDingbats = FALSE)

par(mfrow=c(3,2), mar=c(0,0,0,0), oma=c(3.5,5,2,1), 
    ps=10, mgp=c(3,0.5,0), tcl = -0.25, las=1)

loc.models <- list(genrich.loc.relative,
                   grrich.loc.relative,
                   genus.relative.beta.topdiss,
                   gr.relative.beta.topdiss,
                   genus.relative.beta.seqdiss,
                   gr.relative.beta.seqdiss)

all.models <- list(genrich.all.relative,
                   grrich.all.relative,
                   genus.relative.beta.topdiss.overall,
                   gr.relative.beta.topdiss.overall,
                   genus.relative.beta.seqdiss.overall,
                   gr.relative.beta.seqdiss.overall)

sapply(1:6, function(n){
  
  # set up plot
  if(n %in% c(1,2)){ylims=c(2,7)}
  if(n %in% c(3,4)){ylims=c(0,1.05)}
  if(n %in% c(5,6)){ylims=c(0,0.75)}
  
  plot(x=NULL, y=NULL, xlim=c(0,2200), ylim=ylims, 
       axes=FALSE, yaxs="i", xlab="", ylab="")
  
  
  if(n %in% c(1,3,5)){
    axis(side=2)
    mtext(side=2, line=2, 
          text=c("Rarefied richness (n = 14)",
                 "Dissimilarity to\nfirst time-point",
                 "Dissimilarity to\nprevious time-point")[(n + 1) / 2],
          cex=0.8, las=0)
  } else {axis(side=2, labels=NA)}
  
  
  axis(side=1, at=seq(0,2500,500), labels=NA, tcl= -0.125)
  axis(side=1, at=seq(0,2500,1000), labels=NA)
  
  if(n %in% c(5,6)){
    par(xpd=NA)
    text(x=seq(0,2500,1000),
         y=relative.axis.point(-0.05, "y"),
         labels=c(0, paste0(seq(1,2,1),",000")), 
         srt=30, adj=1)
    par(xpd=FALSE)
    
  }
  
  if(n==5){
    mtext(side=1, at=par("usr")[2], line=2.15,
          text="Years from time-series start", cex=0.8)
  }
  
  if(n %in% c(1:2)){
    mtext(side=3, line=0.2, text=c("Taxonomic (genera)", 
                                   "Functional groups (growth form)")[n], 
          font=2, cex=0.8)
  }
  
  # plot overall trajectory
  temp.all <- all.models[[n]]
  all.preds <- temp.all[[2]]
  
  polygon(x=c(all.preds$raw.date, rev(all.preds$raw.date)),
          y=c(all.preds$upper, rev(all.preds$lower)),
          col="grey85", border=NA)
  
  lines(all.preds$fit ~ all.preds$raw.date, lwd=2)
  
  # plot locality trajectories
  temp.loc <- loc.models[[n]]
  loc.preds <- temp.loc[[2]]
  
  loc.preds$temp.col <- c("blue","red","darkgreen")[loc.preds$site]
  
  sapply(split(loc.preds, f=loc.preds$locality), function(x){
    
    lines(x$fit ~ x$raw.date, col=x$temp.col[1])
    
  })
  
  sapply(split(loc.preds, f=loc.preds$locality), function(x){
    
    text(x=rev(x$raw.date)[1],
         y=rev(x$fit)[1],
         pos=4, offset=0.25,
         labels=locality.number$number[locality.number$locality == x$locality[1]],
         col=x$temp.col[1])
  })
  
  text(x=relative.axis.point(0.025, "x"),
       y=relative.axis.point(0.93, "y"),
       labels=paste0("(", letters[n], ")"), font=2, adj=0)
  
  box()
  
})

dev.off()

#                         compare genus and gr richness trends ####

pdf(date.wrap("./plots/genus & gr richness compared", ".pdf"), height=3.5, width=3.5,
    useDingbats = FALSE)
par(mar=c(3,3,1,1), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(x=NULL, y=NULL, xlim=c(3,6.5), ylim=c(2,6.5), xlab="", ylab="", asp=1,
     xaxt="n")

axis(side=1, at=3:7, mgp=c(3,0.25,0))
axis(side=1,at=seq(3,7,0.5), tcl=-0.125, labels=NA)
axis(side=2,at=seq(3,7,0.5), tcl=-0.125, labels=NA)

mtext(side=1, line=1.25, text="Growth form diversity")
mtext(side=2, line=1.5, text="Genus diversity", las=0)

all.top <- grrich.all.m[[2]]
all.seq <- genrich.all.m[[2]]

# polygon(x=c(all.top$upper, rev(all.top$lower)),
#         y=c(all.seq$upper, rev(all.seq$lower)),
#         border=NA, col="grey85")

library(shape)

loc.color <- data.frame(locality=unique(huon_coral$locality))
loc.color$col <- c("blue","red","darkgreen")[huon_tr$site[match(loc.color$locality,
                                                                huon_tr$locality)]]
# sapply(unique(huon_coral$locality), function(x){
# #   
# #   with(rich.df[rich.df$locality == x,],
# #        points(genus.rarefy ~ gr.rarefy,
# #               pch=16, col=loc.color$col[loc.color$locality==x],
# #               cex=0.4))
# #   
# # })

rect(xleft=par("usr")[1], xright=par("usr")[2],
     ybottom=par("usr")[3], ytop=par("usr")[4],
     border=NA, col=rgb(1,1,1,0.7))

Arrows(x0=all.top$fit[2],
       y0=all.seq$fit[2],
       x1=all.top$fit[1],
       y1=all.seq$fit[1],
       arr.type="triangle", arr.width=0.25, arr.length=0.25)

lines(all.seq$fit ~ all.top$fit, lwd=1.5)

sapply(unique(huon_coral$locality), function(x){
  
  top.preds <- grrich.loc.m[[2]]
  top.preds <- top.preds[top.preds$locality == x, ]
  seq.preds <- genrich.loc.m[[2]]
  seq.preds <- seq.preds[seq.preds$locality == x, ]
  
  lines(seq.preds$fit ~ top.preds$fit,
        col=loc.color$col[loc.color$locality==x])  
  
  Arrows(x0=top.preds$fit[2],
         y0=seq.preds$fit[2],
         x1=top.preds$fit[1],
         y1=seq.preds$fit[1],
         arr.type="triangle", arr.width=0.125, arr.length=0.125,
         col=loc.color$col[loc.color$locality==x])
  
  text(x=top.preds$fit[1],
       y=seq.preds$fit[1],
       pos=4, offset=0.25,
       labels=locality.number$number[locality.number$locality == x],
       col=loc.color$col[loc.color$locality==x],
       cex=0.8)
  
  
})

#abline(a=0,b=1)
dev.off()

#                         compare gen and gr top dissim trends ####

pdf(date.wrap("./plots/genus & gr top dissim compared", ".pdf"), height=3.5, width=3.5,
    useDingbats = FALSE)
par(mar=c(3,3,1,1), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", asp=1,
     xaxt="n")

axis(side=1,mgp=c(3,0.25,0))
axis(side=2)

mtext(side=1, line=1.25, text="Growth form similarity to oldest composition")
mtext(side=2, line=1.5, text="Genus diversity similarity to oldest composition", las=0)

all.top <- gr.top.beta.all.m[[2]]
all.seq <- genus.top.beta.all.m[[2]]

# polygon(x=c(all.top$upper, rev(all.top$lower)),
#         y=c(all.seq$upper, rev(all.seq$lower)),
#         border=NA, col="grey85")

library(shape)

loc.color <- data.frame(locality=unique(huon_coral$locality))
loc.color$col <- c("blue","red","darkgreen")[huon_tr$site[match(loc.color$locality,
                                                                huon_tr$locality)]]
# sapply(unique(huon_coral$locality), function(x){
# #   
# #   with(rich.df[rich.df$locality == x,],
# #        points(genus.rarefy ~ gr.rarefy,
# #               pch=16, col=loc.color$col[loc.color$locality==x],
# #               cex=0.4))
# #   
# # })

Arrows(x0=all.top$fit[2],
       y0=all.seq$fit[2],
       x1=all.top$fit[1],
       y1=all.seq$fit[1],
       arr.type="triangle", arr.width=0.25, arr.length=0.25)

lines(all.seq$fit ~ all.top$fit, lwd=1.5)

sapply(unique(huon_coral$locality), function(x){
  
  top.preds <- gr.top.beta.loc.m[[2]]
  top.preds <- top.preds[top.preds$locality == x, ]
  seq.preds <- genus.top.beta.loc.m[[2]]
  seq.preds <- seq.preds[seq.preds$locality == x, ]
  
  lines(seq.preds$fit ~ top.preds$fit,
        col=loc.color$col[loc.color$locality==x])  
  
  Arrows(x0=top.preds$fit[2],
         y0=seq.preds$fit[2],
         x1=top.preds$fit[1],
         y1=seq.preds$fit[1],
         arr.type="triangle", arr.width=0.125, arr.length=0.125,
         col=loc.color$col[loc.color$locality==x])
  
  text(x=top.preds$fit[1],
       y=seq.preds$fit[1],
       pos=4, offset=0.25,
       labels=locality.number$number[locality.number$locality == x],
       col=loc.color$col[loc.color$locality==x],
       cex=0.8)
  
  
})

#abline(a=0,b=1)
dev.off()


#                         compare gen and gr seq dissim trends ####

pdf(date.wrap("./plots/genus & gr top dissim compared", ".pdf"), height=3.5, width=3.5,
    useDingbats = FALSE)
par(mar=c(3,3,1,1), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", asp=1,
     xaxt="n")

axis(side=1,mgp=c(3,0.25,0))
axis(side=2)

mtext(side=1, line=1.25, text="Growth form similarity to previous composition")
mtext(side=2, line=1.5, text="Genus diversity similarity to previous composition", las=0)

all.top <- gr.region.beta.all.m[[2]]
all.seq <- genus.region.beta.all.m[[2]]

# polygon(x=c(all.top$upper, rev(all.top$lower)),
#         y=c(all.seq$upper, rev(all.seq$lower)),
#         border=NA, col="grey85")

library(shape)

loc.color <- data.frame(locality=unique(huon_coral$locality))
loc.color$col <- c("blue","red","darkgreen")[huon_tr$site[match(loc.color$locality,
                                                                huon_tr$locality)]]
# sapply(unique(huon_coral$locality), function(x){
# #   
# #   with(rich.df[rich.df$locality == x,],
# #        points(genus.rarefy ~ gr.rarefy,
# #               pch=16, col=loc.color$col[loc.color$locality==x],
# #               cex=0.4))
# #   
# # })

Arrows(x0=all.top$fit[2],
       y0=all.seq$fit[2],
       x1=all.top$fit[1],
       y1=all.seq$fit[1],
       arr.type="triangle", arr.width=0.25, arr.length=0.25)

lines(all.seq$fit ~ all.top$fit, lwd=1.5)

sapply(unique(huon_coral$locality), function(x){
  
  top.preds <- gr.top.beta.loc.m[[2]]
  top.preds <- top.preds[top.preds$locality == x, ]
  seq.preds <- genus.top.beta.loc.m[[2]]
  seq.preds <- seq.preds[seq.preds$locality == x, ]
  
  lines(seq.preds$fit ~ top.preds$fit,
        col=loc.color$col[loc.color$locality==x])  
  
  Arrows(x0=top.preds$fit[2],
         y0=seq.preds$fit[2],
         x1=top.preds$fit[1],
         y1=seq.preds$fit[1],
         arr.type="triangle", arr.width=0.125, arr.length=0.125,
         col=loc.color$col[loc.color$locality==x])
  
  text(x=top.preds$fit[1],
       y=seq.preds$fit[1],
       pos=4, offset=0.25,
       labels=locality.number$number[locality.number$locality == x],
       col=loc.color$col[loc.color$locality==x],
       cex=0.8)
  
  
})

#abline(a=0,b=1)
dev.off()

#                         axis rotation ####

library(plotrix)
library(shape)

pdf(date.wrap("./plots/diversity axis rotation plot", ".pdf"), height=7, width=6.5,
    useDingbats = FALSE)

par(mfrow=c(3,3), mar=c(2,0,1,0), oma=c(1,4,1,2), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

# LOCALITY PREDS 

richness.lims <- c(4.25, 7.75)

plot(x=NULL, y=NULL, xlim=richness.lims, ylim=c(-1,1.5), xlab="", ylab="", asp=1, xaxt="n")

axis(side=1, mgp=c(3,0.25,0))
mtext(side=2, line=2, text="Redundancy", las=0, cex=0.8)
mtext(side=3, line=0.2, text="Locality trends", font=2)

loc.color <- data.frame(locality=unique(huon_coral$locality))
loc.color$col <- c("blue","red","darkgreen")[huon_tr$site[match(loc.color$locality,
                                                                huon_tr$locality)]]

sapply(unique(huon_coral$locality), function(x){
  
  gr.pred <- grrich.loc.m[[2]]
  gr.pred <- gr.pred[gr.pred$locality == x, ]
  
  gen.pred <- genrich.loc.m[[2]]
  gen.pred <- gen.pred[gen.pred$locality == x, ]
  
  complexity <- (gr.pred$fit * cos(45 * (pi/180))) + (gen.pred$fit * sin(45 * (pi/180)))
  redundancy <- (-gr.pred$fit * sin(45 * (pi/180))) + (gen.pred$fit * cos(45 * (pi/180)))
  
  lines(redundancy ~ complexity,
        col=loc.color$col[loc.color$locality==x])  
  
  Arrows(x0=complexity[2],
         y0=redundancy[2],
         x1=complexity[1],
         y1=redundancy[1],
         arr.type="triangle", arr.width=0.125, arr.length=0.125,
         col=loc.color$col[loc.color$locality==x])
  
  text(x=complexity[1],
       y=redundancy[1],
       pos=4, offset=0.25,
       labels=locality.number$number[locality.number$locality == x],
       col=loc.color$col[loc.color$locality==x])
  
})

# SITE PREDS

plot(x=NULL, y=NULL, xlim=richness.lims, ylim=c(-1,1.5), xlab="", ylab="", asp=1,
     xaxt="n", yaxt="n")

axis(side=2, labels=NA)
axis(side=1, mgp=c(3,0.25,0))
mtext(side=1, line=1.25, text="Complexity", cex=0.8)
mtext(side=3, line=0.2, text="Site trends", font=2)

sapply(levels(huon_coral$site), function(x){
  
  gr.pred <- grrich.site.m[[2]]
  gr.pred <- gr.pred[gr.pred$site == x, ]
  
  gen.pred <- genrich.site.m[[2]]
  gen.pred <- gen.pred[gen.pred$site == x, ]
  
  complexity <- (gr.pred$fit * cos(45 * (pi/180))) + (gen.pred$fit * sin(45 * (pi/180)))
  redundancy <- (-gr.pred$fit * sin(45 * (pi/180))) + (gen.pred$fit * cos(45 * (pi/180)))
  
  site.col = c("blue", "red", "darkgreen")[levels(gen.pred$site) == x]
  
  lines(redundancy ~ complexity,
        col=site.col)  
  
  Arrows(x0=complexity[2],
         y0=redundancy[2],
         x1=complexity[1],
         y1=redundancy[1],
         arr.type="triangle", arr.width=0.125, arr.length=0.125,
         col=site.col)
  
  text(x=complexity[1],
       y=redundancy[1],
       pos=4, offset=0.25,
       labels=substr(x, 1, 1),
       col=site.col)
  
})

# REGION PREDS 

plot(x=NULL, y=NULL, xlim=richness.lims, ylim=c(-1,1.5), xlab="", ylab="", asp=1,
     xaxt="n", yaxt="n")

axis(side=2, labels=NA)
axis(side=1, mgp=c(3,0.25,0))

gr.pred <- grrich.all.m[[2]]
gen.pred <- genrich.all.m[[2]]

mtext(side=3, line=0.2, text="Region trends", font=2)

text(x=relative.axis.point(0.025, "x"),
     y=relative.axis.point(0.93, "y"),
     labels=paste0("(", letters[c(3,6,9)[n]], ")"), font=2, adj=0)

box()
complexity <- (gr.pred$fit * cos(45 * (pi/180))) + (gen.pred$fit * sin(45 * (pi/180)))
redundancy <- (-gr.pred$fit * sin(45 * (pi/180))) + (gen.pred$fit * cos(45 * (pi/180)))

comp.upper <- (gr.pred$upper * cos(45 * (pi/180))) + (gen.pred$upper * sin(45 * (pi/180)))
comp.lower <- (gr.pred$lower * cos(45 * (pi/180))) + (gen.pred$lower * sin(45 * (pi/180)))
redund.upper <- (-gr.pred$upper * sin(45 * (pi/180))) + (gen.pred$upper * cos(45 * (pi/180)))
redund.lower <- (-gr.pred$lower * sin(45 * (pi/180))) + (gen.pred$lower * cos(45 * (pi/180)))

points(redund.upper ~ comp.upper)
points(redund.lower ~ comp.lower)

Arrows(x0=complexity[2],
       y0=redundancy[2],
       x1=complexity[1],
       y1=redundancy[1],
       arr.type="triangle", arr.width=0.25, arr.length=0.25)

lines(redundancy ~ complexity, lwd=1.5)


mtext(side=4, line=0.4, las=0, text = "Richness", font=2)

# LOCALITY PREDS - top diss

plot(x=NULL, y=NULL, xlim=c(0.35,1.1), ylim=c(-0.375,0.35), xlab="", ylab="", asp=1,
     xaxt="n")

axis(side=1, mgp=c(3,0.25,0))
mtext(side=2, line=2, text="Conservation", las=0, cex=0.8)

loc.color <- data.frame(locality=unique(huon_coral$locality))
loc.color$col <- c("blue","red","darkgreen")[huon_tr$site[match(loc.color$locality,
                                                                huon_tr$locality)]]

sapply(unique(huon_coral$locality), function(x){
  
  gr.pred <- gr.top.beta.loc.m[[2]]
  gr.pred <- gr.pred[gr.pred$locality == x, ]
  
  gen.pred <- genus.top.beta.loc.m[[2]]
  gen.pred <- gen.pred[gen.pred$locality == x, ]
  
  complexity <- (gr.pred$fit * cos(45 * (pi/180))) + (gen.pred$fit * sin(45 * (pi/180)))
  redundancy <- (-gr.pred$fit * sin(45 * (pi/180))) + (gen.pred$fit * cos(45 * (pi/180)))
  
  lines(redundancy ~ complexity,
        col=loc.color$col[loc.color$locality==x])  
  
  Arrows(x0=complexity[2],
         y0=redundancy[2],
         x1=complexity[1],
         y1=redundancy[1],
         arr.type="triangle", arr.width=0.125, arr.length=0.125,
         col=loc.color$col[loc.color$locality==x])
  
  text(x=complexity[1],
       y=redundancy[1],
       pos=4, offset=0.25,
       labels=locality.number$number[locality.number$locality == x],
       col=loc.color$col[loc.color$locality==x])
  
})



# SITE PREDS

plot(x=NULL, y=NULL, xlim=c(0.35,1.1), ylim=c(-0.375,0.35), xlab="", ylab="", asp=1,
     xaxt="n", yaxt="n")

axis(side=2, labels=NA)
axis(side=1, mgp=c(3,0.25,0))
mtext(side=1, line=1.25, text="Dissimilarity", cex=0.8)

sapply(levels(huon_coral$site), function(x){
  
  gr.pred <- gr.top.beta.site.m[[2]]
  gr.pred <- gr.pred[gr.pred$site == x, ]
  
  gen.pred <- genus.top.beta.site.m[[2]]
  gen.pred <- gen.pred[gen.pred$site == x, ]
  
  complexity <- (gr.pred$fit * cos(45 * (pi/180))) + (gen.pred$fit * sin(45 * (pi/180)))
  redundancy <- (-gr.pred$fit * sin(45 * (pi/180))) + (gen.pred$fit * cos(45 * (pi/180)))
  
  site.col = c("blue", "red", "darkgreen")[levels(gen.pred$site) == x]
  
  lines(redundancy ~ complexity,
        col=site.col)  
  
  Arrows(x0=complexity[2],
         y0=redundancy[2],
         x1=complexity[1],
         y1=redundancy[1],
         arr.type="triangle", arr.width=0.125, arr.length=0.125,
         col=site.col)
  
  text(x=complexity[1],
       y=redundancy[1],
       pos=4, offset=0.25,
       labels=substr(x, 1, 1),
       col=site.col)
  
})

# REGION PREDS 

plot(x=NULL, y=NULL, xlim=c(0.35,1.1), ylim=c(-0.375,0.35), xlab="", ylab="", asp=1,
     xaxt="n", yaxt="n")

axis(side=2, labels=NA)
axis(side=1, mgp=c(3,0.25,0))

gr.pred <- gr.top.beta.all.m[[2]]
gen.pred <- genus.top.beta.all.m[[2]]

complexity <- (gr.pred$fit * cos(45 * (pi/180))) + (gen.pred$fit * sin(45 * (pi/180)))
redundancy <- (-gr.pred$fit * sin(45 * (pi/180))) + (gen.pred$fit * cos(45 * (pi/180)))

Arrows(x0=complexity[2],
       y0=redundancy[2],
       x1=complexity[1],
       y1=redundancy[1],
       arr.type="triangle", arr.width=0.25, arr.length=0.25)

lines(redundancy ~ complexity, lwd=1.5)
mtext(side=4, line=0.4, las=0, text = "Turnover (from first time-point)", font=2)

# LOCALITY PREDS - seq diss

plot(x=NULL, y=NULL, xlim=c(0.3,0.85), ylim=c(-0.15,0.25), xlab="", ylab="", asp=1,
     xaxt="n")

axis(side=1, mgp=c(3,0.25,0))
mtext(side=2, line=2, text="Persistence", las=0)

loc.color <- data.frame(locality=unique(huon_coral$locality))
loc.color$col <- c("blue","red","darkgreen")[huon_tr$site[match(loc.color$locality,
                                                                huon_tr$locality)]]

sapply(unique(huon_coral$locality), function(x){
  
  gr.pred <- gr.region.beta.loc.m[[2]]
  gr.pred <- gr.pred[gr.pred$locality == x, ]
  
  gen.pred <- genus.region.beta.loc.m[[2]]
  gen.pred <- gen.pred[gen.pred$locality == x, ]
  
  complexity <- (gr.pred$fit * cos(45 * (pi/180))) + (gen.pred$fit * sin(45 * (pi/180)))
  redundancy <- (-gr.pred$fit * sin(45 * (pi/180))) + (gen.pred$fit * cos(45 * (pi/180)))
  
  lines(redundancy ~ complexity,
        col=loc.color$col[loc.color$locality==x])  
  
  Arrows(x0=complexity[2],
         y0=redundancy[2],
         x1=complexity[1],
         y1=redundancy[1],
         arr.type="triangle", arr.width=0.125, arr.length=0.125,
         col=loc.color$col[loc.color$locality==x])
  
  text(x=complexity[1],
       y=redundancy[1],
       pos=4, offset=0.25,
       labels=locality.number$number[locality.number$locality == x],
       col=loc.color$col[loc.color$locality==x])
  
})

# SITE PREDS

plot(x=NULL, y=NULL, xlim=c(0.3,0.85), ylim=c(-0.15,0.25), xlab="", ylab="", asp=1,
     xaxt="n", yaxt="n")

axis(side=2, labels=NA)
axis(side=1, mgp=c(3,0.25,0))
mtext(side=1, line=1.5, text="Turnover", cex=0.8)

sapply(levels(huon_coral$site), function(x){
  
  gr.pred <- gr.region.beta.site.m[[2]]
  gr.pred <- gr.pred[gr.pred$site == x, ]
  
  gen.pred <- genus.region.beta.site.m[[2]]
  gen.pred <- gen.pred[gen.pred$site == x, ]
  
  complexity <- (gr.pred$fit * cos(45 * (pi/180))) + (gen.pred$fit * sin(45 * (pi/180)))
  redundancy <- (-gr.pred$fit * sin(45 * (pi/180))) + (gen.pred$fit * cos(45 * (pi/180)))
  
  site.col = c("blue", "red", "darkgreen")[levels(gen.pred$site) == x]
  
  lines(redundancy ~ complexity,
        col=site.col)  
  
  Arrows(x0=complexity[2],
         y0=redundancy[2],
         x1=complexity[1],
         y1=redundancy[1],
         arr.type="triangle", arr.width=0.125, arr.length=0.125,
         col=site.col)
  
  text(x=complexity[1],
       y=redundancy[1],
       pos=4, offset=0.25,
       labels=substr(x, 1, 1),
       col=site.col)
  
})

# REGION PREDS 

plot(x=NULL, y=NULL, xlim=c(0.3,0.85), ylim=c(-0.15,0.25), xlab="", ylab="", asp=1,
     xaxt="n", yaxt="n")

axis(side=2, labels=NA)
axis(side=1, mgp=c(3,0.25,0))

gr.pred <- gr.region.beta.all.m[[2]]
gen.pred <- genus.region.beta.all.m[[2]]

complexity <- (gr.pred$fit * cos(45 * (pi/180))) + (gen.pred$fit * sin(45 * (pi/180)))
redundancy <- (-gr.pred$fit * sin(45 * (pi/180))) + (gen.pred$fit * cos(45 * (pi/180)))

Arrows(x0=complexity[2],
       y0=redundancy[2],
       x1=complexity[1],
       y1=redundancy[1],
       arr.type="triangle", arr.width=0.25, arr.length=0.25)

lines(redundancy ~ complexity, lwd=1.5)
mtext(side=4, line=0.4, las=0, text = "Turnover (per 100 years)", font=2)

dev.off()

# ORDINATIONS ####
#                 Genus identity ####

genus.ord <- metaMDS(genusMatProp, index="jaccard")

genus.points <- as.data.frame(genus.ord$points)
genus.points$transect <- rownames(genus.points)
genus.points <- merge(genus.points, huon_site[,c("transect.age", "transect",
                                                 "locality", "site")],
                      all.x=TRUE, all.y=FALSE,
                      by.x="transect", by.y="transect")

col.grad <- data.frame(age=seq(round(min(huon_coral$transect.age),-1),
                               round(max(huon_coral$transect.age),-1),
                               10))

col.grad$Sialum <- colorRampPalette(c("white", "darkgreen"))(dim(col.grad)[1])
col.grad$Hubegong <- colorRampPalette(c("white", "blue"))(dim(col.grad)[1])
col.grad$Kanomi <- colorRampPalette(c("white", "red"))(dim(col.grad)[1])

common.gen <- colSums(ifelse(genusMatProp>0,1,0))/dim(genusMatProp)[1]
common.gen <- names(common.gen)[common.gen>0.4]

library(shape)

date.order <- c("Kwambu Platforms", "Kilasairo River SE", "Kilasairo River NW",
                "Pukau NW", "Paradise Springs", "Bonah River",
                "Midway Cove", "NW Point", "Loto Beach")

pdf(date.wrap("./plots/genus ords", ".pdf"), height=5, width=5)
par(mfrow=c(3,3), mar=c(0,0,0,0), oma=c(3,3,1,1), 
    ps=8, tcl=-0.25, las=1, mgp=c(3,0.5,0))

sapply(1:length(date.order), function(n){
  
  x <- date.order[n]
  
  sub.points <- genus.points[genus.points$locality==x,]
  sub.points <- sub.points[order(sub.points$transect.age, decreasing=TRUE),]
  
  sub.points$color <- as.character(col.grad[match(round(sub.points$transect.age,-1),
                                                  col.grad$age),
                                            as.character(sub.points$site[1])])
  
  sub.species <- genus.ord$species[rownames(genus.ord$species) %in% common.gen,]
  
  plot(genus.ord, type="n", xlim=c(-1,1.2), axes=FALSE, xlab="",
       ylab="")
  
  if(n %in% c(1,4,7)){
    axis(side=2)
    mtext(side=2, line=1.5, text="nMDS2", las=0, cex=0.8)
    } else{axis(side=2, labels=NA)}
  if(n %in% 7:9){
    axis(side=1, mgp=c(3,0,0))
    mtext(side=1, line=1, text="nMDS1", cex=0.8)
    } else{axis(side=1, labels=NA)}

  # other localities
  # sapply(unique(genus.points$locality)[unique(genus.points$locality) !=x],
  #        function(y){
  #          
  #          other.points <- genus.points[genus.points$locality==y,]
  #          other.points <- other.points[order(other.points$transect.age, decreasing=TRUE),]
  #          
  #          # 
  #          # sub.dim <- dim(other.points)[1]
  #          # Arrows(x0=other.points$MDS1[-sub.dim],
  #          #        x1=other.points$MDS1[-1],
  #          #        y0=other.points$MDS2[-sub.dim],
  #          #        y1=other.points$MDS2[-1],
  #          #        lwd=1, arr.type="triangle",
  #          #        arr.length=0.15, arr.width=0.15, arr.adj=1,
  #          #        col="grey70")
  #          
  #          points(other.points$MDS2 ~ other.points$MDS1,
  #                 pch=16, col="grey85", cex=0.75)
  #          
  #        })
  
  sapply(1:nrow(sub.species), function(n1){
  arrows(x0=0, y0=0,
         x1=sub.species[n1,1],
         y1=sub.species[n1,2],
         lwd=1.5, length=0.025, col="grey30")
  
  arrows(x0=0, y0=0,
         x1=sub.species[n1,1],
         y1=sub.species[n1,2],
         lwd=1.25, length=0.025, col="grey70")
  
  shadowtext(x=sub.species[n1,1],
             y=sub.species[n1,2],
             pos = c(2, 3, 2, 4, 4, 1, 4, 2)[n1], offset=0.25,
             labels=rownames(sub.species)[n1], font=1, col="grey70", bg="grey70",
             r=0)
  
  })
  
  sub.dim <- dim(sub.points)[1]
  Arrows(x0=sub.points$MDS1[-sub.dim],
         x1=sub.points$MDS1[-1],
         y0=sub.points$MDS2[-sub.dim],
         y1=sub.points$MDS2[-1],
         lwd=1, arr.type="triangle",
         arr.length=0.10, arr.width=0.10, arr.adj=1)
  
  points(sub.points$MDS2 ~ sub.points$MDS1,
         pch=21, bg=sub.points$color, lwd=0.5, cex=0.75)
  
  loc.number <- locality.number$number[locality.number$locality == x]
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(ifelse(loc.number > 6, 0.95, 0.925), "y"),
       labels = paste0(loc.number, ": ", x), 
       font=2, adj=0, cex=1.1,
       col=c("blue","red","darkgreen")[locality.number$site[locality.number$locality == x]])
box()  
})
dev.off()

#                 Growth form ####

gr.ord <- metaMDS(grMatProp, index="jaccard")

gr.points <- as.data.frame(gr.ord$points)
gr.points$transect <- rownames(gr.points)
gr.points <- merge(gr.points,  huon_site[,c("transect.age", "transect",
                                            "locality", "site")],
                   all.x=TRUE, all.y=FALSE,
                   by.x="transect", by.y="transect")

col.grad <- data.frame(age=seq(round(min(huon_coral$transect.age),-1),
                               round(max(huon_coral$transect.age),-1),
                               10))
col.grad$Sialum <- colorRampPalette(c("white", "darkgreen"))(dim(col.grad)[1])
col.grad$Hubegong <- colorRampPalette(c("white", "blue"))(dim(col.grad)[1])
col.grad$Kanomi <- colorRampPalette(c("white", "red"))(dim(col.grad)[1])

gr.points$color <- col.grad$col[match(round(gr.points$transect.age,-1),
                                      col.grad$age)]

common.gen <- colSums(ifelse(grMatProp>0,1,0))/dim(grMatProp)[1]
common.gen <- names(common.gen)[common.gen>0.4]
common.gen <- colnames(grMatProp)
gr.names <- c("Branch closed", "Branch open", "Columnar", "Corymbose", "Digitate",
              "Encrusting", "Laminar", "Massive", "Solitary", "Submassive", "Table/plate")

library(shape)

pdf(date.wrap("./plots/growth form ords", ".pdf"), height=5, width=5)
par(mfrow=c(3,3), mar=c(0,0,0,0), oma=c(3,3,1,1), 
    ps=8, tcl=-0.25, las=1, mgp=c(3,0.5,0))

sapply(1:length(date.order), function(n){
  
  x <- date.order[n]
  
  sub.points <- gr.points[gr.points$locality==x,]
  sub.points <- sub.points[order(sub.points$transect.age, decreasing=TRUE),]
  
  sub.points$color <- as.character(col.grad[match(round(sub.points$transect.age,-1),
                                                  col.grad$age),
                                            as.character(sub.points$site[1])])
  
  sub.species <- gr.ord$species[rownames(gr.ord$species) %in% common.gen,]
  
  plot(gr.ord, type="n", xlim=c(-1,1.4), axes=FALSE, xlab="",
       ylab="")
  
  if(n %in% c(1,4,7)){
    axis(side=2)
    mtext(side=2, line=1.5, text="nMDS2", las=0, cex=0.8)
  } else{axis(side=2, labels=NA)}
  if(n %in% 7:9){
    axis(side=1, mgp=c(3,0,0))
    mtext(side=1, line=1, text="nMDS1", cex=0.8)
  } else{axis(side=1, labels=NA)}
  
  # other localities
  # sapply(unique(genus.points$locality)[unique(genus.points$locality) !=x],
  #        function(y){
  #          
  #          other.points <- genus.points[genus.points$locality==y,]
  #          other.points <- other.points[order(other.points$transect.age, decreasing=TRUE),]
  #          
  #          # 
  #          # sub.dim <- dim(other.points)[1]
  #          # Arrows(x0=other.points$MDS1[-sub.dim],
  #          #        x1=other.points$MDS1[-1],
  #          #        y0=other.points$MDS2[-sub.dim],
  #          #        y1=other.points$MDS2[-1],
  #          #        lwd=1, arr.type="triangle",
  #          #        arr.length=0.15, arr.width=0.15, arr.adj=1,
  #          #        col="grey70")
  #          
  #          points(other.points$MDS2 ~ other.points$MDS1,
  #                 pch=16, col="grey85", cex=0.75)
  #          
  #        })
  
  sapply(1:nrow(sub.species), function(n1){
    arrows(x0=0, y0=0,
           x1=sub.species[n1,1],
           y1=sub.species[n1,2],
           lwd=1.5, length=0.025, col="grey30")
    
    arrows(x0=0, y0=0,
           x1=sub.species[n1,1],
           y1=sub.species[n1,2],
           lwd=1.25, length=0.025, col="grey70")
    
    shadowtext(x=sub.species[n1,1],
               y=sub.species[n1,2],
               pos = c(1, 2, 3, 2, 4, 1, 1, 3, 2, 4, 2)[n1], offset=0.25,
               labels=gr.names[n1], font=1, col="grey70", bg="grey70",
               r=0)
    
  })
  
  sub.dim <- dim(sub.points)[1]
  Arrows(x0=sub.points$MDS1[-sub.dim],
         x1=sub.points$MDS1[-1],
         y0=sub.points$MDS2[-sub.dim],
         y1=sub.points$MDS2[-1],
         lwd=1, arr.type="triangle",
         arr.length=0.10, arr.width=0.10, arr.adj=1)
  
  points(sub.points$MDS2 ~ sub.points$MDS1,
         pch=21, bg=sub.points$color, lwd=0.5, cex=0.75)
  
  loc.number <- locality.number$number[locality.number$locality == x]
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(ifelse(loc.number > 6, 0.95, 0.925), "y"),
       labels = paste0(loc.number, ": ", x), 
       font=2, adj=0, cex=1.1,
       col=c("blue","red","darkgreen")[locality.number$site[locality.number$locality == x]])
  box()  
})
dev.off()



# test ####

col.ramp <- data.frame(score=seq(-1.25,1.25,0.01))
col.ramp$col <- colorRampPalette(c("white","black"))(nrow(col.ramp))

plot(rich.df$genus.rarefy ~ rich.df$gr.rarefy, pch=21,
     bg=col.ramp$col[match(round(rich.df$redundancy, 2), round(col.ramp$score,2))])

# primary axes are equivalent to the hypotenuse of the right-angled triangle
# formed by the taxonomic and functional ratios
plot(sqrt(rich.df$genus.rarefy^2 + rich.df$gr.rarefy^2) ~ rich.df$complexity)
abline(0,1, lty="31", col="grey80")  

# secondary axes are roughly equivalent to ratios
plot(rich.df$genus.rarefy / rich.df$gr.rarefy ~ rich.df$redundancy, axes=FALSE)
axis(side=1, at= seq(-1,2,0.1))
axis(side=2, at= seq(-1,2,0.1))
box()
abline(v=0, lty="31", col="grey80")  
abline(h=1, lty="31", col="grey80")
abline(h=1/1.2, lty="31", col="red")
abline(v=-0.5, lty="31", col="red")




# how long to add 1 genera & growth forms in fastest trend ####
test <- alpha.model.list[[2]][[2]][[3]]
test <- test[test$locality == "Pukau NW",]
lims <- range(huon_tr$transect.age[huon_tr$locality == "Pukau NW"])
test <- test[test$pred.date >= lims[1] & test$pred.date <= lims[2],]

# now calculate rate of change per 10 years
testDist <- diff(test$fit.1[order(test$pred.date, decreasing=TRUE)])

# collapses to linear fit, slope of 0.02047096 per 10 years
sqrt(2) / (0.02047096 / 10)

# ####

# Site 4 trends
site4Tr <- huon_site[huon_site$locality == "Pukau NW",]
site4Tr <- site4Tr[order(site4Tr$transect.age, decreasing=TRUE),]

site4Gen <- genusMatProp[match(site4Tr$transect, rownames(genusMatProp)),]
site4Gen <- site4Gen[,colSums(site4Gen) > 0.05 ]

site4Gr <- grMatProp[match(site4Tr$transect, rownames(grMatProp)),]
site4Gr <- site4Gr[,colSums(site4Gr) > 0.05 ]

# site 4 gen-int tables
site4Tr36 <- table(huon_coral$genus[huon_coral$transect == 36], 
      huon_coral$growth.form[huon_coral$transect == 36])
site4Tr36 <- site4Tr36[rowSums(site4Tr36) > 0, ]
site4Tr36 <- site4Tr36[,colSums(site4Tr36)>0]

site4Tr26 <- table(huon_coral$genus[huon_coral$transect == 26], 
                   huon_coral$growth.form[huon_coral$transect == 26])
site4Tr26 <- site4Tr26[rowSums(site4Tr26) > 0, ]
site4Tr26 <- site4Tr26[,colSums(site4Tr26)>0]



# site 1 trends
site1Tr <- huon_site[huon_site$locality == "Midway Cove",]
site1Tr <- site1Tr[order(site1Tr$transect.age, decreasing=TRUE),]

site1Gen <- genusMatProp[match(site1Tr$transect, rownames(genusMatProp)),]
site1Gen <- site1Gen[,colSums(site1Gen) > 0.05 ]

site1Gr <- grMatProp[match(site1Tr$transect, rownames(grMatProp)),]
site1Gr <- site1Gr[,colSums(site1Gr) > 0.05 ]

rowSums(site1Gen > 0)
rowSums(site1Gr > 0)

# site 4 gen-int tables
site1Tr64 <- table(huon_coral$genus[huon_coral$transect == 64], 
                   huon_coral$growth.form[huon_coral$transect == 64])
site1Tr64 <- site1Tr64[rowSums(site1Tr64) > 0, ]
site1Tr64 <- site1Tr64[,colSums(site1Tr64)>0]

site1Tr52 <- table(huon_coral$genus[huon_coral$transect == 52], 
                   huon_coral$growth.form[huon_coral$transect == 52])
site1Tr52 <- site1Tr52[rowSums(site1Tr52) > 0, ]
site1Tr52 <- site1Tr52[,colSums(site1Tr52)>0]

# genus vs growth plot for presentation ####

pdf("./plots/richnessComp.pdf", height=4.5, width=4.5, useDingbats = FALSE)
par(mar=c(3,3,1,1), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(turnoverDf$top.gen.diss ~ turnoverDf$top.gr.diss, asp=1,
     xlim=c(0,1), ylim=c(0,1), pch=16, col="grey50", xaxt="n", xlab="", ylab="")

axis(side=1, mgp=c(3,0.2,0))
mtext(side=2, line=2, text="Genus turnover", las=0)
mtext(side=1, line=1.5, text="Growth form turnover")
#abline(a=0,b=1, lty="31", lwd=2)
dev.off()

pdf("./plots/richnessComp2.pdf", height=4.5, width=4.5, useDingbats = FALSE)
par(mar=c(3,3,1,1), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
plot(turnoverDf$conservation ~ turnoverDf$dissimilarity, asp=1, cex=1,
     xlim=c(0.2,1.1), ylim=c(-0.4,0.4), pch=16, col="grey50", xaxt="n", xlab="", ylab="")
axis(side=1, mgp=c(3,0.2,0))
mtext(side=2, line=2, text="Offset (or ratio) of turnover", las=0)
mtext(side=1, line=1.5, text="Coordinated turnover")
abline(h=0, lty="31", lwd=2)
dev.off()


# Examine specific sites ####

temp <- huon_site[huon_site$locality == "Midway Cove",]
temp <- temp[order(temp$transect.age, decreasing=TRUE),]

tempGen <- genusMatProp[match(temp$transect, rownames(genusMatProp)),]
tempGen[,colSums(tempGen) > 0]

tempGr <- grMatProp[match(temp$transect, rownames(genusMatProp)),]
tempGr[,colSums(tempGr) > 0]
