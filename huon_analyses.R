## ################################################# ####
# Diversity dyanmics in Huon Holocence transect data ####
# Author: Timothy L Staples                          ####
# ################################################## ####
# WORKING DIRECTORY ####

rm(list = ls())
setwd("/Users/uqtstapl/Dropbox/Tim/Post-doc/Research projects/Upload_projects/coral_diversity_dynamics")
#setwd("PATH TO THIS FILE")

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

package.loader(c("vegan", "lme4", "MASS", "DHARMa", "shape",
                 "circular", "maptools", "sp", "rgdal", "rgeos", "fields",
                 "grImport2", "mgcv", "plotrix", "iNEXT", "performance",
                 "TeachingDemos", "tidyr", "raster"))

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


# there are a number of field/taxonomic growth form assessments that
# resulted in multiple Veron growth forms. We need to decide how to
# assess these.


huon_coral$binom = ifelse(!is.na(huon_coral$species),
                          paste0(huon_coral$genus, " ", huon_coral$species),
                          NA) 
spTable <- as.data.frame(table(huon_coral$binom))
spTable = spTable[order(spTable$Freq, decreasing=TRUE),]
write.csv(spTable, "./raw.datafiles/spTable.csv")

# Step 1 - get the synonymized growth form assessments.
# Step 2 - if multiple, collapse to most likely of the two for the intercept genus
#          (E.g., if Acropora intercept listed as "massive" and "branching_open", record
#           as "branching_open")
gr.cats <- read.csv("./raw.datafiles/growth.grep.csv")

growth.split <- strsplit(as.character(huon_coral$growth.comb), ":")
growth.split <- do.call("rbind",lapply(growth.split, function(x){
  if(length(x)==2){return(x)}
  return(c(x,""))
}))
growth.mat <- sapply(1:nrow(gr.cats), function(n){
  
  posArg = gr.cats$pos.arg[n]
  negArg = gr.cats$neg.arg[n]
  
  growthCount = grepl(posArg,growth.split[,1]) +
                grepl(posArg, growth.split[,2])
  
  if(!is.na(negArg)){
    growthCount = growthCount - grepl(negArg, growth.split[,1]) -
                                grepl(negArg, growth.split[,2])
  }
  return(growthCount)
})
colnames(growth.mat) = gr.cats$growth.form

singles = ifelse(growth.mat >0 , 1, 0) + matrix(rep(rowSums(growth.mat>0)==1, each=ncol(growth.mat)), ncol=ncol(growth.mat), byrow=TRUE)
singles = which(singles == 2, arr.ind=TRUE)
singles = singles[order(singles[,1]),]

growth.form <- rep(NA, nrow(growth.mat))
growth.form[singles[,1]] = colnames(growth.mat)[singles[,2]]

# for multiples, if ET/ID (coral taxonomists) offer an ID that disagrees with the
# field ID, we go with the actual ID
EDID.favoured <- which((growth.mat>1)==1 & rowSums(growth.mat>1, na.rm=TRUE)==1, arr.ind=TRUE)
EDID.favoured <- EDID.favoured[!EDID.favoured[,1] %in% singles[,1],]
growth.form[EDID.favoured[,1]] = colnames(growth.mat)[EDID.favoured[,2]]

# next are ones with no growth form ID ('unknown' or similar for field and ET/ID)
unknowns <- which(rowSums(growth.mat)==0, arr.ind=TRUE)
growth.form[unknowns] = "unknown"

# next are ones where we have several growth form options, both of which
# are equally weighted. Here is where we need to start looking at taxonomic
# likelihood

# we'll use the known growth forms by genus
gen.gr.mat <- table(huon_coral$genus, growth.form)
genus.dominant <- data.frame(genus = rownames(gen.gr.mat),
                             growth.form = colnames(gen.gr.mat)[apply(gen.gr.mat, 1, which.max)])

# check if the dominant form is one of the max growth form scores
gen.dom.long <- data.frame(growth.form = as.character(genus.dominant$growth.form[match(huon_coral$genus,
                                                                                       genus.dominant$genus)]),
                           stringsAsFactors = FALSE)
gen.dom.long$col.match <- match(gen.dom.long$growth.form, colnames(growth.mat))

gen.dom.long$growth.score <- NA
gen.dom.long$growth.score[!is.na(gen.dom.long$col.match)] = growth.mat[cbind(which(!is.na(gen.dom.long$col.match)),
                                                                             gen.dom.long$col.match[!is.na(gen.dom.long$col.match)])]

gen.dom.long$needed <- is.na(growth.form)
gen.dom.long$row.max <- apply(growth.mat, 1, max)
gen.dom.long$valid <- gen.dom.long$needed & (gen.dom.long$growth.score == gen.dom.long$row.max)
gen.dom.long$valid[is.na(gen.dom.long$valid)] = FALSE

growth.form[gen.dom.long$valid] = gen.dom.long$growth.form[gen.dom.long$valid]

growth.form[is.na(growth.form)] = "unknown"
table(growth.form)

huon_coral$growth.form <- growth.form
huon_coral <- huon_coral[,!colnames(huon_coral) %in% gr.cats$growth.form]

# add in additional species-level traits that are unknown from the intercept
speciesGr <- read.csv("./raw.datafiles/spTableFilled.csv")
spRows <- which(huon_coral$growth.form == "unknown" & !is.na(huon_coral$binom))

growthMatch <- speciesGr$Growth[match(huon_coral$binom[spRows], speciesGr$Var1)]
huon_coral$growth.form[spRows] = growthMatch
huon_coral$growth.form <- as.factor(huon_coral$growth.form)

huon_coral_all <- huon_coral
huon_coral <- droplevels(huon_coral[huon_coral$growth.form != "unknown",])


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
                                method="bray")

# ANALYSES ####
#                     Axis correlations ####
#                       Compare tax and fun alpha richness ####
H0cor <- lmer(grH0 ~ genusH0 + (1|site/locality), data=richDf)
r2(H0cor)
summary(H0cor)
H0corC <- summary(H0cor)$coefficients
cbind(H0corC[,1] - 1.96 * H0corC[,2],
      H0corC[,1] + 1.96 * H0corC[,2])

H1cor <- lmer(grH1 ~ genusH1 + (1|site), data=richDf)
r2(H1cor)
summary(H1cor)
H1corC <- summary(H1cor)$coefficients
cbind(H1corC[,1] - 1.96 * H1corC[,2],
      H1corC[,1] + 1.96 * H1corC[,2])

H2cor <- lmer(grH2 ~ genusH2 + (1|site), data=richDf)
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

#                            Directional versus sequential ####
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

#                                                     PLOT ####

loc.color <- data.frame(locality=unique(huon_coral$locality))
loc.color$col <- c("blue","red","darkgreen")[huon_tr$site[match(loc.color$locality,
                                                                huon_tr$locality)]]

mList = list(H0cor, H1cor, H2cor, dirBcor, seqBcor, taxSeqDir, funSeqDir)
xList = c("genusH0", "genusH1", "genusH2", "top.gen.diss", "seq.gen.diss", "seq.gen.diss", "seq.gr.diss")
yList = c("grH0", "grH1", "grH2", "top.gr.diss", "seq.gr.diss", "top.gen.diss", "top.gr.diss")
xLabs = c(expression("Tax. "*alpha*" diversity"),
          expression("Tax. "*alpha*" diversity"),
          expression("Tax. "*alpha*" diversity"),
          expression("Tax. "*beta*" diversity"),
          expression("Tax. "*beta*" diversity"),
          expression("Directional "*beta),
          expression("Directional "*beta))

yLabs = c(expression("Fun. "*alpha*" diversity"),
          expression("Fun. "*alpha*" diversity"),
          expression("Fun. "*alpha*" diversity"),
          expression("Fun. "*beta*" diversity"),
          expression("Fun. "*beta*" diversity"),
          expression("Sequential "*beta),
          expression("Sequential "*beta))

pdf("./plots/correlationModels.pdf", height=4.5, width=5)
par(mfrow=c(3,3), mar=c(2,2,1,1), oma=c(1,1,0,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)

for(i in 1:length(mList)){
  
  print(i)
  if(i < 4){predData = richDf} else {predData = turnoverDf}
  xT = xList[i]
  yT = yList[i]
  mT = mList[[i]]
  
  predDf = data.frame(X = seq(min(predData[, xT], na.rm=TRUE), max(predData[,xT], na.rm=TRUE), len=200))
  colnames(predDf) = xT
  predDf = cbind(predDf, mer.ci(mT, predDf, 999, 4))
  
  plotRange = range(c(0, predData[,xT], predData[,yT], predDf$lower, predDf$upper), na.rm=TRUE)
  plotRange = plotRange + c(0, 0.05 * diff(plotRange))
  plot(NULL, xlim=plotRange, ylim=plotRange,
       xaxt="n", yaxs="i", xaxs="i", xlab="", ylab="")

  axis(side=1, mgp=c(3,0.2,0))
  mtext(side=1, line=1.5, text=xLabs[i], cex=0.8)
  mtext(side=2, line=1.5, text=yLabs[i], cex=0.8, las=0)
  
  points(jitter(predData[,yT], amount=ifelse(i==1, 0.2,0)) ~ predData[,xT], pch=16,
         col=loc.color$col[match(predData$locality, loc.color$locality)])
  
  polygon(y=c(predDf$lower, rev(predDf$upper)),
           x=c(predDf[,1], rev(predDf[,1])),
           border=NA, col=rgb(0.5,0.5,0.5,0.5))
  lines(x=predDf[,1], y=predDf$fit, lwd=3)
  
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.925, "y"),
       labels=paste0("(",LETTERS[i],")"), adj=0, font=2)
  
  text(x=relative.axis.point(0.15, "x"),
       y=relative.axis.point(0.925, "y"),
       labels=c("Hill q = 0", "Hill q = 1", "Hill q = 2", "Directional", "Sequential",
                "Taxonomic", "Functional")[i], adj=0)
  
  mr2 = performance(mT)
  text(x=relative.axis.point(0.97, "x"),
       y=relative.axis.point(0.06, "y"),
       labels=round(mr2$R2_marginal, 3), adj=1)
  
  text(x=relative.axis.point(0.775, "x"),
       y=relative.axis.point(0.08, "y"),
       labels=expression("R"^2*" = "), adj=1)
  
  
  if(i == 5){plot.new()}
  }

dev.off()

#               Multivariate GAM models ####

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

# model validation
pdf("./plots/model diagnostics.pdf", height=9, width=9, useDingbats = FALSE)
par(mfrow=c(5,4), mar=c(0,3,0,2), oma=c(3,2,2,1), ps=10, las=1, tcl=-0.25)
sapply(1:5, function(n){

x=list(alpha.model.list[[1]]$models[[1]],
       alpha.model.list[[2]]$models[[1]],
       alpha.model.list[[3]]$models[[1]],
       beta.model.list[[1]]$models[[1]],
       beta.model.list[[2]]$models[[1]])[[n]]
xData = x$model

qqnorm(residuals(x)[,1], main="", xlab="", ylab="", xaxt="n")
if(n < 5){axis(side=1, labels=NA)} else {axis(side=1, mgp=c(3,0.5,0)); mtext(side=1, line=2, text="Theoretical quantiles")}
qqline(residuals(x)[,1])

if(n==1){mtext(side=3,  at=par('usr')[2], adj=0.5, line=0.1, text="Taxonomic diversity")}
if(n==3){mtext(side=2, adj=0.5, line=2, text="Sample quantiles", las=0)}

mtext(side=2, line=3, text=c(expression("Hill 0 "*alpha), expression("Hill 1 "*alpha), expression("Hill 2 "*alpha), expression("Directional "*beta), expression("Sequential "*beta))[n], las=0, font=2, cex=1.25)

acf(residuals(x)[,1][order(xData$pred.date, decreasing=TRUE)], xaxt="n")
if(n==3){mtext(side=2, adj=0.5, line=2, text="Autocorrelation", las=0)}
if(n < 5){axis(side=1, labels=NA)} else {axis(side=1, mgp=c(3,0.5,0)); mtext(side=1, line=2, text="Lag")}

qqnorm(residuals(x)[,2], main="", xlab="", ylab="", xaxt="n")
if(n < 5){axis(side=1, labels=NA)} else {axis(side=1, mgp=c(3,0.5,0)); mtext(side=1, line=2, text="Theoretical quantiles")}
qqline(residuals(x)[,2])

if(n==1){mtext(side=3, at=par('usr')[2], adj=0.5, line=0.1, text="Functional diversity")}
if(n==3){mtext(side=2, adj=0.5, line=2, text="Sample quantiles", las=0)}

acf(residuals(x)[,2][order(xData$pred.date, decreasing=TRUE)], xaxt="n")
if(n==3){mtext(side=2, adj=0.5, line=2, text="Autocorrelation", las=0)}
if(n < 5){axis(side=1, labels=NA)} else {axis(side=1, mgp=c(3,0.5,0)); mtext(side=1, line=2, text="Lag")}

})
dev.off()

# NULL MODELS ####
#                      Markov Simulation ####    

simIter = 999

# generate simulated turnover data for each site
print("Generating simulations")
markovSimList <- lapply(unique(huon_site$locality), function(loc){
  
  print(as.character(loc))
  
  # set up intercept data with probabilities of occurrence
  loc_sub <- droplevels(huon_site[huon_site$locality == loc,])
  loc_int <- droplevels(huon_coral[huon_coral$locality == loc,])
  loc_int <- loc_int[complete.cases(loc_int[,c("int.delta", "genus", "growth.form")]),]
  loc_int$gen.gr = paste0(loc_int$genus, ":", loc_int$growth.form)
  
  spTrMat = tapply(loc_int$int.delta,
                   list(loc_int$gen.gr, as.factor(loc_int$transect)),
                   sum, na.rm=TRUE, simplify = FALSE)
  spTrMat = matrix(sapply(spTrMat, function(x){ifelse(is.null(x), 0, x)}), nrow=dim(spTrMat)[1], ncol=dim(spTrMat)[2],
                   dimnames=list(rownames(spTrMat), NULL))
  spTrProp = prop.table(spTrMat,2)
  spTrDiff = cbind(NA, t(apply(spTrProp, 1, diff)))
  spIdentity = do.call("rbind", strsplit(rownames(spTrProp), ":"))
  sampThresh = min(spTrProp[spTrProp > 0])
  
  # this is the likelihood of being selected to re-emergence into transect
  spMeans = rowMeans(spTrProp)
  
  # we also need the relationship between relative abundance and change in relative abundance
  spDiffDf = cbind(ra = spTrProp[spTrProp != 0 | spTrDiff != 0],
                   dra = spTrDiff[spTrProp != 0 | spTrDiff != 0])
  spDiffDf = as.data.frame(spDiffDf[complete.cases(spDiffDf),])
  diffM = lm(abs(dra) ~ ra, data=spDiffDf)
  
  # return models
  return(diffM)
  
  # compositions recorded as list
  locSims = lapply(1:simIter, function(n){
  
    simpComp = list(t1 = data.frame(locality = loc,
                                    transect = 1,
                                    genus = spIdentity[,1],
                                    growth.form = spIdentity[,2],
                                    prop = spTrProp[,1]))
    
  for(tr in 2:length(unique(loc_int$transect))){
    
    print(tr)
    # begin with past transect
    simpComp[[tr]] = simpComp[[tr-1]]
    simpComp[[tr]]$transect = tr
    
    # draw variance and apply as Gaussian noise
    betaVars <- predict(diffM, newdata=data.frame(ra=simpComp[[tr]]$prop), se.fit=TRUE)
    simVar = rnorm(length(betaVars$fit), betaVars$fit, betaVars$se.fit) * sample(c(-1,1), length(betaVars$fit), replace=TRUE)
            
    simpComp[[tr]]$prop = simpComp[[tr]]$prop + simVar
    simpComp[[tr]]$prop = simpComp[[tr]]$prop / sum(simpComp[[tr]]$prop)
    #simpComp[[tr]]$prop[simpComp[[tr]]$prop < 0.0175] = 0
    #simpComp[[tr]]$prop = simpComp[[tr]]$prop / sum(simpComp[[tr]]$prop)
}
  
  return(do.call("rbind", simpComp))
  })
  
  return(locSims)
})
  
# separate sims so we have data for each locality
print("Data processing and diversity calculation")
markovSimDfs <- lapply(1:simIter, function(n){
      do.call("rbind", lapply(markovSimList, function(x){x[[n]]}))
})

# Next: aggregate sims into gen/gr mats like main analysis
markovSimMats <- lapply(markovSimDfs, function(y){
  
  y$locTran <- paste0(y$locality, y$transect)
  y = y[y$prop > 0,]
  y$genus <- as.factor(y$genus)
  y$growth.form <- as.factor(y$growth.form)
  # table genus
  subGenMat <- do.call("rbind", 
                          lapply(split(y, f= ~y$transect + y$locality), 
                                 function(x){
                                   temp <- tapply(x$prop, x$genus, sum, na.rm=TRUE)
                                   ifelse(is.na(temp), 0, temp)
                                 }))
  subGenMat <- prop.table(subGenMat, 1)
  
  subGrMat <- do.call("rbind", 
                       lapply(split(y, f= ~y$transect + y$locality), 
                              function(x){
                                temp <- tapply(x$prop, x$growth.form, sum, na.rm=TRUE)
                                ifelse(is.na(temp), 0, temp)
                              }))
  subGrMat <- prop.table(subGrMat, 1)
  
  return(list(subGenMat, subGrMat))
  
})

# caclulate alpha and beta matrices
markovSimDivs <- lapply(markovSimMats, function(x){
  
  genMat = x[[1]][complete.cases(x[[1]]),]
  grMat = x[[2]][complete.cases(x[[2]]),]
  
  # just need to superimpose transect numbers running oldest to youngest
  tranNums <- do.call("rbind", lapply(split(huon_site, f=huon_site$locality), function(x){
    x <- x[order(x$transect.age, decreasing=TRUE),]
    x$locTran <- paste0(x$locality, 1:nrow(x))
    return(x)
  }))
  
  rownames(genMat) = tranNums$transect
  rownames(grMat) = tranNums$transect
  
  # alpha Mat
 alphaMat = alphaDiversityCalc(siteData = huon_site,
                     genusMatProp = genMat,
                     grMatProp = grMat)
 
 # beta Mat
 betaMat = betaDiversityCalc(siteData = huon_site,
                               genusMatProp = genMat,
                               grMatProp = grMat,
                             method = "bray")
 
 return(list(alphaMat, betaMat))
  
})

# run models
print("Running diversity models")
markovModels <- lapply(1:length(markovSimDivs), function(n){
  print(n)
  x=markovSimDivs[[n]]
  alphaMs = mvGamsFreq(dataList = list(x[[1]], x[[1]], x[[1]]),
                               primaryvarList = c("genusH0", "genusH1", "genusH2"),
                               secondaryvarList = c("grH0", "grH1", "grH2"),
                               dateDiff=FALSE)
  
  betaMs =  mvGamsFreq(dataList = list(x[[2]], x[[2]]),
                        primaryvarList = c("top.gen.diss", "seq.gen.diss"),
                        secondaryvarList = c("top.gr.diss", "seq.gr.diss"),
                        dateDiff=TRUE)
  
  return(list(alphaMs, betaMs))
  
})

#                     Matrix permutation ####

# generate matrices for permutations
print("Generating simulations")
permuteSimList <- lapply(unique(huon_site$locality), function(loc){
  
  print(as.character(loc))
  
  # set up intercept data with probabilities of occurrence
  loc_sub <- droplevels(huon_site[huon_site$locality == loc,])
  loc_int <- droplevels(huon_coral[huon_coral$locality == loc,])
  loc_int <- loc_int[complete.cases(loc_int[,c("int.delta", "genus", "growth.form")]),]
  loc_int$gen.gr = paste0(loc_int$genus, ":", loc_int$growth.form)
  
  spTrMat = tapply(loc_int$int.delta,
                   list(loc_int$gen.gr, as.factor(loc_int$transect)),
                   sum, na.rm=TRUE, simplify = FALSE)
  spTrMat = matrix(sapply(spTrMat, function(x){ifelse(is.null(x), 0, x)}), nrow=dim(spTrMat)[1], ncol=dim(spTrMat)[2],
                   dimnames=list(rownames(spTrMat), NULL))
  spTrProp = t(prop.table(spTrMat,2))
  
  spTrProp = round(t(prop.table(spTrMat,2))*1000)
  
  locNulls <- nullmodel(spTrProp, method="swap_count")
  locSims <- simulate(locNulls, nsim=simIter, thin=100)
  
  table(apply(locSims, c(1,2), function(x){length(unique(x[x>0]))}))
  
  locSims <- lapply(1:dim(locSims)[3], function(n){locSims[,,n]})
  
  locSims <- lapply(locSims, function(x){
    x = as.data.frame(prop.table(x/1000, 1))
    x$transect = 1:nrow(x)
    xLong = pivot_longer(as.data.frame(x), cols=colnames(x)[-ncol(x)],
                         names_to="taxa", values_to="prop")
    xLong <- as.data.frame(xLong)
    xLong <- cbind(xLong,
                   do.call("rbind", strsplit(xLong$taxa, ":")))
    colnames(xLong) = c("transect", "taxa", "prop", "genus", "growth.form")
    xLong$locality = loc
    return(xLong)
    
  })
  
  return(locSims)
})

# separate sims so we have data for each locality
print("Data processing and diversity calculation")
permuteSimDfs <- lapply(1:simIter, function(n){
  do.call("rbind", lapply(permuteSimList, function(x){x[[n]]}))
})

# Next: aggregate sims into gen/gr mats like main analysis
permuteSimMats <- lapply(permuteSimDfs, function(y){
  
  y$locTran <- paste0(y$locality, y$transect)
  y = y[y$prop > 0,]
  y$genus <- as.factor(y$genus)
  y$growth.form <- as.factor(y$growth.form)
  # table genus
  subGenMat <- do.call("rbind", 
                       lapply(split(y, f= ~y$transect + y$locality), 
                              function(x){
                                temp <- tapply(x$prop, x$genus, sum, na.rm=TRUE)
                                ifelse(is.na(temp), 0, temp)
                              }))
  subGenMat <- prop.table(subGenMat, 1)
  
  subGrMat <- do.call("rbind", 
                      lapply(split(y, f= ~y$transect + y$locality), 
                             function(x){
                               temp <- tapply(x$prop, x$growth.form, sum, na.rm=TRUE)
                               ifelse(is.na(temp), 0, temp)
                             }))
  subGrMat <- prop.table(subGrMat, 1)
  
  return(list(subGenMat, subGrMat))
  
})

# caclulate alpha and beta matrices
permuteSimDivs <- lapply(permuteSimMats, function(x){
  
  genMat = x[[1]][complete.cases(x[[1]]),]
  grMat = x[[2]][complete.cases(x[[2]]),]
  
  # just need to superimpose transect numbers running oldest to youngest
  tranNums <- do.call("rbind", lapply(split(huon_site, f=huon_site$locality), function(x){
    x <- x[order(x$transect.age, decreasing=TRUE),]
    x$locTran <- paste0(x$locality, 1:nrow(x))
    return(x)
  }))
  
  rownames(genMat) = tranNums$transect
  rownames(grMat) = tranNums$transect
  
  # alpha Mat
  alphaMat = alphaDiversityCalc(siteData = huon_site,
                                genusMatProp = genMat,
                                grMatProp = grMat)
  
  # beta Mat
  betaMat = betaDiversityCalc(siteData = huon_site,
                              genusMatProp = genMat,
                              grMatProp = grMat,
                              method = "bray")
  
  return(list(alphaMat, betaMat))
  
})

# run models
print("Running diversity models")
permuteModels <- lapply(1:length(permuteSimDivs), function(n){
  print(n)
  x=permuteSimDivs[[n]]
  alphaMs = mvGamsFreqRegion(dataList = list(x[[1]], x[[1]], x[[1]]),
                             primaryvarList = c("genusH0", "genusH1", "genusH2"),
                             secondaryvarList = c("grH0", "grH1", "grH2"),
                             dateDiff=FALSE)
  
  betaMs =  mvGamsFreqRegion(dataList = list(x[[2]], x[[2]]),
                             primaryvarList = c("top.gen.diss", "seq.gen.diss"),
                             secondaryvarList = c("top.gr.diss", "seq.gr.diss"),
                             dateDiff=TRUE)
  
  return(list(alphaMs, betaMs))
  
})

#                             save files ####

saveRDS(markovModels, "./outputs/markovSimModels.rds")
saveRDS(permuteModels, "./outputs/permuteSimModels.rds")

markovModels <- readRDS("./outputs/markovSimModels.rds")
permuteModels <- readRDS("./outputs/permuteSimModels.rds")

#                  aggregate predictions ####

markovAlphaRegH0 <- do.call("rbind", lapply(1:length(markovModels), function(n){
  sub = markovModels[[n]][[1]][[1]]$preds[[1]]
  sub$sim = n
  return(sub)

}))
markovAlphaRegH1 <- do.call("rbind", lapply(1:length(markovModels), function(n){
  sub = markovModels[[n]][[1]][[2]]$preds[[1]]
  sub$sim = n
  return(sub)
}))
markovAlphaRegH2 <- do.call("rbind", lapply(1:length(markovModels), function(n){
  sub = markovModels[[n]][[1]][[3]]$preds[[1]]
  sub$sim = n
  return(sub)
}))
markovBetaRegDir <- do.call("rbind", lapply(1:length(markovModels), function(n){
  sub = markovModels[[n]][[2]][[1]]$preds[[1]]
  sub$sim = n
  return(sub)
}))
markovBetaRegSeq <- do.call("rbind", lapply(1:length(markovModels), function(n){
  sub = markovModels[[n]][[2]][[2]]$preds[[1]]
  sub$sim = n
  return(sub)
}))

# aggregate predictions
permuteAlphaRegH0 <- do.call("rbind", lapply(1:length(permuteModels), function(n){
  sub = permuteModels[[n]][[1]][[1]]$preds[[1]]
  sub$sim = n
  return(sub)
  
}))
permuteAlphaRegH1 <- do.call("rbind", lapply(1:length(permuteModels), function(n){
  sub = permuteModels[[n]][[1]][[2]]$preds[[1]]
  sub$sim = n
  return(sub)
}))
permuteAlphaRegH2 <- do.call("rbind", lapply(1:length(permuteModels), function(n){
  sub = permuteModels[[n]][[1]][[3]]$preds[[1]]
  sub$sim = n
  return(sub)
}))
permuteBetaRegDir <- do.call("rbind", lapply(1:length(permuteModels), function(n){
  sub = permuteModels[[n]][[2]][[1]]$preds[[1]]
  sub$sim = n
  return(sub)
}))
permuteBetaRegSeq <- do.call("rbind", lapply(1:length(permuteModels), function(n){
  sub = permuteModels[[n]][[2]][[2]]$preds[[1]]
  sub$sim = n
  return(sub)
}))

#                                   plot ####

pdf("./simulationModels.pdf", height=10, width=4.5, useDingbats = FALSE)

par(mar=c(2,2,1,1), oma=c(2,2,1,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1,
    mfcol=c(5,2))

sapply(1:10, function(n){
  
print(n)
obsM <- list(alpha.model.list[[1]],
             alpha.model.list[[2]],
             alpha.model.list[[3]],
             beta.model.list[[1]],
             beta.model.list[[2]])[[rep(1:5, 2)[n]]]$preds[[1]]

nullM <- list(markovAlphaRegH0,
              markovAlphaRegH1,
              markovAlphaRegH2,
              markovBetaRegDir,
              markovBetaRegSeq,
              permuteAlphaRegH0,
              permuteAlphaRegH1,
              permuteAlphaRegH2,
              permuteBetaRegDir,
              permuteBetaRegSeq)[[n]]

lims <- rbind(c(4,16,3,11),
              c(0,12,0,8),
              c(0,10,0,7),
              c(0.1,0.7,0.1,0.7),
              c(0,0.8,0,0.6))[rep(1:5, 2)[n],]

nullCols <- list(c("red", rgb(1,0.8,0.8,1)),
                 c("blue", rgb(0.8,0.8,1,1)))[[rep(1:2, each=5)[n]]]


xlabs=list(expression("Tax. "*alpha*" diversity (H0)"),
           expression("Tax. "*alpha*" diversity (H1)"),
           expression("Tax. "*alpha*" diversity (H2)"),
           expression("Tax. "*beta*" diversity (Dir)"),
           expression("Tax. "*beta*" diversity (Seq)"))[[rep(1:5, 2)[n]]]

ylabs=list(expression("Fun. "*alpha*" diversity (H0)"),
           expression("Fun. "*alpha*" diversity (H1)"),
           expression("Fun. "*alpha*" diversity (H2)"),
           expression("Fun. "*beta*" diversity (Dir)"),
           expression("Fun. "*beta*" diversity (Seq)"))[[rep(1:5, 2)[n]]]


plot(NULL, xlim=lims[1:2], ylim=lims[3:4],xlab="", ylab="", xaxt="n")
axis(side=1, mgp=c(3,0.2,0))
mtext(side=1, line=1.25, text=xlabs, cex=0.8)
mtext(side=2, line=1.5, text=ylabs, las=0, cex=0.8)

sapply(split(nullM, f=nullM$sim), function(x){
  lines(x$fit.2 ~ x$fit.1, col=nullCols[2])
  arrows(x0=x$fit.1[10], x1=x$fit.1[1],
         y0=x$fit.2[10], y1=x$fit.2[1], col=nullCols[2], length=0.05)
})
simMeanFit.1 = tapply(nullM$fit.1, nullM$pred.date, mean)
simMeanFit.2 = tapply(nullM$fit.2, nullM$pred.date, mean)
lines(simMeanFit.2 ~ simMeanFit.1, lwd=2, col=nullCols[1], lty="11")
arrows(x0=simMeanFit.1[10], x1=simMeanFit.1[1],
       y0=simMeanFit.2[10], y1=simMeanFit.2[1], length=0.1, lwd=2, col=nullCols[1])
lines(obsM$fit.2 ~ obsM$fit.1, lwd=2)
arrows(x0=obsM$fit.1[2], x1=obsM$fit.1[1],
       y0=obsM$fit.2[2], y1=obsM$fit.2[1], length=0.1, lwd=2)
text(x=relative.axis.point(0.02, 'x'), y=relative.axis.point(0.95, "y"),
     adj=0, labels=paste0("(", LETTERS[n],")"), font=2)

if(n==1){mtext(side=3, line=0.1, text="Markov chains", font=2)}
if(n==6){mtext(side=3, line=0.1, text="Matrix permutations", font=2)}

})
dev.off()

#       markov simulation RA vs dRA plot ####

pdf("./plots/suppMarkovRAplot.pdf", height=4.5, width=4.5, useDingbats = FALSE)
par(mar=c(3,3,1,1), ps=10, tcl=-0.25, las=1, mgp=c(3,0.5,0))
plot(NULL, xlim=c(0,1), ylim=c(0,0.6), xlab="", ylab="", xaxt="n")
abline(a=0,b=1, lty="31")
axis(side=1, mgp=c(3,0.2,0))
mtext(side=1, line=1.5, text="Relative abundance of genus/growth form pair")
mtext(side=2, line=1.75, text=expression("|"*Delta*"Relative abundance|"), las=0)
sapply(1:length(markovSimList), function(n){
x<-markovSimList[[n]]
subDf <- x$model
print(summary(x))

loc <- locality.number[locality.number$locality == unique(huon_site$locality)[n],]
locCol <- loc.color$col[loc.color$locality == loc$locality]
locCol <- locCol[!is.na(locCol)]
locRgb <- col2rgb(locCol)/255
points(subDf$`abs(dra)` ~ subDf$ra, pch=16, cex=0.5, col=colorRampPalette(c(locCol, "white"))(3)[2])
})

sapply(1:length(markovSimList), function(n){
  
  x<-markovSimList[[n]]
  subDf <- x$model
  print(summary(x))

  loc <- locality.number[locality.number$locality == unique(huon_site$locality)[n],]
  locCol <- loc.color$col[loc.color$locality == loc$locality]
  locCol <- locCol[!is.na(locCol)]
  locRgb <- col2rgb(locCol)/255
  
  subPreds <- predict(x, newdata=data.frame(ra = seq(min(subDf$ra), max(subDf$ra), len=200)), se.fit=TRUE)
  
  polygon(x=c(seq(min(subDf$ra), max(subDf$ra), len=200), rev(seq(min(subDf$ra), max(subDf$ra), len=200))),
          y=c(subPreds$fit + 1.96 * subPreds$se.fit, rev(subPreds$fit - 1.96 * subPreds$se.fit)),
          border=NA, col=rgb(locRgb[1],locRgb[2],locRgb[3],0.25))
  lines(subPreds$fit ~ seq(min(subDf$ra), max(subDf$ra), len=200), col=locCol)
  text(x=max(subDf$ra), y=rev(subPreds$fit)[1],
       labels=loc$number, pos=4, offset=0.25, col=locCol, font=2)
})
dev.off()

# OTHER PLOTS ####

#   Dir vs seq beta diversity comparison ####
betaExamplePlot("./Plots/turnover ords.pdf")
#                               Site map ####

siteMapPlot("./Plots/huon map rotated.pdf")

#                        Regional trends ####
#                                                 alpha div ####
# funtion to plot trend arrow in diversity space as well as segments for age
main.trend <- function(preds, col){
  
  fit1 <- preds$fit.1
  fit2 <- preds$fit.2
  pred.date <- preds$pred.date
  
  if(var(fit1) < 1e-3 & var(fit2) < 1e-3){
  
    points(fit2[1] ~ fit1[1], col=col, pch=16, cex=0.85)
  } else {
  
  lines(fit2 ~ fit1, col=col)  
  Arrows(x0=fit1[3], x1=fit1[1],
         y0=fit2[3], y1=fit2[1],
         arr.type="triangle", arr.length=0.1,
         arr.width=0.1, col=col)
  
} 
}

raw.plot <- function(ylims){
  par(mar=c(0,0,0,0), ps=10, tcl=-0.25, mgp=c(3,0,0), las=1)
  plot(x=NULL, y=NULL, xlim=c(10000,6250), ylim=ylims,
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
  
  # sapply(max(preds$pred.date) - seq(0,5000,250), function(seg.n){
  #   
  #   temp.match <- which(preds$pred.date == seg.n)
  #   base.age <- max(preds$pred.date)
  #   diff.age <- base.age - seg.n
  #   
  #   if(length(temp.match)>0){
  #     
  #     target <- preds[temp.match - c(1,0),]
  #     if(nrow(target)<2){return(NULL)}
  #     arrows(x0=target$pred.date[1], x1=target$pred.date[2],
  #            y0=target$trend[2], y1=target$trend[2],
  #            angle=90, col=line.col,
  #            length=ifelse(diff.age %% 1000 == 0, 0.075,0.025),
  #            lwd=ifelse(diff.age %% 1000 == 0, 1,0.5))
  #     
  #   }
  #   
  # })
  
}

pdf(date.wrap("./Plots/alpha region trends",".pdf"), 
    height=8.5, width=5.75, useDingbats = FALSE)

height <- ((0.97-0.06)-0.1) / 3
bots <- c(0.06) + height * c(0:2) + c(0, 0.05, 0.1)
tops <- c(0.06) + height * c(1:3)

xlims=list(c(5,12), c(2,6), c(1.5,4.5))
ylims=list(c(5,12), c(2,6), c(1.5,4.5))

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

sapply(1:3, function(n){
                        
screen(n)
print(n)
x.lims<-xlims[[n]]
y.lims<-ylims[[n]]
x.lab<-xlabs[[n]]
y.lab<-ylabs[[n]]

if(n <= 3){
temp.preds <- alpha.model.list[[n]]$preds[[1]]
temp.data = alpha.model.list[[n]]$models[[1]]$model
} else {
temp.preds <- beta.model.list[[n-3]]$preds[[1]]
temp.data = beta.model.list[[n-3]]$models[[1]]$model
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

raw.plot(range(temp.data$primary) + (diff(range(temp.data$primary)) * 0.05) * c(-1,1))
axis(side=1, at=seq(6000,9000,1000), labels=NA)
axis(side=2, mgp=c(3,0.5,0))
mtext(side=2, line=1.75, text=x.lab, las=0)

points(temp.data$primary ~ temp.data$pred.date, pch=16, cex=0.6, col="grey70")

raw.trend(preds=temp.preds,
          fit.var="fit.1",line.col="black", poly.col=rgb(0.7,0.7,0.7,0.5))
text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.9, "y"),
     labels=paste0("(",LETTERS[c(2,5,8)][n],")"), font=2, adj=0)
box()
close.screen(n+(n-1)+3)

screen(n+(n-1)+4)

fit1.ylims <- c(min(temp.preds$fit.2 - 1.96 * temp.preds$se.fit.2),
                max(temp.preds$fit.2 + 1.96 * temp.preds$se.fit.2)) * c(0.9,1.1)

raw.plot(range(temp.data$secondary) + (diff(range(temp.data$secondary)) * 0.05) * c(-1,1))
axis(side=1, at=seq(6000,9000,1000), labels=format(seq(6000,9000,1000), big.mark=","), mgp=c(3,0.2,0))
mtext(side=1, line=1, text="Time (years before present)")
axis(side=2, mgp=c(3,0.5,0))
mtext(side=2, line=1.75, text=y.lab, las=0)
points(temp.data$secondary ~ temp.data$pred.date, pch=16, cex=0.6, col="grey70")
raw.trend(preds=temp.preds,
          fit.var="fit.2",line.col="black", poly.col=rgb(0.7,0.7,0.7,0.5))
text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.9, "y"),
     labels=paste0("(",LETTERS[c(3,6,9)][n],")"), font=2, adj=0)
box()
close.screen(n+(n-1)+4)

})
close.screen(all.screens=TRUE)
dev.off()

#                                                  beta div ####

pdf(date.wrap("./Plots/beta region trends",".pdf"), 
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
xlims=list(c(0,0.55), c(0,0.55))
ylims=list(c(0,0.55), c(0,0.55))

sapply(1:2, function(n){
  
  screen(n)
  print(n)
  x.lims<-xlims[[n]]
  y.lims<-ylims[[n]]
  x.lab<-xlabs[[n]]
  y.lab<-ylabs[[n]]
  
  temp.preds <- beta.model.list[[n]]$preds[[1]]
  temp.data <- beta.model.list[[n]]$models[[1]]$model

  par(mar=c(0,0,0,0), ps=10, las=1, tcl=-0.25, mgp=c(3,0.5,0))
  
  plot(x=NULL, y=NULL, xlim=x.lims, ylim=y.lims, asp=1,
       xlab="",ylab="", yaxs="i", xaxs="i", axes=FALSE)
  axis(side=1, mgp=c(3,0.2,0))
  axis(side=2, las=1)
  mtext(side=2, line=1.35, text=y.lab, las=0)
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
  
  raw.plot(range(temp.data$primary) + (diff(range(temp.data$primary)) * 0.05) * c(-1,1))
  axis(side=1, at=seq(6000,9000,1000), labels=NA)
  axis(side=2, mgp=c(3,0.5,0))
  mtext(side=2, line=2, text=x.lab, las=0)
  points(temp.data$primary ~ temp.data$pred.date, pch=16, cex=0.6, col="grey70")
  raw.trend(preds=temp.preds,
            fit.var="fit.1",line.col="black", poly.col=rgb(0.7,0.7,0.7,0.5))
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.9, "y"),
       labels=paste0("(",LETTERS[c(2,5,8)][n],")"), font=2, adj=0)
  box()
  close.screen(n+(n-1)+3)
  
  screen(n+(n-1)+4)
  
  fit1.ylims <- c(min(temp.preds$fit.2 - 1.96 * temp.preds$se.fit.2),
                  max(temp.preds$fit.2 + 1.96 * temp.preds$se.fit.2)) * c(0.9,1.1)
  
  raw.plot(range(temp.data$secondary) + (diff(range(temp.data$secondary)) * 0.05) * c(-1,1))
  axis(side=1, at=seq(6000,9000,1000), labels=format(seq(6000,9000,1000), big.mark=","), mgp=c(3,0.2,0))
  mtext(side=1, line=1, text="Time (years before present)")
  axis(side=2, mgp=c(3,0.5,0))
  mtext(side=2, line=2, text=y.lab, las=0)
  points(temp.data$secondary ~ temp.data$pred.date, pch=16, cex=0.6, col="grey70")
  raw.trend(preds=temp.preds,
            fit.var="fit.2",line.col="black", poly.col=rgb(0.7,0.7,0.7,0.5))
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.9, "y"),
       labels=paste0("(",LETTERS[c(3,6,9)][n],")"), font=2, adj=0)
  box()
  close.screen(n+(n-1)+4)
  
})
close.screen(all.screens=TRUE)
dev.off()

#                            Site trends ####

loc.color <- data.frame(locality=unique(huon_coral$locality))
loc.color$col <- c("blue","red","darkgreen")[huon_tr$site[match(loc.color$locality,
                                                                huon_tr$locality)]]
loc.color = loc.color[!is.na(loc.color$locality),]
rich.df <- richDf

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
xlims=list(c(3,16),c(1.5,7),c(1,5.5), c(0,0.75), c(0,0.6))
ylims=list(c(3,16),c(1.5,7),c(1,5.5), c(0,0.7), c(0,0.6))

pdf(date.wrap("./Plots/site trends",".pdf"), height=8.5, width=5.75, bg="white")

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
    site.range <- range(rich.df$pred.date[rich.df$locality==x])
    
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
                      

#                        Trends with CIs ####

loc.color <- data.frame(locality=unique(huon_coral$locality))
loc.color$col <- c("blue","red","darkgreen")[huon_tr$site[match(loc.color$locality,
                                                                huon_tr$locality)]]
loc.color = loc.color[!is.na(loc.color$locality),]


# Taxonomic alpha plots
sapply(1:length(alpha.model.list), function(n1){
  
  locality.preds <- alpha.model.list[[n1]]$preds[[3]]
  raw.data <- richDf
  
  fit1.ylims<-rbind(c(2,18), c(0.5,9), c(0,6.5))[n1,]
  fit1.ylab <- c(expression("Taxonomic "*alpha*" diversity (Hill q = 0)"),
                 expression("Taxonomic "*alpha*" diversity (Hill q = 1)"),
                 expression("Taxonomic "*alpha*" diversity (Hill q = 2)"))[n1]
  fit1.var <- c("genusH0", "genusH1", "genusH2")[n1]
  
  x.marks <- seq(6000,9000,1000)
  
  # PRIMARY AXIS PLOTS
  
  # set up split.screen. Locality is 3x3 panel, with Site and region as extra side plots
  pdf(date.wrap(paste0("./Plots/taxonomic div ",fit1.var),".pdf"), 
      height=5, width=5, useDingbats = FALSE)
  
  panel3.x <- seq(0.1,0.99,len=4)
  panel3.y <- seq(0.08,0.97,len=4)
  panels <- expand.grid(1:3, 1:3)
  panel.mat <- cbind(panel3.x[-4][panels[,1]],
                     panel3.x[-1][panels[,1]],
                     panel3.y[-4][panels[,2]],
                     panel3.y[-1][panels[,2]])
  
  split.screen(panel.mat)
  
 # locality subplots
  sapply(1:9, function(n2){
    screen(n2)
    
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
         labels=paste0("(",LETTERS[(c(9,10,11,6,7,8,3,4,5)-2)[n2]],")"), font=2, adj=0)
    close.screen(n2)
  })
  
  close.screen(all.screens=TRUE)
  dev.off()
  
})

# Functional alpha plots
sapply(1:length(alpha.model.list), function(n1){
  
  locality.preds <- alpha.model.list[[n1]]$preds[[3]]
  raw.data <- richDf
  
  fit1.ylims<-rbind(c(2.75,12), c(1,7), c(0.75,5))[n1,]
  fit1.ylab <- c(expression("Functional "*alpha*" diversity (Hill q = 0)"),
                 expression("Functional "*alpha*" diversity (Hill q = 1)"),
                 expression("Functional "*alpha*" diversity (Hill q = 2)"))[n1]
  fit1.var <- c("grH0", "grH1", "grH2")[n1]
  
  x.marks <- seq(6000,9000,1000)
  
  # PRIMARY AXIS PLOTS
  
  # set up split.screen. Locality is 3x3 panel, with Site and region as extra side plots
  pdf(date.wrap(paste0("./Plots/functional div ",fit1.var),".pdf"), 
      height=5, width=5, useDingbats = FALSE)
  
  panel3.x <- seq(0.1,0.99,len=4)
  panel3.y <- seq(0.08,0.97,len=4)
  panels <- expand.grid(1:3, 1:3)
  panel.mat <- cbind(panel3.x[-4][panels[,1]],
                     panel3.x[-1][panels[,1]],
                     panel3.y[-4][panels[,2]],
                     panel3.y[-1][panels[,2]])
  
  split.screen(panel.mat)
  
  # locality subplots
  sapply(1:9, function(n2){
    screen(n2)
    
    raw.plot(fit1.ylims)
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
         labels=paste0("(",LETTERS[(c(9,10,11,6,7,8,3,4,5)-2)[n2]],")"), font=2, adj=0)
    close.screen(n2)
  })
  
  close.screen(all.screens=TRUE)
  dev.off()
  
})

# Taxonomic beta plots
sapply(1:length(beta.model.list), function(n1){
  
  locality.preds <- beta.model.list[[n1]]$preds[[3]]
  raw.data <- turnoverDf
  
  fit1.ylims<-rbind(c(0,0.775), c(-0.1,0.75))[n1,]
  fit1.ylab <- c(expression("Taxonomic "*beta*" diversity (directional)"),
                 expression("Taxonomic "*beta*" diversity (sequential)"))[n1]
  fit1.var <- c("top.gen.diss", "seq.gen.diss")[n1]
  
  x.marks <- seq(6000,9000,1000)
  
  # PRIMARY AXIS PLOTS
  
  # set up split.screen. Locality is 3x3 panel, with Site and region as extra side plots
  pdf(date.wrap(paste0("./Plots/taxonomic beta ",fit1.var),".pdf"), 
      height=5, width=5, useDingbats = FALSE)
  
  panel3.x <- seq(0.1,0.99,len=4)
  panel3.y <- seq(0.08,0.97,len=4)
  panels <- expand.grid(1:3, 1:3)
  panel.mat <- cbind(panel3.x[-4][panels[,1]],
                     panel3.x[-1][panels[,1]],
                     panel3.y[-4][panels[,2]],
                     panel3.y[-1][panels[,2]])
  
  split.screen(panel.mat)
  
  # locality subplots
  sapply(1:9, function(n2){
    screen(n2)
    
    raw.plot(fit1.ylims)
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
         labels=paste0("(",LETTERS[(c(9,10,11,6,7,8,3,4,5)-2)[n2]],")"), font=2, adj=0)
    close.screen(n2)
  })
  
  close.screen(all.screens=TRUE)
  dev.off()
  
})

# Functional beta plots
sapply(1:length(beta.model.list), function(n1){
  
  locality.preds <- beta.model.list[[n1]]$preds[[3]]
  raw.data <- turnoverDf
  
  fit1.ylims<-rbind(c(0,0.775), c(-0.1,0.45))[n1,]
  fit1.ylab <- c(expression("Functional "*beta*" diversity (directional)"),
                 expression("Functional "*beta*" diversity (sequential)"))[n1]
  fit1.var <- c("top.gr.diss", "seq.gr.diss")[n1]
  
  x.marks <- seq(6000,9000,1000)
  
  # PRIMARY AXIS PLOTS
  
  # set up split.screen. Locality is 3x3 panel, with Site and region as extra side plots
  pdf(date.wrap(paste0("./Plots/functional beta ",fit1.var),".pdf"), 
      height=5, width=5, useDingbats = FALSE)
  
  panel3.x <- seq(0.1,0.99,len=4)
  panel3.y <- seq(0.08,0.97,len=4)
  panels <- expand.grid(1:3, 1:3)
  panel.mat <- cbind(panel3.x[-4][panels[,1]],
                     panel3.x[-1][panels[,1]],
                     panel3.y[-4][panels[,2]],
                     panel3.y[-1][panels[,2]])
  
  split.screen(panel.mat)
  
  # locality subplots
  sapply(1:9, function(n2){
    screen(n2)
    
    raw.plot(fit1.ylims)
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
         labels=paste0("(",LETTERS[(c(9,10,11,6,7,8,3,4,5)-2)[n2]],")"), font=2, adj=0)
    close.screen(n2)
  })
  
  close.screen(all.screens=TRUE)
  dev.off()
  
})

#            Matrix of genus/growth form ####

gr.gen.mat <- table(huon_coral_all$genus, huon_coral_all$growth.form)
gr.gen.mat <- gr.gen.mat[rowSums(gr.gen.mat) > 0, ]
gen.count <- rowSums(gr.gen.mat)

#Other <- colSums(gr.gen.mat[gen.count < 10, ])
#commons <- gr.gen.mat[gen.count >= 10, ]
#final.mat <- rbind(commons, Other)
final.mat <- gr.gen.mat
#final.mat <- final.mat[, c(2,1,3:ncol(final.mat))]

gen.count <- rowSums(final.mat) / sum(final.mat)
gr.count <- colSums(final.mat) / sum(final.mat)

gr.label.mat <- data.frame(col = colnames(final.mat),
                           label = c("Branch (closed)",
                                     "Branch (open)",
                                     "Columnar",
                                     "Corymbose",
                                     "Digitate",
                                     "Encrusting",
                                     "Hispidose",
                                     "Laminar",
                                     "Massive",
                                     "Solitary",
                                     "Submassive",
                                     "Tabular or Plate", 
                                     "Unknown"),
                           pos = 1:ncol(final.mat), stringsAsFactors = FALSE)

final.mat <- final.mat[nrow(final.mat):1, order(gr.label.mat$pos)]
colnames(final.mat) <- gr.label.mat$label[order(gr.label.mat$pos)]
gr.count <- gr.count[order(gr.label.mat$pos)]
#final.mat <- final.mat / rowSums(final.mat)

pdf(date.wrap("./Plots/genus growth form matrix", ".pdf"), height=7.5, width=10)
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

genus.perc <- rev(round(gen.count, 4))*100
gr.perc <- rev(round(gr.count, 4))*100
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
     labels = sprintf("%.2f", round(rev(gen.count), 4)*100), cex=0.75,
     col = ifelse(rev(gen.count) > 0.5, "white", "black"), font=2)

# gr labels
text(x = seq(gr.x, end.x, len=ncol(final.mat)),
     y = gen.y -  0.05,
     labels = colnames(final.mat), srt=20, adj=c(1,1))

# growth prop
text(x = seq(gr.x, end.x, len=ncol(final.mat)),
     y =  gen.y - 3 * y.width,
     labels = sprintf("%.2f", round(gr.count, 4)*100), cex=0.75,
     col = ifelse(gr.count > 0.25, "white", "black"), font=2)

# Totals text

text(x=gr.x - 2.5 * x.width, 
     y = gen.y - 3 * y.width, 
     labels="Totals\n(%)", font=2)

dev.off()

fullTable <- as.data.frame(table(huon_coral_all$growth.comb, huon_coral_all$binom))
fullTable <- merge(fullTable, huon_coral_all[!duplicated(huon_coral_all$growth.comb),c("growth.form", "growth.comb")],
                   by.x="Var1", by.y="growth.comb", all.x=TRUE, all.y=FALSE, sort=FALSE)
fullTable <- fullTable[fullTable$Freq > 0,]
fullTable <- fullTable[order(fullTable$Freq, decreasing=TRUE),]
colnames(fullTable) <- c("growth assess", "species", "count", "growth.form")
write.csv(fullTable, "./outputs/fullGrowthSummaryFormed.csv")

fullTable <- as.data.frame(table(huon_coral_all$growth.form, huon_coral_all$binom))
fullTable <- fullTable[fullTable$Freq > 0,]
fullTable <- fullTable[order(fullTable$Freq, decreasing=TRUE),]
write.csv(fullTable, "./outputs/speciesGrowth.csv")

# SUPP ANALYSES ####
# diversity vs sedimentation rate ####

# we need to pull in the age models and calculate slope using finite differencing
# then compare this to observed diversity

sedDiv <- do.call("rbind", lapply(levels(huon_coral$locality), function(loc){
  
  divLoc <- richDf[richDf$locality == loc,]
  turnLoc <- turnoverDf[turnoverDf$locality == loc,]

  combDiv <- merge(divLoc, turnLoc[,!colnames(turnLoc) %in% c("site", "locality", "date.diff", "pred.date", "top.date.diff")],
                   by.x="transect", by.y="transect")
  
  # we need transect heights as well
  huon_site <- read.csv("./raw.datafiles/huon_transect_data.csv", stringsAsFactors = TRUE)
  
  combDiv$tranHeight <- huon_site$height.from.dist[match(combDiv$transect, huon_site$transect)]
  
  # bring in the bacon preds
  bacon.dirs <- paste0("/Users/uqtstapl/Dropbox/Tim/Post-doc/Research projects/huon_cleaning/Bacon_runs/", loc)
  bacon.files <- list.files(bacon.dirs, pattern="_ages")
  bacon.ages <- read.table(paste0(bacon.dirs, "/", bacon.files), header=TRUE)
  
  heightDepth = list.files(bacon.dirs, pattern="depthsDist.txt")
  bacon.dist.height = unlist(read.table(paste0(bacon.dirs, "/", heightDepth)))
  bacon.ages$height <- (max(bacon.dist.height) - bacon.ages$depth) / 100
  
  tranRows <- match(combDiv$tranHeight,bacon.ages$height)
  tranSlope = cbind((bacon.ages$mean[tranRows+1] - bacon.ages$mean[tranRows]) / abs(bacon.ages$height[tranRows+1] - bacon.ages$height[tranRows]),
                    (bacon.ages$mean[tranRows] - bacon.ages$mean[ifelse(tranRows-1 == 0, 1, tranRows-1)]) / abs(bacon.ages$height[tranRows] - bacon.ages$height[ifelse(tranRows-1 == 0, 1, tranRows-1)]))
  tranSlope[tranSlope == 0 | tranSlope == -Inf] = NA
  combDiv$tranSlope <- rowMeans(tranSlope, na.rm=TRUE)
  return(combDiv)
}))
cor(sedDiv[,c("date.diff", "tranSlope")], use="complete.obs")


# for each diversity type
divVars <- c("genusH0", "genusH1", "genusH2", "top.gen.diss", "seq.gen.diss", "grH0", "grH1", "grH2", "top.gr.diss", "seq.gr.diss")
pdf("./plots/sedimentSlopes.pdf", height=8, width=4.5, useDingbats = FALSE)
par(mfcol=c(5,2), mar=c(0,2.5,0,1), oma=c(3.5,1,0.5,0), ps=10, tcl=-0.25, mgp=c(3,0.5,0), las=1)
sedModels <- lapply(1:length(divVars),
       function(n){
div <- divVars[n]
tempData <- sedDiv
tempData$div = tempData[,div]
tempData <- tempData[!is.na(tempData$div),]

tempData$invSlope <- (1/tempData$tranSlope) * 100
if(grepl("diss", div)){
  divM <- lm(div ~ invSlope*locality + date.diff, data=tempData)
  divReg <- lmer(div ~ invSlope + date.diff + (1|locality), data=tempData)
} else {
  divM <- lm(div ~ invSlope*locality, data=tempData)
  divReg <- lmer(div ~ invSlope + (1|locality), data=tempData)
}  

print(div)
print(summary(divReg)$coefficients)

slopeLims <- tapply(tempData$invSlope, tempData$locality, range)
predDf <- do.call("rbind", 
                  lapply(slopeLims, function(x){data.frame(invSlope=seq(x[1], x[2], len=200),
                                                           date.diff = mean(tempData$date.diff, na.rm=TRUE))}))
predDf$locality <- as.factor(substr(rownames(predDf), 1, regexpr("\\.", rownames(predDf))-1))
predDf <- cbind(predDf, predict(divM, newdata=predDf, se.fit=TRUE))

regPred <- data.frame(invSlope = seq(min(tempData$invSlope), max(tempData$invSlope), len=200),
                      date.diff = mean(tempData$date.diff, na.rm=TRUE))
regPred <- cbind(regPred, mer.ci(divReg, regPred, 99))

loc.color <- loc.color[complete.cases(loc.color),]
plot(tempData$div ~ tempData$invSlope, pch=16, col=loc.color$col[match(tempData$locality, loc.color$locality)],
     xlab="", ylab="", xaxt="n",
     xlim=c(0.1, 1.2), ylim=range(c(tempData$div, predDf$fit), na.rm=TRUE))
rect(xleft=par("usr")[1], xright=par('usr')[2], ybottom=par('usr')[3], ytop=par("usr")[4],
     border=NA, col=rgb(1,1,1,0.5))
mtext(side=2, line=1.75, text=c("Genus H0", "Genus H1", "Genus H2", 
                               expression("Dir Genus "*beta), expression("Seq Genus "*beta),
                               "Growth form H0", "Growth form H1", "Growth form H2",
                               expression("Dir Growth form "*beta), expression("Seq Growth form "*beta))[n], 
      las=0, cex=0.8)

if(n %in% c(5,10)){
axis(side=1, mgp=c(3,0.2,0))
} else {
axis(side=1, labels=NA)
}

if(n == 10){mtext(side=1, line=1.75, at=relative.axis.point(-0.2, "x"), 
                  text=expression("Sediment accumulation rate (cm year"^-1*")"), cex=0.8)
}
  
sapply(split(predDf, f=predDf$locality), function(x){
  tempCol <- loc.color$col[match(x$locality[1], loc.color$locality)]
  tempRgb <- col2rgb(tempCol)/255
  
  # polygon(x=c(x$invSlope, rev(x$invSlope)),
  #         y=c(x$fit + 1.96 * x$se.fit, rev(x$fit - 1.96 * x$se.fit)),
  #         col=rgb(tempRgb[1], tempRgb[2], tempRgb[3], 0.1), border=NA)
  lines(x$fit ~ x$invSlope, col=tempCol, lwd=1)
  text(rev(x$fit)[1] ~ rev(x$invSlope)[1], 
       labels=locality.number$number[match(x$locality[1], locality.number$locality)],
       col=tempCol,
       pos=4, font=2)
})

polygon(x=c(regPred$invSlope, rev(regPred$invSlope)),
        y=c(regPred$lower, rev(regPred$upper)),
        col=rgb(0,0,0, 0.25), border=NA)
lines(regPred$fit ~ regPred$invSlope, lwd=2)

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.95, "y"),
     labels=paste0("(", LETTERS[n], ")"), adj=0, font=2)

return(list(locM = divM,
            regM = divReg))

})
dev.off()

sedReg <- lapply(sedModels, function(x){x[[2]]})
compare_performance(sedReg)
sedLoc <- lapply(sedModels, function(x){x[[1]]})
compare_performance(sedLoc)

# effect of time differences ####

lapply(beta.model.list, function(x){
  lapply(x$models, function(y){
    ptable <- summary(y)$p.table
    ptable[grepl("date\\.diff", rownames(ptable)),]
    })
})
