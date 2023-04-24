markovSimm <- function(N, commMat, tDates, simN){

require(betareg)
  
# get rank abundance distribution of observed comm
commMat <- commMat[,colSums(commMat) > 0]
commRA <- colSums(commMat) / nrow(commMat)
commRAVar <- apply(commMat, 2, sd)
tDateDiff <- abs(diff(tDates))

# set empty matrix
initGamma <- length(commRA)
storeMat <- matrix(0, nrow=N, ncol=initGamma)

# find gamma/beta distribution that represents community
commBeta <- betareg(commRA ~ 1)

# find gamma/beta distribution that reflects variation as a function
# of relative abundance
commBetaVar <- betareg(beta.tr(commRAVar) ~ commRA)

# convert model parameters to beta shape params
betaCoef <- commBeta$coefficients
betaCoef$mean <- plogis(betaCoef$mean)

# shape parameters used in qbeta.
A = betaCoef$mean * betaCoef$precision
B = betaCoef$precision - A

lapply(1:simN, function(n){

# draw equilibrium abundances
tempMat <- storeMat

tempMat[1,] = rbeta(n=length(commRA), A, B) 
tempMat[1,] = tempMat[1,] / sum(tempMat[1,])

# draw variances
betaVars <- predict(commBetaVar, 
                    newdata=data.frame(commRA=tempMat[1,]), se.fit=TRUE)

# run simulation seeding potential occurrence   
for(i in 2:nrow(tempMat)){
  
tempMat[i,] = tempMat[(i-1),] + rnorm(ncol(tempMat), 0, betaVars)
#tempMat[i,] = tempMat[(i-1),] + rnorm(ncol(tempMat), 0, 0.1 * tDateDiff[i-1] / 100)
  tempMat[i,] = ifelse(tempMat[i,] < 0, 0, tempMat[i,])

# sampling threshold so prop of 0s better matches observed data
tempMat[tempMat <= 0.05] = 0

# re proportionalize data
tempMat[i,] = tempMat[i,] / sum(tempMat[i,])

# re-calculate variance based on new relative abundance
# draw variances
betaVars <- predict(commBetaVar, 
                    newdata=data.frame(commRA=c(tempMat[i,])), se.fit=TRUE)
  
}

return(tempMat)

})

  
}