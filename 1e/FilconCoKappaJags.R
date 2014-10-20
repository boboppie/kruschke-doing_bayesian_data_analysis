graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="FilconCoKappaJags" # for constructing output filenames
source("openGraphSaveGraph.R")
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
model {
   for ( subjIdx in 1:nSubj ) {
      # Likelihood:
      z[subjIdx] ~ dbin( theta[subjIdx] , N[subjIdx] )
      # Prior on theta: Notice nested indexing.
      theta[subjIdx] ~ dbeta( a[cond[subjIdx]] , b[cond[subjIdx]] )T(0.001,0.999)
   }
   for ( condIdx in 1:nCond ) {
      a[condIdx] <- mu[condIdx] * kappa[condIdx]
      b[condIdx] <- (1-mu[condIdx]) * kappa[condIdx]
      # Hyperprior on mu and kappa:
      mu[condIdx] ~ dbeta( Amu , Bmu )
      kappa[condIdx] ~ dgamma( Skappa , Rkappa )
   }
   # Constants for hyperprior:
   Amu <- 1
   Bmu <- 1
   Skappa <- pow(meanGamma,2)/pow(sdGamma,2)
   Rkappa <- meanGamma/pow(sdGamma,2)
   meanGamma ~ dunif( 0.01 , 30 )
   sdGamma ~ dunif( 0.01 , 30 )
}
" # close quote for modelstring
# Write model to a file:
writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.

# For each subject, specify the condition s/he was in,
# the number of trials s/he experienced, and the number correct.
# (These lines intentionally exceed the margins so that they don't take up
# excessive space on the printed page.)
cond = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)
N = c(64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64)
z = c(45,63,58,64,58,63,51,60,59,47,63,61,60,51,59,45,61,59,60,58,63,56,63,64,64,60,64,62,49,64,64,58,64,52,64,64,64,62,64,61,59,59,55,62,51,58,55,54,59,57,58,60,54,42,59,57,59,53,53,42,59,57,29,36,51,64,60,54,54,38,61,60,61,60,62,55,38,43,58,60,44,44,32,56,43,36,38,48,32,40,40,34,45,42,41,32,48,36,29,37,53,55,50,47,46,44,50,56,58,42,58,54,57,54,51,49,52,51,49,51,46,46,42,49,46,56,42,53,55,51,55,49,53,55,40,46,56,47,54,54,42,34,35,41,48,46,39,55,30,49,27,51,41,36,45,41,53,32,43,33)
nSubj = length(cond)
nCond = length(unique(cond))

dataList = list(
 nCond = nCond ,
 nSubj = nSubj ,
 cond = cond ,
 N = N ,
 z = z
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

sqzData = .01+.98*dataList$z/dataList$N
mu = aggregate( sqzData , list(dataList$cond) , mean )[,"x"]
sd = aggregate( sqzData , list(dataList$cond) , sd )[,"x"]
kappa = mu*(1-mu)/sd^2 - 1
initsList = list( theta = sqzData , mu = mu , kappa = kappa , 
                  meanGamma = mean( kappa ), sdGamma = sd( kappa ) )

#------------------------------------------------------------------------------
# RUN THE CHAINS.

parameters = c("mu","kappa","theta","meanGamma","sdGamma") 
adaptSteps = 1000             # Number of steps to "tune" the samplers.
burnInSteps = 2000            # Number of steps to "burn-in" the samplers.
nChains = 3                   # Number of chains to run.
numSavedSteps=50000           # Total number of steps in chains to save.
thinSteps=1                   # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
jagsModel = jags.model( "model.txt" , data=dataList , inits=initsList , 
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( jagsModel , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                            n.iter=nPerChain , thin=thinSteps )
# resulting codaSamples object has these indices: 
#   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS.

checkConvergence = FALSE
if ( checkConvergence ) {
  openGraph(width=7,height=7)
  autocorr.plot( codaSamples[[1]][,c("kappa[1]","kappa[2]","kappa[3]","kappa[4]",
                                     "mu[1]","mu[2]","mu[3]","mu[4]")] ,
                 ask=FALSE ) 
  show( gelman.diag( codaSamples ) )
  effectiveChainLength = effectiveSize( codaSamples ) 
  show( effectiveChainLength )
}

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )

# Extract parameter values and save them.
mu = NULL
kappa = NULL
for ( condIdx in 1:nCond ) {
   mu = rbind( mu , mcmcChain[, paste("mu[",condIdx,"]",sep="") ] )
   kappa = rbind( kappa , mcmcChain[, paste("kappa[",condIdx,"]",sep="") ] )
}
save( mu , kappa , mcmcChain , file=paste(fileNameRoot,"MuKappa.Rdata",sep="") )
chainLength = NCOL(mu)

# Histograms of mu differences:
openGraph(width=10,height=2.75)
layout( matrix(1:3,nrow=1) )
source("plotPost.R")
histInfo = plotPost( mu[1,]-mu[2,] , xlab=expression(mu[1]-mu[2]) , main="" ,
          compVal=0 )
histInfo = plotPost( mu[3,]-mu[4,] , xlab=expression(mu[3]-mu[4]) , main="" ,
          compVal=0 )
histInfo = plotPost( (mu[1,]+mu[2,])/2 - (mu[3,]+mu[4,])/2 ,
          xlab=expression( (mu[1]+mu[2])/2 - (mu[3]+mu[4])/2 ) ,
          main="" ,  compVal=0 )
#saveGraph(file=paste(fileNameRoot,"MuDiffs",sep=""),type="eps")

# Scatterplot of mu,kappa in each condition:
windows()
muLim = c(.60,1) ; kappaLim = c( 2.0 , 40 ) ; mainLab="Posterior"
thindex = round( seq( 1 , chainLength , len=300 ) )
plot( mu[1,thindex] , kappa[1,thindex] , main=mainLab ,
      xlab=expression(mu[c]) , ylab=expression(kappa[c]) , cex.lab=1.75 ,
      xlim=muLim , ylim=kappaLim , log="y" , col="red" , pch="1" )
points( mu[2,thindex] , kappa[2,thindex] , col="blue" , pch="2" )
points( mu[3,thindex] , kappa[3,thindex] , col="green" , pch="3" )
points( mu[4,thindex] , kappa[4,thindex] , col="black" , pch="4" )
#saveGraph(file=paste(fileNameRoot,"Scatter",sep=""),type="eps")
