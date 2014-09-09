graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="SystemsJags" # for constructing output filenames
source("openGraphSaveGraph.R")
source("plotPost.R")
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.
modelstring = "
model {
    for( i in 1 : Ndata ) {
        y[i] ~ dnorm( mu[ subj[i] ] , tau[ subj[i] ] )
    }
    for ( j in 1 : Nsubj ) {
        mu[j] ~ dnorm( muG , tauG )
        tau[j] ~ dgamma( sG , rG )
    }
    muG ~ dnorm( 2.3 , 0.1 )
    tauG ~ dgamma( 1 , .5 )
    sG <- pow(m,2) / pow(d,2)
    rG <- m / pow(d,2)
    m ~ dgamma( 1 , .25 )
    d ~ dgamma( 1 , .5 )
}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.

# Load the aircraft data:
load( "Systems.Rdata" ) # loads dataMat
nSubj = length( unique( dataMat[,"Aircraft"] ) )
# Transform the data:
DaysTransf = dataMat[,"Days"]^(1/5)
dataMat = cbind( dataMat , DaysTransf )
colnames( dataMat ) = c( colnames( dataMat )[1:3] , "DaysTransf" )

# Specify data, as a list.
dataList = list(
    y = dataMat[,"DaysTransf"] ,
    subj = dataMat[,"Aircraft"] ,
    Ndata = NROW(dataMat) ,
    Nsubj = nSubj
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

# initsList = list( ** )

#------------------------------------------------------------------------------
# RUN THE CHAINS

parameters = c( "muG" , "tauG" , "mu" , "tau" , "m" , "d") 
adaptSteps = 1000              # Number of steps to "tune" the samplers.
burnInSteps = 5000            # Number of steps to "burn-in" the samplers.
nChains = 10                   # Number of chains to run.
numSavedSteps=50000           # Total number of steps in chains to save.
thinSteps=1                   # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
jagsModel = jags.model( "model.txt" , data=dataList , #inits=initsList , 
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
# EXAMINE THE RESULTS

checkConvergence = FALSE
if ( checkConvergence ) {
  openGraph(width=7,height=7)
  autocorr.plot( codaSamples[[1]] , ask=FALSE )
  show( gelman.diag( codaSamples ) )
  effectiveChainLength = effectiveSize( codaSamples ) 
  show( effectiveChainLength )
}

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )

# Extract chains from BUGS into R:
muGsample = mcmcChain[, "muG" ]
tauGsample = mcmcChain[, "tauG" ]
mSample = mcmcChain[, "m" ]
dSample = mcmcChain[, "d" ]
muSample = NULL
tauSample = NULL
for ( sIdx in 1:nSubj ) {
    muSample = rbind( muSample , mcmcChain[, paste("mu[",sIdx,"]",sep="") ] )
    tauSample = rbind( tauSample , mcmcChain[, paste("tau[",sIdx,"]",sep="") ] )
}

# Plot the aircraft mu:
openGraph(width=15,height=5)
layout( matrix( 1:nSubj , nrow=2 , byrow=T ) )
xLim=quantile(muSample,probs=c(0.001,0.999))
for ( sIdx in 1:nSubj ) {
    histInfo = plotPost( muSample[sIdx,] , xlab=bquote(mu[.(sIdx)]) , xlim=xLim )
}
saveGraph(file=paste(fileNameRoot,"PostMu",sep=""),type="eps")

# Plot the aircraft tau:
openGraph(width=15,height=5)
layout( matrix( 1:nSubj , nrow=2 , byrow=T ) )
xLim=quantile(tauSample,probs=c(0.001,0.999))
for ( sIdx in 1:nSubj ) {
    histInfo = plotPost( tauSample[sIdx,] , xlab=bquote(tau[.(sIdx)]) ,
                         HDItextPlace=.5 , showMode=TRUE , xlim=xLim )
}
saveGraph(file=paste(fileNameRoot,"PostTau",sep=""),type="eps")

# Plot the hyperdistributions:
openGraph(width=15,height=3)
layout( matrix(1:4,ncol=4) )
histInfo = plotPost( muGsample , xlab=expression(mu[G]) )
histInfo = plotPost( tauGsample , xlab=expression(tau[G]) , showMode=TRUE )
histInfo = plotPost( mSample , xlab=expression(m) , showMode=TRUE )
histInfo = plotPost( dSample , xlab=expression(d) , HDItextPlace=.1 , showMode=TRUE )
saveGraph(file=paste(fileNameRoot,"PostHyper",sep=""),type="eps")

#------------------------------------------------------------------------------