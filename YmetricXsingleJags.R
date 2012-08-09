graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="YmetricXsingleJags" # for constructing output filenames
if ( .Platform$OS.type != "windows" ) { 
  windows <- function( ... ) X11( ... ) 
}
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.
modelstring = "
model {
    # Likelihood:
    for( i in 1 : N ) {
        y[i] ~ dnorm( mu , tau ) # tau is precision, not SD
    }
    # Prior:
    tau ~ dgamma( 0.01 , 0.01 )
    mu ~ dnorm( 0 , 1.0E-10 )
}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.

# Generate random data from known parameter values:
set.seed(47405)
trueM = 100
trueSD = 15
y = round( rnorm( n=500 , mean=trueM , sd=trueSD ) ) # R dnorm uses mean and SD

dataList = list(
    y = y ,
    N = length( y )
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

initsList = list( mu = mean( dataList$y ) , tau = 1 / sd( dataList$y )^2 )

#------------------------------------------------------------------------------
# RUN THE CHAINS

parameters = c("mu" , "tau")  # The parameter(s) to be monitored.
adaptSteps = 500              # Number of steps to "tune" the samplers.
burnInSteps = 1000            # Number of steps to "burn-in" the samplers.
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
# EXAMINE THE RESULTS

checkConvergence = T
if ( checkConvergence ) {
  show( summary( codaSamples ) )
  windows()
  plot( codaSamples , ask=F )  
  windows()
  autocorr.plot( codaSamples , ask=F )
}

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )

muSample = mcmcChain[, "mu" ]
tauSample = mcmcChain[, "tau" ]
sigmaSample <- 1 / sqrt( tauSample ) # Convert precision to SD

source("plotPost.R")
windows()
plotPost( muSample , xlab="mu" , breaks=30 , main="Posterior" )
savePlot(file=paste(fileNameRoot,"PostMu.eps",sep=""),type="eps")

nPts = length(muSample) ; nPtsForDisplay = min( nPts , 2000 )
thinIdx = seq( 1 , nPts , nPts / nPtsForDisplay )
windows()
plot( muSample[thinIdx] , sigmaSample[thinIdx] , col="gray" ,
      xlab="mu" , ylab="sigma" , cex.lab=1.5 , main="Posterior" , log="y" )
points( mean(muSample) , mean(sigmaSample) , pch="+" , cex=2 )
text( mean(muSample) , mean(sigmaSample) ,
      bquote( .(round(mean(muSample),1)) *"  "* .(round(mean(sigmaSample),1)) ),
      adj=c(.5,-0.5) )
savePlot(file=paste(fileNameRoot,"PostMuSigma.eps",sep=""),type="eps")

#------------------------------------------------------------------------------