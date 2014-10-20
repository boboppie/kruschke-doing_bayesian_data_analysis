graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="YmetricXsingleJags" # for constructing output filenames
source("openGraphSaveGraph.R")
source("plotPost.R")
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
adaptSteps = 1000              # Number of steps to "tune" the samplers.
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

muSample = mcmcChain[, "mu" ]
tauSample = mcmcChain[, "tau" ]
sigmaSample <- 1 / sqrt( tauSample ) # Convert precision to SD

# specify preferred graph-file type:
graphType="eps" # or "jpg" or whatever you want

# Marginal distributions of mu and sigma:
openGraph(width=7,height=3.5)
layout(matrix(1:2,nrow=1))
plotPost( muSample , xlab=bquote(mu) )
plotPost( sigmaSample , xlab=bquote(sigma) , showMode=TRUE )
saveGraph(file=paste(fileNameRoot,"PostMuSigmaMarg",sep=""),type=graphType)

# Joint distribution of mu and sigma:
nPts = length(muSample) ; nPtsForDisplay = min( nPts , 2000 )
thinIdx = seq( 1 , nPts , nPts / nPtsForDisplay )
openGraph(width=5,height=5)
plot( muSample[thinIdx] , sigmaSample[thinIdx] , col="skyblue" ,
      xlab=bquote(mu) , ylab=bquote(sigma) , cex.lab=1.5 )
points( mean(muSample) , mean(sigmaSample) , pch="+" , cex=2 )
text( mean(muSample) , mean(sigmaSample) ,
      bquote( .(round(mean(muSample),1)) *"  "* .(round(mean(sigmaSample),1)) ),
      adj=c(.5,-0.5) )
saveGraph(file=paste(fileNameRoot,"PostMuSigmaJoint",sep=""),type=graphType)

# Posterior predictive check:
openGraph(width=5,height=4)
histInfo = hist( y , xlab="y (data)" , 
                 main="Data with Posterior Pred. Distrib." , 
                 breaks=20 , col="grey" , border="white" , prob=TRUE )
yLim = range( histInfo$breaks )
yComb = seq( yLim[1] , yLim[2] , length=501 )
chainLength = length(muSample)
chainIdxVec = floor( seq(1,chainLength,length=20) )
for ( i in chainIdxVec ) {
  lines( yComb , dnorm( yComb , mean=muSample[i] , sd=sigmaSample[i] ) ,
         col="skyblue" )
}
saveGraph(file=paste(fileNameRoot,"PostPredict",sep=""),type=graphType)



#------------------------------------------------------------------------------