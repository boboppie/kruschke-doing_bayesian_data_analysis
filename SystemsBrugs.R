graphics.off()
rm(list=ls(all=TRUE))
library(BRugs)         # Kruschke, J. K. (2010). Doing Bayesian data analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.
modelstring = "
# BUGS model specification begins here...
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
# ... end BUGS model specification
" # close quote for modelstring
writeLines(modelstring,con="model.txt")
modelCheck( "model.txt" )

#------------------------------------------------------------------------------
# THE DATA.

# Load the aircraft data:
load( "Systems.Rdata" ) # loads dataMat
Nsubj = length( unique( dataMat[,"Aircraft"] ) )
# Transform the data:
DaysTransf = dataMat[,"Days"]^(1/5)
dataMat = cbind( dataMat , DaysTransf )
colnames( dataMat ) = c( colnames( dataMat )[1:3] , "DaysTransf" )
# Put it into generic variables so easier to change data in other applications:
y = dataMat[,"DaysTransf"] 
subj = dataMat[,"Aircraft"] 
Ndata = NROW(dataMat) 

# Specify data, as a list.
datalist = list(
    y = y ,
    subj = subj ,
    Ndata = Ndata ,
    Nsubj = Nsubj
)
# Get the data into BRugs: (default filename is data.txt).
modelData( bugsData( datalist ) )

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

# First, compile the model:
nchain = 10
modelCompile( numChains = nchain )

modelGenInits() # works when the priors are not too flat

#------------------------------------------------------------------------------
# RUN THE CHAINS

# burn in
BurnInSteps = 1000
modelUpdate( BurnInSteps )
# actual samples
samplesSet( c( "muG" , "tauG" , "mu" , "tau" , "m" , "d" ) )
stepsPerChain = ceiling(10000/nchain)
thinStep = 100
modelUpdate( stepsPerChain , thin=thinStep )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

source("plotChains.R")
source("plotPost.R")
fileNameRoot = "SystemsBrugs"

# Examine chains for convergence and autocorrelation:
muSum = plotChains( "muG" , saveplots=F , filenameroot=fileNameRoot )
tauSum = plotChains( "tauG" , saveplots=F , filenameroot=fileNameRoot )
mSum = plotChains( "m" , saveplots=F , filenameroot=fileNameRoot )
dSum = plotChains( "d" , saveplots=F , filenameroot=fileNameRoot )
mu1Sum = plotChains( "mu[1]" , saveplots=F , filenameroot=fileNameRoot )
tau1Sum = plotChains( "tau[1]" , saveplots=F , filenameroot=fileNameRoot )

# Extract chains from BUGS into R:
muGsample = samplesSample( "muG" )
tauGsample = samplesSample( "tauG" )
mSample = samplesSample( "m" )
dSample = samplesSample( "d" )
muSample = NULL
tauSample = NULL
for ( sIdx in 1:Nsubj ) {
    muSample = rbind( muSample , samplesSample( paste("mu[",sIdx,"]",sep="") ) )
    tauSample = rbind( tauSample , samplesSample( paste("tau[",sIdx,"]",sep="") ) )
}

# Plot the aircraft mu:
windows(15,6)
layout( matrix( 1:Nsubj , nrow=2 , byrow=T ) )
for ( sIdx in 1:Nsubj ) {
    histInfo = plotPost( muSample[sIdx,] , xlab=bquote(mu[.(sIdx)]) )
}
dev.copy2eps(file=paste(fileNameRoot,"PostMu.eps",sep=""))

# Plot the aircraft tau:
windows(15,6)
layout( matrix( 1:Nsubj , nrow=2 , byrow=T ) )
for ( sIdx in 1:Nsubj ) {
    histInfo = plotPost( tauSample[sIdx,] , xlab=bquote(tau[.(sIdx)]) , 
						HDItextPlace=.3 )
}
dev.copy2eps(file=paste(fileNameRoot,"PostTau.eps",sep=""))

# Plot the hyperdistributions:
windows(15,3)
layout( matrix(1:4,ncol=4) )
histInfo = plotPost( muGsample , xlab=expression(mu[G]) , breaks=30 )
histInfo = plotPost( tauGsample , xlab=expression(tau[G]) , breaks=30 )
histInfo = plotPost( mSample , xlab=expression(m) , breaks=30 )
histInfo = plotPost( dSample , xlab=expression(d) , breaks=30 , HDItextPlace=.1 )
dev.copy2eps(file=paste(fileNameRoot,"PostHyper.eps",sep=""))

#------------------------------------------------------------------------------