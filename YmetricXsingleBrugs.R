graphics.off()
rm(list=ls(all=TRUE))
library(BRugs)         # Kruschke, J. K. (2010). Doing Bayesian data analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.
modelstring = "
# BUGS model specification begins here...
model {
    # Likelihood:
    for( i in 1 : N ) {
        y[i] ~ dnorm( mu , tau ) # tau is precision, not SD
    }
    # Prior:
    tau ~ dgamma( 0.01 , 0.01 )
    mu ~ dnorm( 0 , 1.0E-10 )
}
# ... end BUGS model specification
" # close quote for modelstring
writeLines(modelstring,con="model.txt")
modelCheck( "model.txt" )

#------------------------------------------------------------------------------
# THE DATA.

# Generate random data from known parameter values:
set.seed(47405)
trueM = 100
trueSD = 15
y = round( rnorm( n=500 , mean=trueM , sd=trueSD ) ) # R dnorm uses mean and SD

datalist = list(
    y = y ,
    N = length( y )
)

# Get the data into BRugs: (default filename is data.txt).
modelData( bugsData( datalist ) )

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

nchain = 3
modelCompile( numChains = nchain )

automaticInit = F # TRUE or FALSE
if ( automaticInit ) {
    modelGenInits()  # automatically initialize chains from prior
} else {
    genInitList <- function() { # manually initialize chains near the data
        list( mu = mean( datalist$y ) ,
              tau = 1 / sd( datalist$y )^2 )
    }
    for ( chainIdx in 1 : nchain ) {
        modelInits( bugsInits( genInitList ) )
    }
}

#------------------------------------------------------------------------------
# RUN THE CHAINS

# burn in
BurnInSteps = 500
modelUpdate( BurnInSteps )
# actual samples
samplesSet( c( "mu" , "tau" ) )
stepsPerChain = 2000
thinStep = 1
modelUpdate( stepsPerChain , thin=thinStep )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

filenamert = "YmetricXsingleBrugs"

source("plotChains.R")
muSum = plotChains( "mu" , saveplots=F , filenamert )
sigmaSum = plotChains( "tau" , saveplots=F , filenamert )

muSample = samplesSample( "mu" )
tauSample = samplesSample( "tau" )
sigmaSample <- 1 / sqrt( tauSample ) # Convert precision to SD

source("plotPost.R")
windows()
plotPost( muSample , xlab="mu" , breaks=30 , main="Posterior" )
dev.copy2eps(file=paste(filenamert,"PostMu.eps",sep=""))

nPts = length(muSample) ; nPtsForDisplay = min( nPts , 2000 )
thinIdx = seq( 1 , nPts , nPts / nPtsForDisplay )
windows()
plot( muSample[thinIdx] , sigmaSample[thinIdx] , col="gray" ,
      xlab="mu" , ylab="sigma" , cex.lab=1.5 , main="Posterior" , log="y" )
points( mean(muSample) , mean(sigmaSample) , pch="+" , cex=2 )
text( mean(muSample) , mean(sigmaSample) ,
      bquote( .(round(mean(muSample),1)) *"  "* .(round(mean(sigmaSample),1)) ),
      adj=c(.5,-0.5) )
dev.copy2eps(file=paste(filenamert,"PostMuSigma.eps",sep=""))

#------------------------------------------------------------------------------