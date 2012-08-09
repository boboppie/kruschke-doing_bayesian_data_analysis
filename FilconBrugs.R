graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="FilconBrugs" # for constructing output filenames
library(BRugs)         # Kruschke, J. K. (2010). Doing Bayesian data analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
# BUGS model specification begins here...
model {
   for ( subjIdx in 1:nSubj ) {
      # Likelihood:
      z[subjIdx] ~ dbin( theta[subjIdx] , N[subjIdx] )
      # Prior on theta: Notice nested indexing.
      theta[subjIdx] ~ dbeta( a[cond[subjIdx]] , b[cond[subjIdx]] )I(0.001,0.999)
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
   meanGamma <- 10
   sdGamma <- 10
}
# ... end BUGS model specification
" # close quote for modelstring
# Write model to a file:
writeLines(modelstring,con="model.txt")
# Load model file into BRugs and check its syntax:
modelCheck( "model.txt" )

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

# Specify the data in a form that is compatible with BRugs model, as a list:
datalist = list(
 nCond = nCond ,
 nSubj = nSubj ,
 cond = cond ,
 N = N ,
 z = z
)

# Get the data into BRugs:
# Function bugsData stores the data file (default filename is data.txt).
# Function modelData loads data file into BRugs (default filename is data.txt).
modelData( bugsData( datalist ) )

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

nChain = 3
modelCompile( numChains=nChain )

if ( F ) {
   modelGenInits() # often won't work for diffuse prior
} else {
  #  initialization based on data
  genInitList <- function() {
    sqzData = .01+.98*datalist$z/datalist$N
    mu = aggregate( sqzData , list(datalist$cond) , "mean" )[,"x"]
    sd = aggregate( sqzData , list(datalist$cond) , "sd" )[,"x"]
    kappa = mu*(1-mu)/sd^2 - 1
    return(
      list(
        theta = sqzData ,
        mu = mu ,
        kappa = kappa
      )
    )
  }
  for ( chainIdx in 1 : nChain ) {
    modelInits( bugsInits( genInitList ) )
  }
}

#------------------------------------------------------------------------------
# RUN THE CHAINS.

burninSteps = 2000
modelUpdate( burninSteps )
samplesSet( c("mu","kappa","theta","a","b") )
nPerChain = ceiling(5000/nChain)
modelUpdate( nPerChain , thin=10 )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS.

# Check for convergence, mixing and autocorrelation:
source("plotChains.R")
sumInfo = plotChains( "mu" , saveplots=T , fileNameRoot )
sumInfo = plotChains( "kappa" , saveplots=F )
sumInfo = plotChains( "theta[1]" , saveplots=F )

# Extract parameter values and save them.
mu = NULL
kappa = NULL
for ( condIdx in 1:nCond ) {
   mu = rbind( mu , samplesSample( paste("mu[",condIdx,"]",sep="") ) )
   kappa = rbind( kappa , samplesSample( paste("kappa[",condIdx,"]",sep="") ) )
}
save( mu , kappa , file=paste(fileNameRoot,"MuKappa.Rdata",sep="") )
chainLength = NCOL(mu)

# Histograms of mu differences:
windows(10,2.75)
layout( matrix(1:3,nrow=1) )
source("plotPost.R")
plotPost( mu[1,]-mu[2,] , xlab=expression(mu[1]-mu[2]) , main="" ,
          breaks=20 , compVal=0 )
plotPost( mu[3,]-mu[4,] , xlab=expression(mu[3]-mu[4]) , main="" ,
          breaks=20 , compVal=0 )
plotPost( (mu[1,]+mu[2,])/2 - (mu[3,]+mu[4,])/2 ,
          xlab=expression( (mu[1]+mu[2])/2 - (mu[3]+mu[4])/2 ) ,
          main="" , breaks=20 , compVal=0 )
dev.copy2eps(file=paste(fileNameRoot,"MuDiffs.eps",sep=""))

# Scatterplot of mu,kappa in each condition:
windows()
muLim = c(.60,1) ; kappaLim = c( 4.0 , 40 ) ; mainLab="Posterior"
thindex = round( seq( 1 , chainLength , len=300 ) )
plot( mu[1,thindex] , kappa[1,thindex] , main=mainLab ,
      xlab=expression(mu[c]) , ylab=expression(kappa[c]) , cex.lab=1.75 ,
      xlim=muLim , ylim=kappaLim , log="y" , col="red" , pch="1" )
points( mu[2,thindex] , kappa[2,thindex] , col="blue" , pch="2" )
points( mu[3,thindex] , kappa[3,thindex] , col="green" , pch="3" )
points( mu[4,thindex] , kappa[4,thindex] , col="black" , pch="4" )
dev.copy2eps(file=paste(fileNameRoot,"Scatter.eps",sep=""))
