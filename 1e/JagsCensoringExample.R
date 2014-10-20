# This program illustrates estimation of parameters for right-censored data in 
# JAGS.  To produce the graphics at the end, you must have plotPost.R and 
# HDIofMCMC.R in the same folder as this program (and tell R that this is the 
# working directory). Run the program as-is to see that JAGS produces good 
# estimates of the true underlying parameter values, even for censored data. 

# To see what happens when you do *not* tell JAGS that the data are censored, 
# just comment out the follow four lines (and run it again).
# In the model specification:
#   isCensored[i] ~ dinterval( y[i] , censorLimitVec[i] )
# In the dataList:
#   , isCensored = as.numeric(isCensored) # JAGS dinterval needs 0,1 not T,F
#   , censorLimitVec = censorLimitVec
# In the initsList:
#   , y=yInit 
# Notice that the resulting estimated parameters are way too low!

graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="JagsCensoringExample" # for constructing output filenames
source("openGraphSaveGraph.R")
source("plotPost.R")
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.

#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
model {
  for ( i in 1:N ) {
    isCensored[i] ~ dinterval( y[i] , censorLimitVec[i] )
    y[i] ~ dnorm( mu , tau ) 
  }
  tau <- 1/pow(sigma,2)
  sigma ~ dunif(0,100)
  mu ~ dnorm(0,1E-6)
}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.
# For real data, censored values should be recorded as NA, not as the censoring 
# limit, and the indicator value isCensored should be set to TRUE.

# Generate some random data:
set.seed(47405)
y = rnorm( n=500 )
desiredMean = 100 ; desiredSD = 15
actualMean = mean(y) ; actualSD = sd(y)
y = desiredSD*(y-actualMean)/actualSD + desiredMean
N = length(y)
# Now artificially right-censor the data. 
# First, construct the censoring limits. There is a censoring limit defined
# for every observation, put in the vector censorLimitVec.
censorLimit = quantile( y , probs=c(pnorm(1)) ) # cuts at +1SD
censorLimitVec = rep( censorLimit , length(y) )
# Define indicator vector for which observations are censored:
isCensored = ( y > censorLimitVec ) # Note that this has T,F logical values.
# Now cut the data above the censor limits. For real data, censored values
# should be recorded as NA, not as the censor limit, and the indicator value
# isCensored should be set to TRUE.
y[isCensored] = NA

# Specify the data in a form that is compatible with model, as a list:
dataList = list( y = y , N = N 
                 , isCensored = as.numeric(isCensored) # JAGS dinterval needs 0,1 not T,F
                 , censorLimitVec = censorLimitVec
                 )

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

sigmaInit = 15
muInit = 100
# intial values of censored data:
yInit = rep( NA , length(y) )
yInit[isCensored] = censorLimitVec[isCensored]+1

initsList = list( mu=muInit , sigma = sigmaInit 
                  , y=yInit 
                  )

#------------------------------------------------------------------------------
# RUN THE CHAINS

parameters = c( "mu" , "sigma" )
adaptSteps = 1000 # overkill, just shows general form
burnInSteps = 2000 # overkill, just shows general form
nChains = 3 
numSavedSteps=50000
thinSteps=1
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
# Create, initialize, and adapt the model:
jagsModel = jags.model( "model.txt" , data=dataList , inits=initsList , 
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( jagsModel , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                            n.iter=nIter , thin=thinSteps )
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

openGraph(width=10,height=3.0)
layout(matrix(1:3,nrow=1))
par( mar=c(3.5,1.5,3.0,0.5) , mgp=c(2,0.7,0) )
# plot the data:
xLim = desiredMean+c(-3.5*desiredSD,3.5*desiredSD)
xBreaks = seq( xLim[1] , xLim[2] , desiredSD/2 )
xGrid = seq( xLim[1] , xLim[2] , length=201 )  
histInfo = hist( y , main="Data" , xlab="y" , freq=F ,
                 xlim=xLim , breaks=c(min(y,na.rm=T),xBreaks) , 
                 col="grey" , border="white" )
lines( xGrid , dnorm( xGrid , mean=desiredMean , sd=desiredSD )/pnorm(1) , lwd=3 ) 
abline(v=censorLimit,lty="dotted",col="red")
# plot the posterior:
histInfo = plotPost( mcmcChain[,"mu"] , main="Mu" , xlab=expression(mu) , 
                     showMode=FALSE )
histInfo = plotPost( mcmcChain[,"sigma"] , main="Sigma" , xlab=expression(sigma) , 
                     showMode=TRUE )
# save plot:
saveGraph( file=fileNameRoot , type="eps" )
saveGraph( file=fileNameRoot , type="jpg" )
  