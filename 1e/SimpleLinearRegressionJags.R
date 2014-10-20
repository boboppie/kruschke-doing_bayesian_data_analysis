graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="SimpleLinearRegressionJags" # for constructing output filenames
source("openGraphSaveGraph.R")
source("plotPost.R")
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.
modelstring = "
model {
    for( i in 1 : Ndata ) {
        y[i] ~ dnorm( mu[i] , tau )
        mu[i] <- beta0 + beta1 * x[i]
    }
    beta0 ~ dnorm( 0 , 1.0E-12 )
    beta1 ~ dnorm( 0 , 1.0E-12 )
    tau ~ dgamma( 0.001 , 0.001 )
}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.

# Simulated height and weight data:
source("HtWtDataGenerator.R")
HtWtData = HtWtDataGenerator( 30 , rndsd=5678 )
nSubj = NROW( HtWtData )
x = HtWtData[,"height"]
y = HtWtData[,"weight"]

# Re-center data at mean, to reduce autocorrelation in MCMC sampling.
# Standardize (divide by SD) to make initialization easier.
xM = mean( x ) ; xSD = sd( x )
yM = mean( y ) ; ySD = sd( y )
zx = ( x - xM ) / xSD
zy = ( y - yM ) / ySD

# Specify data, as a list.
dataList = list(
  x = zx ,
  y = zy ,
  Ndata = nSubj
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

r = cor(x,y)
initsList = list(
    beta0 = 0 ,    # because data are standardized
    beta1 = r ,        # because data are standardized
    tau = 1 / (1-r^2)  # because data are standardized
)

#------------------------------------------------------------------------------
# RUN THE CHAINS

parameters = c("beta0" , "beta1" , "tau")  # The parameter(s) to be monitored.
adaptSteps = 500              # Number of steps to "tune" the samplers.
burnInSteps = 500            # Number of steps to "burn-in" the samplers.
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

# Extract chain values:
z0 = mcmcChain[, "beta0" ]
z1 = mcmcChain[, "beta1" ]
zTau = mcmcChain[, "tau" ]
zSigma = 1 / sqrt( zTau ) # Convert precision to SD

# Convert to original scale:
b1 = z1 * ySD / xSD
b0 = ( z0 * ySD + yM - z1 * ySD * xM / xSD )
sigma = zSigma * ySD

# Posterior prediction:
# Specify x values for which predicted y's are needed:
xPostPred = seq(55,80,1)
# Define matrix for recording posterior predicted y values at each x value.
# One row per x value, with each row holding random predicted y values.
postSampSize = length(b1)
yPostPred = matrix( 0 , nrow=length(xPostPred) , ncol=postSampSize )
# Define matrix for recording HDI limits of posterior predicted y values:
yHDIlim = matrix( 0 , nrow=length(xPostPred) , ncol=2 )
# Generate posterior predicted y values.
# This gets only one y value, at each x, for each step in the chain.
for ( chainIdx in 1:postSampSize ) {
    yPostPred[,chainIdx] = rnorm( length(xPostPred) ,
                           mean = b0[chainIdx] + b1[chainIdx] * xPostPred ,
                           sd = rep( sigma[chainIdx] , length(xPostPred) ) )
}
source("HDIofMCMC.R")
for ( xIdx in 1:length(xPostPred) ) {
    yHDIlim[xIdx,] = HDIofMCMC( yPostPred[xIdx,] )
}

# Display believable beta0 and b1 values
openGraph(10,5)
par( mar=c(4,4,1,1)+0.1 , mgp=c(2.5,0.8,0) )
layout( matrix(1:2,nrow=1) )
thinIdx = seq(1,length(b0),length=700)
plot( z1[thinIdx] , z0[thinIdx] , cex.lab=1.75 ,
      ylab="Standardized Intercept" , xlab="Standardized Slope" ,
      col="skyblue")
plot( b1[thinIdx] , b0[thinIdx] , cex.lab=1.75 ,
      ylab="Intercept (ht when wt=0)" , xlab="Slope (pounds per inch)" ,
      col="skyblue")
saveGraph(file=paste(fileNameRoot,"SlopeIntercept",sep=""),type="eps")

# Display the posterior of the b1:
openGraph(10,4)
par( mar=c(4,4,1,1)+0.1 , mgp=c(2.5,0.8,0) )
layout( matrix(1:2,nrow=1) )
histInfo = plotPost( z1 , xlab="Standardized slope" , compVal=0.0 )
histInfo = plotPost( b1 , xlab="Slope (pounds per inch)" , compVal=0.0 )
saveGraph(file=paste(fileNameRoot,"PostSlope",sep=""),type="eps")

# Display data with believable regression lines and posterior predictions.
openGraph()
par( mar=c(3,3,2,1)+0.5 , mgp=c(2.1,0.8,0) )
# Plot data values:
xRang = max(x)-min(x)
yRang = max(y)-min(y)
limMult = 0.25
xLim= c( min(x)-limMult*xRang , max(x)+limMult*xRang )
yLim= c( min(y)-limMult*yRang , max(y)+limMult*yRang )
plot( x , y , cex=1.5 , lwd=2 , col="black" , xlim=xLim , ylim=yLim ,
      xlab="X (height in inches)" , ylab="Y (weight in pounds)" , cex.lab=1.5 ,
      main="Data with credible regression lines" , cex.main=1.33  )
# Superimpose a smattering of believable regression lines:
for ( i in seq(from=1,to=length(b0),length=50) ) {
    abline( b0[i] , b1[i] , col="skyblue" )
}
saveGraph(file=paste(fileNameRoot,"DataLines",sep=""),type="eps")

# Display data with HDIs of posterior predictions.
openGraph()
par( mar=c(3,3,2,1)+0.5 , mgp=c(2.1,0.8,0) )
# Plot data values:
yLim= c( min(yHDIlim) , max(yHDIlim) )
plot( x , y , cex=1.5 , lwd=2 , col="black" , xlim=xLim , ylim=yLim ,
      xlab="X (height in inches)" , ylab="Y (weight in pounds)" , cex.lab=1.5 ,
      main="Data with 95% HDI & Mean of Posterior Predictions" , cex.main=1.33  )
# Superimpose posterior predicted 95% HDIs:
segments( xPostPred, yHDIlim[,1] , xPostPred, yHDIlim[,2] , lwd=3, col="skyblue" )
points( xPostPred , rowMeans( yPostPred ) , pch="+" , cex=2 , col="skyblue" )
saveGraph(file=paste(fileNameRoot,"DataPred",sep=""),type="eps")

#------------------------------------------------------------------------------
