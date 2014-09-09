# Example of modifying a program to include a new parameter. This programs
# accompanies a blog entry at http://doingbayesiandataanalysis.blogspot.com/
graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="SimpleLinearRegressionJags-WithQuadTrend" # for output filenames
source("openGraphSaveGraph.R")
source("plotPost.R")   # Adapted by John Kruschke from programs accompanying
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.
modelstring = "
model {
    for( i in 1 : Ndata ) {
        y[i] ~ dnorm( mu[i] , tau )
        mu[i] <- beta0 + beta1 * x[i] + beta2 * pow(x[i],2) # CHANGED
    }
    beta0 ~ dnorm( 0 , 1.0E-12 )
    beta1 ~ dnorm( 0 , 1.0E-12 )
    beta2 ~ dnorm( 0 , 1.0E-12 ) # CHANGED
    tau ~ dgamma( 0.001 , 0.001 )
}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.

# Generate some random data with a quadratic trend:
set.seed(47405)
N=50
x = rnorm( N , mean=0 , sd=3 ) 
y = rnorm( N )*0.4 + ( -4 + 2.6*(x-mean(x)) + 0.23*(x-mean(x))^2 )
x=x+4
y=y+15
Ndata = length(y)

# Specify data, as a list.
dataList = list(
  x = x ,
  y = y ,
  Ndata = Ndata
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

# Use R's built-in least-squares regression to get plausible initial values:
lmInfo = lm( dataList$y ~ dataList$x ) 
b0Init = lmInfo$coef[1]
bInit = lmInfo$coef[2]
tauInit = length(dataList$y) / sum(lmInfo$res^2)
initsList = list(
  beta0 = b0Init ,
  beta1 = bInit ,
  beta2 = 0 , # CHANGED
  tau = tauInit
)

#------------------------------------------------------------------------------
# RUN THE CHAINS

parameters = c("beta0" , "beta1" , "beta2" , "tau")  # CHANGED
adaptSteps = 1000             # Number of steps to "tune" the samplers.
burnInSteps = 2000            # Number of steps to "burn-in" the samplers.
nChains = 3                   # Number of chains to run.
numSavedSteps=100000           # Total number of steps in chains to save.
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
chainLength = NROW(mcmcChain)

# For convenience later, append a column with tau converted to sigma:
sigma = 1 / sqrt( mcmcChain[, "tau" ] ) # Convert precision to SD
mcmcChain = cbind( mcmcChain , sigma )

# Specify preferred file type for saved graphs:
graphFileType = "jpg"

# Posterior prediction:
# Specify x values for which predicted y's are wanted:
extrapolationExtent = 0.25*(range(x)[2]-range(x)[1])
lowX = range(x)[1] - extrapolationExtent
highX = range(x)[2] + extrapolationExtent
xPostPred = seq( lowX , highX , length=20 )
# Define matrix for recording posterior predicted y values at each x value.
# One row per x value, with each row holding random predicted y values.
yPostPred = matrix( 0 , nrow=length(xPostPred) , ncol=chainLength )
# Define matrix for recording HDI limits of posterior predicted y values:
yHDIlim = matrix( 0 , nrow=length(xPostPred) , ncol=2 )
# Generate posterior predicted y values.
# This gets only one y value, at each x, for each step in the chain.
for ( chainIdx in 1:chainLength ) {
    yPostPred[,chainIdx] = rnorm( length(xPostPred) ,
                           mean = ( mcmcChain[chainIdx,"beta0"] 
                                    + mcmcChain[chainIdx,"beta1"] * xPostPred 
                                    + mcmcChain[chainIdx,"beta2"] * xPostPred^2
                                    # CHANGED
                                    ) ,
                           sd = rep( sigma[chainIdx] , length(xPostPred) ) )
}
source("HDIofMCMC.R")
for ( xIdx in 1:length(xPostPred) ) {
    yHDIlim[xIdx,] = HDIofMCMC( yPostPred[xIdx,] )
}

# Display pair plots of parameter values
openGraph()
par( mar=c(4,4,1,1)+0.1 , mgp=c(2.5,0.8,0) )
thinIdx = seq(1,chainLength,length=700)
pairs( mcmcChain[thinIdx,] , col="skyblue" )
saveGraph(file=paste(fileNameRoot,"Pairs",sep=""),type=graphFileType)

# Display the posterior of the beta1:
openGraph()
par( mar=c(4,4,1,1)+0.1 , mgp=c(2.5,0.8,0) )
histInfo = plotPost( mcmcChain[,"beta1"] , xlab="Slope (beta1)" , compVal=0.0 )
saveGraph(file=paste(fileNameRoot,"PostSlope",sep=""),type=graphFileType)

# CHANGED
# Display the posterior of the beta2:
openGraph()
par( mar=c(4,4,1,1)+0.1 , mgp=c(2.5,0.8,0) )
histInfo = plotPost( mcmcChain[,"beta2"] , xlab="Quadratic Coef. (beta2)" ,
                     compVal=0.0 )
saveGraph(file=paste(fileNameRoot,"PostQuad",sep=""),type=graphFileType)

# Display data with believable regression curves and posterior predictions.
openGraph()
par( mar=c(3,3,2,1)+0.5 , mgp=c(2.1,0.8,0) )
# Plot data values:
xRang = max(x)-min(x)
yRang = max(y)-min(y)
limMult = 0.25
xLim= c( min(x)-limMult*xRang , max(x)+limMult*xRang )
yLim= c( min(y)-limMult*yRang , max(y)+limMult*yRang )
plot( x , y , cex=1.5 , lwd=2 , col="black" , xlim=xLim , ylim=yLim ,
      xlab="X" , ylab="Y" , cex.lab=1.5 ,
      main="Data with credible regression lines" , cex.main=1.33  )
# Superimpose a smattering of believable regression lines:
xComb = seq(xLim[1],xLim[2],length=201)
for ( i in round(seq(from=1,to=chainLength,length=50)) ) {
  lines( xComb , 
         mcmcChain[i,"beta0"] + mcmcChain[i,"beta1"]*xComb 
         + mcmcChain[i,"beta2"]*xComb^2 , # CHANGED 
         col="skyblue" )
}
saveGraph(file=paste(fileNameRoot,"DataPostCurves",sep=""),type=graphFileType)

# Display data with HDIs of posterior predictions.
openGraph()
par( mar=c(3,3,2,1)+0.5 , mgp=c(2.1,0.8,0) )
# Plot data values:
yLim= c( min(yHDIlim) , max(yHDIlim) )
plot( x , y , cex=1.5 , lwd=2 , col="black" , xlim=xLim , ylim=yLim ,
      xlab="X" , ylab="Y" , cex.lab=1.5 ,
      main="Data with 95% HDI & Mean of Posterior Predictions" , cex.main=1.33  )
# Superimpose posterior predicted 95% HDIs:
segments( xPostPred, yHDIlim[,1] , xPostPred, yHDIlim[,2] , lwd=3, col="skyblue" )
points( xPostPred , rowMeans( yPostPred ) , pch=19 , cex=2 , col="skyblue" )
saveGraph(file=paste(fileNameRoot,"DataPostPred",sep=""),type=graphFileType)

#------------------------------------------------------------------------------
