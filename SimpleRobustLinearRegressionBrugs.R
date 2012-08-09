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
        y[i] ~ dt( mu[i] , tau , tdf )
        mu[i] <- beta0 + beta1 * x[i]
    }
    beta0 ~ dnorm( 0 , 1.0E-12 )
    beta1 ~ dnorm( 0 , 1.0E-12 )
    tau ~ dgamma( 0.001 , 0.001 )
    udf ~ dunif(0,1)
    tdf <- 1 - tdfGain * log(1-udf) # tdf in [1,Inf).
    # tdfGain specified in data section
}
# ... end BUGS model specification
" # close quote for modelstring
writeLines(modelstring,con="model.txt")
modelCheck( "model.txt" )

#------------------------------------------------------------------------------
# THE DATA.

cigData = read.csv(file="McIntyre1994data.csv")
nSubj = NROW(cigData)
x = cigData[,"Wt"]
xName="Weight"
y = cigData[,"Tar"]
yName="Tar"

# Re-center data at mean, to reduce autocorrelation in MCMC sampling.
# Standardize (divide by SD) to make initialization easier.
xM = mean( x ) ; xSD = sd( x )
yM = mean( y ) ; ySD = sd( y )
zx = ( x - xM ) / xSD
zy = ( y - yM ) / ySD

# Specify data, as a list.
tdfGain = 1 # 1 for low-baised tdf, 100 for high-biased tdf
datalist = list(
  x = zx ,
  y = zy ,
  Ndata = nSubj ,
  tdfGain = tdfGain
)
# Get the data into BRugs:
modelData( bugsData( datalist ) )

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

nchain = 3
modelCompile( numChains = nchain )

genInitList <- function() {
    r = cor(x,y)
    list(
        beta0 = 0 ,    # because data are standardized
        beta1 = r ,        # because data are standardized
        tau = 1 / (1-r^2) , # because data are standardized
        udf = 0.95 # tdf = 4
    )
}
for ( chainIdx in 1 : nchain ) {
    modelInits( bugsInits( genInitList ) )
}

#------------------------------------------------------------------------------
# RUN THE CHAINS

# burn in
BurnInSteps = 100
modelUpdate( BurnInSteps )
# actual samples
samplesSet( c( "beta0" , "beta1" , "tau" , "tdf" ) )
stepsPerChain = ceiling(50000/nchain)
thinStep = 10
modelUpdate( stepsPerChain , thin=thinStep )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

source("plotChains.R")

fname = paste("SimpleRobustLinearRegressionBrugsTdfGain",tdfGain,
              sep="")
#beta0Sum = plotChains( "beta0" , saveplots=F , filenameroot=fname )
#beta1Sum = plotChains( "beta1"     , saveplots=F , filenameroot=fname )
#tauSum   = plotChains( "tau"       , saveplots=F , filenameroot=fname )
#tdfSum   = plotChains( "tdf"       , saveplots=F , filenameroot=fname )

# Extract chain values:
tdfSamp = samplesSample( "tdf" )
tdfM = mean( tdfSamp )
z0 = samplesSample( "beta0" )
z1 = samplesSample( "beta1" )
zTau = samplesSample( "tau" )
zSigma = 1 / sqrt( zTau ) # Convert precision to SD

# Convert to original scale:
b1 = z1 * ySD / xSD
b0 = ( z0 * ySD + yM - z1 * ySD * xM / xSD )
sigma = zSigma * ySD

# Posterior prediction:
# Specify x values for which predicted y's are needed:
xRang = max(x)-min(x)
yRang = max(y)-min(y)
limMult = 0.25
xLim= c( min(x)-limMult*xRang , max(x)+limMult*xRang )
yLim= c( min(y)-limMult*yRang , max(y)+limMult*yRang )
yLim = c(-10,35)
xPostPred = seq( xLim[1] , xLim[2] , length=20 )
# Define matrix for recording posterior predicted y values at each x value.
# One row per x value, with each row holding random predicted y values.
postSampSize = length(b1)
yPostPred = matrix( 0 , nrow=length(xPostPred) , ncol=postSampSize )
# Define matrix for recording HDI limits of posterior predicted y values:
yHDIlim = matrix( 0 , nrow=length(xPostPred) , ncol=2 )
# Generate posterior predicted y values.
# This gets only one y value, at each x, for each step in the chain.
for ( chainIdx in 1:postSampSize ) {
	# p. 398 of book: dt(x|mu,tau,df) in R is sqrt(tau)*dt((x-mu)*sqrt(tau),df).
	# From R's random number generator rt(n,df), need to multiply by SD 
	# and add mean.
    yPostPred[,chainIdx] = ( rt( length(xPostPred) , df=tdfSamp[chainIdx] )
							* sigma[chainIdx] 
							+ b0[chainIdx] + b1[chainIdx] * xPostPred )
}
source("HDIofMCMC.R")
for ( xIdx in 1:length(xPostPred) ) {
    yHDIlim[xIdx,] = HDIofMCMC( yPostPred[xIdx,] )
}

# Display believable beta0 and b1 values
windows()
par( mar=c(4,4,1,1)+0.1 , mgp=c(2.5,0.8,0) )
#layout( matrix(1:2,nrow=1) )
thinIdx = seq(1,length(b0),length=700)
#plot( z1[thinIdx] , z0[thinIdx] , cex.lab=1.75 ,
#      ylab="Standardized Intercept" , xlab="Standardized Slope" )
plot( b1[thinIdx] , b0[thinIdx] , cex.lab=1.75 ,
      ylab="Intercept" , xlab="Slope" )
dev.copy2eps(file=paste(fname,"SlopeIntercept.eps",sep=""))

# Display the posterior of the b1:
source("plotPost.R")
windows(7,4)
par( mar=c(4,4,1,1)+0.1 , mgp=c(2.5,0.8,0) )
#layout( matrix(1:2,nrow=1) )
#histInfo = plotPost( z1 , xlab="Standardized slope" , compVal=0.0 ,
#                     breaks=30  )
histInfo = plotPost( b1 , main=bquote("Mean tdf"==.(signif(tdfM,3))) , cex.main=2 ,
                     xlab=bquote("Slope (" * Delta * .(yName) / Delta * .(xName)
                                        * ")") , compVal=0.0 , breaks=30  )
dev.copy2eps(file=paste(fname,"PostSlope.eps",sep=""))

# Display data with believable regression lines and posterior predictions.
windows()
par( mar=c(3,3,2,1)+0.5 , mgp=c(2.1,0.8,0) )
# Plot data values:
plot( x , y , cex=1.5 , lwd=2 , col="black" , xlim=xLim , ylim=yLim ,
      xlab=xName , ylab=yName , cex.lab=1.5 ,
      main="Data with credible regression lines" , cex.main=1.33  )
# Superimpose a smattering of believable regression lines:
for ( i in seq(from=1,to=length(b0),length=50) ) {
    abline( b0[i] , b1[i] , col="skyblue" )
}
dev.copy2eps(file=paste(fname,"DataLines.eps",sep=""))

# Display data with HDIs of posterior predictions.
windows()
par( mar=c(3,3,2,1)+0.5 , mgp=c(2.1,0.8,0) )
# Plot data values:
#yLim= c( min(c(yHDIlim,y)) , max(c(yHDIlim,y)) )
plot( x , y , cex=1.5 , lwd=2 , col="black" , xlim=xLim , ylim=yLim ,
      xlab=xName , ylab=yName , cex.lab=1.5 ,
      main="Data with 95% HDI & Mean of Posterior Predictions" , cex.main=1.33  )
# Superimpose posterior predicted 95% HDIs:
segments( xPostPred, yHDIlim[,1] , xPostPred, yHDIlim[,2] , lwd=3, col="skyblue" )
points( xPostPred , apply(yPostPred,1,median) , pch="+" , cex=2 , col="skyblue" )
dev.copy2eps(file=paste(fname,"DataPred.eps",sep=""))

#------------------------------------------------------------------------------
