graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="OrdinalProbitRegressionJags" # for constructing output filenames
source("openGraphSaveGraph.R")
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
model {
    for( i in 1 : nData ) {
        y[i] ~ dcat( pr[i,1:nYlevels] )
        pr[i,1] <- phi( ( thresh[1] - mu[i] ) / sigma )
        for ( k in 2:(nYlevels-1) ) {
            pr[i,k] <- max( 0 ,  phi( ( thresh[ k ] - mu[i] ) / sigma )
                               - phi( ( thresh[k-1] - mu[i] ) / sigma ) )
        }
        pr[i,nYlevels] <- 1 - phi( ( thresh[nYlevels-1] - mu[i] ) / sigma )
        mu[i] <- b0 + inprod( b[1:nPredictors] , x[i,1:nPredictors] )
    }
    bPrec <- pow( nYlevels/4 , -2 ) # max plausible slope is 1SD
    for ( j in 1:nPredictors ) {
        b[j] ~ dnorm(0,bPrec) # modest precision because of normalized x,y values
    }
    threshPriorPrec <- 1
    for ( k in 1:(nYlevels-1) ) {
        threshPriorMean[k] <- k+0.5
        thresh[k] ~ dnorm( threshPriorMean[k] , threshPriorPrec )
    }
}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.

dataSource = c("Random","Movies")[2]

# The loading of data must produce a matrix called dataMat that has
# one row per datum, where the first column is the ordinal predicted value
# and the 2nd - last columns are the predictor values. The columns should
# be named.

if ( dataSource=="Random" ) {
   fileNameRoot = paste( fileNameRoot , dataSource , sep="" )
   # Generate some random toy data.
   source( "OrdinalProbitDataGenerator.R" )
   nYlevels = 7
   dataMat = OrdinalProbitDataGenerator( nData = 200 ,
              normPrec=200 , slope=c(-1.0,1.26) , # c(-1.0,1.26) matches Movies
              thresh=c(-Inf,seq(-1.2,1.2,length=nYlevels-1),Inf) ,
              nYlevels=nYlevels , makePlots=F , rndSeed=47405 )
   # Change x values to arbitrary non-standardized scales:
   dataMat[,2] =  1963.64 + 18.13 * dataMat[,2]
   dataMat[,3] =  92.87 + 18.26 * dataMat[,3]
}

if ( dataSource=="Movies" ) {
   fileNameRoot = paste( fileNameRoot , dataSource , sep="" )
   dataFram = read.table( "Moore2006data.txt" , header=T )
   rateVals = sort( unique( dataFram[,"Rating"] ) )
   rankVals = match( dataFram[,"Rating"] , rateVals ) # convert to ranks
   dataMat = cbind( rankVals , dataFram[,"Year"] , dataFram[,"Length"] )
   colnames(dataMat) = c("Rating","Year","Length")
}

# Rename for use by generic processing later:
nData = NROW(dataMat)
x = dataMat[,-1]
predictorNames = colnames(dataMat)[-1]
nPredictors = NCOL(x)
y = as.matrix(dataMat[,1])
predictedName = colnames(dataMat)[1]
nYlevels = max(y)

# Re-center x values at mean, to reduce autocorrelation in MCMC sampling.
# Standardize (divide by SD) to make prior-setting easier.
standardizeCols = function( dataMat ) {
    zDataMat = dataMat
    for ( colIdx in 1:NCOL( dataMat ) ) {
        mCol = mean( dataMat[,colIdx] )
        sdCol = sd( dataMat[,colIdx] )
        zDataMat[,colIdx] = ( dataMat[,colIdx] - mCol ) / sdCol
    }
    return( zDataMat )
}
zx = standardizeCols( x )
# Don't standarize y because they must be integers, 1 to nYlevels

lmInfo = lm( y ~ zx ) # R function returns MLE
b0Init = lmInfo$coef[1]
bInit = lmInfo$coef[-1]
sigmaInit = sqrt(sum(lmInfo$res^2)/nData)

dataList = list(
           x = zx ,
           y = as.vector( y ) , # BUGS does not treat 1-column mat as vector
           nPredictors = nPredictors ,
           nData = nData ,
           nYlevels = nYlevels ,
           sigma = sigmaInit , # fixed, not estimated
           b0 = b0Init  # fixed, not estimated
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

initsList = list(
        b = bInit , # from lm(y~zx), above
        thresh = 1:(nYlevels-1)+.5
    )

#------------------------------------------------------------------------------
# RUN THE CHAINS

parameters = c( "b" , "thresh" )  
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

checkConvergence = F
if ( checkConvergence ) {
  show( summary( codaSamples ) )
  openGraph()
  plot( codaSamples , ask=F )  
  openGraph()
  autocorr.plot( codaSamples , ask=F )
}

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )


# Extract parameter values:
zbSamp = NULL
for ( j in 1:nPredictors ) {
   zbSamp = cbind( zbSamp , mcmcChain[, paste("b[",j,"]",sep="") ] )
}
chainLength = NROW(zbSamp)
zthreshSamp = NULL
for ( j in 1:(nYlevels-1) ) {
   zthreshSamp = cbind( zthreshSamp ,
                        mcmcChain[,paste("thresh[",j,"]",sep="")] )
}

# Convert to original scale:
bSamp = zbSamp * matrix( 1/(sigmaInit*apply(x,2,sd)) , byrow=TRUE ,
                         ncol=nPredictors , nrow=chainLength )
threshSamp = (1/sigmaInit) * ( zthreshSamp - b0Init +
                rowSums( zbSamp * matrix( apply(x,2,mean)/apply(x,2,sd) ,
                                          byrow=TRUE , ncol=nPredictors ,
                                          nrow=chainLength ) ) )
b0 = 0
sigma = 1

source("plotPost.R")

# Scatter plots of parameter values, pairwise:
if ( (nPredictors+nYlevels) <= 10 ) { # don't display if too many
    openGraph()
    thinIdx = ceiling(seq(1,chainLength,length=200))
    pairs( cbind( zbSamp[thinIdx,] , zthreshSamp[thinIdx,] )  ,
           labels=c( paste("zb",predictorNames,sep="") ,
                     paste("zthresh",1:nYlevels,sep="")) )
    openGraph()
    pairs( cbind( bSamp[thinIdx,] , threshSamp[thinIdx,] )  ,
           labels=c( paste("b",predictorNames,sep="") ,
                     paste("thresh",1:nYlevels,sep="")) )
    saveGraph(file=paste(fileNameRoot,"PostPairs.eps",sep=""),type="eps")
}

# Display the posterior:
nPlotPerRow = 5
nPlotRow = ceiling((nPredictors+nYlevels-1)/nPlotPerRow)
nPlotCol = ceiling((nPredictors+nYlevels-1)/nPlotRow)
openGraph(3.5*nPlotCol,2.25*nPlotRow)
layout( matrix(1:(nPlotRow*nPlotCol),nrow=nPlotRow,ncol=nPlotCol,byrow=T) )
par( mar=c(4,3,2.5,0) , mgp=c(2,0.7,0) )
for ( sIdx in 1:nPredictors ) {
histInfo = plotPost( bSamp[,sIdx] , xlab="Slope Value" , compVal=0.0 ,
                     breaks=30 ,
                     main=bquote( b *.(predictorNames[sIdx]) ) ,
                     cex.main=1.67 , cex.lab=1.33 )
}
for ( sIdx in 1:(nYlevels-1) ) {
histInfo = plotPost( threshSamp[,sIdx] , xlab="Thresh Value" , compVal=NULL ,
                     breaks=30 ,
                     main=bquote( theta * .(sIdx) ) ,
                     cex.main=1.67 , cex.lab=1.33 )
}
saveGraph(file=paste(fileNameRoot,"PostHist.eps",sep=""),type="eps")

# Plot the data
if ( nPredictors == 2 ) {
openGraph()
plot( x[,1] , x[,2] , xlab=colnames(x)[1] , ylab=colnames(x)[2] ,
      main=paste( "The Data (" , dataSource , ")" , sep="") ,
      pch=as.character(y) )
for ( chainIdx in round(seq(1,chainLength,len=3)) ) {
  for ( threshIdx in 1:(nYlevels-1) ) {
    abline( threshSamp[chainIdx,threshIdx]/bSamp[chainIdx,2] ,
            -bSamp[chainIdx,1]/bSamp[chainIdx,2] ,
            lwd = 2 , lty=chainIdx , col="grey" )
  }
}
saveGraph(file=paste(fileNameRoot,"Data.eps",sep=""),type="eps")

} # end if nPredictors == 2

# Posterior prediction.
xProbe = c( 1991 , 94 ) # Note order of values: x1 is year and x2 is duration.
# Set up a matrix for storing the values of p(y|xProbe) at each step in chain.
py = matrix( 0 , nrow=chainLength , ncol=nYlevels )
# Step through chain and compute p(y|xProbe) and each step:
for ( chainIdx in 1:chainLength ) {
    yValue = 1
    py[chainIdx,yValue] = (
        pnorm( threshSamp[chainIdx,yValue]
               - sum( bSamp[chainIdx,] * xProbe ) ) )
    for ( yValue in 2:(nYlevels-1) ) {
        py[chainIdx,yValue] = (
            pnorm( threshSamp[chainIdx,yValue]
                   - sum( bSamp[chainIdx,] * xProbe ) )
             - pnorm( threshSamp[chainIdx,yValue-1]
                      - sum( bSamp[chainIdx,] * xProbe ) ) )
    }
    yValue = nYlevels
    py[chainIdx,yValue] = ( 1 -
        pnorm( threshSamp[chainIdx,yValue-1]
               - sum( bSamp[chainIdx,] * xProbe ) ) )
}
# Now average across the chain:
pyAve = colMeans( py )
