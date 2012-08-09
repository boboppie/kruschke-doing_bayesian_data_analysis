graphics.off()
rm(list=ls(all=TRUE))
fname = "OrdinalProbitRegressionBrugs"
library(BRugs)        # Kruschke, J. K. (2010). Doing Bayesian data analysis:
                      # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
# BUGS model specification begins here...
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
# ... end BUGS model specification
" # close quote for modelstring
writeLines(modelstring,con="model.txt")
modelCheck( "model.txt" )

#------------------------------------------------------------------------------
# THE DATA.

dataSource = c("Random","Movies")[2]

# The loading of data must produce a matrix called dataMat that has
# one row per datum, where the first column is the ordinal predicted value
# and the 2nd - last columns are the predictor values. The columns should
# be named.

if ( dataSource=="Random" ) {
   fname = paste( fname , dataSource , sep="" )
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
   fname = paste( fname , dataSource , sep="" )
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

# Get the data into BUGS:
datalist = list(
           x = zx ,
           y = as.vector( y ) , # BUGS does not treat 1-column mat as vector
           nPredictors = nPredictors ,
           nData = nData ,
           nYlevels = nYlevels ,
           sigma = sigmaInit , # fixed, not estimated
           b0 = b0Init  # fixed, not estimated
)
modelData( bugsData( datalist ) )

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

nChain = 3
modelCompile( numChains = nChain )

genInitList <- function() {
    list(
        b = bInit , # from lm(y~zx), above
        thresh = 1:(nYlevels-1)+.5
    )
}
for ( chainIdx in 1 : nChain ) {
    modelInits( bugsInits( genInitList ) )
}

#------------------------------------------------------------------------------
# RUN THE CHAINS

# burn in
BurnInSteps = 2000
modelUpdate( BurnInSteps )
# actual samples
samplesSet( c( "b" , "thresh" ) )
stepsPerChain = ceiling(5000/nChain)
thinStep = 20
modelUpdate( stepsPerChain , thin=thinStep )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

source("plotChains.R")
source("plotPost.R")

checkConvergence = T
if ( checkConvergence ) {
  bSum      = plotChains( "b"      , saveplots=F , filenameroot=fname )
  threshSum = plotChains( "thresh" , saveplots=F , filenameroot=fname )
}

# Extract chain values:
zbSamp = NULL
for ( j in 1:nPredictors ) {
   zbSamp = cbind( zbSamp , samplesSample( paste("b[",j,"]",sep="") ) )
}
chainLength = NROW(zbSamp)
zthreshSamp = NULL
for ( j in 1:(nYlevels-1) ) {
   zthreshSamp = cbind( zthreshSamp ,
                        samplesSample(paste("thresh[",j,"]",sep="")) )
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

# Scatter plots of parameter values, pairwise:
if ( (nPredictors+nYlevels) <= 10 ) { # don't display if too many
    windows()
    thinIdx = ceiling(seq(1,chainLength,length=200))
    pairs( cbind( zbSamp[thinIdx,] , zthreshSamp[thinIdx,] )  ,
           labels=c( paste("zb",predictorNames,sep="") ,
                     paste("zthresh",1:nYlevels,sep="")) )
    windows()
    pairs( cbind( bSamp[thinIdx,] , threshSamp[thinIdx,] )  ,
           labels=c( paste("b",predictorNames,sep="") ,
                     paste("thresh",1:nYlevels,sep="")) )
    dev.copy2eps(file=paste(fname,"PostPairs.eps",sep=""))
}

# Display the posterior:
nPlotPerRow = 5
nPlotRow = ceiling((nPredictors+nYlevels-1)/nPlotPerRow)
nPlotCol = ceiling((nPredictors+nYlevels-1)/nPlotRow)
windows(3.5*nPlotCol,2.25*nPlotRow)
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
dev.copy2eps(file=paste(fname,"PostHist.eps",sep=""))

# Plot the data
if ( nPredictors == 2 ) {
windows()
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
dev.copy2eps(file=paste(fname,"Data.eps",sep=""))

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
