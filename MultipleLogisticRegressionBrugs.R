graphics.off()
rm(list=ls(all=TRUE))
fname = "MultipleLogisticRegressionBrugs"
library(BRugs)         # Kruschke, J. K. (2010). Doing Bayesian data analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.
modelstring = "
# BUGS model specification begins here...
model {
  for( i in 1 : nData ) {
    y[i] ~ dbern( mu[i] )
    mu[i] <- 1/(1+exp(-( b0 + inprod( b[] , x[i,] ))))
  }
  b0 ~ dnorm( 0 , 1.0E-12 )
  for ( j in 1 : nPredictors ) {
    b[j] ~ dnorm( 0 , 1.0E-12 )
  }
}
# ... end BUGS model specification
" # close quote for modelstring
# Write model to a file:
writeLines(modelstring,con="model.txt")
modelCheck( "model.txt" )

#------------------------------------------------------------------------------
# THE DATA.

dataSource = c( "HtWt" , "Cars" , "HeartAttack" )[3]

if ( dataSource == "HtWt" ) {
  fname = paste( fname , dataSource , sep="" )
  # Generate random but realistic data:
  source( "HtWtDataGenerator.R" )
  dataMat = HtWtDataGenerator( nSubj = 70 , rndsd=474 )
  predictedName = "male"
  predictorNames = c( "height" , "weight" )
  nData = NROW( dataMat )
  y = as.matrix( dataMat[,predictedName] )
  x = as.matrix( dataMat[,predictorNames] )
  nPredictors = NCOL( x )
}

if ( dataSource == "Cars" ) {
  fname = paste( fname , dataSource , sep="" )
  dataMat = read.table(file="Lock1993data.txt",header=T,sep=" ")
  predictedName = "AirBag"
  predictorNames = c( "MidPrice" , "RPM" , "Uturn" )
  nData = NROW( dataMat )
  y = as.matrix( as.numeric( dataMat[,predictedName] > 0 ) ) # 0,1,2 to 0,1
  x = as.matrix( dataMat[,predictorNames] )
  nPredictors = NCOL( x )
}

if ( dataSource == "HeartAttack" ) {
  fname = paste( fname , dataSource , sep="" )
  dataMat = read.table(file="BloodDataGeneratorOutput.txt",header=T,sep=" ")
  predictedName = "HeartAttack"
  predictorNames = c( "Systolic", "Diastolic", "Weight", "Cholesterol",
                      "Height", "Age" )
#  predictorNames = c( "Systolic", "Diastolic" )
  nData = NROW( dataMat )
  y = as.matrix( dataMat[,predictedName] )
  x = as.matrix( dataMat[,predictorNames] )
  nPredictors = NCOL( x )
}

# Prepare data for BUGS:
# Re-center data at mean, to reduce autocorrelation in MCMC sampling.
# Standardize (divide by SD) to make initialization easier.
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
zy = y  # y is not standardized; must be 0,1

# Get the data into BUGS:
datalist = list(
           x = zx ,
           y = as.vector( zy ) , # BUGS does not treat 1-column mat as vector
           nPredictors = nPredictors ,
           nData = nData
)
modelData( bugsData( datalist ) )


#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

nchain = 3
modelCompile( numChains = nchain )

genInitList <- function() { 
    glmInfo = glm( datalist$y ~ datalist$x , family=binomial(logit) ) # R func.
    show( glmInfo ) ; flush.console() # display in case glm() has troubles
    b0Init = glmInfo$coef[1] 
    bInit = glmInfo$coef[-1] 
    return( list(
        b0 = b0Init ,
        b = bInit
    ) )
}
for ( chainIdx in 1 : nchain ) {
    modelInits( bugsInits( genInitList ) )
}

#------------------------------------------------------------------------------
# RUN THE CHAINS

# burn in
BurnInSteps = 1000
modelUpdate( BurnInSteps )
# actual samples
samplesSet( c( "b0" , "b" ) )
stepsPerChain = ceiling(5000/nchain)
thinStep = 50  # check autocorrelation and increase as needed
modelUpdate( stepsPerChain , thin=thinStep )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

source("plotChains.R")
source("plotPost.R")

# Check chains for mixing
checkConvergence = T
if ( checkConvergence ) {
  b0Sum = plotChains( "b0" , saveplots=F , filenameroot=fname )
  bSum  = plotChains( "b"  , saveplots=F , filenameroot=fname )
}

# Extract chain values:
zb0Sample = matrix( samplesSample( "b0" ) )
chainLength = length(zb0Sample)
zbSample = NULL
for ( j in 1:nPredictors ) {
   zbSample = cbind( zbSample , samplesSample( paste("b[",j,"]",sep="") ) )
}

# Convert to original scale:
x = dataMat[,predictorNames]
y = dataMat[,predictedName]
My = mean(y)
SDy = sd(y)
Mx = apply(x,2,mean)
SDx = apply(x,2,sd)
b0Sample = 0 * zb0Sample
bSample = 0 * zbSample
for ( stepIdx in 1:chainLength ) {
    b0Sample[stepIdx] = ( zb0Sample[stepIdx]
                          - sum( Mx / SDx * zbSample[stepIdx,] ) )
    for ( j in 1:nPredictors ) {                      
      bSample[stepIdx,j] = zbSample[stepIdx,j] / SDx[j] 
    }
}

# Examine sampled values, z scale:
windows()
thinIdx = ceiling(seq(1,chainLength,length=700))
pairs(  cbind( zb0Sample[thinIdx] , zbSample[thinIdx,] )  ,
       labels=c( "zb0", paste("zb",predictorNames,sep="") ) )
# Examine sampled values, original scale:
windows()
pairs( cbind( b0Sample[thinIdx] , bSample[thinIdx,] ) ,
       labels=c( "b0", paste("b_",predictorNames,sep="") ) )
dev.copy2eps(file=paste(fname,"PostPairs.eps",sep=""))

# Display the posterior :
windows(3.5*(1+nPredictors),2.75)
layout( matrix(1:(1+nPredictors),nrow=1) )
histInfo = plotPost( b0Sample , xlab="b0 Value" , compVal=NULL , breaks=30 , 
                     main=paste( "logit(p(", predictedName ,
                                 "=1)) when predictors = zero" , sep="" ) )
for ( bIdx in 1:nPredictors ) {
histInfo = plotPost( bSample[,bIdx] , xlab=paste("b",bIdx," Value",sep="") , 
                     compVal=0.0 , breaks=30 ,
                     main=paste(predictorNames[bIdx]) )
}
dev.copy2eps(file=paste(fname,"PostHist.eps",sep=""))

# Plot data with .5 level contours of believable logistic surfaces.
# The contour lines are best interpreted when there are only two predictors.
for ( p1idx in 1:(nPredictors-1) ) {
  for ( p2idx in (p1idx+1):nPredictors ) {
    windows()
    xRange = range(x[,p1idx])
    yRange = range(x[,p2idx])
    # make empty plot
    plot( NULL , NULL , main=predictedName , xlim=xRange , ylim=yRange , 
          xlab=predictorNames[p1idx] , ylab=predictorNames[p2idx] )
    # Some of the 50% level contours from the posterior sample.
    for ( chainIdx in ceiling(seq( 1 , chainLength , length=20 )) ) {
      abline( -( b0Sample[chainIdx]
                 + if (nPredictors>2) {
                     bSample[chainIdx,c(-p1idx,-p2idx)]*Mx[c(-p1idx,-p2idx)]
                   } else { 0 } )
                 / bSample[chainIdx,p2idx] ,
              -bSample[chainIdx,p1idx]/bSample[chainIdx,p2idx] ,
              col="grey" , lwd = 2 )
    }
    # The data points:
    for ( yVal in 0:1 ) {
      rowIdx = ( y == yVal )
      points( x[rowIdx,p1idx] , x[rowIdx,p2idx] , pch=as.character(yVal) ,
              cex=1.75 )
    }
    dev.copy2eps(file=paste(fname,"PostContours",p1idx,p2idx,".eps",sep=""))
  }
}

#------------------------------------------------------------------------------

# MLE logistic regression:
glmRes = glm( datalist$y ~ as.matrix(x) , family=binomial(logit) )
show( glmRes )
