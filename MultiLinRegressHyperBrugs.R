graphics.off()
rm(list=ls(all=TRUE))
fname = "MultiLinRegressHyper"
library(BRugs)         # Kruschke, J. K. (2010). Doing Bayesian data analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
# BUGS model specification begins here...
model {
    for( i in 1 : nData ) {
        y[i] ~ dnorm( mu[i] , tau )
        mu[i] <- b0 + inprod( b[] , x[i,] )
    }
    tau ~ dgamma(.01,.01)
    b0 ~ dnorm(0,1.0E-12) 
    for ( j in 1:nPredictors ) {
        b[j] ~ dt( muB , tauB , tdfB )
    }
    muB ~ dnorm( 0 , .100 )
    tauB ~ dgamma(.01,.01)
    udfB ~ dunif(0,1)
    tdfB <- 1 + tdfBgain * ( -log( 1 - udfB ) )
}
# ... end BUGS model specification
" # close quote for modelstring
writeLines(modelstring,con="model.txt")
modelCheck( "model.txt" )

#------------------------------------------------------------------------------
# THE DATA.

tdfBgain = 1

dataSource = c("Guber1999","McIntyre1994","random")[3]

if ( dataSource=="Guber1999" ) {
   fname = paste("Guber1999","tdf",tdfBgain,sep="")
   dataMat = read.table( file="Guber1999data.txt" ,
                         col.names = c( "State","Spend","StuTchRat","Salary",
                                        "PrcntTake","SATV","SATM","SATT") )
   # Specify variables to be used in BUGS analysis:
   predictedName = "SATT"
   predictorNames = c( "Spend" , "PrcntTake" )
   #predictorNames = c( "Spend" , "PrcntTake" , "Salary" , "StuTchRat" )
   nData = NROW( dataMat )
   y = as.matrix( dataMat[,predictedName] )
   x = as.matrix( dataMat[,predictorNames] )
   nPredictors = NCOL( x )
}

if ( dataSource=="McIntyre1994" ) {
   fname = paste("McIntyre1994","tdf",tdfBgain,sep="")
   dataMat = read.csv(file="McIntyre1994data.csv")
   predictedName = "CO"
   predictorNames = c("Tar","Nic","Wt")
   nData = NROW( dataMat )
   y = as.matrix( dataMat[,predictedName] )
   x = as.matrix( dataMat[,predictorNames] )
   nPredictors = NCOL( x )
}

if ( dataSource=="random" ) {
   fname = paste("Random","tdf",tdfBgain,sep="")
   # Generate random data.
   # True parameter values:
   betaTrue = c( 100 , 1 , 2 , rep(0,21) ) # beta0 is first component
   nPredictors = length( betaTrue ) - 1
   sdTrue = 2
   tauTrue = 1/sdTrue^2
   # Random X values:
   set.seed(47405)
   xM = 5 ; xSD = 2
   nData = 100
   x = matrix( rnorm( nPredictors*nData , xM , xSD ) , nrow=nData )
   predictorNames = colnames(x) = paste("X",1:nPredictors,sep="")
   # Random Y values generated from linear model with true parameter values:
   y = x %*% matrix(betaTrue[-1],ncol=1) + betaTrue[1] + rnorm(nData,0,sdTrue)
   predictedName = "Y"
   # Select which predictors to include
   includeOnly = 1:nPredictors # default is to include all
   #includeOnly = 1:6 # subset of predictors overwrites default
   x = x[,includeOnly]
   predictorNames = predictorNames[includeOnly]
   nPredictors = NCOL(x)
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
zy = standardizeCols( y )

# Get the data into BUGS:
datalist = list(
           x = zx ,
           y = as.vector( zy ) , # BUGS does not treat 1-column mat as vector
           nPredictors = nPredictors ,
           nData = nData ,
           tdfBgain = tdfBgain
)
modelData( bugsData( datalist ) )

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

nChain = 3
modelCompile( numChains = nChain )

genInitList <- function(nPred=nPredictors) {
    lmInfo = lm( y ~ x ) # R function returns least-squares (normal MLE) fit.
    bInit = lmInfo$coef[-1] * apply(x,2,sd) / apply(y,2,sd)
    tauInit = (length(y)*apply(y,2,sd)^2)/sum(lmInfo$res^2)
    list(
        b0 = 0 ,
        b = bInit ,
        tau = tauInit ,
        muB =  mean( bInit ) ,
        tauB = 1 / sd( bInit )^2 ,
        udfB = 0.95 # tdfB = 4
    )
}
for ( chainIdx in 1 : nChain ) {
    modelInits( bugsInits( genInitList ) )
}

#------------------------------------------------------------------------------
# RUN THE CHAINS

# burn in
BurnInSteps = 100
modelUpdate( BurnInSteps )
# actual samples
samplesSet( c( "b0" , "b" , "tau" , "muB" , "tauB" , "tdfB" ) )
stepsPerChain = ceiling(10000/nChain)
thinStep = 2
modelUpdate( stepsPerChain , thin=thinStep )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

source("plotChains.R")
source("plotPost.R")

checkConvergence = F
if ( checkConvergence ) {
  b0Sum  = plotChains( "b0"  , saveplots=F , filenameroot=fname )
  bSum   = plotChains( "b"   , saveplots=F , filenameroot=fname )
  tauSum = plotChains( "tau" , saveplots=F , filenameroot=fname )
  muBSum = plotChains( "muB" , saveplots=F , filenameroot=fname )
  tauBSum = plotChains( "tauB" , saveplots=F , filenameroot=fname )
  tdfBSum = plotChains( "tdfB" , saveplots=F , filenameroot=fname )
}

# Extract chain values:
zb0Samp = matrix( samplesSample( "b0" ) )
zbSamp = NULL
for ( j in 1:nPredictors ) {
   zbSamp = cbind( zbSamp , samplesSample( paste("b[",j,"]",sep="") ) )
}
zTauSamp = matrix( samplesSample( "tau" ) )
zSigmaSamp = 1 / sqrt( zTauSamp ) # Convert precision to SD
chainLength = length(zTauSamp)

# Convert to original scale:
bSamp = zbSamp * matrix( apply(y,2,sd)/apply(x,2,sd) , byrow=TRUE ,
                     ncol=nPredictors , nrow=NROW(zbSamp) )
b0Samp = ( zb0Samp * apply(y,2,sd)
          + mean(y)
          - rowSums( zbSamp
          * matrix( apply(y,2,sd)/apply(x,2,sd) , byrow=TRUE ,
                    ncol=nPredictors , nrow=NROW(zbSamp) )
          * matrix( apply(x,2,mean) , byrow=TRUE ,
                    ncol=nPredictors , nrow=NROW(zbSamp) ) ) )
sigmaSamp = zSigmaSamp * apply(y,2,sd)

# Scatter plots of parameter values, pairwise:
if ( nPredictors <= 6 ) { # don't display if too many predictors
    windows()
    thinIdx = round(seq(1,length(zb0Samp),length=200))
    pairs( cbind( zSigmaSamp[thinIdx] , zb0Samp[thinIdx] , zbSamp[thinIdx,] )  ,
      labels=c("Sigma zy","zIntercept",paste("zSlope",predictorNames,sep="")))
    windows()
    thinIdx = seq(1,length(b0Samp),length=700)
    pairs( cbind( sigmaSamp[thinIdx] , b0Samp[thinIdx] , bSamp[thinIdx,] ) ,
      labels=c( "Sigma y" , "Intercept", paste("Slope",predictorNames,sep="")))
    dev.copy2eps(file=paste(fname,"PostPairs.eps",sep=""))
}
# Show correlation matrix on console:
cat("\nCorrlations of posterior sigma, b0, and bs:\n")
show( cor( cbind( sigmaSamp , b0Samp , bSamp ) ) )

# Display the posterior:
nPlotPerRow = 5
nPlotRow = ceiling((2+nPredictors)/nPlotPerRow)
nPlotCol = ceiling((2+nPredictors)/nPlotRow)
windows(3.5*nPlotCol,2.25*nPlotRow)
layout( matrix(1:(nPlotRow*nPlotCol),nrow=nPlotRow,ncol=nPlotCol,byrow=T) )
par( mar=c(4,3,2.5,0) , mgp=c(2,0.7,0) )
histInfo = plotPost( sigmaSamp , xlab="Sigma Value" , compVal=NULL ,
                     breaks=30 , main=bquote(sigma[y]) ,
                     cex.main=1.67 , cex.lab=1.33 )
histInfo = plotPost( b0Samp , xlab="Intercept Value" , compVal=NULL ,
                     breaks=30 , main=bquote(.(predictedName) *" at "* x==0) ,
                     cex.main=1.67 , cex.lab=1.33 )
for ( sIdx in 1:nPredictors ) {
histInfo = plotPost( bSamp[,sIdx] , xlab="Slope Value" , compVal=0.0 ,
                     breaks=30 ,
                     main=bquote( Delta * .(predictedName) /
                                  Delta * .(predictorNames[sIdx]) ) ,
                     cex.main=1.67 , cex.lab=1.33 )
}
dev.copy2eps(file=paste(fname,"PostHist.eps",sep=""))

# Posterior prediction:
# Specify x values for which predicted y's are needed.
# xPostPred is a matrix such that ncol=nPredictors and nrow=nPostPredPts.
xPostPred = rbind(
    apply(x,2,mean)-3*apply(x,2,sd) , # mean of data x minus thrice SD of data x
    apply(x,2,mean)                 , # mean of data x
    apply(x,2,mean)+3*apply(x,2,sd)   # mean of data x plus thrice SD of data x
)
# Define matrix for recording posterior predicted y values for each xPostPred.
# One row per xPostPred value, with each row holding random predicted y values.
postSampSize = chainLength
yPostPred = matrix( 0 , nrow=NROW(xPostPred) , ncol=postSampSize )
# Define matrix for recording HDI limits of posterior predicted y values:
yHDIlim = matrix( 0 , nrow=NROW(xPostPred) , ncol=2 )
# Generate posterior predicted y values.
# This gets only one y value, at each x, for each step in the chain.
for ( chainIdx in 1:chainLength ) {
    yPostPred[,chainIdx] = rnorm( NROW(xPostPred) ,
                           mean = b0Samp[chainIdx]
                                  + xPostPred %*% cbind(bSamp[chainIdx,]) ,
                           sd = rep( sigmaSamp[chainIdx] , NROW(xPostPred) ) )
}
source("HDIofMCMC.R")
for ( xIdx in 1:NROW(xPostPred) ) {
    yHDIlim[xIdx,] = HDIofMCMC( yPostPred[xIdx,] )
}
cat( "\nPosterior predicted y for selected x:\n" )
show( cbind( xPostPred , yPostPredMean=rowMeans(yPostPred) , yHDIlim ) )

#------------------------------------------------------------------------------