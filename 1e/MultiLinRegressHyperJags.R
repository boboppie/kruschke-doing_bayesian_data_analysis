graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="MultiLinRegressHyperJags" # for constructing output filenames
source("openGraphSaveGraph.R")
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
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
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.

tdfBgain = 1

dataSource = c("Guber1999","McIntyre1994","random")[3]

if ( dataSource=="Guber1999" ) {
   fileNameRoot = paste("Guber1999","tdf",tdfBgain,sep="")
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
   fileNameRoot = paste("McIntyre1994","tdf",tdfBgain,sep="")
   dataMat = read.csv(file="McIntyre1994data.csv")
   predictedName = "CO"
   predictorNames = c("Tar","Nic","Wt")
   nData = NROW( dataMat )
   y = as.matrix( dataMat[,predictedName] )
   x = as.matrix( dataMat[,predictorNames] )
   nPredictors = NCOL( x )
}

if ( dataSource=="random" ) {
   fileNameRoot = paste("Random","tdf",tdfBgain,sep="")
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

dataList = list(
           x = zx ,
           y = as.vector( zy ) , # BUGS does not treat 1-column mat as vector
           nPredictors = nPredictors ,
           nData = nData ,
           tdfBgain = tdfBgain
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

lmInfo = lm( y ~ x ) # R function returns least-squares (normal MLE) fit.
bInit = lmInfo$coef[-1] * apply(x,2,sd) / apply(y,2,sd)
tauInit = (length(y)*apply(y,2,sd)^2)/sum(lmInfo$res^2)
initsList = list(
    b0 = 0 ,
    b = bInit ,
    tau = tauInit ,
    muB =  mean( bInit ) ,
    tauB = 1 / sd( bInit )^2 ,
    udfB = 0.95 # tdfB = 4
)

#------------------------------------------------------------------------------
# RUN THE CHAINS

parameters = c( "b0" , "b" , "tau" , "muB" , "tauB" , "tdfB" )  
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

# Extract chain values:
zb0Samp = matrix( mcmcChain[, "b0" ] )
zbSamp = NULL
for ( j in 1:nPredictors ) {
   zbSamp = cbind( zbSamp , mcmcChain[, paste("b[",j,"]",sep="") ] )
}
zTauSamp = matrix( mcmcChain[, "tau" ] )
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

source("plotPost.R")

# Scatter plots of parameter values, pairwise:
if ( nPredictors <= 6 ) { # don't display if too many predictors
    openGraph()
    thinIdx = round(seq(1,length(zb0Samp),length=200))
    pairs( cbind( zSigmaSamp[thinIdx] , zb0Samp[thinIdx] , zbSamp[thinIdx,] )  ,
      labels=c("Sigma zy","zIntercept",paste("zSlope",predictorNames,sep="")))
    openGraph()
    thinIdx = seq(1,length(b0Samp),length=700)
    pairs( cbind( sigmaSamp[thinIdx] , b0Samp[thinIdx] , bSamp[thinIdx,] ) ,
      labels=c( "Sigma y" , "Intercept", paste("Slope",predictorNames,sep="")))
    saveGraph(file=paste(fileNameRoot,"PostPairs.eps",sep=""),type="eps")
}
# Show correlation matrix on console:
cat("\nCorrlations of posterior sigma, b0, and bs:\n")
show( cor( cbind( sigmaSamp , b0Samp , bSamp ) ) )

# Display the posterior:
nPlotPerRow = 5
nPlotRow = ceiling((2+nPredictors)/nPlotPerRow)
nPlotCol = ceiling((2+nPredictors)/nPlotRow)
openGraph(3.5*nPlotCol,2.25*nPlotRow)
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
saveGraph(file=paste(fileNameRoot,"PostHist.eps",sep=""),type="eps")

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