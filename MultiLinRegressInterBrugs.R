graphics.off()
rm(list=ls(all=TRUE))
fname = "MultiLinRegressInterBrugs"
library(BRugs)         # Kruschke, J. K. (2010). Doing Bayesian data analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
# BUGS model specification begins here...
model {
    for( i in 1 : nData ) {
        y[i] ~ dnorm( mu[i] , tau )
        mu[i] <- b0  +  b1 * x[i,1]  +  b2 * x[i,2]  +  b12 * x[i,1] * x[i,2]
    }
    tau ~ dgamma(.001,.001)       
    b0  ~ dnorm(0,1.0E-12)
    b1  ~ dnorm(0,1.0E-12)
    b2  ~ dnorm(0,1.0E-12)
    b12 ~ dnorm(0,1.0E-12)
}
# ... end BUGS model specification
" # close quote for modelstring
writeLines(modelstring,con="model.txt")
modelCheck( "model.txt" )

#------------------------------------------------------------------------------
# THE DATA.

dataSource = c("Guber1999","McIntyre1994","random")[1]

if ( dataSource=="Guber1999" ) {
   fname = paste(fname,"Guber1999",sep="") # file name for saved graphs
   dataMat = read.table( file="Guber1999data.txt" ,
                         col.names = c( "State","Spend","StuTchRat","Salary",
                                        "PrcntTake","SATV","SATM","SATT") )
   # Specify variables to be used in analysis:
   predictedName = "SATT"
   predictorNames = c( "Spend" , "PrcntTake" ) # only two predictors allowed
   nData = NROW( dataMat )
   y = as.matrix( dataMat[,predictedName] )
   x = as.matrix( dataMat[,predictorNames] )
   #nPredictors = NCOL( x )
}

if ( dataSource=="McIntyre1994" ) {
   fname = paste(fname,"McIntyre1994",sep="") # file name for saved graphs
   dataMat = read.csv(file="McIntyre1994data.csv")
   predictedName = "CO"
   predictorNames = c("Tar","Nic") # only two predictors allowed
   nData = NROW( dataMat )
   y = as.matrix( dataMat[,predictedName] )
   x = as.matrix( dataMat[,predictorNames] )
   #nPredictors = NCOL( x )
}

if ( dataSource=="random" ) {
   fname = paste(fname,"Random",sep="")  # file name for saved graphs
   # Generate random data.
   # True parameter values:
   betaTrue = c( 100 , 1 , 1 , -1 )       # b0,b1,b2,b12
   nPredictors = 2
   sdTrue = 1
   tauTrue = 1/sdTrue^2
   # Random X values:
   set.seed(47405)
   xM = 5 ; xSD = 4*sdTrue
   nData = 100
   x = matrix( rnorm( nPredictors*nData , xM , xSD ) , nrow=nData )
   predictorNames = colnames(x) = paste("X",1:nPredictors,sep="")
   # Random Y values generated from linear model with true parameter values:
   y = cbind( betaTrue[1]
              + betaTrue[2] * x[,1]
              + betaTrue[3] * x[,2]
              + betaTrue[4] * x[,1]*x[,2]
              + rnorm(nData,0,sdTrue) )
   predictedName = "Y"
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
           nData = nData
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
        b1 = bInit[1] ,
        b2 = bInit[2] ,
        b12 = 0 ,
        tau = tauInit
    )
}
for ( chainIdx in 1 : nChain ) {
    modelInits( bugsInits( genInitList ) )
}

#------------------------------------------------------------------------------
# RUN THE CHAINS

# burn in
BurnInSteps = 1000
modelUpdate( BurnInSteps )
# actual samples
samplesSet( c( "b0" , "b1" , "b2" , "b12" , "tau" ) )
stepsPerChain = ceiling(100000/nChain)
thinStep = 2
modelUpdate( stepsPerChain , thin=thinStep )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

source("plotChains.R")
source("plotPost.R")

checkConvergence = F
if ( checkConvergence ) {
  b0Sum  = plotChains( "b0"  , saveplots=F , filenameroot=fname )
  b1Sum  = plotChains( "b1"  , saveplots=F , filenameroot=fname )
  b2Sum  = plotChains( "b2"  , saveplots=F , filenameroot=fname )
  b12Sum = plotChains( "b12" , saveplots=F , filenameroot=fname )
  tauSum = plotChains( "tau" , saveplots=F , filenameroot=fname )
}

# Extract chain values:
zb0Samp = matrix( samplesSample( "b0" ) )
zb1Samp = matrix( samplesSample( "b1" ) )
zb2Samp = matrix( samplesSample( "b2" ) )
zb12Samp = matrix( samplesSample( "b12" ) )
zTauSamp = matrix( samplesSample( "tau" ) )
zSigmaSamp = 1 / sqrt( zTauSamp ) # Convert precision to SD

# Convert to original scale:
chainLength = length( zb0Samp )
My = mean(y)          # y is original-scale data
SDy = apply(y,2,sd)
Mx = apply(x,2,mean)  # x is original-scale data
SDx = apply(x,2,sd)
b0Samp = 0 * zb0Samp
b1Samp = 0 * zb1Samp
b2Samp = 0 * zb2Samp
b12Samp = 0 * zb12Samp
for ( stepIdx in 1:chainLength ) {
    b0Samp[stepIdx] = ( My + SDy * zb0Samp[stepIdx]
                       - SDy * (Mx[1]/SDx[1]) * zb1Samp[stepIdx]
                       - SDy * (Mx[2]/SDx[2]) * zb2Samp[stepIdx]
                       + SDy * (Mx[1]/SDx[1]) * (Mx[2]/SDx[2]) * zb12Samp[stepIdx] )
    b1Samp[stepIdx] = ( zb1Samp[stepIdx] * SDy / SDx[1]
                       - zb12Samp[stepIdx] * SDy * (1/SDx[1]) * (Mx[2]/SDx[2]) ) # corrected
    b2Samp[stepIdx] = ( zb2Samp[stepIdx] * SDy / SDx[2]
                       - zb12Samp[stepIdx] * SDy * (Mx[1]/SDx[1]) * (1/SDx[2]) ) # corrected
    b12Samp[stepIdx] = zb12Samp[stepIdx] * SDy * (1/SDx[1]) * (1/SDx[2])
}
sigmaSamp = zSigmaSamp * SDy

save( b0Samp , b1Samp , b2Samp , b12Samp , sigmaSamp ,
      file=paste(fname,".Rdata",sep="") )

# Scatter plots of parameter values, pairwise:
windows()
thinIdx = seq(1,length(b0Samp),length=200)
pairs( cbind( sigmaSamp[thinIdx] , b0Samp[thinIdx] , b1Samp[thinIdx,] ,
                  b2Samp[thinIdx,] , b12Samp[thinIdx,] ) ,
       labels=c( "Sigma y" , "Intercept" ,
                 paste("Beta",predictorNames[1],sep="") ,
                 paste("Beta",predictorNames[2],sep="") ,
                 "Interaction" ) , col="skyblue" )
#dev.copy2eps(file=paste(fname,"PostPairs.eps",sep=""))

# Display the posterior:
windows(3.5*5,2.75)
layout( matrix(1:5,nrow=1) )
par( mar=c(4,3,5,0) , mgp=c(2,0.7,0) )
histInfo = plotPost( sigmaSamp , xlab="Sigma Value" , compVal=NULL ,
                     breaks=30 , main=bquote(sigma[y]) ,
                     cex.main=1.67 , cex.lab=1.33 )
histInfo = plotPost( b0Samp , xlab="Intercept Value" , compVal=NULL ,
                     breaks=30 , main=bquote(.(predictedName) *" at all "* x==0) ,
                     cex.main=1.67 , cex.lab=1.33 )
histInfo = plotPost( b1Samp , xlab="Beta Value" , compVal=NULL ,
                     breaks=30 , cex.main=1.5 , cex.lab=1.33 ,
                     main=bquote( atop( Delta * .(predictedName) /
                                  Delta * .(predictorNames[1])  ,
                                  " at "* .(predictorNames[2])==0 ) ) )
histInfo = plotPost( b2Samp , xlab="Beta Value" , compVal=NULL ,
                     breaks=30 , cex.main=1.5 , cex.lab=1.33 ,
                     main=bquote( atop( Delta * .(predictedName) /
                                  Delta * .(predictorNames[2])  ,
                                  " at "* .(predictorNames[1])==0 ) ) )
histInfo = plotPost( b12Samp , xlab="Interaction Value" , compVal=NULL ,
                     breaks=30 , cex.main=1.67 , cex.lab=1.33 ,
                     main=paste( predictorNames[1],"x",predictorNames[2] ) )
#dev.copy2eps(file=paste(fname,"PostHist.eps",sep=""))

# Credible slopes as function of value of other predictor:
source("HDIofMCMC.R")
#
windows(7,5)
par( mar=c(4,4,3,0) , mgp=c(2,0.7,0) )
x2low = max( min(x[,2]) - 0.1 * ( max(x[,2]) - min(x[,2]) ) , 0 )
x2high = max(x[,2]) + 0.1 * ( max(x[,2]) - min(x[,2]) )
x2comb = seq( x2low , x2high , length=20 )
beta1HDI = matrix( 0 , nrow=3 , ncol=length(x2comb) )
for ( x2idx in 1:length(x2comb) ) {
    slope1Samp = b1Samp + b12Samp * x2comb[x2idx]
    HDIlim = HDIofMCMC( slope1Samp )
    beta1HDI[,x2idx] = c( HDIlim[1] , mean(slope1Samp) , HDIlim[2] )
}
plot( x2comb , beta1HDI[2,] , type="o" , pch="+" , cex=2 , col="skyblue" ,
      ylim=c(min(beta1HDI),max(beta1HDI)) ,
      xlab=bquote("Value of "*.(predictorNames[2])) ,
      ylab=bquote("Slope along "*.(predictorNames[1])) ,
      main="Posterior mean and 95% HDI of slope" ,
      cex.lab=1.5 )
abline( h=0 , lty="dashed" )
segments( x2comb , beta1HDI[1,] , x2comb , beta1HDI[3,] , lwd=4 , col="skyblue" )
#dev.copy2eps(file=paste(fname,"PostSlope1.eps",sep=""))
#
windows(7,5)
par( mar=c(4,4,3,0) , mgp=c(2,0.7,0) )
# x1low = max( min(x[,1]) - 0.1 * ( max(x[,1]) - min(x[,2]) ) , 0 )
# x1high = max(x[,1]) + 0.1 * ( max(x[,1]) - min(x[,1]) )
x1low = 0 ; x1high = 50
x1comb = seq( x1low , x1high , length=20 )
beta2HDI = matrix( 0 , nrow=3 , ncol=length(x1comb) )
for ( x1idx in 1:length(x1comb) ) {
    slope2Samp = b2Samp + b12Samp * x1comb[x1idx]
    HDIlim = HDIofMCMC( slope2Samp )
    beta2HDI[,x1idx] = c( HDIlim[1] , mean(slope2Samp) , HDIlim[2] )
}
plot( x1comb , beta2HDI[2,] , type="o" , pch="+" , cex=2 , col="skyblue" ,
      ylim=c(min(beta2HDI),max(beta2HDI)) ,
      xlab=bquote("Value of "*.(predictorNames[1])) ,
      ylab=bquote("Slope along "*.(predictorNames[2])) ,
      main="Posterior mean and 95% HDI of slope" ,
      cex.lab=1.5 )
abline( h=0 , lty="dashed" )
segments( x1comb , beta2HDI[1,] , x1comb , beta2HDI[3,] , lwd=4 , col="skyblue" )
#dev.copy2eps(file=paste(fname,"PostSlope2.eps",sep=""))

#------------------------------------------------------------------------------