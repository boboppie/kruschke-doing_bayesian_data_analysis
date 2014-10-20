graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="MultiLinRegressInterJags" # for constructing output filenames
if ( .Platform$OS.type != "windows" ) { 
  windows <- function( ... ) X11( ... ) 
}
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
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
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.

dataSource = c("Guber1999","McIntyre1994","random")[1]

if ( dataSource=="Guber1999" ) {
   fileNameRoot = paste(fileNameRoot,"Guber1999",sep="") # file name for saved graphs
   dataMat = read.table( file="Guber1999data.txt" ,
                         col.names = c( "State","Spend","StuTchRat","Salary",
                                        "PrcntTake","SATV","SATM","SATT") )
   # Specify variables to be used in analysis:
   predictedName = "SATT"
   predictorNames = c( "Spend" , "PrcntTake" )
   #predictorNames = c( "Spend" , "PrcntTake" , "Salary" , "StuTchRat" )
   nData = NROW( dataMat )
   y = as.matrix( dataMat[,predictedName] )
   x = as.matrix( dataMat[,predictorNames] )
   nPredictors = NCOL( x )
}

if ( dataSource=="McIntyre1994" ) {
   fileNameRoot = paste(fileNameRoot,"McIntyre1994",sep="") # file name for saved graphs
   dataMat = read.csv(file="McIntyre1994data.csv")
   predictedName = "CO"
   predictorNames = c("Tar","Nic","Wt")
   nData = NROW( dataMat )
   y = as.matrix( dataMat[,predictedName] )
   x = as.matrix( dataMat[,predictorNames] )
   nPredictors = NCOL( x )
}

if ( dataSource=="random" ) {
   fileNameRoot = paste(fileNameRoot,"Random",sep="")  # file name for saved graphs
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
   #includeOnly = 1:10 # subset of predictors overwrites default
   x = x[,includeOnly]
   predictorNames = predictorNames[includeOnly]
   nPredictors = NCOL(x)
}

# Prepare data for JAGS:
# Re-center data at mean, to reduce autocorrelation in MCMC sampling.
# Standardize (divide by SD) to make prior specification easier.
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
           nData = nData
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

lmInfo = lm( y ~ x ) # R function returns least-squares (normal MLE) fit.
bInit = lmInfo$coef[-1] * apply(x,2,sd) / apply(y,2,sd)
tauInit = (length(y)*apply(y,2,sd)^2)/sum(lmInfo$res^2)
initsList = list(
    b0 = 0 ,
    b1 = bInit[1] ,
    b2 = bInit[2] ,
    b12 = 0 ,
    tau = tauInit
)

#------------------------------------------------------------------------------
# RUN THE CHAINS

parameters = c( "b0" , "b1" , "b2" , "b12" , "tau" ) 
adaptSteps = 5000             # Number of steps to "tune" the samplers.
burnInSteps = 10000           # Number of steps to "burn-in" the samplers.
nChains = 3                   # Number of chains to run.
numSavedSteps=100000          # Total number of steps in chains to save.
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
  windows()
  plot( codaSamples , ask=F )  
  windows()
  autocorr.plot( codaSamples , ask=F )
}

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )

# Extract chain values:
zb0Samp = matrix( mcmcChain[, "b0" ] )
zb1Samp = matrix( mcmcChain[, "b1" ] )
zb2Samp = matrix( mcmcChain[, "b2" ] )
zb12Samp = matrix( mcmcChain[, "b12" ] )
zTauSamp = matrix( mcmcChain[, "tau" ] )
zSigmaSamp = 1 / sqrt( zTauSamp ) # Convert precision to SD
chainLength = length(zTauSamp)

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
      file=paste(fileNameRoot,".Rdata",sep="") )

source("plotPost.R")

# Scatter plots of parameter values, pairwise:
windows()
thinIdx = seq(1,length(b0Samp),length=200)
pairs( cbind( sigmaSamp[thinIdx] , b0Samp[thinIdx] , b1Samp[thinIdx,] ,
                  b2Samp[thinIdx,] , b12Samp[thinIdx,] ) ,
       labels=c( "Sigma y" , "Intercept" ,
                 paste("Beta",predictorNames[1],sep="") ,
                 paste("Beta",predictorNames[2],sep="") ,
                 "Interaction" ) , col="skyblue" )
savePlot(file=paste(fileNameRoot,"PostPairs",sep=""),type="eps")
savePlot(file=paste(fileNameRoot,"PostPairs",sep=""),type="jpg")

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
savePlot(file=paste(fileNameRoot,"PostHist",sep=""),type="eps")
savePlot(file=paste(fileNameRoot,"PostHist",sep=""),type="jpg")

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
savePlot(file=paste(fileNameRoot,"PostSlope1",sep=""),type="eps")
savePlot(file=paste(fileNameRoot,"PostSlope1",sep=""),type="jpg")
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
savePlot(file=paste(fileNameRoot,"PostSlope2",sep=""),type="eps")
savePlot(file=paste(fileNameRoot,"PostSlope2",sep=""),type="jpg")

#------------------------------------------------------------------------------