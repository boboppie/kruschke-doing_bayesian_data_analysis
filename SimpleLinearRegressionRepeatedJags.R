graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="SimpleLinearRegressionRepeatedJags" 
source("openGraphSaveGraph.R")
source("plotPost.R")
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.
modelstring = "
model {
    for( r in 1 : Ndata ) {
        y[r] ~ dnorm( mu[r] , tau[ subj[r] ] )
        mu[r] <- b0[ subj[r] ] + b1[ subj[r] ] * x[r]
    }
    for ( s in 1 : Nsubj ) {
        b0[s] ~ dnorm( mu0G , tau0G )
        b1[s] ~ dnorm( mu1G , tau1G )
        tau[s] ~ dgamma( sG , rG )
    }
    mu0G ~ dnorm(0,.01)
    tau0G ~ dgamma(.1,.1)
    mu1G ~ dnorm(0,.01)
    tau1G ~ dgamma(.1,.1)
    sG <- pow(m,2)/pow(d,2)
    rG <- m/pow(d,2)
    m ~ dgamma(1,.1)
    d ~ dgamma(1,.1)
}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.

# Data from H. A. Feldman, 1988, Table 4, p. 1731.
# Columns are "group" , "subjID" , "time" , "retention"
source("Feldman1988Table4data.R")
# Remove missing data:
includeRowVec = is.finite( Feldman1988Table4data[,"retention"] )
dataMat = Feldman1988Table4data[ includeRowVec , ]
# Retain only the Group 1 (lung) data:
dataMat = dataMat[ dataMat[,"group"]==1 , ]
# Convert to log10(retention):
dataMat[,"retention"] = log10( dataMat[,"retention"] )
# Column names and plot labels
yColName = "retention" ; yPlotLab = "log10 Retention"
xColName = "time" ; xPlotLab = "Day"
subjColName = "subjID" ; subjPlotLab = "Subject"
fileNameRoot = "SimpleLinearRegressionRepeatedBrugs"

if ( F ) { # change to T to use income data instead of contam.retention data.
 # Data from http://www.census.gov/hhes/www/income/statemedfaminc.html
 # Downloaded Dec. 06, 2009.
 load("IncomeFamszState.Rdata") # loads IncomeFamszState
 dataMat = IncomeFamszState
 yColName="Income" ; yPlotLab = "Income"
 xColName="Famsz" ; xPlotLab="Family Size"
 subjColName="State" ; subjPlotLab="State"
 fileNameRoot = "IncomeFamszState"
}

# Extract data info:
Ndata = NROW(dataMat)
# To make sure that subj has same order of subjects as dataMat, must use
# levels=unique() argument in factor statement...
subj = as.integer( factor( dataMat[,subjColName] ,
							levels=unique(dataMat[,subjColName] ) ) )
Nsubj = length(unique(subj))
x = as.numeric(dataMat[,xColName])
y = as.numeric(dataMat[,yColName])

# Re-center data at mean, to reduce autocorrelation in MCMC sampling.
# Standardize (divide by SD) to make initialization easier.
xM = mean( x ) ; xSD = sd( x )
yM = mean( y ) ; ySD = sd( y )
zx = ( x - xM ) / xSD
zy = ( y - yM ) / ySD

# Specify data, as a list.
dataList = list(
  Ndata = Ndata ,
  Nsubj = Nsubj ,
  subj = subj ,
  x = zx ,
  y = zy
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

b0 = b1 = tau = rep(0,length=Nsubj)
for ( sIdx in 1:Nsubj ) {
   yVec = dataList$y[dataList$subj==sIdx]
   xVec = dataList$x[dataList$subj==sIdx]
   lmInfo = lm( yVec ~ xVec )
   b0[sIdx] = lmInfo$coef[1]
   b1[sIdx] = lmInfo$coef[2]
   tau[sIdx] = length(yVec) / sum(lmInfo$res^2)
}
mu0G = mean(b0)
tau0G = 1/sd(b0)^2
mu1G = mean(b1)
tau1G = 1/sd(b1)^2
m = mean(tau)
d = sd(tau)
initsList = list( b0=b0 , b1=b1 , tau=tau ,
      mu0G=mu0G , tau0G=tau0G ,
      mu1G=mu1G , tau1G=tau1G ,
      m=m , d=d )

#------------------------------------------------------------------------------
# RUN THE CHAINS

parameters = c("b0","b1","tau" , "mu0G","tau0G", "mu1G","tau1G", "m","d")
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

# Extract chain values for subsequent examination:
zmu0Gsamp = mcmcChain[, "mu0G" ]
zmu1Gsamp = mcmcChain[, "mu1G" ]
zb0samp = NULL
zb1samp = NULL
for ( subjIdx in 1:Nsubj ) {
    zb0samp = rbind( zb0samp , mcmcChain[, paste("b0[",subjIdx,"]",sep="") ])
    zb1samp = rbind( zb1samp , mcmcChain[, paste("b1[",subjIdx,"]",sep="") ])
}

# Convert to original scale:
mu0Gsamp = zmu0Gsamp * ySD + yM - zmu1Gsamp * ySD * xM / xSD
mu1Gsamp = zmu1Gsamp * ySD / xSD
b0samp   = zb0samp   * ySD + yM - zb1samp   * ySD * xM / xSD
b1samp   = zb1samp   * ySD / xSD

# Display believable intercept and slope values
openGraph(10,5.5)
par( mar=c(4,4,1.75,1)+0.1 , mgp=c(2.5,0.8,0) )
layout( matrix(1:2,nrow=1) )
thinIdx = round(seq(1,length(mu0Gsamp),length=700))
plot( zmu1Gsamp[thinIdx] , zmu0Gsamp[thinIdx] , cex.lab=1.75 ,
      ylab="Standardized Intercept" , xlab="Standardized Slope" , 
      col="skyblue")
plot( mu1Gsamp[thinIdx] , mu0Gsamp[thinIdx] , cex.lab=1.0 ,
      ylab=paste("Intercept (",yPlotLab," when ",xPlotLab," =0)",sep="") ,
      xlab=paste("Slope (change in",yPlotLab,"per unit",xPlotLab,")") ,
      col="skyblue")
saveGraph(file=paste(fileNameRoot,"SlopeIntercept",sep=""),type="eps")

# Make graphs of data and corresponding believable slopes:
openGraph(12,6)
par( mar=c(4,4,1.75,1)+0.1 , mgp=c(2.5,0.8,0) )
layout(matrix(c(1:5,1:5,6:10),nrow=3,byrow=T))
xlims = c( min( dataMat[,xColName] ) ,  max( dataMat[,xColName] ) )
ylims = c( min( dataMat[,yColName] ) ,  max( dataMat[,yColName] ) )
sIdVec = unique( dataMat[,subjColName] )
# Plot data of individual subjects:
nSubjPlots = 4 # number of representative subject plots to make
subjIdxVec = round(seq(1,length(sIdVec),length=nSubjPlots))
for ( sIdx in subjIdxVec ) {
    rVec = ( dataMat[,subjColName] == sIdVec[sIdx] )
    plot( dataMat[rVec,xColName] , dataMat[rVec,yColName] , type="o" ,
          ylim=ylims , ylab=yPlotLab , xlab=xPlotLab , cex.lab=1.5 ,
          pch=sIdx%%26 , lty=sIdx , main=bquote(.(subjPlotLab) *" "* .(sIdx)) ,
          cex.main=1.75 )
}
# Plot data of all subjects superimposed
plot( NULL,NULL, xlab=xPlotLab,xlim=xlims , ylab=yPlotLab,ylim=ylims ,
      cex.lab=1.5 , main=paste("All ",subjPlotLab,"s",sep="") , cex.main=1.75 )
for ( sIdx in 1:length(sIdVec) ) {
    rVec = ( dataMat[,subjColName] == sIdVec[sIdx] )
    lines( dataMat[rVec,xColName] , dataMat[rVec,yColName] ,
           lty=sIdx , pch=sIdx%%26 , type="o")
}
# Plot histograms of corresponding posterior slopes:
xLim=quantile(c(b1samp,mu1Gsamp),probs=c(0.001,0.999))
for ( sIdx in subjIdxVec ) {
    histInfo = plotPost( b1samp[sIdx,] , xlab="Slope" , compVal=0.0 , 
                     HDItextPlace=0.9 , xlim=xLim )
}
histInfo = plotPost( mu1Gsamp , xlab="Slope, Group Level" , compVal=0.0 ,
                     HDItextPlace=0.9 , xlim=xLim )
saveGraph(file=paste(fileNameRoot,"Data",sep=""),type="eps")

#------------------------------------------------------------------------------
