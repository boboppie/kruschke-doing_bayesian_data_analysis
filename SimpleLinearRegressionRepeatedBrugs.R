graphics.off()
rm(list=ls(all=TRUE))
library(BRugs)         # Kruschke, J. K. (2010). Doing Bayesian data analysis:
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
modelCheck( "model.txt" )

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
fname = "SimpleLinearRegressionRepeatedBrugs"

if ( F ) { # change to T to use income data instead of contam.retention data.
 # Data from http://www.census.gov/hhes/www/income/statemedfaminc.html
 # Downloaded Dec. 06, 2009.
 load("IncomeFamszState.Rdata") # loads IncomeFamszState
 dataMat = IncomeFamszState
 yColName="Income" ; yPlotLab = "Income"
 xColName="Famsz" ; xPlotLab="Family Size"
 subjColName="State" ; subjPlotLab="State"
 fname = "IncomeFamszState"
}

# Extract data info to pass to BUGS:
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
datalist = list(
  Ndata = Ndata ,
  Nsubj = Nsubj ,
  subj = subj ,
  x = zx ,
  y = zy
)
# Get the data into BRugs:
modelData( bugsData( datalist ) )

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

nchain = 3
modelCompile( numChains = nchain )

genInitList <- function() {
    b0 = b1 = tau = rep(0,length=Nsubj)
    for ( sIdx in 1:Nsubj ) {
       yVec = datalist$y[datalist$subj==sIdx]
       xVec = datalist$x[datalist$subj==sIdx]
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
    list( b0=b0 , b1=b1 , tau=tau ,
          mu0G=mu0G , tau0G=tau0G ,
          mu1G=mu1G , tau1G=tau1G ,
          m=m , d=d )
}
for ( chainIdx in 1 : nchain ) {
    modelInits( bugsInits( genInitList ) )
}

#------------------------------------------------------------------------------
# RUN THE CHAINS

# burn in
BurnInSteps = 500
modelUpdate( BurnInSteps )
# actual samples
samplesSet( c( "b0","b1","tau" , "mu0G","tau0G", "mu1G","tau1G", "m","d" ) )
stepsPerChain = ceiling(5000/nchain)
thinStep = 100 # 40 or more
modelUpdate( stepsPerChain , thin=thinStep )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

source("plotChains.R")
source("plotPost.R")

# Check convergence and autocorrelation:
checkConvergence = T  # check this first time through, examine m,d,tau0G,tau1G
if ( checkConvergence ) {
   # check a few selected chains
   b01Sum   = plotChains( "b0[1]"  , saveplots=F , filenameroot=fname )
   b11Sum   = plotChains( "b1[1]"  , saveplots=F , filenameroot=fname )
   tau1Sum  = plotChains( "tau[1]" , saveplots=F , filenameroot=fname )
   mu0GSum  = plotChains( "mu0G"   , saveplots=F , filenameroot=fname )
   tau0GSum = plotChains( "tau0G"  , saveplots=F , filenameroot=fname )
   mu1GSum  = plotChains( "mu1G"   , saveplots=F , filenameroot=fname )
   tau1GSum = plotChains( "tau1G"  , saveplots=F , filenameroot=fname )
   mSum     = plotChains( "m"      , saveplots=F , filenameroot=fname )
   dSum     = plotChains( "d"      , saveplots=F , filenameroot=fname )
}

# Extract chain values for subsequent examination:
zmu0Gsamp = samplesSample( "mu0G" )
zmu1Gsamp = samplesSample( "mu1G" )
zb0samp = NULL
zb1samp = NULL
for ( subjIdx in 1:Nsubj ) {
    zb0samp = rbind( zb0samp , samplesSample( paste("b0[",subjIdx,"]",sep="") ))
    zb1samp = rbind( zb1samp , samplesSample( paste("b1[",subjIdx,"]",sep="") ))
}

# Convert to original scale:
mu0Gsamp = zmu0Gsamp * ySD + yM - zmu1Gsamp * ySD * xM / xSD
mu1Gsamp = zmu1Gsamp * ySD / xSD
b0samp   = zb0samp   * ySD + yM - zb1samp   * ySD * xM / xSD
b1samp   = zb1samp   * ySD / xSD

# Display believable intercept and slope values
windows(10,5.5)
par( mar=c(4,4,1.75,1)+0.1 , mgp=c(2.5,0.8,0) )
layout( matrix(1:2,nrow=1) )
thinIdx = round(seq(1,length(mu0Gsamp),length=700))
plot( zmu1Gsamp[thinIdx] , zmu0Gsamp[thinIdx] , cex.lab=1.75 ,
      ylab="Standardized Intercept" , xlab="Standardized Slope" )
plot( mu1Gsamp[thinIdx] , mu0Gsamp[thinIdx] , cex.lab=1.0 ,
      ylab=paste("Intercept (",yPlotLab," when ",xPlotLab," =0)",sep="") ,
      xlab=paste("Slope (change in",yPlotLab,"per unit",xPlotLab,")") )
dev.copy2eps(file=paste(fname,"SlopeIntercept.eps",sep=""))

# Make graphs of data and corresponding believable slopes:
windows(12,6)
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
for ( sIdx in subjIdxVec ) {
    histInfo = plotPost( b1samp[sIdx,] , xlab="Slope" , compVal=0.0 , breaks=30 ,
                     HDItextPlace=0.9 )
}
histInfo = plotPost( mu1Gsamp , xlab="Slope, Group Level" , compVal=0.0 ,
                     breaks=30 , HDItextPlace=0.9 )
dev.copy2eps(file=paste(fname,"Data.eps",sep=""))

#------------------------------------------------------------------------------
