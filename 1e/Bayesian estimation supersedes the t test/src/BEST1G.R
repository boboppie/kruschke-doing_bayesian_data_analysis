# MODIFIED FROM BEST.R FOR ONE GROUP INSTEAD OF TWO.
# Version of January 19, 2013.
# John K. Kruschke  
# johnkruschke@gmail.com
# http://www.indiana.edu/~kruschke/BEST/
#
# This program is believed to be free of errors, but it comes with no guarantee! 
# The user bears all responsibility for interpreting the results. 
# Please check the webpage above for updates or corrections.
#
### ***************************************************************
### ******** SEE FILE BEST1Gexample.R FOR INSTRUCTIONS ************
### ***************************************************************

source("openGraphSaveGraph.R") # graphics functions for Windows, Mac OS, Linux

BEST1Gmcmc = function( y , numSavedSteps=100000, thinSteps=1, showMCMC=FALSE) { 
  # This function generates an MCMC sample from the posterior distribution.
  # Description of arguments:
  # showMCMC is a flag for displaying diagnostic graphs of the chains.
  #    If F (the default), no chain graphs are displayed. If T, they are.
  
  require(rjags)

  #------------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  model {
    for ( i in 1:Ntotal ) {
      y[i] ~ dt( mu , tau , nu )
    }
    mu ~ dnorm( muM , muP )
    tau <- 1/pow( sigma , 2 )
    sigma ~ dunif( sigmaLow , sigmaHigh )
    nu <- nuMinusOne+1
    nuMinusOne ~ dexp(1/29)
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="BESTmodel.txt" )
  
  #------------------------------------------------------------------------------
  # THE DATA.
  # Load the data:
  Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    y = y ,
    Ntotal = Ntotal ,
    muM = mean(y) ,
    muP = 0.000001 * 1/sd(y)^2 ,
    sigmaLow = sd(y) / 1000 ,
    sigmaHigh = sd(y) * 1000 
  )
  
  #------------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  mu = mean(y) 
  sigma = sd(y) 
  # Regarding initial values in next line: (1) sigma will tend to be too big if 
  # the data have outliers, and (2) nu starts at 5 as a moderate value. These
  # initial values keep the burn-in period moderate.
  initsList = list( mu = mu , sigma = sigma , nuMinusOne = 4 )
  
  #------------------------------------------------------------------------------
  # RUN THE CHAINS

  parameters = c( "mu" , "sigma" , "nu" )     # The parameters to be monitored
  adaptSteps = 500               # Number of steps to "tune" the samplers
  burnInSteps = 1000
  nChains = 3 
  nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
  # Create, initialize, and adapt the model:
  jagsModel = jags.model( "BESTmodel.txt" , data=dataList , inits=initsList , 
                          n.chains=nChains , n.adapt=adaptSteps )
  # Burn-in:
  cat( "Burning in the MCMC chain...\n" )
  update( jagsModel , n.iter=burnInSteps )
  # The saved MCMC chain:
  cat( "Sampling final MCMC chain...\n" )
  codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                              n.iter=nIter , thin=thinSteps )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  
  #------------------------------------------------------------------------------
  # EXAMINE THE RESULTS
  if ( showMCMC ) {
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
  return( mcmcChain )

} # end function BESTmcmc

#==============================================================================

BEST1Gplot = function( y , mcmcChain , compValm=0.0 , ROPEm=NULL , 
                       compValsd=NULL , ROPEsd=NULL , 
                       ROPEeff=NULL , showCurve=FALSE , pairsPlot=FALSE ) {
  # This function plots the posterior distribution (and data).
  # Description of arguments:
  # y is the data vector.
  # mcmcChain is a list of the type returned by function BEST1G.
  # compValm is hypothesized mean for comparing with estimated mean.
  # ROPEm is a two element vector, such as c(-1,1), specifying the limit
  #   of the ROPE on the mean.
  # ROPEsd is a two element vector, such as c(14,16), specifying the limit
  #   of the ROPE on the difference of standard deviations.
  # ROPEeff is a two element vector, such as c(-1,1), specifying the limit
  #   of the ROPE on the effect size.
  # showCurve is TRUE or FALSE and indicates whether the posterior should
  #   be displayed as a histogram (by default) or by an approximate curve.
  # pairsPlot is TRUE or FALSE and indicates whether scatterplots of pairs
  #   of parameters should be displayed.
  mu = mcmcChain[,"mu"]
  sigma = mcmcChain[,"sigma"]
  nu = mcmcChain[,"nu"]
  
  if ( pairsPlot ) {
    # Plot the parameters pairwise, to see correlations:
    openGraph(width=7,height=7)
    nPtToPlot = 1000
    plotIdx = floor(seq(1,length(mu),by=length(mu)/nPtToPlot))
    panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
      usr = par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r = (cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex=1.25 ) # was cex=cex.cor*r
    }
    pairs( cbind( mu , sigma , log10(nu) )[plotIdx,] ,
           labels=c( expression(mu) ,  
                     expression(sigma) , 
                     expression(log10(nu)) ) , 
           lower.panel=panel.cor , col="skyblue" )
  }

  # Define matrix for storing summary info:
  summaryInfo = NULL
  
  source("plotPost.R")
  # Set up window and layout:
  openGraph(width=6.0,height=5.0)
  layout( matrix( c(3,3,4,4,5,5, 1,1,1,1,2,2) , nrow=6, ncol=2 , byrow=FALSE ) )
  par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  
  # Select thinned steps in chain for plotting of posterior predictive curves:
  chainLength = NROW( mcmcChain )
  nCurvesToPlot = 30
  stepIdxVec = seq( 1 , chainLength , floor(chainLength/nCurvesToPlot) )
  xRange = range( y )
  xLim = c( xRange[1]-0.1*(xRange[2]-xRange[1]) , 
            xRange[2]+0.1*(xRange[2]-xRange[1]) )
  xVec = seq( xLim[1] , xLim[2] , length=200 )
  maxY = max( dt( 0 , df=max(nu[stepIdxVec]) ) / min(sigma[stepIdxVec]) )
  # Plot data y and smattering of posterior predictive curves:
  stepIdx = 1
  plot( xVec , dt( (xVec-mu[stepIdxVec[stepIdx]])/sigma[stepIdxVec[stepIdx]] , 
                   df=nu[stepIdxVec[stepIdx]] )/sigma[stepIdxVec[stepIdx]] , 
        ylim=c(0,maxY) , cex.lab=1.75 ,
        type="l" , col="skyblue" , lwd=1 , xlab="y" , ylab="p(y)" , 
        main="Data w. Post. Pred." )
  for ( stepIdx in 2:length(stepIdxVec) ) {
    lines(xVec, dt( (xVec-mu[stepIdxVec[stepIdx]])/sigma[stepIdxVec[stepIdx]] , 
                    df=nu[stepIdxVec[stepIdx]] )/sigma[stepIdxVec[stepIdx]] , 
          type="l" , col="skyblue" , lwd=1 )
  }
  histBinWd = median(sigma)/2
  histCenter = mean(mu)
  histBreaks = sort( c( seq( histCenter-histBinWd/2 , min(xVec)-histBinWd/2 ,
                             -histBinWd ),
                        seq( histCenter+histBinWd/2 , max(xVec)+histBinWd/2 ,
                             histBinWd ) , xLim ) )
  histInfo = hist( y , plot=FALSE , breaks=histBreaks )
  yPlotVec = histInfo$density 
  yPlotVec[ yPlotVec==0.0 ] = NA
  xPlotVec = histInfo$mids
  xPlotVec[ yPlotVec==0.0 ] = NA
  points( xPlotVec , yPlotVec , type="h" , lwd=3 , col="red" )
  text( max(xVec) , maxY , bquote(N==.(length(y))) , adj=c(1.1,1.1) )
  
  # Plot posterior distribution of parameter nu:
  paramInfo = plotPost( log10(nu) , col="skyblue" , # breaks=30 ,
                       showCurve=showCurve ,
                       xlab=bquote("log10("*nu*")") , cex.lab = 1.75 , showMode=TRUE ,
                       main="Normality" ) #  (<0.7 suggests kurtosis)
  summaryInfo = rbind( summaryInfo , paramInfo )
  rownames(summaryInfo)[NROW(summaryInfo)] = "log10(nu)"
  
  # Plot posterior distribution of parameters mu1, mu2, and their difference:
  xlim = range( c( mu , compValm ) )
  paramInfo = plotPost( mu ,  xlim=xlim , cex.lab = 1.75 , ROPE=ROPEm ,
                       showCurve=showCurve , compVal = compValm ,
                       xlab=bquote(mu) , main=paste("Mean") , 
                       col="skyblue" )
  summaryInfo = rbind( summaryInfo , paramInfo )
  rownames(summaryInfo)[NROW(summaryInfo)] = "mu"
  
  # Plot posterior distribution of sigma:
  xlim=range( sigma )
  paramInfo = plotPost( sigma ,  xlim=xlim , cex.lab = 1.75 , ROPE=ROPEsd ,
                       showCurve=showCurve , compVal=compValsd ,
                       xlab=bquote(sigma) , main=paste("Std. Dev.") , 
                       col="skyblue" , showMode=TRUE )
  summaryInfo = rbind( summaryInfo , paramInfo )
  rownames(summaryInfo)[NROW(summaryInfo)] = "sigma"
  
  # Plot of estimated effect size:
  effectSize = ( mu - compValm ) / sigma
  paramInfo = plotPost( effectSize , compVal=0 ,  ROPE=ROPEeff ,
                       showCurve=showCurve ,
                       xlab=bquote( (mu-.(compValm)) / sigma ),
                       showMode=TRUE , cex.lab=1.75 , main="Effect Size" , col="skyblue" )
  summaryInfo = rbind( summaryInfo , paramInfo )
  rownames(summaryInfo)[NROW(summaryInfo)] = "effSz"
  
  return( summaryInfo )

  } # end of function BEST1Gplot

#==============================================================================

