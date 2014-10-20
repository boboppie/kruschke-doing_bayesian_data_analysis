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
### ******** SEE FILE BESTexample.R FOR INSTRUCTIONS **************
### ***************************************************************

source("openGraphSaveGraph.R") # graphics functions for Windows, Mac OS, Linux

BESTmcmc = function( y1, y2, numSavedSteps=100000, thinSteps=1, showMCMC=FALSE) { 
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
      y[i] ~ dt( mu[x[i]] , tau[x[i]] , nu )
    }
    for ( j in 1:2 ) {
      mu[j] ~ dnorm( muM , muP )
      tau[j] <- 1/pow( sigma[j] , 2 )
      sigma[j] ~ dunif( sigmaLow , sigmaHigh )
    }
    nu <- nuMinusOne+1
    nuMinusOne ~ dexp(1/29)
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="BESTmodel.txt" )
  
  #------------------------------------------------------------------------------
  # THE DATA.
  # Load the data:
  y = c( y1 , y2 ) # combine data into one vector
  x = c( rep(1,length(y1)) , rep(2,length(y2)) ) # create group membership code
  Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    y = y ,
    x = x ,
    Ntotal = Ntotal ,
    muM = mean(y) ,
    muP = 0.000001 * 1/sd(y)^2 ,
    sigmaLow = sd(y) / 1000 ,
    sigmaHigh = sd(y) * 1000 
  )
  
  #------------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  mu = c( mean(y1) , mean(y2) )
  sigma = c( sd(y1) , sd(y2) )
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

BESTsummary = function( y1 , y2 , mcmcChain ) {
  source("HDIofMCMC.R")
  mcmcSummary = function( paramSampleVec , compVal=NULL ) {
    meanParam = mean( paramSampleVec )
    medianParam = median( paramSampleVec )
    dres = density( paramSampleVec )
    modeParam = dres$x[which.max(dres$y)]
    hdiLim = HDIofMCMC( paramSampleVec )
    if ( !is.null(compVal) ) {
      pcgtCompVal = ( 100 * sum( paramSampleVec > compVal ) 
                      / length( paramSampleVec ) )
    } else {
      pcgtCompVal=NA
    }
    return( c( meanParam , medianParam , modeParam , hdiLim , pcgtCompVal ) )
  }
  # Define matrix for storing summary info:
  summaryInfo = matrix( 0 , nrow=9 , ncol=6 , dimnames=list(
    PARAMETER=c( "mu1" , "mu2" , "muDiff" , "sigma1" , "sigma2" , "sigmaDiff" ,
             "nu" , "nuLog10" , "effSz" ),
    SUMMARY.INFO=c( "mean" , "median" , "mode" , "HDIlow" , "HDIhigh" ,
                    "pcgtZero" ) 
    ) )
  summaryInfo[ "mu1" , ] = mcmcSummary( mcmcChain[,"mu[1]"] )
  summaryInfo[ "mu2" , ] = mcmcSummary( mcmcChain[,"mu[2]"] )
  summaryInfo[ "muDiff" , ] = mcmcSummary( mcmcChain[,"mu[1]"]
                                           - mcmcChain[,"mu[2]"] , 
                                           compVal=0 )
  summaryInfo[ "sigma1" , ] = mcmcSummary( mcmcChain[,"sigma[1]"] )
  summaryInfo[ "sigma2" , ] = mcmcSummary( mcmcChain[,"sigma[2]"] )
  summaryInfo[ "sigmaDiff" , ] = mcmcSummary( mcmcChain[,"sigma[1]"]
                                              - mcmcChain[,"sigma[2]"] , 
                                              compVal=0 )
  summaryInfo[ "nu" , ] = mcmcSummary( mcmcChain[,"nu"] )
  summaryInfo[ "nuLog10" , ] = mcmcSummary( log10(mcmcChain[,"nu"]) )
  
  N1 = length(y1)
  N2 = length(y2)
  effSzChain = ( ( mcmcChain[,"mu[1]"] - mcmcChain[,"mu[2]"] ) 
            / sqrt( ( mcmcChain[,"sigma[1]"]^2 + mcmcChain[,"sigma[2]"]^2 ) / 2 ) ) 
  summaryInfo[ "effSz" , ] = mcmcSummary( effSzChain , compVal=0 )
  # Or, use sample-size weighted version:
  # effSz = ( mu1 - mu2 ) / sqrt( ( sigma1^2 *(N1-1) + sigma2^2 *(N2-1) ) 
  #                               / (N1+N2-2) )
  # Be sure also to change plot label in BESTplot function, below.
  return( summaryInfo )
}

#==============================================================================

BESTplot = function( y1 , y2 , mcmcChain , ROPEm=NULL , ROPEsd=NULL , 
                    ROPEeff=NULL , showCurve=FALSE , pairsPlot=FALSE ) {
  # This function plots the posterior distribution (and data).
  # Description of arguments:
  # y1 and y2 are the data vectors.
  # mcmcChain is a list of the type returned by function BTT.
  # ROPEm is a two element vector, such as c(-1,1), specifying the limit
  #   of the ROPE on the difference of means.
  # ROPEsd is a two element vector, such as c(-1,1), specifying the limit
  #   of the ROPE on the difference of standard deviations.
  # ROPEeff is a two element vector, such as c(-1,1), specifying the limit
  #   of the ROPE on the effect size.
  # showCurve is TRUE or FALSE and indicates whether the posterior should
  #   be displayed as a histogram (by default) or by an approximate curve.
  # pairsPlot is TRUE or FALSE and indicates whether scatterplots of pairs
  #   of parameters should be displayed.
  mu1 = mcmcChain[,"mu[1]"]
  mu2 = mcmcChain[,"mu[2]"]
  sigma1 = mcmcChain[,"sigma[1]"]
  sigma2 = mcmcChain[,"sigma[2]"]
  nu = mcmcChain[,"nu"]
  if ( pairsPlot ) {
    # Plot the parameters pairwise, to see correlations:
    openGraph(width=7,height=7)
    nPtToPlot = 1000
    plotIdx = floor(seq(1,length(mu1),by=length(mu1)/nPtToPlot))
    panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
      usr = par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r = (cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex=1.25 ) # was cex=cex.cor*r
    }
    pairs( cbind( mu1 , mu2 , sigma1 , sigma2 , log10(nu) )[plotIdx,] ,
           labels=c( expression(mu[1]) , expression(mu[2]) , 
                     expression(sigma[1]) , expression(sigma[2]) , 
                     expression(log10(nu)) ) , 
           lower.panel=panel.cor , col="skyblue" )
  }
  source("plotPost.R")
  # Set up window and layout:
  openGraph(width=6.0,height=8.0)
  layout( matrix( c(4,5,7,8,3,1,2,6,9,10) , nrow=5, byrow=FALSE ) )
  par( mar=c(3.5,3.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
  
  # Select thinned steps in chain for plotting of posterior predictive curves:
  chainLength = NROW( mcmcChain )
  nCurvesToPlot = 30
  stepIdxVec = seq( 1 , chainLength , floor(chainLength/nCurvesToPlot) )
  xRange = range( c(y1,y2) )
  xLim = c( xRange[1]-0.1*(xRange[2]-xRange[1]) , 
            xRange[2]+0.1*(xRange[2]-xRange[1]) )
  xVec = seq( xLim[1] , xLim[2] , length=200 )
  maxY = max( dt( 0 , df=max(nu[stepIdxVec]) ) /
    min(c(sigma1[stepIdxVec],sigma2[stepIdxVec])) )
  # Plot data y1 and smattering of posterior predictive curves:
  stepIdx = 1
  plot( xVec , dt( (xVec-mu1[stepIdxVec[stepIdx]])/sigma1[stepIdxVec[stepIdx]] , 
                   df=nu[stepIdxVec[stepIdx]] )/sigma1[stepIdxVec[stepIdx]] , 
        ylim=c(0,maxY) , cex.lab=1.75 ,
        type="l" , col="skyblue" , lwd=1 , xlab="y" , ylab="p(y)" , 
        main="Data Group 1 w. Post. Pred." )
  for ( stepIdx in 2:length(stepIdxVec) ) {
    lines(xVec, dt( (xVec-mu1[stepIdxVec[stepIdx]])/sigma1[stepIdxVec[stepIdx]] , 
                      df=nu[stepIdxVec[stepIdx]] )/sigma1[stepIdxVec[stepIdx]] , 
           type="l" , col="skyblue" , lwd=1 )
  }
  histBinWd = median(sigma1)/2
  histCenter = mean(mu1)
  histBreaks = sort( c( seq( histCenter-histBinWd/2 , min(xVec)-histBinWd/2 ,
                             -histBinWd ),
                        seq( histCenter+histBinWd/2 , max(xVec)+histBinWd/2 ,
                             histBinWd ) , xLim ) )
  histInfo = hist( y1 , plot=FALSE , breaks=histBreaks )
  yPlotVec = histInfo$density 
  yPlotVec[ yPlotVec==0.0 ] = NA
  xPlotVec = histInfo$mids
  xPlotVec[ yPlotVec==0.0 ] = NA
  points( xPlotVec , yPlotVec , type="h" , lwd=3 , col="red" )
  text( max(xVec) , maxY , bquote(N[1]==.(length(y1))) , adj=c(1.1,1.1) )
  # Plot data y2 and smattering of posterior predictive curves:
  stepIdx = 1
  plot( xVec , dt( (xVec-mu2[stepIdxVec[stepIdx]])/sigma2[stepIdxVec[stepIdx]] , 
                   df=nu[stepIdxVec[stepIdx]] )/sigma2[stepIdxVec[stepIdx]] , 
        ylim=c(0,maxY) , cex.lab=1.75 , 
        type="l" , col="skyblue" , lwd=1 , xlab="y" , ylab="p(y)" , 
        main="Data Group 2 w. Post. Pred." )
  for ( stepIdx in 2:length(stepIdxVec) ) {
    lines(xVec, dt( (xVec-mu2[stepIdxVec[stepIdx]])/sigma2[stepIdxVec[stepIdx]] , 
                      df=nu[stepIdxVec[stepIdx]] )/sigma2[stepIdxVec[stepIdx]] , 
           type="l" , col="skyblue" , lwd=1 )
  }
  histBinWd = median(sigma2)/2
  histCenter = mean(mu2)
  histBreaks = sort( c( seq( histCenter-histBinWd/2 , min(xVec)-histBinWd/2 ,
                             -histBinWd ),
                        seq( histCenter+histBinWd/2 , max(xVec)+histBinWd/2 ,
                             histBinWd ) , xLim ) )
  histInfo = hist( y2 , plot=FALSE , breaks=histBreaks )
  yPlotVec = histInfo$density 
  yPlotVec[ yPlotVec==0.0 ] = NA
  xPlotVec = histInfo$mids
  xPlotVec[ yPlotVec==0.0 ] = NA
  points( xPlotVec , yPlotVec , type="h" , lwd=3 , col="red" )
  text( max(xVec) , maxY , bquote(N[2]==.(length(y2))) , adj=c(1.1,1.1) )

  # Plot posterior distribution of parameter nu:
  histInfo = plotPost( log10(nu) , col="skyblue" , # breaks=30 ,
                       showCurve=showCurve ,
                  xlab=bquote("log10("*nu*")") , cex.lab = 1.75 , showMode=TRUE ,
                  main="Normality" ) #  (<0.7 suggests kurtosis)

  # Plot posterior distribution of parameters mu1, mu2, and their difference:
  xlim = range( c( mu1 , mu2 ) )
  histInfo = plotPost( mu1 ,  xlim=xlim , cex.lab = 1.75 ,
                       showCurve=showCurve ,
                  xlab=bquote(mu[1]) , main=paste("Group",1,"Mean") , 
                  col="skyblue" )
  histInfo = plotPost( mu2 ,  xlim=xlim , cex.lab = 1.75 ,
                       showCurve=showCurve ,
                  xlab=bquote(mu[2]) , main=paste("Group",2,"Mean") , 
                  col="skyblue" )
  histInfo = plotPost( mu1-mu2 , compVal=0 ,  showCurve=showCurve ,
                  xlab=bquote(mu[1] - mu[2]) , cex.lab = 1.75 , ROPE=ROPEm ,
                  main="Difference of Means" , col="skyblue" )

  # Plot posterior distribution of param's sigma1, sigma2, and their difference:
  xlim=range( c( sigma1 , sigma2 ) )
  histInfo = plotPost( sigma1 ,  xlim=xlim , cex.lab = 1.75 ,
                       showCurve=showCurve ,
                  xlab=bquote(sigma[1]) , main=paste("Group",1,"Std. Dev.") , 
                  col="skyblue" , showMode=TRUE )
  histInfo = plotPost( sigma2 ,  xlim=xlim , cex.lab = 1.75 ,
                       showCurve=showCurve ,
                  xlab=bquote(sigma[2]) , main=paste("Group",2,"Std. Dev.") , 
                  col="skyblue" , showMode=TRUE )
  histInfo = plotPost( sigma1-sigma2 , 
                       compVal=0 ,  showCurve=showCurve ,
                       xlab=bquote(sigma[1] - sigma[2]) , cex.lab = 1.75 , 
                       ROPE=ROPEsd ,
               main="Difference of Std. Dev.s" , col="skyblue" , showMode=TRUE )

  # Plot of estimated effect size. Effect size is d-sub-a from 
  # Macmillan & Creelman, 1991; Simpson & Fitter, 1973; Swets, 1986a, 1986b.
  effectSize = ( mu1 - mu2 ) / sqrt( ( sigma1^2 + sigma2^2 ) / 2 )
  histInfo = plotPost( effectSize , compVal=0 ,  ROPE=ROPEeff ,
                        showCurve=showCurve ,
                  xlab=bquote( (mu[1]-mu[2])
                    /sqrt((sigma[1]^2 +sigma[2]^2 )/2 ) ),
              showMode=TRUE , cex.lab=1.0 , main="Effect Size" , col="skyblue" )
  # Or use sample-size weighted version:
  # Hedges 1981; Wetzels, Raaijmakers, Jakab & Wagenmakers 2009.
  # N1 = length(y1)
  # N2 = length(y2)
  # effectSize = ( mu1 - mu2 ) / sqrt( ( sigma1^2 *(N1-1) + sigma2^2 *(N2-1) )
  #                                    / (N1+N2-2) )
  # Be sure also to change BESTsummary function, above.
  # histInfo = plotPost( effectSize , compVal=0 ,  ROPE=ROPEeff ,
  #          showCurve=showCurve ,
  #          xlab=bquote( (mu[1]-mu[2])
  #          /sqrt((sigma[1]^2 *(N[1]-1)+sigma[2]^2 *(N[2]-1))/(N[1]+N[2]-2)) ),
  #          showMode=TRUE , cex.lab=1.0 , main="Effect Size" , col="skyblue" )
  return( BESTsummary( y1 , y2 , mcmcChain ) )
} # end of function BESTplot

#==============================================================================

BESTpower = function( mcmcChain , N1 , N2 , ROPEm , ROPEsd , ROPEeff ,
                     maxHDIWm , maxHDIWsd , maxHDIWeff , nRep=200 ,
                     mcmcLength=10000 , saveName="BESTpower.Rdata" ,
                     showFirstNrep=0 ) {
  # This function estimates power.
  # Description of arguments:
  # mcmcChain is a matrix of the type returned by function BESTmcmc.
  # N1 and N2 are sample sizes for the two groups.
  # ROPEm is a two element vector, such as c(-1,1), specifying the limit
  #   of the ROPE on the difference of means.
  # ROPEsd is a two element vector, such as c(-1,1), specifying the limit
  #   of the ROPE on the difference of standard deviations.
  # ROPEeff is a two element vector, such as c(-1,1), specifying the limit
  #   of the ROPE on the effect size.
  # maxHDIWm is the maximum desired width of the 95% HDI on the difference 
  #   of means.
  # maxHDIWsd is the maximum desired width of the 95% HDI on the difference 
  #   of standard deviations.
  # maxHDIWeff is the maximum desired width of the 95% HDI on the effect size.
  # nRep is the number of simulated experiments used to estimate the power.
  source("HDIofMCMC.R")
  source("HDIofICDF.R")
  chainLength = NROW( mcmcChain )
  # Select thinned steps in chain for posterior predictions:
  stepIdxVec = seq( 1 , chainLength , floor(chainLength/nRep) )
  goalTally = list( HDIm.gt.ROPE = c(0) ,
                    HDIm.lt.ROPE = c(0) ,
                    HDIm.in.ROPE = c(0) ,
                    HDIm.wd.max = c(0) ,
                    HDIsd.gt.ROPE = c(0) ,
                    HDIsd.lt.ROPE = c(0) ,
                    HDIsd.in.ROPE = c(0) ,
                    HDIsd.wd.max = c(0) ,
                    HDIeff.gt.ROPE = c(0) ,
                    HDIeff.lt.ROPE = c(0) ,
                    HDIeff.in.ROPE = c(0) ,
                    HDIeff.wd.max = c(0) )
  power     = list( HDIm.gt.ROPE = c(0,0,0) ,
                    HDIm.lt.ROPE = c(0,0,0) ,
                    HDIm.in.ROPE = c(0,0,0) ,
                    HDIm.wd.max = c(0,0,0) ,
                    HDIsd.gt.ROPE = c(0,0,0) ,
                    HDIsd.lt.ROPE = c(0,0,0) ,
                    HDIsd.in.ROPE = c(0,0,0) ,
                    HDIsd.wd.max = c(0,0,0) ,
                    HDIeff.gt.ROPE = c(0,0,0) ,
                    HDIeff.lt.ROPE = c(0,0,0) ,
                    HDIeff.in.ROPE = c(0,0,0) ,
                    HDIeff.wd.max = c(0,0,0) )
  nSim = 0
  for ( stepIdx in stepIdxVec ) {
    nSim = nSim + 1
    cat( "\n:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n" )
    cat( paste( "Power computation: Simulated Experiment" , nSim , "of" , 
                length(stepIdxVec) , ":\n\n" ) )
    # Get parameter values for this simulation:
    mu1Val = mcmcChain[stepIdx,"mu[1]"]
    mu2Val = mcmcChain[stepIdx,"mu[2]"]
    sigma1Val = mcmcChain[stepIdx,"sigma[1]"]
    sigma2Val = mcmcChain[stepIdx,"sigma[2]"]
    nuVal = mcmcChain[stepIdx,"nu"]
    # Generate simulated data:
    y1 = rt( N1 , df=nuVal ) * sigma1Val + mu1Val    
    y2 = rt( N2 , df=nuVal ) * sigma2Val + mu2Val    
    # Get posterior for simulated data:
    simChain = BESTmcmc( y1, y2, numSavedSteps=mcmcLength, thinSteps=1, 
                         showMCMC=FALSE )
    if (nSim <= showFirstNrep ) { 
      BESTplot( y1, y2, simChain, ROPEm=ROPEm, ROPEsd=ROPEsd, ROPEeff=ROPEeff ) 
    }
    # Assess which goals were achieved and tally them:
    HDIm = HDIofMCMC( simChain[,"mu[1]"] - simChain[,"mu[2]"] ) 
    if ( HDIm[1] > ROPEm[2] ) { 
      goalTally$HDIm.gt.ROPE = goalTally$HDIm.gt.ROPE + 1 }
    if ( HDIm[2] < ROPEm[1] ) { 
      goalTally$HDIm.lt.ROPE = goalTally$HDIm.lt.ROPE + 1 }
    if ( HDIm[1] > ROPEm[1] &  HDIm[2] < ROPEm[2] ) { 
      goalTally$HDIm.in.ROPE = goalTally$HDIm.in.ROPE + 1 } 
    if ( HDIm[2] - HDIm[1] < maxHDIWm ) { 
      goalTally$HDIm.wd.max = goalTally$HDIm.wd.max + 1 }
    HDIsd = HDIofMCMC( simChain[,"sigma[1]"] - simChain[,"sigma[2]"] ) 
    if ( HDIsd[1] > ROPEsd[2] ) { 
      goalTally$HDIsd.gt.ROPE = goalTally$HDIsd.gt.ROPE + 1 }
    if ( HDIsd[2] < ROPEsd[1] ) { 
      goalTally$HDIsd.lt.ROPE = goalTally$HDIsd.lt.ROPE + 1 }
    if ( HDIsd[1] > ROPEsd[1] &  HDIsd[2] < ROPEsd[2] ) { 
      goalTally$HDIsd.in.ROPE = goalTally$HDIsd.in.ROPE + 1 } 
    if ( HDIsd[2] - HDIsd[1] < maxHDIWsd ) { 
      goalTally$HDIsd.wd.max = goalTally$HDIsd.wd.max + 1 }

    HDIeff = HDIofMCMC( ( simChain[,"mu[1]"] - simChain[,"mu[2]"] ) /
      sqrt( ( simChain[,"sigma[1]"]^2 + simChain[,"sigma[2]"]^2 ) / 2 ) ) 
    if ( HDIeff[1] > ROPEeff[2] ) { 
      goalTally$HDIeff.gt.ROPE = goalTally$HDIeff.gt.ROPE + 1 }
    if ( HDIeff[2] < ROPEeff[1] ) { 
      goalTally$HDIeff.lt.ROPE = goalTally$HDIeff.lt.ROPE + 1 }
    if ( HDIeff[1] > ROPEeff[1] &  HDIeff[2] < ROPEeff[2] ) { 
      goalTally$HDIeff.in.ROPE = goalTally$HDIeff.in.ROPE + 1 } 
    if ( HDIeff[2] - HDIeff[1] < maxHDIWeff ) { 
      goalTally$HDIeff.wd.max = goalTally$HDIeff.wd.max + 1 }

    for ( i in 1:length(goalTally) ) {
      s1 = 1 + goalTally[[i]] 
      s2 = 1 + ( nSim - goalTally[[i]] )
      power[[i]][1] = s1/(s1+s2)
      power[[i]][2:3] = HDIofICDF( qbeta , shape1=s1 , shape2=s2 )
    }
    cat( "\nAfter ", nSim ,
         " Simulated Experiments, Power Is: \n( mean, HDI_low, HDI_high )\n" )
    for ( i in 1:length(power) ) {
      cat( names(power[i]) , round(power[[i]],3) , "\n" )
    }
    save( nSim , power , file=saveName )
  }
  return( power )
} # end of function BESTpower

#==============================================================================

makeData = function( mu1 , sd1 , mu2 , sd2 , nPerGrp , 
                     pcntOut=0 , sdOutMult=2.0 , rnd.seed=NULL ) {
  # Auxilliary function for generating random values from a 
  # mixture of normal distibutions.
  if(!is.null(rnd.seed)){set.seed(rnd.seed)} # Set seed for random values.
  nOut = ceiling(nPerGrp*pcntOut/100)        # Number of outliers.
  nIn = nPerGrp - nOut                       # Number from main distribution.
  #
  y1 = rnorm(n=nIn)                     # Random values in main distribution.
  y1 = ((y1-mean(y1))/sd(y1))*sd1 + mu1 # Translate so exactly realize desired.
  y1Out = rnorm(n=nOut)                 # Random values for outliers
  y1Out=((y1Out-mean(y1Out))/sd(y1Out))*(sd1*sdOutMult)+mu1 # Realize exactly.
  y1 = c(y1,y1Out)                      # Concatenate main with outliers.
  #
  y2 = rnorm(n=nIn)                     # Repeat for second group.
  y2 = ((y2-mean(y2))/sd(y2))*sd2 + mu2
  y2Out = rnorm(n=nOut) 
  y2Out=((y2Out-mean(y2Out))/sd(y2Out))*(sd2*sdOutMult)+mu2
  y2 = c(y2,y2Out)
  #
  # Set up window and layout:
  openGraph(width=7,height=7) # Plot the data.
  layout(matrix(1:2,nrow=2))
  histInfo = hist( y1 , main="Simulated Data" , col="pink2" , border="white" , 
                   xlim=range(c(y1,y2)) , breaks=30 , prob=TRUE )
  text( max(c(y1,y2)) , max(histInfo$density) , 
        bquote(N==.(nPerGrp)) , adj=c(1,1) )
  xlow=min(histInfo$breaks)
  xhi=max(histInfo$breaks)
  xcomb=seq(xlow,xhi,length=1001)
  lines( xcomb , dnorm(xcomb,mean=mu1,sd=sd1)*nIn/(nIn+nOut) +
    dnorm(xcomb,mean=mu1,sd=sd1*sdOutMult)*nOut/(nIn+nOut) , lwd=3 )
  lines( xcomb , dnorm(xcomb,mean=mu1,sd=sd1)*nIn/(nIn+nOut) , 
         lty="dashed" , col="green", lwd=3)
  lines( xcomb , dnorm(xcomb,mean=mu1,sd=sd1*sdOutMult)*nOut/(nIn+nOut) , 
         lty="dashed" , col="red", lwd=3)
  histInfo = hist( y2 , main="" , col="pink2" , border="white" , 
                   xlim=range(c(y1,y2)) , breaks=30 , prob=TRUE )
  text( max(c(y1,y2)) , max(histInfo$density) , 
        bquote(N==.(nPerGrp)) , adj=c(1,1) )
  xlow=min(histInfo$breaks)
  xhi=max(histInfo$breaks)
  xcomb=seq(xlow,xhi,length=1001)
  lines( xcomb , dnorm(xcomb,mean=mu2,sd=sd2)*nIn/(nIn+nOut) +
    dnorm(xcomb,mean=mu2,sd=sd2*sdOutMult)*nOut/(nIn+nOut) , lwd=3)
  lines( xcomb , dnorm(xcomb,mean=mu2,sd=sd2)*nIn/(nIn+nOut) , 
         lty="dashed" , col="green", lwd=3)
  lines( xcomb , dnorm(xcomb,mean=mu2,sd=sd2*sdOutMult)*nOut/(nIn+nOut) , 
         lty="dashed" , col="red", lwd=3)
  #
  return( list( y1=y1 , y2=y2 ) )
}

#==============================================================================
