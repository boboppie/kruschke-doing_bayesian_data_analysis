# Load auxiliary functions for plotting the MCMC chain:
source("plotPost.R")
source("openGraphSaveGraph.R")

BMLRmcmc = function( dataMat , numSavedSteps=250000 , checkConvergence=FALSE ) {
                     # Program written in the style of
  require(rjags)     # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                     # A Tutorial with R and JAGS. Academic Press / Elsevier.
  #------------------------------------------------------------------------------
  # THE MODEL.
  modelstring = "
  # JAGS model specification begins here...
  model {
    # Likelihood:
    for( i in 1:N ) {
      y[i] ~ dnorm( y.hat[i] , tau )
      y.hat[i] <- b0 +  inprod( b[1:nPred] , x[i,1:nPred] )
    }
    # Prior (assumes standardized data):
    tau <- 1/pow(sigma,2)
    sigma ~ dunif( 0 , 10 )
    b0 ~ dnorm( 0 , 1.0E-2 ) 
    for ( j in 1:nPred ) {
      b[j] ~ dnorm( 0 , 1.0E-2 )
    }
  }
  # ... end JAGS model specification
  " # close quote for modelstring
  writeLines(modelstring,con="model.txt")
  
  #------------------------------------------------------------------------------
  # THE DATA.
  
  # dataMat is supplied as an argument to this function.
  # The program assumes that the first column is y, and the remaining columns
  # are x (predictors).
  
  # Now re-name variables from dataMat:
  N = NROW(dataMat)
  y = cbind(as.matrix(dataMat[,1]))
  predictedName = colnames(dataMat)[1]
  x = as.matrix(dataMat[,-1]) 
  predictorNames = colnames(dataMat)[-1]
  nPred = NCOL(x)
  
  # Prepare data for JAGS:
  # Re-center data at mean, to reduce autocorrelation in MCMC sampling.
  # Divide by SD to make prior specification generic.
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
  
  # Get the data into JAGS:
  dataList = list(
             x = zx ,
             y = as.vector( zy ) , # JAGS does not treat 1-column mat as vector
             N = N ,
             nPred = nPred
  )
  
  #------------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  
  # Start the chains at the least-squares fit:
  lmInfo = lm( dataList$y ~ dataList$x )
  initsList = list(
      b0 = lmInfo$coef[1] ,   
      b = lmInfo$coef[-1] ,        
      sigma = sqrt(mean(lmInfo$resid^2)) 
  )
  
  #------------------------------------------------------------------------------
  # RUN THE CHAINS
  
  parameters = c("b0" , "b" , "sigma" )  
  adaptSteps = 1000          # Number of steps to "tune" the samplers.
  burnInSteps = 1000        # Number of steps to "burn-in" the samplers.
  nChains = 3               # Number of chains to run.
  #numSavedSteps=250000       # Total number of steps in chains to save.
  thinSteps=1               # Number of steps to "thin" (1=keep every step).
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
  chainLength = NROW(mcmcChain)
  
  # Rename results, for convenience:
  zsigma = mcmcChain[,"sigma" ]
  zb0 = mcmcChain[,"b0"]
  zb = matrix( 0 , nrow=chainLength , ncol=nPred )
  for ( j in 1:nPred ) {
    zb[,j] = mcmcChain[,paste("b[",j,"]",sep="")]
  }
  colnames(zb) = paste("zb",1:nPred,sep="")
    
  # Convert to original scale (see book Eqn 17.1):
  sigma = zsigma * sd(as.vector(y))
  # for b0:
  subtractTerm = rep(0,chainLength)
  for ( j in 1:nPred ) {
    subtractTerm = subtractTerm + zb[,j] * sd(as.vector(y)) * mean(x[,j]) / sd(x[,j])
  }
  b0 = zb0 * sd(as.vector(y)) + mean(y) - subtractTerm 
  # for b:
  b = 0*zb 
  for ( j in 1:nPred ) {
    b[,j] = zb[,j] * sd(as.vector(y)) / sd(x[,j])
  }
  colnames(b) = paste("b",1:nPred,sep="")
  
  # Compute R^2 for credible parameters:
  Rsq = rep(0,chainLength)
  YcorX = cor( y , x ) # correlation of y with each x predictor
  for ( stepIdx in 1:chainLength ) {
    Rsq[stepIdx] = sum( zb[stepIdx,] * YcorX )
  }
  ## Old version, perhaps useless:
  #Rsq = rep(0,chainLength)
  #for ( stepIdx in 1:chainLength ) {
  #  predY = dataList$x %*% cbind(zb[stepIdx,]) + zb0[stepIdx]
  #  Rsq[stepIdx] = cor(dataList$y,predY)^2
  #}
  
  # Combine results into one big matrix:
  mcmcChain = cbind( zsigma,zb0,zb, sigma,b0,b, Rsq )
  return( mcmcChain )
}


BMLRplot = function( mcmcChain , ROPEbeta=NULL , ROPEbetaDiff=NULL ) {
    
  summaryMatrix = matrix( NA , nrow=0 , ncol=11 , 
                          dimnames=list(ParameterName=NULL,PosteriorInfo=NULL) )
  
  # Compute number of predictors from NCOL of mcmcChain matrix. 
  # Assumes mcmcChain has columns zsigma,zb0,zb,sigma,b0,b,Rsq
  nPred = (NCOL(mcmcChain)-1)/2 - 2
  
  # Display scatter plots of parameter values, pairwise:
  openGraph(width=7,height=7)
  chainLength = NROW(mcmcChain)
  nToPlot = 700
  plotIdx = seq(1,chainLength,by=ceiling(chainLength/nToPlot))
  pairs( mcmcChain[plotIdx,c("Rsq","zsigma","zb0",paste("zb",1:nPred,sep=""))] ,  
         col="skyblue" )
  openGraph(width=7,height=7)
  pairs( mcmcChain[plotIdx,c("Rsq","sigma","b0",paste("b",1:nPred,sep=""))] ,  
         col="skyblue" )
  
  #   # Show correlation matrix on console:
  #   cat("\nCorrlations of posterior sigma, b0, and b's:\n")
  #   show( round( cor( mcmcChain ) , 3) )
  #   cat("\nCovariances of posterior sigma, b0, and b's:\n")
  #   show( round( cov( mcmcChain ) , 5) )

  # Display the marginals of posterior:
  nPostCol = 3 # arbitrary number of columns for display
  nPostRow = ceiling((3+nPred)/nPostCol)
  openGraph(width=nPostCol*2.5,height=nPostRow*2.0)
  layout(matrix(1:(nPostCol*nPostRow),nrow=nPostRow,byrow=TRUE))
  par( mar=c(4,2,3,1) , mgp=c(2.5,0.7,0) )
  histInfo = plotPost( mcmcChain[,"Rsq"] , xlab=bquote(R^2) , main="" ,
                       cex.main=1.67 , cex.lab=1.67 , col="skyblue" )
  summaryMatrix = rbind( summaryMatrix , histInfo )
  histInfo = plotPost( mcmcChain[,"zsigma"] , xlab=bquote(sigma[z]) , main="" ,
                       cex.main=1.67 , cex.lab=1.67 , col="skyblue" )
  summaryMatrix = rbind( summaryMatrix , histInfo )
  histInfo = plotPost( mcmcChain[,"zb0"] , xlab=bquote(beta[0]) , main="" ,
                       cex.main=1.67 , cex.lab=1.67 , col="skyblue" )
  summaryMatrix = rbind( summaryMatrix , histInfo )
  for ( j in 1:nPred ) {
    histInfo = plotPost( mcmcChain[,paste("zb",j,sep="")] , 
                         xlab=bquote(beta[.(j)]) , main="" ,
                         compVal = 0.0 , ROPE=ROPEbeta ,
                         cex.main=1.67 , cex.lab=1.67 , col="skyblue" )
    summaryMatrix = rbind( summaryMatrix , histInfo )
  }
  #
  nPostCol = 3 # arbitrary number of columns for display
  nPostRow = ceiling((3+nPred)/nPostCol)
  openGraph(width=nPostCol*2.5,height=nPostRow*2.0)
  layout(matrix(1:(nPostCol*nPostRow),nrow=nPostRow,byrow=TRUE))
  par( mar=c(4,2,3,1) , mgp=c(2.5,0.7,0) )
  histInfo = plotPost( mcmcChain[,"Rsq"] , xlab=bquote(R^2) , main="" ,
                       cex.main=1.67 , cex.lab=1.67 , col="skyblue" )
  #summaryMatrix = rbind( summaryMatrix , histInfo )
  histInfo = plotPost( mcmcChain[,"sigma"] , xlab=bquote(sigma[y]) , main="" ,
                       cex.main=1.67 , cex.lab=1.67 , col="skyblue" )
  summaryMatrix = rbind( summaryMatrix , histInfo )
  histInfo = plotPost( mcmcChain[,"b0"] , xlab=bquote(b[0]) , main="" ,
                       cex.main=1.67 , cex.lab=1.67 , col="skyblue" )
  summaryMatrix = rbind( summaryMatrix , histInfo )
  for ( j in 1:nPred ) {
  histInfo = plotPost( mcmcChain[,paste("b",j,sep="")] , 
                       xlab=bquote(b[.(j)]) , main="" ,
                       compVal = 0.0 ,
                       cex.main=1.67 , cex.lab=1.67 , col="skyblue" )
  summaryMatrix = rbind( summaryMatrix , histInfo )
  }
  
  # Display differences of standardized regression coefficients:
  nPostCol = 3
  nPostCell = nPred*(nPred-1)/2
  nPostRow = ceiling(nPostCell/nPostCol)
  openGraph(width=nPostCol*2.5,height=nPostRow*2.0)
  layout(matrix(1:(nPostCol*nPostRow),nrow=nPostRow,byrow=TRUE))
  par( mar=c(4,2,3,1) , mgp=c(2.5,0.7,0) )
  for ( j in 1:(nPred-1) ) {
    for ( k in (j+1):nPred ) {
      histInfo = plotPost( mcmcChain[,paste("zb",j,sep="")]
                           - mcmcChain[,paste("zb",k,sep="")], 
                           xlab=bquote(beta[.(j)] - beta[.(k)]) , main="" ,
                           compVal = 0.0 , ROPE=ROPEbetaDiff ,
                           cex.main=1.67 , cex.lab=1.67 , col="skyblue" )
      summaryMatrix = rbind( summaryMatrix , histInfo )
    }
  }
  
  return( summaryMatrix )
  
} # end of function BMLRplot
