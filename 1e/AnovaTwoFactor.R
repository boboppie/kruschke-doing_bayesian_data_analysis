# If you use or modify this program, please retain this header.
# USER BEWARE! This program can yield seemingly strange results because of 
# extreme shrinkage of estimates produced by the hierarchical model. That's 
# not wrong, it's merely what the model implies, given the data. The priors 
# are set to try to prevent implosive shrinkage, but that might not be 
# appropriate for many applications. It is the user's responsibility to 
# interpret the results and adjust the model appropriately.
# John Kruschke, 13 April 2013. Programmed in the style of:
# Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
# A Tutorial with R and BUGS. Academic Press / Elsevier.
# This program is analogous to two-factor between-subjects ANOVA, but:
# * It is "robust" (i.e., accommodates outliers) by using a t distribution
#   instead of normal distribution to describe data in each cell.
# * It uses a distinct variance in each cell of the design, instead of 
#   assuming homogeneous variances.
# * It is hierarchical, using higher-level distributions for the marginal
#   means, the interaction of means, and the cell variances (without factor
#   structure for the cell variances), thereby providing shrinkage on the 
#   estimates of those parameters.
# * Its predecessor, ANOVAtwowayJagsSTZ.R, standaradized the y values.
#   This program does not, instead scaling the prior to be vague on the
#   scale of the data. 

graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="AnovaTwoFactor" # for constructing output filenames
source("openGraphSaveGraph.R")
source("plotPost.R")
require(rjags)         

#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
model {
  for ( i in 1:Ntotal ) {
    y[i] ~ dt( mu[i] , 1/sigma[x1[i],x2[i]]^2 , nu )
    mu[i] <- a0 + a1[x1[i]] + a2[x2[i]] + a1a2[x1[i],x2[i]]
  }
  # 
  # For sparse data with lots of outliers, there can be multimodal small-nu
  # estimates, in which case you may want to change the prior to force a 
  # larger value of nu, such as 
  # nuMinusOne ~ dgamma(5.83,0.0483) # mode 100, sd 50
  nu <- nuMinusOne+1
  nuMinusOne ~  dexp(1/29) 
  #
  # For small samples or near-equal variances, shrinkage can produce 
  # an implosion at sigmaSD=0, yielding strong peaks at sigmaSD=0 or
  # bimodal estimates with one mode at sigmaSD=0. That's merely the model
  # telling you the correct answer under the model assumptions. If the shrinkage
  # is not appropriate for your purposes, try a different prior on sigmaSD,
  # or remove the hierarchical prior altogether.
  for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
    sigma[j1,j2] ~ dgamma( sigmaSh , sigmaRa )
  } }
  sigmaSh <- 1 + sigmaMode * sigmaRa
  sigmaRa <- ( sigmaMode + sqrt( sigmaMode^2 + 4*sigmaSD^2 ) ) /(2*sigmaSD^2)
  sigmaMode ~ dgamma(a1a2sigmaModegammaShRa[1],a1a2sigmaModegammaShRa[2]) 
  sigmaSD ~ dgamma(a1a2sigmaSDgammaShRa[1],a1a2sigmaSDgammaShRa[2]) 
  #
  a0 ~ dnorm(dataM,1/(dataSD*10)^2) 
  #
  for ( j1 in 1:Nx1Lvl ) { a1[j1] ~ dnorm( 0.0 , 1/a1SD^2 ) }
  a1SD ~ dgamma(a1gammaShRa[1],a1gammaShRa[2]) # or try a folded t (Cauchy)
  #
  for ( j2 in 1:Nx2Lvl ) { a2[j2] ~ dnorm( 0.0 , 1/a2SD^2 ) }
  a2SD ~ dgamma(a2gammaShRa[1],a2gammaShRa[2]) # or try a folded t (Cauchy)
  #
  # For small samples or near-equal interaction deflections, shrinkage can 
  # produce an implosion at a1a2SD=0, yielding strong peaks at a1a2SD=0 or
  # bimodal estimates with one mode at a1a2SD=0. That's merely the model
  # telling you the correct answer under the model assumptions. If the shrinkage
  # is not appropriate for your purposes, try a different prior on a1a2SD,
  # or remove the hierarchical prior altogether.
  for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
    a1a2[j1,j2] ~ dnorm( 0.0 , 1/a1a2SD^2 )
  } }
  a1a2SD ~ dgamma(a1a2gammaShRa[1],a1a2gammaShRa[2]) # or try a folded t (Cauchy)

  # Convert a0,a1[],a2[],a1a2[,] to sum-to-zero b0,b1[],b2[],b1b2[,] :
  for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
    m[j1,j2] <- a0 + a1[j1] + a2[j2] + a1a2[j1,j2] # cell means 
  } }
  b0 <- mean( m[1:Nx1Lvl,1:Nx2Lvl] )
  for ( j1 in 1:Nx1Lvl ) { b1[j1] <- mean( m[j1,1:Nx2Lvl] ) - b0 }
  for ( j2 in 1:Nx2Lvl ) { b2[j2] <- mean( m[1:Nx1Lvl,j2] ) - b0 }
  for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
    b1b2[j1,j2] <- m[j1,j2] - ( b0 + b1[j1] + b2[j2] )  
  } }
}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.

# Read in the data file:
datarecord = read.table( "AnovaTwoFactorData.csv" , header=TRUE , sep="," )
fileNameRoot = paste( fileNameRoot , "-Random-" , sep="" ) # for output filenames
# Convert data file columns to generic x,y variable names for model and graphics:
y = as.numeric(datarecord$Score)
x1 = as.numeric(datarecord$FactorA)
x1names = levels(datarecord$FactorA)
x2 = as.numeric(datarecord$FactorB)
x2names = levels(datarecord$FactorB)
Ntotal = length(y)
Nx1Lvl = length(unique(x1))
Nx2Lvl = length(unique(x2))
# Specify desired contrasts to examine:
normalize = function( v ){ return( v / sum(v) ) }
x1contrastList = list( 
  A1vA2 = (x1names=="A1")-(x1names=="A2") ,
  A1vA3 = (x1names=="A1")-(x1names=="A3")
)
x2contrastList = list( 
  B1vB2 = (x2names=="B1")-(x2names=="B2") ,
  B1B2vB3B4 = ( normalize((x2names=="B1")|(x2names=="B2")) -
                normalize((x2names=="B3")|(x2names=="B4") )
  )
)
x1x2contrastList = list(
  A1vA2xB1vB2 = outer( 
    (x1names=="A1")-(x1names=="A2") , 
    (x2names=="B1")-(x2names=="B2") 
  ) 
) 

# Compute scale properties of data, for passing into prior to make the prior
# vague on the scale of the data. 
# First, a utility function for gamma-shaped hyper-priors, that returns 
# shape and rate values for desired mode and SD of the gamma distribution:
GammaShRaFromModeSd = function( mode , stndv ) {
  ra = ( mode + sqrt( mode^2 + 4*stndv^2 ) ) / ( 2 * stndv^2 )
  sh = 1 + mode * ra
  return( c(sh,ra) )
}
# For prior on baseline, etc.:
dataM = mean(y)
dataSD = sd(y)
# For hyper-prior on Factor 1 deflections:
a1eff = aggregate( y , list( x1 ) , mean )[,2] - dataM
a1effSD = sd(a1eff)
a1gammaShRa = GammaShRaFromModeSd( mode=a1effSD , stndv=a1effSD ) 
# For hyper-prior on Factor 2 deflections:
a2eff = aggregate( y , list( x2 ) , mean )[,2] - dataM
a2effSD = sd(a2eff)
a2gammaShRa = GammaShRaFromModeSd( mode=a2effSD , stndv=a2effSD )
# For hyper-prior on interaction deflections:
linpred = as.vector( outer( a1eff , a2eff , "+" ) + dataM )
a1a2eff = aggregate( y, list(x1,x2), mean)[,3] - linpred
a1a2effSD = sd(a1a2eff)
a1a2gammaShRa = GammaShRaFromModeSd( mode=a1a2effSD , stndv=a1a2effSD )
# For hyper-prior on cell sigma's:
a1a2sigma = aggregate( y, list(x1,x2), sd)[,3] 
a1a2sigmaMean = mean(a1a2sigma)
a1a2sigmaSD = sd(a1a2sigma)
a1a2sigmaModegammaShRa = GammaShRaFromModeSd( mode=a1a2sigmaMean , stndv=a1a2sigmaSD )
a1a2sigmaSDgammaShRa = GammaShRaFromModeSd( mode=a1a2sigmaSD , stndv=a1a2sigmaSD )

# Specify the data in a list for sending to JAGS:
dataList = list(
  y = y ,
  x1 = x1 ,
  x2 = x2 ,
  Ntotal = Ntotal ,
  Nx1Lvl = Nx1Lvl ,
  Nx2Lvl = Nx2Lvl ,
  # data properties for scaling the prior:
  dataM = dataM ,
  dataSD = dataSD ,
  a1gammaShRa = a1gammaShRa ,
  a2gammaShRa = a2gammaShRa ,
  a1a2gammaShRa = a1a2gammaShRa ,
  a1a2sigmaModegammaShRa = a1a2sigmaModegammaShRa ,
  a1a2sigmaSDgammaShRa = a1a2sigmaSDgammaShRa
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.
# Let JAGS do it automatically

#------------------------------------------------------------------------------
# RUN THE CHAINS

parameters = c( "a0" ,  "a1" ,  "a2" ,  "a1a2" , # don't really need a0,a1, etc.
                "b0" ,  "b1" ,  "b2" ,  "b1b2" ,
                "sigma" , "sigmaMode" , "sigmaSD" ,
                "nu" , "a1SD" , "a2SD" , "a1a2SD" , "m" )
adaptSteps = 1000 
burnInSteps = 2000 
nChains = 4 
numSavedSteps=50000 
thinSteps=1
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
# Create, initialize, and adapt the model:
jagsModel = jags.model( "model.txt" , data=dataList , #inits=initsList , 
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

checkConvergence = FALSE
if ( checkConvergence ) {
  openGraph(width=7,height=7)
  autocorr.plot( codaSamples[[1]][,c("b1[1]","b2[1]","b1b2[1,1]","sigma[1,1]")] )
  show( gelman.diag( codaSamples ) )
  effectiveChainLength = effectiveSize( codaSamples ) 
  show( effectiveChainLength )
}

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )
chainLength = NROW(mcmcChain)

savePlotImages = TRUE 
savePlotType = "jpg"

# Plot the SDs and nu:
openGraph(width=7,height=7)
layout( matrix(1:6,nrow=3,byrow=TRUE) )
par( mar=c(3,1,2.5,0) , mgp=c(2,0.7,0) )
postInfo = plotPost( mcmcChain[,"sigmaMode"], xlab="sigma mode", 
                     main="Modal Cell sigma's", showMode=TRUE )
postInfo = plotPost( mcmcChain[,"sigmaSD"], xlab="sigma SD", 
                     main="SD of Cell sigma's", showMode=TRUE )
postInfo = plotPost( mcmcChain[,"a1SD"] , xlab="a1SD" , main="a1 SD" , showMode=TRUE )
postInfo = plotPost( mcmcChain[,"a2SD"] , xlab="a2SD" , main="a2 SD" , showMode=TRUE )
postInfo = plotPost( mcmcChain[,"a1a2SD"] , xlab="a1a2SD" , main="Interaction SD" ,
                     showMode=TRUE )
postInfo = plotPost( log10(mcmcChain[,"nu"]) , xlab="log10(nu)" , main="Normality" ,
                     showMode=TRUE )
if ( savePlotImages ) {
  saveGraph( file=paste(fileNameRoot,"SD",sep="") , type=savePlotType )
}

# Plot the cell sigma's:
if ( Nx1Lvl <=6 & Nx2Lvl <=6 ) {
  openGraph(width=(dataList$Nx1Lvl)*2.75,height=(dataList$Nx2Lvl)*2.0)
  layoutMat = matrix( 1:(dataList$Nx1Lvl*dataList$Nx2Lvl) , byrow=TRUE ,
                      nrow=(dataList$Nx2Lvl) , ncol=(dataList$Nx1Lvl) )
  layout( layoutMat )
  par( mar=c(4,0.5,2.5,0.5) , mgp=c(2,0.7,0) )
  for ( x2idx in 1:dataList$Nx2Lvl ) {
    for ( x1idx in 1:dataList$Nx1Lvl ) {
      postInfo = plotPost( mcmcChain[,paste("sigma[",x1idx,",",x2idx,"]",sep="")] , 
                           xlab=bquote(sigma[.(x1idx)*","*.(x2idx)]) , showMode=TRUE ,
                           main=paste("x1:",x1names[x1idx],", x2:",x2names[x2idx])  )
    }
  }
  if ( savePlotImages ) {
    saveGraph( file=paste(fileNameRoot,"sigma",sep="") , type=savePlotType )
  }
}

# Plot the means:
if ( Nx1Lvl <=6 & Nx2Lvl <=6 ) {
  openGraph(width=(dataList$Nx1Lvl+1)*2.75,height=(dataList$Nx2Lvl+1)*2.0)
  layoutMat = matrix( 0 , nrow=(dataList$Nx2Lvl+1) , ncol=(dataList$Nx1Lvl+1) )
  layoutMat[1,1] = 1
  layoutMat[1,2:(dataList$Nx1Lvl+1)] = 1:dataList$Nx1Lvl + 1
  layoutMat[2:(dataList$Nx2Lvl+1),1] = 1:dataList$Nx2Lvl + (dataList$Nx1Lvl + 1)
  layoutMat[2:(dataList$Nx2Lvl+1),2:(dataList$Nx1Lvl+1)] = matrix(
    1:(dataList$Nx1Lvl*dataList$Nx2Lvl) + (dataList$Nx2Lvl+dataList$Nx1Lvl+1) ,
    ncol=dataList$Nx1Lvl , byrow=T )
  layout( layoutMat )
  par( mar=c(4,0.5,2.5,0.5) , mgp=c(2,0.7,0) )
  postInfo = plotPost( mcmcChain[,"b0"] , xlab=expression(beta * 0) , main="Baseline" )
  for ( x1idx in 1:dataList$Nx1Lvl ) {
    postInfo = plotPost( mcmcChain[,"b0"]+mcmcChain[,paste("b1[",x1idx,"]",sep="")] ,
                         xlab=bquote(m*1[.(x1idx)]) ,
                         main=paste("x1:",x1names[x1idx]) )
  }
  for ( x2idx in 1:dataList$Nx2Lvl ) {
    postInfo = plotPost( mcmcChain[,"b0"]+mcmcChain[,paste("b2[",x2idx,"]",sep="")] ,
                         xlab=bquote(m*2[.(x2idx)]) ,
                         main=paste("x2:",x2names[x2idx]) )
  }
  for ( x2idx in 1:dataList$Nx2Lvl ) {
    for ( x1idx in 1:dataList$Nx1Lvl ) {
      postInfo = plotPost( mcmcChain[,paste("m[",x1idx,",",x2idx,"]",sep="")] , 
                           xlab=bquote(m*12[.(x1idx)*","*.(x2idx)]) ,
                           main=paste("x1:",x1names[x1idx],", x2:",x2names[x2idx])  )
    }
  }
  if ( savePlotImages ) {
    saveGraph( file=paste(fileNameRoot,"m",sep="") , type=savePlotType )
  }
}

# Plot b values:
if ( Nx1Lvl <=6 & Nx2Lvl <=6 ) {
  openGraph(width=(dataList$Nx1Lvl+1)*2.75,height=(dataList$Nx2Lvl+1)*2.0)
  layoutMat = matrix( 0 , nrow=(dataList$Nx2Lvl+1) , ncol=(dataList$Nx1Lvl+1) )
  layoutMat[1,1] = 1
  layoutMat[1,2:(dataList$Nx1Lvl+1)] = 1:dataList$Nx1Lvl + 1
  layoutMat[2:(dataList$Nx2Lvl+1),1] = 1:dataList$Nx2Lvl + (dataList$Nx1Lvl + 1)
  layoutMat[2:(dataList$Nx2Lvl+1),2:(dataList$Nx1Lvl+1)] = matrix(
    1:(dataList$Nx1Lvl*dataList$Nx2Lvl) + (dataList$Nx2Lvl+dataList$Nx1Lvl+1) ,
    ncol=dataList$Nx1Lvl , byrow=T )
  layout( layoutMat )
  par( mar=c(4,0.5,2.5,0.5) , mgp=c(2,0.7,0) )
  postInfo = plotPost( mcmcChain[,"b0"] , xlab=expression(beta * 0) , main="Baseline" )
  for ( x1idx in 1:dataList$Nx1Lvl ) {
    postInfo = plotPost( mcmcChain[,paste("b1[",x1idx,"]",sep="")] ,
                         xlab=bquote(beta*1[.(x1idx)]) ,
                         main=paste("x1:",x1names[x1idx]) )
  }
  for ( x2idx in 1:dataList$Nx2Lvl ) {
    postInfo = plotPost( mcmcChain[,paste("b2[",x2idx,"]",sep="")] ,
                         xlab=bquote(beta*2[.(x2idx)]) ,
                         main=paste("x2:",x2names[x2idx]) )
  }
  for ( x2idx in 1:dataList$Nx2Lvl ) {
    for ( x1idx in 1:dataList$Nx1Lvl ) {
      postInfo = plotPost( mcmcChain[,paste("b1b2[",x1idx,",",x2idx,"]",sep="")] , 
                           xlab=bquote(beta*12[.(x1idx)*","*.(x2idx)]) ,
                           main=paste("x1:",x1names[x1idx],", x2:",x2names[x2idx])  )
    }
  }
  if ( savePlotImages ) {
    saveGraph( file=paste(fileNameRoot,"b",sep="") , type=savePlotType )
  }
}

# Display contrast analyses
nContrasts = length( x1contrastList )
if ( nContrasts > 0 ) {
   nPlotPerRow = 5
   nPlotRow = ceiling(nContrasts/nPlotPerRow)
   nPlotCol = ceiling(nContrasts/nPlotRow)
   openGraph(width=3.75*nPlotCol,height=2.5*nPlotRow)
   layout( matrix(1:(nPlotRow*nPlotCol),nrow=nPlotRow,ncol=nPlotCol,byrow=T) )
   par( mar=c(4,0.5,2.5,0.5) , mgp=c(2,0.7,0) )
   for ( cIdx in 1:nContrasts ) {
       contrast = matrix( x1contrastList[[cIdx]],nrow=1) # make it a row matrix
       incIdx = contrast!=0
       postInfo = plotPost( contrast %*%
         t(mcmcChain[,paste("b1[",1:dataList$Nx1Lvl,"]",sep="")]) , 
         xlab=paste( round(contrast[incIdx],2) , x1names[incIdx] ,
                     c(rep("+",sum(incIdx)-1),"") , collapse=" " ) ,
         cex.lab = 1.0 , compVal=0 ,
         main=paste( "X1 Contrast:", names(x1contrastList)[cIdx] )  )
   }
   if ( savePlotImages ) {
     saveGraph( file=paste(fileNameRoot,"x1contrasts",sep="") , type=savePlotType )
   }
}
#
nContrasts = length( x2contrastList )
if ( nContrasts > 0 ) {
  nPlotPerRow = 5
  nPlotRow = ceiling(nContrasts/nPlotPerRow)
  nPlotCol = ceiling(nContrasts/nPlotRow)
  openGraph(width=3.75*nPlotCol,height=2.5*nPlotRow)
  layout( matrix(1:(nPlotRow*nPlotCol),nrow=nPlotRow,ncol=nPlotCol,byrow=T) )
  par( mar=c(4,0.5,2.5,0.5) , mgp=c(2,0.7,0) )
  for ( cIdx in 1:nContrasts ) {
    contrast = matrix( x2contrastList[[cIdx]],nrow=1) # make it a row matrix
    incIdx = contrast!=0
    postInfo = plotPost( contrast %*%
      t(mcmcChain[,paste("b2[",1:dataList$Nx2Lvl,"]",sep="")]) , 
                         xlab=paste( round(contrast[incIdx],2) , x2names[incIdx] ,
                                     c(rep("+",sum(incIdx)-1),"") , collapse=" " ) ,
                         cex.lab = 1.0 , compVal=0 ,
                         main=paste( "X2 Contrast:", names(x2contrastList)[cIdx] )  )
  }
  if ( savePlotImages ) {
    saveGraph( file=paste(fileNameRoot,"x2contrasts",sep="") , type=savePlotType )
  }
}
#
nContrasts = length( x1x2contrastList )
if ( nContrasts > 0 ) {
   nPlotPerRow = 5
   nPlotRow = ceiling(nContrasts/nPlotPerRow)
   nPlotCol = ceiling(nContrasts/nPlotRow)
   openGraph(width=3.75*nPlotCol,height=2.5*nPlotRow)
   layout( matrix(1:(nPlotRow*nPlotCol),nrow=nPlotRow,ncol=nPlotCol,byrow=T) )
   par( mar=c(4,0.5,2.5,0.5) , mgp=c(2,0.7,0) )
   b1b2Sample = array(0, dim=c( dataList$Nx1Lvl , dataList$Nx2Lvl , chainLength ) )
   for ( x1idx in 1:dataList$Nx1Lvl ) {
     for ( x2idx in 1:dataList$Nx2Lvl ) {
       b1b2Sample[x1idx,x2idx,] = mcmcChain[, paste( "b1b2[",x1idx,",",x2idx,"]",
                                                     sep="" ) ]
     }
   }
   for ( cIdx in 1:nContrasts ) {
       contrast = x1x2contrastList[[cIdx]]
       contrastArr = array( rep(contrast,chainLength) ,
                            dim=c(NROW(contrast),NCOL(contrast),chainLength) )
       contrastLab = ""
       for ( x1idx in 1:Nx1Lvl ) {
         for ( x2idx in 1:Nx2Lvl ) {
           if ( contrast[x1idx,x2idx] != 0 ) {
             contrastLab = paste( contrastLab , "+" ,
                                  signif(contrast[x1idx,x2idx],2) ,
                                  x1names[x1idx] , x2names[x2idx] )
           }
         }
       }
       postInfo = plotPost( apply( contrastArr * b1b2Sample , 3 , sum ) ,
                compVal=0 ,  xlab=contrastLab , cex.lab = 0.75 ,
                main=paste( names(x1x2contrastList)[cIdx] ) )
   }
   if ( savePlotImages ) {
     saveGraph( file=paste(fileNameRoot,"x1x2Contrasts",sep="") , type=savePlotType )
   }
}

#==============================================================================
# Do NHST ANOVA:

theData = data.frame( y=y , x1=factor(x1,labels=x1names) ,
                            x2=factor(x2,labels=x2names) )
openGraph(width=7,height=7)
interaction.plot( theData$x1 , theData$x2 , theData$y , type="b" )
#saveGraph( file=paste(fileNameRoot,"DataPlot",sep="") , type="eps" )
aovresult = aov( y ~ x1 * x2 , data = theData )
cat("\n------------------------------------------------------------------\n\n")
print( summary( aovresult ) )
cat("\n------------------------------------------------------------------\n\n")
print( model.tables( aovresult , type = "effects", se = TRUE ) , digits=3 )
cat("\n------------------------------------------------------------------\n\n")

#==============================================================================
