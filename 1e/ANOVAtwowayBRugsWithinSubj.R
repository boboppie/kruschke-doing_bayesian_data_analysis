graphics.off()
rm(list=ls(all=TRUE))
fnroot = "ANOVAtwowayBrugsWithinSubj"
library(BRugs)         # Kruschke, J. K. (2010). Doing Bayesian data analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
# BUGS model specification begins here...
model {
  for ( i in 1:Ntotal ) {
    y[i] ~ dnorm( mu[i] , tau )
    mu[i] <- a0 + a1[x1[i]] + a2[x2[i]] + a1a2[x1[i],x2[i]] + aS[S[i]]
  }
  #
  tau <- pow( sigma , -2 )
  sigma ~ dunif(0,10) # y values are assumed to be standardized
  #
  a0 ~ dnorm(0,0.001) # y values are assumed to be standardized
  #
  for ( j1 in 1:Nx1Lvl ) { a1[j1] ~ dnorm( 0.0 , a1tau ) }
  a1tau <- 1 / pow( a1SD , 2 )
  a1SD <- abs( a1SDunabs ) + .1
  a1SDunabs ~ dt( 0 , 0.001 , 2 )
  #
  for ( j2 in 1:Nx2Lvl ) { a2[j2] ~ dnorm( 0.0 , a2tau ) }
  a2tau <- 1 / pow( a2SD , 2 )
  a2SD <- abs( a2SDunabs ) + .1
  a2SDunabs ~ dt( 0 , 0.001 , 2 )
  #
  for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
    a1a2[j1,j2] ~ dnorm( 0.0 , a1a2tau )
  } }
  a1a2tau <- 1 / pow( a1a2SD , 2 )
  a1a2SD <- abs( a1a2SDunabs ) + .1
  a1a2SDunabs ~ dt( 0 , 0.001 , 2 )
  #
  for ( jS in 1:NSLvl ) { aS[jS] ~ dnorm( 0.0 , aStau ) }
  aStau <- 1 / pow( aSSD , 2 )
  aSSD <- abs( aSSDunabs ) + .1
  aSSDunabs ~ dt( 0 , 0.001 , 2 )
}
# ... end BUGS model specification
" # close quote for modelstring
# Write model to a file, and send to BUGS:
writeLines(modelstring,con="model.txt")
modelCheck( "model.txt" )

#------------------------------------------------------------------------------
# THE DATA.
# Specify data source:
dataSource = c( "Ex19.3" )[1]

# Load the data:
if ( dataSource == "Ex19.3" ) {
  fnroot = paste( fnroot , dataSource , sep="" )
  y = c( 101,102,103,105,104, 104,105,107,106,108, 105,107,106,108,109, 109,108,110,111,112 )
  x1 = c( 1,1,1,1,1, 1,1,1,1,1, 2,2,2,2,2, 2,2,2,2,2 )
  x2 = c( 1,1,1,1,1, 2,2,2,2,2, 1,1,1,1,1, 2,2,2,2,2 )
  S = c( 1,2,3,4,5, 1,2,3,4,5, 1,2,3,4,5, 1,2,3,4,5 )
  x1names = c("x1.1","x1.2")
  x2names = c("x2.1","x2.2")
  Snames = c("S1","S2","S3","S4","S5")
  Ntotal = length(y)
  Nx1Lvl = length(unique(x1))
  Nx2Lvl = length(unique(x2))
  NSLvl = length(unique(S))
  x1contrastList = list( X1.2vX1.1 = c( -1 , 1 ) )
  x2contrastList = list( X2.2vX2.1 = c( -1 , 1 ) )
  x1x2contrastList = NULL # list( matrix( 1:(Nx1Lvl*Nx2Lvl) , nrow=Nx1Lvl ) )
}

# Specify the data in a form that is compatible with BRugs model, as a list:
ySDorig = sd(y)
yMorig = mean(y)
z = ( y - yMorig ) / ySDorig
datalist = list(
  y = z ,
  x1 = x1 ,
  x2 = x2 ,
  S = S ,
  Ntotal = Ntotal ,
  Nx1Lvl = Nx1Lvl ,
  Nx2Lvl = Nx2Lvl ,
  NSLvl = NSLvl
)
# Get the data into BRugs:
modelData( bugsData( datalist ) )

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

# Autocorrelation within chains is large, so use several chains to reduce
# degree of thinning. But we still have to burn-in all the chains, which takes
# more time with more chains.
nchain = 10
modelCompile( numChains = nchain )

if ( F ) {
   modelGenInits() # often won't work for diffuse prior
} else {
  #  initialization based on data
  theData = data.frame( y=datalist$y , x1=factor(x1,labels=x1names) ,
                        x2=factor(x2,labels=x2names) , S=factor(S,labels=Snames) )
  a0 = mean( theData$y )
  a1 = aggregate( theData$y , list( theData$x1 ) , mean )[,2] - a0
  a2 = aggregate( theData$y , list( theData$x2 ) , mean )[,2] - a0
  aS = aggregate( theData$y , list( theData$S ) , mean )[,2] - a0
  linpred = as.vector( outer( a1 , a2 , "+" ) + a0 )
  a1a2 = aggregate( theData$y, list(theData$x1,theData$x2), mean)[,3] - linpred
  genInitList <- function() {
    return(
        list(
            a0 = a0 ,
            a1 = a1 ,
            a2 = a2 ,
            aS = aS ,
            a1a2 = matrix( a1a2 , nrow=Nx1Lvl , ncol=Nx2Lvl ) ,
            sigma = sd(theData$y)/2 , # lazy
            a1SDunabs = sd(a1) ,
            a2SDunabs = sd(a2) ,
            a1a2SDunabs = sd(a1a2) ,
            aSSDunabs = sd(aS)
        )
    )
  }
  for ( chainIdx in 1 : nchain ) {
    modelInits( bugsInits( genInitList ) )
  }
}

#------------------------------------------------------------------------------
# RUN THE CHAINS

# burn in
BurnInSteps = 10000
modelUpdate( BurnInSteps )
# actual samples
samplesSet( c( "a0" ,  "a1" ,  "a2" ,  "a1a2" , "aS" ,
               "sigma" , "a1SD" , "a2SD" , "a1a2SD" , "aSSD" ) )
stepsPerChain = ceiling(2000/nchain)
thinStep = 100 # 500
modelUpdate( stepsPerChain , thin=thinStep )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

source("plotChains.R")
source("plotPost.R")

checkConvergence = F
if ( checkConvergence ) {
   sumInfo = plotChains( "a0" , saveplots=F , filenameroot=fnroot )
   sumInfo = plotChains( "a1" , saveplots=F , filenameroot=fnroot )
   sumInfo = plotChains( "a2" , saveplots=F , filenameroot=fnroot )
   sumInfo = plotChains( "a1a2" , saveplots=F , filenameroot=fnroot )
   sumInfo = plotChains( "aS" , saveplots=F , filenameroot=fnroot )
   readline("Press any key to clear graphics and continue")
   graphics.off()
   sumInfo = plotChains( "sigma" , saveplots=F , filenameroot=fnroot )
   sumInfo = plotChains( "a1SD" , saveplots=F , filenameroot=fnroot )
   sumInfo = plotChains( "a2SD" , saveplots=F , filenameroot=fnroot )
   sumInfo = plotChains( "a1a2SD" , saveplots=F , filenameroot=fnroot )
   sumInfo = plotChains( "aSSD" , saveplots=F , filenameroot=fnroot )
   readline("Press any key to clear graphics and continue")
   graphics.off()
}

# Extract and plot the SDs:
sigmaSample = samplesSample("sigma")
a1SDSample = samplesSample("a1SD")
a2SDSample = samplesSample("a2SD")
a1a2SDSample = samplesSample("a1a2SD")
aSSDSample = samplesSample("aSSD")
windows()
layout( matrix(1:6,nrow=2) )
par( mar=c(3,1,2.5,0) , mgp=c(2,0.7,0) )
plotPost( sigmaSample , xlab="sigma" , main="Cell SD" , breaks=30 )
plotPost( a1SDSample , xlab="a1SD" , main="a1 SD" , breaks=30 )
plotPost( a2SDSample , xlab="a2SD" , main="a2 SD" , breaks=30 )
plotPost( a1a2SDSample , xlab="a1a2SD" , main="Interaction SD" , breaks=30 )
plotPost( aSSDSample , xlab="aSSD" , main="aS SD" , breaks=30 )
dev.copy2eps(file=paste(fnroot,"SD.eps",sep=""))

# Extract a values:
a0Sample = samplesSample( "a0" )
chainLength = length(a0Sample)
a1Sample = array( 0 , dim=c( datalist$Nx1Lvl , chainLength ) )
for ( x1idx in 1:datalist$Nx1Lvl ) {
   a1Sample[x1idx,] = samplesSample( paste("a1[",x1idx,"]",sep="") )
}
a2Sample = array( 0 , dim=c( datalist$Nx2Lvl , chainLength ) )
for ( x2idx in 1:datalist$Nx2Lvl ) {
   a2Sample[x2idx,] = samplesSample( paste("a2[",x2idx,"]",sep="") )
}
a1a2Sample = array(0, dim=c( datalist$Nx1Lvl , datalist$Nx2Lvl , chainLength ) )
for ( x1idx in 1:datalist$Nx1Lvl ) {
  for ( x2idx in 1:datalist$Nx2Lvl ) {
    a1a2Sample[x1idx,x2idx,] = samplesSample( paste( "a1a2[",x1idx,",",x2idx,"]",
                                                     sep="" ) )
  }
}
aSSample = array( 0 , dim=c( datalist$NSLvl , chainLength ) )
for ( Sidx in 1:datalist$NSLvl ) {
   aSSample[Sidx,] = samplesSample( paste("aS[",Sidx,"]",sep="") )
}

# Convert the a values to zero-centered b values.
# m12Sample is predicted cell means at every step in MCMC chain:
m12Sample = array( 0, dim=c( datalist$Nx1Lvl , datalist$Nx2Lvl ,
                             datalist$NSLvl , chainLength ) )
for ( stepIdx in 1:chainLength ) {
  for ( a1idx in 1:Nx1Lvl ) {
    for ( a2idx in 1:Nx2Lvl ) {
      for ( aSidx in 1:NSLvl ) {
        m12Sample[ a1idx , a2idx , aSidx , stepIdx ] = (
           a0Sample[stepIdx]
           + a1Sample[a1idx,stepIdx]
           + a2Sample[a2idx,stepIdx]
           + a1a2Sample[a1idx,a2idx,stepIdx]
           + aSSample[aSidx,stepIdx] )
      }
    }
  }
}

# b0Sample is mean of the cell means at every step in chain:
b0Sample = apply( m12Sample , 4 , mean )
# b1Sample is deflections of factor 1 marginal means from b0Sample:
b1Sample = ( apply( m12Sample , c(1,4) , mean )
           - matrix(rep( b0Sample ,Nx1Lvl),nrow=Nx1Lvl,byrow=T) )
# b2Sample is deflections of factor 2 marginal means from b0Sample:
b2Sample = ( apply( m12Sample , c(2,4) , mean )
           - matrix(rep( b0Sample ,Nx2Lvl),nrow=Nx2Lvl,byrow=T) )
# bSSample is deflections of factor S marginal means from b0Sample:
bSSample = ( apply( m12Sample , c(3,4) , mean )
           - matrix(rep( b0Sample ,NSLvl),nrow=NSLvl,byrow=T) )
# linpredSample is linear combination of the marginal effects:
linpredSample = 0*m12Sample
for ( stepIdx in 1:chainLength ) {
  for ( a1idx in 1:Nx1Lvl ) {
    for ( a2idx in 1:Nx2Lvl ) {
      for ( aSidx in 1:NSLvl ) {
        linpredSample[a1idx,a2idx,aSidx,stepIdx ] = (
           b0Sample[stepIdx]
           + b1Sample[a1idx,stepIdx]
           + b2Sample[a2idx,stepIdx]
           + bSSample[aSidx,stepIdx] )
      }
    }
  }
}
# b1b2Sample is the interaction deflection, i.e., the difference
# between the cell means and the linear combination:
b1b2Sample = apply( m12Sample - linpredSample , c(1,2,4) , mean )

# Convert from standardized b values to original scale b values:
b0Sample = b0Sample * ySDorig + yMorig
b1Sample = b1Sample * ySDorig
b2Sample = b2Sample * ySDorig
bSSample = bSSample * ySDorig
b1b2Sample = b1b2Sample * ySDorig

# Plot b values:
windows((datalist$Nx1Lvl+1)*2.75,(datalist$Nx2Lvl+1)*2.0)
layoutMat = matrix( 0 , nrow=(datalist$Nx2Lvl+1) , ncol=(datalist$Nx1Lvl+1) )
layoutMat[1,1] = 1
layoutMat[1,2:(datalist$Nx1Lvl+1)] = 1:datalist$Nx1Lvl + 1
layoutMat[2:(datalist$Nx2Lvl+1),1] = 1:datalist$Nx2Lvl + (datalist$Nx1Lvl + 1)
layoutMat[2:(datalist$Nx2Lvl+1),2:(datalist$Nx1Lvl+1)] = matrix(
    1:(datalist$Nx1Lvl*datalist$Nx2Lvl) + (datalist$Nx2Lvl+datalist$Nx1Lvl+1) ,
    ncol=datalist$Nx1Lvl , byrow=T )
layout( layoutMat )
par( mar=c(4,0.5,2.5,0.5) , mgp=c(2,0.7,0) )
histinfo = plotPost( b0Sample , xlab=expression(beta * 0) , main="Baseline" ,
                     breaks=30  )
for ( x1idx in 1:datalist$Nx1Lvl ) {
  histinfo = plotPost( b1Sample[x1idx,] , xlab=bquote(beta*1[.(x1idx)]) ,
                       main=paste("x1:",x1names[x1idx]) , breaks=30 )
}
for ( x2idx in 1:datalist$Nx2Lvl ) {
  histinfo = plotPost( b2Sample[x2idx,] , xlab=bquote(beta*2[.(x2idx)]) ,
                       main=paste("x2:",x2names[x2idx]) , breaks=30 )
}
for ( x2idx in 1:datalist$Nx2Lvl ) {
  for ( x1idx in 1:datalist$Nx1Lvl ) {
    histinfo = plotPost( b1b2Sample[x1idx,x2idx,] , breaks=30 ,
              xlab=bquote(beta*12[.(x1idx)*","*.(x2idx)]) ,
              main=paste("x1:",x1names[x1idx],", x2:",x2names[x2idx])  )
  }
}
dev.copy2eps(file=paste(fnroot,"b.eps",sep=""))

# Display contrast analyses
nContrasts = length( x1contrastList )
if ( nContrasts > 0 ) {
   nPlotPerRow = 5
   nPlotRow = ceiling(nContrasts/nPlotPerRow)
   nPlotCol = ceiling(nContrasts/nPlotRow)
   windows(3.75*nPlotCol,2.5*nPlotRow)
   layout( matrix(1:(nPlotRow*nPlotCol),nrow=nPlotRow,ncol=nPlotCol,byrow=T) )
   par( mar=c(4,0.5,2.5,0.5) , mgp=c(2,0.7,0) )
   for ( cIdx in 1:nContrasts ) {
       contrast = matrix( x1contrastList[[cIdx]],nrow=1) # make it a row matrix
       incIdx = contrast!=0
       histInfo = plotPost( contrast %*% b1Sample , compVal=0 , breaks=30 ,
                xlab=paste( round(contrast[incIdx],2) , x1names[incIdx] ,
                            c(rep("+",sum(incIdx)-1),"") , collapse=" " ) ,
                cex.lab = 1.0 ,
                main=paste( "X1 Contrast:", names(x1contrastList)[cIdx] ) )
   }
   dev.copy2eps(file=paste(fnroot,"x1Contrasts.eps",sep=""))
}
#
nContrasts = length( x2contrastList )
if ( nContrasts > 0 ) {
   nPlotPerRow = 5
   nPlotRow = ceiling(nContrasts/nPlotPerRow)
   nPlotCol = ceiling(nContrasts/nPlotRow)
   windows(3.75*nPlotCol,2.5*nPlotRow)
   layout( matrix(1:(nPlotRow*nPlotCol),nrow=nPlotRow,ncol=nPlotCol,byrow=T) )
   par( mar=c(4,0.5,2.5,0.5) , mgp=c(2,0.7,0) )
   for ( cIdx in 1:nContrasts ) {
       contrast = matrix( x2contrastList[[cIdx]],nrow=1) # make it a row matrix
       incIdx = contrast!=0
       histInfo = plotPost( contrast %*% b2Sample , compVal=0 , breaks=30 ,
                xlab=paste( round(contrast[incIdx],2) , x2names[incIdx] ,
                            c(rep("+",sum(incIdx)-1),"") , collapse=" " ) ,
                cex.lab = 1.0 ,
                main=paste( "X2 Contrast:", names(x2contrastList)[cIdx] ) )
   }
   dev.copy2eps(file=paste(fnroot,"x2Contrasts.eps",sep=""))
}
#
nContrasts = length( x1x2contrastList )
if ( nContrasts > 0 ) {
   nPlotPerRow = 5
   nPlotRow = ceiling(nContrasts/nPlotPerRow)
   nPlotCol = ceiling(nContrasts/nPlotRow)
   windows(3.75*nPlotCol,2.5*nPlotRow)
   layout( matrix(1:(nPlotRow*nPlotCol),nrow=nPlotRow,ncol=nPlotCol,byrow=T) )
   par( mar=c(4,0.5,2.5,0.5) , mgp=c(2,0.7,0) )
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
       histInfo = plotPost( apply( contrastArr * b1b2Sample , 3 , sum ) ,
                compVal=0 , breaks=30 , xlab=contrastLab , cex.lab = 0.75 ,
                main=paste( names(x1x2contrastList)[cIdx] ) )
   }
   dev.copy2eps(file=paste(fnroot,"x1x2Contrasts.eps",sep=""))
}

#==============================================================================
# Do NHST ANOVA:

theData = data.frame( y=y , x1=factor(x1,labels=x1names) ,
                            x2=factor(x2,labels=x2names) )
windows()
interaction.plot( theData$x1 , theData$x2 , theData$y , type="b" )
dev.copy2eps(file=paste(fnroot,"DataPlot.eps",sep=""))
aovresult = aov( y ~ x1 * x2 , data = theData )
cat("\n------------------------------------------------------------------\n\n")
print( summary( aovresult ) )
cat("\n------------------------------------------------------------------\n\n")
print( model.tables( aovresult , type = "effects", se = TRUE ) , digits=3 )
cat("\n------------------------------------------------------------------\n\n")

#==============================================================================
