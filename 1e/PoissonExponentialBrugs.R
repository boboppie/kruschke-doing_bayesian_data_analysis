graphics.off()
rm(list=ls(all=TRUE))
fnroot = "PoissonExponentialBrugs"
library(BRugs)         # Kruschke, J. K. (2010). Doing Bayesian data analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
# BUGS model specification begins here...
model {
  for ( i in 1:Ncells ) {
    y[i] ~ dpois( lambda[i] )
    lambda[i] <- exp( a0 + a1[x1[i]] + a2[x2[i]] + a1a2[x1[i],x2[i]] )
  }
  #
  a0 ~ dnorm(10,1.0E-6)
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
}
# ... end BUGS model specification
" # close quote for modelstring
# Write model to a file, and send to BUGS:
writeLines(modelstring,con="model.txt")
modelCheck( "model.txt" )

#------------------------------------------------------------------------------
# THE DATA.
# Specify data source:
dataSource = c( "HairEye" , "CrimeDrink" , "Toy" )[1]

# Load the data:
if ( dataSource == "HairEye" ) {
  fnroot = paste( fnroot , dataSource , sep="" )
  dataFrame = data.frame( # from Snee (1974)
    Freq = c(68,119,26,7,20,84,17,94,15,54,14,10,5,29,14,16) ,
    Eye  = c("Brown","Brown","Brown","Brown","Blue","Blue","Blue","Blue","Hazel","Hazel","Hazel","Hazel","Green","Green","Green","Green"),
    Hair = c("Black","Brunette","Red","Blond","Black","Brunette","Red","Blond","Black","Brunette","Red","Blond","Black","Brunette","Red","Blond") )
  y = as.numeric(dataFrame$Freq)
  x1 = as.numeric(dataFrame$Eye)
  x1names = levels(dataFrame$Eye)
  x2 = as.numeric(dataFrame$Hair)
  x2names = levels(dataFrame$Hair)
  Ncells = length(y)
  Nx1Lvl = length(unique(x1))
  Nx2Lvl = length(unique(x2))
  x1contrastList = list( GREENvHAZEL = c(0,0,1,-1) )
  x2contrastList = list( BLONDvRED = c(0,1,0,-1) )
  x1x2contrastList = list( BLUE.BROWNxBLACK.BLOND
                           = outer(c(-1,1,0,0),c(-1,1,0,0)) )
}

if ( dataSource == "CrimeDrink" ) {
  fnroot = paste( fnroot , dataSource , sep="" )
  dataFrame = data.frame( # from Kendall (1943) via Snee (1974)
    Freq = c(50,88,155,379,18,63,43,62,110,300,14,144) ,
    Crime  = c("Arson","Rape","Violence","Theft","Coining","Fraud","Arson","Rape","Violence","Theft","Coining","Fraud"),
    Drink = c("Drinker","Drinker","Drinker","Drinker","Drinker","Drinker","Nondrink","Nondrink","Nondrink","Nondrink","Nondrink","Nondrink") )
  y = as.numeric(dataFrame$Freq)
  x1 = as.numeric(dataFrame$Crime)
  x1names = levels(dataFrame$Crime)
  x2 = as.numeric(dataFrame$Drink)
  x2names = levels(dataFrame$Drink)
  Ncells = length(y)
  Nx1Lvl = length(unique(x1))
  Nx2Lvl = length(unique(x2))
  x1contrastList = list( FRAUDvOTHER = c(-1/5,-1/5,1,-1/5,-1/5,-1/5) ,
                         FRAUDvRAPE = c(0,0,1,-1,0,0) )
  x2contrastList = list( DRINKERvNON = c(1,-1) )
  x1x2contrastList = list( FRAUD.OTHERxDRINKER.NON
                           = outer(c(-1/5,-1/5,1,-1/5,-1/5,-1/5),c(-1,1)) ,
                           FRAUD.RAPExDRINKER.NON
                           = outer(c(0,0,1,-1,0,0),c(-1,1)) )
}

if ( dataSource == "Toy" ) {
  dataMultiplier = 2 # Try 2 (chi-sq warns) , 6 (p>.05) , 7 (p<.05) , 10
  fnroot = paste( fnroot , dataSource , dataMultiplier , sep="" )
  dataFrame = data.frame(
    Freq = c( 2,2,1,1, 2,2,1,1, 1,1,2,2, 1,1,2,2 ) * dataMultiplier ,
    Col  = c("C1","C2","C3","C4", "C1","C2","C3","C4", "C1","C2","C3","C4", "C1","C2","C3","C4"),
    Row = c("R1","R1","R1","R1", "R2","R2","R2","R2", "R3","R3","R3","R3", "R4","R4","R4","R4") )
  y = as.numeric(dataFrame$Freq)
  x1 = as.numeric(dataFrame$Col)
  x1names = levels(dataFrame$Col)
  x2 = as.numeric(dataFrame$Row)
  x2names = levels(dataFrame$Row)
  Ncells = length(y)
  Nx1Lvl = length(unique(x1))
  Nx2Lvl = length(unique(x2))
  x1contrastList = NULL
  x2contrastList = NULL
  x1x2contrastList = list( R12.R34xC12.C34 = outer(c(1,1,-1,-1)/2,c(1,1,-1,-1)/2) )
}

# Specify the data in a form that is compatible with BRugs model, as a list:
datalist = list(
  y = y ,
  x1 = x1 ,
  x2 = x2 ,
  Ncells = Ncells ,
  Nx1Lvl = Nx1Lvl ,
  Nx2Lvl = Nx2Lvl
)
# Get the data into BRugs:
modelData( bugsData( datalist ) )

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

nchain = 5
modelCompile( numChains = nchain )

if ( F ) {
   modelGenInits() # often won't work for diffuse prior
} else {
  #  initialization based on data
  theData = data.frame( y=log(y) , x1=factor(x1,labels=x1names) ,
                              x2=factor(x2,labels=x2names) )
  a0 = mean( theData$y )
  a1 = aggregate( theData$y , list( theData$x1 ) , mean )[,2] - a0
  a2 = aggregate( theData$y , list( theData$x2 ) , mean )[,2] - a0
  linpred = as.vector( outer( a1 , a2 , "+" ) + a0 )
  a1a2 = aggregate( theData$y, list(theData$x1,theData$x2), mean)[,3] - linpred
  genInitList <- function() {
    return(
        list(
            a0 = a0 ,
            a1 = a1 ,
            a2 = a2 ,
            a1a2 = matrix( a1a2 , nrow=Nx1Lvl , ncol=Nx2Lvl ) ,
            a1SDunabs = sd(a1) ,
            a2SDunabs = sd(a2) ,
            a1a2SDunabs = sd(a1a2)
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
BurnInSteps = 1000
modelUpdate( BurnInSteps )
# actual samples
samplesSet( c( "a0" ,  "a1" ,  "a2" ,  "a1a2" , "a1SD" , "a2SD" , "a1a2SD" ) )
stepsPerChain = ceiling(5000/nchain)
thinStep = 500
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
   readline("Press any key to clear graphics and continue")
   graphics.off()
   sumInfo = plotChains( "a1SD" , saveplots=F , filenameroot=fnroot )
   sumInfo = plotChains( "a2SD" , saveplots=F , filenameroot=fnroot )
   sumInfo = plotChains( "a1a2SD" , saveplots=F , filenameroot=fnroot )
   readline("Press any key to clear graphics and continue")
   graphics.off()
}

# Extract and plot the SDs:
a1SDSample = samplesSample("a1SD")
a2SDSample = samplesSample("a2SD")
a1a2SDSample = samplesSample("a1a2SD")
windows(10,3)
layout( matrix(1:3,nrow=1) )
par( mar=c(3,1,2.5,0) , mgp=c(2,0.7,0) )
histInfo = plotPost( a1SDSample , xlab="a1SD" , main="a1 SD" , breaks=30 )
histInfo = plotPost( a2SDSample , xlab="a2SD" , main="a2 SD" , breaks=30 )
histInfo = plotPost( a1a2SDSample , xlab="a1a2SD" , main="Interaction SD" ,
                     breaks=30 )
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

# Convert to zero-centered b values:
m12Sample = array( 0, dim=c( datalist$Nx1Lvl , datalist$Nx2Lvl , chainLength ) )
for ( stepIdx in 1:chainLength ) {
    m12Sample[,,stepIdx ] = ( a0Sample[stepIdx]
                            + outer( a1Sample[,stepIdx] ,
                                     a2Sample[,stepIdx] , "+" )
                            + a1a2Sample[,,stepIdx] )
}
b0Sample = apply( m12Sample , 3 , mean )
b1Sample = ( apply( m12Sample , c(1,3) , mean )
           - matrix(rep( b0Sample ,Nx1Lvl),nrow=Nx1Lvl,byrow=T) )
b2Sample = ( apply( m12Sample , c(2,3) , mean )
           - matrix(rep( b0Sample ,Nx2Lvl),nrow=Nx2Lvl,byrow=T) )
linpredSample = array(0,dim=c(datalist$Nx1Lvl,datalist$Nx2Lvl,chainLength))
for ( stepIdx in 1:chainLength ) {
    linpredSample[,,stepIdx ] = ( b0Sample[stepIdx]
                                + outer( b1Sample[,stepIdx] ,
                                         b2Sample[,stepIdx] , "+" ) )
}
b1b2Sample = m12Sample - linpredSample

# Plot b values:
windows((datalist$Nx1Lvl+1)*2.75,(datalist$Nx2Lvl+1)*2.25)
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
    hdiLim = HDIofMCMC(b1b2Sample[x1idx,x2idx,])
    if ( hdiLim[1]>0 | hdiLim[2]<0 ) { compVal=0 } else { compVal=NULL }
    histinfo = plotPost( b1b2Sample[x1idx,x2idx,] , breaks=30 , compVal=compVal ,
              xlab=bquote(beta*12[list(x1==.(x1idx),x2==.(x2idx))]) ,
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

# Compute credible cell probability at each step in the MCMC chain
lambda12Sample = 0 * b1b2Sample
for ( chainIdx in 1:chainLength ) {
    lambda12Sample[,,chainIdx] = exp(
        b0Sample[chainIdx]
        + outer( b1Sample[,chainIdx] , b2Sample[,chainIdx] , "+" )
        + b1b2Sample[,,chainIdx] )
}
cellp = 0 * lambda12Sample
for ( chainIdx in 1:chainLength ) {
    cellp[,,chainIdx] = ( lambda12Sample[,,chainIdx]
                          / sum( lambda12Sample[,,chainIdx] ) )
}
# Display credible cell probabilities
windows((datalist$Nx1Lvl)*2.75,(datalist$Nx2Lvl)*2.25)
layoutMat = matrix( 1:(datalist$Nx2Lvl*datalist$Nx1Lvl) ,
                    nrow=(datalist$Nx2Lvl) , ncol=(datalist$Nx1Lvl) , byrow=T )
layout( layoutMat )
par( mar=c(4,1.5,2.5,0.5) , mgp=c(2,0.7,0) )
maxp = max( cellp )
for ( x2idx in 1:datalist$Nx2Lvl ) {
  for ( x1idx in 1:datalist$Nx1Lvl ) {
    histinfo = plotPost( cellp[x1idx,x2idx,] ,
               breaks=seq(0,maxp,length=51) , xlim=c(0,maxp) ,
               xlab=bquote(probability[list(x1==.(x1idx),x2==.(x2idx))]) ,
               main=paste("x1:",x1names[x1idx],", x2:",x2names[x2idx]) ,
               HDItextPlace=0.95 )
  }
}
dev.copy2eps(file=paste(fnroot,"CellP.eps",sep=""))


#==============================================================================
# Conduct NHST Pearson chi-square test of independence.

# Convert dataFrame to frequency table:
obsFreq = matrix( 0 , nrow=Nx1Lvl , ncol=Nx2Lvl )
for ( x1idx in 1:Nx1Lvl ) {
  for ( x2idx in 1:Nx2Lvl ) {
    obsFreq[x1idx,x2idx] = y[ dataFrame[,2]==x1names[x1idx]
                              & dataFrame[,3]==x2names[x2idx] ]
  }
}
obsFreq = t(obsFreq) # merely to match orientation of histogram display
chisqtest = chisq.test( obsFreq )
print( "obs :" )
print( chisqtest$observed )
print( "( obs - exp )^2 / exp :" )
print( ( chisqtest$observed - chisqtest$expected )^2 / chisqtest$expected )
print( chisqtest )

#==============================================================================
