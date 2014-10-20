graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="PoissonExponentialJagsSTZ" # for constructing output filenames
source("openGraphSaveGraph.R")
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
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
  a1SD ~ dgamma(1.01005,0.01005) # mode=1,sd=100
  #
  for ( j2 in 1:Nx2Lvl ) { a2[j2] ~ dnorm( 0.0 , a2tau ) }
  a2tau <- 1 / pow( a2SD , 2 )
  a2SD ~ dgamma(1.01005,0.01005) # mode=1,sd=100
  #
  for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
    a1a2[j1,j2] ~ dnorm( 0.0 , a1a2tau )
  } }
  a1a2tau <- 1 / pow( a1a2SD , 2 )
  a1a2SD ~ dgamma(1.01005,0.01005) # mode=1,sd=100
  # Convert a0,a1[],a2[],a1a2[,] to sum-to-zero b0,b1[],b2[],b1b2[,] :
  for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
    m[j1,j2] <- a0 + a1[j1] + a2[j2] + a1a2[j1,j2]  
  } }
  b0 <- mean( m[1:Nx1Lvl,1:Nx2Lvl] )
  for ( j1 in 1:Nx1Lvl ) { b1[j1] <- mean( m[j1,1:Nx2Lvl] ) - b0 }
  for ( j2 in 1:Nx2Lvl ) { b2[j2] <- mean( m[1:Nx1Lvl,j2] ) - b0 }
  for ( j1 in 1:Nx1Lvl ) { for ( j2 in 1:Nx2Lvl ) {
    b1b2[j1,j2] <- m[j1,j2] - ( b0 + b1[j1] + b2[j2] )  
  } }
}
" # close quote for modelstring
# Write model to a file, and send to JAGS:
writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.
# Specify data source:
dataSource = c( "HairEye" , "CrimeDrink" , "Toy" )[1]

# Load the data:
if ( dataSource == "HairEye" ) {
  fileNameRoot = paste( fileNameRoot , dataSource , sep="" )
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
  normalize = function( v ){ return( v / sum(v) ) }
  x1contrastList = list( GREENvHAZEL = (x1names=="Hazel")-(x1names=="Green") )
  x2contrastList = list( BLONDvRED = (x2names=="Blond")-(x2names=="Red") )
  x1x2contrastList = list( BLUE.BROWNxBLACK.BLOND
                           = outer((x1names=="Blue")-(x1names=="Brown"),
                                   (x2names=="Blond")-(x2names=="Black")) )
}

if ( dataSource == "CrimeDrink" ) {
  fileNameRoot = paste( fileNameRoot , dataSource , sep="" )
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
  normalize = function( v ){ return( v / sum(v) ) }
  x1contrastList = list( FRAUDvOTHER = (x1names=="Fraud")-normalize(x1names!="Fraud") ,
                         FRAUDvRAPE = (x1names=="Fraud")-(x1names=="Rape") )
  x2contrastList = list( DRINKERvNON = (x2names=="Drinker")-normalize(x2names=="Nondrink") )
  x1x2contrastList = list( FRAUD.OTHERxDRINKER.NON
                           = outer((x1names=="Fraud")-normalize(x1names!="Fraud"),
                                   (x2names=="Drinker")-normalize(x2names=="Nondrink")) ,
                           FRAUD.RAPExDRINKER.NON
                           = outer((x1names=="Fraud")-(x1names=="Rape"),
                                   (x2names=="Drinker")-normalize(x2names=="Nondrink")) )
}

if ( dataSource == "Toy" ) {
  dataMultiplier = 2 # Try 2 (chi-sq warns) , 6 (p>.05) , 7 (p<.05) , 10
  fileNameRoot = paste( fileNameRoot , dataSource , dataMultiplier , sep="" )
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
  normalize = function( v ){ return( v / sum(v) ) }
  x1contrastList = NULL
  x2contrastList = NULL
  x1x2contrastList = list( R12.R34xC12.C34 = 
                            outer( normalize(x1names=="C1"|x1names=="C2") -
                                   normalize(x1names=="C3"|x1names=="C4")  ,
                                   normalize(x2names=="R1"|x2names=="R2") -
                                   normalize(x2names=="R3"|x2names=="R4") ) )
}

dataList = list(
  y = y ,
  x1 = x1 ,
  x2 = x2 ,
  Ncells = Ncells ,
  Nx1Lvl = Nx1Lvl ,
  Nx2Lvl = Nx2Lvl
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

theData = data.frame( y=log(y) , x1=factor(x1,labels=x1names) ,
                            x2=factor(x2,labels=x2names) )
a0 = mean( theData$y )
a1 = aggregate( theData$y , list( theData$x1 ) , mean )[,2] - a0
a2 = aggregate( theData$y , list( theData$x2 ) , mean )[,2] - a0
linpred = as.vector( outer( a1 , a2 , "+" ) + a0 )
a1a2 = aggregate( theData$y, list(theData$x1,theData$x2), mean)[,3] - linpred
initsList = list( a0 = a0 , a1 = a1 , a2 = a2 ,
            a1a2 = matrix( a1a2 , nrow=Nx1Lvl , ncol=Nx2Lvl ) ,
            a1SD = max(sd(a1),0.1) , a2SD = max(sd(a2),0.1) , a1a2SD = max(sd(a1a2),0.1) )

#------------------------------------------------------------------------------
# RUN THE CHAINS

parameters = c( "a0" ,  "a1" ,  "a2" ,  "a1a2" , 
                "b0" ,  "b1" ,  "b2" ,  "b1b2" ,
                "a1SD" , "a2SD" , "a1a2SD" )
adaptSteps = 1000 
burnInSteps = 2000 
nChains = 5 
numSavedSteps=50000 
thinSteps=1
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
# Create, initialize, and adapt the model:
jagsModel = jags.model( "model.txt" , data=dataList , inits=initsList , 
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

checkConvergence = F
if ( checkConvergence ) {
  show( summary( codaSamples ) )
  openGraph()
  plot( codaSamples , ask=F )  
  openGraph()
  autocorr.plot( codaSamples , ask=F )
}

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )


# Extract the SDs:
a1SDSample = mcmcChain[,"a1SD"]
a2SDSample = mcmcChain[,"a2SD"]
a1a2SDSample = mcmcChain[,"a1a2SD"]

# Extract b values:
b0Sample = mcmcChain[, "b0" ]
chainLength = length(b0Sample)
b1Sample = array( 0 , dim=c( dataList$Nx1Lvl , chainLength ) )
for ( x1idx in 1:dataList$Nx1Lvl ) {
   b1Sample[x1idx,] = mcmcChain[, paste("b1[",x1idx,"]",sep="") ]
}
b2Sample = array( 0 , dim=c( dataList$Nx2Lvl , chainLength ) )
for ( x2idx in 1:dataList$Nx2Lvl ) {
   b2Sample[x2idx,] = mcmcChain[, paste("b2[",x2idx,"]",sep="") ]
}
b1b2Sample = array(0, dim=c( dataList$Nx1Lvl , dataList$Nx2Lvl , chainLength ) )
for ( x1idx in 1:dataList$Nx1Lvl ) {
  for ( x2idx in 1:dataList$Nx2Lvl ) {
    b1b2Sample[x1idx,x2idx,] = mcmcChain[, paste( "b1b2[",x1idx,",",x2idx,"]",
                                                     sep="" ) ]
  }
}

source("plotPost.R")

openGraph(10,3)
layout( matrix(1:3,nrow=1) )
par( mar=c(3,1,2.5,0) , mgp=c(2,0.7,0) )
histInfo = plotPost( a1SDSample , xlab="a1SD" , main="a1 SD" , showMode=T )
histInfo = plotPost( a2SDSample , xlab="a2SD" , main="a2 SD" , showMode=T )
histInfo = plotPost( a1a2SDSample , xlab="a1a2SD" , main="Interaction SD" ,
                     showMode=T )
#saveGraph(file=paste(fileNameRoot,"SD",sep=""),type="eps")

# Plot b values:
openGraph((dataList$Nx1Lvl+1)*2.75,(dataList$Nx2Lvl+1)*2.25)
layoutMat = matrix( 0 , nrow=(dataList$Nx2Lvl+1) , ncol=(dataList$Nx1Lvl+1) )
layoutMat[1,1] = 1
layoutMat[1,2:(dataList$Nx1Lvl+1)] = 1:dataList$Nx1Lvl + 1
layoutMat[2:(dataList$Nx2Lvl+1),1] = 1:dataList$Nx2Lvl + (dataList$Nx1Lvl + 1)
layoutMat[2:(dataList$Nx2Lvl+1),2:(dataList$Nx1Lvl+1)] = matrix(
    1:(dataList$Nx1Lvl*dataList$Nx2Lvl) + (dataList$Nx2Lvl+dataList$Nx1Lvl+1) ,
    ncol=dataList$Nx1Lvl , byrow=T )
layout( layoutMat )
par( mar=c(4,0.5,2.5,0.5) , mgp=c(2,0.7,0) )
histinfo = plotPost( b0Sample , xlab=expression(beta * 0) , main="Baseline" )
for ( x1idx in 1:dataList$Nx1Lvl ) {
  histinfo = plotPost( b1Sample[x1idx,] , xlab=bquote(beta*1[.(x1idx)]) ,
                       main=paste("x1:",x1names[x1idx]) )
}
for ( x2idx in 1:dataList$Nx2Lvl ) {
  histinfo = plotPost( b2Sample[x2idx,] , xlab=bquote(beta*2[.(x2idx)]) ,
                       main=paste("x2:",x2names[x2idx]) )
}
for ( x2idx in 1:dataList$Nx2Lvl ) {
  for ( x1idx in 1:dataList$Nx1Lvl ) {
    hdiLim = HDIofMCMC(b1b2Sample[x1idx,x2idx,])
    if ( hdiLim[1]>0 | hdiLim[2]<0 ) { compVal=0 } else { compVal=NULL }
    histinfo = plotPost( b1b2Sample[x1idx,x2idx,] , compVal=compVal ,
              xlab=bquote(beta*12[list(x1==.(x1idx),x2==.(x2idx))]) ,
              main=paste("x1:",x1names[x1idx],", x2:",x2names[x2idx])  )
  }
}
#saveGraph(file=paste(fileNameRoot,"b",sep=""),type="eps")

# Display contrast analyses
nContrasts = length( x1contrastList )
if ( nContrasts > 0 ) {
   nPlotPerRow = 5
   nPlotRow = ceiling(nContrasts/nPlotPerRow)
   nPlotCol = ceiling(nContrasts/nPlotRow)
   openGraph(3.75*nPlotCol,2.5*nPlotRow)
   layout( matrix(1:(nPlotRow*nPlotCol),nrow=nPlotRow,ncol=nPlotCol,byrow=T) )
   par( mar=c(4,0.5,2.5,0.5) , mgp=c(2,0.7,0) )
   for ( cIdx in 1:nContrasts ) {
       contrast = matrix( x1contrastList[[cIdx]],nrow=1) # make it a row matrix
       incIdx = contrast!=0
       histInfo = plotPost( contrast %*% b1Sample , compVal=0 ,
                xlab=paste( round(contrast[incIdx],2) , x1names[incIdx] ,
                            c(rep("+",sum(incIdx)-1),"") , collapse=" " ) ,
                cex.lab = 1.0 ,
                main=paste( "X1 Contrast:", names(x1contrastList)[cIdx] ) )
   }
   #saveGraph(file=paste(fileNameRoot,"x1Contrasts",sep=""),type="eps")
}
#
nContrasts = length( x2contrastList )
if ( nContrasts > 0 ) {
   nPlotPerRow = 5
   nPlotRow = ceiling(nContrasts/nPlotPerRow)
   nPlotCol = ceiling(nContrasts/nPlotRow)
   openGraph(3.75*nPlotCol,2.5*nPlotRow)
   layout( matrix(1:(nPlotRow*nPlotCol),nrow=nPlotRow,ncol=nPlotCol,byrow=T) )
   par( mar=c(4,0.5,2.5,0.5) , mgp=c(2,0.7,0) )
   for ( cIdx in 1:nContrasts ) {
       contrast = matrix( x2contrastList[[cIdx]],nrow=1) # make it a row matrix
       incIdx = contrast!=0
       histInfo = plotPost( contrast %*% b2Sample , compVal=0 ,
                xlab=paste( round(contrast[incIdx],2) , x2names[incIdx] ,
                            c(rep("+",sum(incIdx)-1),"") , collapse=" " ) ,
                cex.lab = 1.0 ,
                main=paste( "X2 Contrast:", names(x2contrastList)[cIdx] ) )
   }
   #saveGraph(file=paste(fileNameRoot,"x2Contrasts",sep=""),type="eps")
}
#
nContrasts = length( x1x2contrastList )
if ( nContrasts > 0 ) {
   nPlotPerRow = 5
   nPlotRow = ceiling(nContrasts/nPlotPerRow)
   nPlotCol = ceiling(nContrasts/nPlotRow)
   openGraph(3.75*nPlotCol,2.5*nPlotRow)
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
                compVal=0 , xlab=contrastLab , cex.lab = 0.75 ,
                main=paste( names(x1x2contrastList)[cIdx] ) )
   }
   #saveGraph(file=paste(fileNameRoot,"x1x2Contrasts",sep=""),type="eps")
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
openGraph((dataList$Nx1Lvl)*2.75,(dataList$Nx2Lvl)*2.25)
layoutMat = matrix( 1:(dataList$Nx2Lvl*dataList$Nx1Lvl) ,
                    nrow=(dataList$Nx2Lvl) , ncol=(dataList$Nx1Lvl) , byrow=T )
layout( layoutMat )
par( mar=c(4,1.5,2.5,0.5) , mgp=c(2,0.7,0) )
maxp = max( cellp )
for ( x2idx in 1:dataList$Nx2Lvl ) {
  for ( x1idx in 1:dataList$Nx1Lvl ) {
    histinfo = plotPost( cellp[x1idx,x2idx,] ,
               breaks=seq(0,maxp,length=51) , xlim=c(0,maxp) ,
               xlab=bquote(probability[list(x1==.(x1idx),x2==.(x2idx))]) ,
               main=paste("x1:",x1names[x1idx],", x2:",x2names[x2idx]) ,
               HDItextPlace=0.95 )
  }
}
#saveGraph(file=paste(fileNameRoot,"CellP",sep=""),type="eps")


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
