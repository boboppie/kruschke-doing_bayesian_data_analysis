graphics.off()
rm(list=ls(all=TRUE))
fnroot = "LogisticOnewayAnovaHeteroVarBrugs"
library(BRugs)         # Kruschke, J. K. (2010). Doing Bayesian data analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
# BUGS model specification begins here...
model {
  for ( i in 1:Ntotal ) {
    z[i] ~ dbin( theta[i] , n[i] )  
    theta[i] ~ dbeta( aBeta[x[i]] , bBeta[x[i]] )I(0.001,0.999)
  }
  for ( j in 1:NxLvl ) {
    aBeta[j] <- mu[j] * k[j]
    bBeta[j] <- (1-mu[j]) * k[j]
    mu[j] <- 1 / ( 1 + exp( -( a0 + a[j] ) ) )
    a[j] ~ dnorm( 0.0 , atau )
    k[j] ~ dgamma( skappa , rkappa )
  }
  a0 ~ dnorm( 0 , 0.001 )
  atau <- 1 / pow( aSD , 2 )
  aSD <- abs( aSDunabs ) + .1
  aSDunabs ~ dt( 0 , 0.001 , 2 )
  skappa <- pow(mg,2)/pow(sg,2)
  rkappa <- mg/pow(sg,2)
  mg ~ dunif(0,50)
  sg ~ dunif(0,30)
}
# ... end BUGS model specification
" # close quote for modelstring
# Write model to a file, and send to BUGS:
writeLines(modelstring,con="model.txt")
modelCheck( "model.txt" )

#------------------------------------------------------------------------------
# THE DATA.

# Specify data source:
dataSource = c( "Filcon" , "Relshift" , "Random" )[1]
# Load the data:

sigmoid = function( x ) { return( 1 / ( 1 + exp( -x ) ) ) }
logit = function( y ) { return( log( y / (1-y) ) ) }

if ( dataSource == "Filcon" ) {
  fnroot = paste( fnroot , dataSource , sep="" )
  x = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)
  n = c(64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64)
  z = c(45,63,58,64,58,63,51,60,59,47,63,61,60,51,59,45,61,59,60,58,63,56,63,64,64,60,64,62,49,64,64,58,64,52,64,64,64,62,64,61,59,59,55,62,51,58,55,54,59,57,58,60,54,42,59,57,59,53,53,42,59,57,29,36,51,64,60,54,54,38,61,60,61,60,62,55,38,43,58,60,44,44,32,56,43,36,38,48,32,40,40,34,45,42,41,32,48,36,29,37,53,55,50,47,46,44,50,56,58,42,58,54,57,54,51,49,52,51,49,51,46,46,42,49,46,56,42,53,55,51,55,49,53,55,40,46,56,47,54,54,42,34,35,41,48,46,39,55,30,49,27,51,41,36,45,41,53,32,43,33)
  Ntotal = length(z)
  xnames = c("FiltLR","FiltHt","Condns1","Condns2")
  NxLvl = length(unique(x))
  contrastList = list( FiltLRvFiltHt = c(1,-1,0,0) ,
                       Cond1vCond2 = c(0,0,1,-1) ,
                       FiltvCond = c(1/2,1/2,-1/2,-1/2) )
}

if ( dataSource == "Relshift" ) {
  fnroot = paste( fnroot , dataSource , sep="" )
  #source( "Kruschke1996CSdata.R" ) # if it has not yet been run
  load("Kruschke1996CSdatsum.Rdata") # loads CondOfSubj, nCorrOfSubj, nTrlOfSubj
  x = CondOfSubj
  n = nTrlOfSubj
  z = nCorrOfSubj
  Ntotal = length(z)
  xnames = c("Rev","Rel","Irr","Cmp")
  NxLvl = length(unique(x))
  contrastList = list( REVvREL = c(1,-1,0,0) , RELvIRR = c(0,1,-1,0) ,
                       IRRvCMP = c(0,0,1,-1) , CMPvOneRel = c(0,-1/2,-1/2,1) ,
                       FourExvEightEx = c(-1,1/3,1/3,1/3) ,
                       OneRelvTwoRel = c(-1/2,1/2,1/2,-1/2) )
}

if ( dataSource == "Random" ) {
  fnroot = paste( fnroot , dataSource , sep="" )
  #set.seed(47405)
  a0true = -0.5
  atrue = c( 0.8 , -0.3 , -0.5 ) # sum to zero
  ktrue = 100
  subjPerCell = 50
  nPerSubj = 100
  datarecord = matrix( 0, ncol=3 , nrow=length(atrue)*subjPerCell )
  colnames(datarecord) = c("x","z","n")
  rowidx = 0
  for ( xidx in 1:length(atrue) ) {
    for ( subjidx in 1:subjPerCell ) {
      rowidx = rowidx + 1
      datarecord[rowidx,"x"] = xidx
      mu = sigmoid(a0true+atrue[xidx])
      theta = rbeta( 1 , mu*ktrue , (1-mu)*ktrue )
      datarecord[rowidx,"z"] = rbinom( 1 , prob=theta , size=nPerSubj )
      datarecord[rowidx,"n"] = nPerSubj                                       
    }
  }
  datarecord = data.frame( x=as.factor(datarecord[,"x"]) , z=datarecord[,"z"] , 
                           n=datarecord[,"n"] )
  z = as.numeric(datarecord$z)
  Ntotal = length(z)
  n = as.numeric(datarecord$n)
  x = as.numeric(datarecord$x)
  xnames = levels(datarecord$x)
  NxLvl = length(unique(x))
  # Construct list of all pairwise comparisons, to compare with NHST TukeyHSD:
  contrastList = NULL
  for ( g1idx in 1:(NxLvl-1) ) {
    for ( g2idx in (g1idx+1):NxLvl ) {
      cmpVec = rep(0,NxLvl)
      cmpVec[g1idx] = -1
      cmpVec[g2idx] = 1
      contrastList = c( contrastList , list( cmpVec ) )
    }
  }
}

# Specify the data in a form that is compatible with BRugs model, as a list:
datalist = list(
  z = z ,
  n = n ,
  x = x ,
  Ntotal = Ntotal ,
  NxLvl = NxLvl
)
# Get the data into BRugs:
modelData( bugsData( datalist ) )

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

# Autocorrelation within chains is large, so use several chains to reduce
# degree of thinning. But we still have to burn-in all the chains, which takes
# more time with more chains.
nchain = 3
modelCompile( numChains = nchain )

if ( F ) {
   modelGenInits() # often won't work for diffuse prior
} else {
  #  initialization based on data
  theData = data.frame( pr=.01+.98*datalist$z/datalist$n , 
                        x=factor(x,labels=xnames) )
  a0 = mean( logit(theData$pr) )
  a = aggregate( logit(theData$pr) , list( theData$x ) , mean )[,2] - a0
  mGrp = aggregate( theData$pr , list( theData$x ) , mean )[,2]
  sdGrp = aggregate( theData$pr , list( theData$x ) , sd )[,2]
  kGrp = mGrp*(1-mGrp)/sdGrp^2 - 1
  k = mean(kGrp)
  genInitList <- function() {
    return(
        list(
            a0 = a0 ,
            a = a ,
            aSDunabs = sd(a) ,
            theta = theData$pr ,
            k = kGrp ,
            mg = 10 ,
            sg = 10
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
BurnInSteps = 5000
modelUpdate( BurnInSteps )
# actual samples
samplesSet( c( "a0" ,  "a" , "aSD" , "k" , "mg" , "sg" ) )
stepsPerChain = ceiling(2000/nchain)
thinStep = 500 # 750
modelUpdate( stepsPerChain , thin=thinStep )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

source("plotChains.R")
source("plotPost.R")

checkConvergence = T
if ( checkConvergence ) {
   sumInfo = plotChains( "a0" , saveplots=F , filenameroot=fnroot )
   sumInfo = plotChains( "a" , saveplots=F , filenameroot=fnroot )
   sumInfo = plotChains( "aSD" , saveplots=F , filenameroot=fnroot )
   sumInfo = plotChains( "k" , saveplots=F , filenameroot=fnroot )
   sumInfo = plotChains( "mg" , saveplots=F , filenameroot=fnroot )
   sumInfo = plotChains( "sg" , saveplots=F , filenameroot=fnroot )
}

###### NEEDS MODIFICATION BELOW THIS POINT! #########################

# Extract and plot the SDs:
aSDSample = samplesSample("aSD")
windows()
par( mar=c(3,1,2.5,0) , mgp=c(2,0.7,0) )
histInfo = plotPost( aSDSample , xlab="aSD" , main="a SD" , breaks=30 )
dev.copy2eps(file=paste(fnroot,"SD.eps",sep=""))

# Extract a values:
a0Sample = samplesSample( "a0" )
chainLength = length(a0Sample)
aSample = array( 0 , dim=c( datalist$NxLvl , chainLength ) )
for ( xidx in 1:datalist$NxLvl ) {
   aSample[xidx,] = samplesSample( paste("a[",xidx,"]",sep="") )
}

# Convert to zero-centered b values:
mSample = array( 0, dim=c( datalist$NxLvl , chainLength ) )
for ( stepIdx in 1:chainLength ) {
    mSample[,stepIdx ] = ( a0Sample[stepIdx] + aSample[,stepIdx] )
}
b0Sample = apply( mSample , 2 , mean )
bSample = mSample - matrix(rep( b0Sample ,NxLvl),nrow=NxLvl,byrow=T)

# Plot b values:
windows(datalist$NxLvl*2.75,2.5)
layout( matrix( 1:datalist$NxLvl , nrow=1 ) )
par( mar=c(3,1,2.5,0) , mgp=c(2,0.7,0) )
for ( xidx in 1:datalist$NxLvl ) {
    plotPost( bSample[xidx,] , breaks=30 ,
              xlab=bquote(beta[.(xidx)]) ,
              main=paste(xnames[xidx])  )
}
dev.copy2eps(file=paste(fnroot,"b.eps",sep=""))

# Consider parameter correlations:
kSample = samplesSample("k")
windows()
pairs( cbind( b0Sample , t(bSample) , kSample ) , labels=c("b0",xnames,"k") )

# Display contrast analyses
nContrasts = length( contrastList )
if ( nContrasts > 0 ) {
   nPlotPerRow = 5
   nPlotRow = ceiling(nContrasts/nPlotPerRow)
   nPlotCol = ceiling(nContrasts/nPlotRow)
   windows(3.75*nPlotCol,2.5*nPlotRow)
   layout( matrix(1:(nPlotRow*nPlotCol),nrow=nPlotRow,ncol=nPlotCol,byrow=T) )
   par( mar=c(4,0.5,2.5,0.5) , mgp=c(2,0.7,0) )
   for ( cIdx in 1:nContrasts ) {
       contrast = matrix( contrastList[[cIdx]],nrow=1) # make it a row matrix
       incIdx = contrast!=0
       histInfo = plotPost( contrast %*% bSample , compVal=0 , breaks=30 ,
                xlab=paste( round(contrast[incIdx],2) , xnames[incIdx] ,
                            c(rep("+",sum(incIdx)-1),"") , collapse=" " ) ,
                cex.lab = 1.5 ,
                main=paste( "Contrast:", names(contrastList)[cIdx] ) )
   }
   dev.copy2eps(file=paste(fnroot,"xContrasts.eps",sep=""))
}

#==============================================================================
# Do NHST ANOVA:

theData = data.frame( y=z/n , x=factor(x,labels=xnames) )
aovresult = aov( y ~ x , data = theData )
cat("\n------------------------------------------------------------------\n\n")
print( summary( aovresult ) )
cat("\n------------------------------------------------------------------\n\n")
print( model.tables( aovresult , "means" ) , digits=4 )
windows()
boxplot( y ~ x , data = theData )
cat("\n------------------------------------------------------------------\n\n")
print( TukeyHSD( aovresult , "x" , ordered = FALSE ) )
windows()
plot( TukeyHSD( aovresult , "x" ) )
if ( F ) {
  for ( xIdx1 in 1:(NxLvls-1) ) {
    for ( xIdx2 in (xIdx1+1):NxLvls ) {
      cat("\n----------------------------------------------------------\n\n")
      cat( "xIdx1 = " , xIdx1 , ", xIdx2 = " , xIdx2 ,
           ", M2-M1 = " , mean(score[x==xIdx2])-mean(score[x==xIdx1]) ,
           "\n" )
      print( t.test( score[ x == xIdx2 ] , score[ x == xIdx1 ] ) )
    }
  }
}
cat("\n------------------------------------------------------------------\n\n")

#==============================================================================
