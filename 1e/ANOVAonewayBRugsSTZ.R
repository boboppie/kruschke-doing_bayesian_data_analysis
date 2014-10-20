graphics.off()
rm(list=ls(all=TRUE))
fnroot = "ANOVAonewayBrugsSTZ"
library(BRugs)         # Kruschke, J. K. (2011). Doing Bayesian data analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
# BUGS model specification begins here...
model {
  for ( i in 1:Ntotal ) {
    y[i] ~ dnorm( mu[i] , tau )
    mu[i] <- a0 + a[x[i]]
  }
  #
  tau <- pow( sigma , -2 )
  sigma ~ dunif(0,10) # y values are assumed to be standardized
  #
  a0 ~ dnorm(0,0.001) # y values are assumed to be standardized
  #
  for ( j in 1:NxLvl ) { a[j] ~ dnorm( 0.0 , atau ) }
  atau <- 1 / pow( aSD , 2 )
  aSD <- abs( aSDunabs ) + .1
  aSDunabs ~ dt( 0 , 0.001 , 2 )
  # Convert a0,a[] to sum-to-zero b0,b[] :
  for ( j in 1:NxLvl ) { m[j] <- a0 + a[j] } 
  b0 <- mean( m[1:NxLvl] )
  for ( j in 1:NxLvl ) { b[j] <- m[j] - b0 }
}
# ... end BUGS model specification
" # close quote for modelstring
# Write model to a file, and send to BUGS:
writeLines(modelstring,con="model.txt")
modelCheck( "model.txt" )

#------------------------------------------------------------------------------
# THE DATA.

# Specify data source:
dataSource = c( "McDonaldSK1991" , "SolariLS2008" , "Random" )[1]
# Load the data:

if ( dataSource == "McDonaldSK1991" ) {
  fnroot = paste( fnroot , dataSource , sep="" )
  datarecord = read.table( "McDonaldSK1991data.txt", header=T ,
                           colClasses=c("factor","numeric") )
  y = as.numeric(datarecord$Size)
  Ntotal = length(datarecord$Size)
  x = as.numeric(datarecord$Group)
  xnames = levels(datarecord$Group)
  NxLvl = length(unique(datarecord$Group))
  contrastList = list( BIGvSMALL = c(-1/3,-1/3,1/2,-1/3,1/2) ,
                       ORE1vORE2 = c(1,-1,0,0,0) ,
                       ALAvORE = c(-1/2,-1/2,1,0,0) ,
                       NPACvORE = c(-1/2,-1/2,1/2,1/2,0) ,
                       USAvRUS = c(1/3,1/3,1/3,-1,0) ,
                       FINvPAC = c(-1/4,-1/4,-1/4,-1/4,1) ,
                       ENGvOTH = c(1/3,1/3,1/3,-1/2,-1/2) ,
                       FINvRUS = c(0,0,0,-1,1) )
}

if ( dataSource == "SolariLS2008" ) {
  fnroot = paste( fnroot , dataSource , sep="" )
  datarecord = read.table("SolariLS2008data.txt", header=T ,
                           colClasses=c("factor","numeric") )
  y = as.numeric(datarecord$Acid)
  Ntotal = length(datarecord$Acid)
  x = as.numeric(datarecord$Type)
  xnames = levels(datarecord$Type)
  NxLvl = length(unique(datarecord$Type))
  contrastList = list( G3vOTHER = c(-1/8,-1/8,1,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8) )
}

if ( dataSource == "Random" ) {
  fnroot = paste( fnroot , dataSource , sep="" )
  #set.seed(47405)
  ysdtrue = 4.0
  a0true = 100
  atrue = c( 2 , -2 ) # sum to zero
  npercell = 8
  datarecord = matrix( 0, ncol=2 , nrow=length(atrue)*npercell )
  colnames(datarecord) = c("y","x")
  rowidx = 0
  for ( xidx in 1:length(atrue) ) {
    for ( subjidx in 1:npercell ) {
      rowidx = rowidx + 1
      datarecord[rowidx,"x"] = xidx
      datarecord[rowidx,"y"] = ( a0true + atrue[xidx] + rnorm(1,0,ysdtrue) )
    }
  }
  datarecord = data.frame( y=datarecord[,"y"] , x=as.factor(datarecord[,"x"]) )
  y = as.numeric(datarecord$y)
  Ntotal = length(y)
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
ySDorig = sd(y)
yMorig = mean(y)
z = ( y - yMorig ) / ySDorig
datalist = list(
  y = z ,
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
# more time with more chains (on serial CPUs).
nchain = 5
modelCompile( numChains = nchain )

if ( F ) {
   modelGenInits() # often won't work for diffuse prior
} else {
  #  initialization based on data
  theData = data.frame( y=datalist$y , x=factor(x,labels=xnames) )
  a0 = mean( theData$y )
  a = aggregate( theData$y , list( theData$x ) , mean )[,2] - a0
  ssw = aggregate( theData$y , list( theData$x ) ,
                  function(x){var(x)*(length(x)-1)} )[,2]
  sp = sqrt( sum( ssw ) / length( theData$y ) )
  genInitList <- function() {
    return(
        list(
            a0 = a0 ,
            a = a ,
            sigma = sp ,
            aSDunabs = sd(a)
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
samplesSet( c( "a0" ,  "a" , "b0" , "b" , "sigma" , "aSD" ) )
stepsPerChain = ceiling(5000/nchain)
thinStep = 40
modelUpdate( stepsPerChain , thin=thinStep )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

source("plotChains.R")
source("plotPost.R")

checkConvergence = T
if ( checkConvergence ) {
  # Merely for insight, take a look at the a0,a[] values:
  sumInfo = plotChains( "a0" , saveplots=T , filenameroot=fnroot )
  sumInfo = plotChains( "a" , saveplots=T , filenameroot=fnroot )
   sumInfo = plotChains( "b0" , saveplots=T , filenameroot=fnroot )
   sumInfo = plotChains( "b" , saveplots=T , filenameroot=fnroot )
   readline("Press any key to clear graphics and continue")
   graphics.off()
   sumInfo = plotChains( "sigma" , saveplots=T , filenameroot=fnroot )
   sumInfo = plotChains( "aSD" , saveplots=T , filenameroot=fnroot )
}

# Extract and plot the SDs:
sigmaSample = samplesSample("sigma")
aSDSample = samplesSample("aSD")
windows()
layout( matrix(1:2,nrow=2) )
par( mar=c(3,1,2.5,0) , mgp=c(2,0.7,0) )
histInfo = plotPost( sigmaSample , xlab="sigma" , main="Cell SD" , breaks=30 , col="skyblue" )
histInfo = plotPost( aSDSample , xlab="aSD" , main="a SD" , breaks=30 , col="skyblue" )
savePlot( file=paste(fnroot,"SD.eps",sep="") , type="eps" )

readline("Press any key to clear graphics and continue")
graphics.off()

# Extract b values:
b0Sample = samplesSample( "b0" )
chainLength = length(b0Sample)
bSample = array( 0 , dim=c( datalist$NxLvl , chainLength ) )
for ( xidx in 1:datalist$NxLvl ) {
   bSample[xidx,] = samplesSample( paste("b[",xidx,"]",sep="") )
}
# Convert from standardized b values to original scale b values:
b0Sample = b0Sample * ySDorig + yMorig
bSample = bSample * ySDorig

# Plot b values:
windows(datalist$NxLvl*2.75,2.5)
layout( matrix( 1:datalist$NxLvl , nrow=1 ) )
par( mar=c(3,1,2.5,0) , mgp=c(2,0.7,0) )
for ( xidx in 1:datalist$NxLvl ) {
    histInfo = plotPost( bSample[xidx,] , breaks=30 , col="skyblue" ,
              xlab=bquote(beta*1[.(xidx)]) ,
              main=paste("x:",xnames[xidx])  )
}
savePlot( file=paste(fnroot,"b.eps",sep="") , type="eps" )

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
                cex.lab = 1.0 ,
                main=paste( "X Contrast:", names(contrastList)[cIdx] ), 
                col="skyblue" )
   }
   savePlot( file=paste(fnroot,"xContrasts.eps",sep="") , type="eps" )
}

#==============================================================================
# Do NHST ANOVA and t tests:

theData = data.frame( y=y , x=factor(x,labels=xnames) )
aovresult = aov( y ~ x , data = theData ) # NHST ANOVA
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
if ( T ) {
  for ( xIdx1 in 1:(NxLvl-1) ) {
    for ( xIdx2 in (xIdx1+1):NxLvl ) {
      cat("\n----------------------------------------------------------\n\n")
      cat( "xIdx1 = " , xIdx1 , ", xIdx2 = " , xIdx2 ,
           ", M2-M1 = " , mean(y[x==xIdx2])-mean(y[x==xIdx1]) , "\n" )
      print( t.test( y[x==xIdx2] , y[x==xIdx1] , var.equal=T ) ) # t test
    }
  }
}
cat("\n------------------------------------------------------------------\n\n")

#==============================================================================
