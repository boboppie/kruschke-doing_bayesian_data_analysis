graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="ANOVAonewayNonhomogvarJagsSTZ" # for constructing output filenames
if ( .Platform$OS.type != "windows" ) { 
  windows <- function( ... ) X11( ... ) 
}
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
model {
  for ( i in 1:Ntotal ) {
    y[i] ~ dnorm( mu[i] , tau[x[i]] )
    mu[i] <- a0 + a[x[i]]
  }
  a0 ~ dnorm(0,0.001)
  for ( j in 1:NxLvl ) {
     a[j] ~ dnorm( 0.0 , atau ) 
     tau[j] ~ dgamma( sG , rG ) 
  }
  sG <- pow(m,2)/pow(d,2)
  rG <- m/pow(d,2)
  m ~ dgamma(1,1)
  d ~ dgamma(1,1)
  atau <- 1 / pow( aSD , 2 )
  aSD <- abs( aSDunabs ) + .1
  aSDunabs ~ dt( 0 , 0.001 , 2 )
  # Convert a0,a[] to sum-to-zero b0,b[] :
  for ( j in 1:NxLvl ) { mpred[j] <- a0 + a[j] } 
  b0 <- mean( mpred[1:NxLvl] )
  for ( j in 1:NxLvl ) { b[j] <- mpred[j] - b0 }
}
" # close quote for modelstring
# Write model to a file, and send to BUGS:
writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.

# Specify data source:
dataSource = c( "McDonaldSK1991" , "SolariLS2008" , "Random" )[1]
# Load the data:

if ( dataSource == "McDonaldSK1991" ) {
  fileNameRoot = paste( fileNameRoot , dataSource , sep="" )
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
  fileNameRoot = paste( fileNameRoot , dataSource , sep="" )
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
  fileNameRoot = paste( fileNameRoot , dataSource , sep="" )
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
dataList = list(
  y = z ,
  x = x ,
  Ntotal = Ntotal ,
  NxLvl = NxLvl
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

theData = data.frame( y=dataList$y , x=factor(x,labels=xnames) )
a0 = mean( theData$y )
a = aggregate( theData$y , list( theData$x ) , mean )[,2] - a0
tau = 1/(aggregate( theData$y , list( theData$x ) , sd )[,2])^2
initsList = list( a0 = a0 , a = a , tau = tau , m = mean( tau ) , 
                  d = sd( tau ) , aSDunabs = sd(a) )

#------------------------------------------------------------------------------
# RUN THE CHAINS

parameters = c( "a0" ,  "a" , "b0" , "b" , "tau", "m", "d", "aSD" )  
adaptSteps = 500              # Number of steps to "tune" the samplers.
burnInSteps = 500            # Number of steps to "burn-in" the samplers.
nChains = 3                   # Number of chains to run.
numSavedSteps=50000           # Total number of steps in chains to save.
thinSteps=1                   # Number of steps to "thin" (1=keep every step).
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

checkConvergence = F
if ( checkConvergence ) {
  show( summary( codaSamples ) )
  windows()
  plot( codaSamples , ask=F )  
  windows()
  autocorr.plot( codaSamples , ask=F )
}

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )
chainLength = NROW(mcmcChain)

# Extract parameters:
aSDSample = mcmcChain[,"aSD"]
tauSample = array( 0 , dim=c( dataList$NxLvl , chainLength ) )
for ( xidx in 1:dataList$NxLvl ) {
   tauSample[xidx,] = mcmcChain[, paste("tau[",xidx,"]",sep="") ]
}
a0Sample = mcmcChain[, "a0" ]
aSample = array( 0 , dim=c( dataList$NxLvl , chainLength ) )
for ( xidx in 1:dataList$NxLvl ) {
   aSample[xidx,] = mcmcChain[, paste("a[",xidx,"]",sep="") ]
}
b0Sample = mcmcChain[, "b0" ]
bSample = array( 0 , dim=c( dataList$NxLvl , chainLength ) )
for ( xidx in 1:dataList$NxLvl ) {
   bSample[xidx,] = mcmcChain[, paste("b[",xidx,"]",sep="") ]
}

# Convert from standardized b values to original scale b values:
b0Sample = b0Sample * ySDorig + yMorig
bSample = bSample * ySDorig

source("plotPost.R")

# Plot aSD
windows()
layout( matrix(1:2,nrow=2) )
par( mar=c(3,1,2.5,0) , mgp=c(2,0.7,0) )
plotPost( aSDSample , xlab="aSD" , main="a SD" , breaks=30 , showMode=T )
savePlot(file=paste(fileNameRoot,"SD.eps",sep=""),type="eps")

# Plot b values:
windows(dataList$NxLvl*2.75,2.5)
layout( matrix( 1:dataList$NxLvl , nrow=1 ) )
par( mar=c(3,1,2.5,0) , mgp=c(2,0.7,0) )
for ( xidx in 1:dataList$NxLvl ) {
    plotPost( bSample[xidx,] , breaks=30 ,
              xlab=bquote(beta*1[.(xidx)]) ,
              main=paste("x:",xnames[xidx])  )
}
savePlot(file=paste(fileNameRoot,"b.eps",sep=""),type="eps")

# Plot tau values:
windows(dataList$NxLvl*2.75,2.5)
layout( matrix( 1:dataList$NxLvl , nrow=1 ) )
par( mar=c(3,1,2.5,0) , mgp=c(2,0.7,0) )
for ( xidx in 1:dataList$NxLvl ) {
    plotPost( tauSample[xidx,] , breaks=30 ,
              xlab=bquote(tau[.(xidx)]) ,
              main=paste("x:",xnames[xidx]) , showMode=T )
}
savePlot(file=paste(fileNameRoot,"tau.eps",sep=""),type="eps")

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
                main=paste( "X Contrast:", names(contrastList)[cIdx] ) )
   }
   savePlot(file=paste(fileNameRoot,"xContrasts.eps",sep=""),type="eps")
}
