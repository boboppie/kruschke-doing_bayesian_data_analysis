graphics.off()
rm(list=ls(all=TRUE))
source("openGraphSaveGraph.R")
source("plotPost.R")
fileNameRoot="ANOVAonewayNonhomogvarJagsSTZ" # for constructing output filenames
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
dataSource = c( "McDonaldSK1991" , "SolariLS2008" , "Random" , "Nonhomogvar" )[1]
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

if ( dataSource == "Nonhomogvar" ) {
  fileNameRoot = paste( fileNameRoot , dataSource , sep="" )
  datarecord = read.csv( "NonhomogVarData.csv" )
  y = datarecord$Y
  Ntotal = length(y)
  x = as.numeric(datarecord$Group)
  xnames = levels(datarecord$Group)
  NxLvl = length(levels(datarecord$Group))
  normalize = function( v ){ return( v / sum(v) ) }
  contrastList = list( 
    BvA = (xnames=="B")-(xnames=="A") ,
    CvA = (xnames=="C")-(xnames=="A") ,
    DvA = (xnames=="D")-(xnames=="A") , # !
    CvB = (xnames=="C")-(xnames=="B") , # !
    DvB = (xnames=="D")-(xnames=="B") ,
    DvC = (xnames=="D")-(xnames=="C") 
    )
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
adaptSteps = 1000              # Number of steps to "tune" the samplers.
burnInSteps = 5000            # Number of steps to "burn-in" the samplers.
nChains = 3                   # Number of chains to run.
numSavedSteps=100000           # Total number of steps in chains to save.
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

checkConvergence = FALSE
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
sigmaSample = 1/sqrt(tauSample) * ySDorig

# Plot aSD
openGraph(width=7,height=7)
layout( matrix(1:2,nrow=2) )
par( mar=c(3,1,2.5,0) , mgp=c(2,0.7,0) )
plotPost( aSDSample , xlab="aSD" , main="a SD"  , showMode=T )
saveGraph(file=paste(fileNameRoot,"SD",sep=""),type="eps")

# Plot b values:
openGraph(width=dataList$NxLvl*2.75,height=2.5)
layout( matrix( 1:dataList$NxLvl , nrow=1 ) )
par( mar=c(3,1,2.5,0) , mgp=c(2,0.7,0) )
for ( xidx in 1:dataList$NxLvl ) {
    plotPost( bSample[xidx,]  ,
              xlab=bquote(beta*1[.(xidx)]) ,
              main=paste("x:",xnames[xidx])  )
}
saveGraph(file=paste(fileNameRoot,"b",sep=""),type="eps")

# Plot tau values:
openGraph(width=dataList$NxLvl*2.75,height=2.5)
layout( matrix( 1:dataList$NxLvl , nrow=1 ) )
par( mar=c(3,1,2.5,0) , mgp=c(2,0.7,0) )
for ( xidx in 1:dataList$NxLvl ) {
    plotPost( tauSample[xidx,]  ,
              xlab=bquote(tau[.(xidx)]) ,
              main=paste("x:",xnames[xidx]) , showMode=T )
}
saveGraph(file=paste(fileNameRoot,"tau",sep=""),type="eps")

# Display contrast analyses
nContrasts = length( contrastList )
if ( nContrasts > 0 ) {
   nPlotPerRow = 5
   nPlotRow = ceiling(nContrasts/nPlotPerRow)
   nPlotCol = ceiling(nContrasts/nPlotRow)
   openGraph(width=3.75*nPlotCol,height=2.5*nPlotRow)
   layout( matrix(1:(nPlotRow*nPlotCol),nrow=nPlotRow,ncol=nPlotCol,byrow=T) )
   par( mar=c(4,0.5,2.5,0.5) , mgp=c(2,0.7,0) )
   for ( cIdx in 1:nContrasts ) {
       contrast = matrix( contrastList[[cIdx]],nrow=1) # make it a row matrix
       incIdx = contrast!=0
       histInfo = plotPost( contrast %*% bSample , compVal=0  ,
                xlab=paste( round(contrast[incIdx],2) , xnames[incIdx] ,
                            c(rep("+",sum(incIdx)-1),"") , collapse=" " ) ,
                cex.lab = 1.0 ,
                main=paste( "X Contrast:", names(contrastList)[cIdx] ) )
   }
   saveGraph(file=paste(fileNameRoot,"xContrasts",sep=""),type="eps")
}

# Display data with posterior predictive distributions
openGraph(width=1.5*NxLvl,height=5)
plot(0,0, 
     xlim=c(0.2,NxLvl+0.1) , xlab="X" , 
     xaxt="n" ,
     ylim=c(min(y)-0.2*(max(y)-min(y)),max(y)+0.2*(max(y)-min(y))) , ylab="Y" ,
     main="Data with Posterior Predictive Distrib.")
axis( 1 , at=1:NxLvl , lab=xnames )
for ( j in 1:NxLvl ) {
  yVals = y[x==j]
  points( rep(j,length(yVals))+runif(length(yVals),-0.03,0.03) , 
          yVals , pch=20 , cex=1.5 , col="red" )
  chainSub = round(seq(1,chainLength,length=20))
  for ( chnIdx in chainSub ) {
    m = b0Sample[chnIdx] + bSample[j,chnIdx]
    s = sigmaSample[j,chnIdx]
    yl = m-1.96*s
    yh = m+1.96*s
    ycomb=seq(yl,yh,length=201)
    ynorm = dnorm(ycomb,mean=m,sd=s)
    ynorm = 0.75*ynorm/max(ynorm)
    lines( j-ynorm , ycomb , col="skyblue" ) # col=chnIdx )
  }
}
saveGraph(file=paste(fileNameRoot,"PostPred",sep=""),type="eps")

