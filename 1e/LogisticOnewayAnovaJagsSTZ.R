graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="LogisticOnewayAnovaJagsSTZ" # for constructing output filenames
source("openGraphSaveGraph.R")
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
model {
  for ( i in 1:Ntotal ) {
    z[i] ~ dbin( theta[i] , n[i] )  
    theta[i] ~ dbeta( aBeta[x[i]] , bBeta[x[i]] )T(0.001,0.999)
  }
  for ( j in 1:NxLvl ) {
    aBeta[j] <- mu[j] * k
    bBeta[j] <- (1-mu[j]) * k
    mu[j] <- 1 / ( 1 + exp( -( a0 + a[j] ) ) )
    a[j] ~ dnorm( 0.0 , atau )
  }
  k ~ dgamma( 1.0 , 0.01 )
  a0 ~ dnorm( 0 , 0.001 )
  atau <- 1 / pow( aSD , 2 )
  aSD <- abs( aSDunabs ) + .1
  aSDunabs ~ dt( 0 , 0.001 , 2 )
  # Convert a0,a[] to sum-to-zero b0,b[] :
  for ( j in 1:NxLvl ) { m[j] <- a0 + a[j] } 
  b0 <- mean( m[1:NxLvl] )
  for ( j in 1:NxLvl ) { b[j] <- m[j] - b0 }
}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.

# Specify data source:
dataSource = c( "Filcon" , "Relshift" , "Random" )[1]
# Load the data:

sigmoid = function( x ) { return( 1 / ( 1 + exp( -x ) ) ) }
logit = function( y ) { return( log( y / (1-y) ) ) }

if ( dataSource == "Filcon" ) {
  fileNameRoot = paste( fileNameRoot , dataSource , sep="" )
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
  fileNameRoot = paste( fileNameRoot , dataSource , sep="" )
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
  fileNameRoot = paste( fileNameRoot , dataSource , sep="" )
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

dataList = list(
  z = z ,
  n = n ,
  x = x ,
  Ntotal = Ntotal ,
  NxLvl = NxLvl
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

theData = data.frame( pr=.01+.98*dataList$z/dataList$n , 
                      x=factor(x,labels=xnames) )
a0 = mean( logit(theData$pr) )
a = aggregate( logit(theData$pr) , list( theData$x ) , mean )[,2] - a0
mGrp = aggregate( theData$pr , list( theData$x ) , mean )[,2]
sdGrp = aggregate( theData$pr , list( theData$x ) , sd )[,2]
kGrp = mGrp*(1-mGrp)/sdGrp^2 - 1
k = mean(kGrp)
initsList = list( a0 = a0 , a = a , aSDunabs = sd(a) , theta = theData$pr ,
                  k = k )

#------------------------------------------------------------------------------
# RUN THE CHAINS

parameters = c( "a0" ,  "a" , "b0" , "b" , "aSD" , "k" )  
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
  openGraph()
  plot( codaSamples , ask=F )  
  openGraph()
  autocorr.plot( codaSamples , ask=F )
}

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )


# Extract parameter values:
chainLength = NROW(mcmcChain)
kSample = mcmcChain[,"k"]
aSDSample = mcmcChain[,"aSD"]
b0Sample = mcmcChain[, "b0" ]
bSample = array( 0 , dim=c( dataList$NxLvl , chainLength ) )
for ( xidx in 1:dataList$NxLvl ) {
   bSample[xidx,] = mcmcChain[, paste("b[",xidx,"]",sep="") ]
}

source("plotPost.R")
# plot the SD:
openGraph()
par( mar=c(3,1,2.5,0) , mgp=c(2,0.7,0) )
histInfo = plotPost( aSDSample , xlab="aSD" , main="a SD" , breaks=30 , showMode=T)
saveGraph(file=paste(fileNameRoot,"SD.eps",sep=""),type="eps")
# Plot b values:
openGraph(dataList$NxLvl*2.75,2.5)
layout( matrix( 1:dataList$NxLvl , nrow=1 ) )
par( mar=c(3,1,2.5,0) , mgp=c(2,0.7,0) )
for ( xidx in 1:dataList$NxLvl ) {
    plotPost( bSample[xidx,] , breaks=30 ,
              xlab=bquote(beta[.(xidx)]) ,
              main=paste(xnames[xidx])  )
}
saveGraph(file=paste(fileNameRoot,"b.eps",sep=""),type="eps")

# Consider parameter correlations:
openGraph()
nPts=700 ; thinIdx=round(seq(1,chainLength,length=nPts))
pairs( cbind( b0Sample[thinIdx] , t(bSample[,thinIdx]) , kSample[thinIdx] ) ,
       labels=c("b0",xnames,"k") )

# Display contrast analyses
nContrasts = length( contrastList )
if ( nContrasts > 0 ) {
   nPlotPerRow = 5
   nPlotRow = ceiling(nContrasts/nPlotPerRow)
   nPlotCol = ceiling(nContrasts/nPlotRow)
   openGraph(3.75*nPlotCol,2.5*nPlotRow)
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
   saveGraph(file=paste(fileNameRoot,"xContrasts.eps",sep=""),type="eps")
}

#==============================================================================
# Do NHST ANOVA:

theData = data.frame( y=z/n , x=factor(x,labels=xnames) )
aovresult = aov( y ~ x , data = theData )
cat("\n------------------------------------------------------------------\n\n")
print( summary( aovresult ) )
cat("\n------------------------------------------------------------------\n\n")
print( model.tables( aovresult , "means" ) , digits=4 )
openGraph()
boxplot( y ~ x , data = theData )
cat("\n------------------------------------------------------------------\n\n")
print( TukeyHSD( aovresult , "x" , ordered = FALSE ) )
openGraph()
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
