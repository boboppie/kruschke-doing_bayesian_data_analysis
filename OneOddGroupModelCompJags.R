graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="OneOddGroupModelCompJags" # for constructing output filenames
source("openGraphSaveGraph.R")
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
model {
   for ( i in 1:nSubj ) {
      # Likelihood:
      nCorrOfSubj[i] ~ dbin( theta[i] , nTrlOfSubj[i] )
      # Prior on theta (notice nested indexing):
      theta[i] ~ dbeta( aBeta[ CondOfSubj[i] ] , bBeta[ CondOfSubj[i] ] )T(0.0001,0.9999)
   }
   # Re-parameterization of aBeta[j],bBeta[j] in terms of mu and kappa:
   for ( j in 1:nCond ) {
      # Model 1: Distinct mu[j] each group.  Model 2: Shared mu0 all groups.
      aBeta[j] <-    ( mu[j]*equals(mdlIdx,1) + mu0*equals(mdlIdx,2) ) *kappa[j]
      bBeta[j] <- (1-( mu[j]*equals(mdlIdx,1) + mu0*equals(mdlIdx,2) ))*kappa[j]
   }
   # Hyperpriors for mu and kappa:
   for ( j in 1:nCond ) {
      mu[j] ~ dbeta( a[j,mdlIdx] , b[j,mdlIdx] )
   }
   for ( j in 1:nCond ) {
      kappa[j] ~ dgamma( shk , rak )
   }
   mu0 ~ dbeta( a0[mdlIdx] , b0[mdlIdx] )

   # Constants for hyperprior:
   # (There is no higher-level distribution of across-group relationships,
   # merely to keep the focus here on model comparison.)
   shk <- 1.0
   rak <- 0.1
   aP <- 1
   bP <- 1

   a0[1] <- .53*400     # pseudo
   b0[1] <- (1-.53)*400 # pseudo
   
   a0[2] <- aP # true
   b0[2] <- bP # true

   a[1,1] <- aP # true
   a[2,1] <- aP # true
   a[3,1] <- aP # true
   a[4,1] <- aP # true
   b[1,1] <- bP # true
   b[2,1] <- bP # true
   b[3,1] <- bP # true
   b[4,1] <- bP # true

   a[1,2] <- .61*100     # pseudo
   a[2,2] <- .50*100      # pseudo
   a[3,2] <- .49*100      # pseudo
   a[4,2] <- .51*100      # pseudo
   b[1,2] <- (1-.61)*100 # pseudo
   b[2,2] <- (1-.50)*100  # pseudo
   b[3,2] <- (1-.49)*100  # pseudo
   b[4,2] <- (1-.51)*100  # pseudo

   # Hyperprior on model index:
   mdlIdx ~ dcat( modelProb[] )
   modelProb[1] <- .5
   modelProb[2] <- .5
}
" # close quote for modelstring
writeLines(modelstring,con="model.txt") # Write model to a file.

#------------------------------------------------------------------------------
# THE DATA.

# For each subject, specify the condition s/he was in,
# the number of trials s/he experienced, and the number correct.
# (Randomly generated fictitious data.)
npg = 20  # number of subjects per group
ntrl = 20 # number of trials per subject
CondOfSubj = c( rep(1,npg) , rep(2,npg) , rep(3,npg) , rep(4,npg) )
nTrlOfSubj = rep( ntrl , 4*npg )
set.seed(47401)
nCorrOfSubj = c( rbinom(npg,ntrl,.61) , rbinom(npg,ntrl,.50) ,
                 rbinom(npg,ntrl,.49) , rbinom(npg,ntrl,.51) )
nSubj = length(CondOfSubj)
nCond = length(unique(CondOfSubj))
# Display mean number correct in each group:
for ( condIdx in 1:nCond ) {
    show( mean( nCorrOfSubj[ CondOfSubj==condIdx ] ) )
}

# Specify the data in a form that is compatible with BRugs model, as a list:
dataList = list(
 nCond = nCond ,
 nSubj = nSubj ,
 CondOfSubj = CondOfSubj ,
 nTrlOfSubj = nTrlOfSubj ,
 nCorrOfSubj = nCorrOfSubj
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

# Specific initialization not necessary in this case.
# initsList = list( ** ) 

#------------------------------------------------------------------------------
# RUN THE CHAINS.

parameters = c("mu","kappa","mu0","theta","mdlIdx") 
adaptSteps = 1000           # Number of steps to "tune" the samplers.
burnInSteps = 5000          # Number of steps to "burn-in" the samplers.
nChains = 3                 # Number of chains to run.
numSavedSteps=50000         # Total number of steps in chains to save.
thinSteps=1                 # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
jagsModel = jags.model( "model.txt" , data=dataList , # inits=initsList , 
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
# EXAMINE THE RESULTS.

checkConvergence = F
if ( checkConvergence ) {
  show( summary( codaSamples ) )
  plot( codaSamples , ask=T )  
  autocorr.plot( codaSamples , ask=T )
}

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )

# Display the model index
modelIdxSample = mcmcChain[, "mdlIdx" ]
pM1 = sum( modelIdxSample == 1 ) / length( modelIdxSample )
pM2 = 1 - pM1
string1 =paste("p(DiffMu|D)=",round(pM1,3),sep="")
string2 =paste("p(SameMu|D)=",round(pM2,3),sep="")
openGraph(10,4)
nStepsToPlot = 1000
plot( 1:nStepsToPlot , modelIdxSample[1:nStepsToPlot] , type="l" ,
      xlab="Step in Markov chain" , ylab="Model Index (1, 2)" ,
      main=paste(string1,", ",string2,sep="") )
saveGraph(file=paste(fileNameRoot,"_mdlIdx",".eps",sep="") , type="eps")

# Display the mu0 posterior
mu0sampleM1 = mcmcChain[, "mu0" ][ modelIdxSample == 1 ]
mu0sampleM2 = mcmcChain[, "mu0" ][ modelIdxSample == 2 ]
openGraph()
layout( matrix(1:2,nrow=2) )
hist( mu0sampleM1 , main="Post. mu0 for M = 1 (DiffMu)" ,
      xlab=expression(mu[0]) , freq=F , xlim=c(0,1) ,
      col="grey" , border="white" )
hist( mu0sampleM2 , main="Post. mu0 for M = 2 (SameMu)"  ,
      xlab=expression(mu[0]) , freq=F , xlim=c(0,1) ,
      col="grey" , border="white" )
saveGraph(file=paste(fileNameRoot,"_mu0",".eps",sep="") , type="eps")

# Display the mu[j] posterior
mu1sampleM1 = mcmcChain[, "mu[1]" ][ modelIdxSample == 1 ]
mu2sampleM1 = mcmcChain[, "mu[2]" ][ modelIdxSample == 1 ]
mu3sampleM1 = mcmcChain[, "mu[3]" ][ modelIdxSample == 1 ]
mu4sampleM1 = mcmcChain[, "mu[4]" ][ modelIdxSample == 1 ]
mu1sampleM2 = mcmcChain[, "mu[1]" ][ modelIdxSample == 2 ]
mu2sampleM2 = mcmcChain[, "mu[2]" ][ modelIdxSample == 2 ]
mu3sampleM2 = mcmcChain[, "mu[3]" ][ modelIdxSample == 2 ]
mu4sampleM2 = mcmcChain[, "mu[4]" ][ modelIdxSample == 2 ]
openGraph(10,5)
layout( matrix(1:8,nrow=2,byrow=T) )
hist( mu1sampleM1 , main="Post. mu[1] for M = 1 (DiffMu)" ,
      xlab=expression(mu[1]) , freq=F , xlim=c(0,1) ,
      col="grey" , border="white" )
hist( mu2sampleM1 , main="Post. mu[2] for M = 1 (DiffMu)" ,
      xlab=expression(mu[2]) ,  freq=F , xlim=c(0,1) ,
      col="grey" , border="white" )
hist( mu3sampleM1 , main="Post. mu[3] for M = 1 (DiffMu)" ,
      xlab=expression(mu[3]) , freq=F , xlim=c(0,1) ,
      col="grey" , border="white" )
hist( mu4sampleM1 , main="Post. mu[4] for M = 1 (DiffMu)" ,
      xlab=expression(mu[4]) , freq=F , xlim=c(0,1) ,
      col="grey" , border="white" )
hist( mu1sampleM2 , main="Post. mu[1] for M = 2 (SameMu)" ,
      xlab=expression(mu[1]) , freq=F , xlim=c(0,1) ,
      col="grey" , border="white" )
hist( mu2sampleM2 , main="Post. mu[2] for M = 2 (SameMu)" ,
      xlab=expression(mu[2]) , freq=F , xlim=c(0,1) ,
      col="grey" , border="white" )
hist( mu3sampleM2 , main="Post. mu[3] for M = 2 (SameMu)" ,
      xlab=expression(mu[3]) , freq=F , xlim=c(0,1) ,
      col="grey" , border="white" )
hist( mu4sampleM2 , main="Post. mu[4] for M = 2 (SameMu)" ,
      xlab=expression(mu[4]) , freq=F , xlim=c(0,1) ,
      col="grey" , border="white" )
saveGraph(file=paste(fileNameRoot,"_mucond",".eps",sep="") , type="eps")

# Display the differences of mu[j]'s
muSample = rbind( mu1sampleM1 , mu2sampleM1 , mu3sampleM1 , mu4sampleM1 )
source("plotPost.R")
openGraph(10,5)
layout( matrix(1:6,nrow=2,ncol=3,byrow=T) )
xmin = -0.25
xmax = 0.25
for ( i in 1:3 ) {
    for ( j in (i+1):4 ) {
        plotPost( muSample[i,]-muSample[j,] , compVal=0.0 ,
                  xlab=bquote(mu[.(i)]-mu[.(j)]) ,
                  breaks=unique( c( min(c(xmin,muSample[i,]-muSample[j,])),
                            seq(xmin,xmax,len=20),
                            max(c(xmax,muSample[i,]-muSample[j,])) )) ,
                  main="" , xlim=c(xmin,xmax) )
    }
}
saveGraph(file=paste(fileNameRoot,"_mudiff",".eps",sep="") , type="eps")
