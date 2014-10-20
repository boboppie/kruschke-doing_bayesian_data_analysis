graphics.off()
rm(list=ls(all=TRUE))
library(BRugs)        # Kruschke, J. K. (2010). Doing Bayesian data analysis:
                      # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
# BUGS model specification begins here...
model {
   for ( i in 1:nSubj ) {
      # Likelihood:
      nCorrOfSubj[i] ~ dbin( theta[i] , nTrlOfSubj[i] )
      # Prior on theta (notice nested indexing):
      theta[i] ~ dbeta( aBeta[ CondOfSubj[i] ] , bBeta[ CondOfSubj[i] ] )I(0.0001,0.9999)
   }
   # Re-parameterization of aBeta[j],bBeta[j] in terms of mu and kappa:
   for ( j in 1:nCond ) {
      # Model 1: Distinct mu[j] each group.  Model 2: Shared mu0 all groups.
      aBeta[j] <-       ( mu[j]*(2-mdlIdx) + mu0*(mdlIdx-1) )   * kappa[j]
      bBeta[j] <- ( 1 - ( mu[j]*(2-mdlIdx) + mu0*(mdlIdx-1) ) ) * kappa[j]
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
# ... end BUGS model specification
" # close quote for modelstring
# Write model to a file:
writeLines( text=modelstring , con="model.txt" )
# Load model file into BRugs and check its syntax:
modelCheck( "model.txt" )

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
datalist = list(
 nCond = nCond ,
 nSubj = nSubj ,
 CondOfSubj = CondOfSubj ,
 nTrlOfSubj = nTrlOfSubj ,
 nCorrOfSubj = nCorrOfSubj
)

# Get the data into BRugs:
modelData( bugsData( datalist ) )

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

nchain = 3
modelCompile( numChains=nchain )
modelGenInits()

#------------------------------------------------------------------------------
# RUN THE CHAINS.

burninSteps = 5000
modelUpdate( burninSteps )
samplesSet( c("mu","kappa","mu0","theta","mdlIdx") )
nPerChain = 5000 ; nThin = 10
modelUpdate( nPerChain , thin=nThin )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS.

filenamebase = "OneOddGroupModelComp"

# Check burnin(convergence) and clumpiness(autocorrelation):
source("plotChains.R")
plotChains("mu0")
plotChains("mu")
plotChains("kappa")

# Display the model index
modelIdxSample = samplesSample( "mdlIdx" )
pM1 = sum( modelIdxSample == 1 ) / length( modelIdxSample )
pM2 = 1 - pM1
string1 =paste("p(DiffMu|D)=",round(pM1,3),sep="")
string2 =paste("p(SameMu|D)=",round(pM2,3),sep="")
windows(10,4)
nStepsToPlot = 1000
plot( 1:nStepsToPlot , modelIdxSample[1:nStepsToPlot] , type="l" ,
      xlab="Step in Markov chain" , ylab="Model Index (1, 2)" ,
      main=paste(string1,", ",string2,sep="") )
dev.copy2eps(file=paste(filenamebase,"_mdlIdx",".eps",sep=""))

# Display the mu0 posterior
mu0sampleM1 = samplesSample( "mu0" )[ modelIdxSample == 1 ]
mu0sampleM2 = samplesSample( "mu0" )[ modelIdxSample == 2 ]
windows()
layout( matrix(1:2,nrow=2) )
hist( mu0sampleM1 , main="Post. mu0 for M = 1 (DiffMu)" ,
      xlab=expression(mu[0]) , freq=F , xlim=c(0,1) ,
      col="grey" , border="white" )
hist( mu0sampleM2 , main="Post. mu0 for M = 2 (SameMu)"  ,
      xlab=expression(mu[0]) , freq=F , xlim=c(0,1) ,
      col="grey" , border="white" )
dev.copy2eps(file=paste(filenamebase,"_mu0",".eps",sep=""))

# Display the mu[j] posterior
mu1sampleM1 = samplesSample( "mu[1]" )[ modelIdxSample == 1 ]
mu2sampleM1 = samplesSample( "mu[2]" )[ modelIdxSample == 1 ]
mu3sampleM1 = samplesSample( "mu[3]" )[ modelIdxSample == 1 ]
mu4sampleM1 = samplesSample( "mu[4]" )[ modelIdxSample == 1 ]
mu1sampleM2 = samplesSample( "mu[1]" )[ modelIdxSample == 2 ]
mu2sampleM2 = samplesSample( "mu[2]" )[ modelIdxSample == 2 ]
mu3sampleM2 = samplesSample( "mu[3]" )[ modelIdxSample == 2 ]
mu4sampleM2 = samplesSample( "mu[4]" )[ modelIdxSample == 2 ]
windows(10,5)
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
dev.copy2eps(file=paste(filenamebase,"_mucond",".eps",sep=""))

# Display the differences of mu[j]'s
muSample = rbind( mu1sampleM1 , mu2sampleM1 , mu3sampleM1 , mu4sampleM1 )
source("plotPost.R")
windows(10,5)
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
dev.copy2eps(file=paste(filenamebase,"_mudiff",".eps",sep=""))
