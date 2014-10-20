graphics.off()
rm(list=ls(all=TRUE))
library(BRugs)         # Kruschke, J. K. (2010). Doing Bayesian data analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
# BUGS model specification begins here...
model {
   for ( i in 1:nSubj ) {
      # Likelihood:
      nCorrOfSubj[i] ~ dbin( theta[i] , nTrlOfSubj[i] )
      # Prior on theta: Notice nested indexing.
      theta[i] ~ dbeta( aBeta[ CondOfSubj[i] ] ,
                        bBeta[ CondOfSubj[i] ] )I(0.0001,0.9999)
   }
   # Hyperprior on mu and kappa:
   kappa0 ~ dgamma( shk0[mdlIdx] , rak0[mdlIdx] )
   for ( j in 1:nCond ) {
      mu[j] ~ dbeta( aHyperbeta , bHyperbeta )
      kappa[j] ~ dgamma( shk[j,mdlIdx] , rak[j,mdlIdx] )
   }
   for ( j in 1:nCond ) {
      aBeta[j] <- mu[j]     * ((kappa[j]*(2-mdlIdx))+(kappa0*(mdlIdx-1)))
      bBeta[j] <- (1-mu[j]) * ((kappa[j]*(2-mdlIdx))+(kappa0*(mdlIdx-1)))
      # BUGS equals(,) won't work here, for no apparent reason.
      # Took me hours to isolate this problem (argh!). So, DO NOT use:
      # aBeta[j] <- mu[j]     * (kappa[j]*equals(mdlIdx,1)+kappa0*equals(mdlIdx,2))
      # bBeta[j] <- (1-mu[j]) * (kappa[j]*equals(mdlIdx,1)+kappa0*equals(mdlIdx,2))
   }
   # Constants for hyperprior:
   aHyperbeta <- 1
   bHyperbeta <- 1
   
   # Actual priors:
   shP <- 1.0 # shape for prior
   raP <- 0.1 # rate for prior
   # shape, rate kappa0[ model ]
   shk0[2] <- shP
   rak0[2] <- raP
   # shape kappa[ condition , model ]
   shk[1,1] <- shP
   shk[2,1] <- shP
   shk[3,1] <- shP
   shk[4,1] <- shP
   # rate kappa[ condition , model ]
   rak[1,1] <- raP
   rak[2,1] <- raP
   rak[3,1] <- raP
   rak[4,1] <- raP

   # Pseudo priors:
   # shape, rate kappa0[ model ]
   shk0[1] <- 54.0
   rak0[1] <- 4.35
   # shape kappa[ condition , model ]
   shk[1,2] <- 11.8
   shk[2,2] <- 11.9
   shk[3,2] <- 13.6
   shk[4,2] <- 12.6
   # rate kappa[ condition , model ]
   rak[1,2] <- 1.34
   rak[2,2] <- 1.11
   rak[3,2] <- 0.903
   rak[4,2] <- 0.748

   # Hyperprior on model index:
   mdlIdx ~ dcat( modelProb[] )
   modelProb[1] <- .5
   modelProb[2] <- .5
}
# ... end BUGS model specification
" # close quote for modelstring
# Write model to a file:
.temp = file("model.txt","w") ; writeLines(modelstring,con=.temp) ; close(.temp) 
# Load model file into BRugs and check its syntax:
modelCheck( "model.txt" )

#------------------------------------------------------------------------------
# THE DATA.

# For each subject, specify the condition s/he was in,
# the number of trials s/he experienced, and the number correct.
# (These lines intentionally exceed the margins so that they don't take up
# excessive space on the printed page.)
CondOfSubj = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)
nTrlOfSubj = c(64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64)
nCorrOfSubj = c(45,63,58,64,58,63,51,60,59,47,63,61,60,51,59,45,61,59,60,58,63,56,63,64,64,60,64,62,49,64,64,58,64,52,64,64,64,62,64,61,59,59,55,62,51,58,55,54,59,57,58,60,54,42,59,57,59,53,53,42,59,57,29,36,51,64,60,54,54,38,61,60,61,60,62,55,38,43,58,60,44,44,32,56,43,36,38,48,32,40,40,34,45,42,41,32,48,36,29,37,53,55,50,47,46,44,50,56,58,42,58,54,57,54,51,49,52,51,49,51,46,46,42,49,46,56,42,53,55,51,55,49,53,55,40,46,56,47,54,54,42,34,35,41,48,46,39,55,30,49,27,51,41,36,45,41,53,32,43,33)
nSubj = length(CondOfSubj)
nCond = length(unique(CondOfSubj))

# Specify the data in a form that is compatible with BRugs model, as a list:
datalist = list(
 nCond = nCond ,
 nSubj = nSubj ,
 CondOfSubj = CondOfSubj ,
 nTrlOfSubj = nTrlOfSubj ,
 nCorrOfSubj = nCorrOfSubj
)

# Get the data into BRugs:
# Function bugsData stores the data file (default filename is data.txt).
# Function modelData loads data file into BRugs (default filename is data.txt).
modelData( bugsData( datalist ) )

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

nchain = 1
modelCompile( numChains=nchain )
modelGenInits()

#------------------------------------------------------------------------------
# RUN THE CHAINS.

burninSteps = 1000
modelUpdate( burninSteps )
samplesSet( c("mu","kappa","kappa0","theta","mdlIdx") )
nPerChain = 10000
modelUpdate( nPerChain , thin=1 ) # takes nPerChain * thin steps

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS.

filenamebase = "FilconModelCompPseudoPriorBrugs"

modelIdxSample = samplesSample( "mdlIdx" )
pM1 = sum( modelIdxSample == 1 ) / length( modelIdxSample )
pM2 = 1 - pM1
string1 =paste("p(M1|D)=",round(pM1,3),sep="")
string2 =paste("p(M2|D)=",round(pM2,3),sep="")
windows(10,4)
plot( 1:length(modelIdxSample) , modelIdxSample , type="l" ,
      xlab="Step in Markov chain" , ylab="Model Index (1, 2)" ,
      main=paste(string1,", ",string2,sep="") )
dev.copy2eps(file=paste(filenamebase,"_mdlIdx",".eps",sep=""))

kappa0sampleM1 = samplesSample( "kappa0" )[ modelIdxSample == 1 ]
kappa0sampleM2 = samplesSample( "kappa0" )[ modelIdxSample == 2 ]
windows()
layout( matrix(1:2,nrow=2) )
hist( kappa0sampleM1 , main="Post. kappa0 for M = 1" ,
      xlab=expression(kappa[0]) , xlim=c(0,30) , freq=F , ylab="" ,
      col="grey" , border="white" , cex.lab=1.75 , breaks=c(seq(0,30,len=19),10000) )
lines( seq(0,30,.1) , dgamma( seq(0,30,.1) , 1 , .1 ) )
hist( kappa0sampleM2 , main="Post. kappa0 for M = 2"  ,
      xlab=expression(kappa[0]) , xlim=c(0,30) , freq=F , ylab="" ,
      col="grey" , border="white" , cex.lab=1.75 , breaks=c(seq(0,30,len=19),10000) )
dev.copy2eps(file=paste(filenamebase,"_k0",".eps",sep=""))

kappa1sampleM1 = samplesSample( "kappa[1]" )[ modelIdxSample == 1 ]
kappa2sampleM1 = samplesSample( "kappa[2]" )[ modelIdxSample == 1 ]
kappa3sampleM1 = samplesSample( "kappa[3]" )[ modelIdxSample == 1 ]
kappa4sampleM1 = samplesSample( "kappa[4]" )[ modelIdxSample == 1 ]
kappa1sampleM2 = samplesSample( "kappa[1]" )[ modelIdxSample == 2 ]
kappa2sampleM2 = samplesSample( "kappa[2]" )[ modelIdxSample == 2 ]
kappa3sampleM2 = samplesSample( "kappa[3]" )[ modelIdxSample == 2 ]
kappa4sampleM2 = samplesSample( "kappa[4]" )[ modelIdxSample == 2 ]
windows(10,5)
layout( matrix(1:8,nrow=2,byrow=T) )
hist( kappa1sampleM1 , main="Post. kappa[1] for M = 1" ,
      xlab=expression(kappa[1]) , xlim=c(0,30) , freq=F , ylab="" ,
      col="grey" , border="white" , cex.lab=1.75 , breaks=c(seq(0,30,len=19),10000) )
hist( kappa2sampleM1 , main="Post. kappa[2] for M = 1" ,
      xlab=expression(kappa[2]) , xlim=c(0,30) , freq=F , ylab="" ,
      col="grey" , border="white" , cex.lab=1.75 , breaks=c(seq(0,30,len=19),10000) )
hist( kappa3sampleM1 , main="Post. kappa[3] for M = 1" ,
      xlab=expression(kappa[3]) , xlim=c(0,30) , freq=F , ylab="" ,
      col="grey" , border="white" , cex.lab=1.75 , breaks=c(seq(0,30,len=19),10000) )
hist( kappa4sampleM1 , main="Post. kappa[4] for M = 1" ,
      xlab=expression(kappa[4]) , xlim=c(0,30) , freq=F , ylab="" ,
      col="grey" , border="white" , cex.lab=1.75 , breaks=c(seq(0,30,len=19),10000) )
hist( kappa1sampleM2 , main="Post. kappa[1] for M = 2" ,
      xlab=expression(kappa[1]) , xlim=c(0,30) , freq=F , ylab="" ,
      col="grey" , border="white" , cex.lab=1.75 , breaks=c(seq(0,30,len=19),10000) )
lines( seq(0,30,.1) , dgamma( seq(0,30,.1) , 1 , .1 ) )
hist( kappa2sampleM2 , main="Post. kappa[2] for M = 2" ,
      xlab=expression(kappa[2]) , xlim=c(0,30) , freq=F , ylab="" ,
      col="grey" , border="white" , cex.lab=1.75 , breaks=c(seq(0,30,len=19),10000) )
lines( seq(0,30,.1) , dgamma( seq(0,30,.1) , 1 , .1 ) )
hist( kappa3sampleM2 , main="Post. kappa[3] for M = 2" ,
      xlab=expression(kappa[3]) , xlim=c(0,30) , freq=F , ylab="" ,
      col="grey" , border="white" , cex.lab=1.75 , breaks=c(seq(0,30,len=19),10000) )
lines( seq(0,30,.1) , dgamma( seq(0,30,.1) , 1 , .1 ) )
hist( kappa4sampleM2 , main="Post. kappa[4] for M = 2" ,
      xlab=expression(kappa[4]) , xlim=c(0,30) , freq=F , ylab="" ,
      col="grey" , border="white" , cex.lab=1.75 , breaks=c(seq(0,30,len=19),10000) )
lines( seq(0,30,.1) , dgamma( seq(0,30,.1) , 1 , .1 ) )
dev.copy2eps(file=paste(filenamebase,"_kcond",".eps",sep=""))

source("plotPost.R")
windows(10,5)
layout( matrix(1:6,nrow=2,byrow=T) )
histInfo = plotPost( kappa1sampleM1 - kappa2sampleM1 , cex.lab=2 ,
        xlab=bquote(kappa[1]        -kappa[2]) , compVal=0 , breaks=20 )
histInfo = plotPost( kappa1sampleM1 - kappa3sampleM1 , cex.lab=2 ,
        xlab=bquote(kappa[1]        -kappa[3]) , compVal=0 , breaks=20 )
histInfo = plotPost( kappa1sampleM1 - kappa4sampleM1 , cex.lab=2 ,
        xlab=bquote(kappa[1]        -kappa[4]) , compVal=0 , breaks=20 )
histInfo = plotPost( kappa2sampleM1 - kappa3sampleM1 , cex.lab=2 ,
        xlab=bquote(kappa[2]        -kappa[3]) , compVal=0 , breaks=20 )
histInfo = plotPost( kappa2sampleM1 - kappa4sampleM1 , cex.lab=2 ,
        xlab=bquote(kappa[2]        -kappa[4]) , compVal=0 , breaks=20 )
histInfo = plotPost( kappa3sampleM1 - kappa4sampleM1 , cex.lab=2 ,
        xlab=bquote(kappa[3]        -kappa[4]) , compVal=0 , breaks=20 )
dev.copy2eps(file=paste(filenamebase,"_kdiff",".eps",sep=""))
