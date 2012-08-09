graphics.off()
rm(list=ls(all=TRUE))
library(BRugs)         # Kruschke, J. K. (2010). Doing Bayesian data analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
# BUGS model specification begins here...
model {
   for ( i in 1:nFlip ) {
      # Likelihood:
      y[i] ~ dbern( theta )
   }
   # Prior
   theta <- ( (2-mdlIdx) * 1/(1+exp( -nu )) # theta from model index 1
            + (mdlIdx-1) * exp( -eta ) )    # theta from model index 2
   nu ~ dnorm(0,.1)      #  0,.1  vs  1,1
   eta ~ dgamma(.1,.1)   # .1,.1  vs  1,1
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

# Specify the data in a form that is compatible with BRugs model, as a list:
N = 30
z = 8
datalist = list(
   y = c( rep(1,z) , rep(0,N-z) ) ,
   nFlip = N
)
# Get the data into BRugs:
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
samplesSet( c("theta","nu","eta","mdlIdx") )
nPerChain = 10000
modelUpdate( nPerChain , thin=5 ) # takes nPerChain * thin steps

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS.

filenamebase = "ToyModelComp1"

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

thetaSampleM1 = samplesSample( "theta" )[ modelIdxSample == 1 ]
thetaSampleM2 = samplesSample( "theta" )[ modelIdxSample == 2 ]
source("plotPost.R")
windows()
layout( matrix(1:2,nrow=2) )
h1 = plotPost( thetaSampleM1 , main="Post. theta for M1" , breaks=21 )
h2 = plotPost( thetaSampleM2 , main="Post. theta for M2" , breaks=21 )
dev.copy2eps(file=paste(filenamebase,"_theta",".eps",sep=""))

nuSampleM1 = samplesSample( "nu" )[ modelIdxSample == 1 ]
etaSampleM2 = samplesSample( "eta" )[ modelIdxSample == 2 ]
windows()
layout( matrix(1:2,nrow=2) )
h1 = plotPost( nuSampleM1 ,
               main=bquote("p("*nu*"|D,M1), with p(M1|D)="*.(round(pM1,3))) ,
               breaks=21 , xlab=expression(nu) , xlim=c(-3,4) )
h2 = plotPost( etaSampleM2 ,
               main=bquote("p("*eta*"|D,M2), with p(M2|D)="*.(round(pM2,3))) ,
               breaks=seq(0,50,.25) , xlab=expression(eta) , xlim=c(0,7) )
dev.copy2eps(file=paste(filenamebase,"_nu_eta",".eps",sep=""))