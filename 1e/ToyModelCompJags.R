graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="FilconModelCompJags" # for constructing output filenames
source("openGraphSaveGraph.R")
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
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
" # close quote for modelstring
writeLines(modelstring,con="model.txt") # Write model to a file.

#------------------------------------------------------------------------------
# THE DATA.

N = 30
z = 8
dataList = list(
   y = c( rep(1,z) , rep(0,N-z) ) ,
   nFlip = N
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

# Specific initialization not necessary in this case.
# initsList = list( ** ) 

#------------------------------------------------------------------------------
# RUN THE CHAINS.

parameters = c("theta","nu","eta","mdlIdx") 
adaptSteps = 1000           # Number of steps to "tune" the samplers.
burnInSteps = 1000          # Number of steps to "burn-in" the samplers.
nChains = 1                 # Number of chains to run.
numSavedSteps=10000         # Total number of steps in chains to save.
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

modelIdxSample = mcmcChain[, "mdlIdx" ]
pM1 = sum( modelIdxSample == 1 ) / length( modelIdxSample )
pM2 = 1 - pM1
string1 =paste("p(M1|D)=",round(pM1,3),sep="")
string2 =paste("p(M2|D)=",round(pM2,3),sep="")
openGraph(10,4)
plot( 1:length(modelIdxSample) , modelIdxSample , type="l" ,
      xlab="Step in Markov chain" , ylab="Model Index (1, 2)" ,
      main=paste(string1,", ",string2,sep="") )
saveGraph(file=paste(fileNameRoot,"_mdlIdx",".eps",sep="") , type="eps")

thetaSampleM1 = mcmcChain[, "theta" ][ modelIdxSample == 1 ]
thetaSampleM2 = mcmcChain[, "theta" ][ modelIdxSample == 2 ]
source("plotPost.R")
openGraph()
layout( matrix(1:2,nrow=2) )
h1 = plotPost( thetaSampleM1 , main="Post. theta for M1" , breaks=21 )
h2 = plotPost( thetaSampleM2 , main="Post. theta for M2" , breaks=21 )
saveGraph(file=paste(fileNameRoot,"_theta",".eps",sep="") , type="eps")

nuSampleM1 = mcmcChain[, "nu" ][ modelIdxSample == 1 ]
etaSampleM2 = mcmcChain[, "eta" ][ modelIdxSample == 2 ]
openGraph()
layout( matrix(1:2,nrow=2) )
h1 = plotPost( nuSampleM1 ,
               main=bquote("p("*nu*"|D,M1), with p(M1|D)="*.(round(pM1,3))) ,
               breaks=21 , xlab=expression(nu) , xlim=c(-3,4) )
h2 = plotPost( etaSampleM2 ,
               main=bquote("p("*eta*"|D,M2), with p(M2|D)="*.(round(pM2,3))) ,
               breaks=seq(0,50,.25) , xlab=expression(eta) , xlim=c(0,7) )
saveGraph(file=paste(fileNameRoot,"_nu_eta",".eps",sep="") , type="eps")