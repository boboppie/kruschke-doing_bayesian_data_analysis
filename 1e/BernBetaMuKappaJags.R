rm(list = ls())
graphics.off()
source("openGraphSaveGraph.R")
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

# Specify the model in JAGS language, but save it as a string in R:
modelString = "
# JAGS model specification begins ...
model {
    # Likelihood:
    for ( t in 1:nTrialTotal ) {
        y[t] ~ dbern( theta[ coin[ t ] ] )
    }
    # Prior:
    for ( j in 1:nCoins ) {
        theta[j] ~ dbeta( a , b )T(0.001,0.999)
    }
    a <- mu * kappa
    b <- ( 1.0 - mu ) * kappa
    mu ~ dbeta( Amu , Bmu )
    kappa ~ dgamma( Skappa , Rkappa )
    Amu <- 2.0
    Bmu <- 2.0
    Skappa <- pow(10,2)/pow(10,2)
    Rkappa <- 10/pow(10,2)
}
# ... JAGS model specification ends.
" # close quote to end modelString

# Write the modelString to a file, using R commands:
writeLines(modelString,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.

ncoins = 5 ; nflipspercoin = 50 
muAct = .7 ; kappaAct = 20
thetaAct = rbeta( ncoins ,muAct*kappaAct ,(1-muAct)*kappaAct )
z = rbinom( n=ncoins ,size=nflipspercoin ,prob=thetaAct )
N = rep( nflipspercoin , ncoins )

# Demo data for various figures in the book:
# N =  c( 5, 5, 5, 5, 5 ) # c( 10, 10, 10 )  # c( 15, 5 ) 
# z =  c( 1, 1, 1, 1, 5 ) # c(  1,  5,  9 )  # c(  3, 4 ) 
 
# Therapeutic touch data:
#  z = c(1,2,3,3,3,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,5,5,6,6,7,7,7,8)
#  N = rep(10,length(z))
# Convert z,N to vectors of individual flips.
# coin vector is index of coin for each flip.
# y vector is head or tail for each flip.
# For example,
#  coin = c( 1, 1, 2, 2, 2 )
#  y    = c( 1, 0, 0, 0, 1 )
# means that the first flip was of coin 1 and it was a head, the second flip
# was of coin 1 and it was a tail, the third flip was of coin 2 and it was a
# tail, etc.

coin = NULL ; y = NULL
for ( coinIdx in 1:length(N) ) {
    coin = c( coin , rep(coinIdx,N[coinIdx]) )
    y = c( y , rep(1,z[coinIdx]) , rep(0,N[coinIdx]-z[coinIdx]) )
}
nTrialTotal = length( y )
nCoins = length( unique( coin ) )

dataList = list(
    y = y ,
    coin = coin ,
    nTrialTotal = nTrialTotal ,
    nCoins = nCoins
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAIN.

#------------------------------------------------------------------------------
# RUN THE CHAINS.

parameters = c( "mu" , "kappa" , "theta" )   # The parameter(s) to be monitored.
adaptSteps = 1000              # Number of steps to "tune" the samplers.
burnInSteps = 5000            # Number of steps to "burn-in" the samplers.
nChains = 3                   # Number of chains to run.
numSavedSteps=50000           # Total number of steps in chains to save.
thinSteps=1                   # Number of steps to "thin" (1=keep every step).
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
jagsModel = jags.model( "model.txt" , data=dataList , # inits=initsList , 
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( jagsModel , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                            n.iter=nIter , thin=thinSteps )
# resulting codaSamples object has these indices: 
#   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS.

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

# Extract the posterior sample from JAGS:
muSample = mcmcChain[,"mu"]
kappaSample = mcmcChain[,"kappa"] # BRugs gets sample from JAGS
thetaSample = matrix( 0 , nrow=nCoins , ncol=nChains*nIter )
for ( coinIdx in 1:nCoins ) {
    nodeName = paste( "theta[" , coinIdx , "]" , sep="" )
    thetaSample[coinIdx,] = mcmcChain[,nodeName]
}

# Make a graph using R commands:
source("plotPost.R")
if ( nCoins <= 5 ) { # Only make this figure if there are not too many coins
  openGraph(width=3.2*3,height=2.5*(1+nCoins))
  layout( matrix( 1:(3*(nCoins+1)) , nrow=(nCoins+1) , byrow=T ) )
  par(mar=c(2.95,2.95,1.0,0),mgp=c(1.35,0.35,0),oma=c( 0.1, 0.1, 0.1, 0.1))
  nPtsToPlot = 500
  plotIdx = floor(seq(1,length(muSample),length=nPtsToPlot))
  kPltLim = signif( quantile( kappaSample , p=c(.01,.99) ) , 4 )
  plot( muSample[plotIdx] , kappaSample[plotIdx] , type="p" , ylim=kPltLim ,
        xlim=c(0,1) , xlab=expression(mu) , ylab=expression(kappa) , cex.lab=1.5 , 
        col="skyblue" )
  plotPost( muSample , xlab=bquote(mu) , xlim=c(0,1) , main="")
  plotPost( kappaSample , xlab=bquote(kappa) , main="" , showMode=TRUE ) 
  for ( coinIdx in 1:nCoins ) {
      plotPost( thetaSample[coinIdx,] , xlab=bquote(theta[.(coinIdx)]) ,
                xlim=c(0,1) , main="" )
      plot( thetaSample[coinIdx,plotIdx] , muSample[plotIdx] , type="p" ,
            xlim=c(0,1) , ylim=c(0,1) , cex.lab=1.5 ,
            xlab=bquote(theta[.(coinIdx)]) , ylab=expression(mu) , col="skyblue")
      plot( thetaSample[coinIdx,plotIdx] , kappaSample[plotIdx] , type="p" ,
            xlim=c(0,1) , ylim=kPltLim , cex.lab=1.5 ,
            xlab=bquote(theta[.(coinIdx)]) , ylab=expression(kappa) , col="skyblue")
  }
} # end if ( nCoins <= ...

#
openGraph(width=7,height=5)
layout( matrix( 1:4 , nrow=2 , byrow=T ) )
par(mar=c(2.95,2.95,1.0,0),mgp=c(1.35,0.35,0),oma=c( 0.1, 0.1, 0.1, 0.1) )
plotPost( muSample , xlab=bquote(mu) , main="" , compVal=0.5 )
plotPost( kappaSample , xlab=bquote(kappa) , main="" ,  HDItextPlace=.1 ,
          showMode=TRUE)
plotPost( thetaSample[1,] , xlab=bquote(theta[1]) , main="" , compVal=0.5 )
lastIdx = length(z)
plotPost( thetaSample[lastIdx,] , xlab=bquote(theta[.(lastIdx)]) , 
          main="" , compVal=0.5 )
