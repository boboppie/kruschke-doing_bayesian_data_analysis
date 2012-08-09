graphics.off()
rm(list=ls(all=TRUE))
fileNameRoot="BernBetaModelCompJags" # for constructing output filenames
if ( .Platform$OS.type != "windows" ) { 
  windows <- function( ... ) X11( ... ) 
}
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
model {
    # Likelihood:
    for ( i in 1:nflips ) {
        y[i] ~ dbern( theta )  # y[i] distributed as Bernoulli
    }
    # Prior distribution:
    theta ~ dbeta( aTheta , bTheta ) # theta distributed as beta density
    aTheta <- muTheta * kappaTheta
    bTheta <- (1-muTheta) * kappaTheta
    # Hyperprior:
    muTheta <- muThetaModel[ modelIndex ]
    muThetaModel[1] <- .75
    muThetaModel[2] <- .25
    kappaTheta <- 12
    # Hyperhyperprior:
    modelIndex ~ dcat( modelProb[] )
    modelProb[1] <- .5
    modelProb[2] <- .5
}
" # close quote for modelstring
# Write model to a file:
writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.

# Specify the data in a form that is compatible with JAGS model, as a list:
y = c( rep(0,3) , rep(1,6) )
nflips = length( y )
dataList = list(
    nflips = nflips ,
    y = y
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

# Specific initialization not necessary in this case, but here it is if wanted:
# initsList = list( theta=0.5 , modelIndex=1 ) 

#------------------------------------------------------------------------------
# RUN THE CHAINS.

parameters = c("theta","modelIndex") 
adaptSteps = 1000             # Number of steps to "tune" the samplers.
burnInSteps = 10000           # Number of steps to "burn-in" the samplers.
nChains = 1                   # Number of chains to run.
numSavedSteps=100000          # Total number of steps in chains to save.
thinSteps=1                   # Number of steps to "thin" (1=keep every step).
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
  plot( codaSamples , ask=T )  # not very useful for this application
  autocorr.plot( codaSamples , ask=T )
}

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )

# Get the posterior sample of modelIndex:
modelIdxSample = mcmcChain[,"modelIndex"]
# Compute the proportion of modelIndex at each value:
pM1 = sum( modelIdxSample == 1 ) / length( modelIdxSample )
pM2 = 1 - pM1

# Get the posterior sample of theta:
thetaSample = mcmcChain[,"theta"]
# Extract theta values when modelIndex is 1:
thetaSampleM1 = thetaSample[ modelIdxSample == 1 ]
# Extract theta values when modelIndex is 2:
thetaSampleM2 = thetaSample[ modelIdxSample == 2 ]

# Plot histograms of sampled theta values for each model,
# with pM displayed.
windows()
layout( matrix(1:2,nrow=2) )
hist( thetaSampleM1 , main="Posterior Theta_1 when Model Index = 1" ,
      xlab=expression(theta) , xlim=c(0,1) ,
      col="grey" , border="white" )
text( 0 , 0 , bquote( "p(M1|D)" == .(signif(pM1,3)) ) , adj=c(0,-2) , cex=1.5 )
hist( thetaSampleM2 , main="Posterior Theta_2 when Model Index = 2" ,
      xlab=expression(theta) , xlim=c(0,1) ,
      col="grey" , border="white" )
text( 0 , 0 , bquote( "p(M2|D)" == .(signif(pM2,3)) ) , adj=c(0,-2) , cex=1.5 )

savePlot( file=paste(fileNameRoot,".eps",sep="") , type="eps" )
savePlot( file=paste(fileNameRoot,".jpg",sep="") , type="jpg" )
