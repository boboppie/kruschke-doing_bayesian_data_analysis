library(BRugs)         # Kruschke, J. K. (2010). Doing Bayesian data analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
# BUGS model specification begins here...
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
# ... end BUGS model specification
" # close quote for modelstring
# Write model to a file:
.temp = file("model.txt","w") ; writeLines(modelstring,con=.temp) ; close(.temp) 
# Load model file into BRugs and check its syntax:
modelCheck( "model.txt" )

#------------------------------------------------------------------------------
# THE DATA.

# Specify the data in a form that is compatible with BRugs model, as a list:
y = c( rep(0,3) , rep(1,6) )
nflips = length( y )
datalist = list(
    nflips = nflips ,
    y = y
)

# Get the data into BRugs:
# Function bugsData stores the data file (default filename is data.txt).
# Function modelData loads data file into BRugs (default filename is data.txt).
modelData( bugsData( datalist ) )

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

modelCompile( numChains=1 )
modelGenInits()

#------------------------------------------------------------------------------
# RUN THE CHAINS.

burninSteps = 10000
modelUpdate( burninSteps )
samplesSet( c("theta","modelIndex") )
nPerChain = 100000
modelUpdate( nPerChain , thin=1 ) # takes nPerChain * thin steps

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS.

# Get the posterior sample of modelIndex:
modelIdxSample = samplesSample( "modelIndex" )
# Compute the proportion of modelIndex at each value:
pM1 = sum( modelIdxSample == 1 ) / length( modelIdxSample )
pM2 = 1 - pM1

# Get the posterior sample of theta:
thetaSample = samplesSample( "theta" )
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

dev.copy2eps(file="BernBetaModelCompBrugs.eps")