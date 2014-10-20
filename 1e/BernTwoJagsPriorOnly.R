rm(list = ls())
graphics.off()
source("openGraphSaveGraph.R")
source("plotPost.R")
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelString = "
# JAGS model specification begins here...
model {
    # Likelihood. Each flip is Bernoulli. 
    for ( i in 1 : N1 ) { y1[i] ~ dbern( theta1 ) }
    for ( i in 1 : N2 ) { y2[i] ~ dbern( theta2 ) }
    # Prior. Independent beta distributions.
    theta1 ~ dbeta( 3 , 3 )
    theta2 ~ dbeta( 3 , 3 )
}
# ... end JAGS model specification
" # close quote for modelstring

# Write the modelString to a file, using R commands:
writeLines(modelString,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.

# Specify the data in a form that is compatible with JAGS model, as a list:
dataList = list(
    N1 = 7 ,
    #y1 = c( 1,1,1,1,1,0,0 ) ,
    N2 = 7 #,
    #y2 = c( 1,1,0,0,0,0,0 )
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAIN.

# Can be done automatically in jags.model() by commenting out inits argument.
# Otherwise could be established as:
# initsList = list( theta1 = sum(dataList$y1)/length(dataList$y1) , 
#                   theta2 = sum(dataList$y2)/length(dataList$y2) )

#------------------------------------------------------------------------------
# RUN THE CHAINS.

parameters = c( "theta1" , "theta2" )     # The parameter(s) to be monitored.
adaptSteps = 500              # Number of steps to "tune" the samplers.
burnInSteps = 1000            # Number of steps to "burn-in" the samplers.
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

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )

theta1Sample = mcmcChain[,"theta1"] # Put sampled values in a vector.
theta2Sample = mcmcChain[,"theta2"] # Put sampled values in a vector.

# Plot the trajectory of the last 500 sampled values.
openGraph(width=7,height=7)
par( pty="s" )
chainlength=NROW(mcmcChain)
plot( theta1Sample[(chainlength-500):chainlength] ,
      theta2Sample[(chainlength-500):chainlength] , type = "o" ,
      xlim = c(0,1) , xlab = bquote(theta[1]) , ylim = c(0,1) ,
      ylab = bquote(theta[2]) , main="JAGS Result" , col="skyblue" )
# Display means in plot.
theta1mean = mean(theta1Sample)
theta2mean = mean(theta2Sample)
if (theta1mean > .5) { xpos = 0.0 ; xadj = 0.0
} else { xpos = 1.0 ; xadj = 1.0 }
if (theta2mean > .5) { ypos = 0.0 ; yadj = 0.0
} else { ypos = 1.0 ; yadj = 1.0 }
text( xpos , ypos ,
	bquote(
	"M=" * .(signif(theta1mean,3)) * "," * .(signif(theta2mean,3))
	) ,adj=c(xadj,yadj) ,cex=1.5  )
#saveGraph(file="BernTwoJagsPriorOnly",type="eps")

# Plot a histogram of the posterior differences of theta values.
thetaDiff = theta1Sample - theta2Sample

openGraph(width=7,height=4)
plotPost( thetaDiff , xlab=expression(theta[1]-theta[2]) , compVal=0.0 )
#saveGraph(file="BernTwoJagsPriorOnlyDiff",type="eps")

# For Exercise 8.5:
# Posterior prediction. For each step in the chain, use the posterior thetas 
# to flip the coins.
chainLength = length( theta1Sample )
# Create matrix to hold results of simulated flips:
yPred = matrix( NA , nrow=2 , ncol=chainLength ) 
for ( stepIdx in 1:chainLength ) { # step through the chain
  # flip the first coin:
  pHead1 = theta1Sample[stepIdx]
  yPred[1,stepIdx] = sample( x=c(0,1), prob=c(1-pHead1,pHead1), size=1 )
  # flip the second coin:
  pHead2 = theta2Sample[stepIdx]
  yPred[2,stepIdx] = sample( x=c(0,1), prob=c(1-pHead2,pHead2), size=1 )
}
# Now determine the proportion of times that y1==1 and y2==0
pY1eq1andY2eq0 = sum( yPred[1,]==1 & yPred[2,]==0 ) / chainLength
