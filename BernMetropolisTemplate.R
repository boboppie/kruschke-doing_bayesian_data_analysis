# Use this program as a template for experimenting with the Metropolis
# algorithm applied to a single parameter called theta, defined on the 
# interval [0,1].

# Specify the data, to be used in the likelihood function.
# This is a vector with one component per flip,
# in which 1 means a "head" and 0 means a "tail".
myData = c( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0 )

# Define the Bernoulli likelihood function, p(D|theta).
# The argument theta could be a vector, not just a scalar.
likelihood = function( theta , data ) {
  z = sum( data == 1 )
  N = length( data )
  pDataGivenTheta = theta^z * (1-theta)^(N-z)
  # The theta values passed into this function are generated at random,
  # and therefore might be inadvertently greater than 1 or less than 0.
  # The likelihood for theta > 1 or for theta < 0 is zero:
  pDataGivenTheta[ theta > 1 | theta < 0 ] = 0
  return( pDataGivenTheta )
}

# Define the prior density function. For purposes of computing p(D),
# at the end of this program, we want this prior to be a proper density.
# The argument theta could be a vector, not just a scalar.
prior = function( theta ) {
  prior = rep( 1 , length(theta) ) # uniform density over [0,1]
   # For kicks, here's a bimodal prior. To try it, uncomment the next line.
   #prior = dbeta( pmin(2*theta,2*(1-theta)) ,2,2 )
  # The theta values passed into this function are generated at random,
  # and therefore might be inadvertently greater than 1 or less than 0.
  # The prior for theta > 1 or for theta < 0 is zero:
  prior[ theta > 1 | theta < 0 ] = 0
  return( prior )
}

# Define the relative probability of the target distribution, 
# as a function of vector theta. For our application, this
# target distribution is the unnormalized posterior distribution.
targetRelProb = function( theta , data ) {
  targetRelProb =  likelihood( theta , data ) * prior( theta )
  return( targetRelProb )
}

# Specify the length of the trajectory, i.e., the number of jumps to try:
trajLength = 55556 # arbitrary large number
# Initialize the vector that will store the results:
trajectory = rep( 0 , trajLength )
# Specify where to start the trajectory:
trajectory[1] = 0.50 # arbitrary value
# Specify the burn-in period:
burnIn = ceiling( 0.1 * trajLength ) # arbitrary number, less than trajLength
# Initialize accepted, rejected counters, just to monitor performance:
nAccepted = 0
nRejected = 0
# Specify seed to reproduce same random walk:
set.seed(47405)

# Now generate the random walk. The 't' index is time or trial in the walk.
for ( t in 1:(trajLength-1) ) {
	currentPosition = trajectory[t]
	# Use the proposal distribution to generate a proposed jump.
	# The shape and variance of the proposal distribution can be changed
	# to whatever you think is appropriate for the target distribution.
	proposedJump = rnorm( 1 , mean = 0 , sd = 0.1 )
	# Compute the probability of accepting the proposed jump.
	probAccept = min( 1,
		targetRelProb( currentPosition + proposedJump , myData )
		/ targetRelProb( currentPosition , myData ) )
	# Generate a random uniform value from the interval [0,1] to
	# decide whether or not to accept the proposed jump.
	if ( runif(1) < probAccept ) {
		# accept the proposed jump
		trajectory[ t+1 ] = currentPosition + proposedJump
		# increment the accepted counter, just to monitor performance
		if ( t > burnIn ) { nAccepted = nAccepted + 1 }
	} else {
		# reject the proposed jump, stay at current position
		trajectory[ t+1 ] = currentPosition
		# increment the rejected counter, just to monitor performance
		if ( t > burnIn ) { nRejected = nRejected + 1 }
	}
}

# Extract the post-burnIn portion of the trajectory.
acceptedTraj = trajectory[ (burnIn+1) : length(trajectory) ]

# End of Metropolis algorithm.

#-----------------------------------------------------------------------
# Display the posterior.

source("openGraphSaveGraph.R")
source("plotPost.R")
openGraph(width=5,height=5)
mcmcInfo = plotPost( acceptedTraj , xlim=c(0,1) , xlab=bquote(theta) )

# Display rejected/accepted ratio in the plot.
# Get the highest point and mean of the plot for subsequent text positioning:
densMax = max( density( acceptedTraj )$y )
meanTraj = mean( acceptedTraj )
sdTraj = sd( acceptedTraj )
if ( meanTraj > .5 ) {
  xpos = 0.0 ; xadj = 0.0
} else {
  xpos = 1.0 ; xadj = 1.0
}
text( xpos , 0.75*densMax ,
	bquote(	N[pro] * "=" * .(length(acceptedTraj)) * "  " *
	frac(N[acc],N[pro]) * "=" * .(signif( nAccepted/length(acceptedTraj) , 3 ))
	) , adj=c(xadj,0)  )

#------------------------------------------------------------------------
# Evidence for model, p(D).

# Compute a,b parameters for beta distribution that has the same mean
# and stdev as the sample from the posterior. This is a useful choice
# when the likelihood function is Bernoulli.
a =   meanTraj   * ( (meanTraj*(1-meanTraj)/sdTraj^2) - 1 )
b = (1-meanTraj) * ( (meanTraj*(1-meanTraj)/sdTraj^2) - 1 )

# For every theta value in the posterior sample, compute 
# dbeta(theta,a,b) / likelihood(theta)*prior(theta)
# This computation assumes that likelihood and prior are proper densities,
# i.e., not just relative probabilities. This computation also assumes that
# the likelihood and prior functions were defined to accept a vector argument,
# not just a single-component scalar argument.
wtdEvid = dbeta( acceptedTraj , a , b ) / (
		likelihood( acceptedTraj , myData ) * prior( acceptedTraj ) )
pData = 1 / mean( wtdEvid )

# Display p(D) in the graph
if ( meanTraj > .5 ) { xpos = 0.0 ; xadj = 0.0
} else { xpos = 1.0 ; xadj = 1.0 }
text( xpos , 0.9*densMax , bquote( p(D)==.( signif(pData,3) ) ) ,
      adj=c(xadj,0) , cex=1.5 )

# Uncomment next line if you want to save the graph.
saveGraph(file="BernMetropolisTemplate",type="eps")