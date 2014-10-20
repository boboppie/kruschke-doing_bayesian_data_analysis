# Specify known values of prior and actual data.
priorA = 100
priorB = 1
actualDataZ  = 8
actualDataN  = 12
# Compute posterior parameter values.
postA = priorA + actualDataZ
postB = priorB + actualDataN - actualDataZ
# Number of flips in a simulated sample should match the actual sample size:
simSampleSize = actualDataN
# Designate an arbitrarily large number of simulated samples.
nSimSamples = 10000
# Set aside a vector in which to store the simulation results.
simSampleZrecord = vector( length=nSimSamples )
# Now generate samples from the posterior.
for ( sampleIdx in 1:nSimSamples ) {
	# Generate a theta value for the new sample from the posterior.
	sampleTheta = rbeta( 1 , postA , postB )
	# Generate a sample, using sampleTheta.
	sampleData = sample( x=c(0,1) , prob=c( 1-sampleTheta , sampleTheta ) ,
                          size=simSampleSize , replace=TRUE )
	# Store the number of heads in sampleData.
	simSampleZrecord[ sampleIdx ] = sum( sampleData )
}
# Make a histogram of the number of heads in the samples.
hist( simSampleZrecord )
                   # Kruschke, J. K. (2011). Doing Bayesian data analysis: A
                   # Tutorial with R and BUGS. Academic Press / Elsevier.