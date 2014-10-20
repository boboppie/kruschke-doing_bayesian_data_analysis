graphics.off()          # Clears all graphical displays.
rm(list=ls(all=TRUE))   # Removes all variables from memory!

# Get the Bayesian functions into R's working memory:
source("BMLR.R")

# Load the data into R:
dataMat = read.csv( file="BMLRexampleData.csv" )
# Important: The matrix dataMat must have the criterion y in its first column, 
# and the predictors in the subsequent columns!

# Run the Bayesian analysis and put the results in variable named mcmcChain:
mcmcChain = BMLRmcmc( dataMat )

# Plot the posterior distribution and put summary in variable named postInfo:
postInfo = BMLRplot( mcmcChain ) 
# Display the summary on the console:
show(postInfo)


# Another example, using a large N data set:
dataMat = read.csv( file="BMLRexampleLargeNdata.csv" )
mcmcChain = BMLRmcmc( dataMat , numSavedSteps=100000 )
postInfo = BMLRplot( mcmcChain , 
                     ROPEbeta=c(-0.05,0.05) , ROPEbetaDiff=c(-0.05,0.05) ) 
show(postInfo)
