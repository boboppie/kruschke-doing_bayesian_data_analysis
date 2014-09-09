
# OPTIONAL: Clear R's memory and graphics:
rm(list=ls())  # Careful! This clears all of R's memory!
graphics.off() # This closes all of R's graphics windows.

# Specify the data
y = c(101,100,102,104,102,97,105,105,98,101,100,123,105,103,100,95,102,106,
       109,102,82,102,100,102,102,101,102,102,103,103,97,97,103,101,97,104,
       96,103,124,101,101,100,101,101,104,100,101)

# Run the Bayesian analysis:
source("BEST1G.R")
mcmcChain = BEST1Gmcmc( y )

# Display the results:
BEST1Gplot( y , mcmcChain , compValm=100 , ROPEeff=c(-0.1,0.1) , pairsPlot=TRUE )
