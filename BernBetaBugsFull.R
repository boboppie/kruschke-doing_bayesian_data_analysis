library(BRugs)         # Kruschke, J. K. (2010). Doing Bayesian data analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

# Specify the model in BUGS language, but save it as a string in R:
modelString = "
# BUGS model specification begins ...
model {
    # Likelihood:
    for ( i in 1:nFlips ) {
        y[i] ~ dbern( theta )
    }
    # Prior distribution:
    theta ~ dbeta( priorA , priorB )
    priorA <- 1
    priorB <- 1
}
# ... BUGS model specification ends.
" # close quote to end modelString

# Write the modelString to a file, using R commands:
writeLines(modelString,con="model.txt")
# Use BRugs to send the model.txt file to BUGS, which checks the model syntax:
modelCheck( "model.txt" )

#------------------------------------------------------------------------------
# THE DATA.

# Specify the data in R, using a list format compatible with BUGS:
dataList = list(
    nFlips = 14 ,
    y = c( 1,1,1,1,1,1,1,1,1,1,1,0,0,0 )
)

# Use BRugs commands to put the data into a file and ship the file to BUGS:
modelData( bugsData( dataList ) )

#------------------------------------------------------------------------------
# INTIALIZE THE CHAIN.

modelCompile()  # BRugs command tells BUGS to compile the model.
modelGenInits() # BRugs command tells BUGS to randomly initialize a chain.

#------------------------------------------------------------------------------
# RUN THE CHAINS.

# BRugs tells BUGS to keep a record of the sampled "theta" values:
samplesSet( "theta" )
# R command defines a new variable that specifies an arbitrary chain length:
chainLength = 10000
# BRugs tells BUGS to generate a MCMC chain:
modelUpdate( chainLength )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS.

thetaSample = samplesSample( "theta" ) # BRugs asks BUGS for the sample values.
thetaSummary = samplesStats( "theta" ) # BRugs asks BUGS for summary statistics.

# Make a graph using R commands:
windows(10,6)
layout( matrix( c(1,2) , nrow=1 ) )
plot( thetaSample[1:500] , 1:length(thetaSample[1:500]) , type="o" ,
      xlim=c(0,1) , xlab=bquote(theta) , ylab="Position in Chain" ,
      cex.lab=1.25 , main="BUGS Results" )
source("plotPost.R")
histInfo = plotPost( thetaSample , xlim=c(0,1) )
dev.copy2eps(file="BernBetaBugsFull.eps")

# Posterior prediction:
# For each step in the chain, use posterior theta to flip a coin:
chainLength = length( thetaSample )
yPred = rep( NULL , chainLength )  # define placeholder for flip results
for ( stepIdx in 1:chainLength ) {
  pHead = thetaSample[stepIdx]
  yPred[stepIdx] = sample( x=c(0,1), prob=c(1-pHead,pHead), size=1 )
}
# Jitter the 0,1 y values for plotting purposes:
yPredJittered = yPred + runif( length(yPred) , -.05 , +.05 )
# Now plot the jittered values:
windows(5,5.5)
plot( thetaSample[1:500] , yPredJittered[1:500] , xlim=c(0,1) ,
      main="posterior predictive sample" ,
      xlab=expression(theta) , ylab="y (jittered)" )
points( mean(thetaSample) , mean(yPred) , pch="+" , cex=2 )
text( mean(thetaSample) , mean(yPred) ,
      bquote( mean(y) == .(signif(mean(yPred),2)) ) ,
      adj=c(1.2,.5) )
text( mean(thetaSample) , mean(yPred) , srt=90 ,
      bquote( mean(theta) == .(signif(mean(thetaSample),2)) ) ,
      adj=c(1.2,.5) )
abline( 0 , 1 , lty="dashed" , lwd=2 )
dev.copy2eps(file="BernBetaBugsPost.eps")