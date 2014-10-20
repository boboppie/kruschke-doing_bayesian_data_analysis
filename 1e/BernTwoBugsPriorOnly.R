library(BRugs)         # Kruschke, J. K. (2010). Doing Bayesian data analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
# BUGS model specification begins here...
model {
    # Likelihood. Each flip is Bernoulli. 
    for ( i in 1 : N1 ) { y1[i] ~ dbern( theta1 ) }
    for ( i in 1 : N2 ) { y2[i] ~ dbern( theta2 ) }
    # Prior. Independent beta distributions.
    theta1 ~ dbeta( 3 , 3 )
    theta2 ~ dbeta( 3 , 3 )
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
datalist = list(
    N1 = 7 ,
#    y1 = c( 1,1,1,1,1,0,0 ) ,
    N2 = 7 #,
#    y2 = c( 1,1,0,0,0,0,0 )
)
# Get the data into BRugs:
modelData( bugsData( datalist ) )  # commented out

#------------------------------------------------------------------------------
# INTIALIZE THE CHAIN.

modelCompile()
modelGenInits()

#------------------------------------------------------------------------------
# RUN THE CHAINS.

samplesSet( c( "theta1" , "theta2" ) ) # Keep a record of sampled "theta" values
chainlength = 10000                     # Arbitrary length of chain to generate.
modelUpdate( chainlength )             # Actually generate the chain.

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS.

theta1Sample = samplesSample( "theta1" ) # Put sampled values in a vector.
theta2Sample = samplesSample( "theta2" ) # Put sampled values in a vector.

# Plot the trajectory of the last 500 sampled values.
windows()
par( pty="s" )
plot( theta1Sample[(chainlength-500):chainlength] ,
      theta2Sample[(chainlength-500):chainlength] , type = "o" ,
      xlim = c(0,1) , xlab = bquote(theta[1]) , ylim = c(0,1) ,
      ylab = bquote(theta[2]) , main="BUGS Result" )
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
dev.copy2eps(file="BernTwoBugsPriorOnly.eps")

# Plot a histogram of the posterior differences of theta values.
thetaDiff = theta1Sample - theta2Sample
windows(7,4)
source("plotPost.R")
plotPost( thetaDiff , xlab=expression(theta[1]-theta[2]) ,
          breaks=20 , main="" )
dev.copy2eps(file="BernTwoBugsPriorOnlyDiff.eps")