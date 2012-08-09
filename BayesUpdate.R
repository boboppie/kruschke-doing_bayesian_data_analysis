# Theta is the vector of candidate values for the parameter theta.
# nThetaVals is the number of candidate theta values.
# To produce the examples in the book, set nThetaVals to either 3 or 63.
nThetaVals = 3
# Now make the vector of theta values:
Theta = seq( from = 1/(nThetaVals+1) , to = nThetaVals/(nThetaVals+1) ,
             by = 1/(nThetaVals+1) )

# pTheta is the vector of prior probabilities on the theta values.
pTheta = pmin( Theta , 1-Theta ) # Makes a triangular belief distribution.
pTheta = pTheta / sum( pTheta )  # Makes sure that beliefs sum to 1.

# Specify the data. To produce the examples in the book, use either
# Data = c(1,1,1,0,0,0,0,0,0,0,0,0) or Data = c(1,0,0,0,0,0,0,0,0,0,0,0).
Data = c(1,1,1,0,0,0,0,0,0,0,0,0)
nHeads = sum( Data == 1 )
nTails = sum( Data == 0 )

# Compute the likelihood of the data for each value of theta:
pDataGivenTheta = Theta^nHeads * (1-Theta)^nTails

# Compute the posterior:
pData = sum( pDataGivenTheta * pTheta )
pThetaGivenData = pDataGivenTheta * pTheta / pData   # This is Bayes' rule!

# Plot the results.
windows(7,10) # create window of specified width,height inches.
layout( matrix( c( 1,2,3 ) ,nrow=3 ,ncol=1 ,byrow=FALSE ) ) # 3x1 panels
par(mar=c(3,3,1,0))         # number of margin lines: bottom,left,top,right
par(mgp=c(2,1,0))           # which margin lines to use for labels
par(mai=c(0.5,0.5,0.3,0.1)) # margin size in inches: bottom,left,top,right

# Plot the prior:
plot( Theta , pTheta , type="h" , lwd=3 , main="Prior" ,
      xlim=c(0,1) , xlab=bquote(theta) ,
      ylim=c(0,1.1*max(pThetaGivenData)) , ylab=bquote(p(theta)) ,
      cex.axis=1.2 , cex.lab=1.5 , cex.main=1.5 )

# Plot the likelihood:
plot( Theta , pDataGivenTheta , type="h" , lwd=3 , main="Likelihood" ,
      xlim=c(0,1) , xlab=bquote(theta) ,
      ylim=c(0,1.1*max(pDataGivenTheta)) , ylab=bquote(paste("p(D|",theta,")")),
      cex.axis=1.2 , cex.lab=1.5 , cex.main=1.5 )
text( .55 , .85*max(pDataGivenTheta) , cex=2.0 ,
      bquote( "D=" * .(nHeads) * "H," * .(nTails) * "T" ) , adj=c(0,.5) )

# Plot the posterior:
plot( Theta , pThetaGivenData , type="h" , lwd=3 , main="Posterior" ,
      xlim=c(0,1) , xlab=bquote(theta) ,
      ylim=c(0,1.1*max(pThetaGivenData)) , ylab=bquote(paste("p(",theta,"|D)")),
      cex.axis=1.2 , cex.lab=1.5 , cex.main=1.5 )
text( .55 , .85*max(pThetaGivenData) , cex=2.0 ,
      bquote( "p(D)=" * .(signif(pData,3)) ) , adj=c(0,.5) )

# Save the plot as an EPS file.
if ( nThetaVals == 3 ) { modeltype = "simpleModel" }
if ( nThetaVals == 63 ) { modeltype = "complexModel" }
if ( nHeads == 3 & nTails == 9 ) { datatype = "simpleData" }
if ( nHeads == 1 & nTails == 11 ) { datatype = "complexData" }
filename = paste( "BayesUpdate_" ,modeltype ,"_" ,datatype ,".eps" ,sep="" )
# The command dev.copy2eps, used below, doesn't work on all systems.
# Try help("dev.copy2eps") for info about saving graphs in other file formats.
dev.copy2eps( file=filename )
