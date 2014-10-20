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
# Data = c(rep(1,3),rep(0,9)) or Data = c(rep(1,1),rep(0,11)).
Data = c(rep(1,3),rep(0,9))
nHeads = sum( Data )
nTails = length( Data ) - nHeads

# Compute the likelihood of the data for each value of theta:
pDataGivenTheta = Theta^nHeads * (1-Theta)^nTails

# Compute the posterior:
pData = sum( pDataGivenTheta * pTheta )
pThetaGivenData = pDataGivenTheta * pTheta / pData   # This is Bayes' rule!

# Plot the results.
source("openGraphSaveGraph.R") # read in graph functions
openGraph(width=7,height=10,mag=0.7) # open a window for the graph
layout( matrix( c( 1,2,3 ) ,nrow=3 ,ncol=1 ,byrow=FALSE ) ) # 3x1 panels
par(mar=c(3,3,1,0))         # number of margin lines: bottom,left,top,right
par(mgp=c(2,1,0))           # which margin lines to use for labels
par(mai=c(0.5,0.5,0.3,0.1)) # margin size in inches: bottom,left,top,right

# Plot the prior:
plot( Theta , pTheta , type="h" , lwd=3 , main="Prior" ,
      xlim=c(0,1) , xlab=bquote(theta) ,
      ylim=c(0,1.1*max(pThetaGivenData)) , ylab=bquote(p(theta)) ,
      cex.axis=1.2 , cex.lab=1.5 , cex.main=1.5 , col="skyblue" )

# Plot the likelihood:
plot( Theta , pDataGivenTheta , type="h" , lwd=3 , main="Likelihood" ,
      xlim=c(0,1) , xlab=bquote(theta) ,
      ylim=c(0,1.1*max(pDataGivenTheta)) , ylab=bquote(paste("p(D|",theta,")")),
      cex.axis=1.2 , cex.lab=1.5 , cex.main=1.5 , col="skyblue" )
text( .55 , .85*max(pDataGivenTheta) , cex=2.0 ,
      bquote( "D=" * .(nHeads) * "H," * .(nTails) * "T" ) , adj=c(0,.5) )

# Plot the posterior:
plot( Theta , pThetaGivenData , type="h" , lwd=3 , main="Posterior" ,
      xlim=c(0,1) , xlab=bquote(theta) ,
      ylim=c(0,1.1*max(pThetaGivenData)) , ylab=bquote(paste("p(",theta,"|D)")),
      cex.axis=1.2 , cex.lab=1.5 , cex.main=1.5 , col="skyblue" )
text( .55 , .85*max(pThetaGivenData) , cex=2.0 ,
      bquote( "p(D)=" * .(signif(pData,3)) ) , adj=c(0,.5) )

# Save the plot 
if ( nThetaVals == 3 ) { modeltype = "simpleModel" }
if ( nThetaVals == 63 ) { modeltype = "complexModel" }
if ( nHeads == 3 & nTails == 9 ) { datatype = "simpleData" }
if ( nHeads == 1 & nTails == 11 ) { datatype = "complexData" }
filename = paste( "BayesUpdate_" ,modeltype ,"_" ,datatype,sep="" )
saveGraph( file=filename , type="eps" )
