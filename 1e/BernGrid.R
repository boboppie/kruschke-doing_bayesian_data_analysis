BernGrid = function( Theta , pTheta , Data ,
                     credib=.95 , nToPlot=length(Theta) ) {
# Bayesian updating for Bernoulli likelihood and prior specified on a grid.
# Input arguments:
#  Theta is a vector of theta values, all between 0 and 1.
#  pTheta is a vector of corresponding probability _masses_.
#  Data is a vector of 1's and 0's, where 1 corresponds to a and 0 to b.
#  credib is the probability mass of the credible interval, default is 0.95.
#  nToPlot is the number of grid points to plot; defaults to all of them.
# Output:
#  pThetaGivenData is a vector of posterior probability masses over Theta.
#  Also creates a three-panel graph of prior, likelihood, and posterior 
#  probability masses with credible interval.
# Example of use:
#  # Create vector of theta values.
#  > binwidth = 1/1000
#  > thetagrid = seq( from=binwidth/2 , to=1-binwidth/2 , by=binwidth )
#  # Specify probability mass at each theta value.
#  > relprob = pmin(thetagrid,1-thetagrid) # relative prob at each theta
#  > prior = relprob / sum(relprob) # probability mass at each theta
#  # Specify the data vector.
#  > datavec = c( rep(1,3) , rep(0,1) ) # 3 heads, 1 tail
#  # Call the function.
#  > posterior = BernGrid( Theta=thetagrid , pTheta=prior , Data=datavec )
# Hints:
#  You will need to "source" this function before calling it.
#  You may want to define a tall narrow window before using it; e.g.,
#  > source("openGraphSaveGraph.R")
#  > openGraph(width=7,height=10,mag=0.7)

# Create summary values of Data
z = sum( Data==1 ) # number of 1's in Data
N = length( Data ) # number of flips in Data
# Compute the likelihood of the Data for each value of Theta.
pDataGivenTheta = Theta^z * (1-Theta)^(N-z)
# Compute the evidence and the posterior.
pData = sum( pDataGivenTheta * pTheta )
pThetaGivenData = pDataGivenTheta * pTheta / pData

# Plot the results.
layout( matrix( c( 1,2,3 ) ,nrow=3 ,ncol=1 ,byrow=FALSE ) ) # 3x1 panels
par( mar=c(3,3,1,0) , mgp=c(2,1,0) , mai=c(0.5,0.5,0.3,0.1) ) # margin settings
dotsize = 5 # how big to make the plotted dots
# If the comb has a zillion teeth, it's too many to plot, so plot only a
# thinned out subset of the teeth.
nteeth = length(Theta)
if ( nteeth > nToPlot ) {
  thinIdx = seq( 1, nteeth , round( nteeth / nToPlot ) )
  if ( length(thinIdx) < length(Theta) ) {
    thinIdx = c( thinIdx , nteeth ) # makes sure last tooth is included
  }
} else { thinIdx = 1:nteeth }
# Plot the prior.
meanTheta = sum( Theta * pTheta ) # mean of prior, for plotting
plot( Theta[thinIdx] , pTheta[thinIdx] , type="p" , pch="." , cex=dotsize ,
      xlim=c(0,1) , ylim=c(0,1.1*max(pThetaGivenData)) , cex.axis=1.2 ,
      xlab=bquote(theta) , ylab=bquote(p(theta)) , cex.lab=1.5 ,
      main="Prior" , cex.main=1.5 , col="skyblue" )
if ( meanTheta > .5 ) {
   textx = 0 ; textadj = c(0,1)
} else {
  textx = 1 ; textadj = c(1,1)
}
text( textx , 1.0*max(pThetaGivenData) ,
      bquote( "mean(" * theta * ")=" * .(signif(meanTheta,3)) ) ,
      cex=2.0 , adj=textadj )
# Plot the likelihood: p(Data|Theta)
plot(Theta[thinIdx] ,pDataGivenTheta[thinIdx] ,type="p" ,pch="." ,cex=dotsize
	,xlim=c(0,1) ,cex.axis=1.2 ,xlab=bquote(theta) 
	,ylim=c(0,1.1*max(pDataGivenTheta)) 
	,ylab=bquote( "p(D|" * theta * ")" )  
	,cex.lab=1.5 ,main="Likelihood" ,cex.main=1.5 , col="skyblue" )
if ( z > .5*N ) { textx = 0 ; textadj = c(0,1) }
else { textx = 1 ; textadj = c(1,1) }
text( textx ,1.0*max(pDataGivenTheta) ,cex=2.0
	,bquote( "Data: z=" * .(z) * ",N=" * .(N) ) ,adj=textadj )
# Plot the posterior.
meanThetaGivenData = sum( Theta * pThetaGivenData )
plot(Theta[thinIdx] ,pThetaGivenData[thinIdx] ,type="p" ,pch="." ,cex=dotsize
	,xlim=c(0,1) ,ylim=c(0,1.1*max(pThetaGivenData)) ,cex.axis=1.2 
	,xlab=bquote(theta) ,ylab=bquote( "p(" * theta * "|D)" )
	,cex.lab=1.5 ,main="Posterior" ,cex.main=1.5 , col="skyblue" )
if ( meanThetaGivenData > .5 ) { textx = 0 ; textadj = c(0,1) } 
else { textx = 1 ; textadj = c(1,1) }
text(textx ,1.00*max(pThetaGivenData) ,cex=2.0
	,bquote( "mean(" * theta * "|D)=" * .(signif(meanThetaGivenData,3)) ) 
	,adj=textadj )
text(textx ,0.75*max(pThetaGivenData) ,cex=2.0
	,bquote( "p(D)=" * .(signif(pData,3)) ) ,adj=textadj )
# Mark the highest density interval. HDI points are not thinned in the plot.
source("HDIofGrid.R")
HDIinfo = HDIofGrid( pThetaGivenData , credMass=credib )
points( Theta[ HDIinfo$indices ] ,
       rep( HDIinfo$height , length( HDIinfo$indices ) ) , pch="-" , cex=1.0 )
text( mean( Theta[ HDIinfo$indices ] ) , HDIinfo$height ,
         bquote( .(100*signif(HDIinfo$mass,3)) * "% HDI" ) ,
         adj=c(0.5,-1.5) , cex=1.5 )
# Mark the left and right ends of the waterline. 
# Find indices at ends of sub-intervals:
inLim = HDIinfo$indices[1] # first point
for ( idx in 2:(length(HDIinfo$indices)-1) ) {
  if ( ( HDIinfo$indices[idx] != HDIinfo$indices[idx-1]+1 ) | # jumps on left, OR
    ( HDIinfo$indices[idx] != HDIinfo$indices[idx+1]-1 ) ) { # jumps on right
    inLim = c(inLim,HDIinfo$indices[idx]) # include idx
  }
}
inLim = c(inLim,HDIinfo$indices[length(HDIinfo$indices)]) # last point
# Mark vertical lines at ends of sub-intervals:
for ( idx in inLim ) {
  lines( c(Theta[idx],Theta[idx]) , c(-0.5,HDIinfo$height) , type="l" , lty=2 , 
         lwd=1.5 )
  text( Theta[idx] , HDIinfo$height , bquote(.(round(Theta[idx],3))) ,
        adj=c(0.5,-0.1) , cex=1.2 )
}

return( pThetaGivenData )
} # end of function