OrdinalProbitDataGenerator = function( nData , nYlevels=5 ,
                              normPrec=1 , slope=c(.5,.5) ,
                              thresh=c(-Inf,seq(-1.25,1.25,length=nYlevels-1),Inf) ,
                              makePlots = FALSE , rndSeed=NULL ) {

if ( !is.null(rndSeed) ) { set.seed( rndSeed ) }

# Generate random _standardized_ X values.
nPredictors = length(slope)
Xdata = matrix( rnorm( nPredictors*nData , 0 , 1 ) ,
                nrow=nData, ncol=nPredictors)

# Standardize the X values:
for ( colIdx in 1:NCOL(Xdata) ) {
    mX = mean( Xdata[,colIdx] )
    sdX = sd( Xdata[,colIdx] )
    Xdata[,colIdx] =  ( Xdata[,colIdx] - mX ) / sdX
}

# Generate continuous mu values as linear function of X.
bias = 0
slope = slope / sum(abs(slope))
dim(slope) = c(nPredictors,1) # make it a column vector
mu = bias + Xdata %*% slope

# Convert continuous mu values to discrete (ordinal) Y values:
# Utility function: Cumulative normal, parameterized with precision not SD.
cumnorm = function( x , prec=1 ) {
	y = pnorm( x , mean=0 , sd=(1/sqrt(prec)) )
	return( y )
}
# Randomly generate discrete Y values based on proximity to bin thresholds.
Ydata = matrix( 0 , nrow=nData , ncol=1 )
for ( subjIdx in 1:nData ) {
    Yprob = ( cumnorm( thresh[2:(nYlevels+1)] - mu[subjIdx] , normPrec )
            - cumnorm( thresh[1:(nYlevels)]   - mu[subjIdx] , normPrec ) )
    Ydata[subjIdx] = sample( 1:nYlevels , prob = Yprob ,
                             size=1 , replace=T )
}

# Combine X and Y into a data matrix.
datamatrix = cbind( Ydata , Xdata )
colnames( datamatrix ) = c( "Y" , paste("X",1:nPredictors,sep="") )

# Plot the data.
if ( makePlots ) {
  colorlist = rep( c("black","red","blue","green","gold","purple","cyan","brown"),
                   length=nYlevels )
  # 2D scatterplot
  xrange = range(datamatrix[,2])
  yrange = range(datamatrix[,3])
  rowIdx = ( datamatrix[,1] == 1 )
  plot( datamatrix[rowIdx,2] , datamatrix[rowIdx,3] , pch=as.character(1) , 
        col=colorlist[respIdx=1] ,
        main=paste("Ordinal Values (1-",nYlevels,")",sep="") ,
        #xlim=c(10000,90000) , ylim=c(0,10) ,
        xlim=xrange , ylim=yrange ,
        xlab="X1" , ylab="X2" )
  for ( respIdx in 2:nYlevels ) {
    rowIdx = ( datamatrix[,1] == respIdx )
    points( datamatrix[rowIdx,2] , datamatrix[rowIdx,3] , col=colorlist[respIdx] , 
		pch=as.character(respIdx) )
  }

} # end if makePlots

return( datamatrix )

} # end function   # Kruschke, J. K. (2011). Doing Bayesian data analysis: A
                   # Tutorial with R and BUGS. Academic Press / Elsevier.