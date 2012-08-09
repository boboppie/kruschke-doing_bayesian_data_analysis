graphics.off()
rm(list=ls(all=TRUE))
filenamebase = "NormalModelCompBrugs"
library(BRugs)  # John K. Kruschke, 2011. For details on how to use this
                # program, see similar programs in the book: Doing Bayesian
                # Data Analysis: A Tutorial with R and BUGS. ISBN 9780123814852 

#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
# BUGS model specification begins here...
model {
  # Likelihood:
  for ( i in 1:N ) {
    y[i] ~ dnorm( mu , tau )
  }
  # Prior:
  mu ~ dnorm( M[mIdx] , T[mIdx] )
  tau <- pow(sigma,-2)
  sigma ~ dunif( L[mIdx] , H[mIdx] )
  M[1] <- 0
  T[1] <- pow(nullPriorSD,-2)
  L[1] <- 0
  H[1] <- 10
  M[2] <- 0
  T[2] <- pow(altPriorSD,-2) 
  L[2] <- 0
  H[2] <- 10
  # Hyperprior:
  mIdx ~ dcat( mProb[] )
  mProb[1] <- 0.5
  mProb[2] <- 0.5
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
N = 40
set.seed(47405)
#SD = 2 ; M = 0.6*SD # HDI excludes zero, and alt prior wins
SD = 2 ; M = 0.4*SD # HDI excludes zero, but null prior wins
#SD = 2 ; M = 0.3*SD # HDI includes zero, and null prior wins
y = rnorm( N ) 
y = (y-mean(y))/sd(y) * SD + M

altPriorSD = c(1.5,20,50)[2]
dataType = c("Prior","Post")[2]
nullPriorSD = 0.01
if ( dataType=="Post" ) {
  datalist = list(
     y = y ,
     N = N ,
     altPriorSD = altPriorSD ,
     nullPriorSD = nullPriorSD
  )
} else {
  datalist = list(
     N = N ,
     altPriorSD = altPriorSD ,
     nullPriorSD = nullPriorSD
  )
}
# Get the data into BRugs:
modelData( bugsData( datalist ) )

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

nchain = 1
modelCompile( numChains=nchain )
modelGenInits()

#------------------------------------------------------------------------------
# RUN THE CHAINS.

burninSteps = 1000
modelUpdate( burninSteps )
samplesSet( c("mu","sigma","mIdx") )
nPerChain = 20000
modelUpdate( nPerChain , thin=100 ) # takes nPerChain * thin steps

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS.

source("plotChains.R")
source("plotPost.R")

#plotChains("muM1")

mIdxCh = samplesSample( "mIdx" )
mu1Ch = samplesSample( "mu" )[ mIdxCh == 1 ]
mu2Ch = samplesSample( "mu" )[ mIdxCh == 2 ]
muRange = range( c(mu1Ch,mu2Ch) )
sigma1Ch = samplesSample( "sigma" )[ mIdxCh == 1 ]
sigma2Ch = samplesSample( "sigma" )[ mIdxCh == 2 ]

windows(10,7)
layout( matrix( c(1,1,2,2, 3,4,4,4, 5,5,6,6 ) , nrow=4 ) , 
        heights=1+c(1,1,1,1) , widths=1+c(2,1,2)  )

# mu1
hi = plotPost( mu1Ch , xlab=bquote(mu) , xlim=muRange ,  
              main=paste("Model 1: Prior SD on mu =",nullPriorSD ) , 
              breaks=30 , col="skyblue" , cex.lab=1.75 , border="skyblue" )
# sigma1
hi = plotPost( sigma1Ch , xlab=bquote(sigma) , main="Model 1" , 
              breaks=30 , col="skyblue" , cex.lab=1.75 )

# data boxplot
par(xpd=NA)
if ( !is.null(datalist$y) ) {
  boxplot( datalist$y , horizontal=T , main="Data" ) 
  text( mean(datalist$y) , 1.5 , adj=c(0.5,1) , cex=1.5 ,
      bquote( "N=" * .(datalist$N) * 
              ", m=" * .(round(mean(datalist$y),2)) *
              ", sd=" * .(round(sd(datalist$y),2)) ) )
} else  {
  plot( 0,0,main="Empty Data for Prior" )
}

# model index
pM1 = sum( mIdxCh == 1 ) / length( mIdxCh )
pM2 = 1 - pM1
string1 =paste("p(M1|D)=",round(pM1,3),sep="")
string2 =paste("p(M2|D)=",round(pM2,3),sep="")
plot( mIdxCh[1:min(2000,length(mIdxCh))] , 1:min(2000,length(mIdxCh)) , type="l" ,
      ylab="Step in Markov chain" , xlab="Model Index (1, 2)" ,
      main=paste(string1,", ",string2,sep="") , cex.lab=1.5 , col="skyblue" )

# mu2
hi = plotPost( mu2Ch , xlab=bquote(mu) , xlim=muRange , 
              main=paste("Model 2: Prior SD on mu =",altPriorSD ) , 
              breaks=30 , col="skyblue" , cex.lab=1.75 , compVal=0.0 )

# sigma2
hi = plotPost( sigma2Ch , xlab=bquote(sigma) , main="Model 2" , 
              breaks=30 , col="skyblue" , cex.lab=1.75 )

savePlot(file=paste( filenamebase,dataType,M,altPriorSD,".eps",sep="") , type="eps" )
savePlot(file=paste( filenamebase,dataType,M,altPriorSD,".jpg",sep="") , type="jpg" )


show( t.test(y) )
# http://pcl.missouri.edu/bf-one-sample
