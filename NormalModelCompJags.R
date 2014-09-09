graphics.off()
rm(list=ls(all=TRUE))
filenamebase = "NormalModelCompJags"
source("openGraphSaveGraph.R")
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.

#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
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
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.

# Specify the data:
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
  dataList = list(
     y = y ,
     N = N ,
     altPriorSD = altPriorSD ,
     nullPriorSD = nullPriorSD
  )
} else {
  dataList = list(
     N = N ,
     altPriorSD = altPriorSD ,
     nullPriorSD = nullPriorSD
  )
}

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

# Let JAGS do it randomly.

#------------------------------------------------------------------------------
# RUN THE CHAINS.


parameters = c("mu","sigma","mIdx")
adaptSteps = 500              # Number of steps to "tune" the samplers.
burnInSteps = 1000            # Number of steps to "burn-in" the samplers.
nChains = 3                   # Number of chains to run.
numSavedSteps=100000         # Total number of steps in chains to save.
thinSteps=10                   # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
jagsModel = jags.model( "model.txt" , data=dataList , # inits=initsList , 
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( jagsModel , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                            n.iter=nPerChain , thin=thinSteps )
# resulting codaSamples object has these indices: 
#   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]


#------------------------------------------------------------------------------
# EXAMINE THE RESULTS.

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )

source("plotPost.R")


mIdxCh = mcmcChain[, "mIdx" ]
mu1Ch = mcmcChain[, "mu" ][ mIdxCh == 1 ]
mu2Ch = mcmcChain[, "mu" ][ mIdxCh == 2 ]
muRange = range( c(mu1Ch,mu2Ch) )
sigma1Ch = mcmcChain[, "sigma" ][ mIdxCh == 1 ]
sigma2Ch = mcmcChain[, "sigma" ][ mIdxCh == 2 ]
sigmaRange = range( c(sigma1Ch,sigma2Ch) )

openGraph(10,7)
layout( matrix( c(1,1,2,2, 3,4,4,4, 5,5,6,6 ) , nrow=4 ) , 
        heights=1+c(1,1,1,1) , widths=1+c(2,1,2)  )

# mu1
hi = plotPost( mu1Ch , xlab=bquote(mu) , xlim=muRange ,  
              main=paste("Model 1: Prior SD on mu =",nullPriorSD ) , 
               cex.lab=1.75 , border="skyblue" )
# sigma1
hi = plotPost( sigma1Ch , xlab=bquote(sigma) , main="Model 1" , xlim=sigmaRange , 
               showMode=TRUE , cex.lab=1.75 )

# data boxplot
par(xpd=NA)
if ( !is.null(dataList$y) ) {
  boxplot( dataList$y , horizontal=T , main="Data" ) 
  text( mean(dataList$y) , 1.5 , adj=c(0.5,1) , cex=1.5 ,
      bquote( "N=" * .(dataList$N) * 
              ", m=" * .(round(mean(dataList$y),2)) *
              ", sd=" * .(round(sd(dataList$y),2)) ) )
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
              cex.lab=1.75 , compVal=0.0 )

# sigma2
hi = plotPost( sigma2Ch , xlab=bquote(sigma) , main="Model 2" , xlim=sigmaRange , 
              showMode=TRUE , cex.lab=1.75 )

#saveGraph(file=paste( filenamebase,dataType,M,altPriorSD,sep="") , type="eps" )
#saveGraph(file=paste( filenamebase,dataType,M,altPriorSD,sep="") , type="jpg" )

show( t.test(y) )
# http://pcl.missouri.edu/bf-one-sample
