#graphics.off()
rm(list=ls(all=TRUE))
library(BRugs)        # Kruschke, J. K. (2010). Doing Bayesian data analysis:
                      # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

# Specify the model in BUGS language, but save it as a string in R:
modelString = "
# BUGS model specification begins ...
model {
    # Likelihood:
    for ( t in 1:nTrialTotal ) {
        y[t] ~ dbern( theta[ coin[ t ] ] )
    }
    # Prior:
    for ( j in 1:nCoins ) {
        theta[j] ~ dbeta( a , b )I(0.0001,0.9999)
    }
    a <- mu * kappa
    b <- ( 1.0 - mu ) * kappa
    mu ~ dbeta( Amu , Bmu )
    kappa ~ dgamma( Skappa , Rkappa )
    Amu <- 2.0
    Bmu <- 2.0
    Skappa <- pow(10,2)/pow(10,2)
    Rkappa <- 10/pow(10,2)
}
# ... BUGS model specification ends.
" # close quote to end modelString

# Write the modelString to a file, using R commands:
.temp = file("model.txt","w") ; writeLines(modelString,con=.temp) ; close(.temp)
# Use BRugs to send the model.txt file to BUGS, which checks the model syntax:
modelCheck( "model.txt" )

#------------------------------------------------------------------------------
# THE DATA.

# Demo data for various figures in the book:
# N =  c( 5, 5, 5, 5, 5 ) # c( 10, 10, 10 )  # c( 15, 5 ) # c( 5, 5, 5, 5, 5 )
# z =  c( 1, 1, 1, 1, 5 ) # c(  1,  5,  9 )  # c(  3, 4 ) # c( 1, 1, 1, 1, 5 )

ncoins = 5 ; nflipspercoin = 50 
muAct = .7 ; kappaAct = 20
thetaAct = rbeta( ncoins ,muAct*kappaAct ,(1-muAct)*kappaAct )
z = rbinom( n=ncoins ,size=nflipspercoin ,prob=thetaAct )
N = rep( nflipspercoin , ncoins )
 
# Therapeutic touch data:
#  z = c(1,2,3,3,3,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,5,5,6,6,7,7,7,8)
#  N = rep(10,length(z))
# Convert z,N to vectors of individual flips.
# coin vector is index of coin for each flip.
# y vector is head or tail for each flip.
# For example,
#  coin = c( 1, 1, 2, 2, 2 )
#  y    = c( 1, 0, 0, 0, 1 )
# means that the first flip was of coin 1 and it was a head, the second flip
# was of coin 1 and it was a tail, the third flip was of coin 2 and it was a
# tail, etc.
coin = NULL ; y = NULL
for ( coinIdx in 1:length(N) ) {
    coin = c( coin , rep(coinIdx,N[coinIdx]) )
    y = c( y , rep(1,z[coinIdx]) , rep(0,N[coinIdx]-z[coinIdx]) )
}
nTrialTotal = length( y )
nCoins = length( unique( coin ) )
dataList = list(
    y = y ,
    coin = coin ,
    nTrialTotal = nTrialTotal ,
    nCoins = nCoins
)

# Use BRugs commands to put the data into a file and ship the file to BUGS:
modelData( bugsData( dataList ) )

#------------------------------------------------------------------------------
# INTIALIZE THE CHAIN.

nChains = 3
modelCompile( numChains = nChains )  # BRugs tells BUGS to compile the model.
modelGenInits() # BRugs tells BUGS to randomly initialize the chains.

#------------------------------------------------------------------------------
# RUN THE CHAINS.

# Run some initial steps without recording them, to burn-in the chains:
burninSteps = 1000
modelUpdate( burninSteps )
# BRugs tells BUGS to keep a record of the sampled values:
samplesSet( c( "mu" , "kappa" , "theta" ) )
nPerChain = 1000
modelUpdate( nPerChain , thin=10 )

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS.

# Check for mixing and autocorrelation:
source("plotChains.R")
plotChains( "mu" , saveplots=F )
plotChains( "kappa" , saveplots=F )
plotChains( "theta[1]" , saveplots=F )

# Extract the posterior sample from BUGS:
muSample = samplesSample( "mu" ) # BRugs gets sample from BUGS
kappaSample = samplesSample( "kappa" ) # BRugs gets sample from BUGS
thetaSample = matrix( 0 , nrow=nCoins , ncol=nChains*nPerChain )
for ( coinIdx in 1:nCoins ) {
    nodeName = paste( "theta[" , coinIdx , "]" , sep="" )
    thetaSample[coinIdx,] = samplesSample( nodeName )
}

# Make a graph using R commands:
source("plotPost.R")
if ( nCoins <= 5 ) { # Only make this figure if there are not too many coins
windows(3.2*3,2.5*(1+nCoins))
layout( matrix( 1:(3*(nCoins+1)) , nrow=(nCoins+1) , byrow=T ) )
par(mar=c(2.95,2.95,1.0,0),mgp=c(1.35,0.35,0),oma=c( 0.1, 0.1, 0.1, 0.1))
nPtsToPlot = 500
plotIdx = floor(seq(1,length(muSample),length=nPtsToPlot))
kPltLim = signif( quantile( kappaSample , p=c(.01,.99) ) , 4 )
plot( muSample[plotIdx] , kappaSample[plotIdx] , type="p" , ylim=kPltLim ,
      xlim=c(0,1) , xlab=expression(mu) , ylab=expression(kappa) , cex.lab=1.5 )
plotPost( muSample , xlab="mu" , xlim=c(0,1) , main="" , breaks=20 )
plotPost( kappaSample , xlab="kappa" , main="" , breaks=20 , HDItextPlace=.6 )
for ( coinIdx in 1:nCoins ) {
    plotPost( thetaSample[coinIdx,] , xlab=paste("theta",coinIdx,sep="") ,
              xlim=c(0,1) , main="" , breaks=20 , HDItextPlace=.3 )
    plot( thetaSample[coinIdx,plotIdx] , muSample[plotIdx] , type="p" ,
          xlim=c(0,1) , ylim=c(0,1) , cex.lab=1.5 ,
          xlab=bquote(theta[.(coinIdx)]) , ylab=expression(mu) )
    plot( thetaSample[coinIdx,plotIdx] , kappaSample[plotIdx] , type="p" ,
          xlim=c(0,1) , ylim=kPltLim , cex.lab=1.5 ,
          xlab=bquote(theta[.(coinIdx)]) , ylab=expression(kappa) )
}
#dev.copy2eps(file=paste("BernBetaMuKappaBugs",paste(z,collapse=""),".eps",sep=""))
} # end if ( nCoins <= ...

## Uncomment the following if using therapeutic touch data:
windows(7,5)
layout( matrix( 1:4 , nrow=2 , byrow=T ) )
par(mar=c(2.95,2.95,1.0,0),mgp=c(1.35,0.35,0),oma=c( 0.1, 0.1, 0.1, 0.1) )
plotPost( muSample , xlab="mu" , main="" , breaks=20 , compVal=0.5 )
plotPost( kappaSample , xlab="kappa" , main="" , breaks=20 , HDItextPlace=.1 )
plotPost( thetaSample[1,] , xlab="theta1" , main="" , breaks=20 , compVal=0.5 )
#plotPost( thetaSample[28,] , xlab="theta28" , main="" , breaks=20 , compVal=0.5 )
#dev.copy2eps(file="BernBetaMuKappaBugsTT.eps")