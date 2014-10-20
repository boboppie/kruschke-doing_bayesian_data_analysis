# MCMC demo. Simple Metropolis algorithm applied to 1D discrete values, with
# proposal distribution being merely to try one value up or one value down.
graphics.off()

# Specify probability mass in each interval
# Comment out all but one of spec's below
#
relprob = c(0,1,2,3,4,5,6,7,0) # MUST HAVE ZERO AT EACH END!
distribname = "incline"
#
#relprob = c(0,1,2,3,4,5,0) # MUST HAVE ZERO AT EACH END!
#distribname = "shortincline"
#
#relprob = c(0,1,2,3,4,5,1,2,3,4,0) # MUST HAVE ZERO AT EACH END!
#distribname = "sawtooth"

############################################################################
# Part I: Plot the probability of being in each position as a function of 
# update number. This is not an individual trajectory, but a distribution
# of probabilities that the particle is in each possible position.
############################################################################

# pin is the target relative probability distribution with matrix attributes.
pin = matrix( relprob , nrow=1 , ncol=length(relprob) , byrow=TRUE )

# Specify transition matrix, tm. Rows are from, columns are to.
# This assumes that the proposal distribution is either choosing one interval 
# up or one interval down.
pup = .5 # MUST BE <=.5, prob. of choosing interval above current interval
pdown = pup # probability of choosing interval below current interval
tm = matrix( 0 , nrow=length(relprob) , ncol=length(relprob) )
# NEXT STATEMENT RELIES ON ZEROS AT ENDS OF relprob
for ( rowint in (1+1) : (length(relprob)-1) ) { 
	# determine p(rowint-1)
	tm[rowint,rowint-1] = pdown * min( pin[1,rowint-1]/pin[1,rowint] , 1 )
	# determine p(rowint+1)
	tm[rowint,rowint+1] = pup * min( pin[1,rowint+1]/pin[1,rowint] , 1 )
	# determine p(rowint)
	tm[rowint,rowint] = 1 - tm[rowint,rowint-1] - tm[rowint,rowint+1]	
} 

# Specify initial probability vector. Must have length same as relprob.
pvec = matrix( 0 , nrow=1 , ncol=length(relprob) )
pvec[1,ceiling(median(1:length(pvec)))] = 1

# Specify number of updates to execute
nupdate = 99
# Specify graphics window details.
nr = 4  # 3 for shortincline, 5 otherwise
nc = 4
windows(7,7)
layout( matrix(1:nupdate,nrow=nr,ncol=nc,byrow=FALSE) )
par(mar=c(2,2,0,0)) # number of margin lines: bottom,left,top,right
par(mgp=c(1,1,0)) # which margin lines to use for labels
par(mai=c(0.25,0.25,0.05,0.05)) # margin size in inches: bottom,left,top,right

# Plot the updated distributions
for ( t in 1:nupdate ) {
	if ( t <= (nr*nc - 2) ) {
		plot( 1:(length(pvec)-2) , pvec[2:(length(pvec)-1)] 
			, type="h" , lwd=4 , xlab="" , ylab="" 
			,ylim=c(0,max(pvec[2:(length(pvec)-1)])) )
		text( 1,max(pvec),bquote(t==.(t)) ,adj=c(0,1) ,cex=1.5)
	}
	pvec = pvec %*% tm
}

# Plot the final updated distribution
plot( 1:(length(pvec)-2) , pvec[2:(length(pvec)-1)] 
	, type="h" , lwd=4 , xlab="" , ylab="" 
	,ylim=c(0,max(pvec[2:(length(pvec)-1)])) )
text( 1,max(pvec),bquote(t==.(t)) ,adj=c(0,1) ,cex=1.5)

# Plot the target distribution.
prob = relprob 
plot( 1:(length(pvec)-2) , prob[2:(length(prob)-1)] 
		, type="h" , lwd=4 , xlab="" , ylab="" 
		,ylim=c(0,max(prob[2:(length(prob)-1)])) )
text( 1,max(prob),bquote(target) ,adj=c(0,1) ,cex=1.5)

# Save the plot as an EPS file.
filename = paste( "BinomMCMCdemo_", distribname ,".eps",sep="")
dev.copy2eps(file=filename)


############################################################################
# Part II: Plot some individual trajectories.
############################################################################

set.seed(47)

trajlength = 2000
trajectory = rep( 0 , trajlength )

# Start in middle of range
trajectory[1] = ceiling(median(1:length(relprob)))

for ( t in 1:(trajlength-1) ) {
	currentposition = trajectory[t]
	proposedjump = sample( c(-1,1) , 1 )
	probaccept = min( 1, 
		relprob[ currentposition + proposedjump ]
		/ relprob[ currentposition ] )
	if ( runif(1) < probaccept ) {
		trajectory[ t+1 ] = trajectory[t] + proposedjump
	} else {
		trajectory[ t+1 ] = trajectory[t]
	}
}


trajectory = trajectory-1
trajectoryX = 1:(length(relprob)-2)

windows(7,10)
layout( matrix( c( 1,1,1,1, 2,2, 3,3,3,3,3,3, 4,4, 5,5,5,5 )
	, nrow=9 , ncol=2 , byrow=TRUE ) )

# plot 1 is histogram
par(mar=c(2,2,0,0)) # number of margin lines: bottom,left,top,right
par(mgp=c(2,1,0)) # which margin lines to use for labels
par(mai=c(0.4,0.45,0.05,0.05)) # margin size in inches: bottom,left,top,right
burnin=200
hist( trajectory[burnin:length(trajectory)] 
	, breaks = c(0,trajectoryX)+.5
	, xlab=bquote(theta) ,main="" ,cex.lab=1.5 )


# make upward arrow 
# Specify margin measurements for subplot.
par(mai=c(0.0,0.4,0.0,0.0)) # margin size in inches: bottom,left,top,right
plot( 0 , 0 ,type="n" ,bty="n" ,xaxt="n" ,yaxt="n" ,xlab="" ,ylab="" )
text(0,0 ,bquote( " " %dblup% " " ) ,cex=4 )


# plot 3 is trajectory
par(mar=c(3,2,0,0)) # number of margin lines: bottom,left,top,right
par(mgp=c(2,1,0)) # which margin lines to use for labels
par(mai=c(0.4,0.45,0.05,0.05)) # margin size in inches: bottom,left,top,right
yvec = seq( trajlength , 1 , by=(-1) )
plot( trajectory , 1:trajlength , type='o' 
	, xlab=bquote(theta)
	,log='y' , ylab="Time" , main="" ,cex.lab=1.5)

# plot 4 is upward arrow
# Specify margin measurements for subplot.
par(mai=c(0.0,0.4,0.0,0.0)) # margin size in inches: bottom,left,top,right
plot( 0 , 0 ,type="n" ,bty="n" ,xaxt="n" ,yaxt="n" ,xlab="" ,ylab="" )
text(0,0 ,bquote( " " %dblup% " " ) ,cex=4 )

# plot 5 is target relprob
par(mar=c(2,2,0,0)) # number of margin lines: bottom,left,top,right
par(mgp=c(2,1,0)) # which margin lines to use for labels
par(mai=c(0.4,0.45,0.05,0.05)) # margin size in inches: bottom,left,top,right
plot( trajectoryX , relprob[2:(length(relprob)-1)] 
	, type='h' ,lwd=4
      , xlab=bquote(theta)
	, ylim = c(0,max(relprob)) , ylab=expression(P(theta))  ,cex.lab=1.5)

# Save the plot as an EPS file.
# To save graphs, please see update at
# http://doingbayesiandataanalysis.blogspot.com/2013/01/uniform-r-code-for-opening-saving.html
#filename = paste( "BinomMCMCdemowalk_", distribname ,".eps",sep="")
#dev.copy2eps(file=filename)

