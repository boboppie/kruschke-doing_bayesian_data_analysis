# Grid on 1st-level parameter
b=7 ; nbin = 2*(2*b+1) ; binwidth = 1/nbin
theta = seq( from=binwidth/2 , to=1.0-binwidth/2 , by=binwidth ) 
# Grid on hyperparameter
mu = seq( from=binwidth/2 , to=1.0-binwidth/2 , by=binwidth ) 

#############################################################################
## UNCOMMENT JUST ONE OF THE FOLLOWING SETS OF SPECIFICATIONS

setSpec = c("twoval","betaSmallK","betaLargeK")[3]
if ( setSpec=="twoval" ) {
 # Two-value mu for model-comparison hyperprior
 hypertype = "twoval"
 priorMuN = 12
 datay = 6 ; dataN = 9
}
if ( setSpec=="betaSmallK" ) {
 # Beta hyperprior, low dependence
 hypertype = "beta"
 priorMuN = 6
 hyper_a = 20 ; hyper_b = 20
 datay = 9 ; dataN = 12
}
if ( setSpec=="betaLargeK" ) {
 # Beta hyperprior, high dependence
 hypertype = "beta"
 priorMuN = 100
 hyper_a = 2 ; hyper_b = 2
 datay = 9 ; dataN = 12
}

#############################################################################

# Specify prior
# Prior at 1st level
# Small priorMuN implies low certainty about dependence of theta on mu.
# Large priorMuN implies high certainty about dependence of theta on mu.
pThetaAtMu = function( theta , mu , n=priorMuN ) {
	dbeta(theta ,n*mu ,(1-mu)*n)
}
# Hyperprior
if ( hypertype == "twoval" ) {
	pMu = function( mu ) { 
		valone = .25 ; valtwo = .75 ; slicewidth = .001 
		( abs(mu-valone)<slicewidth ) | ( abs(mu-valtwo)<slicewidth ) 
	}
}
if ( hypertype == "beta" ) {
	# Small values for hyper_a,hyper_b imply low certainty about mu.
	# Large values for hyper_a,hyper_b imply high certainty about mu.
	pMu = function( mu ) { dbeta( mu , hyper_a , hyper_b ) }
}
# Joint prior
pThetaAndMu = function( theta , mu ) { pThetaAtMu(theta,mu) * pMu(mu) }
jointprior = outer( theta , mu , pThetaAndMu )
jointprior = jointprior / sum( jointprior ) # mass at discrete points
jointprior = jointprior / (binwidth^2) # density
# Marginal of joint prior
priorMuMarg = colSums( jointprior ) * binwidth # density
priorThetaMarg = rowSums( jointprior ) * binwidth # density

# Specify likelihood

likeatjointpoint = function( t , m , y=datay , N=dataN ) { 
	t^y * (1-t)^(N-y) # notice this depends only on t, not on m
}
likelihood = outer( theta , mu , likeatjointpoint )

# Determine posterior
pData = sum( jointprior * likelihood )
posterior = ( jointprior * likelihood ) / pData
posterior = posterior / (binwidth^2) # density
# Marginals of joint posterior
postMuMarg = colSums(posterior) * binwidth # density
postThetaMarg = rowSums(posterior) * binwidth # density

#------ plot ----------------------------

ipp = 3 # inches per plot
windows(3.2*ipp,5*ipp)
layout( 
	matrix( c( 	1,1,   2,2,   3,3,
			1,1,   2,2,   3,3,
			4,4,   5,5,   6,6,
			4,4,   5,5,   7,7,
			8,8,   9,9,   10,10,
			8,8,   9,9,   10,10,
			11,11, 12,12, 13,13,
			11,11, 12,12, 13,13,
			14,14, 15,15, 16,16,
			14,14, 15,15, 17,17
			) 
	,nrow=10,ncol=6,byrow=TRUE ) 
) 

cexfac = (0.75)
par(cex=cexfac,cex.axis=cexfac,cex.lab=cexfac,cex.main=cexfac*1.2,cex.sub=cexfac)
par(mex=0.9*cexfac)
par(mar=c(2.95,2.95,1.0,0)) # number of margin lines: bottom,left,top,right
par(tcl=-0.25) # tick length as proportion of character ht
par(mgp=c(1.35,0.35,0)) # which margin lines title,label,line
par( oma = c( 0.1, 0.1, 0.1, 0.1) ) # outer margin

rotate = (-25)
tilt = 25
parallelness = 5.0
shadeval = 0.05
perspcex = 0.7
ncontours = 9
zmax = max(c(max(posterior),max(jointprior)))
muMargMax = max(c(max(priorMuMarg),max(postMuMarg)))
thetaMargMax = max(c(max(priorThetaMarg),max(postThetaMarg)))

# 1
par(pty="m") 
persp( theta ,mu ,jointprior ,main="Prior" ,cex=perspcex ,zlab="prior"
	,xlim=c(0,1) ,ylim=c(0,1) ,zlim=c(0,zmax) ,lwd=0.1
	,theta=rotate, phi=tilt ,d=parallelness ,shade=shadeval )

#2
par(pty="m") 
contour( theta ,mu ,jointprior ,main=bquote(" ") 
	,levels=signif(seq(0,zmax,length=ncontours),3) ,drawlabels=FALSE
	,xlab=bquote(theta) ,ylab=bquote(mu) )

#3
par(pty="m") 
plot( priorMuMarg ,mu ,type="l"
	,ylab=bquote(mu)
	,xlab=bquote("Marginal p("*mu*")") ,xlim=c(0,muMargMax) )

#4
plot( c(0,1),c(0,1) ,type="n" ,bty="n" ,xaxt="n" ,yaxt="n" ,xlab="" ,ylab="" )
text( 1,.8,"Prior",adj=c(1,-0.1),cex=1.5)
if ( hypertype == "beta" ) {
text( 1,0.8,bquote( list( A[mu]==.(hyper_a) , B[mu]==.(hyper_b) ) ),adj=c(1,1.2),cex=1.0)
}
text( 1,0.8,bquote( list( K==.(priorMuN) ) ),adj=c(1,2.6),cex=1.0)



#5
par(pty="m") 
plot( theta, priorThetaMarg ,type="l"
	,xlab=bquote(theta)
	,ylab=bquote("Marginal p("*theta*")") ,ylim=c(0,thetaMargMax) )

#6
plotmuval=.75
plotmuidx = which.min( abs(mu-plotmuval) )
plotmuval = mu[plotmuidx]
ymax = max( c( max( jointprior[,plotmuidx]/sum(jointprior[,plotmuidx]) ) 
	, max( posterior[,plotmuidx]/sum(posterior[,plotmuidx]) )))*(1/binwidth)
par(pty="m") 
plot( theta , jointprior[,plotmuidx]/sum(jointprior[,plotmuidx])*(1/binwidth)
	, type="l" ,ylim=c(0,ymax) ,xlab=bquote(theta)
	,ylab=bquote("p("*theta*"|"*mu*"=.75)")
	,cex.lab=cexfac*.8 ,cex.axis=cexfac*.8 )

#7
plotmuval=.25
plotmuidx = which.min( abs(mu-plotmuval) )
plotmuval = mu[plotmuidx]
ymax = max( c( max( jointprior[,plotmuidx]/sum(jointprior[,plotmuidx]) ) 
	, max( posterior[,plotmuidx]/sum(posterior[,plotmuidx]) )))*(1/binwidth)
par(pty="m") 
plot( theta , jointprior[,plotmuidx]/sum(jointprior[,plotmuidx])*(1/binwidth)
	, type="l" ,ylim=c(0,ymax) ,xlab=bquote(theta)
	,ylab=bquote("p("*theta*"|"*mu*"=.25)")
	,cex.lab=cexfac*.8 ,cex.axis=cexfac*.8 )

#8
par(pty="m") 
persp( theta,mu, likelihood ,main="Likelihood" ,cex=perspcex ,lwd=0.1
	,theta=rotate, phi=tilt ,d=parallelness ,shade=shadeval )

#9
par(pty="m") 
contour( theta ,mu ,likelihood ,main=bquote(" ") ,nlevels=(ncontours-1)
	,xlab=bquote(theta) ,ylab=bquote(mu) ,drawlabels=FALSE )

#10
plot( c(0,1),c(0,1) ,type="n" ,bty="n" ,xaxt="n" ,yaxt="n" ,xlab="" ,ylab="" )
text( 0.0,0.5,"Likelihood",adj=c(0,-0.2),cex=1.5)
text( 0.0,0.5
	,bquote( "D = " * .(datay) *" heads, "* .(dataN-datay) *" tails" )
	,adj=c(0,1.2),cex=1.0)

#11
par(pty="m") 
persp( theta,mu, posterior ,main="Posterior" ,cex=perspcex
	,xlim=c(0,1) ,ylim=c(0,1) ,zlim=c(0,zmax) ,lwd=0.1
	,theta=rotate, phi=tilt ,d=parallelness ,shade=shadeval )

#12
par(pty="m") 
contour( theta ,mu ,posterior ,main=bquote(" ") 
	,levels=signif(seq(0,zmax,length=ncontours),3) ,drawlabels=FALSE
	,xlab=bquote(theta) ,ylab=bquote(mu) )

#13
par(pty="m") 
plot( postMuMarg ,mu ,type="l"
	,ylab=bquote(mu)
	,xlab=bquote("Marginal p("*mu*"|D)") ,xlim=c(0,muMargMax) )

#14
plot( c(0,1),c(0,1) ,type="n" ,bty="n" ,xaxt="n" ,yaxt="n" ,xlab="" ,ylab="" )
text( 1,.8,"Posterior",adj=c(1,-0.1),cex=1.5)

#15
par(pty="m") 
plot( theta, postThetaMarg ,type="l"
	,xlab=bquote(theta)
	,ylab=bquote("Marginal p("*theta*"|D)") ,ylim=c(0,thetaMargMax) )

#16
plotmuval=.75
plotmuidx = which.min( abs(mu-plotmuval) )
plotmuval = mu[plotmuidx]
ymax = max( c( max( jointprior[,plotmuidx]/sum(jointprior[,plotmuidx]) ) 
	, max( posterior[,plotmuidx]/sum(posterior[,plotmuidx]))))*(1/binwidth)
par(pty="m") 
plot( theta , posterior[,plotmuidx]/sum(posterior[,plotmuidx])*(1/binwidth) 
	, type="l" ,ylim=c(0,ymax) ,xlab=bquote(theta)
	,ylab=bquote("p("*theta*"|"*mu*"=.75,D)")
	,cex.lab=cexfac*.8 ,cex.axis=cexfac*.8 )

#17
plotmuval=.25
plotmuidx = which.min( abs(mu-plotmuval) )
plotmuval = mu[plotmuidx]
ymax = max( c( max( jointprior[,plotmuidx]/sum(jointprior[,plotmuidx]) ) 
	, max( posterior[,plotmuidx]/sum(posterior[,plotmuidx]))))*(1/binwidth)
par(pty="m") 
plot( theta , posterior[,plotmuidx]/sum(posterior[,plotmuidx])*(1/binwidth) 
	, type="l" ,ylim=c(0,ymax) ,xlab=bquote(theta)
	,ylab=bquote("p("*theta*"|"*mu*"=.25,D)") 
	,cex.lab=cexfac*.8 ,cex.axis=cexfac*.8 )

want_eps_file = FALSE
# To save graphs, please see update at
# http://doingbayesiandataanalysis.blogspot.com/2013/01/uniform-r-code-for-opening-saving.html
if ( want_eps_file ) {
	if ( hypertype == "beta" ) {
		epsfilename = paste( "BinomHierGrid" 
			,"_",bquote(.(priorMuN))
			,"_",bquote(.(hyper_a))
			,"_",bquote(.(hyper_b))
			,"_",bquote(.(datay))
			,"_",bquote(.(dataN))
			,".eps" ,sep="" ) 
	}
	if ( hypertype == "twoval" ) {
		epsfilename = "BinomHierGridTwoval.eps"
	}
	dev.copy2eps( file = epsfilename )
}