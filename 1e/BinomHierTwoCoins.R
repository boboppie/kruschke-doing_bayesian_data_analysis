graphics.off()

# Grids
b=12 ; nbin = 2*(2*b+1) ; binwidth = 1/nbin
theta1 = seq( from=binwidth/2 , to=1.0-binwidth/2 , by=binwidth ) 
theta2 = seq( from=binwidth/2 , to=1.0-binwidth/2 , by=binwidth ) 
mu = seq( from=binwidth/2 , to=1.0-binwidth/2 , by=binwidth ) 

# Specify prior
n = 75 ; c = 2 ; d = 2
# n = 5 ; c = 2 ; d = 2
pt1t2mu = function( mu , theta1 , theta2 ) {
	pmu = dbeta( mu , c , d )
	pt1gmu = dbeta( theta1 , n*mu , (1-mu)*n )
	pt2gmu = dbeta( theta2 , n*mu , (1-mu)*n )
	pt1t2mu = pmu * pt1gmu * pt2gmu
}
prior = array( 0 , dim = c( length(mu) , length(theta1) , length(theta2 ) ) )
for ( t2 in 1:length(theta2) ) {
	for ( t1 in 1:length(theta1) ) {
		for ( m in 1:length(mu) ) {
			prior[m,t1,t2] = pt1t2mu( mu[m],theta1[t1],theta2[t2] )
		}
	}
}
prior = prior / sum(prior)

# Marginals of joint prior
priorMuMarg = apply( prior , 1 , sum )
priorTheta1Marg = apply( prior , 2 , sum )
priorTheta2Marg = apply( prior , 3 , sum )
priorMuTheta1Marg = apply( prior , c(1,2) , sum )
priorMuTheta2Marg = apply( prior , c(1,3) , sum )

# Specify likelihood
y1 = 3 ; N1 = 15 ; y2 = 4 ; N2 = 5
likelihood = array( 0 , dim = c( length(mu) , length(theta1) , length(theta2 ) ) )
for ( t2 in 1:length(theta2) ) {
	for ( t1 in 1:length(theta1) ) {
		for ( m in 1:length(mu) ) {
			likelihood[m,t1,t2] = (
                theta1[t1]^y1 * (1-theta1[t1])^(N1-y1) *
                theta2[t2]^y2 * (1-theta2[t2])^(N2-y2) )
		}
	}
}
# Marginals of likelihood
likelihoodMuTheta1Marg = apply( likelihood , c(1,2) , sum )
likelihoodMuTheta2Marg = apply( likelihood , c(1,3) , sum )


# Determine posterior
pData = sum( prior * likelihood )
posterior = ( prior * likelihood ) / pData

# Marginals of joint posterior
posteriorMuMarg = apply( posterior , 1 , sum )
posteriorTheta1Marg = apply( posterior , 2 , sum )
posteriorTheta2Marg = apply( posterior , 3 , sum )
posteriorMuTheta1Marg = apply( posterior , c(1,2) , sum )
posteriorMuTheta2Marg = apply( posterior , c(1,3) , sum )

#------ plot ----------------------------

MuThetaMax = max(c(
	max(priorMuTheta1Marg),
	max(priorMuTheta2Marg),
	max(posteriorMuTheta1Marg),
	max(posteriorMuTheta2Marg)))

MargMax = max(c(
	max(priorMuMarg), 
	max(priorTheta1Marg) ,
	max(priorTheta2Marg) ,
	max(posteriorMuMarg), 
	max(posteriorTheta1Marg) ,
	max(posteriorTheta2Marg) ))

windows(4.5,6.5)
layout( matrix( 1:15 ,nrow=5,ncol=3,byrow=TRUE ) ) 

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

# 1
par(pty="m") 
contour( theta1 ,mu ,t(priorMuTheta1Marg) ,main=bquote(" ") 
	,drawlabels=FALSE
	,levels=signif(seq(0,MuThetaMax,length=ncontours),3)
	,xlab=bquote(theta[1]) ,ylab=bquote(mu) )
text( 0,1,bquote(p(theta[1],mu)),adj=c(0,1),cex=0.75)

#2
par(pty="m") 
contour( theta2 ,mu ,t(priorMuTheta2Marg) ,main=bquote(" ") 
	,drawlabels=FALSE
	,levels=signif(seq(0,MuThetaMax,length=ncontours),3)
	,xlab=bquote(theta[2]) ,ylab=bquote(mu) )
text( 0,1,bquote(p(theta[2],mu)),adj=c(0,1),cex=0.75)

#3
par(pty="m") 
plot( priorMuMarg ,mu ,type="l"
	,ylab=bquote(mu)
	,xlab=bquote("p("*mu*")") ,xlim=c(0,MargMax)  )

#4
par(pty="m") 
plot( theta1, priorTheta1Marg ,type="l"
	,xlab=bquote(theta[1]) ,ylim=c(0,MargMax)
	,ylab=bquote("p("*theta[1]*")")  )


#5
par(pty="m") 
plot( theta2, priorTheta2Marg ,type="l"
	,xlab=bquote(theta[2]) ,ylim=c(0,MargMax)
	,ylab=bquote("p("*theta[2]*")")  )

#6
plot( c(0,1),c(0,1) ,type="n" ,bty="n" ,xaxt="n" ,yaxt="n" ,xlab="" ,ylab="" )
text( 0,.8,"Prior",adj=c(0,-0.1),cex=1.5)
text( 0.0,0.8,bquote( list( A[mu]==.(c) , B[mu]==.(d) ) ),adj=c(0,1.2),cex=1.0)
text( 0.0,0.8,bquote( list( K==.(n) ) ),adj=c(0,2.6),cex=1.0)

likelihoodMax = max( c( max(likelihoodMuTheta1Marg) 
			, max(likelihoodMuTheta2Marg) ) )

# 7
par(pty="m") 
contour( theta1 ,mu ,t(likelihoodMuTheta1Marg) ,main=bquote(" ") 
	,drawlabels=FALSE
	,levels=signif(seq(0,likelihoodMax,length=ncontours),3)
	,xlab=bquote(theta[1]) ,ylab=bquote(mu) )

#8
par(pty="m") 
contour( theta2 ,mu ,t(likelihoodMuTheta2Marg) ,main=bquote(" ") 
	,drawlabels=FALSE
	,levels=signif(seq(0,likelihoodMax,length=ncontours),3)
	,xlab=bquote(theta[2]) ,ylab=bquote(mu) )

#9
plot( c(0,1),c(0,1) ,type="n" ,bty="n" ,xaxt="n" ,yaxt="n" ,xlab="" ,ylab="" )
text( 0.0,0.5,"Likelihood",adj=c(0,-0.2),cex=1.5)
text( 0.0,0.5
	,bquote("D1: "* .(y1) *" heads, "* .(N1-y1) *" tails" )
	,adj=c(0,1.2),cex=1.0)
text( 0.0,0.5
	,bquote("D2: "* .(y2) *" heads, "* .(N2-y2) *" tail" )
	,adj=c(0,2.4),cex=1.0)

#10
par(pty="m") 
contour( theta1 ,mu ,t(posteriorMuTheta1Marg) ,main=bquote(" ") 
	,drawlabels=FALSE
	,levels=signif(seq(0,MuThetaMax,length=ncontours),3)
	,xlab=bquote(theta[1]) ,ylab=bquote(mu) )
text( 0,1,bquote(p(theta[1],mu*"|"*D)),adj=c(0,1),cex=0.75)

#11
par(pty="m") 
contour( theta2 ,mu ,t(posteriorMuTheta2Marg) ,main=bquote(" ") 
	,drawlabels=FALSE
	,levels=signif(seq(0,MuThetaMax,length=ncontours),3)
	,xlab=bquote(theta[2]) ,ylab=bquote(mu) )
text( 0,1,bquote(p(theta[2],mu*"|"*D)),adj=c(0,1),cex=0.75)

#12
par(pty="m") 
plot( posteriorMuMarg ,mu ,type="l"
	,ylab=bquote(mu) ,xlim=c(0,MargMax)
	,xlab=bquote("p("*mu*"|D)")  )

#13
par(pty="m") 
plot( theta1, posteriorTheta1Marg ,type="l"
	,xlab=bquote(theta[1]) ,ylim=c(0,MargMax)
	,ylab=bquote("p("*theta[1]*"|D)")  )


#14
par(pty="m") 
plot( theta2, posteriorTheta2Marg ,type="l"
	,xlab=bquote(theta[2]) ,ylim=c(0,MargMax)
	,ylab=bquote("p("*theta[2]*"|D)")  )

#15
plot( c(0,1),c(0,1) ,type="n" ,bty="n" ,xaxt="n" ,yaxt="n" ,xlab="" ,ylab="" )
text( 0,1,"Posterior",adj=c(0,1),cex=1.5)

want_eps_file = FALSE
# To save graphs, please see update at
# http://doingbayesiandataanalysis.blogspot.com/2013/01/uniform-r-code-for-opening-saving.html
if ( want_eps_file ) {
	epsfilename = paste( "BinomHierTwoCoins" 
		,"_",bquote(.(n))
		,"_",bquote(.(c))
		,"_",bquote(.(d))
		,"_",bquote(.(y1))
		,"_",bquote(.(N1))
		,"_",bquote(.(y2))
		,"_",bquote(.(N2))
		,".eps" ,sep="" ) 
	dev.copy2eps( file = epsfilename )
}