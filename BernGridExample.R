graphics.off()
source("openGraphSaveGraph.R") # for openGraph() function, used below.
source("BernGrid.R")

# For Figure 6.4:
# Specify theta values.
thetagrid = seq(0,1,length=1001)
# Specify probability mass at each theta value.
relprob = sin( 2*pi*thetagrid )^6
prior = relprob / sum(relprob) # probability mass at each theta
# Specify the data vector.
datavec = c( rep(1,2) , rep(0,1) ) 
# Open a window.
openGraph(width=7,height=10,mag=0.7)
# Call the function.
posterior = BernGrid( Theta=thetagrid , pTheta=prior , Data=datavec )
saveGraph(file="Fig.6.4",type="jpg")

# For Figure 6.5:
pTheta = c( 50:1 , rep(1,50) , 1:50 , 50:1 , rep(1,50) , 1:50 )
pTheta = pTheta / sum( pTheta )
width = 1 / length(pTheta)
Theta = seq( from = width/2 , to = 1-width/2 , by = width )
dataVec = c( rep(1,3) , rep(0,1) )
openGraph(width=7,height=10,mag=0.7)
posterior = BernGrid( Theta=Theta , pTheta=pTheta , Data=dataVec )
saveGraph(file="Fig.6.5left",type="jpg")
dataVec = c( rep(1,12) , rep(0,4) )
openGraph(width=7,height=10,mag=0.7)
posterior = BernGrid( Theta=Theta , pTheta=posterior , Data=dataVec )
saveGraph(file="Fig.6.5right",type="jpg")
