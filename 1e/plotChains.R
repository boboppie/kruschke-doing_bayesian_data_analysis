plotChains = function( nodename , saveplots=F , filenameroot="DeleteMe" ) {
    summarytable = samplesStats(nodename)
    show( summarytable )
    nCompon = NROW(summarytable)
    nPlotPerRow = 5
    nPlotRow = ceiling(nCompon/nPlotPerRow)
    nPlotCol = ceiling(nCompon/nPlotRow)
    windows(3.75*nPlotCol,3.5*nPlotRow)
    par( mar=c(4,4,3,1) , mgp=c(2,0.7,0) )
    samplesHistory( nodename , ask=F , mfrow=c(nPlotRow,nPlotCol) ,
                    cex.lab=1.5 , cex.main=1.5 )
    if ( saveplots ) {
       dev.copy2eps( file=paste( filenameroot , toupper(nodename) ,
                                 "history.eps" , sep="" )) }
    windows(3.75*nPlotCol,3.5*nPlotRow)
    par( mar=c(4,4,3,1) , mgp=c(2,0.7,0) )
    samplesAutoC( nodename , chain=1 , ask=F , mfrow=c(nPlotRow,nPlotCol) ,
                  cex.lab=1.5 , cex.main=1.5 )
    if ( saveplots ) {
       dev.copy2eps( file=paste( filenameroot , toupper(nodename) ,
                     "autocorr.eps" , sep="" )) }
    windows(3.75*nPlotCol,3.5*nPlotRow)
    par( mar=c(4,4,3,1) , mgp=c(2,0.7,0) )
    samplesBgr( nodename , ask=F , mfrow=c(nPlotRow,nPlotCol) ,
                cex.lab=1.5 , cex.main=1.5 )
    if ( saveplots ) {
       dev.copy2eps( file=paste( filenameroot , toupper(nodename) ,
                     "bgr.eps" , sep="" )) }
    return( summarytable )
}
