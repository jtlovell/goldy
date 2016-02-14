#' export
autoGoldy<-function(filename, trait1="Photo", trait2="Cond", window=6, verbose=TRUE,weights=c(1,1,1,1)){

  if(verbose) cat("parsing LICOR datafile ...\n")
  dat<-parseLICOR(filename)
  par(mfrow=c(3,1), mar=c(5,5,5,5))

  if(verbose) cat("standardizing traits: ", trait1, trait2, "...\n")
  stdplot<-twoPlot(dat=dat, trait1=trait1, trait2=trait2)

  if(verbose) cat("running sliding window for convergence ...\n")
  sw.t1<-cvSlopeSW(stdplot$std.y1, window=window)
  sw.t2<-cvSlopeSW(stdplot$std.y2, window=window)
  lines(stdplot$x, sw.t1$mean, col="red")
  lines(stdplot$x, sw.t2$mean, col="red")

  if(verbose) cat("defining optimal window ...\n")
  ranks<-weightedRank(x=stdplot$x, slope1=sw.t1$slope, slope2=sw.t2$slope,
                      cv1=sw.t1$cv, cv2=sw.t2$cv, window=window,weights=weights)
  plotDiagnotics(x=stdplot$x, slope1=sw.t1$slope, slope2=sw.t2$slope,
                 cv1=sw.t1$cv, cv2=sw.t2$cv, weightedRankOutput=ranks,
                 trait1=trait1, trait2=trait2,
                 inputFile=filename)
  out<-getBV(dat, ranks, filename=filename, window=window)
  return(out)
}
