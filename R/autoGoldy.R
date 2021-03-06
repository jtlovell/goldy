#' @title Automated processing of LICOR / IRGA Data
#'
#' @description
#' \code{autoGoldy} Define the most stable set of observations in time-series infrared
#' gas exchange data. Currently optimized for LICOR 6400 .rtf data output.
#'
#' @param filename The path and name of a file to parse and optimize
#' @param trait1 Character string that matches the name of a column in the input dataset. This
#' data will be used for optimization computation
#' @param trait2 A (optional) second trait to assist with optimization
#' @param window The size of the window in the sliding window optimization protocol
#' @param verbose Logical, should updates be printed
#' @param weights A numeric vector of length 4. These values indicate the weighting of the
#' slope [1] and coeficient of variation (cv) [2] of trait1 and trait2 [3,4]
#' @param ... additional arguments, not currently in use.
#'
#' @details This function runs the following pipeline:
##' \itemize{
##'  \item{1. }{Import and parse rtf LICOR datasets}
##'  \item{2. }{Standardize (0-1) and plot the raw data}
##'  \item{3. }{Run sliding window analysis to calculate the slope and CV of standardized data}
##'  \item{4. }{Find the most stable window}
##'  \item{5. }{Plot diagnostics}
##'  \item{6. }{Output the optimal data}
##' }
##'
#' @return a single row dataframe that gives the most stable data from the time series
#' data input originally. A single page figure with some information is also printed.
#' @examples
#' ## Not Run out<-autoGoldy(filename="./data/RO_052714_rr11_h57s", window=4, weights=c(1,1,.5,.5))
#' @import  zoo
#' @export
autoGoldy<-function(filename=NULL,dirname=NULL, trait1="Photo", trait2="Cond",
                    window=6, verbose=TRUE,weights=c(1,1,1,1), output.dir=NULL){
  if(is.null(output.dir)) {
    output.dir<-getwd()
  }
  if(length(filename)==1){
    if(verbose) cat("parsing LICOR datafile ...\n")
    dat<-parseLICOR(filename)
    par(mfrow=c(3,1), mar=c(5,5,5,5))

    if(verbose) cat("standardizing traits: ", trait1, trait2, "...\n")
    stdplot<-twoPlot(dat=dat, trait1=trait1, trait2=trait2, title=filename)
    dat$time.seconds<-stdplot$x

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
    out<-getBV(dat=dat, ranks=ranks, filename=filename, window=window)
    return(out)

  }else{
    if(!is.null(dirname)){
      setwd(dirname)
      filename<-dir(dirname)
      print(filename)
    }
    if(verbose) cat("checking column names\n")
    cnames<-unique(unlist(lapply(filename, function(x) colnames(parseLICOR(x)))))

    if(verbose) cat("parsing", length(filename), "irga files\n")
    if(verbose) cat("writing pdf output to file: diagnosticPlots.pdf\n")
    pdf(paste(output.dir, "/diagnosticPlots.pdf", sep=""))
    all.out<-data.frame()
    for(i in filename){
      if(verbose) cat("parsing and running sliding window for convergence\t",i,"\n")
      dat<-parseLICOR(i)

      par(mfrow=c(3,1), mar=c(5,5,5,5))

      stdplot<-twoPlot(dat=dat, trait1=trait1, trait2=trait2, title=i)
      dat$time.seconds<-stdplot$x

      sw.t1<-cvSlopeSW(stdplot$std.y1, window=window)
      sw.t2<-cvSlopeSW(stdplot$std.y2, window=window)
      lines(stdplot$x, sw.t1$mean, col="red")
      lines(stdplot$x, sw.t2$mean, col="red")

      ranks<-weightedRank(x=stdplot$x, slope1=sw.t1$slope, slope2=sw.t2$slope,
                          cv1=sw.t1$cv, cv2=sw.t2$cv, window=window,weights=weights)
      plotDiagnotics(x=stdplot$x, slope1=sw.t1$slope, slope2=sw.t2$slope,
                     cv1=sw.t1$cv, cv2=sw.t2$cv, weightedRankOutput=ranks,
                     trait1=trait1, trait2=trait2,
                     inputFile=i)
      out<-getBV(dat=dat, ranks=ranks, filename=i, window=window)
      nanames<-cnames[!cnames %in% colnames(out)]
      if(length(nanames)>0){
        for(j in nanames) out[,j]<-NA
      }
      out<-out[,c(cnames, "time.seconds","filename","window.size.observations")]
      all.out<-rbind(all.out, out)
    }
    dev.off()
    if(verbose) cat("writing dataset of values to file: meanValues.csv")
    write.csv(all.out, file=paste(output.dir, "/meanValues.csv", sep=""), row.names=F)
    return(all.out)
  }
}
