#' Import zoo
#' Export
parseLICOR<-function(x, use.rtf=TRUE){
  datStart<-readLines(x)
  datStart<-grep("STARTOFDATA",datStart)
  dat<-read.delim(x, sep="\t",  header=T, skip=datStart, stringsAsFactors=F)
  dat<-dat[dat$HHMMSS!="",]
  info<-dat$Obs[5]
  return(dat)
}
hhmmss2Time<-function(x){
  x<-as.character(x)
  h<-as.numeric(sapply(x, function(y) strsplit(y, ":")[[1]][1]))
  m<-as.numeric(sapply(x, function(y) strsplit(y, ":")[[1]][2]))
  s<-as.numeric(sapply(x, function(y) strsplit(y, ":")[[1]][3]))
  d<-data.frame(h,m,s)
  d$time<-with(d, ((h*3600)+(m*60)+s))
  d$time<-round(d$time-(min(d$time)),2)
  return(d$time)
}
std.01<-function(x) {
  (x-min(x))/max(x-min(x))
}
twoPlot<-function(timename="HHMMSS",dat,trait1,trait2, title){
  y1<-dat[,trait1]
  y2<-dat[,trait2]
  x<-dat[,timename]
  x<-hhmmss2Time(x)
  std.y1<-std.01(y1)
  std.y2<-std.01(y2)
  par(mar=c(5,5,5,5))
  plot(x, std.y1, type="n",bty="n", yaxt="n", ylab="", xaxt="n", xlab="time (minutes since recording began)",
       xlim=c(0, ceiling(max(x)/60)*60), main=title)
  axis(side=1, at=seq(from=0, to=ceiling(max(x)/60)*60,by=60), labels=seq(from=0, to=ceiling(max(x)/60),by=1))
  axis(side=2, at=c(0,.2,.4,.6,.8,1),
       labels=round(as.numeric(sapply(c(0,.2,.4,.6,.8,1), function(x) quantile(y1,x))),1),
       las=1)
  mtext(trait1, side=2, line=3)
  axis(side=4, at=c(0,.2,.4,.6,.8,1),
       labels=round(as.numeric(sapply(c(0,.2,.4,.6,.8,1), function(x) quantile(y2,x))),2),
       las=1)
  mtext(trait2, side=4, line=3)
  points(x,std.y1, cex=.5)
  points(x,std.y2, cex=.8, pch=2, col="blue")
  lines(x,std.y1);lines(x,std.y2, col="blue")
  return(data.frame(x=x, std.y1=std.y1, std.y2=std.y2))
}
cvSlopeSW<-function(x, window){
  temp<-rollapply(x, width=window, fill=NA, partial=FALSE, align="center", function(y) {
    cv<-sd(y)
    ind<-1:length(y)
    mean<-mean(y)
    slope<-abs(as.numeric(coef(lm(y~ind))[2]))
    n<-length(y)
    return(c(cv,mean,slope,n))
  })
  return(data.frame(mean=temp[,2], slope=temp[,3], cv=temp[,1], n.obs=temp[,4]))
}

weightedRank<-function(x, slope1,cv1,slope2=NULL,cv2=NULL, weights=NULL,window){
  if(is.null(weights)){
    if(is.null(slope2)){
      weights<-rep(1,2)
    }else{
      weights<-rep(1,4)
    }
  }
  if(is.null(slope2)){
    summedRank<-
      (rank(slope1)*weights[1]) +
      (rank(cv1)*weights[2])
  }else{
    summedRank<-
      (rank(slope1)*weights[1]) +
      (rank(cv1)*weights[2]) +
      (rank(slope2)*weights[3]) +
      (rank(cv2)*weights[4])
  }
  summedQuant<-summedRank/max(summedRank, na.rm=TRUE)
  best.index<-which(summedQuant==min(summedQuant, na.rm=TRUE))
  best.point<-x[best.index]
  if((best.index-floor((window)/2)) != 0){
    best.low<-x[best.index-floor((window)/2)]
  }else{
    best.low<-1
  }
  if((best.index+ceiling((window)/2))<length(x)){
    best.high<-x[best.index+ceiling((window)/2)]
  }else{
    best.high<-x[length(x)]
  }

  wind.size<-x[which(summedQuant==min(summedQuant,na.rm=T))]
  rect(xleft=best.low, xright=best.high, ytop=1, ybottom=0,
       border="forestgreen", col=rgb(0,1,0,.2))

  return(list(summedRank=summedRank, summedQuant=summedQuant, wind=c(best.low,best.high),
              weights=weights, window.size=window))
}

plotDiagnotics<-function(x, slope1, slope2=NULL, cv1, cv2=NULL,weightedRankOutput, inputFile,
                         trait1, trait2){
  good<-which(!is.na(slope1))
  x<-x[good]
  slope1<-slope1[good]
  cv1<-cv1[good]
  if(!is.null(cv2)){
    slope2<-slope2[good]
    cv2<-cv2[good]
  }
  plot(x, slope1, ylim=c(0,1), type="n", bty="n",yaxt="n",xaxt="n", ylab="", xlab="time (minutes since recording began)",
       xlim=c(0, ceiling(max(x)/60)*60))
  axis(side=1, at=seq(from=0, to=ceiling(max(x)/60)*60,by=60), labels=seq(from=0, to=ceiling(max(x)/60),by=1))

  lines(x, weightedRankOutput$summedQuant[good], col=rgb(1,.5,0,1), lwd=4)
  y1<-rank(slope1)/length(slope1)
  y2<-rank(cv1)/length(cv1)
  axis(side=2, at=c(0,.2,.4,.6,.8,1),
       labels=c(0,.2,.4,.6,.8,1),
       las=1)
  mtext("Quantile of Slope/CV", side=2, line=3)
  lines(x, y1, lty=1)
  lines(x, y2, lty=2)
  if(!is.null(slope2)){
    y1<-rank(slope2)/length(slope2)
    y2<-rank(cv2)/length(cv2)
    lines(x, y1, lty=1, col="blue")
    lines(x, y2, lty=2, col="blue")
  }

  rect(xleft=weightedRankOutput$wind[1], xright=weightedRankOutput$wind[2], ytop=1, ybottom=0,
       border="forestgreen", col=rgb(0,1,0,.2))

  plot(1:6,1:6, type="n",xaxt="n",yaxt="n", ylab="",xlab="", bty="n")
  text(x=1,y=6,adj=c(0,1),"upper plot presents raw and smoothed physiology data highlighting the best window")
  text(x=1,y=5,adj=c(0,1),"lower plot presents quantiles of cv (dashed) and slope (solid).")
  text(x=1,y=4,adj=c(0,1), paste("window size =", weightedRankOutput$window.size))
  text(x=1,y=3,adj=c(0,1), paste("weights (slope1, cv1, slope2, cv2) = ", paste(weightedRankOutput$weights, collapse=", ")))
  text(x=1,y=2,adj=c(0,1), paste("input file:", inputFile))
  title(paste(trait1,": Black ...",trait2,":Blue"))
  legend("bottomright", c("mean trait",trait1, trait2,"summed, weighted quantile"), col=c("red","black","blue",rgb(1,.5,0,1)), lwd=c(1,1,1,4),
         pt.cex=c(0,1,.8,0), pch=c(1,1,2,1))
}

getBV<-function(dat, ranks, window,filename){
  temp<-dat[dat$time.seconds>=ranks$wind[1] & dat$time.seconds<=ranks$wind[2],]
  out1<-temp[ceiling(nrow(temp)/2),1:2]
  temp2<-temp[,-c(1:2)]
  out2<-c(out1,colMeans(temp2))
  out<-data.frame(t(unlist(out2)))
  out$filename=filename
  out$window.size.observations=window
  return(out)
}
