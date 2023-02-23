plot.rect<-function(xc, xr, yc, yr){
  n = length(xc)
  xl <- xc-xr/2
  xu <- xc+xr/2
  yu <- yc+yr/2
  yl <- yc-yr/2
  
  plot(0,0,type="n",main="Rectangle Plot",xlab="X",ylab="Y",
       cex.lab=0.8, cex.axis=0.8, cex.main=0.8, 
       xlim=c(min(xl)-(max(xu)-min(xl))/8, max(xu)+(max(xu)-min(xl))/8), 
       ylim=c(min(yl)-(max(yu)-min(yl))/8, max(yu)+(max(yu)-min(yl))/8), 
       )
  grid()
  segments(x0=xl,y0=yl,x1=xu,y1=yl,col=c(rep("black",(n))))
  segments(x0=xl,y0=yu,x1=xu,y1=yu,col=c(rep("black",(n))))
  segments(x0=xl,y0=yl,x1=xl,y1=yu,col=c(rep("black",(n))))
  segments(x0=xu,y0=yl,x1=xu,y1=yu,col=c(rep("black",(n))))
}


plot.cross<-function(xc, xr, yc, yr){
  n = length(xc)
  xl <- xc-xr/2
  xu <- xc+xr/2
  yu <- yc+yr/2
  yl <- yc-yr/2
  
  plot(0,0,type="n",main="Cross Plot",xlab="X",ylab="Y",
       cex.lab=0.8, cex.axis=0.8, cex.main=0.8, 
       xlim=c(min(xl)-(max(xu)-min(xl))/8, max(xu)+(max(xu)-min(xl))/8), 
       ylim=c(min(yl)-(max(yu)-min(yl))/8, max(yu)+(max(yu)-min(yl))/8), 
       )
  
  grid()
  points(xc,yc,pch=16,col=c(rep("black",(length(xl)))))
  segments(x0=xl,y0=yc,x1=xu,y1=yc,lwd=1,col=c(rep("black",(length(xl)))))
  segments(x0=xc,y0=yl,x1=xc,y1=yu,lwd=1,col=c(rep("black",(length(xl)))))
}


plot.segment<-function(xc, xr, yc, yr){
  ### Standardize the data
  xc1 <- (xc - mean(xc))/sd(xc)
  yc1 <- (yc - mean(yc))/sd(yc)
  xr1 <- (xr - mean(xr))/sd(xr)
  yr1 <- (yr - mean(yr))/sd(yr)
  plot.segs <- data.frame(x1 = xc1, y1=xr1, xend =yc1, yend =yr1,colors="black ")
  plot.inter2 <- data.frame(x = c(xc1,yc1), y = c(xr1,yr1),
                            colors = c(rep("black",length(xc1)),rep("black",length(yc))),
                            shapes = c(rep(21,length(xc1)),rep(24,length(yc))),
                            bgs = c(rep("black",length(xc1)),rep("white",length(yc))))
  
  
  par(mar = c(5,5,5,8))
  par(pin=c(2.8,2.8))
  
  a0=min(plot.segs$x1)-(max(plot.segs$x1)-min(plot.segs$x1))/8
  b0=max(plot.segs$x1)+(max(plot.segs$x1)-min(plot.segs$x1))/8
  c0=min(plot.segs$y1)-(max(plot.segs$y1)-min(plot.segs$y1))/8
  d0=max(plot.segs$y1)+(max(plot.segs$y1)-min(plot.segs$y1))/8
  a1 = max(abs(c(a0,b0,c0,d0)))
  a=-1*a1
  b=a1
  c=-1*a1
  d=a1
  
  plot(0,0,col="transparent",xlim=c(a,b),ylim=c(c,d),
       main="Interval Distance",xlab="Center",ylab="Range",cex.main=1, cex.lab=1, cex.axis=1)
  grid()
  abline(v=0,lty=2)
  abline(h=0,lty=2)
  abline(0,1,lty=2)
  abline(0,-1,lty=2)
  segments(plot.segs$x1,plot.segs$y1,plot.segs$xend,plot.segs$yend,col="black",lwd=1)
  points(plot.inter2$x,plot.inter2$y,col=as.character(plot.inter2$colors),
         xlab="Center",ylab="Range",pch=plot.inter2$shapes,cex=1,bg=plot.inter2$bgs)
}

plot.dandelion<-function(xc, xr, yc, yr, polygon=TRUE){
  ### Standardize the data
  xc1 <- (xc - mean(xc))/sd(xc)
  yc1 <- (yc - mean(yc))/sd(yc)
  xr1 <- (xr - mean(xr))/sd(xr)
  yr1 <- (yr - mean(yr))/sd(yr)
  
  
  ### Move y to the original point (0,0)
  xcstar <- xc1 - yc1
  xrstar <- xr1 - yr1
  ycstar <- yc1 - yc1
  yrstar <- yr1 - yr1
  
  plot.segs <- data.frame(x1 = xcstar, y1=xrstar, xend =ycstar, yend =yrstar, colors="black")
  plot.guide<-data.frame(x1=c(-2/sqrt(pi),-2/sqrt(pi),-2/sqrt(pi),0,2/sqrt(pi),2/sqrt(pi),2/sqrt(pi),0),
                         xend=c(-2/sqrt(pi),-2/sqrt(pi),0,2/sqrt(pi),2/sqrt(pi),2/sqrt(pi),0,-2/sqrt(pi)),
                         y1=c(-2/sqrt(pi),0,2/sqrt(pi),2/sqrt(pi),2/sqrt(pi),0,-2/sqrt(pi),-2/sqrt(pi)),
                         yend=c(0,2/sqrt(pi),2/sqrt(pi),2/sqrt(pi),0,-2/sqrt(pi),-2/sqrt(pi),-2/sqrt(pi)))
  plot.true<-data.frame(x1=c(-mean(abs(xcstar+xrstar))*sqrt(2)/2,-mean(abs(xcstar)),-mean(abs(xcstar-xrstar))*sqrt(2)/2,0,mean(abs(xcstar+xrstar))*sqrt(2)/2,mean(abs(xcstar)),mean(abs(xcstar-xrstar))*sqrt(2)/2,0),
                        xend=c(-mean(abs(xcstar)),-mean(abs(xcstar-xrstar))*sqrt(2)/2,0,mean(abs(xcstar+xrstar))*sqrt(2)/2,mean(abs(xcstar)),mean(abs(xcstar-xrstar))*sqrt(2)/2,0,-mean(abs(xcstar+xrstar))*sqrt(2)/2),
                        y1=c(-mean(abs(xcstar+xrstar))*sqrt(2)/2,0,mean(abs(xcstar-xrstar))*sqrt(2)/2,mean(abs(xrstar)),mean(abs(xcstar+xrstar))*sqrt(2)/2,0,-mean(abs(xcstar-xrstar))*sqrt(2)/2,-mean(abs(xrstar))),
                        yend=c(0,mean(abs(xcstar-xrstar))*sqrt(2)/2,mean(abs(xrstar)),mean(abs(xcstar+xrstar))*sqrt(2)/2,0,-mean(abs(xcstar-xrstar))*sqrt(2)/2,-mean(abs(xrstar)),-mean(abs(xcstar+xrstar))*sqrt(2)/2))
  plot.inter2 <- data.frame(x = c(xcstar,ycstar), y = c(xrstar,yrstar),
                            colors = c(rep("black",length(xc1)),rep("black",length(yc))),
                            shapes = c(rep(21,length(xc1)),rep(24,length(yc))),
                            bgs = c(rep("black",length(xc1)),rep("white",length(yc))))
  
  
  # par(mar = c(5,5,5,8))
  par(pin=c(2.8,2.8))

  a0=min(plot.segs$x1)-(max(plot.segs$x1)-min(plot.segs$x1))/8
  b0=max(plot.segs$x1)+(max(plot.segs$x1)-min(plot.segs$x1))/8
  c0=min(plot.segs$y1)-(max(plot.segs$y1)-min(plot.segs$y1))/8
  d0=max(plot.segs$y1)+(max(plot.segs$y1)-min(plot.segs$y1))/8
  a1 = max(abs(c(a0,b0,c0,d0)))
  a=-1*a1
  b=a1
  c=-1*a1
  d=a1
  
  plot(0,0,col="transparent",xlim=c(a,b),ylim=c(c,d),main="Relative Interval Distance",xlab="Center",ylab="Range",
       cex.main=1, cex.lab=1, cex.axis=1)
  grid()
  abline(v=0,lty=2)
  abline(h=0,lty=2)
  abline(0,1,lty=2)
  abline(0,-1,lty=2)
  segments(plot.segs$x1,plot.segs$y1,plot.segs$xend,plot.segs$yend,col="black",lwd=1)
  
  points(plot.inter2$x,plot.inter2$y,col=as.character(plot.inter2$colors),
         xlab="Center",ylab="Range",pch=plot.inter2$shapes,cex=1,bg=plot.inter2$bgs)
  if(polygon==T)
  {
    segments(plot.guide$x1,plot.guide$y1,plot.guide$xend,plot.guide$yend,col="black",lwd=1,lty=2)
    segments(plot.true$x1,plot.true$y1,plot.true$xend,plot.true$yend,col="red",lwd=1)
  }
}
