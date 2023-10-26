#!/bin/sh

#  kriging.sh
#  
#
#  Created by Martin Kapun on 18.04.18.
#  

echo '''

legend.col <- function(col, lev){
opar <- par
n <- length(col)
bx <- par("usr")
box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
box.cy <- c(bx[3], bx[3])
box.sy <- (bx[4] - bx[3]) / n
xx <- rep(box.cx, each = 2)
par(xpd = TRUE)
for(i in 1:n){
yy <- c(box.cy[1] + (box.sy * (i - 1)),
box.cy[1] + (box.sy * (i)),
box.cy[1] + (box.sy * (i)),
box.cy[1] + (box.sy * (i - 1)))
polygon(xx, yy, col = col[i], border = col[i])
}
par(new = TRUE)
plot(0, 0, type = "n",
ylim = c(min(lev), max(lev)),
yaxt = "n", ylab = "",
xaxt = "n", xlab = "",
frame.plot = FALSE)
map.axes(side = 4, las = 2, tick = FALSE, line = .25)
par <- opar
}

library(rworldmap)
library(kriging)
library(maps)
dataset=read.table("/Users/martinkapun/Documents/Work/PostDoc/manuscripts/Inversion_review_Dmel/InvFreq/worldInv.txt",header=T,na.string="NA")
#dataset=read.table("/Users/martinkapun/Documents/Work/PostDoc/manuscripts/Inversion_review_Dmel/InvFreq/EurInv.txt",header=T,na.string="NA")

dataset=na.omit(dataset)
attach(dataset)

newmap<-getMap(resolution="low")

X=c(min(Longitude),max(Longitude))
Y=c(min(Latitude),max(Latitude))
color=colorRampPalette(c("white", "yellow","red","purple","blue","darkgreen","black"),alpha=0.5)
par(cex=1.25)

pdf("/Users/martinkapun/Documents/Work/PostDoc/manuscripts/Inversion_review_Dmel/InvFreq/worldInv_In(2L)t.pdf",width=20,height=12)
Zi=-0.1#min(In.2L.t)-0.1
Za=1#max(In.2L.t)+0.1
K=kriging(Longitude,Latitude,response=In.2L.t, pixels = 500)
plot(newmap,xlim=X,ylim=Y,col=rgb(0,0,0,0.2),ylab="Latitude",xlab="Longitude")
image(K,zlim=c(Zi,Za),add=T,col=color(100))
plot(newmap,add=T,col=rgb(0,0,0,0.2))
points(Longitude,Latitude,pch=21,cex=2)
legend.col(color(100),seq(Zi,Za,0.01))
map.axes(1)
map.axes(2)
box()
dev.off()

pdf("/Users/martinkapun/Documents/Work/PostDoc/manuscripts/Inversion_review_Dmel/InvFreq/worldInv_In(3L)P.pdf",width=20,height=12)
Zi=-0.1#min(In.3L.P)-0.1
Za=1#max(In.3L.P)+0.1
K=kriging(Longitude,Latitude,response=In.3L.P, pixels = 500)
plot(newmap,xlim=X,ylim=Y,col=rgb(0,0,0,0.2),ylab="Latitude",xlab="Longitude")
image(K,zlim=c(Zi,Za),add=T,col=color(100))
plot(newmap,add=T,col=rgb(0,0,0,0.2))
points(Longitude,Latitude,pch=21,cex=2)
legend.col(color(100),seq(Zi,Za,0.01))
map.axes(1)
map.axes(2)
box()
dev.off()

pdf("/Users/martinkapun/Documents/Work/PostDoc/manuscripts/Inversion_review_Dmel/InvFreq/worldInv_In(2R)NS.pdf",width=20,height=12)
Zi=-0.1#min(In.2R.NS)-0.1
Za=1#max(In.2R.NS)+0.1
K=kriging(Longitude,Latitude,response=In.2R.NS, pixels = 500)
plot(newmap,xlim=X,ylim=Y,col=rgb(0,0,0,0.2),ylab="Latitude",xlab="Longitude")
image(K,zlim=c(Zi,Za),add=T,col=color(100))
plot(newmap,add=T,col=rgb(0,0,0,0.2))
points(Longitude,Latitude,pch=21,cex=2)
legend.col(color(100),seq(Zi,Za,0.01))
map.axes(1)
map.axes(2)
box()
dev.off()


pdf("/Users/martinkapun/Documents/Work/PostDoc/manuscripts/Inversion_review_Dmel/InvFreq/worldInv_In(3R)P.pdf",width=20,height=12)
Zi=-0.1#min(In.3R.Payne)-0.1
Za=1#max(In.3R.Payne)+0.1
K=kriging(Longitude,Latitude,response=In.3R.Payne, pixels = 500)
plot(newmap,xlim=X,ylim=Y,col=rgb(0,0,0,0.2),ylab="Latitude",xlab="Longitude")
map.axes()
image(K,zlim=c(Zi,Za),add=T,col=color(100))
plot(newmap,add=T,col=rgb(0,0,0,0.2))
points(Longitude,Latitude,pch=21,cex=2)
legend.col(color(100),seq(Zi,Za,0.01))
dev.off()

pdf("/Users/martinkapun/Documents/Work/PostDoc/manuscripts/Inversion_review_Dmel/InvFreq/worldInv_In(2L)t_dots.pdf",width=20,height=12)
Zi=-0.1#min(In.2L.t)-0.1
Za=1#max(In.2L.t)+0.1
plot(newmap,xlim=X,ylim=Y,col=rgb(0,0,0,0.2),ylab="Latitude",xlab="Longitude")
map.axes()
points(Longitude,Latitude,pch=21,cex=4,col="black",bg=color(100)[as.numeric(cut(In.2L.t,breaks = 100))])
legend.col(color(100),seq(Zi,Za,0.01))
dev.off()

pdf("/Users/martinkapun/Documents/Work/PostDoc/manuscripts/Inversion_review_Dmel/InvFreq/worldInv_In(2R)NS_dots.pdf",width=20,height=12)
Zi=-0.1#min(In.2L.t)-0.1
Za=1#max(In.2L.t)+0.1
plot(newmap,xlim=X,ylim=Y,col=rgb(0,0,0,0.2),ylab="Latitude",xlab="Longitude")
map.axes()
points(Longitude,Latitude,pch=21,cex=4,col="black",bg=color(100)[as.numeric(cut(In.2R.NS,breaks = 100))])
legend.col(color(100),seq(Zi,Za,0.01))
dev.off()

pdf("/Users/martinkapun/Documents/Work/PostDoc/manuscripts/Inversion_review_Dmel/InvFreq/worldInv_In(3L)P_dots.pdf",width=20,height=12)
Zi=-0.1#min(In.2L.t)-0.1
Za=1#max(In.2L.t)+0.1
plot(newmap,xlim=X,ylim=Y,col=rgb(0,0,0,0.2),ylab="Latitude",xlab="Longitude")
map.axes()
points(Longitude,Latitude,pch=21,cex=4,col="black",bg=color(100)[as.numeric(cut(In.3L.P,breaks = 100))])
legend.col(color(100),seq(Zi,Za,0.01))
dev.off()

pdf("/Users/martinkapun/Documents/Work/PostDoc/manuscripts/Inversion_review_Dmel/InvFreq/worldInv_In(3R)Payne_dots.pdf",width=20,height=12)
Zi=-0.1#min(In.2L.t)-0.1
Za=1#max(In.2L.t)+0.1
plot(newmap,xlim=X,ylim=Y,col=rgb(0,0,0,0.2),ylab="Latitude",xlab="Longitude")
map.axes()
points(Longitude,Latitude,pch=21,cex=4,col="black",bg=color(100)[as.numeric(cut(In.3R.Payne,breaks = 100))])
legend.col(color(100),seq(Zi,Za,0.01))
dev.off()



detach(dataset)

