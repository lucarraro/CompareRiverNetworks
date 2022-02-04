rm(list=ls())

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(OCNet)
library(fields)
library(rnaturalearth)
library(sp)
library(rgdal)

source('support_functions.R')

load('results/data20.rda')
load('results/data100.rda')
load('results/data500.rda')
load('results/HortonLaws.rda'); for(i in 1:length(ll)) {assign(names(ll)[i], ll[[i]])}
load('results/scaling.rda'); for(i in 1:length(scaling)) {assign(names(scaling)[i], scaling[[i]])}

rivers <- read.csv("rivers/real_rivers.csv")
riverNames <- rivers$Name

riversFull <- read.csv('results/TableS1_draft.csv')

# Fig. 1c ####
pdf(file="Fig1c_draft.pdf",width=9/2.54, height=5/2.54)
set.seed(7); OCN <- create_OCN(20,10)
draw_simple_OCN(OCN, thrADraw=10)
OCN <- landscape_OCN(OCN)
OCN <- aggregate_OCN(OCN,10)
abline(h=0:10); abline(v=0:20)
points(OCN$RN$X,OCN$RN$Y,cex=2)
title(sprintf("N = %d - p = %.2f",OCN$RN$nNodes,OCN$AG$nNodes/OCN$RN$nNodes))
dev.off()

# Fig. 1dfgi ####
pdf(file="Fig1dfgi_draft.pdf", width=17/2.54, height=9/2.54)
par(mfrow=c(2,2))
load("rivers/Adige.rda")
Adige20 <- extract_river(Adige, 20)
Adige20_short <- extract_river(Adige, 20, maxReachLength = 1.99*Adige20$cellsize) 

node <- which(Adige20_short$AG$A==min(Adige20_short$AG$A))[1]
nodes <- Adige20_short$SC$toFD[node][[1]]

subset_X <- Adige20_short$FD$X[nodes]
subset_Y <- Adige20_short$FD$Y[nodes]
X_mesh <- c(min(subset_X)-Adige20$cellsize, sort(unique(subset_X)), max(subset_X)+Adige20$cellsize) 
Y_mesh <- c(min(subset_Y)-Adige20$cellsize, sort(unique(subset_Y)), max(subset_Y)+Adige20$cellsize) 
mesh <- matrix(data=0,nrow=length(Y_mesh),ncol=length(X_mesh))
for (i in 1:length(subset_X)){
  ind_x <- which(X_mesh==subset_X[i])
  ind_y <- which(Y_mesh==subset_Y[i])
  mesh[ind_y,ind_x] <- 1
}
count <- contourLines(X_mesh,Y_mesh,t(mesh),levels=1)
X_contour <- vector("list", length = length(count))
Y_contour <- vector("list", length = length(count))
for (k in 1:length(count)){
  X_contour[[k]] <- count[[k]]$x
  Y_contour[[k]] <- count[[k]]$y
} 

draw_thematic_catch(Adige20, backgroundColor = "white", addLegend=F)
for (k in 1:length(X_contour)){
  lines(X_contour[[k]],Y_contour[[k]],lwd=2,col="black")
}
title(sprintf('thrA = %.2f km2  -  N = %d  -  p = %.2f',
              Adige20$thrA*Adige$cellsize^2/1e6, 
              Adige20$RN$nNodes, Adige20$AG$nNodes/Adige20$RN$nNodes))

##
Adige500 <- extract_river(Adige, 500)
Adige500_short <- extract_river(Adige, 500, maxReachLength = 1.99*Adige500$cellsize) 

node <- which(Adige500_short$AG$A==min(Adige500_short$AG$A))[1]
nodes <- Adige500_short$SC$toFD[node][[1]]

subset_X <- Adige500_short$FD$X[nodes]
subset_Y <- Adige500_short$FD$Y[nodes]
X_mesh <- c(min(subset_X)-Adige500$cellsize, sort(unique(subset_X)), max(subset_X)+Adige500$cellsize)
Y_mesh <- c(min(subset_Y)-Adige500$cellsize, sort(unique(subset_Y)), max(subset_Y)+Adige500$cellsize) 
mesh <- matrix(data=0,nrow=length(Y_mesh),ncol=length(X_mesh))
for (i in 1:length(subset_X)){
  ind_x <- which(X_mesh==subset_X[i])
  ind_y <- which(Y_mesh==subset_Y[i])
  mesh[ind_y,ind_x] <- 1
}
count <- contourLines(X_mesh,Y_mesh,t(mesh),levels=1)
X_contour <- vector("list", length = length(count))
Y_contour <- vector("list", length = length(count))
for (k in 1:length(count)){
  X_contour[[k]] <- count[[k]]$x
  Y_contour[[k]] <- count[[k]]$y
} 

draw_thematic_catch(Adige500, backgroundColor = "white", addLegend=F)
for (k in 1:length(X_contour)){
  lines(X_contour[[k]],Y_contour[[k]],lwd=2,col="black")
}
title(sprintf('thrA = %.2f km2  -  N = %d  -  p = %.2f',
              Adige500$thrA*Adige$cellsize^2/1e6, 
              Adige500$RN$nNodes, Adige500$AG$nNodes/Adige500$RN$nNodes))

##
load("OCN/OCN31.rda")
OCN <- landscape_OCN(OCN31)
OCN_20 <- aggregate_OCN(OCN, thrA=20)
OCN_20_short <- aggregate_OCN(OCN, thrA=20, maxReachLength = 1.99)

node <- which(OCN_20_short$AG$A==min(OCN_20_short$AG$A))[1]
nodes <- OCN_20_short$SC$toFD[node][[1]]
subset_X <- OCN_20_short$FD$X[nodes]
subset_Y <- OCN_20_short$FD$Y[nodes]
X_mesh <- c(min(subset_X)-OCN$cellsize, sort(unique(subset_X)), max(subset_X)+OCN$cellsize) 
Y_mesh <- c(min(subset_Y)-OCN$cellsize, sort(unique(subset_Y)), max(subset_Y)+OCN$cellsize) 
mesh <- matrix(data=0,nrow=length(Y_mesh),ncol=length(X_mesh))
for (i in 1:length(subset_X)){
  ind_x <- which(X_mesh==subset_X[i])
  ind_y <- which(Y_mesh==subset_Y[i])
  mesh[ind_y,ind_x] <- 1
}
count <- contourLines(X_mesh,Y_mesh,t(mesh),levels=1)
X_contour <- vector("list", length = length(count))
Y_contour <- vector("list", length = length(count))
for (k in 1:length(count)){
  X_contour[[k]] <- count[[k]]$x
  Y_contour[[k]] <- count[[k]]$y
} 

draw_simple_OCN(OCN, thrA=20, easyDraw=T)
lines(OCN$CM$XContour[[1]][[1]],OCN$CM$YContour[[1]][[1]])
for (k in 1:length(X_contour)){
  lines(X_contour[[k]],Y_contour[[k]],lwd=2,col="black")
}
title(sprintf('thrA = %d cells  -  N = %d  -  p = %.2f',
              OCN_20$thrA, 
              OCN_20$RN$nNodes, OCN_20$AG$nNodes/OCN_20$RN$nNodes))

##
OCN_500 <- aggregate_OCN(OCN, thrA=500)
OCN_500_short <- aggregate_OCN(OCN, thrA=500, maxReachLength = 1.99)

node <- which(OCN_500_short$AG$A==min(OCN_500_short$AG$A))[1]
nodes <- OCN_500_short$SC$toFD[node][[1]]
subset_X <- OCN_500_short$FD$X[nodes]
subset_Y <- OCN_500_short$FD$Y[nodes]
X_mesh <- c(min(subset_X)-OCN$cellsize, sort(unique(subset_X)), max(subset_X)+OCN$cellsize) 
Y_mesh <- c(min(subset_Y)-OCN$cellsize, sort(unique(subset_Y)), max(subset_Y)+OCN$cellsize) 
mesh <- matrix(data=0,nrow=length(Y_mesh),ncol=length(X_mesh))
for (i in 1:length(subset_X)){
  ind_x <- which(X_mesh==subset_X[i])
  ind_y <- which(Y_mesh==subset_Y[i])
  mesh[ind_y,ind_x] <- 1
}
count <- contourLines(X_mesh,Y_mesh,t(mesh),levels=1)
X_contour <- vector("list", length = length(count))
Y_contour <- vector("list", length = length(count))
for (k in 1:length(count)){
  X_contour[[k]] <- count[[k]]$x
  Y_contour[[k]] <- count[[k]]$y
} 
draw_simple_OCN(OCN, thrA=500, easyDraw=T)
lines(OCN$CM$XContour[[1]][[1]],OCN$CM$YContour[[1]][[1]])
for (k in 1:length(X_contour)){
  lines(X_contour[[k]],Y_contour[[k]],lwd=2,col="black")
}
title(sprintf('thrA = %d cells  -  N = %d  -  p = %.2f',
              OCN_500$thrA, 
              OCN_500$RN$nNodes, OCN_500$AG$nNodes/OCN_500$RN$nNodes))
dev.off()

# Fig. 1eh ####
pdf(file="Fig1eh_draft.pdf", width=17/2.54, height=9/2.54)
par(mfrow=c(1,2))

Adige100 <- extract_river(Adige, 100)
Adige100_short <- extract_river(Adige, 100, maxReachLength = 1.99*Adige100$cellsize) 

node <- which(Adige100_short$AG$A==min(Adige100_short$AG$A))[1]
nodes <- Adige100_short$SC$toFD[node][[1]]

subset_X <- Adige100_short$FD$X[nodes]
subset_Y <- Adige100_short$FD$Y[nodes]
X_mesh <- c(min(subset_X)-Adige100$cellsize, sort(unique(subset_X)), max(subset_X)+Adige100$cellsize)
Y_mesh <- c(min(subset_Y)-Adige100$cellsize, sort(unique(subset_Y)), max(subset_Y)+Adige100$cellsize) 
mesh <- matrix(data=0,nrow=length(Y_mesh),ncol=length(X_mesh))
for (i in 1:length(subset_X)){
  ind_x <- which(X_mesh==subset_X[i])
  ind_y <- which(Y_mesh==subset_Y[i])
  mesh[ind_y,ind_x] <- 1
}
count <- contourLines(X_mesh,Y_mesh,t(mesh),levels=1)
X_contour <- vector("list", length = length(count))
Y_contour <- vector("list", length = length(count))
for (k in 1:length(count)){
  X_contour[[k]] <- count[[k]]$x
  Y_contour[[k]] <- count[[k]]$y
} 

draw_thematic_catch(Adige100, backgroundColor = "white", addLegend=F)
for (k in 1:length(X_contour)){
  lines(X_contour[[k]],Y_contour[[k]],lwd=2,col="black")
}

OCN_100_short <- aggregate_OCN(OCN, thrA=100, maxReachLength = 1.99)

node <- which(OCN_100_short$AG$A==min(OCN_100_short$AG$A))[1]
nodes <- OCN_100_short$SC$toFD[node][[1]]
subset_X <- OCN_100_short$FD$X[nodes]
subset_Y <- OCN_100_short$FD$Y[nodes]
X_mesh <- c(min(subset_X)-OCN$cellsize, sort(unique(subset_X)), max(subset_X)+OCN$cellsize) 
Y_mesh <- c(min(subset_Y)-OCN$cellsize, sort(unique(subset_Y)), max(subset_Y)+OCN$cellsize) 
mesh <- matrix(data=0,nrow=length(Y_mesh),ncol=length(X_mesh))
for (i in 1:length(subset_X)){
  ind_x <- which(X_mesh==subset_X[i])
  ind_y <- which(Y_mesh==subset_Y[i])
  mesh[ind_y,ind_x] <- 1
}
count <- contourLines(X_mesh,Y_mesh,t(mesh),levels=1)
X_contour <- vector("list", length = length(count))
Y_contour <- vector("list", length = length(count))
for (k in 1:length(count)){
  X_contour[[k]] <- count[[k]]$x
  Y_contour[[k]] <- count[[k]]$y
} 

draw_simple_OCN(OCN, thrA=100, easyDraw=T)
lines(OCN$CM$XContour[[1]][[1]],OCN$CM$YContour[[1]][[1]])
for (k in 1:length(X_contour)){
  lines(X_contour[[k]],Y_contour[[k]],lwd=2,col="black")
}

dev.off()

# Fig. 2 ####
pdf('Fig2_draft.pdf',width=17/2.54,height=10/2.54)
par(mfrow=c(1,2))
thrA <- exp(log(10)*seq(0, 4, length.out=100))
A <- exp(log(10)*seq(4, 7, length.out=100))
thr_a <- thrA %*% (t(A)^(-1))

N <- 0.435*thrA^(-0.446) %*% t(A)
p <- 1.531*thrA^(-0.523) %*% t(A^(-0.032))
p[p>1] <- 1

image.plot(log10(thrA),log10(A),log10(N),
           xlab='AT',ylab='A',yaxt="n",col=tim.colors(1000),
           breaks=seq(1,7,length.out=1001))
axis(side=2,at=4:7)
contour(log10(thrA),log10(A),log10(thr_a),add=T, levels=-(6:1))
title('N')
points(log10(c(20,100,500)),log10(rep(4e4,3)),pch=19)


image.plot(log10(thrA),log10(A),p,
           xlab='AT',ylab='A',yaxt="n",col=tim.colors(1000),
           breaks=seq(0,1,length.out=1001))
axis(side=2,at=4:7)
contour(log10(thrA),log10(A),log10(thr_a),add=T, levels=-(6:1))
title('p')
points(log10(c(20,100,500)),log10(rep(4e4,3)),pch=19)
dev.off()

# Fig. 3 ####
pdf(file="Fig3_draft.pdf",width=17/2.54, height=6/2.54)
p20 <- data20$nLinks$Rivers/data20$nNodes$Rivers
p100 <- data100$nLinks$Rivers/data100$nNodes$Rivers
p500 <- data500$nLinks$Rivers/data500$nNodes$Rivers
aa <- cbind(p20, p100, p500)
p20_n <- (p20-mean(p20))/sd(p20)
p100_n <- (p100-mean(p100))/sd(p100)
p500_n <- (p500-mean(p500))/sd(p500)
aa_n <- cbind(p20_n, p100_n, p500_n)
specialRivers <- c(6,7,13,14,22,40,44,45,46,48,49)
cols <- rainbow(11)
par(mfrow=c(1,2),mai=c(0,0,0,0))
plot(c(20,100,500),aa[1,],type="o",col="gray",log="xy",pch=19,xlim=c(100/6,500),ylim=c(0.01,0.5),
     bty="n",xaxt="n",yaxt="n",xlab="A_T [no. pixels]",ylab="p_r")
axis(1,pos=0.01,at=c(20,100,500)); axis(2,pos=100/6)
for(i in 2:50){if (!(i %in% specialRivers)){lines(c(20,100,500),aa[i,],type="o",pch=19,col="gray")}}
for (i in 1:length(specialRivers)){lines(c(20,100,500),aa[specialRivers[i],],lwd=1,type="o",pch=19,col=cols[i])}
legend(200,0.5,legend=riverNames[specialRivers],col=cols,lty=1+numeric(11))

plot(c(20,100,500),aa_n[1,],type="o",col="gray",log="x",pch=19,ylim=c(-3,4.5),xlim=c(100/6,500),
     bty="n",xaxt="n",yaxt="n",xlab="A_T [no. pixels]",ylab="Normalized p_r")
axis(1,pos=-3,at=c(20,100,500)); axis(2,pos=100/6,at=seq(-3,4.5,1.5))
for(i in 2:50){if (!(i %in% specialRivers)){lines(c(20,100,500),aa_n[i,],type="o",pch=19,col="gray")}}
for (i in 1:length(specialRivers)){lines(c(20,100,500),aa_n[specialRivers[i],],lwd=1,type="o",pch=19,col=cols[i])}
#legend(20,4.5,legend=riverNames[specialRivers],col=cols,lty=1+numeric(11))
dev.off()

# Fig. 4a ####
pdf(file="Fig4a_draft.pdf",width=17/2.54, height=10/2.54)
plot(ne_countries(scale="medium",continent="Europe"),xlim=c(3,18),ylim=c(43,48.6),col="#EEEEEE")

riverContinent <- rivers$continent

for (i in which(riverContinent=="Europe")){
  load(paste0('rivers/',riverNames[i],'.rda'))
  eval(parse(text=paste0('river <- ',riverNames[i])))
  eval(parse(text=paste0('rm(',riverNames[i],')')))
  p = Polygon(cbind(river$CM$XContour,river$CM$YContour))
  ps = Polygons(list(p),1)
  xy = SpatialPolygons(list(ps), proj4string=CRS(paste0("+init=epsg:",rivers$EPSG[i])))
  XY <- spTransform(xy, CRS("+init=epsg:4326"))
  plot(XY,add=T,col="#AAAAAA",lwd=0.1)
  text(XY@bbox[1], XY@bbox[2], labels=sprintf('%d. %s',i,riverNames[i]))
}

lines(c(5,5),c(30,60),col="gray")
lines(c(10,10),c(30,60),col="gray")
lines(c(15,15),c(30,60),col="gray")
lines(c(0,20),c(45,45),col="gray")
lines(c(0,20),c(43,43),col="gray")
lines(c(0,20),c(47,47),col="gray")
lines(c(0,20),c(49,49),col="gray")
dev.off()

# Fig. 4b ####
pdf(file="Fig4b_draft.pdf",width=17/2.54, height=10/2.54)
plot(ne_countries(scale="medium",continent="North America"),xlim=c(-130,-70),ylim=c(30,70),col="#EEEEEE")
for (i in which(riverContinent=="NorthAmerica")){
  load(paste0('rivers/',riverNames[i],'.rda'))
  eval(parse(text=paste0('river <- ',riverNames[i])))
  eval(parse(text=paste0('rm(',riverNames[i],')')))
  p = Polygon(cbind(river$CM$XContour,river$CM$YContour))
  ps = Polygons(list(p),1)
  xy = SpatialPolygons(list(ps), proj4string=CRS(paste0("+init=epsg:",rivers$EPSG[i])))
  XY <- spTransform(xy, CRS("+init=epsg:4326"))
  plot(XY,add=T,col="#AAAAAA",lwd=0.1)
  text(XY@bbox[1], XY@bbox[2], labels=sprintf('%d. %s',i,riverNames[i]))
}
lines(c(5,5),c(30,60),col="gray")
lines(c(-80,-80),c(30,70),col="gray")
lines(c(-100,-100),c(30,70),col="gray")
lines(c(-120,-120),c(30,70),col="gray")
lines(c(-150,-70),c(50,50),col="gray")
lines(c(-150,-70),c(40,40),col="gray")
lines(c(-150,-70),c(60,60),col="gray")
lines(c(-150,-70),c(70,70),col="gray")
dev.off()


# Fig. 5ab ####
pdf(file="Fig5ab_draft.pdf", width=21/2.54, height=9/2.54)

cols <- hcl.colors(4,"Spectral",rev=T)
par(mfrow=c(1,3))
# N_omega
SS <- c( sO_thrA20$BBT, sO_thrA20$RBN, sO_thrA20$OCN,  sO_thrA20$Rivers)
NN <- c(N_omega_thrA20$BBT, N_omega_thrA20$RBN, N_omega_thrA20$OCN, N_omega_thrA20$Rivers)
VV <- c(rep("A_BBT",length(sO_thrA20$BBT)), rep("B_RBN",length(sO_thrA20$RBN)), 
        rep("C_OCN",length(sO_thrA20$OCN)), rep("D_Rivers",length(sO_thrA20$Rivers)) )
boxplot(NN ~ VV*SS,log="y", col=cols, frame=F,  xlab="Stream Order",ylab="N_omega",ylim=c(1,1000),
        xaxt="n",yaxt="n"); title("thrA = 20")
axis(1, pos=1, at = seq(2.5, 31, 4),labels=1:8)
axis(2, pos=0.5, at=c(1,10,100,1000))
for(i in seq(4.5 , 36 , 4)){ 
  abline(v=i,lty=1, col="grey")
}
lmod <- lm(log(N_omega_thrA20$BBT) ~ sO_thrA20$BBT)
text(28,240,labels=sprintf("BBT: RB = %.2f + %.2f",mean(RB_thrA20$BBT),sd(RB_thrA20$BBT)))
lines(seq(1,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[1])
lmod <- lm(log(N_omega_thrA20$RBN) ~ sO_thrA20$RBN)
text(28,300,labels=sprintf("RBN: RB = %.2f + %.2f",mean(RB_thrA20$RBN),sd(RB_thrA20$RBN)))
lines(seq(2,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[2])
lmod <- lm(log(N_omega_thrA20$OCN) ~ sO_thrA20$OCN)
text(28,380,labels=sprintf("OCN: RB = %.2f + %.2f",mean(RB_thrA20$OCN),sd(RB_thrA20$OCN)))
lines(seq(3,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[3])
lmod <- lm(log(N_omega_thrA20$Rivers) ~ sO_thrA20$Rivers)
text(28,500,labels=sprintf("Rivers: RB = %.2f  + %.2f",mean(RB_thrA20$Rivers),sd(RB_thrA20$Rivers))) # exp(-lmod$coefficients[2])
lines(seq(4,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[4])

# L_omega
SS <- c(sO_thrA20$BBT, sO_thrA20$RBN, sO_thrA20$OCN, sO_thrA20$Rivers)
LL <- c(L_omega_thrA20$BBT, L_omega_thrA20$RBN,  L_omega_thrA20$OCN, L_omega_thrA20$Rivers)
LL[LL==0] <- NA
VV <- c(rep("A_BBT",length(sO_thrA20$BBT)), rep("B_RBN",length(sO_thrA20$RBN)),  
        rep("C_OCN",length(sO_thrA20$OCN)), rep("D_Rivers",length(sO_thrA20$Rivers)) )
boxplot(LL ~ VV*SS,log="y", col=cols,ylim=c(1,1000), xlab="Stream Order",ylab="L_omega",
        frame=F, xaxt="n", yaxt="n"); title("thrA = 20")
axis(1,pos=1,  at = seq(2.5, 31, 4),labels=1:8)
axis(2, pos=0.5, at=c(1,10,100,1000))
for(i in seq(4.5 , 36 , 4)){ 
  abline(v=i,lty=1, col="grey")
}
tmp <- which(sO_thrA20$BBT==max(sO_thrA20$BBT))
lmod <- lm(log(L_omega_thrA20$BBT[-tmp]) ~ sO_thrA20$BBT[-tmp])
text(5,240,labels=sprintf("BBT: RL = %.2f + %.2f",mean(RL_thrA20$BBT),sd(RL_thrA20$BBT)))
lines(seq(1,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[1])
tmp <- which(sO_thrA20$RBN==max(sO_thrA20$RBN))
lmod <- lm(log(L_omega_thrA20$RBN[-tmp]) ~ sO_thrA20$RBN[-tmp])
text(5,300,labels=sprintf("RBN: RL = %.2f + %.2f",mean(RL_thrA20$RBN),sd(RL_thrA20$RBN)))
lines(seq(2,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[2])
tmp <- which(sO_thrA20$OCN==max(sO_thrA20$OCN))
lmod <- lm(log(L_omega_thrA20$OCN[-tmp]) ~ sO_thrA20$OCN[-tmp])
text(5,380,labels=sprintf("OCN: RL = %.2f + %.2f",mean(RL_thrA20$OCN),sd(RL_thrA20$OCN)))
lines(seq(3,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[3])
tmp <- which(sO_thrA20$Rivers==max(sO_thrA20$Rivers))
lmod <- lm(log(L_omega_thrA20$Rivers[-tmp]) ~ sO_thrA20$Rivers[-tmp])
text(5,500,labels=sprintf("Rivers: RL = %.2f + %.2f",mean(RL_thrA20$Rivers),sd(RL_thrA20$Rivers)))
lines(seq(4,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[4])
dev.off()


# Fig. 5c ####
pdf(file="Fig5c_draft.pdf", width=21/2.54, height=9/2.54)
cols <- hcl.colors(4,"Spectral")
par(mfrow=c(1,3))
thr <- ceiling(0.05*length(P_A_OCN))
vv <- c(1:1000,seq(1010,1e4,10),seq(10100,sum(P_A_river>1e-3),100))
plot(vv,P_A_river[vv],pch=1,col=cols[1],cex=1,log="xy",ylim=c(1e-3,1),xaxt="n",yaxt="n",bty="n",
     xlab="a [no. pixels]",ylab="P [A > a]"); #title('Scaling of areas')
axis(1,pos=1e-3, at=c(1,10,1e2,1e3,1e4,1e5))
axis(2,pos=1, at=c(1,1e-3,1e-2,1e-1,1))
vv <- c(1:1000,seq(1010,1e4,10),seq(10100,sum(P_A_OCN>1e-3),100))
points(vv,P_A_OCN[vv],pch=1,col=cols[2],cex=0.75)
vv <- c(1:1000,seq(1010,1e4,10),seq(10100,sum(P_A_RBN>1e-3),100))
points(vv,P_A_RBN[vv],pch=1,col=cols[3],cex=0.5)
vv <- c(1:1000,seq(1010,sum(P_A_BBT>1e-3),10))
points(vv,P_A_BBT[vv],pch=1,col=cols[4],cex=0.25)
lmod_BBT_1 <- lm(log(P_A_BBT[11:thr]) ~ log(11:thr))
lmod_RBN_1 <- lm(log(P_A_RBN[11:thr]) ~ log(11:thr))
lmod_OCN_1 <- lm(log(P_A_OCN[11:thr]) ~ log(11:thr))
lmod_river_1 <- lm(log(P_A_river[11:thr]) ~ log(11:thr))
summary(lmod_BBT_1); summary(lmod_RBN_1); summary(lmod_OCN_1); summary(lmod_river_1)
lines(c(1,4e4),exp(lmod_river_1$coefficients[1])*c(1,4e4)^lmod_river_1$coefficients[2], col=cols[1])
lines(c(1,4e4),exp(lmod_OCN_1$coefficients[1])*c(1,4e4)^lmod_OCN_1$coefficients[2], col=cols[2])
lines(c(1,4e4),exp(lmod_RBN_1$coefficients[1])*c(1,4e4)^lmod_RBN_1$coefficients[2], col=cols[3])
lines(c(2,4e4),exp(lmod_BBT_1$coefficients[1])*c(2,4e4)^lmod_BBT_1$coefficients[2], col=cols[4], ylim=c(1e-3, 1))
lines(c(thr,thr),c(1e-3,1))
text(1e4,0.9,labels = paste0(sprintf('River: beta = %.2f',lmod_river_1$coefficients[2])))
text(1e4,0.7,labels = paste0(sprintf('OCN: beta = %.2f',lmod_OCN_1$coefficients[2])))
text(1e4,0.55,labels = paste0(sprintf('RBN: beta = %.2f',lmod_RBN_1$coefficients[2])))
text(1e4,0.45,labels = paste0(sprintf('BBT: beta = %.2f',lmod_BBT_1$coefficients[2])))
dev.off()

# Fig. 6 ####
pdf(file="Fig6_draft.pdf",width=17/2.54, height=24/2.54)
par(mfrow=c(4,3))
boxplot(data20$CVm_100[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T), ylim=c(0.37,0.8)); title("CVm")
boxplot(data100$CVm_100[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T), ylim=c(0.37,0.8)); title("CVm")
boxplot(data500$CVm_100[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T), ylim=c(0.37,0.8)); title("CVm")
boxplot(data20$lambdaM_100[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity") 
boxplot(data100$lambdaM_100[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity") 
boxplot(data500$lambdaM_100[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity")

boxplot(data20$CVm_100_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T), ylim=c(0.4,0.85)); title("CVm H")
boxplot(data100$CVm_100_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T), ylim=c(0.4,0.85)); title("CVm H")
boxplot(data500$CVm_100_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T), ylim=c(0.4,0.85)); title("CVm H")
boxplot(data20$lambdaM_100_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity H") 
boxplot(data100$lambdaM_100_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity H") 
boxplot(data500$lambdaM_100_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity H")
dev.off()


# Table S2 ####
# scaling for single real rivers
thr <- 2000
beta_exp <- numeric(length(riverNames))

for (i in 1:length(riverNames)){
  fnam <- riverNames[i]
  eval(parse(text=paste0('load("rivers/',fnam,'.rda")')))
  eval(parse(text=paste0('river <- ',fnam)))
  eval(parse(text=paste0('rm(',fnam,')')))
  dA_river <- river$FD$A[!is.na(river$FD$A)]
  P_A_river <- numeric(max(dA_river))
  for (j in 1:length(P_A_river)){
    P_A_river[j] <- sum(dA_river >= j)/length(dA_river)
  }
  lmod_river_1 <- lm(log(P_A_river[11:thr]) ~ log(11:thr)); 
  lmod_river_1 <- summary(lmod_river_1);
  beta_exp[i] <- -lmod_river_1$coefficients[2,1]
}

# Horton for single rivers
j <- 0; RB <- numeric(length(riverNames))
RL <- numeric(length(riverNames))
for (i in 1:length(riverNames)){
  j <- j+1; tmp <- j
  while (sO_thrA20$Rivers[j+1]>sO_thrA20$Rivers[j]){
    j <- j+1
    tmp <- c(tmp,j)
    if (j==length(sO_thrA20$Rivers)){
      break
    }
  }
  N_omega <- N_omega_thrA20$Rivers[tmp]
  L_omega <- L_omega_thrA20$Rivers[tmp[1:(length(tmp)-1)]]
  sO <- sO_thrA20$Rivers[tmp]
  sO_L <- sO_thrA20$Rivers[tmp[1:(length(tmp)-1)]]
  lmod_RB <- lm(log(N_omega) ~ sO)
  lmod_RL <- lm(log(L_omega) ~ sO_L)
  RB[i] <- exp(-lmod_RB$coefficients[2])
  RL[i] <- exp(lmod_RL$coefficients[2])
}

df <- data.frame(river=riverNames, RB=RB,RL=RL,beta_exp=beta_exp)
write.csv(df,file='results/TableS2_draft.csv')

# Fig. S1 ####
pdf(file="FigS1_draft.pdf", width=21/2.54, height=9/2.54)
par(mfrow=c(1,3))
# N_omega
SS <- c(sO_thrA100$BBT, sO_thrA100$RBN, sO_thrA100$OCN, sO_thrA100$Rivers)
NN <- c(N_omega_thrA100$BBT, N_omega_thrA100$RBN,   N_omega_thrA100$OCN, N_omega_thrA100$Rivers)
VV <- c(rep("A_BBT",length(sO_thrA100$BBT)),  rep("B_RBN",length(sO_thrA100$RBN)), 
        rep("C_OCN",length(sO_thrA100$OCN)), rep("D_Rivers",length(sO_thrA100$Rivers)))
boxplot(NN ~ VV*SS,log="y", col=cols,  xlab="Stream Order",ylab="N_omega",
        frame=F, xaxt="n",yaxt="n"); title("thrA = 100")
axis(1, pos=1, at = seq(2.5, 31, 4),labels=1:8)
axis(2, pos=0.5, at=c(1,10,100,1000))
for(i in seq(4.5 , 36 , 4)){ 
  abline(v=i,lty=1, col="grey")
}
lmod <- lm(log(N_omega_thrA100$BBT) ~ sO_thrA100$BBT)
text(20,100,labels=sprintf("BBT: RB = %.2f + %.2f",mean(RB_thrA100$BBT),sd(RB_thrA100$BBT)))
lines(seq(1,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[1])
lmod <- lm(log(N_omega_thrA100$RBN) ~ sO_thrA100$RBN)
text(20,85,labels=sprintf("RBN: RB = %.2f + %.2f",mean(RB_thrA100$RBN),sd(RB_thrA100$RBN)))
lines(seq(2,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[2])
lmod <- lm(log(N_omega_thrA100$OCN) ~ sO_thrA100$OCN)
text(20,72,labels=sprintf("OCN: RB = %.2f + %.2f",mean(RB_thrA100$OCN),sd(RB_thrA100$OCN)))
lines(seq(3,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[3])
lmod <- lm(log(N_omega_thrA100$Rivers) ~ sO_thrA100$Rivers)
text(20,60,labels=sprintf("Rivers: RB = %.2f + %.2f",mean(RB_thrA100$Rivers),sd(RB_thrA100$Rivers)))
lines(seq(4,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[4])

# L_omega
SS <- c(sO_thrA100$BBT, sO_thrA100$RBN, sO_thrA100$OCN, sO_thrA100$Rivers)
LL <- c(L_omega_thrA100$BBT, L_omega_thrA100$RBN, L_omega_thrA100$OCN, L_omega_thrA100$Rivers)
LL[LL==0] <- NA
VV <- c(rep("A_BBT",length(sO_thrA100$BBT)),  rep("B_RBN",length(sO_thrA100$RBN)),
        rep("C_OCN",length(sO_thrA100$OCN)), rep("D_Rivers",length(sO_thrA100$Rivers)))
boxplot(LL ~ VV*SS,log="y", col=cols,ylim=c(1,500), xlab="Stream Order",ylab="L_omega",
        frame=F,xaxt="n",yaxt="n"); title("thrA = 100")
axis(1, pos=1, at = seq(2.5, 31, 4),labels=1:8)
axis(2, pos=0.5, at=c(1,10,100,1000))
for(i in seq(4.5 , 36 , 4)){ 
  abline(v=i,lty=1, col="grey")
}
tmp <- which(sO_thrA100$BBT==max(sO_thrA100$BBT))
lmod <- lm(log(L_omega_thrA100$BBT[-tmp]) ~ sO_thrA100$BBT[-tmp])
text(5,500,labels=sprintf("BBT: RL = %.2f + %.2f",mean(RL_thrA100$BBT),sd(RL_thrA100$BBT)))
lines(seq(1,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[1])
tmp <- which(sO_thrA100$RBN==max(sO_thrA100$RBN))
lmod <- lm(log(L_omega_thrA100$RBN[-tmp]) ~ sO_thrA100$RBN[-tmp])
text(5,380,labels=sprintf("RBN: RL = %.2f + %.2f",mean(RL_thrA100$RBN),sd(RL_thrA100$RBN)))
lines(seq(2,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[2])
tmp <- which(sO_thrA100$OCN==max(sO_thrA100$OCN))
lmod <- lm(log(L_omega_thrA100$OCN[-tmp]) ~ sO_thrA100$OCN[-tmp])
text(5,300,labels=sprintf("OCN: RL = %.2f + %.2f",mean(RL_thrA100$OCN),sd(RB_thrA100$OCN)))
lines(seq(3,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[3])
tmp <- which(sO_thrA100$Rivers==max(sO_thrA100$Rivers))
lmod <- lm(log(L_omega_thrA100$Rivers[-tmp]) ~ sO_thrA100$Rivers[-tmp])
text(5,240,labels=sprintf("Rivers: RL = %.2f + %.2f",mean(RL_thrA100$Rivers),sd(RL_thrA100$Rivers)))
lines(seq(4,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[4])
dev.off()


# Fig. S2 ####
pdf(file="FigS2_draft.pdf", width=21/2.54, height=9/2.54)
par(mfrow=c(1,3))
# N_omega
SS <- c(sO_thrA500$BBT, sO_thrA500$RBN, sO_thrA500$OCN, sO_thrA500$Rivers)
NN <- c(N_omega_thrA500$BBT, N_omega_thrA500$RBN, N_omega_thrA500$OCN,  N_omega_thrA500$Rivers)
VV <- c(rep("A_BBT",length(sO_thrA500$BBT)), rep("B_RBN",length(sO_thrA500$RBN)), 
        rep("C_OCN",length(sO_thrA500$OCN)), rep("D_Rivers",length(sO_thrA500$Rivers)))
boxplot(NN ~ VV*SS,log="y", col=cols, ylim=c(1,100), xlab="Stream Order",ylab="N_omega",
        frame=F, xaxt="n",yaxt="n"); title("thrA = 500")
axis(1, pos=1, at = seq(2.5, 31, 4),labels=1:8)
axis(2, pos=0.5, at=c(1,10,100,1000))
for(i in seq(4.5 , 36 , 4)){ 
  abline(v=i,lty=1, col="grey")
}
lmod <- lm(log(N_omega_thrA500$BBT) ~ sO_thrA500$BBT)
text(14,40,labels=sprintf("BBT: RB = %.2f + %.2f",mean(RB_thrA500$BBT),sd(RB_thrA500$BBT)))
lines(seq(1,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[1])
lmod <- lm(log(N_omega_thrA500$RBN) ~ sO_thrA500$RBN)
text(14,34.5,labels=sprintf("RBN: RB = %.2f + %.2f",mean(RB_thrA500$RBN),sd(RB_thrA500$RBN)))
lines(seq(2,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[2])
lmod <- lm(log(N_omega_thrA500$OCN) ~ sO_thrA500$OCN)
text(14,30,labels=sprintf("OCN: RB = %.2f + %.2f",mean(RB_thrA500$OCN),sd(RB_thrA500$OCN)))
lines(seq(3,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[3])
lmod <- lm(log(N_omega_thrA500$Rivers) ~ sO_thrA500$Rivers)
text(14,26,labels=sprintf("Rivers: RB = %.2f + %.2f",mean(RB_thrA500$Rivers),sd(RB_thrA500$Rivers)))
lines(seq(4,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[4])

# L_omega
SS <- c(sO_thrA500$BBT, sO_thrA500$RBN, sO_thrA500$OCN,  sO_thrA500$Rivers)
LL <- c(L_omega_thrA500$BBT,L_omega_thrA500$RBN, L_omega_thrA500$OCN, L_omega_thrA500$Rivers)
LL[LL==0] <- NA
VV <- c(rep("A_BBT",length(sO_thrA500$BBT)), rep("B_RBN",length(sO_thrA500$RBN)),  
        rep("C_OCN",length(sO_thrA500$OCN)), rep("D_Rivers",length(sO_thrA500$Rivers)) )
boxplot(LL ~ VV*SS,log="y", col=cols,ylim=c(1,1000),  xlab="Stream Order",ylab="L_omega",
        frame=F,xaxt="n",yaxt="n"); title("thrA = 500")
axis(1, pos=1, at = seq(2.5, 31, 4),labels=1:8)
axis(2, pos=0.5, at=c(1,10,100,1000))
for(i in seq(4.5 , 36 , 4)){ 
  abline(v=i,lty=1, col="grey")
}
tmp <- which(sO_thrA500$BBT==max(sO_thrA500$BBT))
lmod <- lm(log(L_omega_thrA500$BBT[-tmp]) ~ sO_thrA500$BBT[-tmp])
text(5,500,labels=sprintf("BBT: RL = %.2f + %.2f",mean(RL_thrA500$BBT),sd(RL_thrA500$BBT)))
lines(seq(1,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[1])
tmp <- which(sO_thrA500$RBN==max(sO_thrA500$RBN))
lmod <- lm(log(L_omega_thrA500$RBN[-tmp]) ~ sO_thrA500$RBN[-tmp])
text(5,380,labels=sprintf("RBN: RL = %.2f + %.2f",mean(RL_thrA500$RBN),sd(RL_thrA500$RBN)))
lines(seq(2,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[2])
tmp <- which(sO_thrA500$OCN==max(sO_thrA500$OCN))
lmod <- lm(log(L_omega_thrA500$OCN[-tmp]) ~ sO_thrA500$OCN[-tmp])
text(5,300,labels=sprintf("OCN: RL = %.2f + %.2f",mean(RL_thrA500$OCN),sd(RL_thrA500$OCN)))
lines(seq(3,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[3])
tmp <- which(sO_thrA500$Rivers==max(sO_thrA500$Rivers))
lmod <- lm(log(L_omega_thrA500$Rivers[-tmp]) ~ sO_thrA500$Rivers[-tmp])
text(5,240,labels=sprintf("Rivers: RL = %.2f + %.2f",mean(RL_thrA500$Rivers),sd(RL_thrA500$Rivers)))
lines(seq(4,40,4), exp(lmod$coefficients[1])*exp((1:10)*lmod$coefficients[2]), col=cols[4])
dev.off()

# Fig. S3 ####
pdf(file="FigS3_draft.pdf",width=9/2.54, height=8/2.54)
plot(riversFull$A_cells, df$beta_exp, pch=19, cex=0.5, xlim=c(32000,44000),ylim=c(0.35,0.6),
     xaxt="n", yaxt="n", bty="n",xlab="A [no. pixels]", ylab="beta")
axis(1,pos=0.35); axis(2,pos=32000)
lmod <- lm(df$beta_exp ~ riversFull$A_cells); lmod <- summary(lmod)
lines(c(32000,44000), lmod$coefficients[1]+lmod$coefficients[2]*c(32000,44000))
text(42000,0.59,sprintf('p: %.3f',lmod$coefficients[2,4]))
text(42000,0.57,sprintf('R^2: %.3f',lmod$adj.r.squared))
dev.off()


# Fig. S4 ####
pdf(file="FigS4_draft.pdf",width=17/2.54, height=24/2.54)
par(mfrow=c(4,3))
boxplot(data20$CVm_10[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.07,0.17)); title("CVm")
boxplot(data100$CVm_10[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.07,0.17)); title("CVm")
boxplot(data500$CVm_10[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.07,0.17)); title("CVm")
boxplot(data20$lambdaM_10[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity") 
boxplot(data100$lambdaM_10[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity") 
boxplot(data500$lambdaM_10[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity") 
boxplot(data20$CVm_10_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.11,0.24)); title("CVm H")
boxplot(data100$CVm_10_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.11,0.24)); title("CVm H")
boxplot(data500$CVm_10_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.11,0.24)); title("CVm H")
boxplot(data20$lambdaM_10_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity H") 
boxplot(data100$lambdaM_10_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity H") 
boxplot(data500$lambdaM_10_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity H") 
dev.off()

# Fig. S5 ####
pdf(file="FigS5_draft.pdf",width=17/2.54, height=24/2.54)
par(mfrow=c(4,3))
boxplot(data20$CVm_20[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.12,0.34)); title("CVm")
boxplot(data100$CVm_20[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.12,0.34)); title("CVm")
boxplot(data500$CVm_20[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.12,0.34)); title("CVm")
boxplot(data20$lambdaM_20[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity") 
boxplot(data100$lambdaM_20[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity") 
boxplot(data500$lambdaM_20[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity") 
boxplot(data20$CVm_20_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.17,0.42)); title("CVm H")
boxplot(data100$CVm_20_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.17,0.42)); title("CVm H")
boxplot(data500$CVm_20_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.17,0.42)); title("CVm H")
boxplot(data20$lambdaM_20_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity H") 
boxplot(data100$lambdaM_20_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity H") 
boxplot(data500$lambdaM_20_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity H")
dev.off()

# Fig. S6 ####
pdf(file="FigS6_draft.pdf",width=17/2.54, height=24/2.54)
par(mfrow=c(4,3))
boxplot(data20$CVm_200[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.50,0.89)); title("CVm")
boxplot(data100$CVm_200[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.50,0.89)); title("CVm")
boxplot(data500$CVm_200[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.50,0.89)); title("CVm")
boxplot(data20$lambdaM_200[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral")); title("metapop capacity") 
boxplot(data100$lambdaM_200[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity") 
boxplot(data500$lambdaM_200[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity")
boxplot(data20$CVm_200_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.59,0.91)); title("CVm H")
boxplot(data100$CVm_200_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.59,0.91)); title("CVm H")
boxplot(data500$CVm_200_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.59,0.91)); title("CVm H")
boxplot(data20$lambdaM_200_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity H") 
boxplot(data100$lambdaM_200_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity H") 
boxplot(data500$lambdaM_200_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity H")
dev.off()

# Fig. S7 ####
pdf(file="FigS7_draft.pdf",width=17/2.54, height=24/2.54)
par(mfrow=c(4,3))
boxplot(data20$CVm_1000[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.85,0.98)); title("CVm")
boxplot(data100$CVm_1000[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.85,0.98)); title("CVm")
boxplot(data500$CVm_1000[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.85,0.98)); title("CVm")
boxplot(data20$lambdaM_1000[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity")
boxplot(data100$lambdaM_1000[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity")
boxplot(data500$lambdaM_1000[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity")
boxplot(data20$CVm_1000_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.88,1)); title("CVm H")
boxplot(data100$CVm_1000_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.88,1)); title("CVm H")
boxplot(data500$CVm_1000_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(0.88,1)); title("CVm H")
boxplot(data20$lambdaM_1000_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity H") 
boxplot(data100$lambdaM_1000_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity H") 
boxplot(data500$lambdaM_1000_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("metapop capacity H")
dev.off()



# Fig. S8 ####
pdf(file="FigS8_draft.pdf",width=17/2.54, height=12/2.54)
par(mfrow=c(2,3))
boxplot(data20$nNodes[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("N")
abline(h=0.435*20^(-0.446) * 40000)
boxplot(data100$nNodes[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("N")
abline(h=0.435*100^(-0.446) * 40000)
boxplot(data500$nNodes[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("N")
abline(h=0.435*500^(-0.446) * 40000)
boxplot(data20$nLinks[,c(4,3,2,1)]/data20$nNodes[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("p")
abline(h=1.531*20^(-0.523) * 40000^(-0.032))
boxplot(data100$nLinks[,c(4,3,2,1)]/data100$nNodes[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("p")
abline(h=1.531*100^(-0.523) * 40000^(-0.032))
boxplot(data500$nLinks[,c(4,3,2,1)]/data500$nNodes[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("p")
abline(h=1.531*500^(-0.523) * 40000^(-0.032))
dev.off()

# Fig. S9 ####
pdf(file="FigS9_draft.pdf",width=17/2.54, height=12/2.54)
par(mfrow=c(2,3))
boxplot(data20$mean_dist[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(50,350)); title("mean distance")
boxplot(data100$mean_dist[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(50,350)); title("mean distance")
boxplot(data500$mean_dist[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T),ylim=c(50,350)); title("mean distance")
boxplot(data20$CVm_0_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("CVm 0")
boxplot(data100$CVm_0_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("CVm 0")
boxplot(data500$CVm_0_H[,c(4,3,2,1)], notch=T, col=hcl.colors(4,"spectral",rev=T)); title("CVm 0")
dev.off()


# Fig. S10 ####
pdf("FigS10_draft.pdf",width=17/2.54,height=7/2.54)
df <- read.csv("results/TableS1_draft.csv")
Acell <- df$A_cells; Akm2 <- df$A_km2
par(mfrow=c(1,2))
hist(Acell[1:50],ylim=c(0,15))
hist(log10(Akm2[1:50]),breaks=seq(2,5,0.5),ylim=c(0,20))
dev.off()





