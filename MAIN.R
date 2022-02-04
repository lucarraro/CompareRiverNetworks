rm(list=ls())

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(elevatr)
library(sp)
library(raster)
library(rgdal)
library(shapefiles)
library(spam)
library(fields)
library(Matrix)
library(OCNet)

# follow instructions on https://hydrology.usu.edu/taudem/taudem5/TauDEMRScript.txt to install TauDEM

source("support_functions.R")

rivers <- read.csv("rivers/real_rivers.csv")
riverNames <- rivers$Name

# Generate rivers ####
if (!file.exists("results/TableS1_draft.csv")){
  Acell <- Akm2 <- cellsize <- zoom <- lat <- long <- numeric(length(rivers$Name))
  dA <- NULL
  for (i in 1:50){
    nam <- rivers$Name[i]
    cat(sprintf("%s \n",nam))
    outlet_pos = c(rivers$x_outlet[i], rivers$y_outlet[i])
    domain = c(rivers$x_ll[i], rivers$x_tr[i], rivers$y_ll[i], rivers$y_tr[i])
    EPSG=rivers$EPSG[i]; 
    zoom[i] <- rivers$zoom[i] 
    river <- river_from_DEM(outlet_pos=outlet_pos, domain=domain, EPSG=EPSG,zoom=zoom[i], showFig=F)
    Akm2[i] <- river$CM$A/1e6
    Acell[i] <- river$CM$A/river$cellsize^2
    cellsize[i] <- river$cellsize
    # convert outlet coordinates into WGS84
    cc <- spTransform(SpatialPoints(list(outlet_pos[1],outlet_pos[2]),proj4string = CRS(SRS_string = paste0("EPSG:",EPSG))),CRS("+init=epsg:4326"))
    lat[i] <- cc@coords[2]; long[i] <- cc@coords[1]
    eval(parse(text=paste0(nam,' <- river')))
    eval(parse(text=paste0('save(',nam,', file="rivers/',nam,'.rda")')))
    cat("  \n")
  }
  df <- data.frame(name=rivers$Name, A_km2=Akm2, A_cells=Acell, l_m=cellsize, zoom=zoom, country=rivers$country, lat=lat, long=long)
  write.csv(df,file="results/TableS1_draft.csv")
}

# Generate OCN ####
for (i in 1:50){
  fname <- paste0('OCN/OCN',i,'.rda')
  if (!file.exists(fname)){
    cat(sprintf('i: %d \n',i))
    set.seed(i)
    OCN <- NULL
    eval(parse(text=paste0('OCN',i,' <- OCN')))
    eval(parse(text=paste0('save(OCN',i,',file="OCN/OCN',i,'.rda")')))
    
    OCN <- create_OCN(200, 200, outletPos = 4*i,
                      nIter = 80*200*200,
                      coolingRate = 0.5,
                      displayUpdates=2, typeInitialState = c("T","V","I")[(i %% 3) + 1]) # alternate I,T,V as initial states
    
    
    eval(parse(text=paste0('OCN',i,' <- OCN')))
    eval(parse(text=paste0('save(OCN',i,',file="OCN/OCN',i,'.rda")')))
  }
}

# Evaluate CVm and lambdaM (+ create RBN & BBT) ####
thrA_vec <- c(500,100,20) 
for (thrA in thrA_vec){
  fname <- paste0('results/data_new',thrA,'.rda')
  if (!file.exists(fname)){
    
    mean_dist <- data.frame(matrix(0,50,4))
    CVm_1000 <- lambdaM_1000 <- data.frame(matrix(0,50,4))
    CVm_200 <- lambdaM_200 <-  data.frame(matrix(0,50,4))
    CVm_100 <- lambdaM_100 <-  data.frame(matrix(0,50,4))
    CVm_20 <- lambdaM_20 <-  data.frame(matrix(0,50,4))
    CVm_10 <- lambdaM_10 <-  data.frame(matrix(0,50,4))
    nNodes <- nLinks <- data.frame(matrix(0,50,4))
    
    CVm_1000_H <- lambdaM_1000_H <- data.frame(matrix(0,50,4))
    CVm_200_H <- lambdaM_200_H <- data.frame(matrix(0,50,4))
    CVm_100_H <- lambdaM_100_H <- data.frame(matrix(0,50,4))
    CVm_20_H <- lambdaM_20_H <- data.frame(matrix(0,50,4))
    CVm_10_H <- lambdaM_10_H <- CVm_0_H <- data.frame(matrix(0,50,4))
    
    names(CVm_1000) <- names(lambdaM_1000)  <- c("Rivers","OCN","RBN","BBT")
    names(CVm_200) <- names(lambdaM_200)  <- c("Rivers","OCN","RBN","BBT")
    names(CVm_100) <- names(lambdaM_100)  <- c("Rivers","OCN","RBN","BBT")
    names(CVm_20) <- names(lambdaM_20)  <- c("Rivers","OCN","RBN","BBT")
    names(CVm_10) <- names(lambdaM_10)  <- c("Rivers","OCN","RBN","BBT")
    names(mean_dist)  <-  c("Rivers","OCN","RBN","BBT")
    names(nNodes) <- names(nLinks) <- c("Rivers","OCN","RBN","BBT")
    
    names(CVm_1000_H) <- names(lambdaM_1000_H)  <- c("Rivers","OCN","RBN","BBT")
    names(CVm_200_H) <- names(lambdaM_200_H)  <- c("Rivers","OCN","RBN","BBT")
    names(CVm_100_H) <- names(lambdaM_100_H)  <- c("Rivers","OCN","RBN","BBT")
    names(CVm_20_H) <- names(lambdaM_20_H)  <- c("Rivers","OCN","RBN","BBT")
    names(CVm_10_H) <- names(lambdaM_10_H)  <- names(CVm_0_H) <- c("Rivers","OCN","RBN","BBT")
    
    for (i in 1:50){
      cat(sprintf("i: %d \n",i))
      for (j in 1:4){
        if (j==1){ # real river
          fnam <- riverNames[i]
          cat(sprintf('   extract river %s ... \n',fnam))
          fname <- paste0('rivers/',fnam,'.rda')
          load(fname)
          eval(parse(text=paste0('river <- ',fnam)))
          eval(parse(text=paste0('rm(',fnam,')')))
          river <- extract_river(river, thrA)
          river <- paths_OCN(river, includeUnconnectedPaths=T)
        }
        if (j==2){ # OCN
          cat('   OCN... \n')
          fnam <- paste0('OCN',i)
          fname <- paste0('OCN/',fnam,'.rda')
          load(fname)
          eval(parse(text=paste0('OCN <- ',fnam)))
          eval(parse(text=paste0('rm(',fnam,')')))
          OCN <- landscape_OCN(OCN, displayUpdates = 2)
          OCN <- aggregate_OCN(OCN, thrA=thrA) 
          RN_nodes <- OCN$RN$nNodes; AG_nodes <- OCN$AG$nNodes
          river <- paths_OCN(OCN, includeUnconnectedPaths=T)
          rm(OCN)
        } else if (j==3){ # RBN
          set.seed(10*i)
          cat('   create RBN... \n')
          # create RBN
          river <- create_RBN(RN_nodes, AG_nodes/RN_nodes)
          fnam <- paste0('RBN',i,'_',thrA)
          eval(parse(text=paste0(fnam, '<- river')))
          eval(parse(text=paste0('save(',fnam,', file="RBN/',fnam,'.rda" )')))
          eval(parse(text=paste0('rm(',fnam,')')))
          cat('   RBN paths... \n')
          river <- paths_OCN(river, includeUnconnectedPaths=T)
          river[["cellsize"]] <- 1
        } else if (j==4){ # BBT
          set.seed(10*i+1)
          cat('   create BBT... \n')
          # create BBT
          river <- create_BBT(RN_nodes, AG_nodes/RN_nodes)
          fnam <- paste0('BBT',i,'_',thrA)
          eval(parse(text=paste0(fnam, '<- river')))
          eval(parse(text=paste0('save(',fnam,', file="BBT/',fnam,'.rda" )')))
          eval(parse(text=paste0('rm(',fnam,')')))
          cat('   BBT paths... \n')
          river <- paths_OCN(river, includeUnconnectedPaths=T)
          river[["cellsize"]] <- 1
        }
        gc(verbose=F) # clean some RAM
        cat('   calculate RN distances... \n')
        distances <- river$RN$downstreamLengthUnconnected + t(river$RN$downstreamLengthUnconnected) + 
          as.matrix.spam(river$RN$downstreamPathLength + t(river$RN$downstreamPathLength))
        
        patchSize <- river$RN$nUpstream^0.5/sum(river$RN$nUpstream^0.5) 
        HS_mat <- matrix(patchSize,river$RN$nNodes,1) %*% matrix(patchSize,1,river$RN$nNodes)
        mean_dist[i,j] <- mean(distances)/river$cellsize
        CVm_0_H[i,j] <- (sum(patchSize^2))^0.5/sum(patchSize)
        
        cat('   alpha=1000... \n')
        expdist <- exp(-distances/river$cellsize/1000)
        diag(expdist) <- numeric(river$RN$nNodes)
        expdist_HS <- expdist*HS_mat
        diag(expdist_HS) <- numeric(river$RN$nNodes)
        CVm_1000[i,j] <- ((1+1/river$RN$nNodes*sum(expdist))/river$RN$nNodes)^0.5
        EE <- eigen(expdist, only.values = T)
        lambdaM_1000[i,j] <- max(EE$values)
        CVm_1000_H[i,j] <- (sum(patchSize^2) + sum(expdist_HS))^0.5/sum(patchSize)
        EEH <- eigen(expdist_HS,only.values=T)
        lambdaM_1000_H[i,j] <- max(EEH$values)
        rm(expdist,expdist_HS)
        gc(verbose=F)
        
        cat('   alpha=200... \n')
        expdist <- exp(-distances/river$cellsize/200)
        diag(expdist) <- numeric(river$RN$nNodes)
        expdist_HS <- expdist*HS_mat
        diag(expdist_HS) <- numeric(river$RN$nNodes)
        CVm_200[i,j] <- ((1+1/river$RN$nNodes*sum(expdist))/river$RN$nNodes)^0.5
        EE <- eigen(expdist, only.values = T)
        lambdaM_200[i,j] <- max(EE$values)
        CVm_200_H[i,j] <- (sum((patchSize)^2) + sum(expdist_HS))^0.5/sum(patchSize)
        EEH <- eigen(expdist_HS,only.values=T)
        lambdaM_200_H[i,j] <- max(EEH$values)
        rm(expdist,expdist_HS)
        gc(verbose=F)
        
        cat('   alpha=100... \n')
        expdist <- exp(-distances/river$cellsize/100)
        diag(expdist) <- numeric(river$RN$nNodes)
        expdist_HS <- expdist*HS_mat
        diag(expdist_HS) <- numeric(river$RN$nNodes)
        CVm_100[i,j] <- ((1+1/river$RN$nNodes*sum(expdist))/river$RN$nNodes)^0.5
        EE <- eigen(expdist, only.values = T)
        lambdaM_100[i,j] <- max(EE$values)
        CVm_100_H[i,j] <- (sum((patchSize)^2) + sum(expdist_HS))^0.5/sum(patchSize)
        EEH <- eigen(expdist_HS,only.values=T)
        lambdaM_100_H[i,j] <- max(EEH$values)
        rm(expdist,expdist_HS)
        gc(verbose=F)
        
        cat('   alpha=20... \n')
        expdist <- exp(-distances/river$cellsize/20)
        diag(expdist) <- numeric(river$RN$nNodes)
        expdist_HS <- expdist*HS_mat
        diag(expdist_HS) <- numeric(river$RN$nNodes)
        CVm_20[i,j] <- ((1+1/river$RN$nNodes*sum(expdist))/river$RN$nNodes)^0.5
        EE <- eigen(expdist, only.values = T)
        lambdaM_20[i,j] <- max(EE$values)
        CVm_20_H[i,j] <- (sum((patchSize)^2) + sum(expdist_HS))^0.5/sum(patchSize)
        EEH <- eigen(expdist_HS,only.values=T)
        lambdaM_20_H[i,j] <- max(EEH$values)
        rm(expdist,expdist_HS)
        gc(verbose=F)
        
        cat('   alpha=10... \n')
        expdist <- exp(-distances/river$cellsize/10)
        diag(expdist) <- numeric(river$RN$nNodes)
        expdist_HS <- expdist*HS_mat
        diag(expdist_HS) <- numeric(river$RN$nNodes)
        CVm_10[i,j] <- ((1+1/river$RN$nNodes*sum(expdist))/river$RN$nNodes)^0.5
        EE <- eigen(expdist, only.values = T)
        lambdaM_10[i,j] <- max(EE$values)
        CVm_10_H[i,j] <- (sum((patchSize)^2) + sum(expdist_HS))^0.5/sum(patchSize)
        EEH <- eigen(expdist_HS,only.values=T)
        lambdaM_10_H[i,j] <- max(EEH$values)
        rm(expdist,expdist_HS)
        gc(verbose=F)
        
        rm(distances,HS_mat)
        gc(verbose=F)
        
        nNodes[i,j] <- river$RN$nNodes
        nLinks[i,j] <- river$AG$nNodes
        rm(river)
        gc(verbose = F)
        
        data <- list(CVm_1000=CVm_1000, lambdaM_1000=lambdaM_1000,
                     CVm_200=CVm_200, lambdaM_200=lambdaM_200,
                     CVm_100=CVm_100, lambdaM_100=lambdaM_100,
                     CVm_20=CVm_20, lambdaM_20=lambdaM_20,
                     CVm_10=CVm_10, lambdaM_10=lambdaM_10,
                     CVm_1000_H=CVm_1000_H, lambdaM_1000_H=lambdaM_1000_H,
                     CVm_200_H=CVm_200_H, lambdaM_200_H=lambdaM_200_H,
                     CVm_100_H=CVm_100_H, lambdaM_100_H=lambdaM_100_H,
                     CVm_20_H=CVm_20_H, lambdaM_20_H=lambdaM_20_H,
                     CVm_10_H=CVm_10_H, lambdaM_10_H=lambdaM_10_H,
                     CVm_0_H=CVm_0_H, mean_dist=mean_dist,  
                     nNodes=nNodes,  nLinks=nLinks)
        eval(parse(text=paste0('data',thrA, '<- data')))
        eval(parse(text=paste0('save(data',thrA, ', file="results/data_new',thrA,'.rda")')))
      }
    }
  }
}
# Evaluate Horton Laws ####
if (!file.exists("results/HortonLaws.rda")){
  N_omega_thrA20 <- L_omega_thrA20  <- sO_thrA20 <- vector("list",4)
  names(N_omega_thrA20) <- names(L_omega_thrA20)  <- names(sO_thrA20) <- c("Rivers","OCN","RBN","BBT")
  N_omega_thrA100 <- L_omega_thrA100  <- sO_thrA100 <- vector("list",4)
  names(N_omega_thrA100) <- names(L_omega_thrA100)  <- names(sO_thrA100) <- c("Rivers","OCN","RBN","BBT")
  N_omega_thrA500 <- L_omega_thrA500  <- sO_thrA500 <- vector("list",4)
  names(N_omega_thrA500) <- names(L_omega_thrA500)  <- names(sO_thrA500) <- c("Rivers","OCN","RBN","BBT")
  
  RB_thrA20 <- RL_thrA20  <- data.frame(matrix(0,50,4))
  names(RB_thrA20) <- names(RL_thrA20)  <-  c("Rivers","OCN","RBN","BBT")
  RB_thrA100 <- RL_thrA100  <- data.frame(matrix(0,50,4))
  names(RB_thrA100) <- names(RL_thrA100)  <-  c("Rivers","OCN","RBN","BBT")
  RB_thrA500 <- RL_thrA500  <- data.frame(matrix(0,50,4))
  names(RB_thrA500) <- names(RL_thrA500) <-  c("Rivers","OCN","RBN","BBT")
  
  for (thrA in c(500,100,20)){
    for (ind in (1:50)){
      cat(sprintf('thrA: %d  -  i: %d \n',thrA,ind))
      for (type in c("Rivers","OCN","RBN","BBT")){
        cat(sprintf('   %s \n',type))
        if (type == "Rivers"){
          load(paste0('Rivers/',riverNames[ind],'.rda'))
          eval(parse(text=paste0('river <- ',riverNames[ind])))
          eval(parse(text=paste0('rm(',riverNames[ind],')')))
          river <- extract_river(river, thrA)
        } else if (type == "OCN"){
          load(paste0('OCN/OCN',ind,'.rda'))
          eval(parse(text=paste0('OCN <- OCN',ind)))  
          OCN <- landscape_OCN(OCN)
          OCN <- aggregate_OCN(OCN,thrA=thrA)
          river <- OCN
          rm(OCN)
          eval(parse(text=paste0('rm(OCN',ind,')')))
        } else if (type == "RBN"){
          load(paste0('RBN/RBN',ind,'_',thrA,'.rda'))
          eval(parse(text=paste0('river <- RBN',ind,'_',thrA)))
          eval(parse(text=paste0('rm(RBN',ind,'_',thrA,')')))
          river[["cellsize"]] <- 1
        } else if (type == "BBT"){
          load(paste0('BBT/BBT',ind,'_',thrA,'.rda'))
          eval(parse(text=paste0('river <- BBT',ind,'_',thrA)))
          eval(parse(text=paste0('rm(BBT',ind,'_',thrA,')')))
          river[["cellsize"]] <- 1
        }
        
        # evaluate streams
        AGtoStream <- numeric(river$AG$nNodes)
        
        stream <- 1:(sum(river$AG$streamOrder==1))
        streamOrder <- 1+numeric(length(stream))
        AGtoStream[which(river$AG$streamOrder==1)] <- 1:(sum(river$AG$streamOrder==1))
        
        for (ss in 2:max(river$AG$streamOrder)){
          tmp <- which(river$AG$streamOrder==ss)
          tmp2 <- tmp
          for (i in 1:length(tmp)){
            link <- tmp[i]
            if (link %in% tmp2 & max(river$AG$streamOrder[which(river$AG$downNode==link)])<ss ){ # if not processed yet and is the most upstream link, it's a new stream
              idStream <- length(stream)+1
              stream <- c(stream,idStream)
              streamOrder <- c(streamOrder,ss)
              while (isTRUE(river$AG$streamOrder[link] == ss)){
                AGtoStream[link] <- idStream
                tmp2 <- tmp2[-(which(tmp2==link))]
                link <- river$AG$downNode[link]
              }
            }
          }
        }
        
        # eval lengths
        leng <- numeric(length(stream))
        for (i in 1:length(stream)){
          leng[i] <- sum(river$AG$leng[which(AGtoStream==i)])/river$cellsize
        }
        
        N_omega <- L_omega <-  numeric(max(river$AG$streamOrder))
        sO <- 1:length(N_omega)
        for (ss in sO){
          N_omega[ss] <- sum(streamOrder==ss)
          L_omega[ss] <- mean(leng[streamOrder==ss])
        }
        L_omega[L_omega==0] <- NA
        lmod_RB <- lm(log(N_omega) ~ sO)
        lmod_RL <- lm(log(L_omega[1:(length(L_omega)-1)]) ~ sO[1:(length(L_omega)-1)]); 
        
        if (thrA == 20){
          N_omega_thrA20[[type]] <- c(N_omega_thrA20[[type]], N_omega)
          L_omega_thrA20[[type]] <- c(L_omega_thrA20[[type]], L_omega)
          sO_thrA20[[type]] <- c(sO_thrA20[[type]], sO)
          RB_thrA20[[type]][ind] <- exp(-lmod_RB$coefficients[2])
          RL_thrA20[[type]][ind] <- exp(lmod_RL$coefficients[2])
        } else if (thrA==100){
          N_omega_thrA100[[type]] <- c(N_omega_thrA100[[type]], N_omega)
          L_omega_thrA100[[type]] <- c(L_omega_thrA100[[type]], L_omega)
          sO_thrA100[[type]] <- c(sO_thrA100[[type]], sO)
          RB_thrA100[[type]][ind] <- exp(-lmod_RB$coefficients[2])
          RL_thrA100[[type]][ind] <- exp(lmod_RL$coefficients[2])
        } else if (thrA==500){
          N_omega_thrA500[[type]] <- c(N_omega_thrA500[[type]], N_omega)
          L_omega_thrA500[[type]] <- c(L_omega_thrA500[[type]], L_omega)
          sO_thrA500[[type]] <- c(sO_thrA500[[type]], sO)
          RB_thrA500[[type]][ind] <- exp(-lmod_RB$coefficients[2])
          RL_thrA500[[type]][ind] <- exp(lmod_RL$coefficients[2])
        }
      }
    }
  }
  ll <- list(sO_thrA500=sO_thrA500, N_omega_thrA500=N_omega_thrA500, 
             L_omega_thrA500=L_omega_thrA500, 
             sO_thrA100=sO_thrA100, N_omega_thrA100=N_omega_thrA100, 
             L_omega_thrA100=L_omega_thrA100, 
             sO_thrA20=sO_thrA20, N_omega_thrA20=N_omega_thrA20, 
             L_omega_thrA20=L_omega_thrA20, 
             RB_thrA20=RB_thrA20, RL_thrA20=RL_thrA20, 
             RB_thrA100=RB_thrA100, RL_thrA100=RL_thrA100, 
             RB_thrA500=RB_thrA500, RL_thrA500=RL_thrA500)
  
  save(ll,file="results/HortonLaws.rda")   
}
# Evaluate scaling laws ####
if (!file.exists("results/scaling.rda")){
  
  # generate space-filling random networks
  # this takes a lot of time (> 1 day) and writes ~ 3GB of data 
  thrA <- 1
  for (i in 1:50){
    cat(sprintf("i: %d \n",i))
    cat('   OCN... \n')
    fnam <- paste0('OCN',i)
    fname <- paste0('OCN/',fnam,'.rda')
    load(fname)
    eval(parse(text=paste0('OCN <- ',fnam)))
    eval(parse(text=paste0('rm(',fnam,')')))
    
    OCN <- landscape_OCN(OCN, displayUpdates = 2)
    OCN <- aggregate_OCN(OCN, thrA=thrA) 
    RN_nodes <- OCN$RN$nNodes; AG_nodes <- OCN$AG$nNodes
    rm(OCN)
    gc(verbose=F)
    
    set.seed(10*i)
    cat('   create RBN... \n')
    river <- create_RBN(RN_nodes, AG_nodes/RN_nodes)
    fnam <- paste0('RBN',i,'_',thrA)
    eval(parse(text=paste0(fnam, '<- river')))
    eval(parse(text=paste0('save(',fnam,', file="RBN/',fnam,'.rda" )')))
    eval(parse(text=paste0('rm(',fnam,')')))
    rm(river)
    gc(verbose=F)
    
    set.seed(10*i+1)
    cat('   create BBT... \n')
    river <- create_BBT(RN_nodes, AG_nodes/RN_nodes)
    fnam <- paste0('BBT',i,'_',thrA)
    eval(parse(text=paste0(fnam, '<- river')))
    eval(parse(text=paste0('save(',fnam,', file="BBT/',fnam,'.rda" )')))
    eval(parse(text=paste0('rm(',fnam,')')))
    rm(river)
    gc(verbose=F)
  }
  
  
  dA_river <- NULL
  for (i in 1:length(riverNames)){
    fnam <- riverNames[i]
    eval(parse(text=paste0('load("rivers/',fnam,'.rda")')))
    eval(parse(text=paste0('river <- ',fnam)))
    eval(parse(text=paste0('rm(',fnam,')')))
    dA_river <- c(dA_river, river$FD$A[!is.na(river$FD$A)])
  }
  P_A_river <- numeric(max(dA_river))
  for (i in 1:length(P_A_river)){
    P_A_river[i] <- sum(dA_river >= i)/length(dA_river)
  }
  
  
  dA_OCN <- NULL
  for (i in 1:length(riverNames)){
    fnam <- paste0('OCN',i)
    eval(parse(text=paste0('load("OCN/',fnam,'.rda")')))
    eval(parse(text=paste0('river <- ',fnam)))
    eval(parse(text=paste0('rm(',fnam,')')))
    dA_OCN <- c(dA_OCN, river$FD$A[!is.na(river$FD$A)])
  }
  P_A_OCN <- numeric(max(dA_OCN))
  for (i in 1:length(P_A_OCN)){
    P_A_OCN[i] <- sum(dA_OCN >= i)/length(dA_OCN)
  }
  
  
  
  dA_RBN <- NULL
  for (i in 1:50){
    fname <- paste0('RBN/RBN',i,'_',thrA,'.rda')
    load(fname)
    eval(parse(text=paste0('RBN <- RBN',i,'_',thrA)))
    eval(parse(text=paste0('rm(RBN',i,'_',thrA,')')))
    dA_RBN <- c(dA_RBN, RBN$RN$nUpstream)
  }
  
  P_A_RBN <- numeric(max(dA_RBN))
  for (i in 1:length(P_A_RBN)){
    P_A_RBN[i] <- sum(dA_RBN >= i)/length(dA_RBN)
  }
  
  dA_BBT <- NULL
  for (i in 1:50){  
    cat(sprintf('BBT  -  i: %d   -  thrA: %d \n',i,thrA))
    fname <- paste0('BBT/BBT',i,'_',thrA,'.rda')
    load(fname)
    eval(parse(text=paste0('BBT <- BBT',i,'_',thrA)))
    eval(parse(text=paste0('rm(BBT',i,'_',thrA,')')))
    dA_BBT <- c(dA_BBT, BBT$RN$nUpstream)
  }
  
  P_A_BBT <- numeric(max(dA_BBT))
  for (i in 1:length(P_A_BBT)){
    P_A_BBT[i] <- sum(dA_BBT >= i)/length(dA_BBT)
  }
  
  
  scaling <- list(P_A_river=P_A_river, P_A_OCN=P_A_OCN, P_A_RBN=P_A_RBN, P_A_BBT=P_A_BBT, 
                  dA_river=dA_river, dA_OCN=dA_OCN, dA_RBN=dA_RBN, dA_BBT=dA_BBT)
  save(scaling,file="results/scaling.rda")
  
} 



