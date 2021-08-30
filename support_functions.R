
# download real river from DEM
river_from_DEM <- function(outlet_pos,domain,EPSG,zoom,
                           showFig=F){
  
  x_outlet <- outlet_pos[1]; y_outlet <- outlet_pos[2]
  x_ll <- domain[1]; x_tr <- domain[2]; y_ll <- domain[3]; y_tr <- domain[4]
  
  crs_str <- CRS(SRS_string = paste0("EPSG:",EPSG)) # 21781 is EPSG code for CRS CH1903/LV03. See also http://www.epsg-registry.org 
  # 32632 is UTM 32N
  
  loc.df <- data.frame(x=c(x_ll, x_tr), y=c(y_ll,y_tr))
  
  cat("Downloading DEM... ")
  z <- get_elev_raster(locations = loc.df, prj = crs_str, z=zoom, verbose=F) # change to z=12 for small Thur
  
  cellsizeX <- (z@extent@xmax - z@extent@xmin)/z@ncols
  cellsizeY <- (z@extent@ymax - z@extent@ymin)/z@nrows
  cellsize <- sqrt(cellsizeX*cellsizeY)
  
  z <- crop(z, domain)
  z <- reclassify(z, cbind(NA,0)) # all pixels with elev=NA are set to 0. Then the pit remove algortihm will take care of them
  
  writeRaster(z,filename="DEM.tif",format="GTiff")
  
  cat("Pit remove... \n")
  system("mpiexec -n 8 pitremove -z DEM.tif -fel DEMfel.tif",show.output.on.console=F,invisible=F)
  fel=raster("DEMfel.tif")
  
  # D8 flow directions
  cat("D8 flow directions... \n")
  system("mpiexec -n 8 D8Flowdir -p DEMp.tif -sd8 DEMsd8.tif -fel DEMfel.tif",show.output.on.console=F,invisible=F)  
  p=raster("DEMp.tif")
  sd8=raster("DEMsd8.tif")
  
  # Contributing area
  cat("Find contributing area... \n")
  system("mpiexec -n 8 AreaD8 -p DEMp.tif -ad8 DEMad8.tif",show.output.on.console=F,invisible=F)
  ad8=raster("DEMad8.tif")
  
  # Threshold
  cat("Apply area threshold... \n")
  #tA <- thrA/cellsizeX/cellsizeY
  tA <- 1000 # increase threshold area to be sure to get the right outlet position
  system(sprintf("mpiexec -n 8 Threshold -ssa DEMad8.tif -src DEMsrc.tif -thresh %.2f",tA),show.output.on.console=F,invisible=F)
  src=raster("DEMsrc.tif")
  
  shp.point(x_outlet,y_outlet,"ApproxOutlet")
  
  # Move Outlets
  cat("Move outlet to stream... \n")
  system("mpiexec -n 8 moveoutletstostreams -p DEMp.tif -src DEMsrc.tif -o ApproxOutlet.shp -om Outlet.shp",show.output.on.console=F,invisible=F)
  
  # Contributing area upstream of outlet
  cat("Contributing area upstream of outlet... \n")
  system("mpiexec -n 8 Aread8 -p DEMp.tif -o Outlet.shp -ad8 DEMssa.tif",show.output.on.console=F,invisible=F)
  ssa=raster("DEMssa.tif")
  
  cont <- rasterToContour(reclassify(ssa,cbind(NA,0)),levels=0.5)
  ext <- extent(cont)
  ext@xmin <- ext@xmin - cellsize; ext@xmax <- ext@xmax + cellsize 
  ext@ymin <- ext@ymin - cellsize; ext@ymax <- ext@ymax + cellsize 
  
  ssa <- crop(ssa,ext)
  p <- crop(p,ext)
  fel <- crop(fel,ext)
  
  if (showFig){
    old.par <- par()
    par(mfrow=c(1,2))
    
    brks=seq(0,max(values(z)),100)
    plot(z,col=terrain.colors(length(brks)),breaks=brks,axes=FALSE)
    plot(cont,add=T)
    
    plot(ad8,col=bpy.colors(255),axes=FALSE)
    # ,xlim=c(x_outlet-50*cellsize,x_outlet+50*cellsize),ylim=c(y_outlet-50*cellsize,y_outlet+50*cellsize)
    plot(cont,col="green",add=T)
    par(old.par)
  }
  
  file.remove("DEM.tif","DEMfel.tif","DEMp.tif","DEMsd8.tif","DEMad8.tif","DEMsrc.tif","DEMssa.tif")
  file.remove("Outlet.shp","Outlet.shx","Outlet.dbf","Outlet.prj")
  file.remove("ApproxOutlet.dbf","ApproxOutlet.shp","ApproxOutlet.shx")
  
  tmp <- coordinates(ssa)
  X_FD <- tmp[,1]
  Y_FD <- tmp[,2]
  rm(tmp)
  Z_FD <- values(fel)
  A_FD <- values(ssa)
  
  ncols <- fel@ncols
  nrows <- fel@nrows
  
  tmp <- values(p)
  Length_FD <- cellsize*(1 + as.numeric(tmp %% 2 ==0)*(sqrt(2)-1)) 
  
  
  FD <- list(X=X_FD, Y=Y_FD, Z=Z_FD, A=A_FD, leng=Length_FD, flowDir=values(p))
  XContour <-  cont@lines[[1]]@Lines[[1]]@coords[,1]
  YContour <-  cont@lines[[1]]@Lines[[1]]@coords[,2]
  CM <- list(XContour=XContour, YContour=YContour, A=max(A_FD,na.rm=T)*cellsize^2)
  river <- list(FD=FD, CM=CM, nrows=nrows, ncols=ncols, cellsize=cellsize)
  
  invisible(river)
}

# extract river with given thrA value
extract_river <- function(river, thrA, maxReachLength=Inf){
  
  cat("Connectivity at RN level... \n")
  indexRNNodes <- which(river$FD$A >=thrA)
  
  X_RN <- river$FD$X[indexRNNodes]
  Y_RN <- river$FD$Y[indexRNNodes]
  A_RN <- river$FD$A[indexRNNodes]
  Z_RN <- river$FD$Z[indexRNNodes]
  Length_RN <- river$FD$leng[indexRNNodes]
  
  ## W, downNode at RN level
  
  nNodes_RN <- length(indexRNNodes)
  W_RN <- spam(0,nNodes_RN,nNodes_RN)
  ind <- matrix(0,nNodes_RN,2)
  downNode_RN <- numeric(nNodes_RN)
  Slope_RN <-  numeric(nNodes_RN)
  
  k <- 1
  for (i in 1:nNodes_RN){
    mov <- neigh(river$FD$flowDir[indexRNNodes[i]])
    d <- which(indexRNNodes==(indexRNNodes[i]+mov[1]+mov[2]*river$ncols)) # indices from top-left corner to the right, then next row...
    if (length(d)!=0){
      ind[k, ] <- c(i,d)
      k <- k + 1
      Slope_RN[i] <- (Z_RN[i]-Z_RN[d])/Length_RN[i]
    } 
  }
  ind <- ind[-k, ]
  downNode_RN[ind[,1]] <- ind[,2]
  W_RN[ind] <- 1
  rm(ind)
  Outlet_RN <- which(downNode_RN==0)
  
  # find AG nodes ####
  cat("Connectivity at AG level... \n")
  DegreeIn <- colSums(W_RN)
  DegreeOut <- rowSums(W_RN)
  Confluence <- DegreeIn>1
  Source <- DegreeIn==0
  SourceOrConfluence <- Source|Confluence
  ConfluenceNotOutlet <- Confluence&(downNode_RN!=0)
  ChannelHeads <- SourceOrConfluence  #Source|ConfluenceNotOutlet
  
  OutletNotChannelHead <- (downNode_RN==0)&(!ChannelHeads)
  IsNodeAG <- SourceOrConfluence|OutletNotChannelHead
  whichNodeAG <- which(IsNodeAG)
  
  nNodes_AG <- sum(IsNodeAG)
  Length_AG <- numeric(nNodes_AG)
  RN_to_AG <- numeric(nNodes_RN)
  reachID <- 1
  X_AG <- NaN*numeric(nNodes_AG)
  Y_AG <- NaN*numeric(nNodes_AG)
  Z_AG <- NaN*numeric(nNodes_AG)
  A_AG <- NaN*numeric(nNodes_AG)
  while (length(whichNodeAG) != 0){ # explore all AG Nodes
    i <- whichNodeAG[1] # select the first
    RN_to_AG[i] <- reachID 
    j <- downNode_RN[i] 
    X_AG[reachID] <- X_RN[i]
    Y_AG[reachID] <- Y_RN[i]
    Z_AG[reachID] <- Z_RN[i]
    A_AG[reachID] <- A_RN[i]
    Length_AG[reachID] <- Length_RN[i]
    tmp_length <- Length_RN[i]
    tmp <- NULL
    j0 <- j
    while (!IsNodeAG[j] && j!=0) {
      tmp <- c(tmp, j)
      tmp_length <-  tmp_length + Length_RN[j]
      j_old <- j
      j <- downNode_RN[j]} 
    
    if (tmp_length > maxReachLength){
      n_splits <- ceiling(tmp_length/maxReachLength)
      new_maxLength <- tmp_length/n_splits
      j <- j0
      while (!IsNodeAG[j] && j!=0 && Length_AG[reachID] <= new_maxLength) {
        RN_to_AG[j] <- reachID 
        Length_AG[reachID] <-  Length_AG[reachID] + Length_RN[j]
        j_old <- j
        j <- downNode_RN[j]}
      if (Length_AG[reachID] > new_maxLength){
        j <- j_old
        Length_AG[reachID] <-  Length_AG[reachID] - Length_RN[j]
        ChannelHeads[j] <- 1
        whichNodeAG <- c(whichNodeAG,j)}
      
    } else {
      RN_to_AG[tmp] <- reachID
      Length_AG[reachID] <- tmp_length
    }
    
    reachID <- reachID + 1
    whichNodeAG <- whichNodeAG[-1]
  }
  nNodes_AG <- length(X_AG)
  
  # W, downNode at AG level ####
  
  downNode_AG <- numeric(nNodes_AG)
  W_AG <- spam(0,nNodes_AG,nNodes_AG)
  ind <- matrix(0,nNodes_AG,2)
  reachID <- sum(ChannelHeads) + 1
  for (i in 1:nNodes_RN){ 
    if (downNode_RN[i] != 0 && RN_to_AG[downNode_RN[i]] != RN_to_AG[i]) {
      downNode_AG[RN_to_AG[i]] <- RN_to_AG[downNode_RN[i]]
      ind[RN_to_AG[i],] <- c(RN_to_AG[i],downNode_AG[RN_to_AG[i]])
    }
  }
  ind <- ind[-which(ind[,1]==0),]
  W_AG[ind] <- 1
  Outlet_AG <- RN_to_AG[Outlet_RN]
  
  AG_to_RN <- vector("list", nNodes_AG)
  for(i in 1:nNodes_AG) { # attribute river network pixels to fields of the AG_to_FD list 
    AG_to_RN[[i]] <- which(RN_to_AG==i) 
  }
  
  FD_to_SC <- NA*numeric(length(river$FD$flowDir))
  SC_to_FD <- vector("list",nNodes_AG)
  FD_to_SC[indexRNNodes] <- RN_to_AG
  for (i in 1:nNodes_AG){
    SC_to_FD[[i]] <- indexRNNodes[which(RN_to_AG==i)]
  }
  
  
  # find FD_to_SC
  drainageArea <- river$FD$A
  indexFDNodes <- which(drainageArea>0)
  ind_head <- indexFDNodes[which(drainageArea[indexFDNodes]==1)]
  
  for (i in 1:length(ind_head)){
    d <- ind_head[i]
    k <- NA; d_new <- d; sub_d <- numeric(0)
    while (is.na(k)){
      k <- FD_to_SC[d_new]
      if (is.na(k)){
        sub_d <- c(sub_d, d_new)
        mov <- neigh(river$FD$flowDir[d_new])
        d_new <- d_new + mov[1] + mov[2]*river$ncols # indices from top-left corner to the right, then next row...
      }
    }
    FD_to_SC[sub_d] <- k
    SC_to_FD[[k]] <- c(SC_to_FD[[k]], sub_d)
  }
  
  
  # Upstream_RN : list containing IDs of all reaches upstream of each reach (plus reach itself)
  Upstream_RN <- vector("list",nNodes_RN)
  nUpstream_RN <- numeric(nNodes_RN)
  for (i in 1:nNodes_RN){
    UpOneLevel <- which(downNode_RN==i) # find reaches at one level upstream
    Upstream_RN[[i]] <- UpOneLevel      # add them to the list
    while (length(UpOneLevel)!=0) { # continue until there are no more reaches upstream
      ContinuePath <- UpOneLevel # jump 1 level above
      UpOneLevel <- which(downNode_RN %in% ContinuePath) # find reaches at one level upstream
      Upstream_RN[[i]] <- c(Upstream_RN[[i]],UpOneLevel) # add them to the list
    }
    Upstream_RN[[i]] <- c(Upstream_RN[[i]],i)
    nUpstream_RN[i] <- length(Upstream_RN[[i]])
    if ((i %% 100)==0){
      message(sprintf("%.2f%% done\r",100*i/nNodes_RN),appendLF = F)
    }
  }
  message('100.00% done \n',appendLF = F)
  
  # Upstream_AG : list containing IDs of all reaches upstream of each reach (plus reach itself)
  Upstream_AG <- vector("list",nNodes_AG)
  nUpstream_AG <- numeric(nNodes_AG)
  for (i in 1:nNodes_AG){
    UpOneLevel <- which(downNode_AG==i) # find reaches at one level upstream
    Upstream_AG[[i]] <- UpOneLevel      # add them to the list
    while (length(UpOneLevel)!=0) { # continue until there are no more reaches upstream
      ContinuePath <- UpOneLevel # jump 1 level above
      UpOneLevel <- which(downNode_AG %in% ContinuePath) # find reaches at one level upstream
      Upstream_AG[[i]] <- c(Upstream_AG[[i]],UpOneLevel) # add them to the list
    }
    Upstream_AG[[i]] <- c(Upstream_AG[[i]],i)
    nUpstream_AG[i] <- length(Upstream_AG[[i]])
  }
  
  
  # calculate Strahler stream order
  StreamOrder_AG <- numeric(nNodes_AG)
  for (i in 1:nNodes_AG){
    j <- order(nUpstream_AG)[i] # index that explores reaches in a downstream direction
    tmp <- which(downNode_AG==j) # set of reaches draining into j
    if (length(tmp)>0){
      IncreaseOrder <- sum(StreamOrder_AG[tmp]==max(StreamOrder_AG[tmp])) # check whether tmp reaches have the same stream order
      if (IncreaseOrder > 1) {
        StreamOrder_AG[j] <- 1 + max(StreamOrder_AG[tmp]) # if so, increase stream order
      } else {StreamOrder_AG[j] <- max(StreamOrder_AG[tmp])} # otherwise, keep previous stream order
    } else {StreamOrder_AG[j] <- 1} # if j is an headwater, impose StreamOrder = 1
  }
  
  
  #print('Length and slope at AG level...',quote=FALSE) 
  # Calculate length and slopes of reaches
  #Length_AG <- rep(0,Nnodes_AG)
  Slope_AG <- numeric(nNodes_AG)
  for (i in 1:nNodes_AG){
    #Length_AG[i] <- sum(OCN$FD$leng[AG_to_FD[[i]]])
    Slope_AG[i] <- (Slope_RN[RN_to_AG==i] %*% Length_RN[RN_to_AG==i])/Length_AG[i] # scalar product between vector of slopes and lengths of nodes at RN level belonging to reach i 
  }
  
  
  nNodes_SC <- nNodes_AG
  Z_SC <- numeric(nNodes_SC)
  Alocal_SC <- numeric(nNodes_SC)
  for (i in 1:nNodes_SC) {
    Z_SC[i] <- mean(river$FD$Z[which(FD_to_SC==i)])
    Alocal_SC[i] <- sum(FD_to_SC==i,na.rm=T)
  }
  
  Areach_AG <- numeric(nNodes_AG)
  for (i in 1:nNodes_AG) {
    Areach_AG[i] <- sum(Alocal_SC[Upstream_AG[[i]]])  
  }
  
  # coordinates of AG nodes considered at the downstream end of the respective edge
  XReach <- numeric(nNodes_AG)
  YReach <- numeric(nNodes_AG)
  ZReach <- numeric(nNodes_AG)
  for (i in 1:nNodes_AG){
    tmp <- AG_to_RN[[i]]
    ind <- which(A_RN[tmp]==max(A_RN[tmp]))
    node <- tmp[ind]
    XReach[i] <- X_RN[node]
    YReach[i] <- Y_RN[node]
    ZReach[i] <- Z_RN[node]
  }
  XReach[Outlet_AG] <- NaN
  YReach[Outlet_AG] <- NaN
  ZReach[Outlet_AG] <- NaN
  
  river$FD[["toSC"]] <- FD_to_SC
  
  river$RN[["A"]] <- A_RN*river$cellsize^2
  river$RN[["downNode"]] <- downNode_RN
  river$RN[["outlet"]] <- which(downNode_RN==0)
  river$RN[["leng"]] <- Length_RN
  river$RN[["upstream"]] <- Upstream_RN
  river$RN[["nUpstream"]] <- nUpstream_RN
  river$RN[["X"]] <- X_RN
  river$RN[["Y"]] <- Y_RN
  river$RN[["toAGReach"]] <- RN_to_AG
  river$RN[["nNodes"]] <- nNodes_RN
  
  river$AG[["A"]] <- A_AG*river$cellsize^2
  river$AG[["AReach"]] <- Areach_AG*river$cellsize^2
  river$AG[["downNode"]] <- downNode_AG
  river$AG[["upstream"]] <- Upstream_AG
  river$AG[["nUpstream"]] <- nUpstream_AG
  river$AG[["leng"]] <- Length_AG
  river$AG[["outlet"]] <- Outlet_AG
  river$AG[["slope"]] <- Slope_AG
  river$AG[["streamOrder"]] <- StreamOrder_AG
  river$AG[["Upstream"]] <- Upstream_AG
  river$AG[["X"]] <- X_AG
  river$AG[["Y"]] <- Y_AG
  river$AG[["Z"]] <- Z_AG
  river$AG[["W"]] <- W_AG
  river$AG[["nNodes"]] <- nNodes_AG
  pl <- initial_permutation(downNode_AG)
  river$AG[["perm"]] <- pl$perm
  
  river$SC[["toFD"]] <- SC_to_FD
  river$SC[["A"]] <- Alocal_SC*river$cellsize^2
  
  river[["thrA"]] <- thrA
  
  cat(sprintf("Total catchment area: %d cells  -   %.2f km2   -  nNodes: %d   -  nLinks: %d \n",
              max(A_RN),river$CM$A/1e6,river$RN$nNodes,river$AG$nNodes))
  
  invisible(river)
}


# create random networks
create_BBT <- function(nNodes_RN, p_branching){
  
  p_BBT <- p_branching/(2-p_branching) # for nNodes_RN >> 1
  
  downNode_RN <- level <- numeric(nNodes_RN)
  level[1] <- 1
  ind_childnode <- which(level==0)[1]
  
  while (sum(level==0)>0 & ind_childnode<=nNodes_RN){
    ind_parentnode <- ind_childnode - 1
    parent_level <- level[ind_parentnode]
    open_branches <- which(level==parent_level)
    for (i in open_branches){
      if (runif(1) <= p_BBT){ # branching
        downNode_RN[ind_childnode] <- i
        downNode_RN[ind_childnode + 1] <- i
        level[ind_childnode] <- level[i] + 1
        level[ind_childnode + 1] <- level[i] + 1
      } else { # non-branching
        downNode_RN[ind_childnode] <- i
        level[ind_childnode] <- level[i] + 1
      }
      ind_childnode <- which(level==0)[1]
    }
  }
  downNode_RN <- downNode_RN[1:nNodes_RN]
  level <- level[1:nNodes_RN]
  
  tmp <- numeric(max(level))
  for (i in 1:max(level)){
    tmp[i] <- sum(level==i)
  }
  
  W_RN <- spam(0,nNodes_RN,nNodes_RN)
  ind <- matrix(0,nNodes_RN*10000,2)
  k <- 1
  for (i in 1:nNodes_RN){
    if (downNode_RN[i] != 0){
      ind[k,] <- c(i, downNode_RN[i]) 
      k <- k+1
    }
  }
  ind <- ind[1:(k-1),]
  W_RN[ind] <- 1
  rm(ind)
  gc(verbose=F)
  
  Length_RN <- 1+numeric(nNodes_RN)
  
  # Upstream_RN : list containing IDs of all nodes upstream of each node (plus node itself)
  Upstream_RN <- vector("list",nNodes_RN)
  Nupstream_RN <- numeric(nNodes_RN)
  for (i in 1:nNodes_RN){
    UpOneLevel <- which(downNode_RN==i) # find reaches at one level upstream
    Upstream_RN[[i]] <- UpOneLevel      # add them to the list
    while (length(UpOneLevel)!=0) { # continue until there are no more reaches upstream
      ContinuePath <- UpOneLevel # jump 1 level above
      UpOneLevel <- which(downNode_RN %in% ContinuePath) # find reaches at one level upstream
      Upstream_RN[[i]] <- c(Upstream_RN[[i]],UpOneLevel) # add them to the list
    }
    Upstream_RN[[i]] <- c(Upstream_RN[[i]],i)
    Nupstream_RN[i] <- length(Upstream_RN[[i]])
  }
  
  # create AG network
  DegreeIn <- colSums(W_RN)
  DegreeOut <- rowSums(W_RN)
  Confluence <- DegreeIn>1
  Source <- DegreeIn==0
  SourceOrConfluence <- Source|Confluence
  ConfluenceNotOutlet <- Confluence&(downNode_RN!=0)
  ChannelHeads <- SourceOrConfluence  #Source|ConfluenceNotOutlet
  
  OutletNotChannelHead <- (downNode_RN==0)&(!ChannelHeads)
  IsNodeAG <- SourceOrConfluence|OutletNotChannelHead
  whichNodeAG <- which(IsNodeAG)
  
  nNodes_AG <- sum(IsNodeAG)
  Length_AG <- numeric(nNodes_AG)
  RN_to_AG <- numeric(nNodes_RN)
  reachID <- 1
  while (length(whichNodeAG) != 0){ # explore all AG Nodes
    i <- whichNodeAG[1] # select the first
    RN_to_AG[i] <- reachID 
    j <- downNode_RN[i] 
    Length_AG[reachID] <- Length_RN[i]
    tmp_length <- Length_RN[i]
    tmp <- NULL
    j0 <- j
    while (!IsNodeAG[j] && j!=0) {
      tmp <- c(tmp, j)
      tmp_length <-  tmp_length + Length_RN[j]
      j_old <- j
      j <- downNode_RN[j]} 
    RN_to_AG[tmp] <- reachID
    Length_AG[reachID] <- tmp_length
    reachID <- reachID + 1
    whichNodeAG <- whichNodeAG[-1]
  }
  
  
  AG_to_RN <- vector("list", nNodes_AG)
  for(i in 1:nNodes_AG) { 
    AG_to_RN[[i]] <- which(RN_to_AG==i) 
  }
  
  downNode_AG <- numeric(nNodes_AG)
  # W_AG <- sparseMatrix(i=1,j=1,x=0,dims=c(Nnodes_AG,Nnodes_AG))
  W_AG <- spam(0,nNodes_AG,nNodes_AG)
  ind <- matrix(0,nNodes_AG,2)
  for (i in 1:nNodes_RN){ 
    if (downNode_RN[i] != 0 && RN_to_AG[downNode_RN[i]] != RN_to_AG[i]) {
      downNode_AG[RN_to_AG[i]] <- RN_to_AG[downNode_RN[i]]
      #W_AG[RN_to_AG[i],DownNode_AG[RN_to_AG[i]]] <- 1
      ind[RN_to_AG[i],] <- c(RN_to_AG[i],downNode_AG[RN_to_AG[i]])
    }
    # contributing area of nodes at AG level
    # if (ChannelHeads[i]){
    #   A_AG[RN_to_AG[i]] <- A_RN[i]
    # } 
  }
  ind <- ind[-which(ind[,1]==0),]
  W_AG[ind] <- 1
  rm(ind)
  gc(verbose=F)
  
  Upstream_AG <- vector("list",nNodes_AG)
  Nupstream_AG <- numeric(nNodes_AG)
  for (i in 1:nNodes_AG){
    UpOneLevel <- which(downNode_AG==i) # find reaches at one level upstream
    Upstream_AG[[i]] <- UpOneLevel      # add them to the list
    while (length(UpOneLevel)!=0) { # continue until there are no more reaches upstream
      ContinuePath <- UpOneLevel # jump 1 level above
      UpOneLevel <- which(downNode_AG %in% ContinuePath) # find reaches at one level upstream
      Upstream_AG[[i]] <- c(Upstream_AG[[i]],UpOneLevel) # add them to the list
    }
    Upstream_AG[[i]] <- c(Upstream_AG[[i]],i)
    Nupstream_AG[i] <- length(Upstream_AG[[i]])
  }
  
  StreamOrder_AG <- numeric(nNodes_AG)
  for (i in 1:nNodes_AG){
    j <- order(Nupstream_AG)[i] # index that explores reaches in a downstream direction
    tmp <- which(downNode_AG==j) # set of reaches draining into j
    if (length(tmp)>0){
      IncreaseOrder <- sum(StreamOrder_AG[tmp]==max(StreamOrder_AG[tmp])) # check whether tmp reaches have the same stream order
      if (IncreaseOrder > 1) {
        StreamOrder_AG[j] <- 1 + max(StreamOrder_AG[tmp]) # if so, increase stream order
      } else {StreamOrder_AG[j] <- max(StreamOrder_AG[tmp])} # otherwise, keep previous stream order
    } else {StreamOrder_AG[j] <- 1} # if j is an headwater, impose StreamOrder = 1
  }
  
  Outlet_RN <- which(downNode_RN==0)
  Outlet_AG <- which(downNode_AG==0)
  
  RN <- list(nNodes=nNodes_RN, downNode=downNode_RN, nUpstream=Nupstream_RN, toAG=RN_to_AG, 
             upstream=Upstream_RN, W=W_RN, leng=Length_RN, outlet=Outlet_RN)
  AG <- list(nNodes=nNodes_AG, downNode=downNode_AG, leng=Length_AG, streamOrder=StreamOrder_AG, toRN=AG_to_RN, W=W_AG, 
             upstream=Upstream_AG, nUpstream=Nupstream_AG, toCM=1+numeric(nNodes_AG), outlet=Outlet_AG)
  BBT <- list(RN=RN,AG=AG)
  
  invisible(BBT)
}

create_RBN <- function(nNodes_RN, p_branch){
  
  Length_AG <- rgeom(2*nNodes_RN,p_branch)+1
  exitflag <- 1
  
  while (exitflag==1){
    foo <- cumsum(Length_AG) 
    thr <- which(foo >= nNodes_RN)[1]
    
    if (foo[thr]==nNodes_RN & (thr %% 2)==1){ # if the no. patches is exactly nNodes_RN and thr is uneven, then no correction needed
      nNodes_AG <- thr
      Length_AG <- Length_AG[1:thr]
      exitflag <- 0
    } else if (foo[thr]>nNodes_RN & (thr %% 2)==1){ # if too many patches, but thr is uneven -> reduce size of last patch
      # note that last patch cannot be 0, otherwise the previous case would have occurred
      nNodes_AG <- thr
      Length_AG <- Length_AG[1:thr]
      Length_AG[thr] <- nNodes_RN - sum(Length_AG[1:(thr-1)])
      exitflag <- 0
    } else if ((thr %% 2)==0){ # if thr is even, remove two last elements from Length_AG and try again with the spare values
      Length_AG <- Length_AG[-c(thr-1,thr)]
    }
  }
  
  downNode_AG <- numeric(nNodes_AG)
  W_AG <- spam(0,nNodes_AG,nNodes_AG)
  ind <- matrix(0,nNodes_AG*10000,2)
  to_attribute_up <- 1:nNodes_AG
  to_attribute_down <- 1:nNodes_AG
  upstream <- vector("list",nNodes_AG)
  for (i in 1:nNodes_AG){
    upstream[[i]] <- i
  }
  k <- 1
  while (length(to_attribute_up)>1){ #
    fromNode <- sample(to_attribute_up,2) # choose 2 random nodes that haven't been attributed a downstream direction
    sset <- setdiff(to_attribute_down,c(upstream[[fromNode[1]]],upstream[[fromNode[2]]]) )
    if (length(sset)==1){sset <- c(sset,sset)} # patch for length(sset)=1 (which would cause problems to sample)
    toNode <- sample(sset,1) 
    # possible downstream direction: one that is not already located upstream of a fromNode
    downNode_AG[fromNode] <- toNode
    ind[k,] <- c(fromNode[[1]],toNode)
    ind[k+1,] <- c(fromNode[[2]],toNode)
    k <- k+2
    upstream[[toNode]] <- c(upstream[[toNode]], upstream[[fromNode[1]]], upstream[[fromNode[2]]]) # add all upstream nodes 
    
    j <- downNode_AG[toNode]
    while (j != 0){ # add forbidden connections also to all catchments downstream!
      upstream[[j]] <- c(upstream[[j]], upstream[[toNode]])
      j <- downNode_AG[j]
    }
    
    to_attribute_up <- setdiff(to_attribute_up,fromNode)
    to_attribute_down <- setdiff(to_attribute_down,toNode)
    # the last node remaining in to_attribute is the outlet
  }
  ind <- ind[1:(k-1),]
  W_AG[ind] <- 1
  
  maxrow <- max(rowSums(as.dgCMatrix.spam(W_AG + t(W_AG))))
  isDAG <- det(Diagonal(nNodes_AG) - as.dgCMatrix.spam( t(W_AG)))
  rm(ind)
  gc(verbose=F)
  
  RN_to_AG <- NULL 
  for (i in 1:nNodes_AG){
    RN_to_AG <- c(RN_to_AG,rep(i,Length_AG[i]))
  } 
  AG_to_RN <- vector("list", nNodes_AG)
  for(i in 1:nNodes_AG) { # attribute river network pixels to fields of the AG_to_FD list 
    AG_to_RN[[i]] <- which(RN_to_AG==i) 
  }
  
  downNode_RN <- numeric(nNodes_RN)
  ind <- 1
  for (i in 1:nNodes_AG){
    if (Length_AG[i]>1){
      downNode_RN[ind:(ind+Length_AG[i]-2)] <- (ind+1):(ind+Length_AG[i]-1)
      downNode_RN[ind+Length_AG[i]-1] <- which(RN_to_AG==downNode_AG[i])[1]
    } else {
      downNode_RN[ind] <- which(RN_to_AG==downNode_AG[i])[1]
    }
    ind <- ind + Length_AG[i]
  }
  downNode_RN[is.na(downNode_RN)] <- 0
  
  W_RN <- spam(0,nNodes_RN,nNodes_RN)
  ind <- matrix(0,nNodes_RN*10000,2)
  k <- 1
  for (i in 1:nNodes_RN){
    if (downNode_RN[i] != 0){
      ind[k,] <- c(i, downNode_RN[i]) 
      k <- k+1
    }
  }
  ind <- ind[1:(k-1),]
  W_RN[ind] <- 1
  rm(ind)
  gc(verbose=F)
  
  Length_RN <- 1+numeric(nNodes_RN)
  #det(diag(nNodes_RN) - t(W_RN))
  
  # Upstream_RN : list containing IDs of all nodes upstream of each node (plus node itself)
  Upstream_RN <- vector("list",nNodes_RN)
  Nupstream_RN <- numeric(nNodes_RN)
  for (i in 1:nNodes_RN){
    UpOneLevel <- which(downNode_RN==i) # find reaches at one level upstream
    Upstream_RN[[i]] <- UpOneLevel      # add them to the list
    while (length(UpOneLevel)!=0) { # continue until there are no more reaches upstream
      ContinuePath <- UpOneLevel # jump 1 level above
      UpOneLevel <- which(downNode_RN %in% ContinuePath) # find reaches at one level upstream
      Upstream_RN[[i]] <- c(Upstream_RN[[i]],UpOneLevel) # add them to the list
    }
    Upstream_RN[[i]] <- c(Upstream_RN[[i]],i)
    Nupstream_RN[i] <- length(Upstream_RN[[i]])
  }
  
  
  Upstream_AG <- vector("list",nNodes_AG)
  Nupstream_AG <- numeric(nNodes_AG)
  for (i in 1:nNodes_AG){
    UpOneLevel <- which(downNode_AG==i) # find reaches at one level upstream
    Upstream_AG[[i]] <- UpOneLevel      # add them to the list
    while (length(UpOneLevel)!=0) { # continue until there are no more reaches upstream
      ContinuePath <- UpOneLevel # jump 1 level above
      UpOneLevel <- which(downNode_AG %in% ContinuePath) # find reaches at one level upstream
      Upstream_AG[[i]] <- c(Upstream_AG[[i]],UpOneLevel) # add them to the list
    }
    Upstream_AG[[i]] <- c(Upstream_AG[[i]],i)
    Nupstream_AG[i] <- length(Upstream_AG[[i]])
  }
  
  StreamOrder_AG <- numeric(nNodes_AG)
  for (i in 1:nNodes_AG){
    j <- order(Nupstream_AG)[i] # index that explores reaches in a downstream direction
    tmp <- which(downNode_AG==j) # set of reaches draining into j
    if (length(tmp)>0){
      IncreaseOrder <- sum(StreamOrder_AG[tmp]==max(StreamOrder_AG[tmp])) # check whether tmp reaches have the same stream order
      if (IncreaseOrder > 1) {
        StreamOrder_AG[j] <- 1 + max(StreamOrder_AG[tmp]) # if so, increase stream order
      } else {StreamOrder_AG[j] <- max(StreamOrder_AG[tmp])} # otherwise, keep previous stream order
    } else {StreamOrder_AG[j] <- 1} # if j is an headwater, impose StreamOrder = 1
  }
  
  Outlet_RN <- which(downNode_RN==0)
  Outlet_AG <- which(downNode_AG==0)
  
  RN <- list(nNodes=nNodes_RN, downNode=downNode_RN, nUpstream=Nupstream_RN, toAG=RN_to_AG, 
             upstream=Upstream_RN, W=W_RN, leng=Length_RN, outlet=Outlet_RN)
  AG <- list(nNodes=nNodes_AG, downNode=downNode_AG, leng=Length_AG, streamOrder=StreamOrder_AG, toRN=AG_to_RN, W=W_AG, 
             upstream=Upstream_AG, nUpstream=Nupstream_AG, toCM=1+numeric(nNodes_AG), outlet=Outlet_AG)
  brn <- list(RN=RN,AG=AG)
  
  invisible(brn)
}


# function to create point shapefile given coordinates
shp.point <- function(x,y,sname="shape"){
  n <- length(x)
  dd <- data.frame(Id=1:n,X=x,Y=y)
  ddTable <- data.frame(Id=c(1),Name=paste("Outlet",1:n,sep=""))
  ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
  write.shapefile(ddShapefile, sname, arcgis=T)
}

# a quick R function to write a shapefile
makeshape.r=function(sname="shape",n=1){
  xy=locator(n=n)
  points(xy)
  
  #Point
  dd <- data.frame(Id=1:n,X=xy$x,Y=xy$y)
  ddTable <- data.frame(Id=c(1),Name=paste("Outlet",1:n,sep=""))
  ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
  write.shapefile(ddShapefile, sname, arcgis=T)
}

# inherited from OCNet
initial_permutation <- function(DownNode){
  
  Outlet <- which(DownNode==0)
  NodesToExplore <- Outlet # start from outlets
  reverse_perm <- numeric(length(DownNode)) # build permutation vector from outlets to headwaters, then flip it
  
  k <- 0
  while (length(NodesToExplore)>0){ # continue until all the network has been explored
    k <- k + 1
    node <- NodesToExplore[1] # explore a node
    reverse_perm[k] <- node # assign position in the permutation vector
    NodesToExplore <- NodesToExplore[-1] # remove explored node
    UpNodes <- which(DownNode==node) # find nodes upstream of node
    while (length(UpNodes)>0){ # continue upstream until a headwater is found
      k <- k + 1
      node <- UpNodes[1] # explore first upstream node
      reverse_perm[k] <- node
      if (length(UpNodes)>1){ # if there is a bifurcation upstream, add the other upstream connections at the top of NodesToExplore
        NodesToExplore <- c(UpNodes[2:length(UpNodes)],NodesToExplore)
      }
      UpNodes <- which(DownNode==node)
    }
  }
  
  perm <- reverse_perm[length(DownNode):1] # flip permutation
  
  OutList = list(perm=perm,noDAG=0)
  
  invisible(OutList)
}

neigh <- function(dir) {
  mov <- c(0,0)
  switch(dir,
         {mov[1] <- 1; mov[2] <- 0},   # case 1 (E)
         {mov[1] <- 1; mov[2] <- -1},   # case 2 (NE) 
         {mov[1] <- 0; mov[2] <- -1},   # case 3 (N)
         {mov[1] <- -1; mov[2] <- -1},  # case 4 (NW)
         {mov[1] <- -1; mov[2] <- 0},  # case 5 (W)
         {mov[1] <- -1; mov[2] <- 1}, # case 6 (SW)
         {mov[1] <- 0; mov[2] <- 1},  # case 7 (S)
         {mov[1] <- 1; mov[2] <- 1})  # case 8 (SE)
  return(mov)
}

# variation of draw_thematic_OCN of OCNet
draw_thematic_catch <- function(catch,theme=NA*numeric(catch$AG$nNodes),
                                chooseAggregation=NULL,
                                discreteLevels=FALSE,
                                colLevels=NULL,
                                cutoff=FALSE,
                                colPalette=colorRampPalette(c("yellow","red","black")),
                                drawNodes=FALSE,
                                nodeType="upstream",
                                cex=2,
                                pch=21,
                                nanColor="#0099FF",
                                riverColor="#0099FF",
                                backgroundColor="#999999",
                                addLegend=TRUE,
                                zoom_xy_lim=NA){
  
  # initialization
  if (discreteLevels == FALSE) {
    if (is.null(colLevels)){
      
      minval <- min(theme[!(is.nan(theme))])
      maxval <- max(theme[!(is.nan(theme))])
      if (is.na(minval) & is.na(maxval)){
        minval <- 0; maxval <- 0;
      }
      N_colLevels <- 1000
      colLevels <- c(minval,maxval,N_colLevels)
    }
    minval <- colLevels[1]; maxval <- colLevels[2]; N_colLevels <- colLevels[3]
    
    if (minval==maxval) {maxval <- minval + 1}
    Breakpoints <- seq(minval,maxval,len = N_colLevels+1)
  } else if (discreteLevels == TRUE) {
    if (is.null(colLevels)){
      N_colLevels <- length(unique(theme[!is.nan(theme)]))
      Breakpoints <- c(sort(unique(theme[!is.nan(theme)])),2*max(theme[!is.nan(theme)]))
    } else {N_colLevels <- length(colLevels) - 1
    Breakpoints <- colLevels}}
  
  if (typeof(colPalette)=="closure") {
    colPalette <- colPalette(N_colLevels)
  } else if (typeof(colPalette)=="character") {
    colPalette <- colPalette[1:N_colLevels] }
  
  if (length(theme)!=catch$AG$nNodes){
    stop('theme has invalid length')
  }
  
  
  if (length(cex)>1 && length(cex) != length(theme)){
    stop('cex has invalid length')
  }
  
  if (length(pch)>1 && length(pch) != length(theme)){
    stop('pch has invalid length')
  }
  
  
  X <- catch$RN$X#[which( catch$FD$toCM %in% chooseCM )]
  Y <- catch$RN$Y#[which( catch$FD$toCM %in% chooseCM )]
  Xc <- catch$CM$XContour
  Yc <- catch$CM$YContour
  
  
  if (length(cex)==1){
    cex_vec <- cex*rep(1,length(theme))
  } else {cex_vec <- cex}
  
  if (length(pch)==1){
    pch_vec <- pch*rep(1,length(theme))
  } else {pch_vec <- pch}
  
  AvailableNodes <- setdiff(1:catch$RN$nNodes,catch$RN$outlet)
  
  if (all(is.na(zoom_xy_lim))){
    zoom_xy_lim <- numeric(4)
    zoom_xy_lim[1] <- min(Xc); zoom_xy_lim[2] <- max(Xc) 
    zoom_xy_lim[3] <- min(Yc); zoom_xy_lim[4] <- max(Yc) 
    showAxes <- FALSE; labX <- ""; labY  <- ""
  } else {showAxes <- TRUE; labX <- "x"; labY <- "y"}
  
  plot(c(zoom_xy_lim[1],zoom_xy_lim[2],zoom_xy_lim[2],zoom_xy_lim[1],zoom_xy_lim[1]), 
       c(zoom_xy_lim[3],zoom_xy_lim[3],zoom_xy_lim[4],zoom_xy_lim[4],zoom_xy_lim[3]),
       type="n",xlab=labX,ylab=labY,asp=1,axes=showAxes)
  
  xy_lim <- par("usr")
  
  if (!is.null(backgroundColor)){ # don't show background when zooming
    polygon(Xc,Yc,col=backgroundColor,lty=0)    
  }
  
  maxA <- max(catch$RN$A[which(X >= xy_lim[1]& X <= xy_lim[2] & Y >= xy_lim[3] & Y <= xy_lim[4])])
  
  for (i in AvailableNodes){
    
    reach <- catch$RN$toAGReach[i]
    if (catch$RN$A[i]>=catch$thrA & X[i] >= xy_lim[1]-1*catch$cellsize & X[i] <= xy_lim[2]+1*catch$cellsize &
        Y[i] >= xy_lim[3]-1*catch$cellsize & Y[i] <= xy_lim[4]+1*catch$cellsize ) {
      if (  ((is.nan(theme[reach])==TRUE | is.na(theme[reach])==TRUE)) ||
            ( cutoff==TRUE && (theme[reach] < min(Breakpoints) || theme[reach] > max(Breakpoints))) )  {
        hexcolor <- nanColor
      } else {
        val <- theme[reach]
        colvalue <- which(Breakpoints > val)[1] - 1  
        if (isTRUE(colvalue==0)) {colvalue <- 1}
        if (is.na(colvalue)) {colvalue <- N_colLevels}
        
        hexcolor <- colPalette[colvalue]
      }
      if (drawNodes==TRUE){
        lines(c(X[i],X[catch$RN$downNode[i]]),c(Y[i],Y[catch$RN$downNode[i]]),
              lwd=0.1+2.9*(catch$RN$A[i]/maxA)^0.5,col=riverColor)
      } else {
        lines(c(X[i],X[catch$RN$downNode[i]]),c(Y[i],Y[catch$RN$downNode[i]]),
              lwd=0.1+2.9*(catch$RN$A[i]/maxA)^0.5,col=hexcolor)}
    }
  }
  
  if (drawNodes==TRUE){
    for (i in 1:catch$AG$nNodes){
      if (is.nan(theme[i]) || (cutoff==TRUE && (theme[i] < min(Breakpoints) || theme[i] > max(Breakpoints)) )){
        hexcolor <- nanColor
      } else {
        colvalue <- which(Breakpoints > theme[i])[1] - 1
        if (isTRUE(colvalue==0)) {colvalue <- 1}
        if (is.na(colvalue)) {colvalue <- N_colLevels}
        hexcolor <- colPalette[colvalue]
      }
      if (nodeType=="upstream") {
        points(catch$AG$X[i],catch$AG$Y[i],bg=hexcolor,pch=pch_vec[i],cex=cex_vec[i])  
      } else if (nodeType=="downstream") {
        points(catch$AG$XReach[i],catch$AG$YReach[i],bg=hexcolor,pch=pch_vec[i],cex=cex_vec[i])  
      }
    }
  }
  
  if (addLegend) {
    if (discreteLevels==FALSE){
      tmp <- par()$plt[2]
      pr <- 1 - tmp
      image.plot(col=colPalette,legend.only=TRUE,zlim=c(minval,maxval),
                 smallplot=c(0.88, 0.9,par()$plt[3],par()$plt[4]))
    } else {
      if (is.null(colLevels)){
        str <- NULL
        for (level in 1:N_colLevels){
          str <- c(str, as.character(round(1000*Breakpoints[level])/1000) )}
      } else { 
        str <- vector(mode="character", N_colLevels)
        for (level in 1:(N_colLevels-1)){
          str[level] <- paste("[",as.character(round(1000*Breakpoints[level])/1000),"; ",
                              as.character(round(1000*Breakpoints[level+1])/1000),")",sep="")
        } 
        str[N_colLevels] <- paste("[",as.character(round(1000*Breakpoints[N_colLevels])/1000),"; ",
                                  as.character(round(1000*Breakpoints[N_colLevels+1])/1000),"]",sep="")
      }
      
      legend(x=1.01*max(Xc),y=max(Yc),
             str,fill=colPalette,ncol=ceiling(N_colLevels/20), xpd=TRUE, cex=0.8, bty="n")
    }
  }
  invisible(xy_lim)
  
}

