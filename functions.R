# ================ Functions ================

# draw boxplot over all types of data over entire time length ("grouped")
# type in ["m", "c", "r", "mue"]
bpEntirely <- function(df, type, selectedTraj = selectTrajValuesByType(type, df), mfrow = c(1,1)){
  ### Boxplot grouped Data types over time
  par(mfrow=mfrow)
  boxplot.matrix(as.matrix.data.frame(selectedTraj))
}

# plot boxplot for all slices of a type ("groups" combined)  
bpSlices <- function(df, type, selectedTraj = selectTrajValuesByType(type, df), nrTs=5){
  par(mfrow=c(1,nrTs))
  for (time in timeSlicing(df, nrTs = nrTs)){
    boxplot(selectedTraj[time,])
  }
}

# plot boxplot at certain time of certain type ("groups" combined)  
bpAtTime <- function(df, type, at, selectedTraj = selectTrajValuesByType(type, df)){
  selectedTraj <- selectedTraj[cvTime2Idx(df, at), ]
  par(mfrow=c(1,1))
  boxplot(selectedTraj)
}

# create static temporary network to check for cycles
createNet <- function(edges, nodes, trajectory, time){
  net <- network(edges, matrix.type="edgelist", vertex.attr = nodes[order(nodes$id),], 
                 loops=T, multiple=T, ignore.eval = F, directed = TRUE)
  ### Set network attributes
  #### Set persistent ids for vertices and edges
  set.network.attribute(netStat, 'vertex.pid', 'id')
  #set.network.attribute(netStat, 'edge.pid', 'reacID')
  #### Set Metadata information about network --> FIXME: ADD!
  #nameNetwork <- "Pulsed Feed Example"
  #subNetwork <- "1 organism, 3 compounds"
  #set.network.attribute(netStat, 'net.name', nameNetwork)
  #set.network.attribute(netStat, 'net.sub', subNetwork)
  
  ### Creating parameters for plot function
  net <- colorfyNet(net)
}

# 




# Cytoscape filter
filterCytoscape <- function(){
  
}




# filter edges with fluxxes in a certain range
filterFlux <- function(edgeDF, LT=NULL, UT=NULL, rmZeros=TRUE){
  if(is.null(UT)){tbDeletedUpper <- NULL}
  else{tbDeletedUpper <- which(edgeDF$flux > UT)}
  if(is.null(LT)){tbDeletedLower <- NULL}  
  else{tbDeletedLower <- which(edgeDF$flux < LT)}
  tbDeletedZeros <- NULL
  if(rmZeros){
    tbDeletedZeros <- which(edgeDF$flux == 0)
  }
  tbDeleted <- union(tbDeletedUpper, union(tbDeletedLower, tbDeletedZeros))
  if(length(tbDeleted) == 0){
    cat("No reduction, adapt thresholds if necessary!")
    return(edgeDF)
  }
  else{
    cat(sprintf("Reduction of %s edges, which were below %s or above %s", length(tbDeleted), LT, UT))
    return(edgeDF[-tbDeleted, ])
  }
}

# reduce edgeSet at a certain timestep() by providing a lower and upper threshold for flux-values
# return reduced edgeSet with corresponding flux-values
reduceEdgeSet <- function(traj=NULL, edges, time, lowerThresh=NULL, upperThresh=NULL, rmBackEdges = FALSE, 
                          firstCall = FALSE){

  
  #1 At the first call, the algorithm has to get the flux values from the trajectory and
  # implement this values into the edges-df and create new edges, if they are negative fluxes
  if(firstCall && !is.null(traj)){
    timeStep <- cvTime2Idx(traj, time)
    initialTraj <- selectTrajValuesByType("flux", traj, edgesDf = edges)[timeStep,]
    if(is.null(upperThresh)){tbDeletedUpper <- NULL}
    else{tbDeletedUpper <- which(abs(initialTraj) > upperThresh)}
    if(is.null(lowerThresh)){tbDeletedLower <- NULL}  
    else{tbDeletedLower <- which(abs(initialTraj) < lowerThresh)}
    tbDeleted <- c(tbDeletedUpper, tbDeletedLower)
    
    # if(length(tbDeleted) == 0){cat("No reduction took place, please adapt the thresholds!\n")}
    # # At the first call, the algorithm has to check for a cycle and create negFluxEdges
    if(length(tbDeleted) == 0){
      reducedTraj <- initialTraj
    }
    else{
      reducedTraj <- initialTraj[,-tbDeleted] 
    }
    reducedEdges <- edges[which(edges$reacID %in% names(reducedTraj)), ]
    
    # Find out the negativ pointing ones
    negFluxIndexes <- integer(0)
    for (i in 1:nrow(reducedEdges)){
      # reducedEdges$flux[i] <- reducedTraj[[reducedEdges$reacID[i]]][1] # given at the end
      # find "negative" fluxes and interchange $from and $to
      if(reducedEdges$sense[[i]] * reducedTraj[[reducedEdges$reacID[i]]][1] < 0){
        negFluxIndexes <- append(negFluxIndexes, i)
        oldFrom <- reducedEdges$from[[i]]
        oldTo <- reducedEdges$to[[i]]
        reducedEdges$from[[i]] <- oldTo
        reducedEdges$to[[i]] <- oldFrom
      }
      # make all flux values positive, as direction info is included in ($from $to) -pair
      reducedEdges$flux[i] <- abs(reducedTraj[[reducedEdges$reacID[i]]][1])
    }
    return(reducedEdges)
  }
  
  #2 Normal Iteration Step, NO FIRST CALL
  if(is.null(upperThresh)){tbDeletedUpper <- NULL}
  else{tbDeletedUpper <- which(edges$flux > upperThresh)}
  if(is.null(lowerThresh)){tbDeletedLower <- NULL}  
  else{tbDeletedLower <- which(edges$flux < lowerThresh)}
  tbDeleted <- c(tbDeletedUpper, tbDeletedLower)
  
  if(length(tbDeleted) == 0){
    cat("No reduction took place, please adapt the thresholds!\n")
    return(edges)  # in this case there is no checking for cyclic property
  }
  else if(length(tbDeleted) == nrow(edges)){
    cat("Complete reduction, no edge present anymore. Adapt the thresholds!\n")
    return(edges[-c(1:nrow(edges)), ])  # returns an empty df with just the colNames counted as 1 row
  }
  else{ # partly reduction
    reducedEdges <- edges[-tbDeleted,]
    return(reducedEdges)
  }
  
  
}

# quick static network from edgelist without optional vertex attributes
pidNetQuick <- function(edges){
  if(is_empty(edges)){
    return(network.initialize(0))
  }
  net <- network(edges, matrix.type="edgelist",  
                 loops=T, multiple=F, ignore.eval = F, directed = TRUE)
  # vertex.attr = nodes[order(nodes$id),],  --> might include in network
  # in diesem Netzwerk sind folglich die anderen nodes Attribute noch nicht vorhanden
  set.network.attribute(net, 'vertex.pid', 'vertex.names')  
  # ODER set.network.attribute(net, 'vertex.pid', 'id')
  set.network.attribute(net, 'edge.pid', 'reacID')
  return(net)
}

# basic static network from edgelist WITH optional vertex attributes
pidNet <- function(edges, nodes){
  if(is_empty(edges)){
    return(network.initialize(0))
  }
  net <- network(edges, matrix.type="edgelist", vertex.attr = nodes[order(nodes$id),], 
                 loops=T, multiple=T, ignore.eval = F, directed = TRUE)
  set.network.attribute(net, 'vertex.pid', 'id')
  set.network.attribute(net, 'edge.pid', 'reacID')
  return(net)
}

# One Iteration including reducing the edgeSet and then checking for cycle
# returns TRUE, if reduction was succesful (no more loops), or FALSE, if it still contains loops
iterationTopDown <- function(traj = NULL, edges, time, lowerThresh, upperThresh, rmBackEdges = FALSE, firstCall = FALSE){
  redEdges <- reduceEdgeSet(traj, edges, time, lowerThresh, upperThresh, rmBackEdges, firstCall)
  #1 No reduction occured AND not first Call
  if(!firstCall && nrow(redEdges) == nrow(edges) ){  
    
    return(list(hasCycle=TRUE, edges=redEdges))
  }
  #2 Complete reduction, all edges removed
  else if(is_empty(redEdges)){  
    return(list(hasCycle=TRUE, edges=redEdges)) # return empty edges-df
  }
  #3 Partly reduction OR no reduction at FIRST CALL
  isCyclic <- hasCycle(pidNetQuick(edges = redEdges))
  return(list(hasCycle=isCyclic, edges=redEdges))
}

# takes the entire edgeSet and removes edges according to flux-value thresholds
# INPUT:
# Set initial lower and upper thresholds. These will be adapted during iteration, unitl the is no more cycle in the network
# Set certain thresStep, to define the lowering of the upper threshold
# different modes define the stepwise mechanism to reduce edgeset
## modes in ["quantiles", "oneByOne"]
# return final network
topDownReduction <- function(traj=NULL, edges, time, initLT = -Inf, initUT = Inf, mode = "quantiles", nrQuant = NULL){
  initialSize <- nrow(edges)
  quantiles <- trajStats(df=traj, type = "flux", time = time, quantFrac = getQuantFrac(nrQuant))
  mean <- mean(unlist(selectTrajValuesByType(type = "flux", traj = traj, edgesDf = edges), use.names = FALSE))
  #1 first Iteration, also kind of prefiltering step with thresholds initLT and initUT
  firstIteration <- iterationTopDown(edges = edges, traj = traj, time = time, lowerThresh = initLT, upperThresh = initUT,
                   firstCall = TRUE) # first Call: firstCall = True
  isCyclic <- firstIteration$hasCycle 
  if(!isCyclic){
    cat(sprintf("Network at time = %s with lowerThreshold %s and upperThreshold %s has no cycles!\n", time, initLT, initUT))
    return(firstIteration$edges)
  }
  
  #2 Loop through quantiles and set them as lowerThreshold
  edges <- firstIteration$edges
  if(mode == "quantiles"){
    counter = 1 # counts the quantiles
    while (isCyclic == TRUE && counter <= length(quantiles)) {
      currentThres <- quantiles[[counter]]
      currentIterationResult <- iterationTopDown(edges = edges, time = time, 
                                                 lowerThresh = currentThres, upperThresh = Inf)
      isCyclic <- currentIterationResult$hasCycle
      edges <- currentIterationResult$edges
      counter <- counter + 1
    }
    breakQuant <- quantiles[[counter-1]]
    finalEdgeSet <- currentIterationResult$edges
    finalSize <- nrow(finalEdgeSet)
    sizeReduction <- initialSize - finalSize
    cat(sprintf("Network got acyclic at a flux threshold of %s", breakQuant))
  }
  return(list(edges = finalEdgeSet, threshold = breakQuant, sizeReduction = sizeReduction))
}

# input: nr of wanted quantiles 
# return:  customized vector of evenely distributed quantile fractions
getQuantFrac <- function(nrQuant){
  if(is.null(nrQuant)){return(NULL)}
  return(c(0:nrQuant)/nrQuant)
}

# time is really a time, eg. 0.102
bottomUpReduction <- function(traj, edges, time, initLT = -Inf, initUT = Inf){
  timeStep <- cvTime2Idx(traj, time)
  initialTraj <- selectTrajValuesByType("flux", traj, edgesDf = edges)[timeStep,]
  initialSize <- nrow(edges)
  
  if(is.null(initUT)){tbDeletedUpper <- NULL}
  else{tbDeletedUpper <- which(abs(initialTraj) > initUT)}
  if(is.null(initLT)){tbDeletedLower <- NULL}  
  else{tbDeletedLower <- which(abs(initialTraj) < initLT)}
  tbDeleted <- c(tbDeletedUpper, tbDeletedLower)
  if(length(tbDeleted) == 0){
    reducedTraj <- initialTraj
  }
  else{
    reducedTraj <- initialTraj[, -tbDeleted] 
    # FIXME: For modTraj* müsste man hier noch alle NA rausfiltern
  }
  filteredEdges <- edges[which(edges$reacID %in% names(reducedTraj)), ]
  # FIXME: For modTraj* müsste man hier noch alle NA rausfiltern
  
  # FIXME: for Schleife bräuchte man alles nicht mehr!
  # Find out the negativ pointing ones
  negFluxIndexes <- integer(0)
  for (i in 1:nrow(filteredEdges)){
    # reducedEdges$flux[i] <- reducedTraj[[reducedEdges$reacID[i]]][1] # given at the end
    # find "negative" fluxes and interchange $from and $to
    if(filteredEdges$sense[[i]] * initialTraj[[filteredEdges$reacID[i]]][1] < 0){
      negFluxIndexes <- append(negFluxIndexes, i)
      oldFrom <- filteredEdges$from[[i]]
      oldTo <- filteredEdges$to[[i]]
      filteredEdges$from[[i]] <- oldTo
      filteredEdges$to[[i]] <- oldFrom
    }
    # make all flux values positive, as direction info is included in ($from $to) -pair
    filteredEdges$flux[i] <- abs(initialTraj[[filteredEdges$reacID[i]]][1])
    # add id
  }
  # Order by flux value in reverse
  filteredEdges <- filteredEdges[order(filteredEdges$flux, decreasing = TRUE), ]
  
  
  # Start building the network
  edgesNew <- data.frame(character(0), character(0), character(0), numeric(0), numeric(0), stringsAsFactors = FALSE)
  names(edgesNew) <- names(filteredEdges)
  for(idx in 1:nrow(filteredEdges)){
    edgesNew[idx,] <- filteredEdges[idx,]
    isCyclic <- hasCycle(pidNetQuick(edges = edgesNew))
    if(isCyclic){
      edgesNew <- edgesNew[-idx,]
      break
    }
  }
  finalSize <- nrow(edgesNew)
  sizeReduction <- initialSize - finalSize
  return(list(hasCycle=isCyclic, edges=edgesNew, sizeReduction=sizeReduction))
}

# Analyse a trajectory/network, by taking timeslices and reduce the network at the 
# specific timepoints to acyclic graphs. The function returns all acyclic edge-Sets at
# the given timepoints, for further analysis
# possible mode: 1 out of c("bottomUp", "topDown")
# custom in paraTimeSteps is a vector which contains row indexes of traj!, not time values itself
netStructAnalysis <- function(traj=trajectory, edges=edges, nodes=nodes, 
                              paraTimeSteps = list(nrSlices=NULL, custom=NULL), 
                              name="", 
                              plot = FALSE, mode = "bottomUp", nrQuant = NULL){
  timeSteps = timeSlices(traj=traj, nrSlices=paraTimeSteps$nrSlices, custom = paraTimeSteps$custom)
  timeSlices <- sapply(timeSteps, cvIdx2Time, traj=traj)
  modes = c("bottomUp", "topDown")
  resultList <- vector(mode = "list", length = length(timeSlices))
  sizeReductionList <- vector(mode = "list", length = length(timeSlices))
  resultDF <- data.frame(time=timeSlices, edges=rep(NA, length(timeSlices)))
  for(idx in 1:length(timeSlices)){
    if(!mode %in% modes){
      stop(cat(sprintf("Wrong mode type provided. Must be within %s!\n", modes))) # FIXME
    }
    else if(mode == "bottomUp"){
      currentEdgeSet <- bottomUpReduction(traj = traj, edges = edges, time = timeSlices[idx])
      #resultDF$edges[[idx]] <- list(currentEdgeSet$edges)
      resultList[[idx]] <- currentEdgeSet$edges
      sizeReductionList[[idx]] <- currentEdgeSet$sizeReduction
    }
    else if(mode == "topDown"){
      currentEdgeSet <- topDownReduction(traj = traj, edges = edges, time = timeSlices[idx], 
                                       nrQuant = nrQuant)
      #resultDF$edges[[idx]] <- list(currentEdgeSet$edges)
      resultList[[idx]] <- currentEdgeSet$edges
      sizeReductionList[[idx]] <- currentEdgeSet$sizeReduction
    }
  }
  names(resultList) <- timeSlices
  names(sizeReductionList) <- timeSlices
  return(list(edgeList = resultList, sizeReduction = sizeReductionList))
}

# calculate betweenes of all nodes by package RBGL
# 1.step: importing graph(s) from cytoscape
# 2.step: calculate betweenes
# 3.step: export results to graph
calcBetweenes <- function(cytoNetSuid){
  
}



# slicing trajectory rows into equal chunks
# return: vector of time steps
timeSlicing <- function(traj=trajectory, nrTs = 5){
  timeSteps <- (nrow(trajectory) - nrow(trajectory) %% nrTs) / nrTs
  timeSlices <- c(1, (timeSteps * seq(1, nrTs-1))+1)
}

# set custom time slices for layout analysis or get automated ones by timeSlicing()
timeSlices <- function(traj, nrSlices=NULL, custom=NULL){
  if(is.null(nrSlices) && is.null(custom)){
    return(timeSlicing(traj = traj))
  }
  else if(is.null(nrSlices)){
    return(custom)
  }
  else if(is.null(custom)){
    return(timeSlicing(traj = traj, nrTs = nrSlices))
  }
  else{
    error("Either provide number of timeSlices or custom timepoints, not both at the same time!")
  }
  
}
  
  
  
# type must be in  ["m", "c", "r", "flux", "ie", "mue"] 
# if time == NULL: alle Werte
trajStats <- function(df, type, selectedTraj = selectTrajValuesByType(type, df), time=NULL, 
                      quantFrac = NULL){
  if(!is.null(time)){
    currentTraj <- selectedTraj[cvTime2Idx(selectedTraj, time),-1]
  }
  else{
    currentTraj <- selectedTraj[,-1]
  }
  if(is.null(quantFrac)){
    quantiles <- quantile(abs(unlist(currentTraj, use.names = FALSE))) # use standard qunatiles 1/4, 2/4, 3/4
  }
  else {
    quantiles <- quantile(abs(unlist(currentTraj, use.names = FALSE)), quantFrac)
  }
  return(quantiles)  
  
}

# # type must be in  ["m", "c", "r", "flux", "ie", "mue"] 
# trajStatsOld <- function(df, type, selectedTraj = selectTrajValuesByType(type, df), 
#                       quantFrac = c(1, 3)/4){
#   for (name in selectedTraj){
#     min <- min(selectedTraj[[name]])
#     quantiles <- quantile(selectedTraj[[name]], quantFrac)
#     mean <- mean(selectedTraj[[name]])
#     median <- median(selectedTraj[[name]])
#     max <- max(selectedTraj[[name]])
#     
#   }
# }

# recursive utility function for hasCycle() function
# Depth first search into the tree
# marks nodes as visited as going through
# after going through all out-neighbors of a node, this node will be removed from the recStack
recursiveCheck <- function(net, id, env){
  #print(parent.env(environment()))
  #print(sys.nframe())
  #print(names(environment()))
  #cat('visited_before_', sys.nframe(), ': ', visited, '\n')
  visTmp <- get("visited", envir = env)
  recStTmp <- get("recStack", envir = env)
  env$visited[id] <- TRUE
  env$recStack[id] <- TRUE
  #cat('visited_before_', sys.nframe(), ': ', env$visited, '\n')
  #cat('recStack_before_', sys.nframe(), ': ', env$recStack, '\n')
  for(neighborID in get.neighborhood(net, id, type = 'out')){
    if(env$visited[neighborID] == FALSE){
      if(recursiveCheck(net, neighborID, env) == TRUE){
        return(TRUE)
      }
    }
    else if(env$recStack[neighborID] == TRUE){  # the condition visited[neighborID] == TRUE is implicitly included
      return(TRUE)
    }
  }  
  # Current node needs to be poppen from recursion Stack before function ends
  env$recStack[id] <- FALSE
  #cat('visited_after_', sys.nframe(), ': ', env$visited, '\n')
  #cat('recStack_after_', sys.nframe(), ': ', env$recStack, '\n')
  
  return(FALSE)

}

# checks if network has cycles
# returns True, if it is cyclic
hasCycle <- function(net){
  cycleEnv <- new_environment()
  cycleEnv$visited <- rep(FALSE, network.size(net))
  cycleEnv$recStack <- rep(FALSE, network.size(net))
  
  storageEnv <- new.env()
  assign("visited", rep(FALSE, network.size(net)), envir = storageEnv)
  assign("recStack", rep(FALSE, network.size(net)), envir = storageEnv)
  
  
  #print(parent.env(environment()))
  #print(sys.nframe())
  #print(names(environment()))
  visited <- rep(FALSE, network.size(net))
  recStack <- rep(FALSE, network.size(net))
  # Perform a depth first search
  for(pid in get.vertex.pid(net, 1:network.size(net))){
    id <- get.vertex.id(net, pid) # necessary conversion to internal ids (basically numbered vertecis)
    if(get("visited", envir = storageEnv)[id] == FALSE){
      if(recursiveCheck(net, id, storageEnv) == TRUE){
        return(TRUE)
      }
      else{
        return(FALSE)
      }
    }
    
  }
}




colorfy <- function(column){
  distinctFeatures = unique(column)
  returnVec <- integer(length(column))
  for(i in 1:length(column)){
    for(j in 1:length(distinctFeatures)){
      if(column[i]==distinctFeatures[j]){
        returnVec[i] <- j
      }
    }
  }
  return(returnVec)
}

colorfySense <- function(reacIDs){
  returnVec <- integer(length(reacIDs))
  vecIdx = 1
  for(id in reacIDs){
    if(endsWith(id, "_c")){
      returnVec[vecIdx] <- 2
      vecIdx <- vecIdx + 1
    }
    else{
      returnVec[vecIdx] <- 1
      vecIdx <- vecIdx + 1
    }
  }
  return(returnVec)
}

# set "colNode" and "colEdge" als attribute in vertex/edge
colorfyNet <- function(net, vertexAttr = "type"){
  uniqueVertexTypes <- length(unique(get.vertex.attribute(net, vertexAttr)))
  colPalNodesType <- viridis_pal(begin = 0, end = 0.6)(uniqueVertexTypes)
  show_col(colPalNodesType)
  set.vertex.attribute(net, "colNode", colPalNodesType[colorfy(get.vertex.attribute(net, vertexAttr))])
  get.vertex.attribute(net, "colNode")
  
  colPalEdgeSense <- viridis_pal(begin = 0.8, end = 1)(2)
  show_col(colPalEdgeSense)
  set.edge.attribute(net, "colEdge", colPalEdgeSense[colorfySense(get.edge.attribute(net, "reacID"))])
  get.edge.attribute(net, "colEdge")
  return(net)
}

# Color node by type
colorNode <- function(typeVec){
  moCol = "green"
  compoundCol = "grey"
  colVec <- vector("character", length(typeVec))
  for(i in 1:length(typeVec)){
    if(typeVec[i]=="m"){
      colVec[i] <- moCol
    }
    else if(typeVec[i]=="c"){
      colVec[i] <- compoundCol
    }
    else{
      colVec[i] <- NA
    }
  }
  return(colVec)
}

#convert time to row-Index
cvTime2Idx <- function(traj=trajectory, time){
  which(near(traj$time, time, tol = 2000*(.Machine$double.eps^0.5)))
}

# convert row-Index to time
cvIdx2Time <- function(traj=trajectory, idx){
  traj$time[idx]
}

# Go one time steo further and adjust all evolutionary values in nodes and edges
getValuesAtTs <- function(nodes, edges, timestep){
  for(nInd in 1:nrow(nodes)){
    nodes$value[nInd] <- trajectory[timestep, nodes$id[nInd]]
  }
  for(eInd in 1:nrow(edges)){
    edges$flux[eInd] <- trajectory[timestep, edges$reacID[eInd]]
  }
  
  return(trajectory[timestep, "time"])
}



# Deactivate certain nodes

# Set dynamic edge Attributes which change over time
updateEdgeAttributes <- function(dynNet, traj = trajectory, edgeAttr='flux', defaultVal=0)
{
  # Initialisation, cb every attribute has to be defined for each tim eperiod
  # see ndtv.pdf, 7.1 / page 18
  activate.edge.attribute(dynNet, edgeAttr, defaultVal, onset=-Inf, terminus=Inf)
  reacIDs <- dynNet %e% 'reacID'
  for(ts in 1:nrow(traj)){
    currentVal <- traj[ts, reacIDs]
    if(ts != nrow(traj)){
      activate.edge.attribute(dynNet, edgeAttr, currentVal, onset = traj$time[ts], terminus = traj$time[ts+1])
    }
    else {
      activate.edge.attribute(dynNet, edgeAttr, currentVal, onset = traj$time[ts], terminus = traj$time[ts])
    }
  }
  return(dynNet)
}

# Set dynamic node Attributes which change over time
updateNodeAttributes <- function(dynNet, traj = trajectory, nodeAttr='value', defaultVal=0)
{
  # Initialisation, cb every attribute has to be defined for each tim eperiod
  # see ndtv.pdf, 7.1 / page 18
  activate.vertex.attribute(dynNet, nodeAttr, defaultVal, onset=-Inf, terminus=Inf)
  nodeIDs <- dynNet %v% 'id'
  for(ts in 1:nrow(traj)){
    currentVal <- traj[ts, nodeIDs]
    if(ts != nrow(traj)){
      activate.vertex.attribute(dynNet, nodeAttr, currentVal, onset = traj$time[ts], terminus = traj$time[ts+1])
    }
    else {
      activate.vertex.attribute(dynNet, nodeAttr, currentVal, onset = traj$time[ts], terminus = traj$time[ts])
    }
  }
  return(dynNet)
}

network.layout.animate.custom2 <- function(net, dist.mat = NULL, default.dist = NULL, 
                                           seed.coords = NULL, layout.par = list(), verbose=FALSE){
  fixedCoords <- plot(net1, interactive=T)
  return(cbind(x,y))
}

network.layout.animate.custom <- function(net, dist.mat = NULL, default.dist = NULL, 
                                          seed.coords = NULL, layout.par = list(), verbose=FALSE){
  return(customCoords)
}

# Select all columns of values of the trajectory data of a certain type and provide its prefix
# eg: for all microorganisms: pref= "m_"
# DOES NOT ALTER THE TRAJECTORY DATA ITSELF!!!
# If I would want to, I need assign("trajectory", alteredTraj, envor = .GlobalEnv)
selectTrajValuesByPref <- function(pref, traj = trajectory){
  return(traj[names(traj)[startsWith(names(traj), pref)]])
}

# select all values from trajectory file by one type of ["m", "c", "mc", "r", "flux", "ie", "mue"]
# FIXME: uses nodes from Global environment
selectTrajValuesByType <- function(type, traj = trajectory, edgesDf = edges, nodesDf = nodes){
  if(type == "m" || type == "c"){
    return(traj[nodesDf$id[startsWith(nodesDf$id, paste(type, "_", sep = ""))]])
  }
  else if(type == "mc"){
    boolVec <- !(startsWith(nodesDf$id, "i") | startsWith(nodesDf$id, "e"))
    return(traj[c("time", nodesDf$id[boolVec])])
  }
  
  else if(type == "r"){
    return(traj[c("time", edgesDf$reacID)])
  }
  else if(type == "mue"){
    return(traj[names(traj)[startsWith(names(traj), type)]])
  }
  else if(type == "flux"){
    boolVec <- !(startsWith(edgesDf$reacID, "i") | startsWith(edgesDf$reacID, "e"))
    return(traj[c("time", edgesDf$reacID[boolVec])])
  }
  else if(type == "ie"){
    boolVec <- startsWith(edgesDf$reacID, "i") | startsWith(edgesDf$reacID, "e")
    return(traj[edges$reacID[boolVec]])
  }
  else{
    stop("Wrong type provided. Must be within [\"m\", \"c\", \"r\", \"flux\", \"ie\",\"mue\"]")
  }
}

scaleValuesType <- function(type, traj = trajectory, minCut = min(traj), maxCut = max(traj), 
                            scRange = c(0,1), func = NULL){
  selectedVals <- selectTrajValuesByType(type, traj)
  selectedVals[selectedVals < minCut] <- NA
  selectedVals[selectedVals > maxCut] <- NA
  normalize(selectedVals, method = "range", range = scRange, margin = 2L)
}

# # scales node values; SAME SCALE for microorganisms and compounds!!! # FIXME: NOT CLEVER!
# scaleNodeValues <- function(traj=trajectory, minCut = NULL, maxCut = NULL, scRange= c(0.5, 5)){
#   trajCol2ScaleM <- selectTrajValuesByType("m", traj = traj)
#   trajCol2ScaleC <- selectTrajValuesByType("c", traj = traj)
#   if(is.null(minCut) == TRUE){
#     minCut <- min(min(trajCol2ScaleM), min(trajCol2ScaleC))
#   }
#   if(is.null(maxCut) == TRUE){
#     maxCut <- max(max(trajCol2ScaleM), max(trajCol2ScaleC))
#   }
#   trajCol2ScaleM[trajCol2ScaleM < minCut] <- NA # Wahrscheinlich besser auf 0 setze
#   trajCol2ScaleM[trajCol2ScaleM > maxCut] <- NA # Auf Inf setzen???! aber kacke für visualisierung
#   trajCol2ScaleC[trajCol2ScaleC < minCut] <- NA 
#   trajCol2ScaleC[trajCol2ScaleC > maxCut] <- NA
#   # join both horizontally --> Aufteilung vorher eigentlich SINNNNNNNLOS!!!
#   namesM <- names(trajCol2ScaleM)
#   namesC <- names(trajCol2ScaleC)
#   joinedNames = c(namesM, namesC)
#   traj[namesM] <- trajCol2ScaleM
#   traj[namesC] <- trajCol2ScaleC
#   # Normalize
#   traj[joinedNames] <- normalizeMe(traj[joinedNames], minCut= minCut, maxCut=maxCut, scRange = scRange)
#   print(paste("Wertebereich von [", minCut, ", ", maxCut, "]", 
#               " auf Wertebereich [", scRange[1], ", ", scRange[2], "] skaliert!", sep = ""))
#   return(traj)
# }

# minCut and maxCut can contain custom values; min/maxCut[1] refers to microorganisms ...[2] to compounds
# scRange can contain custom values; scRange[1,2] refers to microorganisms ...[3,4] to compounds
# VALUES NEED TO BE ABSOLUTE; NOT NEGATIVE ONES!!!
scaleNodeValues <- function(traj=trajectory, minCut = NULL, maxCut = NULL, scRange= c(0.5, 5, 0.5, 5)){
  trajCol2ScaleM <- selectTrajValuesByType("m", traj = traj)
  trajCol2ScaleC <- selectTrajValuesByType("c", traj = traj)
  if(is.null(minCut) == TRUE){
    minCutM <- min(trajCol2ScaleM)
    minCutC <- min(trajCol2ScaleC)
  }
  else {
    minCutM <- minCut[1]
    minCutC <- minCut[2]
  }
  if(is.null(maxCut) == TRUE){
    maxCutM <- max(trajCol2ScaleM)
    maxCutC <- max(trajCol2ScaleC)
  }
  else {
    maxCutM <- maxCut[1]
    maxCutC <- maxCut[2]
  }
  
  namesM <- names(trajCol2ScaleM)
  namesC <- names(trajCol2ScaleC)
  # Normalize
  traj[namesM] <- normalizeMe(traj[namesM], minCut= minCutM, maxCut=maxCutM, scRange = scRange[1:2])
  print(paste("MIKROORGANISMS: Wertebereich von [", minCutM, ", ", maxCutM, "]", 
              " auf Wertebereich [", scRange[1], ", ", scRange[2], "] skaliert!", sep = ""))
  traj[namesC] <- normalizeMe(traj[namesC], minCut= minCutC, maxCut=maxCutC, scRange = scRange[3:4])
  print(paste("COMPOUNDS: Wertebereich von [", minCutC, ", ", maxCutC, "]", 
              " auf Wertebereich [", scRange[3], ", ", scRange[4], "] skaliert!", sep = ""))
  return(traj)
}




normalizeMe <- function(dataFrame, minCut, maxCut, scRange){
  for (colName in names(dataFrame)){
    for (rowIdx in 1:length(dataFrame[[colName]])){
      if (dataFrame[[colName]][rowIdx] < minCut){
        dataFrame[[colName]][rowIdx] <- 0  # FIXME: How to deal with values smaller than minCut
      }
      else if (dataFrame[[colName]][rowIdx] > maxCut){
        dataFrame[[colName]][rowIdx] <- maxCut  # # FIXME: How to deal with values greater than max Cut
      }
      else{
        dataFrame[[colName]][rowIdx] <- ((dataFrame[[colName]][rowIdx] -minCut) /maxCut) *(scRange[2] - scRange[1]) +scRange[1] 
      }
    }
  }
  return(dataFrame)
  # dataFrame <- ((dataFrame-minCut)/maxCut)*(scRange[2] - scRange[1])+scRange[1]
}

scaleEdgeFlux <- function(traj=trajectory, minCut = NULL, maxCut = NULL, scRange= c(0.5, 5)){
  trajCol2ScaleR <- selectTrajValuesByType("r", traj = traj)
  if(is.null(minCut) == TRUE){
    minCut <- min(trajCol2ScaleR)
  }
  if(is.null(maxCut) == TRUE){
    maxCut <- max(trajCol2ScaleR)
  }
  
  trajCol2ScaleR[trajCol2ScaleR < minCut] <- NA # Wahrscheinlich besser auf 0 setze
  trajCol2ScaleR[trajCol2ScaleR > maxCut] <- NA # Auf Inf setzen???! aber kacke für visualisierung
  names2Scale <- names(trajCol2ScaleR)
  # Normalize
  traj[names2Scale] <- normalizeMe(traj[names2Scale], minCut= minCut, maxCut=maxCut, scRange = scRange)
  print(paste("Wertebereich von [", minCut, ", ", maxCut, "]", 
              " auf Wertebereich [", scRange[1], ", ", scRange[2], "] skaliert!", sep = ""))
  return(traj)
}

getFluxValSmallerThan <- function(id, val=0, traj=trajectory){
  foundRows <- which(traj[id] < val)
  return(traj[foundsRows, traj[c("time", id)]])
}

# find continuos stretchtes of time periods with negative values of value*sense 
# and save them with corresponding reacIDs
getInverseFlowPeriods <- function(traj=trajectory){
  edgeActList <- vector(mode = "list", length = length(edges$reacID))
  names(edgeActList) <- edges$reacID
  tsRows <- length(traj$time) # nr. of timesteps in trajectory data
  i <- 1  # reacID counter 
  while (i <= length(edges$reacID)){  
    currReacID <- colnames(traj[edges$reacID])[i] # current reacID as string
    currSense <- edges$sense[which(edges$reacID == currReacID)]
    negRows <- which(traj[currReacID]*currSense < 0) # all negative timesteps for each reaction
    #FIXME: | (traj[currReacID] <0 & currSense <0) --> falls Definition nicht immer RIchtung mo -> comp 
    if(length(negRows) == 0){
      edgeActList[[currReacID]] <- NULL  # Deletes the corresponding reacID -tag from the List, 
      # as there are no negative values
    }
    else if(length(negRows) != 0){ # There are negative Values for the i-th reacID ...
      currentVec <- c() # will contain all c(onset, terminus) pairs for every reacID
      j <- 1
      while(j <= length(negRows)){
        onset = negRows[j] # Zeitpunkt, bei dem negativer Stretch startet (ORDINAL!, nicht exakte Zeit)
        while(j < length(negRows) && (negRows[j] +1) == negRows[j+1]){
          j = j + 1
          #if(j == length(edges$reacID)){break} # letztes ELement in einem längeren Stretch enthalten
        }
        # hier kann es zu einer Überschreitung von terminus über maxIndex traj$time um +1 kommen, 
        # das ist aber so gewollt und wird in den Folgefunktionen entsprechend behandelt
        terminus <- negRows[j] + 1 
        j <- j+1
        currentVec <- c(currentVec, c(onset, terminus))
      }
      edgeActList[[currReacID]] <- currentVec
    }
    i <- i+1  # next reacID
  }
  return(edgeActList)
}

incorpNewNodesInTraj <- function(traj = trajectory, nodesName = c("i", "e"), defVal = 0){
  for(node in nodesName){
    traj <- insertColumn(traj, node, rep(defVal, nrow(traj)))
  }
  return(traj)
}

# Create inverse due to negative values in trajectory
# Parallel to creating the new edges the will be activated in net which is a network-Object
# Input: getInverseFlowPeriods(trajectory)
incorpInverseEdges <- function(traj = trajectory, invFp = getInverseFlowPeriods(traj), edgeL = edges){
  # replace $sense now with $excretion, as all values in Trajectory are positive now, so the meaning of this column changes slightly
  names(edgeL) <- gsub("sense", "excretion", names(edgeL), fixed = TRUE)
  edgeL$excretion <- rep(TRUE, nrow(edgeL))
  
  for(i in 1:length(invFp)){  # loop through all edges which have to be duplicated
    currReacID <- names(invFp)[i]
    # create inverse edges
    # 1. add to edges
    edgeIdx <- which(edgeL$reacID==currReacID)
    newEdge <- edgeL[edgeIdx, ]
    oldFrom <- newEdge$from
    newEdge$from <- newEdge$to
    newEdge$to <- oldFrom
    newReacID <- paste(currReacID, "c", sep = "_")
    newEdge$reacID <- newReacID
    newEdge$excretion <- FALSE
    edgeL <- insertRow(edgeL, newEdge, edgeIdx +1)
    # 2. add to trajectory
    newCol <- traj[[currReacID]]
    traj <- insertColumn(traj, newReacID, newCol, currReacID)
    
    # 3. set values = 0 for forward and inverse edges respectively in the calculated intervals from invFp
    # oldTerminus: basically the onset of each forward flux period, so end to each reverse flow period
    oldTerminus <- NA # Initialization: wenn reverse edge direkt bei timeIdx == 1 startet
    firstOnsetTimeIdx <- invFp[[currReacID]][1]
    # condition for first iteration, to find out if reverse flux already an timeIdx = 1
    if(firstOnsetTimeIdx > 1){
      oldTerminus <- 1
    }
    
    runIdx <- 1
    resetValue <- NA
    while (runIdx < length(invFp[[currReacID]])) {
      onset <- invFp[[currReacID]][runIdx]
      terminus <- invFp[[currReacID]][runIdx+1]
      
      # set flux of all forward edges to zero in [onset, terminus-1]
      traj[onset:(terminus-1), which(names(traj) == currReacID)] <- resetValue
      
      # set flux of all forward edges to zero in [oldTerminus, onset -1)
      # check first if(is.na(terminusOld)), um zu überprüfen, ob im ersten Durchgang eine Löschung erfolgen muss
      if(!is.na(oldTerminus)){
        traj[oldTerminus:(onset-1), which(names(traj) == newReacID)] <- resetValue
      }
      # Set oldTerminus to current terminus
      oldTerminus <- terminus
      # Jump to next onset|terminus pair in invFp[[currReacID]]
      runIdx <- runIdx + 2
    }
    if(oldTerminus<=nrow(traj)){ # bedeute es kommt noch mal eine "forward period" ...
      traj[onset:(terminus-1), which(names(traj) == currReacID)] <- resetValue
      traj[oldTerminus:nrow(traj), which(names(traj) == newReacID)] <- resetValue
    } 
  # make all flux values positive
  fluxNames <- edgeL$reacID
  #backwardFluxNames <- paste(names(invFp), "c", sep = "_")
  #fluxNames <- c(forwardFluxNames, backwardFluxNames)
  traj[fluxNames] <- abs(traj[fluxNames])
  }
  return(list("edges" = edgeL, "trajectory" = traj))
}

# Activation of edges according to return of getInverseFlowPeriods()
edgeActivation <- function(net, traj = trajectory, invFp){
  # All edges which do not change direction during evolution
  #edgeNamesConst <- setdiff(names(traj), c(names(invFp), "time"))
  #edgeNamesVarying <- names(invFp)
  #startTime = traj$time[1]
  #endTime = traj$time[nrow(traj)]
  for (i in 1:length(invFp)) {
    currReacID <- names(invFp)[i]
    forwarEdgeID <- currReacID
    inverseEdgeID <- paste(currReacID, "c", sep = "_")
    runIdx <- 1
    lalala <- length(invFp[[currReacID]])
    while (runIdx < length(invFp[[currReacID]])) {
      # Deactivate Forward facing Edges and activate inverse onces in the calculated intervals from invFp
      onsetIdx <- invFp[[currReacID]][runIdx]
      terminusIdx <- invFp[[currReacID]][runIdx+1]
      
      onset <- traj$time[onsetIdx]
      if(terminusIdx <= length(traj$time)){
        terminus <- traj$time[terminusIdx]
        deactivate.edges(net, onset = onset, terminus = terminus, e = get.edge.id(net, forwarEdgeID))
        activate.edges(net, onset = onset, terminus = terminus, e = get.edge.id(net, inverseEdgeID))
      }
      else{ # terminusIdx überschreitet den timeIdx: onset --> terminus Schreibweise kann nichtmehr verwendet werden
        lastTimePoint <- traj$time[nrow(traj)]
        deactivate.edges(net, onset = onset, length = lastTimePoint - onset, e = get.edge.id(net, forwarEdgeID))
        activate.edges(net, onset = onset, length = lastTimePoint - onset, e = get.edge.id(net, inverseEdgeID))
      }
        
      runIdx <- runIdx + 2
    }
  }
  return(net)
}

#converts the \s to /.
repath <- function(promptMessage = NA) {
  if (is.na(promptMessage)){
    x <- readline(prompt = "Enter Path to directory containing one Dataset: ")
  }
  else{
    x <- readline(prompt = promptMessage)
  }
  xa <- gsub('\\\\', '/', x)
  writeClipboard(paste(xa, collapse=""))
  cat('Here\'s your de-windowsified path. (It\'s also on the clipboard.)\n', xa, '\n')
  return(xa)
}

# inserts row at certain position in a data.frame
insertRow <- function(existingDF, newrow, r) {
  if (r<=nrow(existingDF)){
    existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  }
  existingDF[r,] <- newrow
  existingDF
}

# inserts Column with specified colName and optional position (after can either be an colIdx or colName) in a data.frame
insertColumn <- function(existingDF, colName, newCol, after = length(existingDF)){
  if (is.character(after) && length(after) == 1){ # checks iff after is a String
    nameIdx <- which(names(existingDF) == after)
  }  
  else if(is.numeric(after)){
    nameIdx = after
  }
  
  if (nameIdx != length(existingDF)){
    newOrder <- c(names(existingDF)[1:nameIdx], colName, names(existingDF)[(nameIdx+1):length(existingDF)])
  }
  else {
    newOrder <- c(names(existingDF), colName)
  }
  existingDF[[colName]] <- newCol
  existingDF <- existingDF[, newOrder]
  return(existingDF)
}

#plot Static Network with Standard parameters
plotStatNet <- function(net = netStat, paraList = statPlotPara){
  par(mfrow=c(1,1))
  plot(net, vertex.col= paraList$vertexCol, vertex.cex = paraList$vertexCex, 
       label = get.vertex.attribute(net, paraList$vertexLabel), edge.col = paraList$edgeCol, 
       edge.lwd = paraList$edgeLwd, arrowhead.cex =paraList$arrowheadCex, 
       edge.curve = paraList$edgeCurve, xlab=paraList$xLab, ylab=paraList$yLab, 
       edge.label=paraList$edgeLabel,
       interactive=paraList$interAct, displaylabels=paraList$displayLab, 
       main = get.network.attribute(net, "net.name"), sub = get.network.attribute(net, "net.sub"))
}

#plot Several Static networks
plotStatNets <- function(){}

# create backbone nodeSet from all nodeSets
createBackBoneNodeIds <- function(listNodes){
  unique(unlist(sapply(listNodes, function(x){x$id})))
}

# create backbone edgeSet from all edgeSets
createBackBoneEdgeIDs <- function(listEdges){
  allEdgeIDs <- sapply(listEdges, function(x){x$reacID})
  unique(unlist(allEdgeIDs))
}


#create multiple networks from List of EdgeSets
#FIXME: Might add plot Function
# FIXME: Do not include an additional nodes filtering step
createStatNets <- function(listSets, nodes, title, plot=FALSE, remNodes=FALSE, traj=trajectory, LTnodes = 0.1){
  netList <- vector(mode = "list", length = length(listSets))
  names(netList) <- names(listSets)
  for(idx in 1:length(listSets)){
    # currentTime <- attr(listSets[[idx]], "time")
    currentTime <- as.numeric(names(listSets)[idx])
    currentNodes <- nodes
    currentNodes$value <- unlist(sapply(currentNodes$id, function(x, time, traj){
      if(x %in% c("i", "e")){ # workaround, da in trajectory keine values für "i" und "e" hinterlegt sind
        return(0)
      }
      traj[cvTime2Idx(traj = traj, time = time),][[x]]
      }, time=currentTime, traj=traj, USE.NAMES = FALSE, simplify = TRUE)) # currentNodes$value wird mit tatsächlichen Konzentrationswerten gefüllt
    # FIXME: WARUM NULL BEI INFLUENT UND EFFLUENT node, und warum muss ich unlist machen????!!!
    # SHouldn't occur wenn nicht pulsedfeed
    
    if(remNodes){ # do not add nodes to the relative network, if they are not within the underlying edgeSet
      uniqueNodeIDs <- union(listSets[[idx]]$from, listSets[[idx]]$to)
      currentNodes <- redNodeSet(nodesDF = currentNodes, ids = uniqueNodeIDs)
    }
    
    # The following LTnodes filtering step is useless here, as the will be included from the edgelist by using the
    # network constructor with matrix.type = "edgelist" anyway
    # currentNodes <- currentNodes[which(currentNodes$value > LTnodes), ]
    
    currentNet <- network(listSets[[idx]], matrix.type="edgelist", vertex.attr = currentNodes[order(currentNodes$id),], 
                       loops=T, multiple=T, ignore.eval = F, directed = TRUE)
    # Create pids
    set.network.attribute(currentNet, 'vertex.pid', 'id')
    set.network.attribute(currentNet, 'edge.pid', 'reacID')
    #### Set Metadata information about network
    set.network.attribute(currentNet, 'net.sub', sprintf("t = %.2f", currentTime))
    currentNet <- colorfyNet(currentNet)
    netList[[idx]] <- currentNet
  }
  return(netList)
}

# Create nodeSets from edgeSets, so that nodes are kicked out, which are not present within the 
# relative edgeSet
createNodeSets <- function(listSets, nodes, traj=trajectory, exportCyt = FALSE){
  nodeList <- vector(mode = "list", length = length(listSets))
  names(nodeList) <- names(nodes)
  
  for(idx in 1:length(listSets)){
    # currentTime <- attr(listSets[[idx]], "time")
    currentTime <- as.numeric(names(listSets)[idx])
    currentNodes <- nodes
    currentNodes$value <- unlist(sapply(currentNodes$id, function(x, time, traj){
      if(x %in% c("i", "e")){ # workaround, da in trajectory keine values für "i" und "e" hinterlegt sind
        return(0)
      }
      traj[cvTime2Idx(traj = traj, time = time),][[x]]
    }, time=currentTime, traj=traj, USE.NAMES = FALSE, simplify = TRUE)) # currentNodes$value wird mit tatsächlichen Konzentrationswerten gefüllt
    # FIXME: WARUM NULL BEI INFLUENT UND EFFLUENT node, und warum muss ich unlist machen????!!!
    # SHouldn't occur wenn nicht pulsedfeed
    if("from" %in% names(listSets[[idx]]) && "to" %in% names(listSets[[idx]])){
      uniqueNodeIDs <- union(listSets[[idx]]$from, listSets[[idx]]$to)
    }
    else if("source" %in% names(listSets[[idx]]) && "target" %in% names(listSets[[idx]])){
      uniqueNodeIDs <- unique(union(listSets[[idx]]$source, listSets[[idx]]$target))
    }
    else{
      error("Column Identifier of start and end node id of an edge has to be either 
            'from', 'to' or 'source', 'target'. Please adapt to that effect!")
    }
      
    currentNodes <- redNodeSet(nodesDF = currentNodes, ids = uniqueNodeIDs)
    nodeList[[idx]] <- currentNodes
  }
  return(nodeList)
}

# replace names(edges) from "from" -> "source" AND "to" -> "target"
replEdgeNames <- function(charVec){
  charVec <- gsub("from", "source", charVec, fixed = TRUE)
  charVec <- gsub("to", "target", charVec, fixed = TRUE)
  return(charVec)
}

# from normal edge data.frame format to iGraph format
replEdgeNamesToI <- function(charVec){
  if(any(c("from", "to", "flux") %in% charVec)){
    charVec <- gsub("from", "source", charVec, fixed = TRUE)
    charVec <- gsub("to", "target", charVec, fixed = TRUE)
    charVec <- gsub("flux", "weight", charVec, fixed = TRUE)
  }  
  return(charVec)
}

# from iGraph edge data.frame format to normal
replEdgeNamesFromI <- function(charVec){
  if(any(c("source", "target", "weight") %in% charVec)){
    charVec <- gsub("source", "from", charVec, fixed = TRUE)
    charVec <- gsub("target", "to", charVec, fixed = TRUE)
    charVec <- gsub("weight", "flux", charVec, fixed = TRUE)
  }
  return(charVec)
}

# from normal OR iGraph edge data.frame format to Cyto format
replEdgeNamesToCyto <- function(charVec){
  if(any(c("from", "to", "weight") %in% charVec)){
    charVec <- gsub("from", "source", charVec, fixed = TRUE)
    charVec <- gsub("to", "target", charVec, fixed = TRUE)
    charVec <- gsub("weight", "flux", charVec, fixed = TRUE)
  }  
  return(charVec)
}

# from normal node data.frame format to iGraph format
# "id" als 'code-name'
replNodeNamesToI <- function(charVec){
  if("id" %in% charVec){
    charVec <- gsub("name", "shortName", charVec, fixed = TRUE)
    charVec <- gsub("id", "name", charVec, fixed = TRUE)
  }
  return(charVec)
}

# from iGraph node data.frame format to normal
# "shortName" als 'code-name'
replNodeNamesFromI <- function(charVec){
  if("shortName" %in% charVec){
    charVec <- gsub("name", "id", charVec, fixed = TRUE)
    charVec <- gsub("shortName", "name", charVec, fixed = TRUE)
  }
  return(charVec)
}

# replace names(edges) of edgesDF "from" -> "source" AND "to" -> "target"
replEdgeNamesOfEdgesDF <- function(edgeDF){
  names(edgeDF) <- gsub("from", "source", names(edgeDF), fixed = TRUE)
  names(edgeDF) <- gsub("to", "target", names(edgeDF), fixed = TRUE)
  return(edgeDF)
}

# reduce nodes by removing certain ones
# tbr: nodes-ids, which should be removed
# returns reduced nodeSet
redNodeSet <- function(nodesDF, ids){
  redNodes <- nodesDF[which(nodes$id %in% ids), ]
}

# filter out certain nodes by node ids or by type
# optional filter out certain edges by id
# excludeIE: exclude Inflow and Outflow nodes and corresponding edges
# in Default mode, just IE nodes and edges will be removed
# returns list of new Dataframes and of the ids of the removed nodes and edges
excludeNodesAndEdges <- function(nodesDF, edgesDF, nodeIDs=NULL, edgeIDs=NULL, nodeType = NULL, 
                                excludeIE = TRUE, output = "n"){
  # Change headers if needed
  names(edgesDF) <- replEdgeNamesFromI(names(edgesDF))
  names(nodesDF) <- replNodeNamesFromI(names(nodesDF))
  outputNormal <- c("n", "normal")
  outputIgraph <- c("i", "igraph", "iGraph")
  
  # Identify Rows in nodesDF to be excluded
  nodesRowNrTbEx <- integer(0)
  if(excludeIE){
    nodesRowNrTbEx<- union(nodesRowNrTbEx, which(nodesDF$type %in% "ie"))
  }
  if(!is.null(nodeType)){ # FIXME: Should be possible on top of excludeIE = TRUE, so rather APPEND nodesRowNrTbEx
    nodesRowNrTbEx <- union(nodesRowNrTbEx, which(nodesDF$type %in% nodeType))
  }
  if(!is.null(nodeIDs)){
    nodesRowNrTbEx <- union(nodesRowNrTbEx, which(nodesDF$id %in% nodeIDs))
  }
  nodeIDsTbEx <- nodesDF$id[nodesRowNrTbEx]
  
  # Identify Rows in edgesDF to be excluded
  edgesRowNrTbEx <- integer(0)
  if(excludeIE){  # FIXME: Should not be necessary, if corresponding "i" and "e" nodes exist in nodesDF
    edgesRowNrTbEx <- union(edgesRowNrTbEx, which(edgesDF$from %in% "ie" | edgesDF$to %in% "ie"))
  }
  if(!is.null(edgeIDs)){
    edgesRowNrTbEx <- union(edgesRowNrTbEx, which(edgesDF$reacID %in% edgeIDs))
  }
  if(!length(nodesRowNrTbEx) == 0){
    edgesRowNrTbEx <- union(edgesRowNrTbEx, which(edgesDF$from %in% nodeIDsTbEx | edgesDF$to %in% nodeIDsTbEx))
  }
  edgesReacIDsTbEx <- edgesDF$reacID[edgesRowNrTbEx]
  
  newNodes <- nodesDF[-nodesRowNrTbEx,]
  newEdges <- edgesDF[-edgesRowNrTbEx,]
  if(output %in% outputIgraph){
    names(newNodes) <- replEdgeNamesToI(names(newNodes))
    names(newEdges) <- replEdgeNamesToI(names(newEdges))
    message("Output in iGraph compatible format!")
  }
  else if(output %in% outputNormal){
    message("Output in normal (function.R compatible) format!")
  }
  else{
    error(paste0(output, " is not a valid value for the function parameter 'output'."))
  }
  cat(sprintf("Excluded %s/%s nodes and %s/%s edges from the initial dataframes!\n", 
              length(nodeIDsTbEx), nrow(nodesDF), length(edgesReacIDsTbEx), nrow(edgesDF)))
  return(list(nodes = newNodes, edges = newEdges, rmNodeIds = nodeIDsTbEx, rmEdgeReacIds = edgesReacIDsTbEx))
}


# import network data into cytoscape for one network, adapt layout and vis-props, calc betweenness ...
importNetworkToCytoscape <- function(cytoNodeSet, cytoEdgeSet, title = names(cytoEdgeSet), overallTitle){
  cytoNet <- createNetworkFromDataFrames(cytoNodeSet, cytoEdgeSet, 
                                                    title=title, 
                                                    collection=overallTitle)
  getNetworkSuid()
  # createGroup(group.name = 'Microorganism', nodes = 'm', nodes.by.col = 'type')
  # createGroup(group.name = 'Compounds', nodes = 'c', nodes.by.col = 'type')
  
  # Define Visual Properties
  setVisualStyle('Custom')
  setNodeShapeMapping(table.column = 'type', table.column.values = c('m', 'c'),
                      c('ELLIPSE', 'DIAMOND'), style.name = 'Custom')
  setNodeColorMapping(table.column = 'type', table.column.values = c('m', 'c'), mapping.type = 'd',
                      c('#42C77B', '#D469F3'), style.name = 'Custom')
  setNodeTooltipMapping(table.column = 'value', style.name = 'Custom')
  setEdgeLineWidthMapping(table.column = 'flux', mapping.type = 'c', 
                          table.column.values = c(min(cytoEdgeSets[[idx]]$flux), max(cytoEdgeSets[[idx]]$flux)),
                          widths = c(0.5,5),
                          style.name = 'Custom')
  setEdgeTooltipMapping(table.column = 'flux', style.name = 'Custom')
  layoutNetwork('hierarchical')
  # setNodeSizeMapping(table.column = 'value', mapping.type = 'c', 
  #                    table.column.values = c(min(cytoNodeSets[[idx]]$value), max(cytoNodeSets[[idx]]$value)),
  #                    sizes = c(5,30),
  #                    style.name = 'Custom')
  
  # setNodeBorderColorMapping(table.column = 'type', table.column.values = c('m', 'c'), mapping.type = 'd',
  #                           +                           c('#42C77B', '#D469F3'), style.name = 'Custom')
  
  
  # Set additional network attributes
  networkAttr <- data.frame(size = NULL, stringsAsFactors = FALSE)
  
  # longest and shortest path
  # longest: take leaves() function of graph and provide parameter degree.dir = "in" for starting layer and degree.dir = "out" as end layer
  ## take longest path
  ## take shortest path
  
  # set additional node attributes
  
  # betweenes with RBGL
  nelNet <- createGraphFromNetwork.adapted()
  betweennessBrandes <- brandes.betweenness.centrality(nelNet)
  ## Betweenness on nodes
  absBetweennessCentralityVertices <- t(betweennessBrandes$betweenness.centrality.vertices)
  dfBetweennessNodes <- data.frame(absBetweennessCentralityVertices, stringsAsFactors = FALSE)
  # FIXME: Is that right with stringAsFactors = FALSE ?!
  # add other betweenness parameters
  nVert <- ncol(betweennessBrandes$relative.betweenness.centrality.vertices)
  dfBetweennessNodes$relBetweennessCentralityVertices <- betweennessBrandes$relative.betweenness.centrality.vertices[1:nVert]
  ### load node betweenness into cytoscape with shortNames as keys on both sides
  loadTableData(dfBetweennessNodes, data.key.column = "row.names", 
                table = "node", table.key.column = "name")
  
  ## Betweenness on edges
  betweennessCentralityEdges <- t(betweennessBrandes$betweenness.centrality.edges)
  nEdges <- length(edgeData(nelNet))
  namesEdges <- vector(length = nEdges)
  for (idx in 1:nEdges){
    namesEdges[idx] <- edgeData(nelNet)[[idx]]$reacID
  }
  dfBetweennessEdges <- data.frame(betweennessCentralityEdges)
  rownames(dfBetweennessEdges) <- namesEdges
  ### load edge betweenness into cytoscape with shortNames as keys on both sides
  loadTableData(dfBetweennessEdges, data.key.column = "row.names", 
                table = "edge", table.key.column = "reacID")
}




# set visualization intervalls automatically
getInterval <- function(traj = trajectory, n = 20){
  nrRows <- nrow(traj)
  
}



# set slice.par
setSlicePar <- function(traj = trajectory, startV = traj$time[1], interV = getInterval(traj),
                        aggregateDur = 1, ruleV = 'any'){
  
}

# sliceFunction node-value
sliceFuncNode <- function(slc = slice, attr = "value", scFactor = 1){
  get.vertex.attribute(slice, attr) * scFactor
}

# Create lookup data frame, which contains name, longname, compartment and formula of every compound
importLookupDF <- function(){
  # Reactivate row below, if user should choose the path dynamically
  # pathLookupTxt = repath(promptMessage = "Please enter you absolute path to the lookup *.csv file: \n")
  pathLookupTxt = "C:/Users/rauj/Documents/microbialsim/simulatedTrajectory_20-Dec-2019 11_57_36_Case_3_0.01_init_parallel/cLookup_part.txt"
  lookupDF <- read.csv(pathLookupTxt, header=T, as.is=T)
}

# categorize with Rodrigo's lookup *.txt
# Two Step approach: First using longnames, then using shortnames
categorizeLookupDF <- function(lookupDF){
  folderDir <- "C:/Users/rauj/Documents/microbialsim/lookupData"
  carbonDF <- read.table(paste(folderDir, "carbonSources.txt", sep = '/'), header=F, as.is=T, sep = "\n")
  nitrogenDF <- read.csv(paste(folderDir, "nitrogenSources.txt", sep = '/'), header=F, as.is=T, sep = "\n")
  phosphorusDF <- read.csv(paste(folderDir, "phosphorusSources.txt", sep = '/'), header=F, as.is=T, sep = "\n")
  sulfurDF <- read.csv(paste(folderDir, "sulfurSources.txt", sep = '/'), header=F, as.is=T, sep = "\n")
  lookupDF$category <- rep(NA, nrow(lookupDF))
  # using longName
  commonGroundCarbon <- intersect(lookupDF$longName, carbonDF[,1])
  commonGroundNitrogen <- intersect(lookupDF$longName, nitrogenDF[,1])
  commonGroundPhosphorus <- intersect(lookupDF$longName, phosphorusDF[,1])
  commonGroundSulfur <- intersect(lookupDF$longName, sulfurDF[,1])
  for (idx in 1:nrow(lookupDF)){
    # FIXME: Decide on one if the two options
    # A: using %in%
    # B: using str_contains(x, pattern)
    # if (lookupDF$longName[idx] %in% commonGroundCarbon){
    #   lookupDF$category[idx] <- "carbon"
    # }
    # else if (lookupDF$longName[idx] %in% commonGroundNitrogen){
    #   lookupDF$category[idx] <- "nitrogen"
    # }
    # else if (lookupDF$longName[idx] %in% commonGroundPhosphorus){
    #   lookupDF$category[idx] <- "phosphorus"
    # }
    # else if (lookupDF$longName[idx] %in% commonGroundSulfur){
    #   lookupDF$category[idx] <- "sulfur"
    # }
    
    if (str_contains(commonGroundCarbon, lookupDF$longName[idx], ignore.case = TRUE)){
      lookupDF$category[idx] <- "carbon"
    }
    else if (str_contains(commonGroundNitrogen, lookupDF$longName[idx], ignore.case = TRUE)){
      lookupDF$category[idx] <- "nitrogen"
    }
    else if (str_contains(commonGroundPhosphorus, lookupDF$longName[idx], ignore.case = TRUE)){
      lookupDF$category[idx] <- "phosphorus"
    }
    else if (str_contains(commonGroundSulfur, lookupDF$longName[idx], ignore.case = TRUE)){
      lookupDF$category[idx] <- "sulfur"
    }
  }
  
  # using shortname (which is $name in lookupDF)
  # first create copy of shortnames whith "[e]" removed
  modShortNames <- sapply(lookupDF$name, gsub, pattern = "\\[e\\]", replacement = "")
  for (idx in 1:nrow(lookupDF)){
    if (str_contains(commonGroundCarbon, modShortNames[idx], ignore.case = TRUE)){
      lookupDF$category[idx] <- "carbon"
    }
    else if (str_contains(commonGroundNitrogen, modShortNames[idx], ignore.case = TRUE)){
      lookupDF$category[idx] <- "nitrogen"
    }
    else if (str_contains(commonGroundPhosphorus, modShortNames[idx], ignore.case = TRUE)){
      lookupDF$category[idx] <- "phosphorus"
    }
    else if (str_contains(commonGroundSulfur, modShortNames[idx], ignore.case = TRUE)){
      lookupDF$category[idx] <- "sulfur"
    }
  }
  return(lookupDF)
  
}

# Create lookup list which maps names (short names) to longnames
listLongNames <- function(lookupCompoundsDF = missing_arg()){
  if (is_missing(lookupCompoundsDF)){
    lookupCompoundsDF <- importLookupDF()
  }
  # get name (shortnames) as keys
  keys <- lookupCompoundsDF$name
  values <- lookupCompoundsDF$longName
  lookUpLongnames <- as.list(values)
  names(lookUpLongnames) <- keys
  return(lookUpLongnames)
}

# get Longname from shortname, eg: h[e] -> Proton
getLongName <- function(shortName, listWithLongNames = listLongNames()){
  listWithLongNames[[shortName]]
}

# get shortname from id, eg: c_64 -> h[e]
getShortName <- function(id, nodesDF = nodes){
  nodesDF$name[which(nodesDF$id == id)]
} 

# get id from shortname, h[e] -> c_64
getId <- function(shortName, nodesDF = nodes){
  nodesDF$id[which(nodesDF$name == shortName)]
}
# adaption of Rcy3 function createIgraphFromNetwork(), which uses colum "id" instead of "name" as identifier for nodes
createIgraphFromNetwork.adapted <- function (network = NULL){
  suid = getNetworkSuid(title = network)
  cyedges <- getTableColumns("edge", network = suid)
  cynodes <- getTableColumns("node", network = suid)
  if (!"source" %in% colnames(cyedges) || (!"target" %in% 
                                           colnames(cyedges))) {
    st = data.frame(do.call("rbind", strsplit(cyedges$name, 
                                              " \\(.*\\) ")))
    colnames(st) <- c("source", "target")
    cyedges <- cbind(st, cyedges)
  }
  colnames(cyedges)[colnames(cyedges) == "source"] <- "from"
  colnames(cyedges)[colnames(cyedges) == "target"] <- "to"
  cyedges2 = cbind(cyedges[c("from", "to")], cyedges[, 
                                                     !(names(cyedges) %in% c("from", "to"))])
  # ADAPTION: replaced ["name"] to ["id"]
  cynodes2 = cbind(cynodes["id"], cynodes[, !(names(cynodes) == 
                                                "id")])
  igraph::graph_from_data_frame(cyedges2, directed = TRUE, 
                                vertices = cynodes2)
}

# include adapted createIgraphFromNetwork() function call
createGraphFromNetwork.adapted <- function (network = NULL){
  ig <- createIgraphFromNetwork.adapted(network)
  g <- igraph::igraph.to.graphNEL(ig)
  return(g)
}

# sums trajectory values of certain type over entire time
sumTrajValuesByType <- function(traj = trajectory, type, edgesDF = edges, ordered = TRUE){
  selectedTraj <- selectTrajValuesByType(type = type, traj = traj, edgesDf = edgesDF)
  summedTraj <- colSums(selectedTraj)
  if (ordered){
    summedTraj <- summedTraj[order(summedTraj, decreasing = TRUE)]
  }
  return(summedTraj)
}

# get biomass name from id
getBiomassName <- function(id, nodesDF = nodes){
  nodesDF$name[which(nodes$id == id)]
}

# calculate discrete derivative of trajectory
# absDerivTraj: contains abs. derivations: delta(c) = c(k) - c(k-1)
# relDerivTraj: contains rel. derivations: delta(c) / delta(t)
# --> here $time will show delta(t) #FIXME: Maybe change?!
# return will include absolute derivation and relative derivation
derivTraj <- function(traj){
  absDerivTraj <- traj
  relDerivTraj <- traj 
  # Idee: Veränderung ist mit der jeweiligen Zeitdifferenz gewichtet, 
  # also Änderung/delta(t), da nicht alle delta(t) gleich groß 
  for (idx in 2:nrow(traj)){
    absDerivTraj[idx-1, -1] <- traj[idx, -1] - traj[idx-1, -1]
    deltaT <- traj[idx, 1] - traj[idx-1, 1] # first calc delta(t)
    relDerivTraj[idx-1, -1] <- (traj[idx, -1] - traj[idx-1, -1]) / deltaT
  }
  # Delete last timestamp, as no change can be calculated here
  absDerivTraj <- absDerivTraj[-nrow(traj),]
  relDerivTraj <- relDerivTraj[-nrow(traj),]
  return(list(abs = absDerivTraj, rel = relDerivTraj))
}

# sum up trajectory consumption and production rates per time step
# IDEE: Man könnte sogar noch sd oder quantile ausrechnen, und so eine besser Aussage über die Streuung zulassen ..
# return: data.fram(time=*, production=*, consumption=*)
integrateConsumpRates <- function(derivTraj){
  integDeriv <- data.frame(time = derivTraj$time, production = NA, consumption = NA, overall = NA)
  # get rid of time column
  workingTraj <- derivTraj[,-which(colnames(derivTraj) == "time")]
  result <- apply(workingTraj, 1, function(x){
    production <- sum(x[which(x >= 0)])
    consumption <- sum(x[which(x <= 0)])
    overall <- production + consumption
    return(c(production, consumption, overall))
  })
  resultTransp <- t(result)
  integDeriv$production <- resultTransp[, 1]
  integDeriv$consumption <- resultTransp[, 2]
  integDeriv$overall <- resultTransp[, 3]
  return(integDeriv)
}

# integrate fluxes over time, distinguishing ecretion and uptake fluxes
integrateFluxes <- function(traj, edgeDF){
  resultDF <- data.frame(time = traj$time, excretion.sum = NA, uptake.sum = NA,
                         excretion.mean = NA, uptake.mean = NA)
  # get rid of time column
  workingTraj <- traj[,-which(colnames(traj) == "time")]
  result <- apply(workingTraj, 1, function(x){
    excretion.sum <- sum(x[edgeDF$reacID[which(edgeDF$excretion)]], na.rm = TRUE)
    uptake.sum <- sum(x[edgeDF$reacID[which(!edgeDF$excretion)]], na.rm = TRUE)
    excretion.mean <- mean(x[edgeDF$reacID[which(edgeDF$excretion)]], na.rm = TRUE)
    uptake.mean <- mean(x[edgeDF$reacID[which(!edgeDF$excretion)]], na.rm = TRUE)
    return(c(excretion.sum, uptake.sum, excretion.mean, uptake.mean))
  })
  resultTransp <- t(result)
  resultDF$excretion.sum <- resultTransp[, 1]
  resultDF$uptake.sum <- resultTransp[, 2]
  resultDF$excretion.mean <- resultTransp[, 3]
  resultDF$uptake.mean <- resultTransp[, 4]
  return(resultDF)
}

# checks if choosen dataset is %in% c("b", "batch")
isBatch <- function(string){
  if (string %in% c("b", "batch")){
    TRUE
  }
  else if(string %in% c("c", "chemostat", "cstr")){
    FALSE
  }
  else{
    stop(paste0(string, " is not a valid descriptor", collapse = ""))
  }
}

# applies stype custom to current graph in cytoscape
# FIXME: needs edgesDF right now, because Linewidthmapping checks edges$flux for suitable mapping range
# --> should be possible from within the cytoNet, so that there is no need for further function parameter
applyCustomCytoStyle <- function(edgesDF){
  setVisualStyle('Custom')
  setNodeShapeMapping(table.column = 'type', table.column.values = c('m', 'c'),
                      c('ELLIPSE', 'DIAMOND'), style.name = 'Custom')
  setNodeColorMapping(table.column = 'type', table.column.values = c('m', 'c'), mapping.type = 'd',
                      c('#42C77B', '#D469F3'), style.name = 'Custom')
  setNodeTooltipMapping(table.column = 'value', style.name = 'Custom')
  setEdgeLineWidthMapping(table.column = 'flux', mapping.type = 'c',
                          table.column.values = c(min(edgesDF$flux), max(edgesDF$flux)),
                          widths = c(0.5,5),
                          style.name = 'Custom')
  setEdgeTooltipMapping(table.column = 'flux', style.name = 'Custom')
  layoutNetwork('hierarchical')
}

# fill in Values from trajectory at specific timeIdx
# NOTE: Excludes IE values by design!
fillInValues <- function(traj=trajectory, timeIdx, edgesDF = edges, nodesDF = nodes, 
                         rmEdges = TRUE, rmNodes = TRUE, addMue = FALSE, rmZero=FALSE){
  # slice trajectory at timeIdx
  slicedTrajConc <- selectTrajValuesByType("mc", traj=traj, edgesDf = edgesDF, nodesDf = nodesDF)[timeIdx,]
  slicedTrajFlux <- selectTrajValuesByType("flux", traj=traj, edgesDf = edgesDF, nodesDf = nodesDF)[timeIdx,]
  if(addMue){
    slicedTrajMue <- selectTrajValuesByType("mue", traj=traj, edgesDf = edgesDF, nodesDf = nodesDF)[timeIdx,]
  }
  # fill in values in edge and node data.frames
  nodesDF$value <- unlist(sapply(nodesDF$id, function(x, traj){
    traj[[x]]
  }, traj=slicedTrajConc, USE.NAMES = FALSE, simplify = TRUE)) 
  edgesDF$flux <- unlist(sapply(edgesDF$reacID, function(x, traj){
    traj[[x]]
  }, traj=slicedTrajFlux, USE.NAMES = FALSE, simplify = TRUE))
  if(addMue){
    rowNrsMo <- which(nodesDF$type == "m")
    nodesDF$mue <- rep(NA, nrow(nodesDF))
    nodesDF$mue[rowNrsMo] <- unlist(sapply(nodesDF$id[rowNrsMo], function(x, traj){
      moNumber <- gsub("m_", "", x, fixed = TRUE)
      traj[[paste0("mue_", moNumber, collapse = "")]]
    }, traj=slicedTrajMue, USE.NAMES = FALSE, simplify = TRUE)) 
  }
  # Remove all Edges, if their value is NA
  if(rmEdges){
    edgesDF <- edgesDF[which(!is.na(edgesDF$flux)), ]
  }
  if(rmNodes){
    nodesDF <- nodesDF[which(!is.na(nodesDF$value)), ]
  }
  
  # Remove all nodes and edges with value == 0
  if(rmZero){
    edgesDF <- edgesDF[which(edgesDF$flux>0), ]
  }
  if(rmZero){
    nodesDF <- nodesDF[which(nodesDF$value>0), ]
  }
  
  return(list(nodes = nodesDF, edges = edgesDF))
}


# function to extract flux and concentration information into nodes and edges data.frames
# names(nodes) = c("name", "type", "shortName", "value")
# names(edges) = c("from"/"source", "to"/"target", reacID, "type", "weight")
# return edges in igraph-Format with "source", "target", "weight"
# retrun nodes in igraph-Format with $name and $shortName instead of $id and $name
# output in c("i", "igraph", "iGraph", "n", "normal")
# --> defines output format
extractFullSlices <- function(traj, edgesDF, nodesDF, timeSlices = NA, output = "i", addMue = FALSE, rmZero = FALSE){
  outputNormal <- c("n", "normal")
  outputIgraph <- c("i", "igraph", "iGraph")
  # Adjust names of edgesDF and nodesDF to normal format, in case the input is partly
  # already in iGraphFormat
  names(edgesDF) <- replEdgeNamesFromI(names(edgesDF))
  names(nodesDF) <- replNodeNamesFromI(names(nodesDF))
  
  # check if timeSlices are explicitly given or have to automatically be calculated
  if(!is.na(timeSlices)){
    # evaluate timeSlices to timeIdxs
    timeIdxs <- unlist(sapply(timeSlices, cvTime2Idx, traj=traj))
  }
  else{
    # FIXME: Insert automatically generation of timeslices likt in netStructAnalysis
  }
  
  # initialize list with resulting data.frame
  slicedGraphInfo <- lapply(vector(mode = 'list', length(timeSlices)), function(x){
    x <- vector(mode='list', 2)})
  names(slicedGraphInfo) <- timeSlices # FIXME: Sollen da wirklich die Zeiten als keys rein?
  for (i in 1:length(timeIdxs)){
    modDFs <- fillInValues(traj = traj, timeIdx = timeIdxs[i], 
                           nodesDF = nodesDF, edgesDF = edgesDF, addMue = addMue, rmZero = rmZero)
    # check output format
    if(output %in% outputIgraph){
      names(modDFs$edges) <- replEdgeNamesToI(names(modDFs$edges))
      names(modDFs$nodes) <- replEdgeNamesToI(names(modDFs$nodes))
      message("Output in normal (function.R compatible) format!")
    }
    else if(output %in% outputNormal){
      message("Output in iGraph compatible format!")
    }
    else{
      error(paste0(output, " is not a valid value for the function parameter 'output'."))
    }
    # fill list with current result at timeIdxs[i]
    slicedGraphInfo[[i]] <- modDFs
  }
  return(slicedGraphInfo)
}  
  




# time is really a time, eg. 0.102
bottomUpReductionSlice <- function(edges, initLT = -Inf, initUT = Inf){
  initialSize <- nrow(edges)
  # In case of empty edgeSet
  if(initialSize == 0){
    return(list(hasCycle=FALSE, edges=edges, sizeReduction=0, rmEdgeIDs = vector(mode = "character", 0L)))
  }
  
  # Initialize edgesNew
  edgesOldModes <- unlist(sapply(names(edges), function(x){
    mode(edges[[x]])
  }))
  listNamesWithModes <- sapply(edgesOldModes, function(x){
    vector(mode = x)
  })
  edgesNew <- data.frame(listNamesWithModes, stringsAsFactors = FALSE)
  rmEdges <- vector(mode = "character", 0L)
  
  # Order by flux value in reverse
  edgesOrdered <- edges[order(edges$flux, decreasing = TRUE), ]
  
  # Start building the network and LINEARIZING the network
  if(nrow(edges)>0){
    currentLengthEdgesNew <- 1  # Initialization
    for(idx in 1:nrow(edgesOrdered)){
      edgesNew[currentLengthEdgesNew,] <- edgesOrdered[idx,]
      #isCyclic <- hasCycle(pidNetQuick(edges = edgesNew))
      #IG <- graph_from_data_frame(edgesNew)
      isDag <- is.dag(graph_from_data_frame(edgesNew))
      #if(isCyclic){
      if(!isDag){
        currentEdgeID <- edgesNew$reacID[currentLengthEdgesNew]
        edgesNew <- edgesNew[-currentLengthEdgesNew,]
        rmEdges[idx - nrow(edgesNew)] <- currentEdgeID
        
        # break # FIXME: Shouldn't I rather continue and try to add more?!
      }
      else{
        currentLengthEdgesNew <- currentLengthEdgesNew +1
      }
    }
  }
  
  finalSize <- nrow(edgesNew)
  sizeReduction <- initialSize - finalSize
  return(list(hasCycle=!isDag, edges=edgesNew, sizeReduction=sizeReduction, rmEdgeIDs = rmEdges))
}


# Analyse a trajectory/network, by taking timeslices and reduce the network at the 
# specific timepoints to acyclic graphs. The function returns all acyclic edge-Sets at
# the given timepoints, for further analysis
# possible mode: 1 out of c("bottomUp", "topDown")
# custom in paraTimeSteps is a vector which contains row indexes of traj!, not time values itself
linearizeSlices <- function(slices, mode = "bottomUp", output = "n"){
  outputNormal <- c("n", "normal")
  outputIgraph <- c("i", "igraph", "iGraph")
  
  timeSlices <- names(slices)
  modes = c("bottomUp", "topDown")

  resultList <- lapply(vector(mode = 'list', length(timeSlices)), function(x){
    x <- vector(mode='list', 3)})
  names(resultList) <- timeSlices
  
  
  sizeReductionList <- vector(mode = "list", length = length(timeSlices))
  rmEdgeIDs <- vector(mode = "list", length = length(timeSlices))
  names(sizeReductionList) <- timeSlices
  resultDF <- data.frame(time=timeSlices, edges=rep(NA, length(timeSlices)))
  for(idx in 1:length(timeSlices)){
    if(!mode %in% modes){
      stop(cat(sprintf("Wrong mode type provided. Must be within %s!\n", modes))) # FIXME
    }
    else if(mode == "bottomUp"){
      # replace edge headers to normal, in case they are in igraph format
      inputEdges <- slices[[timeSlices[idx]]]$edges 
      # FIXME: Why does this work? iSlices[[meta$slice[idx]]] does not work, only if 
      # meta$slice[idx] is turned to character like: iSlices[[as.character(meta$slice[idx])]]
      names(inputEdges) <- replEdgeNamesFromI(names(inputEdges))
      # LINEARIZATION STEP
      currentEdgeSet <- bottomUpReductionSlice(inputEdges)
      if(output %in% outputIgraph){
        names(currentEdgeSet$edges) <- replEdgeNamesToI(names(currentEdgeSet$edges))
        message("Output in iGraph compatible format!")
      }
      else if(output %in% outputNormal){
        message("Output in normal (function.R compatible) format!")
      }
      else{
        error(paste0(output, " is not a valid value for the function parameter 'output'."))
      }
      resultList[[idx]] <- list(edges=currentEdgeSet$edges, sizeReduction=currentEdgeSet$sizeReduction, rmEdgeIDs = currentEdgeSet$rmEdgeIDs)
    }
    else if(mode == "topDown"){
      # FIXME: Not yet implemented for slices
      # currentEdgeSet <- topDownReduction(traj = traj, edges = edges, time = timeSlices[idx], 
      #                                    nrQuant = nrQuant)
      # #resultDF$edges[[idx]] <- list(currentEdgeSet$edges)
      # resultList[[idx]] <- currentEdgeSet$edges
      # sizeReductionList[[idx]] <- currentEdgeSet$sizeReduction
    }
  }
  return(linSlices = resultList)
}





# Copy Betweenness to try some stuff
copyBetweenness <- function (graph, v = V(graph), directed = TRUE, weights = NULL, 
                             nobigint = TRUE, normalized = FALSE) 
{
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  v <- as.igraph.vs(graph, v)
  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) {
    weights <- E(graph)$weight
  }
  if (!is.null(weights) && any(!is.na(weights))) {
    weights <- as.numeric(weights)
  }
  else {
    weights <- NULL
  }
  on.exit(.Call(C_R_igraph_finalizer))
  res <- .Call(C_R_igraph_betweenness, graph, v - 1, as.logical(directed), 
               weights, as.logical(nobigint))
  if (normalized) {
    vc <- vcount(graph)
    if (is_directed(graph) && directed) {
      res <- res/(vc * vc - 3 * vc + 2)
    }
    else {
      res <- 2 * res/(vc * vc - 3 * vc + 2)
    }
  }
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    names(res) <- V(graph)$name[v]
  }
  res
}

# ommit certain compounds, definded by a char. vector of short naems
omitComps <- function(compNames, toBeOmitted){
  compNames[!compNames %in% toBeOmitted]
}


# iGraph Stats: Calculate graph related stats on slices, such as diamter, centrality, ...
# input either linearized slices (linearized = TRUE) or non-linear, which is the default
iGraphStats <- function(iSlices, iSlicesStats, iNodes, linearized = FALSE){
  sumOverallEdgesFull <- 0
  for(idx in 1:length(iSlices)){
    cEdges <- iSlices[[idx]]$edges
    names(cEdges) <- replEdgeNamesToI(names(cEdges))
    cGraph <- graph_from_data_frame(d = cEdges, vertices = iNodes)
    
    # Standard graph Parameters for further use
    # FIXME: Does it make sense to collect these here? maybe not, and rather separate
    sumOverallEdgesFull <- sumOverallEdgesFull + nrow(cEdges)
    #names(cEdges) <- replEdgeNamesFromI(names(cEdges))
    flux <- E(cGraph)$weight
    excretion <- E(cGraph)$excretion
      
    # Diameters weighted and not (flux value as edge weight already implemented in data structure)
    dia <- diameter(cGraph, weights = NA)
    diaPath <- get_diameter(cGraph, weights = NA)
    wDia <- diameter(cGraph)
    wDiaPath <- get_diameter(cGraph)
    
    # eigen_centrality
    #eigCentr <- eigen_centrality(cGraph, directed = TRUE, weights = NA)
    #wEigCentr <- eigen_centrality(cGraph, directed = TRUE)
    # --> yields warning from C Script: At centrality.c:362 :Weighted directed graph in eigenvector centrality
    
    # page rank
    pageRank <- page_rank(cGraph, directed = TRUE, weights = NA)
    wPageRank <- page_rank(cGraph, directed = TRUE)
    
    # Betweenness and edge_betweenness
    vertBetween <- igraph::betweenness(graph = cGraph, weights = NA, normalized = TRUE)
    #wVertBetween <- betweenness(graph = cGraph, normalized = TRUE)
    # --> Error in betweenness(graph = cGraph, normalized = TRUE) : At centrality.c:1604 : Weight vector must be positive, Invalid value
    edgeBetween <- edge_betweenness(graph = cGraph, weights = NA)
    #wEdgeBetween <- edge_betweenness(graph = cGraph)
    # --> this is were the crash appears
    
    degreeCentrIn <- centr_degree(cGraph, mode="in")
    degreeCentrOut <- centr_degree(cGraph, mode="out")
    degreeCentrTot <- centr_degree(cGraph, mode="total")

    if(linearized){
      # FIRST: longest weighted path, meaning the path with the highest sum of fluxes
      # therefore all fluxes have to be turned negative
      if(gsize(cGraph) > 0){
        E(cGraph)$weight <- -E(cGraph)$weight
      }
      # FIXME: just mode "in" works, for the others it detects negative loops???!?!
      maxFluxPath <- distances(cGraph, mode = "in")
      # SECOND just distance; to do so, each weight has to have the same value
      E(cGraph)$weight <- rep(-1L, gsize(cGraph))
      # FIXME: this just works with "out" and "in"
      distOut <- distances(cGraph, mode = "out")
      distIn <- distances(cGraph, mode = "in")
      # fill in results
      iSlicesStats[[idx]] <- list(flux=flux, excretion=excretion, dia=dia, diaPath=diaPath, wDia=wDia, 
                                  wDiaPath=wDiaPath, pageRank=pageRank, wPageRank=wPageRank, 
                                  vertBetween=vertBetween, edgeBetween=edgeBetween, maxFluxPath=maxFluxPath, 
                                  degreeCentrIn=degreeCentrIn, degreeCentrOut=degreeCentrOut,
                                  degreeCentrTot=degreeCentrTot, distIn=distIn, distOut=distOut)
    }
    else{
      # fill in results
      iSlicesStats[[idx]] <- list(flux=flux, excretion=excretion, dia=dia, diaPath=diaPath, wDia=wDia, 
                                  wDiaPath=wDiaPath, pageRank=pageRank, wPageRank=wPageRank, 
                                  vertBetween=vertBetween, edgeBetween=edgeBetween, degreeCentrIn=degreeCentrIn, 
                                  degreeCentrOut=degreeCentrOut, degreeCentrTot=degreeCentrTot)
    }
    rm(cEdges)
    rm(cGraph)
  }
  return(iSlicesStats)
}


# compare the slices BEFORE AFTER
compareSlices <- function(statsBefore, statsAfter){
  nrowsDF <- 2*length(statsBefore)
  sliceComparison <- data.frame(timepoints = numeric(nrowsDF),
                                fluxSum = numeric(nrowsDF),
                                fluxSumEx = numeric(nrowsDF),
                                fluxSumUp = numeric(nrowsDF),
                                dia = numeric(nrowsDF),
                                diaPath = I(vector("list", nrowsDF)),
                                wDia = numeric(nrowsDF),
                                wDiaPath = I(vector("list", nrowsDF)),
                                pageRank = I(vector("list", nrowsDF)),
                                wPageRank = I(vector("list", nrowsDF)),
                                vertBetween = I(vector("list", nrowsDF)),
                                edgeBetween = I(vector("list", nrowsDF)),
                                degreeCentrIn = I(vector("list", nrowsDF)),
                                degreeCentrOut = I(vector("list", nrowsDF)),
                                degreeCentrTot = I(vector("list", nrowsDF)),
                                maxFluxPath = I(vector("list", nrowsDF)),
                                dist = I(vector("list", nrowsDF)))
  rowIdx <- 1
  for (timeIdx in 1:length(statsBefore)){
    sliceComparison$timepoints[c(rowIdx, rowIdx+1)] <- names(statsBefore)[timeIdx]
    # Fill in BEFORE values
    sliceComparison$fluxSum[rowIdx] <- sum(statsBefore[[timeIdx]][["flux"]])
    sliceComparison$fluxSumEx[rowIdx] <- sum(statsBefore[[timeIdx]][["flux"]][statsBefore[[timeIdx]][["excretion"]]])
    sliceComparison$fluxSumUp[rowIdx] <- sum(statsBefore[[timeIdx]][["flux"]][!statsBefore[[timeIdx]][["excretion"]]])
    sliceComparison$dia[rowIdx] <- statsBefore[[timeIdx]][["dia"]]
    sliceComparison$diaPath[rowIdx] <- list(names(statsBefore[[timeIdx]][["diaPath"]]))
    sliceComparison$wDia[rowIdx] <- statsBefore[[timeIdx]][["wDia"]]
    sliceComparison$wDiaPath[rowIdx] <- list(names(statsBefore[[timeIdx]][["wDiaPath"]]))
    sliceComparison$pageRank[rowIdx] <- statsBefore[[timeIdx]][["pageRank"]]["vector"]
    sliceComparison$wPageRank[rowIdx] <- statsBefore[[timeIdx]][["wPageRank"]]["vector"]
    sliceComparison$vertBetween[rowIdx] <- statsBefore[[timeIdx]]["vertBetween"]
    sliceComparison$edgeBetween[rowIdx] <- statsBefore[[timeIdx]]["edgeBetween"]
    sliceComparison$degreeCentrIn[rowIdx] <- statsBefore[[timeIdx]][["degreeCentrIn"]]["res"]
    sliceComparison$degreeCentrOut[rowIdx] <- statsBefore[[timeIdx]][["degreeCentrOut"]]["res"]
    sliceComparison$degreeCentrTot[rowIdx] <- statsBefore[[timeIdx]][["degreeCentrTot"]]["res"]
    rowIdx <- rowIdx +1
    
    # Fill in After values
    sliceComparison$fluxSum[rowIdx] <- sum(statsAfter[[timeIdx]][["flux"]])
    sliceComparison$fluxSumEx[rowIdx] <- sum(statsAfter[[timeIdx]][["flux"]][statsAfter[[timeIdx]][["excretion"]]])
    sliceComparison$fluxSumUp[rowIdx] <- sum(statsAfter[[timeIdx]][["flux"]][!statsAfter[[timeIdx]][["excretion"]]])
    sliceComparison$dia[rowIdx] <- statsAfter[[timeIdx]][["dia"]]
    sliceComparison$diaPath[rowIdx] <- list(names(statsAfter[[timeIdx]][["diaPath"]]))
    sliceComparison$wDia[rowIdx] <- statsAfter[[timeIdx]][["wDia"]]
    sliceComparison$wDiaPath[rowIdx] <- list(names(statsAfter[[timeIdx]][["wDiaPath"]]))
    sliceComparison$pageRank[rowIdx] <- statsAfter[[timeIdx]][["pageRank"]]["vector"]
    sliceComparison$wPageRank[rowIdx] <- statsAfter[[timeIdx]][["wPageRank"]]["vector"]
    sliceComparison$vertBetween[rowIdx] <- statsAfter[[timeIdx]]["vertBetween"]
    sliceComparison$edgeBetween[rowIdx] <- statsAfter[[timeIdx]]["edgeBetween"]
    sliceComparison$degreeCentrIn[rowIdx] <- statsAfter[[timeIdx]][["degreeCentrIn"]]["res"]
    sliceComparison$degreeCentrOut[rowIdx] <- statsAfter[[timeIdx]][["degreeCentrOut"]]["res"]
    sliceComparison$degreeCentrTot[rowIdx] <- statsAfter[[timeIdx]][["degreeCentrTot"]]["res"]
    sliceComparison$maxFluxPath[rowIdx] <- statsAfter[[timeIdx]]["maxFluxPath"]
    sliceComparison$dist[rowIdx] <- statsAfter[[timeIdx]]["dist"] 
    rowIdx <- rowIdx +1
  }
  return(sliceComparison)
}


# extract unique union of removed nodes from linearized network slices
uniqueRmEdges <- function(iSlicesLin, from=1, to=length(iSlicesLin)){
  idxs <- c(from:to)
  unionIDs <- vector("character", 0L)
  for (idx in idxs){
    unionIDs <- c(unionIDs, iSlicesLin[[idx]][["rmEdgeIDs"]])
  }
  unique(unionIDs)
}


# count frequency of edge removal
freqRmEdges <- function(iSlicesLin, uniques = uniqueRmEdges(iSlicesLin)){
  nrEdges <- length(uniques)
  freqDF <- data.frame(matrix(NA, nrow = 1, ncol = nrEdges))
  names(freqDF) <- uniques
  frequency <- rep(0, nrEdges)
  occurrence <- vector(mode="list", length=nrEdges)
  names(occurrence) <- uniques
  for(idx in 1:length(iSlicesLin)){
    isPartOfSlice <- uniques %in% iSlicesLin[[idx]][["rmEdgeIDs"]]
    currentFreq <- sapply(isPartOfSlice, function(x){
      if(x){1}
      else{0}
    })
    frequency <- frequency + currentFreq
    currentOccurrences <- uniques[isPartOfSlice]
    for(edge in currentOccurrences){
      occurrence[[edge]] <- c(occurrence[[edge]], names(iSlicesLin)[idx])
    }
  }
  freqDF[1,] <- frequency
  return(list(freq=freqDF, occ=occurrence))
}

# 
checkEdges <- function(full, lin){
  dirT <- "C:/Users/rauj/Documents/00_Masterarbeit/1_Plots/3_5_linearization"
  t1 <- read.csv(paste(dirT, full, sep = '/'), header=T, as.is=T)
  t2 <- read.csv(paste(dirT, lin, sep = '/'), header=T, as.is=T)
  setdiff(t1$reacID,t2$reacID)
}





#slCompMVP[!slCompMVP %in% c("h2o[e]", "h[e]", "time")]
# predefine source for linear represenation of graph

# predefine sinks/Zielvariablen, zb max flux
