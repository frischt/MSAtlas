#######################################################
# load background data
#######################################################

#' WM_IL
#' WM_NAWM
#' WM_AL
#' WM_CA
#' WM_RL
load(file =  "data/WM_LESION.RData")

#######################################################
# network files
#######################################################
mapping_file = read.table("data/ntrez_to_gene_symbol_bioc.txt", header = T,
                          stringsAsFactors = F)

network_file = read.table("data/I2D_network.sif", stringsAsFactors = F,
                          header = F, sep = "\t", comment.char = "")

possibleID = unique(c(WM_AL$other_ID, WM_CA$other_ID, WM_IL$other_ID, 
                      WM_NAWM$other_ID, WM_RL$other_ID))

myID = unique(c(network_file$V1, network_file$V2))

myNodes = data.frame(
  name = myID,
  label= myID,
  color = rep("white", length(myID)),
  size = rep(3, length(myID)), 
  stringsAsFactors = F
)

myEdges = data.frame(
  from = network_file$V1,
  to = network_file$V2, 
  stringsAsFactors = F
)

generateGraph <- function(namesList, k, c, tmpDir, keyPathwayMiner){
  print("function: generateGraph")
  print(sprintf("%d genes are selected and the graph contains %d nodes.", length(namesList), length(myID)))
  nameInID = namesList %in% myID
  namesList = namesList[nameInID]
  
  if(keyPathwayMiner){
    print("generateGraph: keyPathwayMiner")
    subnet = myID %in% runKeyPathwayMiner(namesList, k, tmpDir)
  }else{
    print("generateGraph: bigest connected component")
    subnet = myID %in% namesList
  }
  
  print(sprintf("There are %d nodes!", sum(subnet)))
  if(sum(subnet) == 0){
    return(NULL)
  }
  myNodes$color[subnet] = "black"
  myNodes$color[!subnet] = "white"
  
  #build graph object
  myGraph = graph_from_data_frame(myEdges, directed = F, vertices = myNodes)
  
  #subset graph
  myGraph = induced_subgraph(myGraph, V(myGraph)[color == "black"])

  if(k == 0){
    print("calculating connected component")
    #get the biggest connected component
    myComp = igraph::decompose(myGraph, mode = "weak")
    
    a = unlist(lapply(myComp, gorder))
    
    myOrder = order(a)
    
    print(sprintf("Selected subgraph: %d", max(length(myOrder)-c,1)))
    
    if(length(myOrder) == 0){
      return(NULL)
    }
    
    finGraph = myComp[[(myOrder[max(length(myOrder)-c,1)])]]
    
    return(finGraph)
  } else{
    return(myGraph)
  }
}

#######################################################
# Server background functions
#######################################################

#'This function uses the default input data and joines them
#'based on the ensemble ID. Be aware that it is only selecting
#'the unique ensemble_ID's. Multiple ones might be lost!
generateJoinedSet <- function(){
  print("function: generateJoinedSet")
  
  # figure all possible rownames
  joined_dataSet = unique(c(WM_AL$ensemble_ID, WM_CA$ensemble_ID, WM_IL$ensemble_ID, WM_RL$ensemble_ID, WM_NAWM$ensemble_ID))
  
  # build a dataset for joing
  joined_dataSet = data.frame(ensemble_ID = joined_dataSet, stringsAsFactors = F)
  
  # joining all datasets
  joined_dataSet = join_all(list(joined_dataSet, WM_AL, WM_IL, WM_CA, WM_NAWM, WM_RL), by = "ensemble_ID")
  
  # TODO before subset check if all other_ID are the same!
  
  # subsetting only columns of interest
  myColnames = colnames(joined_dataSet)
  myCols = which(myColnames == "logFC" | myColnames == "FDR" 
                 | myColnames =="PValue")
  joined_dataSet = joined_dataSet[,c(1,2,myCols)]
  
  # renaming the columns
  myNames = c("Ensemble", "Gene Symbol", 
              "AL logFC", "AL PValue", "AL FDR",
              "IL logFC", "IL PValue", "IL FDR",
              "CA logFC", "CA PValue", "CA FDR",
              "NAWM logFC", "NAWM PValue", "NAWM FDR",
              "RL logFC", "RL PValue", "RL FDR")
  colnames(joined_dataSet) = myNames
  
  return(joined_dataSet)
  
}


#'This function is able to decide which columns should be removed based on which
#'lesion types have been selected.
#'@return The dataset only containing columns for specific lesions
subsetData <- function(dataset, al, il, ca, nawm, rl, include_pValue = F, logFC = 0, 
                       fdr = 1000, allFDR, allFC)
{
  print("function: subsetData")
  # 1) select the right columns
  cols       = c(T, T, al, al, al, il, il, il, ca, ca, ca, nawm, nawm, nawm, rl, rl, rl)
  dataset = dataset[,cols]
  myColnames = colnames(dataset)
  
  # 2) Select columns and rows based on the fdr criteria
  fdr_cols = grepl("FDR", myColnames)
  
  if(sum(fdr_cols) > 1){
    if(allFDR){
      fdr_rows = apply(dataset[,fdr_cols], 1, function(x){all(x<fdr, na.rm = T)})
    }else{
      fdr_rows = apply(dataset[,fdr_cols], 1, function(x){any(x<fdr, na.rm = T)})
    }
  }else{
    fdr_rows  = dataset[,fdr_cols] < fdr
  }
  
  fdr_cols[c(1,2)] = T
  
  # 3) Select columns and rows based on the logFC criteria
  logFC_cols = grepl("logFC", myColnames)
  if(sum(logFC_cols) > 1){
    if(allFC){
      logFC_rows = apply(dataset[,logFC_cols], 1, function(x){all(abs(x) > abs(logFC))})
    }else{
      logFC_rows = apply(dataset[,logFC_cols], 1, function(x){any(abs(x) > abs(logFC))})
    }
  }else{
    logFC_rows = abs(dataset[,logFC_cols]) > abs(logFC)
  }
  
  logFC_cols[c(1,2)] = T
  
  
  if(include_pValue){
    pVal_cols = grepl("PValue", myColnames)
    pVal_cols[c(1,2)] = T
  }else{
    pVal_cols = c(T, T, rep(F, length(myColnames) - 2))
  }
  
  # 4) return selected columns and rows
  return(dataset[logFC_rows & fdr_rows, fdr_cols | logFC_cols | pVal_cols])
}


generateHeatMap <- function(dataset,
                            al = T,
                            il = T,
                            ca = T,
                            nawm = T,
                            rl = T,
                            al_reg = "BOTH",
                            il_reg = "BOTH",
                            ca_reg = "BOTH",
                            nawm_reg = "BOTH",
                            rl_reg = "BOTH",
                            fdr = 0.05,
                            logFC = 0,
                            allFC = F,
                            allFDR = F,
                            n_genes = 1/0){
  print("function: generateHeatMap")
  
  newSubset = subsetData(dataset = dataset, al = al, il = il, ca = ca,
                         nawm = nawm, rl = rl, logFC = logFC, fdr = fdr,
                         allFDR = allFDR, allFC = allFC)
  
  myColumns = colnames(newSubset)
  tmp = grepl("logFC", myColumns)
  
  # remove all rows that consist either only of NA or of 0
  if(sum(tmp) > 1){
    newSubset = newSubset[apply(newSubset[,grepl("logFC", myColumns)], 1,
                                function(x){any(x != 0, na.rm = T)}),]
  }else{
    newSubset = newSubset[which(newSubset[,tmp] != 0),]
  }
  
  print(sprintf("Generated a subset with dimensions: (%d, %d)", nrow(newSubset),
                ncol(newSubset)))
  
  rows = rep(T, nrow(newSubset))
  
  
  if(al && al_reg != "BOTH"){
    if(al_reg == "UP"){
      rows = rows & (newSubset[, which(myColumns == "AL logFC")] > 0)
    }else if(al_reg == "DOWN"){
      rows = rows & (newSubset[, which(myColumns == "AL logFC")] < 0)
    }
  }
  
  if(il && il_reg != "BOTH"){
    if(il_reg == "UP"){
      rows = rows & (newSubset[,which(myColumns == "IL logFC")] > 0)
    }else if(il_reg == "DOWN"){
      rows = rows & (newSubset[,which(myColumns == "IL logFC")] < 0)
    }
  }
  
  if(ca && ca_reg != "BOTH"){
    if(ca_reg == "UP"){
      rows = rows & (newSubset[,which(myColumns == "CA logFC")] > 0)
    }else if(ca_reg == "DOWN"){
      rows = rows & (newSubset[,which(myColumns == "CA logFC")] < 0)
    }
  }
  
  if(nawm && nawm_reg != "BOTH"){
    if(nawm_reg == "UP"){
      rows = rows & (newSubset[,which(myColumns == "NAWM logFC")] > 0)
    }else if(nawm_reg == "DOWN"){
      rows = rows & (newSubset[,which(myColumns == "NAWM logFC")] < 0)
    }
  }
  
  if(rl && rl_reg != "BOTH"){
    if(rl_reg == "UP"){
      rows = rows & (newSubset[,which(myColumns == "RL logFC")] > 0)
    }else if(rl_reg == "DOWN"){
      rows = rows & (newSubset[,which(myColumns == "RL logFC")] < 0)
    }
  }
  
  newSubset = newSubset[rows,]
  
  if(is.finite(n_genes) && n_genes < sum(rows)){
    maxLogTwo = apply(newSubset, 1, function(x){max(x)})
    newSubset = newSubset[order(maxLogTwo, decreasing = T)[1:n_genes],]
  }
  
  return(newSubset)
}


generateVolcano <- function(dataset,
                            al = T,
                            il = T,
                            ca = T,
                            nawm = T,
                            rl = T){
  print("function: generateVulcano")
  
  newSubset = subsetData(dataset = dataset, al = al, il = il, ca = ca,
                         nawm = nawm, rl = rl, include_pValue = T, allFDR = F, allFC = F)
  if(sum(al,il,ca,nawm,rl) != 1){
    return(NULL)
  }
  
  colnames(newSubset) = c("Ensemble", "GeneSymbol", "logFC", "FDR", "PValue")
  return(newSubset)
  
  
}


symbol_to_entrez <- function(symbolInput){
  symbolInput = symbolInput[!duplicated(symbolInput[,2]),]
  symbolInput = symbolInput[(symbolInput[,2] %in% mapping_file$SYMBOL),]
  
  
  res = as.character(mapvalues(symbolInput[,2], from = mapping_file$SYMBOL,
                               to = mapping_file$ENTREZID))
  
  return(list(res, symbolInput))
}


runKeyPathwayMiner <- function(namesVector, myK, tmpDir){
  print("runKeyPathwayMiner")
  
  unlink(paste(tmpDir, "/results", sep = ""), recursive = T)
  
  dataset = data.frame(namesVector, rep(1,length(namesVector)))
  
  #write files necessary for running
  write.table(file = paste(tmpDir, "selected_network.txt", sep = "/"),dataset, sep = "\t",
              row.names = F, col.names = F, append = F, quote = F)
  
  myCommand = paste("java -Xmx10G -jar KPM-4.0.jar -numProc=10", 
                    " -graphFile=I2D_network_interaction.sif", 
                    " -matrix1=selected_network.txt -L1=0 -K=", myK, 
                    sep = "")
  
  system(paste("cp KPM-4.0.jar kpm.properties data/I2D_network_interaction.sif",
               tmpDir))
  
  print(getwd())
  lastWorkingDIR = getwd()
  setwd(tmpDir)
  print(myCommand)
  javaOutput = system(myCommand)
  setwd(lastWorkingDIR)
  print(getwd())
  
  res = read.table(paste(tmpDir, "results/Pathway-01-NODES-.txt", sep = "/"),
                   stringsAsFactors = F)
  print("doneKeyPathwayMiner")
  res$V1
}
