#' @author Tobias Frisch 

source("data_processing.R")

#######################################################
# loading the basic dataset
#######################################################

myJoinedData <- generateJoinedSet()

# setting FDR to 1 if it is NA
tmpCols = grepl("FDR", colnames(myJoinedData))
myJoinedData[,tmpCols][is.na(myJoinedData[, tmpCols])] = 1

# setting the logFC to 0 if it is NA
tmpCols = grepl("logFC", colnames(myJoinedData))
myJoinedData[,tmpCols][is.na(myJoinedData[, tmpCols])] = 0

myJoinedData$`Gene Symbol` <- ifelse(is.na(myJoinedData$`Gene Symbol`), myJoinedData$Ensemble, myJoinedData$`Gene Symbol`)
myJoinedData$`Gene Symbol` <- ifelse(myJoinedData$`Gene Symbol` == "NA", myJoinedData$Ensemble, myJoinedData$`Gene Symbol`)

global_limits = range(myJoinedData[,tmpCols], na.rm = T)
#######################################################
# solving the directory problem
#######################################################
makeTmpDir = function(uniqueID){
  newDir = paste0("/home/atlas/", uniqueID)
  
  if(dir.exists(newDir))
    stop("Directory should not exist!")
  
  dir.create(newDir, showWarnings = T);
  
  return(newDir)
  

}



#######################################################
# server implementation
#######################################################
server <- function(input, output, session) {
 
  selectedSet <- NULL
  
  selectedNetwork <- NULL
  
  clusterOrder <- NULL
  
  specifiedGenes <- NULL
  
  ##########
  # initialize and destroy tmp directory
  ##########
  sessionDir = makeTmpDir(session$token)
  
  session$onSessionEnded(function() {
    if(dir.exists(sessionDir))
      unlink(sessionDir, recursive = T)
  })
  
  
  
  
  ########
  generateSelectedSet <-function(){
    
    print("function:generateSelectedSet")
    showNotification("Calculating Dataset", type = "message")
    a = generateHeatMap(dataset = myJoinedData, al = input$AL, il = input$IL, 
                        ca = input$CA, rl = input$RL, nawm = input$NAWM,
                        al_reg = input$AL_REG, il_reg = input$IL_REG,
                        ca_reg = input$CA_REG, nawm_reg = input$NAWM_REG,
                        rl_reg = input$RL_REG, fdr = as.numeric(input$FDR),
                        logFC = as.numeric(input$logFC),
                        allFC = (input$LOGFC_ALL == "all"),
                        allFDR = ((input$FDR_ALL) == "all"),
                        n_genes = as.numeric(input$n_genes))
    
    #deal with NA and duplicated row names!
    rownames(a) = make.names(a[,2], unique = T)
    
    b = a[grepl("logFC", colnames(a))]
    b[is.na(b)] = 0
    
    
    # Checking for specific genes
    genes = trimws(strsplit(input$geneInput, ",")[[1]])
    if(length(genes) > 0){
      
      
      # exact pattern matching 
      #b = b[rownames(b) %in% genes,]
      
      # relative pattern matching
      myIndexes = lapply(genes, function(x){grep(x, rownames(b))})
      myIndexes = unlist(myIndexes, recursive = T)
      b = b[myIndexes,]
    }
 
    selectedSet <<- b
    
    print(sprintf("generateSelectedSet: (%d,%d)", nrow(selectedSet),
                  ncol(selectedSet)))
    
    # Clustering
    print("Distance and Clustering for heatmap")
    if(nrow(selectedSet) >= 2){
      myDis = dist(selectedSet)
      myClus = hclust(myDis)
      clusterOrder <<- myClus$order
    }else{
      clusterOrder <<- NULL
    }
  }
  
  generateselectedNetwork <- function(){
    print("function:generateSelectedNetwork")
    showNotification("Generating Network If you have selected KeyPathwayMiner this can take several minutes!", type = "message")
    localNetwork <- generateGraph(rownames(selectedSet), isolate(as.numeric(input$K)),
                                  isolate(as.numeric(input$C)), 
                                  sessionDir, isolate(input$enrichment == "KeyPathwayMiner"))
  
    if(is.null(localNetwork)){
      selectedNetwork <<- NULL
    }else{
      if(vcount(localNetwork)> 10){
        fine = 500
        pal = colorRampPalette(c("blue", "red"))
        V(localNetwork)$color = 
          pal(fine)[as.numeric(cut(betweenness(localNetwork), breaks = fine))]
        
      }
      
      if(vcount(localNetwork) <= 10){
        V(localNetwork)$size = 5
      }
      
      if(vcount(localNetwork) <= 0){
        showNotification("The network you have selected is empty. There are two possible reasons: (1) The genes that you selected are not part of our network or (2) we where 
                         not able to find connected components within the network. Please consider to change parameters for the network enrichment process to solve the problem.", type = "warning")
      }
      selectedNetwork <<- localNetwork
    }
  }
  
  readFile <- function(fileName){
    df = read.csv(fileName, header = F, sep = ";", quote = "", stringsAsFactors = F)
    df[1,]
  }
  
  observeEvent(input$applyToNetwork, ignoreInit = T, ignoreNULL = T, {
    print("observeEvent input applyToNetwork")
    
    output$network <- renderVisNetwork({
      generateselectedNetwork()
      print("output$network")
      
      print(sprintf("Number of nodes: %d", vcount(selectedNetwork)))
      
      # layout the graph
      l = layout_with_fr(selectedNetwork, dim = 3)
      
      #graphjs object
      #graphjs(g = selectedNetwork, layout = l, showLabels = T)
      data = toVisNetworkData(selectedNetwork)
      visNetwork(nodes = data$nodes, edges = data$edges)
    })
  })
  
  output$volcano <- renderPlotly({
    print("output$volcano")
    
    vs = input$volcanoSlection;
    
    localRes = generateVolcano(myJoinedData, vs == "active", vs == "inactive",
                               vs == "chronic active", vs == "NAWM",
                               vs == "remyelinating")
    
    if(is.null(localRes)){
      showNotification("At the moment you can only select one lesion type for",
                       duration = 10, type = "error")
      return(NULL)
    }
    
    #point labels
    GeneSymbol = localRes$GeneSymbol
    GeneSymbol[is.na(GeneSymbol)] = "Unknown"
    GeneSymbol = make.unique(GeneSymbol)
    
    #point colors
    myColors = rep("brown", nrow(localRes))
    myFDR= localRes$FDR
    myFC = abs(localRes$logFC)
    myColors[which(myFDR < 0.05)] = "red"
    myColors[which(myFC > 1)] = "grey"
    myColors[which(myFC > 1 & myFDR < 0.05)] = "orange"
    #print(table(myColors))

    p <- ggplot(data = localRes, aes(x = logFC, y = -log10(PValue),
                                     label=GeneSymbol)) +
      geom_point(colour = myColors) +
      xlab("log2 fold change") +
      ylab("-log10 p-value")
    return(ggplotly(p))
  })
  
  observeEvent(input$geneList,{
    print("observing: input$geneList")
    
    print(sprintf("Trying to read file: %s", input$geneList$datapath))
    
    specifiedGenes <<- readFile(input$geneList$datapath)
    specifiedGenes <<- as.character(specifiedGenes)
    
    specifiedGenes <<- unlist(strsplit(specifiedGenes, split = ","))
    
    a = generateHeatMap(dataset = myJoinedData, al = input$AL, il = input$IL, 
                        ca = input$CA, rl = input$RL, nawm = input$NAWM,
                        al_reg = input$AL_REG, il_reg = input$IL_REG,
                        ca_reg = input$CA_REG, nawm_reg = input$NAWM_REG,
                        rl_reg = input$RL_REG, fdr = as.numeric(input$FDR),
                        logFC = as.numeric(input$logFC),
                        allFC = (input$LOGFC_ALL == "all"),
                        allFDR = (input$FDR_ALL) == "all")
    
    mySubset = a[,2] %in% specifiedGenes
    
    #print(a)
    #print(mySubset)
    
    a = a[mySubset,]
    
    #remove duplicated row names
    rownames(a) = make.unique(a[,2])
    
    b = a[grepl("logFC", colnames(a))]
    b[is.na(b)] = 0
    
    selectedSet <<- b
    
    # Clustering
    print("Distance and Clustering for heatmap")
    if(nrow(selectedSet) >= 2){
      myDis = dist(selectedSet)
      myClus = hclust(myDis)
      clusterOrder <<- myClus$order
    }else{
      clusterOrder <<- NULL
    }
    
    myLimits=NULL
    if(isolate(input$adjust_colors)){
      myLimits = rep(max(abs(range(selectedSet, na.rm = T))), 2)
      minValue = as.numeric(isolate(input$range_min))
      maxValue = as.numeric(isolate(input$range_max))
      if(!is.na(minValue))
        myLimits[1] = minValue
      
      if(!is.na(maxValue))
        myLimits[2] = maxValue
    }
    
    generateselectedNetwork()
    
    output$heatMap <- renderPlotly({
      print("function: output$heatMap 2")
      
      print(selectedSet)
      
      if(nrow(selectedSet) == 0){
        showNotification("There are no genes that fullfill your requirements.", type = "warning")
        return(NULL)
      }
      
      if(nrow(selectedSet) < 2 && ncol(selectedSet) < 2){
        print("1")
        heatmaply(selectedSet, Rowv = F, Colv = F,
                  c("blue", "grey", "red"), 
                  labRow = rownames(selectedSet), 
                  labCol = colnames(selectedSet), limits = myLimits)
      }else if(nrow(selectedSet) < 2){
        print("2")
        heatmaply(selectedSet, Rowv = F, Colv = T,
                  c("blue", "grey", "red"), 
                  labRow = rownames(selectedSet), 
                  labCol = colnames(selectedSet), limits = myLimits)
      }else if(ncol(selectedSet) < 2){
        print("3")
        heatmaply(as.data.frame(selectedSet[clusterOrder,]), Rowv = isolate(input$row_dendogram), Colv = F,
                  c("blue", "grey", "red"), 
                  labRow = rownames(selectedSet)[clusterOrder], 
                  labCol = colnames(selectedSet), limits = myLimits)
        
      }else{
        print("4")
        heatmaply(selectedSet[clusterOrder,], Rowv = isolate(input$row_dendogram), Colv = T,
                  c("blue", "grey", "red"), 
                  labRow = rownames(selectedSet[clusterOrder,]), 
                  labCol = colnames(selectedSet), limits = myLimits)
      }
    })
    
    output$network <- renderVisNetwork({
      
      print("output$network")
      
      print(sprintf("Number of nodes: %d", vcount(selectedNetwork)))
      
      # layout the graph
      l = layout_with_fr(selectedNetwork, dim = 3)
      
      #graphjs object
      #graphjs(g = selectedNetwork, layout = l, showLabels = T)
      data = toVisNetworkData(selectedNetwork)
      visNetwork(nodes = data$nodes, edges = data$edges)
    })
  })
  
 
  
  observeEvent(input$start, ignoreNULL = F, ignoreInit = F, {
    print("observeEvent input$start")
    
    updateSelectInput(session, "C", selected = "1")
    updateSelectInput(session, "K", selected = "0")
    
    print("Changing input:")
    print(input$C)
    
    generateSelectedSet()
    if(nrow(selectedSet) == 0){
      showNotification("There are no genes that fullfill your requirements.", type = "warning")
      return(NULL)
    }
    generateselectedNetwork()
    
    output$heatMap <- renderPlotly({
      print("function: output$heatMap 2")
      
      if(nrow(selectedSet) == 0){
        showNotification("There are no genes that fullfill your requirements.", type = "warning")
        return(NULL)
      }
      
      myLimits=NULL
      if(isolate(input$adjust_colors)){
        myLimits = rep(max(abs(range(selectedSet, na.rm = T))), 2)
        myLimits[1] = -myLimits[1]
        minValue = as.numeric(isolate(input$range_min))
        maxValue = as.numeric(isolate(input$range_max))
        if(!is.na(minValue))
          myLimits[1] = minValue
        
        if(!is.na(maxValue))
          myLimits[2] = maxValue
      }
      
      if(nrow(selectedSet) < 2 && ncol(selectedSet) < 2){
        print("1")
        heatmaply(selectedSet, Rowv = F, Colv = F,
                  c("blue", "grey", "red"), 
                  labRow = rownames(selectedSet), 
                  labCol = colnames(selectedSet), limits = myLimits)
      }else if(nrow(selectedSet) < 2){
        print("2")
        heatmaply(selectedSet, Rowv = F, Colv = T,
                  c("blue", "grey", "red"), 
                  labRow = rownames(selectedSet), 
                  labCol = colnames(selectedSet), limits = myLimits)
      }else if(ncol(selectedSet) < 2){
        print("3")
        heatmaply(as.data.frame(selectedSet[clusterOrder,]), Rowv = isolate(input$row_dendogram), Colv = F,
                  c("blue", "grey", "red"), 
                  labRow = rownames(selectedSet)[clusterOrder], 
                  labCol = colnames(selectedSet), limits = myLimits)
        
      }else{
        print("4")
        heatmaply(selectedSet[clusterOrder,], Rowv = isolate(input$row_dendogram), Colv = T,
                  c("blue", "grey", "red"), 
                  labRow = rownames(selectedSet[clusterOrder,]), 
                  labCol = colnames(selectedSet), limits = myLimits)
      }
    })
    
    output$network <- renderVisNetwork({
      
      print("output$network")
      
      print(sprintf("Number of nodes: %d", vcount(selectedNetwork)))
      
      # layout the graph
      l = layout_with_fr(selectedNetwork, dim = 3)
      
      #graphjs object
      #graphjs(g = selectedNetwork, layout = l, showLabels = T)
      data = toVisNetworkData(selectedNetwork)
      visNetwork(nodes = data$nodes, edges = data$edges)
    })
  })
  
  output$export_genes <- downloadHandler(
    filename = function(){
      "geneList.csv"
    },
    
    content = function(file){
      write.csv2(selectedSet, file, quote = F, append = F)
    }
  )
  
  output$export_network <- downloadHandler(
    filename = function(){
      "network.csv"
    },
    
    content = function(file){
      myGraph = as_edgelist(selectedNetwork, names = T)
      
      myGraph = myGraph[myGraph[,1] != myGraph[,2],]
      
      write.table(myGraph, file, append = F, quote = F, row.names = F,
                  col.names = F, sep = "\t")
    }
  )
}