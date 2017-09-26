PCA_genes_All_based <- function(dataset_id, dataset, CellType_genes){ 
  print("########################################")
  print("##### Scales based on all samples ####")
  print("########################################")
  
  
  dataset <- dataset[,sapply(dataset[1,], is.numeric)] # Exclude the annotation columns
  groups <- sapply(grep("_", colnames(dataset), value=T),
                   function(x) gsub("_.*", "", x)) %>% table %>% names
  for(grp in groups){
    assign(grp, grep(grp, colnames(dataset)))
  }
  
  PCAresults <- list()
  PCAresults$dataset_id <- dataset_id
  PCAresults$ControlOnly <- list()
  PCAresults$All <- list()
  PCAresults$modified <- list()
  
  for(i in 1:c(length(names(CellType_genes)))){
    print(names(CellType_genes)[i])
    if(length(which(rownames(dataset) %in% CellType_genes[[i]])) > 2){
      data <- dataset %>% mutate(GeneSymbol = rownames(.)) %>%
        filter(GeneSymbol %in% CellType_genes[[i]])
      PCAresults$ControlOnly[[i]] <- data %>%
        select(Cont) %>% t %>% prcomp(scale=T) #The assumption here is that there is a control group
      rownames(PCAresults$ControlOnly[[i]]$rotation) <- data$GeneSymbol
      PCAresults$ControlOnly[[i]] <- correct_sign(names(CellType_genes)[i], PCAresults$ControlOnly[[i]])
      PCAresults$All[[i]] <- data %>% select(-GeneSymbol) %>% t %>% prcomp(scale=TRUE) 
      PCAresults$All[[i]] <- correct_sign(names(CellType_genes)[i], PCAresults$All[[i]])
      rownames(PCAresults$All[[i]]$rotation) <- data$GeneSymbol
      while (sum(PCAresults$All[[i]]$rotation[,1] > 0) < nrow(PCAresults$All[[i]]$rotation)){
        minorGene <- rownames(PCAresults$All[[i]]$rotation)[PCAresults$All[[i]]$rotation[,1] < 0]
        data %<>% filter(!GeneSymbol %in% minorGene)
        rownames(data) <- data$GeneSymbol
        PCAresults$All[[i]] <- data %>% select(-GeneSymbol) %>% t %>% prcomp(scale=TRUE) 
        PCAresults$All[[i]] <- correct_sign(names(CellType_genes)[i], PCAresults$All[[i]])
      }
      
      PCAresults$modified[[i]] <- PCAresults$All[[i]]$x %>% list
      PCAresults$modified[[i]] <- apply(PCAresults$modified[[i]][[1]],2,
                                               function(x) rescale(x,c(0,1))) %>% list
      names(PCAresults$modified[[i]]) <- "x"
      
    } else {
      print(paste("Looki here!!! No genes for", names(CellType_genes)[i]))
      x <- matrix(nrow=ncol(dataset), ncol=1)
      rownames(x) <- names(dataset)[sapply(groups, function(grp){
        eval(as.name(grp))
      }) %>% unlist] ;colnames(x)="x"
      PCAresults$ControlOnly[[i]] <- list(x)
      names(PCAresults$ControlOnly[[i]]) <- "x"
      PCAresults$All[[i]] <- list(x)
      names(PCAresults$All[[i]]) <- "x"
      PCAresults$modified[[i]] <- list(x)
      names(PCAresults$modified[[i]]) <- "x"
    }
    
  }
  for(i in 2:4){
    names(PCAresults[[i]]) <- names(CellType_genes)
  }
  return(PCAresults)
}
