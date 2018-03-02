load("DarmanisHumanExp.rda")
DarmanisExp <- DarmanisHumanExp %>% data.frame()
DarmanisExp$GeneSymbol <- rownames(DarmanisExp)

load("DarmanisHumanMeta.rda")
DarmanisMeta <- DarmanisHumanMeta

GetHumanExp <- function(genes, markerInfo = NULL, corInfo = NULL, CellType = NULL,
                        txtSize = 16, ptSize = 0.4, title = NULL,
                        colors = c("darkorange", "firebrick1", "darkorchid1",
                                   "darkolivegreen4", "darkseagreen", "dodgerblue3")){
  genes = toupper(genes)
  #this part is here because of Gemma annotations
  genes <- sapply(genes, function(x) strsplit(x, "\\.")[[1]][1])
  tempGene <- DarmanisExp %>% filter(GeneSymbol %in% genes)
  if(is.null(markerInfo)){
    marker = ""
  } else {
    marker <- sapply(tempGene$GeneSymbol, function(gene){
      geneOrg = names(genes)[genes %in% gene]
      if(is.na(markerInfo %>% filter(GeneSymbol == geneOrg) %>% .$GeneType)){
        ""
      } else if(markerInfo %>% filter(GeneSymbol == geneOrg) %>% .$GeneType == CellType){
        "*"
      } else {
        "#"
      }
    })
  }
  if(is.null(corInfo)){
    Cor = ""
  } else {
    Cor <- sapply(tempGene$GeneSymbol, function(gene){
      geneOrg = names(genes)[genes %in% gene]
      corInfo %>% filter(GeneSymbol == geneOrg) %>% .[2] %>% signif(digits = 2)
      })
  }
  tempGene %<>% mutate(GeneSymbol = factor(GeneSymbol, levels = genes))
  tempGene %<>% droplevels()
  tempGene %<>% mutate(GeneName = paste0(tempGene$GeneSymbol, marker, " (r =", Cor, ")"))
  tempGene$GeneName <- factor(tempGene$GeneName,
                              levels = tempGene$GeneName[match(levels(tempGene$GeneSymbol),
                                                               tempGene$GeneSymbol)])
  tempGeneMelt <- melt(tempGene, id.vars = c("GeneSymbol", "GeneName"), variable.name = "GSM", value.name = "Expression")
  tempGeneMelt$CellType <- DarmanisMeta$cellType[match(tempGeneMelt$GSM, DarmanisMeta$GSM)]
  tempGeneMelt %<>% filter(!CellType %in% c("hybrid", "fetal_quiescent","fetal_replicating"))
  tempGeneMelt$CellType <- factor(tempGeneMelt$CellType, levels = c("astrocytes", "endothelial",
                                                                    "microglia", "oligodendrocytes",
                                                                    "OPC", "neurons"))
  tempGeneMelt %<>% mutate(Expression = log2(Expression + 1))
  if(is.null(title)){
    title = CellType
  }
  p <- ggplot(tempGeneMelt, aes(CellType, Expression))
  plot <- p + labs(title = title, x = "", y = "log2(RPKM +1) expression") +
    theme_bw(base_size = txtSize) +
    theme(axis.text.y = element_text(size = rel(1.1)),
          axis.text.x = element_text(size = rel(1.1), angle = 50, hjust = 1),
          legend.position = "none",
          panel.grid = element_blank()) +
    geom_boxplot(outlier.color = NA, width=0.8) + 
    geom_jitter(width = 0.2, size = ptSize, aes(color = CellType)) +
    scale_color_manual(values = colors, name = "") +
    facet_wrap(~GeneName, scales = "free_y")
  return(plot)
}


MouseGenesProp <- function(genes, MGPused, ExpData = DarmanisExp, MetaData = DarmanisMeta,
                           ExcludeCell = c("hybrid", "fetal_quiescent","fetal_replicating"),
                           cellOrder = c("astrocytes", "endothelial",
                                         "microglia", "oligodendrocytes",
                                         "OPC", "neurons"),
                           cellColor = c("red", "orange", "grey", "green", "darkgreen", "blue")){
  genes = toupper(genes)
  #this part is here because of Gemma annotations
  genes <- sapply(genes, function(x) strsplit(x, "\\.")[[1]][1])
  tempGene <- ExpData[ExpData$GeneSymbol %in% genes,]
  tempGene %<>% mutate(GeneSymbol = factor(GeneSymbol, levels = genes)) %>% droplevels()
  tempGeneMelt <- melt(tempGene, id.vars = c("GeneSymbol"), variable.name = "GSM", value.name = "Expression")
  tempGeneMelt$CellType <- MetaData$cellType[match(tempGeneMelt$GSM, MetaData$GSM)]
  tempGeneMelt %<>% filter(!CellType %in% ExcludeCell)

  tempGeneMelt %<>% mutate(Expression = log2(Expression + 1))
  tempGeneMelt$IsExpressed <- "NO"
  tempGeneMelt$IsExpressed[which(tempGeneMelt$Expression != 0)] <- "YES"

  
  CellNumber <- group_by(MetaData %>% filter(!cellType %in% ExcludeCell),
                         cellType) %>% summarise(n = n()) %>% data.frame()
  
  CellExpr <- sapply(CellNumber$cellType, function(cell){
    NumCell = CellNumber %>% filter(cellType == cell) %>% .$n
    sapply(unique(tempGeneMelt$GeneSymbol), function(gene){
      PropExp <- nrow(tempGeneMelt %>% filter(GeneSymbol == gene,
                                               CellType == cell,
                                               IsExpressed == "YES"))/NumCell
      MedianExp <- median(tempGeneMelt %>% filter(GeneSymbol == gene,
                                                  CellType == cell) %>% .$Expression)
      data.frame(CellType = cell,
                 GeneSymbol = gene,
                 PropExp = PropExp,
                 MedianExp = MedianExp)
    }, simplify = FALSE) %>% rbindlist()
  }, simplify = FALSE) %>% rbindlist() %>% data.frame
  
  CellExpr$MGPused <- "NO"
  CellExpr$MGPused[which(CellExpr$GeneSymbol %in% MGPused)] <- "YES"
  CellExpr$MGPused <- factor(CellExpr$MGPused, levels = c("YES", "NO"))
  CellExpr$CellType <- factor(CellExpr$CellType, levels = cellOrder)
  
  p <- ggplot(CellExpr, aes(CellType, PropExp)) +
    theme(axis.text.x = element_blank(),
          legend.position = "none") +
    geom_boxplot(outlier.shape = NA, aes(fill = CellType), alpha = 0.8) +
    scale_fill_manual(values = cellColor) +
    geom_jitter(size = 0.5, alpha = 0.2, width = 0.2) + facet_wrap(~MGPused)
  
  p2 <- ggplot(CellExpr, aes(CellType, MedianExp)) +
    theme(axis.text.x = element_blank(),
          legend.position = "none") +
    labs(y = "Median expression (log2(fpkm + 1))") + 
    geom_boxplot(outlier.shape = NA, aes(fill = CellType), alpha = 0.8) +
    scale_fill_manual(values = cellColor) +
    geom_jitter(size = 0.5, alpha = 0.2, width = 0.2) + facet_wrap(~MGPused)
  
  #Normalize the expression 0-1 for heatmap
  ScaledExp <- sapply(tempGeneMelt$GeneSymbol, function(gene){
    subData <- tempGeneMelt[tempGeneMelt$GeneSymbol == gene,]
    subData$NormExp <- rescale(subData$Expression, c(0,1))
    subData
  }, simplify = FALSE) %>% rbindlist %<>% arrange(CellType)
  
  ScaledExp$GSM <- factor(ScaledExp$GSM, levels = unique(ScaledExp$GSM))
  
  ScaledExpCellType <- group_by(ScaledExp, CellType, GeneSymbol) %>%
    summarise(MeanExp = mean(NormExp), MedianExp = median(NormExp)) %>% data.frame() %>%
    arrange(CellType, MedianExp, MeanExp) %>% droplevels()
  
  GeneOrder <- rev(unique(ScaledExpCellType$GeneSymbol))
  ScaledExp$GeneSymbol <- factor(ScaledExp$GeneSymbol,
                                 levels = GeneOrder)
  
  ScaledExp$MGPused <- "NO"
  ScaledExp$MGPused[which(ScaledExp$GeneSymbol %in% MGPused)] <- "YES"
  
  ExpMatrix <- tempGene %>% select(-GeneSymbol) %>% as.matrix
  rownames(ExpMatrix) <- tempGene$GeneSymbol
  ExpMatrix <- apply(ExpMatrix, c(1,2), function(x) log2(x+1))
  ExpMatrixNorm <- apply(ExpMatrix, 1, function(x) rescale(x, c(0,1))) %>% t %>% data.frame
  ExpMatrixNorm %<>% select_(.dots = levels(ScaledExp$GSM))
  ExpMatrixNorm %<>% as.matrix
  ExpMatrixNorm <- ExpMatrixNorm[match(GeneOrder, rownames(ExpMatrixNorm)),]
  
  ann_column = data.frame(CellType = factor(MetaData$cellType[match(colnames(ExpMatrixNorm), DarmanisMeta$GSM)]))
  rownames(ann_column) <- colnames(ExpMatrixNorm)
  ann_colors = list(
    CellType = cellColor)
  names(ann_colors$CellType) <- cellOrder
  PheatData <- list(Mat = ExpMatrixNorm,
                    ann_column = ann_column,
                    ann_colors = ann_colors)
  
  return(list(AllCellExpr = tempGeneMelt,
              CellTypeExr = CellExpr,
              Boxplot1 = p,
              Boxplot2 = p2,
              PheatData = PheatData))
}

