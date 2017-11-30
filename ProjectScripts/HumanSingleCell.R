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
  tempGene %<>% mutate(GeneName = paste0(tempGene$GeneSymbol, marker, " (rho=", Cor, ")"))
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

