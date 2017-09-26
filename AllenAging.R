AllenAgingExp <- read.table("DevelopmentData/AllenAging/fpkm_table_normalized.csv", header = T, sep = ",")
names(AllenAgingExp)[1] <- "gene_id"
names(AllenAgingExp) <- sapply(names(AllenAgingExp), function(x) gsub("X", "", x))
AllenAgingGenes <- read.table("DevelopmentData/AllenAging/rows-genes.csv", header = T, sep = ",")
AllenAgingExp <- merge(AllenAgingGenes %>% select(-chromosome, -gene_entrez_id), AllenAgingExp, by="gene_id")
names(AllenAgingExp)[1:3] <- c("GeneID", "GeneSymbol", "GeneName")
AllenAgingSamples <- read.table("DevelopmentData/AllenAging/columns-samples.csv", header = T, sep = ",")
AllenAgingMeta <- read.table("DevelopmentData/AllenAging/DonorInformation2.csv", header = T, sep = ";", na.strings = "N/A")


GabaPVex <- AllenAgingExpOrg %>% filter(GeneSymbol %in% markerGenesHumanUpdate$GabaPV)
GabaPVexpMean <- sapply(GabaPVex$GeneSymbol, function(x){
  GabaPVex %>% filter(GeneSymbol == x) %>% .[-c(1:3)] %>% unlist %>% mean
}) %>% data.frame(GeneSymbol = GabaPVex$GeneSymbol, MeanExp = .) %>% arrange(MeanExp)

GabaPVexpMean$PVALBcor <- sapply(GabaPVexpMean$GeneSymbol, function(x){
  geneExp <- GabaPVex %>% filter(GeneSymbol == x) %>% .[-c(1:3)] %>% unlist
  pvalbExp <- GabaPVex %>% filter(GeneSymbol == "PVALB") %>% .[-c(1:3)] %>% unlist
  cor(geneExp, pvalbExp, method = "spearman")
})

GabaPVexpMean$OldMark <- sapply(GabaPVexpMean$GeneSymbol, function(x){
  if(x %in% markerGenesHuman$GabaPV){
    "Yes"
  } else {
    "No"
  }
}) %>% factor

GabaPVexpMean$PVALBcor <- sapply(GabaPVexpMean$GeneSymbol, function(x){
  geneExp <- GabaPVex %>% filter(GeneSymbol == x) %>% .[-c(1:3)] %>% unlist
  pvalbExp <- GabaPVex %>% filter(GeneSymbol == "PVALB") %>% .[-c(1:3)] %>% unlist
  cor(geneExp, pvalbExp, method = "spearman")
})

GabaPVexpMean$OldMark <- sapply(GabaPVexpMean$GeneSymbol, function(x){
  if(x %in% markerGenesHuman$GabaPV){
    "Yes"
  } else {
    "No"
  }
}) %>% factor
