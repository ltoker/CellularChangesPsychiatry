# for(study in list.files(pattern="studyFinal.*rda")){
#   load(study)
#   if(grepl("StanleyArray", study)){
#     for(std in names(studyFinal)){
#       assign(paste0("Cortex", std), studyFinal[[std]])
#     }
#   } else {
#     assign(paste0("Cortex", gsub("studyFinal|\\.rda", "", study)), studyFinal[["Cortex"]])
#   }
# }

#rm(studyFinal, std, study)


#PlotAllStudyOneGeneMGPcor(gene = "PVALB", MGPname = "GabaVIPReln_Genes")
#PlotAllStudyOneGeneMGPcor(gene = "TAC1", MGPname = "GabaVIPReln_Genes")
#PlotAllStudyOneGeneMGPcor(gene = "KCNK1", MGPname = "GabaVIPReln_Genes")
#PlotAllStudyOneGeneMGPcor(gene = "TAC3", MGPname = "GabaVIPReln_Genes")
#PlotAllStudyOneGeneMGPcor(gene = "SNCG", MGPname = "GabaVIPReln_Genes")

DarmanisExp <- neuroExpressoAnalysis::DarmanisHumanExp %>% data.frame()
DarmanisExp$GeneSymbol <- rownames(DarmanisExp)
DarmanisMeta <- neuroExpressoAnalysis::DarmanisHumanMeta

GetHumanExp <- function(gene){
  gene = toupper(gene)
  temp <- data.frame(GSM = names(DarmanisExp), value = DarmanisExp %>% filter(GeneSymbol == gene) %>% as.numeric)
  temp$CellType <- DarmanisMeta$cellType[match(temp$GSM, DarmanisMeta$GSM)]
  temp <- temp[-nrow(temp), ]
  temp <- temp[!grepl("fetal|hybrid", temp$CellType),]
  temp %<>% mutate(value = log2(value + 1))
  group_by(temp, CellType) %>% summarise(mean(value))
  p <- ggplot(temp, aes(CellType, value))
  plot <- p + labs(x = "", y = paste0(gene, " expression(log2)")) +
    theme_grey(base_size = 16) +
    theme(axis.text.y = element_text(size = rel(1.2)),
          axis.text.x = element_text(size = rel(1.2), angle = -45, hjust = 0),
          legend.text=element_text(size=rel(1)),
          panel.grid = element_blank()) +
    geom_boxplot(outlier.size = 0, width=0.2) + 
    geom_jitter(width = 0.2, size=1, aes(color = CellType, size = 12))
  print(plot)
}


temp <- DarmanisExp %>% filter(GeneSymbol %in% c("PVALB", "SST", "VIP", "CBLN4", "LRRC38", "BTN2A2", "PIP5KL1", "COX6A2", "LHX6", "SOX6"))  %>% t %>% data.frame()
temp2 <- temp[-nrow(temp),]
temp2 <- apply(temp2, 2, as.numeric) %>% data.frame
names(temp2) <-t(temp["GeneSymbol",]) %>% unlist
temp2 %>% filter(LHX6 > 0)

CortexCells <- neuroExpressoAnalysis::mouseMarkerGenes$Cortex %>% names

allStudyCor <- sapply(ls(pat = "_Cortex$"), function(study){
  data = eval(as.name(study))
  data2 <- sapply(names(data), function(celltype){
    cellData <- data[[celltype]]
    DF <- do.call(rbind, cellData)
    DF$CellType <- celltype
    DF$Profile <- sapply(rownames(DF), function(x) strsplit(x, "\\.")[[1]][1])
    DF$Study <- study
    DF
  }, simplify = FALSE)
  do.call(rbind, data2)
}, simplify = FALSE) %>% do.call(rbind, .) %>% data.frame()
allStudyCor %<>% mutate_each(funs(as.factor), -Cor)

markerGenesHuman <- lapply(markerGenes, function(cell){
  mouse2human(cell)$humanGene
})


##### Mean correlation all groups combined ######
allStudyCorMean <- group_by(allStudyCor, CellType, GeneSymbol) %>%
  summarise(meanCor = mean(Cor), median = median(Cor), min = min(Cor), max = max(Cor)) %>%
  data.frame %>% arrange(CellType, desc(meanCor))

allStudyCorMean$PVdevelop <- sapply(allStudyCorMean$GeneSymbol, function(gene){
  if(gene %in% PVdevelopUP){
    "PV developUP"
  } else if (gene %in% PVdevelopDown){
    "PV developDown"
  } else {
    "NS"
  }
})

allStudyCorMean$Marker <- "NO"

for(cell in names(markerGenesHuman)){
  genes <- markerGenesHuman[[cell]]
  cell = paste0(cell, "_Genes")
  allStudyCorMean$Marker[allStudyCorMean$CellType == cell & allStudyGroupCorMean$GeneSymbol %in% genes] <- "YES"
}

allStudyCorMean %<>% mutate(Marker = as.factor(Marker))

allStudyCorMean$CellType <- sapply(allStudyCorMean$CellType, function(x){
  gsub("_Genes", "_MGP", x)
}) %>% factor(levels = c("Astrocyte_MGP", "Microglia_activation_MGP", "Microglia_deactivation_MGP","Microglia_MGP", "Oligo_MGP",
                         "GabaPV_MGP", "GabaRelnCalb_MGP", "GabaVIPReln_MGP",
                         "Pyramidal_Glt_25d2_MGP","Pyramidal_S100a10_MGP", "Pyramidal_Thy1_MGP","PyramidalCorticoThalam_MGP"))

allStudyCorMean$Profile <- relevel(allStudyCorMean$Profile, ref="Cont")

##### Mean correlation groups separated ######
allStudyGroupCorMean <- group_by(allStudyCor, Profile, CellType, GeneSymbol) %>%
  summarise(meanCor = mean(Cor), median = median(Cor), min = min(Cor), max = max(Cor)) %>%
  data.frame %>% arrange(Profile, CellType, desc(meanCor))

allStudyGroupCorMean$PVdevelop <- sapply(allStudyGroupCorMean$GeneSymbol, function(gene){
  if(gene %in% PVdevelopUP){
    "PV developUP"
  } else if (gene %in% PVdevelopDown){
    "PV developDown"
  } else {
    "NS"
  }
})


allStudyGroupCorMean$Marker <- "NO"

for(cell in names(markerGenesHuman)){
  genes <- markerGenesHuman[[cell]]
  cell = paste0(cell, "_Genes")
  allStudyGroupCorMean$Marker[allStudyGroupCorMean$CellType == cell & allStudyGroupCorMean$GeneSymbol %in% genes] <- "YES"
}

allStudyGroupCorMean %<>% mutate(Marker = as.factor(Marker))

allStudyGroupCorMean$CellType <- sapply(allStudyGroupCorMean$CellType, function(x){
  gsub("_Genes", "_MGP", x)
}) %>% factor(levels = c("Astrocyte_MGP", "Microglia_activation_MGP", "Microglia_deactivation_MGP","Microglia_MGP", "Oligo_MGP",
                         "GabaPV_MGP", "GabaRelnCalb_MGP", "GabaVIPReln_MGP",
                         "Pyramidal_Glt_25d2_MGP","Pyramidal_S100a10_MGP", "Pyramidal_Thy1_MGP","PyramidalCorticoThalam_MGP"))

allStudyGroupCorMean$Profile <- relevel(allStudyGroupCorMean$Profile, ref="Cont")

#############################################################################################

ggplot(allStudyGroupCorMean %>% filter(GeneSymbol %in% MistryGenesDown$GeneSymbol,
                                       !CellType %in% c("GabaRelnCalb_MGP",
                                                        "Microglia_activation_MGP",
                                                        "Microglia_deactivation_MGP")), aes(CellType, meanCor))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap(~Profile, ncol=1) +
  geom_violin(width = 0.2) +
  #geom_boxplot(width = 0.1, outlier.size = 0) +
  geom_abline(slope = 0, intercept = 0, color ="red")

ggplot(allStudyGroupCorMean %>% filter(GeneSymbol %in% MistryGenesUP$GeneSymbol,
                                       !CellType %in% c("GabaRelnCalb_MGP",
                                                        "Microglia_activation_MGP",
                                                        "Microglia_deactivation_MGP")), aes(CellType, meanCor))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap(~Profile, ncol=1) +
  geom_violin(width = 0.2) +
  #geom_boxplot(width = 0.1, outlier.size = 0) +
  geom_abline(slope = 0, intercept = 0, color ="red")

ggplot(allStudyGroupCorMean %>% filter(GeneSymbol %in% markerGenesHuman$GabaVIPReln, #!GeneSymbol %in% c("PIP5KL1", "COX6A2", "ADAMTS15"),
                                       !CellType %in% c("GabaRelnCalb_MGP",
                                                        "Microglia_activation_MGP",
                                                        "Microglia_deactivation_MGP")), aes(CellType, meanCor))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~Profile, ncol=1) +
  geom_violin(width = 0.2) +
  #geom_boxplot(width = 0.1, outlier.size = 0) +
  geom_abline(slope = 0, intercept = 0, color ="red")


