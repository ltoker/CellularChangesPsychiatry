source("SetUp.R")
packageF("parallel")

PreProcces2 <- function(Data, meta){
  aned <- Data$aned
  Metadata <- meta
  source(paste0(GenScriptPath,"pre-proccess.R"), local = TRUE)
  output <- list(aned_high, aned_good, aned_low, MaxNoise,Z_scores_High,
                 exclude_samples_low, exclude_samples_high, exclude_samples, Metadata)
  names(output) <- c("aned_high", "aned_good", "aned_low", "NoiseThreshold", "Z_scores_High",
                     "exclude_samples_low", "exclude_samples_high", "exclude_samples", "Metadata")
  return(output)
}

shortPCA <- function(data, celltype){
  PCAresult <- prcomp(t(data %>% dplyr::select(matches("GSM"))), scale = TRUE)
  PCAresult <- correct_sign("GabaPV", PCAresult)
  while (sum(PCAresult$rotation[,1] > 0) < nrow(PCAresult$rotation)){
    minorGene <- rownames(PCAresult$rotation)[PCAresult$rotation[,1] < 0]
    data %<>% filter(!GeneSymbol %in% minorGene)
    rownames(data) <- data$GeneSymbol
    PCAresult <- data %>% select(matches("_|GSM")) %>% t %>% prcomp(scale=TRUE) 
    PCAresult <- correct_sign("GabaPV", PCAresult)
  }
  return(PCAresult)
}

plotDevelopMGP <- function(data, xVal, yVal, title, xlab = "", ylab = "Gene-MGP correlation", txtSize = 12,
                           colors = c("darkgreen", "darkred")){
  p <- ggplot(data,
               aes_string(x = xVal, y=yVal))
  plot <- p +
    labs(title = title, x = xlab, y = ylab)+
    theme_classic(base_size = txtSize) +
    theme(axis.text.y = element_text(size = rel(1.4)),
          axis.text.x = element_text(size = rel(1.4)),
          panel.grid = element_blank(),
          legend.position = "none") +
    scale_fill_manual(values = colors, name="") +
    geom_violin(aes_string(fill = xVal), alpha = 0.6) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_abline(intercept = 0, slope = 0, color = "red", linetype="dashed")
  return(plot)
}

GeneMgpCorDev <- function(dataExp = OkatyExpMelt,
                          dataMGP = GabaPVdevelopMGP$x[,1],
                          matchCol  = "SampleName",
                          gene1 = "Gm17750", gene2 = "Sgpp2",
                          colors = c("darkgreen", "darkred"),
                          title = "Okaty 2009 (GSE17806)\nmouse development data",
                          txtSize = 12){
  dataSub <- dataExp %>% filter(GeneSymbol %in% c(gene1, gene2))
  dataSub$MGP <- dataMGP[match(dataSub[[matchCol]], names(dataMGP))] 
  p1 <- ggplot(dataSub, aes(Age, Exp, color = GeneSymbol)) +
    theme_classic(base_size = txtSize) +
    theme(axis.text.y = element_text(size = rel(1.4)),
          axis.text.x = element_text(size = rel(1.4)),
          legend.position = c(0.15, 1), legend.justification = c("left", "top"),
          legend.direction = "horizontal", legend.title = element_blank(),
          legend.background = element_rect(fill = NA)) +
    scale_y_continuous(limits = c(5, 11.5)) +
    geom_point() +
    labs(title = title, x = "Age", y = "Expression (log2)") +
    scale_color_manual(values = colors) +
    stat_summary(aes(group=GeneSymbol), fun.y=mean, geom="line")
  
  MgpGeneCor <- data.frame(GeneSymbol = c(gene1, gene2),
                           Cor = sapply(c(gene1, gene2), function(gene){
                             cor.test(~MGP+Exp, data = dataSub %>% filter(GeneSymbol == gene),
                                      method = "spearman")$estimate %>%
                               round(digits = 2) %>%
                               paste0('r["S"] == ', .)
                           }))
  p2 <- ggplot(dataSub, aes(MGP, Exp, color = GeneSymbol)) +
    theme_classic(base_size = txtSize) +
    theme(axis.text.y = element_text(size = rel(1.4)),
          axis.text.x = element_text(size = rel(1.4)),
          legend.position = c(0.15, 1), legend.justification = c("left", "top"),
          legend.direction = "horizontal", legend.title = element_blank(),
          legend.background = element_rect(fill = NA)) +
    scale_y_continuous(limits = c(5, 11.5)) +
    geom_point() +
    labs(title = title, x = "GabaPV MGP (AU)" , y = "Expression (log2)") +
    scale_color_manual(values = c("darkgreen", "darkred")) +
    geom_smooth(method='lm', size = 0.5, fullrange=TRUE,  fill = "lightgrey") +
    annotate("text", x=c(-4, -4), y=c(8.5, 7.5), label =  MgpGeneCor$Cor, color = colors, parse = TRUE)
  plot_grid(p1, p2, nrow = 1, align = "vertical")
  return(list(AgeExpPlot = p1, MgpExpPlot = p2))
}


DevDataPath = "DevelopmentData/GSE17806/data/"
resultsPath = "DevelopmentData/GSE17806/"
if(!resultsPath %in% list.dirs(full.names = FALSE)){
  dir.create(resultsPath)
  dir.create(DevDataPath)
  }
softDown("GSE17806", paste0(DevDataPath,"GSE17806.soft"))

if(length(list.files(pattern = "GSE17806meta.tsv", path = DevDataPath)) == 0){
  OkatyMeta <- ReadSoft(paste0(DevDataPath,"GSE17806.soft"))
  OkatyMeta %<>% mutate(Gender = "Male")
  names(OkatyMeta)[4] <- "CellNum"
  write.table(OkatyMeta, file = paste0(DevDataPath,"GSE17806meta.tsv"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
} else {
  OkatyMeta <- read.table(paste0(DevDataPath, "GSE17806meta.tsv"),
                             header = TRUE, sep = "\t")
}

OkatyMeta$GSM <- OkatyMeta$Series_sample_id
OkatyMeta$Profile <- OkatyMeta$age

DownloadCel(OkatyMeta$GSM, path = DevDataPath)
Data <- ReadCell(path = DevDataPath, CelFiles=list.files(path = DevDataPath,pattern="\\.cel",
                                                         ignore.case = TRUE),
                 QC=0, platform=OkatyMeta$Sample_platform_id %>% unique)

OkatyMeta$ScanDate <- Data$scanDate[match(OkatyMeta$GSM, names(Data$scanDate))] 

OkatyStudy <- PreProcces2(Data, OkatyMeta)

OkatyMeta <- OkatyStudy$Metadata 

OkatyMeta$AgeNumeric <- sapply(OkatyMeta$age, function(x) gsub("P", "", x) %>% as.numeric )
OkatyMeta %<>% arrange(AgeNumeric) 
OkatyMeta$age <- factor(OkatyMeta$age, levels = unique(OkatyMeta$age))

OkatyExp <- OkatyStudy$aned_high
rownames(OkatyExp) <- OkatyExp$Probe

OkatyExpMelt <- melt(OkatyExp, id.vars = c("Probe", "GeneSymbol"),
                     measure.vars = grep("GSM", names(OkatyExp), value=T),
                     variable.name = "SampleName", value.name = "Exp")
OkatyExpMelt$Age <- OkatyMeta$age[match(OkatyExpMelt$SampleName, OkatyMeta$Series_sample_id)] %>% ordered()

OkatyExpMelt$AgeDays <- sapply(OkatyExpMelt$Age, function(x){
  as.numeric(as.character(gsub("P", "", x)))
})

temp <- group_by(OkatyExpMelt, Age, Probe) %>% summarise(Mean = mean(Exp)) %>% data.frame
temp <- group_by(temp, Probe) %>% summarise(Max = max(Mean))

lmResults <- mclapply(unique(OkatyExpMelt$Probe),function(probe){ 
  lm.mod <- lm(Exp ~ 1 + Age, data = OkatyExpMelt %>% filter(Probe == probe))
  anv <- anova(lm.mod)
  cbind(t(lm.mod$coefficients), anv["Age",])
},mc.cores = detectCores())

VarExp <- mclapply(unique(OkatyExpMelt$Probe),function(probe){ 
  lm.mod <- lm(Exp ~ 1 + Age, data = OkatyExpMelt %>% filter(Probe == probe))
  temp <- summary(lm.mod)
  temp$r.squared
},mc.cores = 0.5*detectCores()) %>% unlist

VarExp2 <- mclapply(unique(OkatyExpMelt$Probe),function(probe){ 
  data = OkatyExpMelt %>% filter(Probe == probe)
  loess.mod <- loess(Exp ~ 1 + AgeDays, data = data)
  temp <- predict(loess.mod)
  (cor(temp, data$Exp))^2
},mc.cores = 0.5*detectCores()) %>% unlist

names(lmResults) <- OkatyExpMelt$Probe %>% unique()
lmResults <-  rbindlist(lmResults, use.names = TRUE, idcol="Probe")
lmResults$BH <- p.adjust(lmResults$`Pr(>F)`, method = "BH")
lmResults %<>% arrange(BH)

lmResults$GeneSymbol <- OkatyExp$GeneSymbol[match(lmResults$Probe, OkatyExp$Probe)]
lmResults %<>% arrange(Age.L, BH)

PVdevelopUP <- lmResults %>% 
  filter(BH < 0.01, Age.L > 0) %>%
  .$GeneSymbol %>% unique %>% as.character

PVdevelopDown <- lmResults %>% 
  filter(BH < 0.01, Age.L < 0) %>%
  .$GeneSymbol %>% unique %>% as.character()

commonGenes <- intersect(PVdevelopDown, PVdevelopUP)
PVdevelopDown <- PVdevelopDown[!PVdevelopDown %in% commonGenes]
PVdevelopUP <- PVdevelopUP[!PVdevelopUP %in% commonGenes]

markerGabaPV <- mouseMarkerGenesCombined$Cortex$GabaPV

rownames(OkatyExp) <- OkatyExp$GeneSymbol
OkatyGabaPV <- OkatyExp %>% filter(GeneSymbol %in% markerGabaPV)
rownames(OkatyGabaPV) <- OkatyGabaPV$GeneSymbol

GabaPVdevelopMGP <- shortPCA(OkatyGabaPV, "GabaPV_Genes")
genes <- c(PVdevelopUP, PVdevelopDown)

DevelopCor <- sapply(genes[genes %in% OkatyExp$GeneSymbol], function(gene){
  cor(OkatyExp %>% filter(GeneSymbol == gene) %>% select(matches("GSM")) %>% unlist, GabaPVdevelopMGP$x[,1], method = "spearman")
}) %>% data.frame(row.names = genes[genes %in% OkatyExp$GeneSymbol])

names(DevelopCor) <- "Cor"
DevelopCor %<>% mutate(GeneSymbol = rownames(.))
DevelopCor$PVdevelop <- sapply(DevelopCor$GeneSymbol, function(gene){
  if(gene %in% PVdevelopUP){
    "PV developUP"
  } else if(gene %in% PVdevelopDown){
    "PV developDown"
  } else {
    NA
  }
})

#Plot example for genes down and upregulated during development
GeneMGPcorDevPlots <- GeneMgpCorDev()

#Plot Okaty mouse GabaPV development plot
GeneMGPcorDevPlots$OkatyPlot <- plotDevelopMGP(data = DevelopCor, xVal = "PVdevelop", yVal = "Cor",
                            title = "Okaty 2009 (GSE17806)\nmouse development data")




####### Human development #######
markerGabaPVhuman <- sapply(markerGabaPV, function(x){
  mouse2human(x)$humanGene
})

markerGabaPVhuman <- markerGabaPVhuman[!markerGabaPVhuman %in% c("WIF1", "TMEM132C", "BTN2A2")]

PVdevelopUPhuman <- sapply(PVdevelopUP, function(gene) {
  mouse2human(gene)$humanGene
}) %>% unlist()

PVdevelopDownhuman <- sapply(PVdevelopDown, function(gene) {
  mouse2human(gene)$humanGene
}) %>% unlist()

DevDataPath = "DevelopmentData/GSE25219/data/"
resultsPath = "DevelopmentData/GSE25219/"
softDown("GSE25219", paste0(DevDataPath,"GSE25219.soft"))
if(length(list.files(pattern = "GSE25219meta.tsv", path = DevDataPath)) == 0){
  GSE25219meta <- ReadSoft(paste0(DevDataPath,"GSE25219.soft"))
  write.table(GSE25219meta, file = paste0(DevDataPath,"GSE25219meta.tsv"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
} else {
  GSE25219meta <- read.table(paste0(DevDataPath, "GSE25219meta.tsv"),
                             header = TRUE, sep = "\t")
}

GSE25219meta$GSM <- GSE25219meta$Series_sample_id
GSE25219meta %<>% filter(Sample_platform_id == "GPL5175")
GSE25219meta$Age_years <- sapply(GSE25219meta$age, function(x){
  if(grepl("PCW", x)){
    x = strsplit(as.character(x), " ")[[1]][1] %>% as.numeric()
    signif(-1*x/52, digits = 2) 
  } else if(grepl("PMonth", x)){
    x = strsplit(as.character(x), " ")[[1]][1] %>% as.numeric()
    signif(x/12, digits = 2)
  } else {
    strsplit(as.character(x), " ")[[1]][1] %>% as.numeric()
  }
})

GSE25219meta %<>% filter(Age_years < 50, region == "DFC") %>% droplevels()
GSE25219meta <- GSE25219meta[!duplicated(GSE25219meta$brain.code),]

names(GSE25219meta) <- sapply(names(GSE25219meta), function(x) gsub("\\.", "_", x))

DownloadCel(GSE25219meta$GSM, path = DevDataPath)
Data <- ReadCell(path = DevDataPath, CelFiles=list.files(path = DevDataPath,pattern="\\.cel",
                                                         ignore.case = TRUE),
                 QC=0, platform=GSE25219meta$Sample_platform_id %>% unique %>% .[1])

GSE25219meta$ScanDate <- Data$scanDate[match(GSE25219meta$GSM, names(Data$scanDate))] 

GSE25219data <- PreProcces2(Data, GSE25219meta)

GSE25219GabaPV <- GSE25219data$aned_high %>% filter(GeneSymbol %in% markerGabaPVhuman)
rownames(GSE25219GabaPV) <- GSE25219GabaPV$GeneSymbol

GSE25219GabaPVMGP <- shortPCA(GSE25219GabaPV, "GabaPV_Genes")


genesHuman <- c(PVdevelopUPhuman, PVdevelopDownhuman) 

GSE25219DevelopCor <- sapply(genesHuman[genesHuman %in% GSE25219data$aned_high$GeneSymbol], function(gene){
  cor(GSE25219data$aned_high %>% filter(GeneSymbol == gene) %>% select(matches("GSM")) %>% unlist, GSE25219GabaPVMGP$x[,1], method = "spearman")
}) %>% data.frame(row.names = genesHuman[genesHuman %in% GSE25219data$aned_high$GeneSymbol])

names(GSE25219DevelopCor) <- "Cor"
GSE25219DevelopCor %<>% mutate(GeneSymbol = rownames(.))
GSE25219DevelopCor$PVdevelop <- sapply(GSE25219DevelopCor$GeneSymbol, function(gene){
  if(gene %in% PVdevelopUPhuman){
    "PV developUP"
  } else if(gene %in% PVdevelopDownhuman){
    "PV developDown"
  } else {
    NA
  }
})

GeneMGPcorDevPlots$KangPlot <- plotDevelopMGP(data = GSE25219DevelopCor, xVal = "PVdevelop", yVal = "Cor",
                           title = "Kang 2011 (GSE25219)\nhuman development data")


###### GSE13564 #########################
DevDataPath = "DevelopmentData/GSE13564/data/"
resultsPath = "DevelopmentData/GSE13564/"
softDown("GSE13564", paste0(DevDataPath,"GSE13564.soft"))
if(length(list.files(pattern = "GSE13564meta.tsv", path = DevDataPath)) == 0){
  GSE13564meta <- ReadSoft(paste0(DevDataPath,"GSE13564.soft"))
  GSE13564meta %<>% mutate(Age_years = sapply(as.character(.$V2), function(x) strsplit(x, ",")[[1]][2] %>% as.numeric),
                           PMI = sapply(as.character(.$V2), function(x) strsplit(x, ",")[[1]][4] %>% as.numeric),
                           Gender = sapply(as.character(.$V2), function(x) strsplit(x, ",")[[1]][6]))
  write.table(GSE13564meta, file = paste0(DevDataPath,"GSE13564meta.tsv"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
} else {
  GSE13564meta <- read.table(paste0(DevDataPath, "GSE13564meta.tsv"),
                             header = TRUE, sep = "\t")
}

GSE13564meta$GSM <- GSE13564meta$Series_sample_id

DownloadCel(GSE13564meta$GSM, path = DevDataPath)
Data <- ReadCell(path = DevDataPath, CelFiles=list.files(path = DevDataPath,pattern="\\.cel",
                                                         ignore.case = TRUE),
                 QC=0, platform=GSE13564meta$Sample_platform_id %>% unique)

GSE13564meta$ScanDate <- Data$scanDate[match(GSE13564meta$GSM, names(Data$scanDate))] 

GSE13564data <- PreProcces2(Data, GSE13564meta)

GSE13564GabaPV <- GSE13564data$aned_high %>% filter(GeneSymbol %in% markerGabaPVhuman)

geneRM <- sapply(GSE13564GabaPV$GeneSymbol, function(gene){
  if(GSE13564GabaPV %>% filter(GeneSymbol == gene) %>% select(matches("GSM")) %>% unlist %>% quantile( 0.67) < GSE13564data$NoiseThreshold){
    TRUE
  } else {
    FALSE
  }
})

GSE13564GabaPV <- GSE13564GabaPV[!geneRM,]

rownames(GSE13564GabaPV) <- GSE13564GabaPV$GeneSymbol

GSE13564GabaPVMGP <- shortPCA(GSE13564GabaPV, "GabaPV_Genes")

GSE13564DevelopCor <- sapply(genesHuman[genesHuman %in% GSE13564data$aned_high$GeneSymbol], function(gene){
  cor(GSE13564data$aned_high %>% filter(GeneSymbol == gene) %>% select(matches("GSM")) %>% unlist, GSE13564GabaPVMGP$x[,1], method = "spearman")
}) %>% data.frame(row.names = genesHuman[genesHuman %in% GSE13564data$aned_high$GeneSymbol])

names(GSE13564DevelopCor) <- "Cor"
GSE13564DevelopCor %<>% mutate(GeneSymbol = rownames(.))
GSE13564DevelopCor$PVdevelop <- sapply(GSE13564DevelopCor$GeneSymbol, function(gene){
  if(gene %in% PVdevelopUPhuman){
    "PV developUP"
  } else if(gene %in% PVdevelopDownhuman){
    "PV developDown"
  } else {
    NA
  }
})

GeneMGPcorDevPlots$HarrisPlot <- plotDevelopMGP(data = GSE13564DevelopCor, xVal = "PVdevelop", yVal = "Cor",
                             title = "Harris 2009 (GSE13564)\nhuman development data")

ggsave("GeneralResults/HarrisDevelopment.pdf", plot = GeneMGPcorDevPlots$HarrisPlot, width = 6, height = 6, units = "in",
       dpi=300)


temp <- merge(GSE25219DevelopCor, GSE13564DevelopCor, by="GeneSymbol", suffixes = c("GSE25219", "GSE13564"))
temp %<>% mutate(PVdevelopGSE25219 = as.factor(PVdevelopGSE25219), PVdevelopGSE13564 = as.factor(PVdevelopGSE13564))
ggplot(temp, aes(CorGSE25219,CorGSE13564, color = PVdevelopGSE13564)) + geom_point()

## Psychiatry data ######
#Combine MGP correlations from all studies
load("AnalysisStudies.rda")
StanleyStudies <- read.table("StanleyStudies.txt", sep="\t", header=T)

for(name in AnalysisStudies){
  load(paste0("GeneralResults/MGPcorAllgenes_", name,".rda"))
  rm(list = ls(pat=paste0(name, "_Cortex$")))
  MGPcorAllgenes <- lapply(MGPcorAllgenes, function(std){
    lapply(std, function(cellType){
      lapply(cellType, function(grp){
        grp$GeneSymbol <- rownames(grp)
        grp %>% arrange(desc(Cor))
      })
    })
  })
  
  if(grepl("Stanley", name)){
    for(study in names(MGPcorAllgenes)){
      investigator = gsub("Study.[0-9]?", "", study)
      region = StanleyStudies %>%
        filter(Investigator == investigator) %>%
        select(NeuExpRegion) %>%
        unlist %>% as.character
      assign(paste0(study, "_", region), MGPcorAllgenes[[study]])
    }
  } else {
    assign(paste0(name, "_Cortex"), MGPcorAllgenes$Cortex)
  }
}

allStudyCor <- sapply(ls(pat = "_Cortex$"), function(study){
  data = eval(as.name(study))
  data2 <- sapply(names(data), function(celltype){
    cellData <- data[[celltype]]
    DF <- do.call(rbind, cellData)
    DF$CellType <- celltype
    DF$Profile <- sapply(rownames(DF), function(x) strsplit(x, "\\.")[[1]][1]) %>% factor()
    DF$Study <- study
    DF
  }, simplify = FALSE)
  do.call(rbind, data2)
}, simplify = FALSE) %>% do.call(rbind, .) %>% data.frame()
allStudyCor %<>% mutate_each(funs(as.factor), -Cor)

allStudyCor$Profile <- relevel(allStudyCor$Profile, ref = "Cont")

allStudyCor$CellType <- sapply(allStudyCor$CellType, function(x){
  gsub("_Genes", "_MGP", x)
})

allStudyCorMean <- group_by(allStudyCor, CellType, GeneSymbol) %>%
  summarise(meanCor = mean(Cor)) %>%
  data.frame %>% arrange(CellType, desc(meanCor))

allStudyCorMean$PVdevelop <- sapply(allStudyCorMean$GeneSymbol, function(gene){
  if(gene %in% PVdevelopUPhuman){
    "PV developUP"
  } else if (gene %in% PVdevelopDownhuman){
    "PV developDown"
  } else {
    "NS"
  }
})


GeneMGPcorDevPlots$PsychiatryMeanPlot <- plotDevelopMGP(data = allStudyCorMean %>% filter(PVdevelop != "NS") %>% droplevels(),
                                     xVal = "PVdevelop", yVal = "meanCor",
                                     title = "Mean Psychiatry data")

allGroupCorMean <- group_by(allStudyCor, Study, CellType, GeneSymbol) %>%
  summarise(meanCor = mean(Cor)) %>%
  data.frame %>% arrange(CellType, desc(meanCor))

allGroupCorMean$PVdevelop <- sapply(allGroupCorMean$GeneSymbol, function(gene){
  if(gene %in% PVdevelopUPhuman){
    "PV developUP"
  } else if (gene %in% PVdevelopDownhuman){
    "PV developDown"
  } else {
    "NS"
  }
})

p <- ggplot(allGroupCorMean %>% filter(PVdevelop != "NS") %>% droplevels(),
             aes(x = PVdevelop, y=meanCor))
PsychiatrySinglePlot <- p + theme_bw(base_size = 10) +
  labs(title = "Psychiatry data", x = "", y = "Gene-MGP correlation") +
  theme(axis.text.y = element_text(size = rel(1.1)),
        axis.text.x = element_text(size = rel(1.1)),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("darkgreen", "darkred"), name="") +
  geom_violin(aes(fill = PVdevelop), alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.size = 0) +
  geom_abline(intercept = 0, slope = 0, color = "red", linetype="dashed") +
  facet_wrap(~Study, nrow = 3)
ggsave("GeneralResults/PsychiatrySingle.pdf", plot = PsychiatrySinglePlot, width = 10,  height = 6, units = "in",
       dpi=300)

pdf("GeneralResults/DevelopmentFigure.pdf", width = 12, height = 8,  pointsize = 12, useDingbats = FALSE)
plot_grid(GeneMGPcorDevPlots$AgeExpPlot, GeneMGPcorDevPlots$MgpExpPlot, GeneMGPcorDevPlots$OkatyPlot,
          GeneMGPcorDevPlots$KangPlot, GeneMGPcorDevPlots$HarrisPlot, GeneMGPcorDevPlots$PsychiatryMeanPlot,
          nrow = 2, labels = c("A", "B", "C", "D", "E", "F"))
dev.off()

ggsave("GeneralResults/DevelopmentFigure.pdf", width = 10,  height = 6, units = "in",
       dpi=300)

save(allGroupCorMean, allStudyCorMean, allStudyCor, file = "GeneralResults/GeneMGPcor.rda")
save.image(file = "GeneralResults/Development.rda")
