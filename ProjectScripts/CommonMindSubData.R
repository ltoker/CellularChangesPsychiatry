source("SetUp.R")
packageF("parallel")
source(paste0(GenScriptPath, "Cell_type_PCA.R"))

load("GeneralResults/studyFinalCommonMind.rda")
name = "CommonMind"

studyFinal$BA9$Metadata$Study <- "CommonMind"
PC1allGenes <- prcomp(as.matrix(studyFinal$BA9$aned_high %>% select(matches("_")) %>% t),
                      scale = T)
studyFinal$BA9$Metadata$PC1allGenes <- PC1allGenes$x[match(studyFinal$BA9$Metadata$CommonName, rownames(PC1allGenes$x)),1]
studyFinal$BA9$Metadata$PC1allGenes <- rescale(studyFinal$BA9$Metadata$PC1allGenes, c(0,1))

AstroMGPallData <- studyFinal$BA9$Metadata$Astrocyte_Genes[match(names(studyFinal$BA9$aned_high)[-c(1:3)],
                                                                 studyFinal$BA9$Metadata$CommonName)]
GeneMGPcorCommonMindAll <- apply(studyFinal$BA9$aned_high %>% select(matches("_")), 1, function(Gene){
  cor(Gene, AstroMGPallData, method = "spearman")
})
names(GeneMGPcorCommonMindAll) <- studyFinal$BA9$aned_high$GeneSymbol

HighRINdata <- studyFinal$BA9$Metadata %>% filter(RIN > 8) %>% droplevels()
HighRINexp <- studyFinal$BA9$aned_high %>% select_(.dots = c("GeneSymbol", "Probe", "ensemblID",
                                                             as.character(HighRINdata$CommonName)))
HighRIN_MGPresults <- PCA_genes_All_based(dataset_id="HighRIN",
                                  dataset=HighRINexp,
                                  CellType_genes=CellType_genes[c("Astrocyte_Genes", "GabaPV_Genes")],
                                  NoiseThershold = studyFinal$BA9$NoiseThreshold)
HighRINdata$HighRIN_AstroMGP <- HighRIN_MGPresults$All$Astrocyte_Genes$x[match(HighRINdata$CommonName,
                                                                               rownames(HighRIN_MGPresults$All$Astrocyte_Genes$x)),1]
HighRINdata$HighRIN_GabaPVMGP <- HighRIN_MGPresults$All$GabaPV_Genes$x[match(HighRINdata$CommonName,
                                                                               rownames(HighRIN_MGPresults$All$GabaPV_Genes$x)),1]
HighRINdata$HighRIN_AstroMGP <- rescale(HighRINdata$HighRIN_AstroMGP, c(0,1))
HighRINdata$HighRIN_GabaPVMGP <- rescale(HighRINdata$HighRIN_GabaPVMGP, c(0,1))

# Run MGP analysis on subsampled data 
#####################################
ContSamp <- grep("Cont", studyFinal$BA9$Metadata$CommonName, value = T)
SczSamp <- grep("SCZ", studyFinal$BA9$Metadata$CommonName, value = T)

CellType_genes <- GetMarkers("Cortex")
#Exclude GabaPV genes which are not neuron specific in human (Darmanis) data
CellType_genes$GabaPV_Genes <- CellType_genes$GabaPV_Genes[!CellType_genes$GabaPV_Genes %in% c("WIF1", "TMEM132C", "BTN2A2")]

set.seed(16)
Rot <- as.list(1:10000)
names(Rot) <- paste0("Rot_", 1:10000)

#Subsamppled metadata
SubData <- mclapply(Rot, function(x){
  samples = c(sample(ContSamp, 30, replace = FALSE), sample(SczSamp, 30, replace = FALSE))
  studyFinal$BA9$Metadata[studyFinal$BA9$Metadata$CommonName %in% samples,] %>% droplevels()
}, mc.cores = 0.8*detectCores())

#Run MGP analisys
SubDataMGP <- mclapply(SubData, function(x){
  aned_high <- studyFinal$BA9$aned_high %>% select_(.dots =  c("GeneSymbol", "Probe", "ensemblID", as.character(x$CommonName)))
  MGPresults <- PCA_genes_All_based(dataset_id=x,
                                    dataset=aned_high,
                                    CellType_genes=CellType_genes[c("Astrocyte_Genes", "GabaPV_Genes")],
                                    NoiseThershold = studyFinal$BA9$NoiseThreshold)
  
  PC1allGenes <- prcomp(as.matrix(aned_high %>% select(matches("_")) %>% t),
                        scale = T)
  PC1allGenes <- PC1allGenes$x[,1] %>% rescale(c(0,1))
  list(MGP = MGPresults,
       PC1all = PC1allGenes)
},  mc.cores = 0.8*detectCores())

load(paste0(GeneralResultsPath,"GeneMGPcor.rda"))

SubData2 <- sapply(names(Rot), function(RotNum){
  aned_high <- studyFinal$BA9$aned_high %>% select_(.dots =  c("GeneSymbol", "Probe", "ensemblID", as.character(SubData[[RotNum]]$CommonName)))
  Data = SubData[[RotNum]]
  Data$SubAstrocyte_Genes <- SubDataMGP[[RotNum]]$MGP$All$Astrocyte_Genes$x[match(SubData[[RotNum]]$CommonName,
                                                                          rownames(SubDataMGP[[RotNum]]$MGP$All$Astrocyte_Genes$x)),1]
  Data$SubGabaPV_Genes <- SubDataMGP[[RotNum]]$MGP$All$GabaPV_Genes$x[match(SubData[[RotNum]]$CommonName,
                                                                              rownames(SubDataMGP[[RotNum]]$MGP$All$GabaPV_Genes$x)),1]
  Data$SubAstrocyte_Genes <- rescale(Data$SubAstrocyte_Genes, c(0,1))
  Data$SubGabaPV_Genes <- rescale(Data$SubGabaPV_Genes, c(0,1))
  Data$SubPC1allGenes <- SubDataMGP[[RotNum]]$PC1all[match(SubData[[RotNum]]$CommonName,
                                                        names(SubDataMGP[[RotNum]]$PC1all))]
  
  GeneMGPcor <- apply(aned_high %>% select(matches("_")), 1, function(Gene){
    cor(Gene, Data$SubAstrocyte_Genes, method = "spearman")
  })
  
  names(GeneMGPcor) <- aned_high$GeneSymbol
  list(GeneMGPcor = GeneMGPcor, Data = Data)
}, simplify = FALSE)

SubPC1fsPVcor <- mclapply(SubData2, function(RotNum){
  cor(RotNum$Data$SubGabaPV_Genes, RotNum$Data$SubPC1allGenes, method = "spearman")
}, mc.cores = 0.5*detectCores()) %>% unlist

SubPC1Astrocor <- mclapply(SubData2, function(RotNum){
  cor(RotNum$Data$SubAstrocyte_Genes, RotNum$Data$SubPC1allGenes, method = "spearman")
}, mc.cores = 0.5*detectCores()) %>% unlist

SubRINfsPVcor <- mclapply(SubData2, function(RotNum){
  cor(RotNum$Data$SubGabaPV_Genes, RotNum$Data$RIN, method = "spearman")
}, mc.cores = 0.5*detectCores()) %>% unlist

SubRINAstrocor <- mclapply(SubData2, function(RotNum){
  cor(RotNum$Data$SubAstrocyte_Genes, RotNum$Data$RIN, method = "spearman")
}, mc.cores = 0.5*detectCores()) %>% unlist

SubPC1allRINcor <- mclapply(SubData2, function(RotNum){
  cor(RotNum$Data$SubPC1allGenes, RotNum$Data$RIN, method = "spearman")
}, mc.cores = 0.5*detectCores()) %>% unlist

SubCorDF <- data.frame(Study = "CommonMind",
                       RotNum = names(SubPC1fsPVcor),
                       PC1fsPVcor = abs(SubPC1fsPVcor),
                       RINfsPVcor = abs(SubRINfsPVcor),
                       PC1Astrocor = abs(SubPC1Astrocor),
                       RINAstrocor = abs(SubRINAstrocor),
                       PC1allRINcor = abs(SubPC1allRINcor))

WilcoxShiftAstro <- mclapply(SubData2, function(Rot){
  x = Rot$Data
  Stat <- wilcox.test(SubAstrocyte_Genes~Profile, data = x, conf.int=T)
  -1*Stat$estimate
}, mc.cores = 0.5*detectCores()) %>% unlist

WilcoxShiftPV <- mclapply(SubData2, function(Rot){
  x = Rot$Data
  Stat <- wilcox.test(SubGabaPV_Genes~Profile, data = x, conf.int=T)
  -1*Stat$estimate
}, mc.cores = 0.5*detectCores()) %>% unlist

WilcoxShiftRIN <- mclapply(SubData2, function(Rot){
  x = Rot$Data
  Stat <- wilcox.test(RIN~Profile, data = x, conf.int=T)
  -1*Stat$estimate
}, mc.cores = 0.5*detectCores()) %>% unlist


VarExplainedAstro <- mclapply(SubDataMGP, function(Rot){
  x = Rot$MGP
  x$All$Astrocyte_Genes %>% summary() %>% .$importance %>% data.frame() %>% .[2,1:3] %>% unlist
}, mc.cores = 0.3*detectCores()) %>% do.call(rbind, .) %>% data.frame()

VarExplainedAstro %<>% mutate(Rot = rownames(.))


VarExplainedPV <- mclapply(SubDataMGP, function(Rot){
  x = Rot$MGP
  x$All$GabaPV_Genes %>% summary() %>% .$importance %>% data.frame() %>% .[2,1:3] %>% unlist
}, mc.cores = 0.3*detectCores()) %>% do.call(rbind, .) %>% data.frame()

VarExplainedPV %<>% mutate(Rot = rownames(.))

VarExplainedAstro$ShiftAstro <- WilcoxShiftAstro
VarExplainedAstro$ShiftPV <- WilcoxShiftPV
VarExplainedAstro$ShiftRIN <- WilcoxShiftRIN

VarExplainedPV$ShiftAstro <- WilcoxShiftAstro
VarExplainedPV$ShiftPV <- WilcoxShiftPV
VarExplainedPV$ShiftRIN <- WilcoxShiftRIN


packageF("zoo")
temp <- VarExplainedAstro %>% arrange(PC1)
SliWindAstro <- data.frame(PC1 = rollapply(temp$PC1, width = 10, by = 5, FUN = mean, align = "left"),
                           WlCxDelta = rollapply(temp$ShiftAstro, width = 10, by = 5, FUN = mean, align = "left"),
                           Study = "CommonMind")
ggplot(SliWindAstro, aes(PC1, WlCxDelta)) + geom_point()

temp2 <- VarExplainedPV %>% arrange(PC1)
PC1slidePV <- rollapply(temp2$PC1, width = 10, by = 5, FUN = mean, align = "left")
SliWindPV <- data.frame(PC1 = rollapply(temp2$PC1, width = 10, by = 5, FUN = mean, align = "left"),
                           WlCxDelta = rollapply(temp$ShiftPV, width = 10, by = 5, FUN = mean, align = "left"),
                           Study = "CommonMind")

load("GeneralResults/PCAresultsList.rda")
load("GeneralResults/DataCombined.rda")

names(PCAresults) <- sapply(names(PCAresults), function(x){
  gsub("PCAresults|Cortex", "", x)
})

VarExplainedAstroAll <- mclapply(PCAresults[names(PCAresults) != "Study4Chen"], function(study){
  VarExplainedAstro <- mclapply(study, function(x){
    x$All$Astrocyte_Genes %>% summary() %>% .$importance %>% data.frame() %>% .[2,1:3] %>% unlist
    }, mc.cores = 0.3*detectCores()) %>% do.call(rbind, .) %>% data.frame()
  }, mc.cores = length(names(PCAresults)))
  
WilcoxShiftAstro_All <- mclapply(PCAresults[names(PCAresults) != "Study4Chen"], function(study){
  mclapply(study, function(Rot){
    RotData = data.frame(CommonName = rownames(Rot$All$Astrocyte_Genes$x),
                         Profile = sapply(rownames(Rot$All$Astrocyte_Genes$x), function(Name){
                           strsplit(Name, "_")[[1]][1]
                         }),
                         Astrocyte_Genes = Rot$All$Astrocyte_Genes$x[,1])
    RotData$Astrocyte_Genes <-  rescale(RotData$Astrocyte_Genes, c(0,1))
    RotData %<>% filter(Profile != "BP") %>% droplevels()
    Stat <- wilcox.test(Astrocyte_Genes~Profile, data = RotData, conf.int=T)
    -1*Stat$estimate
  }, mc.cores = 0.5*detectCores()) %>% unlist
}, mc.cores = length(PCAresults)) 

for(study in grep("Study4Chen", names(PCAresults), value = T, invert = T)){
  VarExplainedAstroAll[[study]]$ShiftAstro <- WilcoxShiftAstro_All[[study]]
}

SliWindAstroAll <- sapply(names(VarExplainedAstroAll), function(study){
  temp <- VarExplainedAstroAll[[study]] %>% arrange(PC1)
  data.frame(PC1 = rollapply(temp$PC1, width = 10, by = 5, FUN = mean, align = "left"),
             WlCxDelta = rollapply(temp$ShiftAstro, width = 10, by = 5, FUN = mean, align = "left"),
             Study = study)
}, simplify = FALSE) %>% do.call(rbind, .)

SliWindAstroAll <- rbind(SliWindAstro, SliWindAstroAll)
SliWindAstroAll %<>% mutate(PC1 = 100*PC1)

ggplot(SliWindAstroAll %>% filter(Study != "Study2AltarC") , aes(PC1, WlCxDelta, color = Study)) +
  theme_classic(base_size = 16) +
  labs(y = "Change in Astrocyte MGP", x = "MGP variance explained (%)") +
  geom_point() +
  scale_color_manual(values=c("black", "burlywood4", "deeppink", "darkorchid1", "chocolate1", "cadetblue",
                              "cornflowerblue", "brown1", "darkgoldenrod1", "darkgreen"), name = "Dataset")
ggsave("AstroPC1deltaShift.pdf", path = GeneralResultsPath, width = 8, height = 4, units = "in", dpi=300)

VarExplainedPVAll <- mclapply(PCAresults[names(PCAresults) != "Study4Chen"], function(study){
  VarExplainedPV <- mclapply(study, function(x){
    x$All$GabaPV_Genes %>% summary() %>% .$importance %>% data.frame() %>% .[2,1:3] %>% unlist
  }, mc.cores = 0.3*detectCores()) %>% do.call(rbind, .) %>% data.frame()
}, mc.cores = length(names(PCAresults)))

WilcoxShiftPV_All <- mclapply(PCAresults[names(PCAresults) != "Study4Chen"], function(study){
  mclapply(study, function(Rot){
    RotData = data.frame(CommonName = rownames(Rot$All$GabaPV_Genes$x),
                         Profile = sapply(rownames(Rot$All$GabaPV_Genes$x), function(Name){
                           strsplit(Name, "_")[[1]][1]
                         }),
                         GabaPV_Genes = Rot$All$GabaPV_Genes$x[,1])
    RotData$GabaPV_Genes <-  rescale(RotData$GabaPV_Genes, c(0,1))
    RotData %<>% filter(Profile != "BP") %>% droplevels()
    Stat <- wilcox.test(GabaPV_Genes~Profile, data = RotData, conf.int=T)
    -1*Stat$estimate
  }, mc.cores = 0.5*detectCores()) %>% unlist
}, mc.cores = length(PCAresults)) 

for(study in grep("Study4Chen", names(PCAresults), value = T, invert = T)){
  VarExplainedPVAll[[study]]$ShiftPV <- WilcoxShiftPV_All[[study]]
}

SliWindPVAll <- sapply(names(VarExplainedPVAll), function(study){
  temp <- VarExplainedPVAll[[study]] %>% arrange(PC1)
  data.frame(PC1 = rollapply(temp$PC1, width = 10, by = 5, FUN = mean, align = "left"),
             WlCxDelta = rollapply(temp$ShiftPV, width = 10, by = 5, FUN = mean, align = "left"),
             Study = study)
}, simplify = FALSE) %>% do.call(rbind, .)

SliWindPVAll <- rbind(SliWindPV, SliWindPVAll)
SliWindPVAll %<>% mutate(PC1 = 100*PC1)

ggplot(SliWindPVAll %>% filter(Study != "Study2AltarC") , aes(PC1, WlCxDelta, color = Study)) +
  theme_classic(base_size = 16) +
  labs(y = "Change in fsPV MGP", x = "MGP variance explained (%)") +
  geom_point() +
  scale_color_manual(values=c("black", "burlywood4", "deeppink", "darkorchid1", "chocolate1", "cadetblue",
                              "cornflowerblue", "brown1", "darkgoldenrod1", "darkgreen"), name = "Dataset")
ggsave("GabaPVPC1deltaShift.pdf", path = GeneralResultsPath, width = 8, height = 4, units = "in", dpi=300)

#Get PC1 based on all genes for each subsampling rotation
PC1allDataRot <- sapply(names(PCAresults), function(study){
  data = ExpAll[[study]]
  meta = metaCombined %<>% filter(Study == study)
  PCAresultsStudy = PCAresults[[study]]
  names(PCAresultsStudy) <- paste0("Rot_", 1:length(PCAresultsStudy))
  PC1allDataRot <- mclapply(PCAresultsStudy, function(MGPstudy){
    samples = rownames(MGPstudy$All$Astrocyte_Genes$x)
    RotData = data %>% select_(.dots = samples)
    RotMeta = meta %>% filter(CommonName %in% samples)
    PCA <- prcomp(RotData %>% t(), scale = T)
    PC1 <- PCA$x[,1]
    PC1 <- rescale(PC1, c(0,1))
    temp <- data.frame(CommonName = names(RotData),
                       PC1allGenes = PC1,
                       Study = study)
    temp$Profile <- sapply(as.character(temp$CommonName), function(subj){
      strsplit(subj, "_")[[1]][1]
    })
    temp$Astrocyte_Genes <- RotMeta$Astrocyte_Genes[match(temp$CommonName, RotMeta$CommonName)]
    temp$GabaPV_Genes <- RotMeta$GabaPV_Genes[match(temp$CommonName, RotMeta$CommonName)]
    temp$RIN <- RotMeta$RIN[match(temp$CommonName, RotMeta$CommonName)]
    temp$Profile <- relevel(as.factor(temp$Profile), ref = "Cont")
    temp
  }, mc.cores = 0.5*detectCores())
  PC1allDataRot
}, simplify = FALSE)

RotCorDF <- lapply(PC1allDataRot, function(study){
  SubPC1fsPVcor <- mclapply(study, function(RotNum){
    cor(RotNum$GabaPV_Genes, RotNum$PC1allGenes, method = "spearman")
  }, mc.cores = 0.5*detectCores()) %>% unlist
  
  SubPC1Astrocor <- mclapply(study, function(RotNum){
    cor(RotNum$Astrocyte_Genes, RotNum$PC1allGenes, method = "spearman")
  }, mc.cores = 0.5*detectCores()) %>% unlist
  
  SubRINfsPVcor <- mclapply(study, function(RotNum){
    cor(RotNum$GabaPV_Genes, RotNum$RIN, method = "spearman")
  }, mc.cores = 0.5*detectCores()) %>% unlist
  
  SubRINAstrocor <- mclapply(study, function(RotNum){
    cor(RotNum$Astrocyte_Genes, RotNum$RIN, method = "spearman")
  }, mc.cores = 0.5*detectCores()) %>% unlist
  
  SubPC1allRINcor <- mclapply(study, function(RotNum){
    cor(RotNum$PC1allGenes, RotNum$RIN, method = "spearman")
  }, mc.cores = 0.5*detectCores()) %>% unlist
  
  SubCorDF <- data.frame(RotNum = names(SubPC1fsPVcor),
                         PC1fsPVcor = abs(SubPC1fsPVcor),
                         RINfsPVcor = abs(SubRINfsPVcor),
                         PC1Astrocor = abs(SubPC1Astrocor),
                         RINAstrocor = abs(SubRINAstrocor),
                         PC1allRINcor = abs(SubPC1allRINcor))
  SubCorDF
})

for(study in names(RotCorDF)){
  RotCorDF[[study]]$Study <- study
  RotCorDF[[study]] <- RotCorDF[[study]] %>% select_(.dots = names(SubCorDF))
  }

RotCorDF <- do.call(rbind, RotCorDF)

RotCorDFcombined <- rbind(SubCorDF, RotCorDF)
RotCorDFcombined$Study <- relevel(RotCorDFcombined$Study, ref = "CommonMind")

MedianCor <- group_by(RotCorDFcombined, Study) %>%
  summarise(MedianPC1fsPV = median(PC1fsPVcor),
            MedianPC1Astro = median(PC1Astrocor)) %>%
  data.frame() %>% arrange(desc(MedianPC1fsPV))

RotCorDFcombined$Study <- factor(RotCorDFcombined$Study, levels = MedianCor$Study)

ggplot(RotCorDFcombined, aes(Study, PC1fsPVcor)) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Study", y = "PC1 fsPV MGP cor") +
  geom_boxplot(outlier.size = 1)
ggsave("PC1allGabaPVcorRot.pdf", path = GeneralResultsPath, width = 8, height = 4, units = "in", dpi=300)

ggplot(RotCorDFcombined %>% filter(!is.na(PC1allRINcor)), aes(Study, PC1allRINcor)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Study", y = "PC1 RIN cor") +
  geom_boxplot(outlier.size = 1)
ggsave("PC1allRINcorRot.pdf", path = GeneralResultsPath, width = 8, height = 4, units = "in", dpi=300)

ggplot(RotCorDFcombined %>% filter(!is.na(RINfsPVcor)), aes(Study, RINfsPVcor)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Study", y = "RIN fsPV cor") +
  geom_boxplot(outlier.size = 1)
ggsave("RINGabaPVcorRot.pdf", path = GeneralResultsPath, width = 8, height = 2, units = "in", dpi=300)

ggplot(RotCorDFcombined %>% filter(!is.na(RINfsPVcor)), aes(Study, RINAstrocor)) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Study", y = "RIN Astro cor") +
  geom_boxplot(outlier.size = 1)

MedianCor %<>% arrange(desc(MedianPC1Astro))
RotCorDFcombined$Study <- factor(RotCorDFcombined$Study, levels = MedianCor$Study)

ggplot(RotCorDFcombined, aes(Study, PC1Astrocor)) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Study", y = "PC1 Astro MGP cor") +
  geom_boxplot(outlier.size = 1)

#Get PC1 based on all genes on the aggregated subsampling
PC1allData <- sapply(names(PCAresults), function(study){
  data = ExpAll[[study]]
  meta = metaCombined %<>% filter(Study == study)
  PCA <- prcomp(data %>% select(matches("GSM|_")) %>% t(), scale = T)
  PC1 <- PCA$x[,1]
  PC1 <- rescale(PC1, c(0,1))
  temp <- data.frame(CommonName = names(PC1),
                     PC1allGenes = PC1,
                     Study = study)
  temp$Profile <- sapply(as.character(temp$CommonName), function(subj){
    strsplit(subj, "_")[[1]][1]
  })
  temp$Astrocyte_Genes <- meta$Astrocyte_Genes[match(temp$CommonName, meta$CommonName)]
  temp$GabaPV_Genes <- meta$GabaPV_Genes[match(temp$CommonName, meta$CommonName)]
  temp$RIN <- meta$RIN[match(temp$CommonName, meta$CommonName)]
  temp$Profile <- relevel(as.factor(temp$Profile), ref = "Cont")
  temp
}, simplify = FALSE) %>% do.call(rbind, .)
 
PC1allData <- rbind(PC1allData,
                    studyFinal$BA9$Metadata %>% select(CommonName, PC1allGenes, Study, Profile, Astrocyte_Genes, GabaPV_Genes, RIN))

PC1allData$Study <- relevel(PC1allData$Study, ref = "CommonMind")

corFunc <- function(cor1, cor2, xloc, yloc, PC1Data = PC1allData ){
  CorResult <- sapply(levels(PC1Data$Study), function(study){
    data = PC1Data %>% filter(Study == study) %>% droplevels()
    Stat <- cor.test(data[[cor1]], data[[cor2]], method = "spearman")
    paste0("r = ",signif(Stat$estimate, digits = 2)) 
  })
  CorDF <- data.frame(Study = sapply(names(CorResult), function(x){
    strsplit(x, "\\.")[[1]][1]
    }),
    Cor = CorResult,
    x = xloc,
    y = yloc)
  return(CorDF)
}

PC1_AstroCor <- corFunc("PC1allGenes", "Astrocyte_Genes", 0.2, 0)
PC1_PVCor <- corFunc("PC1allGenes", "GabaPV_Genes", 0.2, 0)
PV_AstroCor <- corFunc("Astrocyte_Genes", "GabaPV_Genes", 0.2, 0)
PC1_RINCor <- corFunc("RIN","PC1allGenes", c(9, 8.5),0.8, 
                       PC1Data = PC1allData %>% filter(Study %in% c("CommonMind", "GSE53987")) %>% droplevels )
fsPV_MGP_RINCor <- corFunc("RIN", "GabaPV_Genes", c(9, 8.5),0.15,
                      PC1Data = PC1allData %>% filter(Study %in% c("CommonMind", "GSE53987")) %>% droplevels )


ggplot(PC1allData, aes(Profile, PC1allGenes)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "PC1 score") +
  geom_boxplot(outlier.shape = NA, aes(fill = Profile), alpha = 0.6) +
  scale_fill_manual(values = c("white", "dodgerblue2", "darkolivegreen4"), name = "Group") +
  geom_jitter(size = 0.6, width = 0.2, alpha = 0.4) +
  facet_wrap(~Study, scales = "free_x") 
ggsave("GroupPC1.pdf", path = GeneralResultsPath, width = 8, height = 6, units = "in", dpi=300)

ggplot(PC1allData, aes(PC1allGenes, Astrocyte_Genes)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "PC1 score", y = "Astrocyte MGP") +
  geom_point(size = 0.6, alpha = 0.4) +
  facet_wrap(~Study, scales = "free_x") +
  geom_text(data = PC1_AstroCor,  mapping = aes(x = x, y = y, label = Cor), inherit.aes = FALSE, color = "red")
ggsave("PC1_AstroCor.pdf", path = GeneralResultsPath, width = 8, height = 6, units = "in", dpi=300)

ggplot(PC1allData, aes(PC1allGenes, GabaPV_Genes)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "PC1 score", y = "fsPV MGP") +
  geom_point(size = 0.6, alpha = 0.4) +
  facet_wrap(~Study, scales = "free_x") +
  geom_text(data = PC1_PVCor,  mapping = aes(x = x, y = y, label = Cor), inherit.aes = FALSE, color = "red")
ggsave("PC1_PVCor.pdf", path = GeneralResultsPath, width = 8, height = 6, units = "in", dpi=300)

ggplot(PC1allData %>% filter(Study %in% c("CommonMind", "GSE53987")), aes(RIN, PC1allGenes)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "RIN", y = "PC1 score") +
  geom_point(size = 0.6, alpha = 0.4) +
  facet_wrap(~Study, scales = "free") +
  geom_text(data = PC1_RINCor,  mapping = aes(x = x, y = y, label = Cor), inherit.aes = FALSE, color = "red")
ggsave("PC1_RINCor.pdf", path = GeneralResultsPath, width = 6, height = 3, units = "in", dpi=300)

ggplot(PC1allData %>% filter(Study %in% c("CommonMind", "GSE53987")), aes(RIN, GabaPV_Genes)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "RIN", y = "fsPVMGP") +
  geom_point(size = 0.6, alpha = 0.4) +
  facet_wrap(~Study, scales = "free") +
  geom_text(data = fsPV_MGP_RINCor,  mapping = aes(x = x, y = y, label = Cor), inherit.aes = FALSE, color = "red")
ggsave("fsPV_RINCor.pdf", path = GeneralResultsPath, width = 6, height = 3, units = "in", dpi=300)

ggplot(PC1allData, aes(Astrocyte_Genes, GabaPV_Genes)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Astrocyte MGP", y = "fsPV MGP") +
  geom_point(size = 0.6, alpha = 0.4) +
  facet_wrap(~Study, scales = "free_x") +
  geom_text(data = PV_AstroCor,  mapping = aes(x = x, y = y, label = Cor), inherit.aes = FALSE, color = "red")

save.image(paste0(GeneralResultsPath,"CommonMindSubSample.RData"))
