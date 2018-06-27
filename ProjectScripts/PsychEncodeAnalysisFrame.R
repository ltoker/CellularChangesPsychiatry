source("SetUp.R")
packageF("biomaRt")
packageF("parallel")
packageF("pheatmap")
name = "PsychEncode"

Count2CPM <- function(countData){
  apply(countData, 2, function(smp){
    TotalCount = sum(smp)
    (10^6)*smp/TotalCount
  })
}

####### Preparing analysis directories #########
path = paste0(getwd(), "/")

if(length(ls(pat = "^name$")) == 0){
  print("no study name")
  browser()
}

if(!name %in% list.dirs(full.names = FALSE)){
  dir.create(name)
  dir.create(paste0(name, "/data"))
}

if(!paste0(name, "/data") %in% list.dirs(full.names = FALSE)){
  dir.create(paste0(name, "/data"))
}

resultsPath = paste0(path, name, "/")

subCol = "NeuExpRegion"

########## BioMart #########
ensembl <- useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl")

geneNames <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), mart = ensembl)
MitoGenes <- geneNames[grepl("MT-", geneNames$hgnc_symbol),]
###### Getting the metadata#########
Meta1 <- read.table("PEC_BrainGVEX_ATACseq_metadata.tsv", header = T, sep = "\t") %>% filter(!is.na(PMI)) %>% droplevels()
Meta1Stanley <- Meta1 %>% filter(BrainBank != "BSHRI") %>% droplevels() 

write.table(Meta1Stanley$RNASeq.BID, "PsychEncodeIDsStanley.tsv", sep = "\t", row.names = FALSE, col.names = FALSE)

Meta2 <- read.table("UIC-UChicago-U01MH103340_BrainGVEX_ClinicalMetadata_AgeCensored_August2016Release.tsv", header = T, sep = "\t")

Meta1Stanley$pH <- Meta2$pH[match(Meta1Stanley$RNASeq.BID, Meta2$Individual_ID..RNAseq.library.BID..Primary.ID.)]

RNAseqMeta <- read.table("UIC-UChicago-U01MH103340_BrainGVEX_RNAseqMetadata_August2016Release.txt", header = T, sep = "\t")
RNAseqMeta$Individual_ID..RNAseq.Library.BID..Synapse. <- sapply(RNAseqMeta$File_Name..on.Synapse..if.available, function(x){
  strsplit(as.character(x), "_")[[1]][7]
})

RNAseqMeta <- merge(RNAseqMeta, Meta2 %>% select(Individual_ID..RNAseq.library.BID..Primary.ID.,BrainBank, Diagnosis), by.x = "Individual_ID..RNAseq.Library.BID..Synapse.",
                    by.y = "Individual_ID..RNAseq.library.BID..Primary.ID.") 

RNAseqMetaStanley <- RNAseqMeta %>%
  filter(Individual_ID..RNAseq.Library.BID..Synapse. %in% as.character(Meta1Stanley$RNASeq.BID)) %>%
  droplevels()

RNAseqMetaStanley %<>% mutate_at(grep("Reads", names(RNAseqMetaStanley), value = T), as.character)
RNAseqMetaStanley %<>% mutate_at(grep("Reads", names(RNAseqMetaStanley), value = T), as.numeric)
RNAseqMetaStanley$rRNARate <- sapply(RNAseqMetaStanley$rRNARate, function(x){
  x <- gsub("%", "", as.character(x))
  as.numeric(x)
})

Meta1Stanley$SequencingPlatform <- RNAseqMetaStanley$SequencingPlatform[match(Meta1Stanley$RNASeq.BID,
                                                                              RNAseqMetaStanley$Individual_ID..RNAseq.Library.BID..Synapse.)]
Meta1Stanley$rRNArate <- RNAseqMetaStanley$rRNARate[match(Meta1Stanley$RNASeq.BID,
                                                          RNAseqMetaStanley$Individual_ID..RNAseq.Library.BID..Synapse.)]
Meta1Stanley$ERCC <- RNAseqMetaStanley$ERCC_Added[match(Meta1Stanley$RNASeq.BID,
                                                          RNAseqMetaStanley$Individual_ID..RNAseq.Library.BID..Synapse.)]


Metadata <- Meta1Stanley %>% select(-ATACSeq.BID)
Metadata %<>% mutate(NeuExpRegion = "Cortex",
                     OrgRegion = as.factor("Cortex"),
                     Series_sample_id  = RNASeq.BID,
                     Filename = RNASeq.BID)

Metadata$Series_sample_id <- sapply(Metadata$Series_sample_id, function(x){
  gsub("PEC_BrainGVEX_UIC.UChicago_FC_mRNA_HiSeq2.00_", "", x)
})

Metadata$Profile <- sapply(as.character(Metadata$Diagnosis), function(x){
  if(x == "BP"){
    "BP"
  } else if (x == "SCZ"){
    "SCZ"
  } else if(x == "Control"){
    "Cont"
  } else {
    NA
  }
}) %>% factor(levels = c("Cont", "BP", "SCZ"))
Metadata$ScanDateOrg <- RNAseqMetaStanley$LibraryBatch[match(Metadata$RNASeq.BID, RNAseqMetaStanley$Individual_ID..RNAseq.Library.BID..Synapse.)]
DateMatch <- data.frame(OrgDate = unique(Metadata$ScanDateOrg))
DateMatch$ConvertDate <- sapply(as.character(DateMatch$OrgDate), function(Date){
  if(grepl("/", Date)){
    Date
  } else {
    Year = strsplit(Date, "-")[[1]][3]
    Month = strsplit(Date, "-")[[1]][1]
    Month = rev(strsplit(Month, "0")[[1]])[1]
    Day = strsplit(Date, "-")[[1]][2]
    Day = gsub("0","", Day)
    paste(Month, Day, Year, sep = "/")
  }
})

Metadata$ScanDate <- DateMatch$ConvertDate[match(Metadata$ScanDateOrg, DateMatch$OrgDate)]


SampleNum <- Metadata %>% group_by(Profile) %>% summarise(n = n()) %>% data.frame()

Metadata %<>% arrange(Profile) %<>% droplevels()

Metadata$CommonName <- apply(SampleNum, 1, function(x){
  paste0(x[1], "_", 1:x[2])
}) %>% unlist

Metadata <- merge(Metadata, Meta2 %>%
                    select(Individual_ID..RNAseq.library.BID..Primary.ID., AgeOnset, DurationIllness, CauseDeath, Height, Weight, EtohResults),
                  by.x = "RNASeq.BID", by.y = "Individual_ID..RNAseq.library.BID..Primary.ID.", all = FALSE)
Metadata %<>% mutate(Weight = as.numeric(as.character(Weight)),
                     Height = as.numeric(as.character(Height)))
Metadata %<>% mutate(BMI = 703*Weight/(Height)^2)

####### Getting expression data ###############
#Get normalized count matrix
countMatrix <- read.table("PsychEncode/data/countMatrix.genes", header = T, sep = "\t", quote = "")
names(countMatrix) <- sapply(names(countMatrix), function(x) gsub("PEC_BrainGVEX_UIC.UChicago_FC_mRNA_HiSeq2.00_", "", x))
names(countMatrix) <- sapply(names(countMatrix), function(x) gsub("\\.", "-", x))

corSample <- cor(countMatrix[-1])
diag(corSample) <- NA
MedianCor <- apply(corSample, 1, function(x) median(x, na.rm=T)) %>% sort
MedianCorLow <- MedianCor[MedianCor < 0.95] %>% names(.) %>% sapply(function(x) strsplit(x, "_")[[1]][1]) 

corSample <- corSample[,match(names(MedianCor), colnames(corSample))]
par(mar = c(10, 4.1, 4.1,2.1))
boxplot(corSample, outline = FALSE, las = 2)
pheatmap(as.matrix(corSample))

ReadTotal <- apply(countMatrix[-1], 2, sum)
ReadMT <- apply(countMatrix[countMatrix$genes %in% MitoGenes$ensembl_gene_id, -1], 2, sum)

ReadsDF <- data.frame(All = ReadTotal, MT = ReadMT)
ReadsDF %<>% mutate(Ratio = MT/All,
                    SampleID = names(countMatrix)[-1])
ReadsDF$SubjectID <- sapply(ReadsDF$SampleID,function(x) strsplit(x, "_")[[1]][1]) %>% as.character()
ReadsDF$MedianCor <- MedianCor[match(ReadsDF$SampleID, names(MedianCor))]
ReadsDF$CommonName <- Metadata$CommonName[match(ReadsDF$SubjectID, Metadata$RNASeq.BID)]
ReadsDF$ScanDate <- Metadata$ScanDate[match(ReadsDF$SubjectID, Metadata$RNASeq.BID)]

ReadsDF$MappedReadSynapse <- Metadata$num_unique_reads_mapped[match(ReadsDF$SubjectID, Metadata$RNASeq.BID)]
ReadsDF$MTratioSynapse <- Metadata$Mitochondrial.contamination[match(ReadsDF$SubjectID, Metadata$RNASeq.BID)]

MitoCountFiltered <- countMatrix[!countMatrix$genes %in% MitoGenes$ensembl_gene_id,]

CPMmatrix <- Count2CPM(MitoCountFiltered[-1])

#log2 transformation
CPMmatrix <- apply(CPMmatrix, c(1,2), function(x) {log2(x+1)})
CPMmatrixCor <- cor(CPMmatrix)
diag(CPMmatrixCor) <- NA
pheatmap(CPMmatrixCor)

#For samples ran in replicates - keep the replicate with the highest read count
RepSubj <- table(ReadsDF$SubjectID)[table(ReadsDF$SubjectID) > 1]
SamplRM <- vector()
for(subj in names(RepSubj)){
  temp <- ReadsDF %>% filter(SubjectID == subj) %>% arrange(desc(All))
  temp <- temp$SampleID[duplicated(temp$SubjectID)] %>% as.character()
  SamplRM <- c(SamplRM, temp)
}

CPMmatrix <- CPMmatrix[,!colnames(CPMmatrix) %in% SamplRM]

#Remove genes which are all zeros (no variance)
VarGenes <- apply(CPMmatrix, 1, sd)
CPMmatrix <- CPMmatrix[VarGenes > 0.01,]

CPMmatrixCor <- cor(CPMmatrix)
diag(CPMmatrixCor) <- NA
pheatmap(CPMmatrixCor)

PCAallGenes <- prcomp(t(CPMmatrix), scale = T)
rownames(PCAallGenes$x) <- sapply(rownames(PCAallGenes$x), function(x) strsplit(x, "_")[[1]][1])
Metadata$PC1allGenes <- PCAallGenes$x[match(Metadata$RNASeq.BID, rownames(PCAallGenes$x)),1]

GeneSymbol <- geneNames$hgnc_symbol[match(MitoCountFiltered$genes, geneNames$ensembl_gene_id)]

ExpData <- data.frame(GeneSymbol = GeneSymbol,
                      Probe = MitoCountFiltered$genes,
                      ensemblID = MitoCountFiltered$genes)

ExpData <- cbind(ExpData[VarGenes > 0.01,], CPMmatrix)
names(ExpData) <- sapply(names(ExpData), function(x){
  strsplit(x, "_")[[1]][1]
})

ExpData <- ExpData[!is.na(ExpData$GeneSymbol),]
ExpData <- ExpData[!ExpData$GeneSymbol == "",]


RegionData <- sapply(levels(Metadata$OrgRegion), function(region){
  subMeta = Metadata %>% filter(OrgRegion == region)
  subExp = ExpData[,c("GeneSymbol", "Probe", "ensemblID", as.character(subMeta$Series_sample_id))]
  names(subExp)[-c(1:3)] <- Metadata$CommonName[match(names(subExp)[-c(1:3)], Metadata$Series_sample_id)] %>% as.character()
  list(Metadata = subMeta, aned = subExp)
}, simplify = FALSE)

studyFinal <- lapply(RegionData, function(region) {
  Metadata = region$Metadata
  aned = region$aned
  source(paste0(GenScriptPath, "pre-proccess.R"), local=T)
  output <- list(aned_high, aned_good, aned_low, MaxNoise,Z_scores_High,
                 exclude_samples_low, exclude_samples_high, exclude_samples, Metadata)
  names(output) <- c("aned_high", "aned_good", "aned_low", "NoiseThreshold", "Z_scores_High",
                     "exclude_samples_low", "exclude_samples_high", "exclude_samples", "Metadata")
  output
})

studyFinal <- lapply(studyFinal, function(regData) {
  regData$Metadata <- GeneSex(regData$aned_good, Metadata=regData$Metadata)
  return(regData)
})

PCAallGenesCombat <- prcomp(t(studyFinal$Cortex$aned_high[-c(1:3)]), scale = TRUE)
studyFinal$Cortex$Metadata$PC1allGenesCombat <- PCAallGenesCombat$x[match(studyFinal$Cortex$Metadata$CommonName, rownames(PCAallGenesCombat$x)),1]

missmatched <- lapply(studyFinal, function(x){
  meta <- x$Metadata
  names(meta) <- tolower(names(meta))
  meta$sex <- sapply(meta$sex, function(sex){
    if(grepl("^male|^man|^m$", tolower(sex))){
      "M"
    } else if(grepl("female|^wom|w|^f$", tolower(sex))){
      "F"
    }
  })
  meta$commonname[meta$sex != meta$biogender]
})

#Create gender HeatMaps
datas <- datasGenerate(c("XIST", "KDM5D", "RPS4Y1"))


HeatMapGen(datas=datas, Meta=Metadata[!duplicated(Metadata$CommonName),], path = resultsPath,  save=1)
HeatMapGen2(datas=datas, Meta=Metadata[!duplicated(Metadata$CommonName),], path = resultsPath, missmatched=missmatched, save=1)

#Estimating cell type proportions
source(paste0(GenScriptPath, "Cell_type_PCA.R"))

#PCA analysis without the missmatched samples
PCA_results <- sapply(names(studyFinal), function(x){
  region = studyFinal[[x]]$Metadata$NeuExpRegion %>% unique
  CellType_genes <- GetMarkers(region)
  
  #Exclude GabaPV genes which are not neuron specific in human (Darmanis) data
  if(region == "Cortex"){
    CellType_genes$GabaPV_Genes <- CellType_genes$GabaPV_Genes[!CellType_genes$GabaPV_Genes %in% c("WIF1", "TMEM132C", "BTN2A2")]
  }
  
  aned_high <- studyFinal[[x]]$aned_high
  aned_high <- aned_high[,!names(aned_high) %in% missmatched[[x]]]
  
  #bootstrap with replacement the samples in each group to ensure equal number of samples/group (90% of the samples in the smaller group)
  groups <- sapply(grep("_", names(aned_high), value=T),
                   function(x) gsub("_.*", "", x)) %>% table
  MinGrp <- round(0.9*min(groups))
  AllSamples <- sapply(names(groups), function(grp){
    grep(grp, names(aned_high), value = T)
  }, simplify = FALSE)
  results <- list()
  for(i in c(1:100)){
    BootSamples <- lapply(AllSamples, function(grp){
      grp[sample(1:length(grp), MinGrp,replace = FALSE)]
    }) %>% unlist %>% as.character
    
    aned_highSub <- aned_high %>% select_(.dots = c("GeneSymbol", BootSamples))
    results[[i]] <- PCA_genes_All_based(dataset_id=x,
                                        dataset=aned_highSub,
                                        CellType_genes=CellType_genes,
                                        NoiseThershold = studyFinal[[x]]$NoiseThreshold)
  }
  
  return(results)
}, simplify = F)

PCA_resultsMean <- lapply(PCA_results, function(region){
  MeanPCA <- sapply(names(region[[1]]$modified), function(celltype){
    temp <- data.frame(CommonName = names(region[[1]]$modified[[celltype]]$x[,1]),
                       Rot = region[[1]]$modified[[celltype]]$x[,1])
    for(i in 2:length(region)){
      temp <- merge(temp, data.frame(CommonName = names(region[[i]]$modified[[celltype]]$x[,1]),
                                     Rot = region[[i]]$modified[[celltype]]$x[,1]), by = "CommonName", all = TRUE)
    }
    names(temp)[2:ncol(temp)] <- paste0("Rot", c(1:c(ncol(temp)-1)))
    temp$MeanRot <- rowMeans(temp[-1], na.rm = T)
    temp
  }, simplify=FALSE)
})

#Add estimation to Metadata 
for(study in names(studyFinal)){
  studyFinal[[study]]$Metadata %<>% select(-matches("_Genes"))
  estimates <- lapply(PCA_resultsMean[[study]], function(cells){
    temp <- cells$MeanRot
    names(temp) <- cells$CommonName
    temp
  }) %>% List2df
  studyFinal[[study]]$Metadata <- merge(studyFinal[[study]]$Metadata,
                                        estimates,
                                        by.x="CommonName",
                                        by.y="row.names",
                                        all.x=TRUE)
  studyFinal[[study]]$Metadata$Profile <- as.factor(studyFinal[[study]]$Metadata$Profile)
  studyFinal[[study]]$Metadata$Profile <- relevel(studyFinal[[study]]$Metadata$Profile, ref="Cont")
}

#Plot results
sapply(names(studyFinal), function(stdName){
  if(!stdName %in% list.dirs(resultsPath,full.names = FALSE)){
    dir.create(paste0(resultsPath, stdName))
  }
  meta <- studyFinal[[stdName]]$Metadata
  meta <- meta[apply(meta, 2, function(x) sum(!is.na(x))) > 0]
  names(meta) <- sapply(names(meta), function(x) gsub(" ", "_", x))
  sapply(grep("_Genes", names(meta), value = TRUE), function(mgp){
    temp <- PlotPCggplot(data=meta, CellVar = mgp,
                         name = paste(name, stdName, sep = "-"), txtSize = 16, pValSize = )
    ggsave(paste0(resultsPath, stdName, "/", gsub("Genes", "MGP", mgp), ".pdf"),
           plot = temp, width = 12, height = 8, units = "in",
           dpi=300)
  })
})

#Get variance explained
VarExplainedPV <- mclapply(PCA_results$Cortex, function(x){
  x$All$GabaPV_Genes %>% summary() %>% .$importance %>% data.frame() %>% .[2,1:3] %>% unlist
  }, mc.cores = 0.5*detectCores()) %>% do.call(rbind, .) %>% data.frame()

VarExplainedPV$Rot <- paste0("Rot_", 1:100)

VarExplainedPV$GabaPVshift <- mclapply(PCA_results$Cortex, function(Rot){
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

VarExplainedAstro <- mclapply(PCA_results$Cortex, function(x){
  x$All$Astrocyte_Genes %>% summary() %>% .$importance %>% data.frame() %>% .[2,1:3] %>% unlist
}, mc.cores = 0.5*detectCores()) %>% do.call(rbind, .) %>% data.frame()

VarExplainedAstro$Rot <- paste0("Rot_", 1:100)

VarExplainedAstro$Astroshift <- mclapply(PCA_results$Cortex, function(Rot){
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

lmAstro <- lm(Astrocyte_Genes~Profile + pH + PMI + RIN + AgeDeath, data = studyFinal$Cortex$Metadata)
lmPV <- lm(GabaPV_Genes~Profile + pH + PMI + RIN + AgeDeath, data = studyFinal$Cortex$Metadata)

MetaMGP <- studyFinal$Cortex$Metadata %>% select(CommonName, Profile, BrainBank, Astrocyte_Genes, GabaPV_Genes)
MetaMGP$AstroCovAdj <- ModelAdj(lmAstro, adj = data.frame(effect = c("pH", "PMI", "RIN", "AgeDeath"), adjValue= ""))
MetaMGP$PVCovAdj <- ModelAdj(lmPV, adj = data.frame(effect = c("pH", "PMI", "RIN", "AgeDeath"), adjValue= ""))
MetaMGPmelt <- melt(MetaMGP, id.vars = c("CommonName", "Profile", "BrainBank"), variable.name = "CellType", value.name = "MGP")
MetaMGPmelt$CovarAdjust <- sapply(as.character(MetaMGPmelt$CellType), function(x){
  if(grepl("Adj", x)){
    "Yes"
  } else {
    "No"
  }
})

ContMedian <- MetaMGPmelt %>% filter(CovarAdjust == "Yes", Profile == "Cont") %>% group_by(CellType) %>%
  summarise(median = median(MGP, na.rm = T)) %>% data.frame

StatAstro <- lmAstro %>% summary %>% .$coefficients %>% data.frame
StatPV <- lmPV %>% summary %>% .$coefficients %>% data.frame

StatDF <- data.frame(Profile = rep(c("BP", "SCZ"), 2),
                     CellType = c(rep("AstroCovAdj",2), rep("PVCovAdj", 2)),
                     Pval = paste0("p < ", c(StatAstro[c("ProfileBP", "ProfileSC"),]$P,
                              StatPV[c("ProfileBP", "ProfileSC"),]$P) %>% scientific(digits = 2)))
StatDF %<>% mutate(x = rep(c(2,3), 2),
                   y = 1.03)

ggplot(MetaMGPmelt %>% filter(CovarAdjust == "Yes"), aes(Profile, MGP)) +
  theme_classic(base_size = 16) +
  labs(x = "", y = "MGP (covariate adjusteed)") + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(size = rel(1.2)),
        axis.text.y = element_text(size = rel(0.9))) +
  ylim(0, 1.13) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size=2)+
  facet_wrap(~CellType, scale = "free") +
  geom_text(data = StatDF, aes(x, y, label = Pval), size=4, inherit.aes = FALSE, color = "red") +
  geom_hline(data = ContMedian, aes(yintercept = median), color = "red", linetype="dashed")
ggsave("GeneralResults/MGPresultsPsychEncode.pdf", width = 12, height = 5, units = "in",
       dpi=300)

ggplot(MetaMGPmelt %>% filter(CovarAdjust == "Yes"), aes(Profile, MGP)) +
  theme_classic(base_size = 16) +
  labs(x = "", y = "MGP (covariate adjusteed)") + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(size = rel(1.2)),
        axis.text.y = element_text(size = rel(0.9))) +
  ylim(0, 1.13) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size=2)+
  facet_wrap(~CellType + BrainBank, scale = "free") 

save.image(paste0(GeneralResultsPath, name, ".RData"))
save(studyFinal, file = paste0(GeneralResultsPath, "studyFinal", name, ".rda"))
save(PCA_results, file = paste0(GeneralResultsPath, "PCAresults", name, ".Rda"))


