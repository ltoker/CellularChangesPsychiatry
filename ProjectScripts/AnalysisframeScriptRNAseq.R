source("SetUp.R")
packageF("biomaRt")

Count2CPM <- function(countData){
  apply(countData, 2, function(smp){
    TotalCount = sum(smp)
    (10^6)*smp/TotalCount
  })
}


ensembl <- useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl")

geneNames <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), mart = ensembl)

MitoGenes <- geneNames[grepl("MT-", geneNames$hgnc_symbol),]

#Reading tables
######################################################################################################
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


Metadata <- GetMeta(name, subCol = subCol)
Metadata %<>% filter(Profile != "MD")

#Remove brain regions not represented in NeuroExpresso
Metadata <- Metadata[!is.na(Metadata$NeuExpRegion),] %>% droplevels()
Metadata %<>% mutate(Filename = Series_sample_id,
                     Study = name)

#Get the count matrix and filter mitochondrial genes

countMatrix <- read.csv(paste0(name, "/data/countMatrix.genes"), header=TRUE, sep = "\t", quote = "")
CountSum <- apply(countMatrix[grepl("GSM", colnames(countMatrix))], 2, sum)

MitoCountSum <- apply(countMatrix[countMatrix$genes %in% MitoGenes$ensembl_gene_id,
                                  grepl("GSM", colnames(countMatrix))], 2, sum)

MitoCountFiltered <- countMatrix[!countMatrix$genes %in% MitoGenes$ensembl_gene_id,]

CPMmatrix <- Count2CPM(MitoCountFiltered[-1])

#log2 transformation
CPMmatrix <- apply(CPMmatrix, c(1,2), function(x) {log2(x+1)})

GeneSymbol <- geneNames$hgnc_symbol[match(MitoCountFiltered$genes, geneNames$ensembl_gene_id)]

ExpData <- data.frame(GeneSymbol = GeneSymbol,
                      Probe = GeneSymbol,
                      ensemblID = MitoCountFiltered$genes)

ExpData <- cbind(ExpData, CPMmatrix) 

RegionData <- sapply(levels(Metadata$OrgRegion), function(region){
  subMeta = Metadata %>% filter(OrgRegion == region)
  subExp = ExpData %>% select_(.dots = c("GeneSymbol", "Probe", "ensemblID", as.character(subMeta$Series_sample_id)))
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

#Print MGP plots
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

save.image(paste0(GeneralResultsPath, name, ".RData"))
save(studyFinal, file = paste0(GeneralResultsPath, "studyFinal", name, ".rda"))
save(PCA_results, file = paste0(GeneralResultsPath, "PCAresults", name, ".Rda"))



