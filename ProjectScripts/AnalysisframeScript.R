source("SetUp.R")

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

if(grepl("Stanley", name)) {
  subCol = "Study"
} else {
  subCol = "NeuExpRegion"
}

Metadata <- GetMeta(name, subCol = subCol)
Metadata %<>% filter(Profile != "MD")

#Download cell files
DownloadCel(Metadata$Filename, paste0(name, "/data/"))


#read cell files
if(grepl("Stanley", name)){
  for(study in unique(Metadata$Study)){
    print(study)
    regex = paste0(Metadata %>%
                     filter(Study == study) %>% .$Filename, collapse="|")
    CelFiles=list.files(paste0(name, "/data/"),
                        pattern=regex)
    Data <- ReadCell(path = paste0(name,"/data/"),CelFiles=CelFiles,
                     QC=0, platform=Metadata %>% filter(Study == study) %>% .$Platform %>% unique)
    names(Data$aned)[-c(1:3)] <- Metadata$CommonName[pmatch(names(Data$aned),
                                                            Metadata$Filename)] %>% na.omit %>% as.character
    
    assign(paste0("Data_", study), Data)
    rm(Data, study, regex)
  }
  #Create eSet objects
  for(regData in ls(pat="Data_")){
    Data <- eval(as.name(regData))
    study = gsub("Data_", "", regData)
    Data_eSet <- makeSet(data = Data, meta=Metadata %>% filter(Study == study), name=name)
    assign(paste0(study, "_eSet"), Data_eSet)
    rm(Data_eSet, study)
  }
} else if(name == "GSE92538"){
  for(platform in unique(Metadata$Platform)){
    print(platform)
    if(platform == "GPL10526"){
      platform2 = "GPL570"
    } else if (platform == "GPL17027"){
      platform2 = "GPL96"
    }
    regex = paste0(Metadata %>%
                     filter(Platform == platform) %>% .$Filename, collapse="|")
    CelFiles=list.files(paste0(name, "/data/"),
                        pattern=regex)
    Data <- ReadCell(path = paste0(name,"/data/"),CelFiles=CelFiles,
                     QC=0, platform=platform2)
    
    names(Data$aned)[-c(1:3)] <- Metadata$CommonName[pmatch(names(Data$aned),
                                                            Metadata$Filename)] %>% na.omit %>% as.character
    
    assign(paste0("Data_", platform), Data)
    rm(Data, region, regex)
  }
  #Create eSet objects
  for(regData in ls(pat="Data_")){
    Data <- eval(as.name(regData))
    platform = strsplit(regData, "_")[[1]][2]
    Data_eSet <- makeSet(data = Data, meta=Metadata %>% filter(Platform == platform), name=name)
    Study = strsplit(regData, "Data_")[[1]][2]
    assign(paste0(Study, "_eSet"), Data_eSet)
    rm(Data_eSet, region)
  }
  } else {
  for(region in unique(Metadata$NeuExpRegion)){
    print(region)
    regex = paste0(Metadata %>%
                     filter(NeuExpRegion == region) %>% .$Filename, collapse="|")
    CelFiles=list.files(paste0(name, "/data/"),
                        pattern=regex)
    Data <- ReadCell(path = paste0(name,"/data/"),CelFiles=CelFiles,
                     QC=0, platform=Metadata$Platform %>% as.character %>% unique)
    
    names(Data$aned)[-c(1:3)] <- Metadata$CommonName[pmatch(names(Data$aned),
                                                            Metadata$Filename)] %>% na.omit %>% as.character
    
    assign(paste0("Data_", region), Data)
    rm(Data, region, regex)
  }
  #Create eSet objects
  for(regData in ls(pat="Data_")){
    Data <- eval(as.name(regData))
    region = strsplit(regData, "_")[[1]][2]
    Data_eSet <- makeSet(data = Data, meta=Metadata %>% filter(NeuExpRegion == region), name=name)
    Study = strsplit(regData, "Data_")[[1]][2]
    assign(paste0(Study, "_eSet"), Data_eSet)
    rm(Data_eSet, region)
  }
}

studyFinal <- sapply(ls(pat="_eSet"), function(x) {
  study = strsplit(x, "_")[[1]][1]
  eSet = eval(as.name(x))
  PreProcces(eSet=eSet, study=study)
}, simplify=FALSE)

names(studyFinal) <- sapply(ls(pat="_eSet"), function(x) strsplit(x, "_eS")[[1]][1]) 

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

rm(list=ls(pat="exp"))
#Create gender HeatMaps
datas <- datasGenerate(c("XIST", "KDM5D", "RPS4Y1|RPS4Y2", "RPS4Y1"))


HeatMapGen(datas=datas, Meta=Metadata[!duplicated(Metadata$CommonName),], path = resultsPath,  save=1)
HeatMapGen2(datas=datas, Meta=Metadata[!duplicated(Metadata$CommonName),], path = resultsPath, missmatched=missmatched, save=1)

#Estimating cell type proportions
source(paste0(GenScriptPath, "Cell_type_PCA.R"))

studyFinal$GPL10526$Metadata %<>% filter(pH >= 6.5)
studyFinal$GPL17027$Metadata %<>% filter(pH >= 6.5)

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


