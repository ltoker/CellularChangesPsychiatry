GenScriptPath = "/home/ltoker/Rscripts/"
ProjScriptPath = "ProjectScripts/"
if(!"GeneralResults" %in% list.dirs(full.names = FALSE)){
  dir.create("GeneralResults")
}
GeneralResultsPath = 'GeneralResults/'

source(paste0(GenScriptPath,"general_functions.R"))
packageF("sva")
packageF("gplots")
packageF("scales")
select = dplyr::select
filter = dplyr::filter
mutate = dplyr::mutate
source("projectFunc.R")

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

Metadata %<>% filter(brain.region %in% c("frontal cortex", "white matter"), Platform == "GPL5175") %>% droplevels()

Metadata$NeuExpRegion <- "Cortex"
Metadata$NeuExpRegion <- factor(Metadata$NeuExpRegion)
Metadata$CommonName <- Metadata$Filename

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
} else {
  for(region in unique(Metadata$NeuExpRegion)){
    print(region)
    regex = paste0(Metadata %>%
                     filter(NeuExpRegion == region) %>% .$Filename, collapse="|")
    CelFiles=list.files(paste0(name, "/data/"),
                        pattern=regex)
    Data <- ReadCell(path = paste0(name,"/data/"),CelFiles=CelFiles,
                     QC=0, platform=Metadata$Platform %>% unique)
    
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
datas <- datasGenerate(c("XIST", "KDM5D", "RPS4Y1|RPS4Y2|RPS4Y1"))


HeatMapGen(datas=datas, Meta=Metadata[!duplicated(Metadata$CommonName),], path = resultsPath,  save=1)
HeatMapGen2(datas=datas, Meta=Metadata[!duplicated(Metadata$CommonName),], path = resultsPath, missmatched=missmatched, save=1)

#Estimating cell type proportions

source(paste0(ProjScriptPath, "Cell_type_PCA.R"))

#PCA analysis without the missmatched samples
PCA_results <- sapply(names(studyFinal), function(x){
  region = studyFinal[[x]]$Metadata$NeuExpRegion %>% as.character %>% unique
  browser()
  CellType_genes <- GetMarkers(region)
  
  #Exclude GabaPV genes which are not neuron specific in human (Darmanis) data
  if(region == "Cortex"){
    CellType_genes$GabaPV_Genes <- CellType_genes$GabaPV_Genes[!CellType_genes$GabaPV_Genes %in% c("WIF1", "TMEM132C", "BTN2A2")]
  }
  
  aned_high <- studyFinal[[x]]$aned_high
  aned_high <- aned_high[,!names(aned_high) %in% missmatched[[x]]]
  
  results <- PCA_genes_All_based(dataset_id=paste0(name,x),
                                 dataset=aned_high,
                                 CellType_genes=CellType_genes,
                                 NoiseThershold = studyFinal[[x]]$NoiseThreshold,
                                 contName = "GSM")
  
  return(results)
}, simplify = F)

#Add estimation to Metadata 
for(study in names(studyFinal)){
  studyFinal[[study]]$Metadata %<>% select(-matches("_Genes"))
  estimates <- lapply(PCA_results[[study]]$modified, function(cells){
    temp <- cells$x[,1]
    names(temp) <- rownames(cells$x)
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


temp <- studyFinal$Cortex$Metadata %>% select(matches("Filename|_Genes")) %>% filter(!is.na(.$Astrocyte_Genes))
cellRM <- apply(temp,2, function(x) sum(is.na(x)))
temp <- temp[!cellRM > 0]
rownames(temp) <- temp$Filename
PCAcells <- prcomp(temp[-1])
temp2 <- data.frame(GSM = rownames(PCAcells$x),
                                PC1 = PCAcells$x[,1],
                                region =studyFinal$Cortex$Metadata$brain.region[match(rownames(PCAcells$x),
                                                                                               studyFinal$Cortex$Metadata$Filename)])
temp2 %<>% arrange(PC1)      
MMsamples <- read.csv("GSE46706/RegionMislabeledSamples.txt", header = FALSE) %>% unlist() %>% as.character()
MMsamples <- c(MMsamples,
               temp2 %>% filter(region == "frontal cortex", PC1 > 0) %>% .$GSM %>% as.character(),
               temp2 %>% filter(region == "white matter", PC1 < 0) %>% .$GSM %>% as.character())

ggplot(temp2, aes(x = region, y = PC1, color = region)) + geom_jitter()

ggplot(studyFinal$Cortex$Metadata, aes(Age, Oligo_Genes, color = brain.region)) + geom_point()
ggplot(studyFinal$Cortex$Metadata, aes(GabaPV_Genes, Oligo_Genes, color = brain.region)) + geom_point()


ggplot(studyFinal$Cortex$Metadata %>% filter(!Filename %in% MMsamples), aes(Age, Oligo_Genes, color = brain.region)) + geom_point()

save.image(paste0(GeneralResultsPath, name, ".RData"))
save(studyFinal, file = paste0(GeneralResultsPath, "studyFinal", name, ".rda"))