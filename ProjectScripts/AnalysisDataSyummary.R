GenScriptPath = "/home/ltoker/Rscripts/"
ProjScriptPath = "ProjectScripts/"
ResultsPath = "GeneralResults/"
source(paste0(GenScriptPath,"general_functions.R"))
source(paste0(ProjScriptPath, "projectFunc.R"))
packageF("sva")
packageF("gplots")
packageF("scales")


StanleySamples <- read.table("MetadataStanley.csv", header = T, sep = "\t")
StanleyStudies <- read.table("StanleyStudies.txt", header = TRUE, sep = "\t")
StudySamplesAll <- c("Dataset", "BrainRegion", "BrainRegion(Org)",  "Platform", "NoiseThreshold", "Samples(Org)", "SingleBatchExcluded",
                               "OutlierExcluded", "MisannotatedExcluded", "Samples(Final)")

for(DataSet in list.files(path = ResultsPath, pattern = "studyFinal")){
  load(paste0(ResultsPath, DataSet))
  DataSet = gsub("studyFinal|\\.rda", "", DataSet)
  if(!grepl("Stanley", DataSet)){
    MetaOrg <- read.table(paste(DataSet, "Metadata.tsv", sep = "/"), header = TRUE, sep = "\t")
    names(MetaOrg)[grep("disease|Diagnosis|profile|phenotype", names(MetaOrg), ignore.case=TRUE)] <- "Profile"
    names(MetaOrg)[grep("(tissue$|region|brain.?region)", names(MetaOrg),ignore.case = TRUE)] <- "Region"
    MetaOrg$NeuExpRegion <- sapply(MetaOrg$Region, function(x) {
      if (grepl("cerebe",tolower(x))){
        "Cerebellum"
      } else if (grepl("cortex|pfc|dlpfc|frontalba|^ba|gyrus",tolower(x))){
        "Cortex"
      } else if (grepl("hippocampus|hip|hpc",tolower(x))){
        "Hippocampus"
      } else if (grepl("thalamus",tolower(x))) {
        "Thalamus"
      } else if (grepl("str",tolower(x))){
        "Striatum"
      } else if (grepl("putamen",tolower(x))){
        "Putamen"
      } else if (grepl("nigra",tolower(x))){
        "SubstantiaNigra"
      } else if (grepl("brain",tolower(x))){
        "Cortex"
      } else {
        NA
      }
    }, simplify=TRUE)
    
    MetaOrg <- MetaOrg[!grepl("MD|depres", MetaOrg$Profile),] %>% filter(!is.na(.$Profile))
  }
  for(study in names(studyFinal)){
    Region = studyFinal[[study]]$Metadata$NeuExpRegion %>% unique  %>% as.character()
    Platform = studyFinal[[study]]$Metadata$Platform %>% unique %>% as.character()
    Outlier <- length(studyFinal[[study]]$exclude_samples) 
    Noise <- studyFinal[[study]]$NoiseThreshold %>% round(digits = 2) %>% as.numeric
    MM <- sum(is.na(studyFinal[[study]]$Metadata$Oligo_Genes))
    if(grepl("Stanley", DataSet)){
      MetaName <- paste0(study, ".csv")
      MetaOrg <- read.table(paste(DataSet, MetaName, sep = "/"), header = T, sep = "\t")
      MDsamples <- StanleySamples %>% filter(Profile... == "Depression" ) %>% .$Stanley.ID... %>% as.character()
      MetaOrg %<>% filter(!MetaOrg$Stanley.ID... %in% MDsamples)
      OrgSamplNum <- nrow(MetaOrg)
      invest = gsub("Study.[1-9]?", "", study)
      RegionOrg = StanleyStudies %>% filter(Investigator == invest) %>% .$Region %>% as.character
    } else {
      MetaOrg2 <- MetaOrg[grepl(Region, MetaOrg$NeuExpRegion),]
      OrgSamplNum <- nrow(MetaOrg2)
      RegionOrg = MetaOrg2$Region %>% unique  %>% as.character()
    }
    
    Batch <- OrgSamplNum - (Outlier + nrow(studyFinal[[study]]$Metadata))
    
    FinalSamp = group_by(studyFinal[[study]]$Metadata %>% filter(!is.na(.$Oligo_Genes)), Profile) %>% summarise(n = n()) %>%
      data.frame %>%
      apply(1, function(x) paste(x[2], x[1])) %>%  paste0(collapse = ", ")
    FinalSamp = paste0(nrow(studyFinal[[study]]$Metadata %>% filter(!is.na(.$Oligo_Genes))), " (", FinalSamp, ")")
    if(grepl("Stanley", DataSet)){
      DataSet2 = study
    } else {
      DataSet2 = DataSet
    }
    studySamples <- c(DataSet2, Region, RegionOrg,  Platform,  Noise, OrgSamplNum, Batch, Outlier,
                               MM, FinalSamp)
    assign(paste0("studyFinal", DataSet, study), studyFinal[[study]])
    assign(paste0("StudySamples", DataSet, study), studySamples)
    StudySamplesAll <- rbind(StudySamplesAll, studySamples)
  }
}

StudySamplesAllDF <- data.frame(StudySamplesAll[-1,])
names(StudySamplesAllDF) <- StudySamplesAll[1,]
write.table(StudySamplesAllDF, "GeneralResults/StudySamples.tsv", row.names = FALSE, col.names = TRUE, sep = "\t")
