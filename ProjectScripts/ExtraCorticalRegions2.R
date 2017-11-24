GenScriptPath = "/home/ltoker/Rscripts/"
ProjScriptPath = "ProjectScripts/"
source(paste0(GenScriptPath,"general_functions.R"))
source(paste0(ProjScriptPath, "projectFunc.R"))
packageF("sva")
packageF("gplots")
packageF("scales")
source(paste0(GenScriptPath, "Cell_type_PCA.R"))

GetStar <-  function(pvals){
  Star = sapply(pvals, function(p){
    if(is.na(p) | p > 0.05){
      ""
    } else if(p < 0.05 & p > 0.01){
      "*"
    } else if (p < 0.01 & p > 0.001){
      "**"
    } else if ( p < 0.001){
      "***"
    }
  })
  return(Star)
}


MetadataAllStanley <- read.table("MetadataStanley.csv", sep="\t", header=T)
names(MetadataAllStanley) <- gsub("\\.", "", names(MetadataAllStanley))
MetadataAllStanley$Profile <- sapply(MetadataAllStanley$Profile, function(x) {
  if(grepl("bipolar", tolower(x))){
    "BP"
  } else if (grepl("control|unaffected", tolower(x))){
    "Cont"
  } else if (grepl("schizo", tolower(x))){
    "SCZ"
  } else if (grepl("depres", tolower(x))){
    "MD"
  } else {
    NA
  }
}, simplify=TRUE)

MetadataAllStanley$CharVec2 <- apply(MetadataAllStanley %>%
                                       select(Sex, BrainPH, PMIh, Profile), 1, function(x) paste0(x, collapse="_"))

StanleyGSE35978 <- new.env()
load("GeneralResults/studyFinalGSE35978.rda", envir = StanleyGSE35978)
studyFinalGSE35978 <- get("studyFinal", envir = StanleyGSE35978)

studyFinalGSE35978$Cortex$Metadata$CharVec <- apply(studyFinalGSE35978$Cortex$Metadata %>% select(Sex, pH, PMI, Profile), 1,
                                                    function(x) paste0(x, collapse="_"))
studyFinalGSE35978$Cerebellum$Metadata$CharVec <- apply(studyFinalGSE35978$Cerebellum$Metadata %>% select(Sex, pH, PMI, Profile), 1,
                                                        function(x) paste0(x, collapse="_"))

Matched <- data.frame(GSE35978 = studyFinalGSE35978$Cortex$Metadata$CharVec,
                      Stanley = NA)
Matched$Stanley <- MetadataAllStanley$CharVec2[match(studyFinalGSE35978$Cortex$Metadata$CharVec, MetadataAllStanley$CharVec2)]


studyFinalGSE35978$Cortex$Metadata <- cbind(studyFinalGSE35978$Cortex$Metadata,
                                            MetadataAllStanley[match(Matched$Stanley,
                                                                     MetadataAllStanley$CharVec2),
                                                               c("CharVec2","StanleyID", "CollectionType")])
Matched <- data.frame(GSE35978 = studyFinalGSE35978$Cerebellum$Metadata$CharVec,
                      Stanley = NA)
Matched$Stanley <- MetadataAllStanley$CharVec2[match(studyFinalGSE35978$Cerebellum$Metadata$CharVec, MetadataAllStanley$CharVec2)]

studyFinalGSE35978$Cerebellum$Metadata <- cbind(studyFinalGSE35978$Cerebellum$Metadata,
                                            MetadataAllStanley[match(Matched$Stanley,
                                                                     MetadataAllStanley$CharVec2),
                                                               c("CharVec2","StanleyID", "CollectionType")])

GSE53987 <- new.env()
load("GeneralResults/studyFinalGSE53987.rda", envir = GSE53987)
studyFinalGSE53987 <- get("studyFinal", envir = GSE53987)
studyFinalGSE53987 <- lapply(studyFinalGSE53987, function(region){
  region$Metadata$CollectionType <- "Pittsburgh" %>% factor
  region
})

StanleyArray <- new.env()
load("GeneralResults/studyFinalStanleyArray.rda", envir = StanleyArray)
studyFinalSA <- get("studyFinal", envir = StanleyArray)

StanleyConsort <- new.env()
load("GeneralResults/studyFinalStanleyConsortium.rda", envir = StanleyConsort)
studyFinalSC <- get("studyFinal", envir = StanleyConsort)

MetadataGSE53987_hpc <- studyFinalGSE53987$Hippocampus$Metadata %>% mutate(Study2 = "Hippocampus_GSE53987")
MetadataGSE53987_str <- studyFinalGSE53987$Striatum$Metadata %>% mutate(Study2 = "Striatum_GSE53987")
MetadataGSE53987_pfc <- studyFinalGSE53987$Cortex$Metadata 
MetadataGSE35978_pfc <- studyFinalGSE35978$Cortex$Metadata %>% .[!duplicated(.$StanleyID),]
MetadataGSE35978_crb <- studyFinalGSE35978$Cerebellum$Metadata %>% mutate(Study2 = "Cerebellum_GSE35978") %>% .[!duplicated(.$StanleyID),]
MetadataSC_thlm <- studyFinalSC$Study16Kemether$Metadata %>% mutate(Study2 = "Thalamus_Study16Kemether")
MetadataSA_hpc <- studyFinalSA$Study17Laeng$Metadata %>% mutate(Study2 = "Hippocampus_Study17Laeng")

MetadataGSE53987_pfc_hpc <- MetadataGSE53987_pfc %>% filter(CommonName %in% MetadataGSE53987_hpc$CommonName) %>% mutate(Study2 = "Hippocampus_GSE53987")
MetadataGSE53987_pfc_str <- MetadataGSE53987_pfc %>% filter(CommonName %in% MetadataGSE53987_str$CommonName) %>% mutate(Study2 = "Striatum_GSE53987")
MetadataGSE35978_pfc_crb <- MetadataGSE35978_pfc %>% filter(StanleyID %in% MetadataGSE35978_crb$StanleyID) %>% mutate(Study2 = "Cerebellum_GSE35978")
MetadataSC_pfc_thlm <- MetadataGSE35978_pfc %>% filter(StanleyID %in% MetadataSC_thlm$StanleyID) %>% mutate(Study2 = "Thalamus_Study16Kemether")
MetadataSA_pfc_hpc <- MetadataGSE35978_pfc %>% filter(StanleyID %in% MetadataSA_hpc$StanleyID) %>% mutate(Study2 = "Hippocampus_Study17Laeng")

rm(list = ls(pat="_pfc$"))

for(study in ls(pat = "Metadata.*_")[!grepl("pfc", ls(pat = "Metadata.*_"))]){
  meta <- eval(as.name(study))
  metaCombinedMelt <- melt(meta, id.vars = grep("_Genes", names(meta), invert = TRUE, value = TRUE),
                           variable.name = "CellType",
                           value.name = "MGP")
  metaCombinedMelt$CellType <- sapply(metaCombinedMelt$CellType, function(x) gsub("_Genes", "", x))
  metaCombinedMelt$StudyRegion <- paste(metaCombinedMelt$Study %>% sapply(function(x) strsplit(as.character(x), "_")[[1]][1]),
                                        metaCombinedMelt$NeuExpRegion, sep = "_")
  EstimateDeltaCell <- sapply(unique(metaCombinedMelt$StudyRegion), function(study){
    studyData <- sapply(metaCombinedMelt %>% filter(StudyRegion == study) %>% .$CellType %>% unique, function(celltype){
      data = metaCombinedMelt %>% filter(StudyRegion == study, CellType == celltype) %>% droplevels()
      if(sum(is.na(data$MGP)) == nrow(data)){
        BP = NA
        BPpval = NA
        SCZ = NA
        SCZpval = NA
      } else {
        groups = levels(data$Profile)
        if("BP" %in% groups){
          data2  = data %>% filter(Profile %in% c("Cont", "BP")) %>% droplevels()
          data2$Profile <- relevel(data2$Profile, ref = "BP")
          wcx <- wilcox.test(MGP~Profile, data = data2, conf.int = TRUE)
          BP = wcx$estimate
          BPpval = wcx$p.value
        } else {
          BP = NA
          BPpval = NA
        }
        if ("SCZ" %in% groups){
          data2  = data %>% filter(Profile %in% c("Cont", "SCZ")) %>% droplevels()
          data2$Profile <- relevel(data2$Profile, ref = "SCZ")
          wcx = wilcox.test(MGP~Profile, data = data2, conf.int = TRUE)
          SCZ = wcx$estimate
          SCZpval = wcx$p.value
        } else {
          SCZ = NA
          SCZpval = NA
        }
      }
      data.frame(Study = data$Study2 %>% unique,
                 Region = data$NeuExpRegion %>% unique,
                 CellType = celltype,
                 BP = BP,
                 SCZ = SCZ,
                 BPpval = BPpval,
                 SCZpval = SCZpval)
    }, simplify = FALSE) %>% do.call(rbind, .)
  }, simplify = FALSE) %>% do.call(rbind, .)
  EstimateDeltaCell[grep("BP$|SCZ$", names(EstimateDeltaCell))] <- apply(EstimateDeltaCell %>% select(matches("BP$|SCZ$")), c(1,2), function(x){
    if(is.na(x)){
      NA
    } else{
      round(x, digits=2)
    }
  }) %>% data.frame
  
  EstimateDeltaCell[grep("pval", names(EstimateDeltaCell))] <- apply(EstimateDeltaCell %>% select(matches("pval")), c(1,2), function(x){
    if(is.na(x)){
      NA
    } else if (x > 0.05){
      signif(x, digits=2)
    } else {
      sigStar = GetStar(x)
      paste0(scientific(x, digits=2), sigStar)
    }
  }) %>% data.frame
  
  EstimateDeltaCell %<>% select(Region, CellType, BP, BPpval, SCZ, SCZpval)
  names(EstimateDeltaCell)[names(EstimateDeltaCell) == "BP"] <- "BPshift"
  names(EstimateDeltaCell)[names(EstimateDeltaCell) == "SCZ"] <- "SCZshift"
  write.table(EstimateDeltaCell, paste0("GeneralResults/MGPshifts", study, ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)
}

for(meta in ls(pat = "Metadata.*_")){
  metaData <- eval(as.name(meta))
  metaData$SubjectID <- if("StanleyID" %in% names(metaData)){
    metaData$StanleyID
    } else {
      metaData$Filename
      }
  assign(meta, metaData)
}

metaCombined <- sapply(ls(pat = "Metadata.*_"), function(studyMeta){
  meta <- eval(as.name(studyMeta)) %>% select(matches("SubjectID|CommonName|BioGender|age$|_Genes|ph$|CollectionType|NeuExp|Profile|Study2", ignore.case = TRUE))
  names(meta)[grepl("BioGender", names(meta), ignore.case = TRUE)] <- "Sex"
  names(meta)[grepl("age", names(meta), ignore.case = TRUE)] <- "Age"
  names(meta)[grepl("ph$", names(meta), ignore.case = TRUE)] <- "pH"
  meta$Study <- gsub("Metadata", "", studyMeta) %>% factor
  meta$NeuExpRegion <- factor(meta$NeuExpRegion)
  meta %<>% select(Profile, Age, pH, Sex, CommonName, SubjectID, Study, Study2, CollectionType, NeuExpRegion, matches("Astrocyte|Bergmann"))
  if(grepl("pfc", studyMeta)){
    meta %<>% mutate(AstrocyteCortex_Genes = Astrocyte_Genes)
  } else if (grepl("crb", studyMeta)){
    names(meta)[grepl("Bergm", names(meta))] <- "Astrocyte_Genes"
  }
  meta
}, simplify = FALSE) %>% do.call(rbind, .) %>% droplevels()
metaCombined$dataset <- sapply(metaCombined$Study, function(x) strsplit(as.character(x), "_")[[1]][1]) %>% factor


metaCombinedMelt <- melt(metaCombined, id.vars = grep("_Genes", names(metaCombined), invert = TRUE, value = TRUE),
                         variable.name = "CellType",
                         value.name = "MGP")
metaCombinedMelt$CellType <- sapply(metaCombinedMelt$CellType, function(x) gsub("_Genes", "", x))
metaCombinedMelt$StudyRegion <- paste(metaCombinedMelt$Study %>% sapply(function(x) strsplit(as.character(x), "_")[[1]][1]),
                                        metaCombinedMelt$NeuExpRegion, sep = "_")

EstimateDeltaCell <- sapply(unique(metaCombinedMelt$StudyRegion), function(study){
  studyData <- sapply(metaCombinedMelt %>% filter(StudyRegion == study) %>% .$CellType %>% unique, function(celltype){
    data = metaCombinedMelt %>% filter(StudyRegion == study, CellType == celltype) %>% droplevels()
    if(sum(is.na(data$MGP)) == nrow(data)){
      BP = NA
      BPpval = NA
      SCZ = NA
      SCZpval = NA
    } else {
      groups = levels(data$Profile)
      if("BP" %in% groups){
        data2  = data %>% filter(Profile %in% c("Cont", "BP")) %>% droplevels()
        data2$Profile <- relevel(data2$Profile, ref = "BP")
        wcx <- wilcox.test(MGP~Profile, data = data2, conf.int = TRUE)
        BP = wcx$estimate
        BPpval = wcx$p.value
      } else {
        BP = NA
        BPpval = NA
      }
      if ("SCZ" %in% groups){
        data2  = data %>% filter(Profile %in% c("Cont", "SCZ")) %>% droplevels()
        data2$Profile <- relevel(data2$Profile, ref = "SCZ")
        wcx = wilcox.test(MGP~Profile, data = data2, conf.int = TRUE)
        SCZ = wcx$estimate
        SCZpval = wcx$p.value
      } else {
        SCZ = NA
        SCZpval = NA
      }
    }
    data.frame(Study = data$Study2 %>% unique,
               Region = data$NeuExpRegion %>% unique,
               CellType = celltype,
               BP = BP,
               SCZ = SCZ,
               BPpval = BPpval,
               SCZpval = SCZpval)
  }, simplify = FALSE) %>% do.call(rbind, .)
}, simplify = FALSE) %>% do.call(rbind, .)

EstimateDeltaCell$BPstar <- GetStar(EstimateDeltaCell$BPpval)

EstimateDeltaCell$SCZstar <- GetStar(EstimateDeltaCell$SCZpval)
rownames(EstimateDeltaCell) <- c(1:nrow(EstimateDeltaCell))
EstimateDeltaCell$Region %<>% relevel(EstimateDeltaCell$Study, ref = "Cortex")
EstimateDeltaCell$RegionCellType <- paste(EstimateDeltaCell$Region, EstimateDeltaCell$CellType, sep="_")
EstimateDeltaCell %<>% filter(!(Region == "Cortex" & CellType == "AstrocyteCortex"))
                          
EstimateDeltaCell %<>% mutate(CellType2 = as.character(CellType))
EstimateDeltaCell$CellType2 <- apply(EstimateDeltaCell %>% select(CellType2, Region), 1, function(x){
  if(x[1] == "Astrocyte" & x[2] == "Cortex"){
    "Astrocytes \n(Cortex)" 
  } else if (x[1] == "Astrocyte" & x[2] != "Cortex"){
    "Astrocytes \n(region)"
  } else if (x[1] == "AstrocyteCortex" ){
    "Cortical astrocyte MGP \n(region)"
  }
})
EstimateDeltaCell$CellType2 <- factor(EstimateDeltaCell$CellType2, levels = c("Astrocytes \n(Cortex)",
                                                                          "Astrocytes \n(region)",
                                                                          "Cortical astrocyte MGP \n(region)"))
EstimateDeltaCell %<>% droplevels()
EstimateDeltaCell$Study <- factor(as.character(EstimateDeltaCell$Study), levels = c("Striatum_GSE53987",
                                                                        "Hippocampus_GSE53987",
                                                                        "Hippocampus_Study17Laeng",
                                                                        "Cerebellum_GSE35978",
                                                                        "Thalamus_Study16Kemether"))
EstimateDeltaCellMelt <- melt(data.table(EstimateDeltaCell), id.vars = c("Study", "CellType2"), measure.vars = patterns("BP$|SCZ$", "pval", "star", cols = names(EstimateDeltaCell)),
                              variable.name = "Profile",
                              value.name = c("Estimate", "pVal", "Star")) %>% data.frame

EstimateDeltaCellMelt$Profile <- sapply(EstimateDeltaCellMelt$Profile, function(x) {
  switch(as.character(x), "1"="BP", "2"="SCZ")
}) %>% factor

EstimateDeltaCellMelt$Region <- sapply(EstimateDeltaCellMelt$CellType2, function(x){
  if(grepl("Cortex", x)){
    "Cortex"
  } else {
    "Extra-cortical regions"
  }
}) %>% factor(levels = c("Extra-cortical regions", "Cortex"))

EstimatePlot <- ggplot(EstimateDeltaCellMelt, aes(CellType2, Study)) +
  theme_light(base_size = 16) +
  theme(axis.text.y = element_text(size = rel(1.4)),
        axis.text.x = element_text(size = rel(1.4)),
        strip.background = element_rect(fill="grey40")) +
  labs(x="", y="") +
  facet_grid(Profile~Region, scales = "free_x", space = "free_x") +
  geom_tile(aes(fill = Estimate), colour = "white") +
  theme(axis.text.x = element_text(angle = 40,hjust = 1)) +
  scale_fill_gradient2(low = "darkolivegreen", high = "darkorange3", mid = "grey90",na.value = "white",
                       midpoint = 0, limit = c(-max(abs(c(EstimateDeltaCellMelt$Estimate)), na.rm=TRUE),
                                               max(abs(c(EstimateDeltaCellMelt$Estimate)), na.rm=TRUE)), name = "Effect")+
  geom_text(aes(label=Star), color="black", size=7)
ggsave("GeneralResults/EstimatePlot_extracortical2.pdf", plot = EstimatePlot, width = 12, height = 6, units = "in",
       dpi=300)

#SSTReln Hippocampus

metaCombinedSSTReln <- sapply(ls(pat = "Metadata.*_hpc")[!grepl("pfc", ls(pat = "Metadata.*_hpc") )], function(studyMeta){
  meta <- eval(as.name(studyMeta)) %>% select(matches("SubjectID|CommonName|BioGender|age$|GabaSSTReln_Genes|ph$|CollectionType|NeuExp|Profile|Study2", ignore.case = TRUE))
  names(meta)[grepl("BioGender", names(meta), ignore.case = TRUE)] <- "Sex"
  names(meta)[grepl("age", names(meta), ignore.case = TRUE)] <- "Age"
  names(meta)[grepl("ph$", names(meta), ignore.case = TRUE)] <- "pH"
  meta$Study <- gsub("Metadata", "", studyMeta) %>% factor
  meta$NeuExpRegion <- factor(meta$NeuExpRegion)
  meta %<>% select(Profile, Age, pH, Sex, CommonName, SubjectID, Study, Study2, CollectionType, NeuExpRegion, GabaSSTReln_Genes)
  meta
}, simplify = FALSE) %>% do.call(rbind, .) %>% droplevels()
metaCombinedSSTReln$dataset <- sapply(metaCombinedSSTReln$Study, function(x) strsplit(as.character(x), "_")[[1]][1]) %>% factor
levels(metaCombinedSSTReln$Study) <- c("GSE53987", "Study17Laeng")

names(metaCombinedSSTReln)[names(metaCombinedSSTReln) == "GabaSSTReln_Genes"] <- "GabaSSTReln"

EstimateDeltaCellhpc <- sapply(unique(metaCombinedSSTReln$Study) %>% as.character(), function(study){
  data <- metaCombinedSSTReln %>% filter(Study == study) %>% droplevels()
  for(grp in c("BP", "SCZ")){
    data2  = data %>% filter(Profile %in% c("Cont", grp)) %>% droplevels()
    data2$Profile <- relevel(data2$Profile, ref = grp)
    wcx <- wilcox.test(GabaSSTReln~Profile, data = data2, conf.int = TRUE)
    assign(grp, wcx$estimate)
    assign(paste0(grp, "pval"), wcx$p.value)
    }
  data.frame(Study = study,
             Study2 = data$Study2 %>% unique,
             CellType = "GabaSSTReln",
             BP = BP,
             SCZ = SCZ,
             BPpval = BPpval,
             SCZpval = SCZpval)
}, simplify = FALSE) %>% do.call(rbind, .)
EstimateDeltaCellhpc$BPstar <- GetStar(EstimateDeltaCellhpc$BPpval)
EstimateDeltaCellhpc$SCZstar <- GetStar(EstimateDeltaCellhpc$SCZpval)
EstimateDeltaCellhpc$Study <- factor(EstimateDeltaCellhpc$Study, levels = EstimateDeltaCellhpc$Study)
EstimateDeltaCellhpcMelt <- melt(data.table(EstimateDeltaCellhpc), id.vars = c("Study", "Study2", "CellType"), measure.vars = patterns("BP$|SCZ$", "pval", "star", cols = names(EstimateDeltaCell)),
                              variable.name = "Profile",
                              value.name = c("Estimate", "pVal", "Star")) %>% data.frame
EstimateDeltaCellhpcMelt$x <- EstimateDeltaCellhpcMelt$Profile %>% as.numeric() %>% +1
EstimateDeltaCellhpcMelt$y <- 1.01
EstimateDeltaCellhpcMelt %<>% droplevels

SampleNum <- metaCombinedSSTReln  %>% filter(!is.na(.$GabaSSTReln))  %>% 
  group_by(Study, Profile) %>%
  summarise(n=n()) %>% data.frame()
ContMedian <- metaCombinedSSTReln %>% filter(Profile == "Cont") %>% group_by(Study) %>%
  summarise(median = median(GabaSSTReln, na.rm = T)) %>% data.frame

metaCombinedSSTReln$Profile2 <- apply(metaCombinedSSTReln %>% select(Study, Profile), 1, function(x){
  SampleNum %>% filter(Study == x[1], Profile == x[2]) %>% .$n
})

metaCombinedSSTReln %<>% mutate(Profile2 = paste0(Profile, "\n(n=", Profile2, ")"))
metaCombinedSSTReln$Profile2 <- factor(metaCombinedSSTReln$Profile2,
                                      levels = c("Cont\n(n=17)", "BP\n(n=16)", "SCZ\n(n=13)",
                                                 "Cont\n(n=19)", "BP\n(n=20)", "SCZ\n(n=20)" ))
SSTRelnPlot <- ggplot(metaCombinedSSTReln, aes(Profile2, GabaSSTReln)) +
  theme_bw(base_size = 16) +
  labs(x = "", y = "GabaSSTReln MGP") +
  ylim(0,1.1) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(size = rel(0.9)),
        axis.text.y = element_text(size = rel(0.9))) +
  #geom_violin(width = 0.8)+
  geom_boxplot(width = 0.5, outlier.size = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~Study, scales = "free_x") +
  geom_hline(data = ContMedian, aes(yintercept = median), color = "red", linetype="dashed") +
  geom_text(data = EstimateDeltaCellhpcMelt,
            aes(x, y, label = Star), size=6, inherit.aes = FALSE)
ggsave("GeneralResults/SSTRelnHPC.pdf", plot = SSTRelnPlot, width = 8, height = 4, units = "in",
       dpi=300)
