GenScriptPath = "/home/ltoker/Rscripts/"
source(paste0(GenScriptPath,"general_functions.R"))
source("projectFunc.R")
packageF("sva")
packageF("gplots")
packageF("scales")
source("ProjectScripts/Cell_type_PCA.R")

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

studyFinalGSE35978$Cortex$Metadata$CharVec <- apply(studyFinalGSE35978$Cortex$Metadata %>% select(Sex, ph, PMI, Profile), 1,
                                                    function(x) paste0(x, collapse="_"))
studyFinalGSE35978$Cerebellum$Metadata$CharVec <- apply(studyFinalGSE35978$Cerebellum$Metadata %>% select(Sex, ph, PMI, Profile), 1,
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
rm(MetadataAllStanley)
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

MetadataGSE53987_pfc <- studyFinalGSE53987$Cortex$Metadata
MetadataGSE53987_hpc <- studyFinalGSE53987$Hippocampus$Metadata
MetadataGSE53987_str <- studyFinalGSE53987$Striatum$Metadata
MetadataGSE35978_pfc <- studyFinalGSE35978$Cortex$Metadata
MetadataGSE35978_crb <- studyFinalGSE35978$Cerebellum$Metadata
MetadataSC_thlm <- studyFinalSC$Study16Kemether$Metadata
MetadataSA_hpc <- studyFinalSA$Study17Laeng$Metadata

for(meta in ls(pat = "Metadata")){
  metaData <- eval(as.name(meta))
  metaData$SubjectID <- if("StanleyID" %in% names(metaData)){
    metaData$StanleyID
    } else {
      metaData$Filename
      }
  assign(meta, metaData)
}

metaCombined <- sapply(ls(pat = "Metadata"), function(studyMeta){
  meta <- eval(as.name(studyMeta)) %>% select(matches("SubjectID|BioGender|age$|_Genes|ph$|CollectionType|NeuExp|Profile", ignore.case = TRUE))
  names(meta)[grepl("BioGender", names(meta), ignore.case = TRUE)] <- "Sex"
  names(meta)[grepl("age", names(meta), ignore.case = TRUE)] <- "Age"
  names(meta)[grepl("ph$", names(meta), ignore.case = TRUE)] <- "pH"
  meta$Study <- gsub("Metadata", "", studyMeta) %>% factor
  meta$NeuExpRegion <- factor(meta$NeuExpRegion)
  meta %<>% select(Profile, Age, pH, Sex, SubjectID, Study, CollectionType, NeuExpRegion, matches("Astrocyte|Bergmann"))
  if(grepl("pfc", studyMeta)){
    meta %<>% mutate(AstrocyteCortex_Genes = Astrocyte_Genes)
  } else if (grepl("crb", studyMeta)){
    names(meta)[grepl("Bergm", names(meta))] <- "Astrocyte_Genes"
  }
  meta
}, simplify = FALSE) %>% do.call(rbind, .) %>% droplevels()

metaCombinedMelt <- melt(metaCombined, id.vars = grep("_Genes", names(metaCombined), invert = TRUE, value = TRUE),
                         variable.name = "CellType",
                         value.name = "MGP")
metaCombinedMelt$CellType <- sapply(metaCombinedMelt$CellType, function(x) gsub("_Genes", "", x))
metaCombinedMelt$CohortRegion <- paste(metaCombinedMelt$CollectionType,
                                        metaCombinedMelt$Study,
                                        metaCombinedMelt$NeuExpRegion, sep = "_")

EstimateDeltaCell <- sapply(unique(metaCombinedMelt$CohortRegion), function(study){
  studyData <- sapply(metaCombinedMelt %>% filter(CohortRegion == study) %>% .$CellType %>% unique, function(celltype){
    data = metaCombinedMelt %>% filter(CohortRegion == study, CellType == celltype) %>% droplevels()
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
    data.frame(Study = study,
               Collection = data$Collection %>% unique,
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
EstimateDeltaCell$Study2 <- paste(EstimateDeltaCell$Collection,
                                  sapply(EstimateDeltaCell$Study,
                                         function(x) strsplit(as.character(x), "_")[[1]][2]),
                                  EstimateDeltaCell$Region, sep="_")
EstimateDeltaCell %<>% filter(!(Region == "Cortex" & CellType == "AstrocyteCortex"))
EstimateDeltaCell <- rbind(EstimateDeltaCell,
                           EstimateDeltaCell %>%
                             filter(RegionCellType == "Cortex_Astrocyte", Study2 == "Array_GSE35978_Cortex") %<>% mutate(Study2 = "Array_SA_Hippocampus"),
                           EstimateDeltaCell %>%
                             filter(RegionCellType == "Cortex_Astrocyte", Study2 == "Consortium_GSE35978_Cortex") %<>% mutate(Study2 = "Consortium_SC_Thalamus"),
                           EstimateDeltaCell %>%
                             filter(RegionCellType == "Cortex_Astrocyte", Study2 == "Array_GSE35978_Cortex") %<>% mutate(Study2 = "Array_GSE35978_Cerebellum"),
                           EstimateDeltaCell %>%
                             filter(RegionCellType == "Cortex_Astrocyte", Study2 == "Consortium_GSE35978_Cortex") %<>% mutate(Study2 = "Consortium_GSE35978_Cerebellum"),
                           EstimateDeltaCell %>%
                             filter(RegionCellType == "Cortex_Astrocyte", Study2 == "Pittsburgh_GSE53987_Cortex") %<>% mutate(Study2 = "Pittsburgh_GSE53987_Hippocampus"),
                           EstimateDeltaCell %>%
                             filter(RegionCellType == "Cortex_Astrocyte", Study2 == "Pittsburgh_GSE53987_Cortex") %<>% mutate(Study2 = "Pittsburgh_GSE53987_Striatum"))
EstimateDeltaCell <- EstimateDeltaCell[!grepl("Cortex", EstimateDeltaCell$Study2),]

EstimateDeltaCell %<>% mutate(CellType2 = as.character(CellType))
EstimateDeltaCell$CellType2[EstimateDeltaCell$Region == "Cortex"] <- paste("Cortical", EstimateDeltaCell$CellType2[EstimateDeltaCell$Region == "Cortex"], sep = "_")
EstimateDeltaCell$CellType2 <- factor(EstimateDeltaCell$CellType2, levels = c("Cortical_Astrocyte",
                                                                          "Astrocyte",
                                                                          "AstrocyteCortex"))
EstimateDeltaCell$Study2 <- factor(EstimateDeltaCell$Study2, levels = c("Pittsburgh_GSE53987_Striatum",
                                                                        "Pittsburgh_GSE53987_Hippocampus",
                                                                        "Array_SA_Hippocampus",
                                                                        "Array_GSE35978_Cerebellum",
                                                                        "Consortium_GSE35978_Cerebellum",
                                                                        "Consortium_SC_Thalamus"))
EstimateDeltaCellMelt <- melt(data.table(EstimateDeltaCell), id.vars = c("Study2", "CellType2"), measure.vars = patterns("BP$|SCZ$", "pval", "star", cols = names(EstimateDeltaCell)),
                              variable.name = "Profile",
                              value.name = c("Estimate", "pVal", "Star")) %>% data.frame

EstimateDeltaCellMelt$Profile <- sapply(EstimateDeltaCellMelt$Profile, function(x) {
  switch(as.character(x), "1"="BP", "2"="SCZ")
}) %>% factor

EstimatePlot <- ggplot(EstimateDeltaCellMelt, aes(CellType2, Study2)) +
  theme_light(base_size = 16) +
  theme(axis.text.y = element_text(size = rel(1)),
        axis.text.x = element_text(size = rel(1.3)),
        strip.background = element_rect(fill="grey40")) +
  facet_wrap(~Profile) +
  geom_tile(aes(fill = Estimate), colour = "white") +
  theme(axis.text.x = element_text(angle = 40,hjust = 1)) +
  scale_fill_gradient2(low = "darkolivegreen", high = "darkorange3", mid = "grey90",na.value = "white",
                       midpoint = 0, limit = c(-0.4, 0.4), name = "Effect")+
  geom_text(aes(label=Star), color="black", size=7)
ggsave("GeneralResults/EstimatePlot_extracortical.pdf", plot = EstimatePlot, width = 12, height = 6, units = "in",
       dpi=300)
