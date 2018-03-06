source("SetUp.R")
source(paste0(GenScriptPath, "Cell_type_PCA.R"))


PlotStudies <- function(data, cell, txtSize = 18, pvalTxtSize = 5, ylab = "Relative MGP"){
  cellGenes = paste0(cell, "_Genes")
  data <- melt(data, id.var = c("CommonName", "Profile", "Profile2", "Study"),
               measure.vars = cellGenes, variable.name = "CellType", value.name = "MGP")
  pvalAll <- sapply(cellGenes, function(cell){
    sapply(levels(data$Profile)[-1], function(grp){
      data = data %>% filter(Profile %in% c("Cont", grp), CellType == cell)
      WlCx <- wilcox.test(formula(paste0("MGP", "~Profile")), data=data)
      if(WlCx$p.value > 0.05){
        pval=paste0("p = ", round(WlCx$p.value, digits=2))
      } else{
        pval = paste0("p = ", scientific(WlCx$p.value))
      }
      char=GetSigChar(WlCx$p.value)
      paste0(char, pval)
    }, simplify = TRUE)
  }) %>% data.frame() %>%
    mutate(Profile = rownames(.)) %>%
    melt(id.vars = "Profile", variable.name = "CellType", value.name= "pVal") %>%
    mutate(x = rep(c(2:3), 2), y = 1.1)
  
  studyMedians =  data %>% filter(!is.na(MGP)) %>% 
    group_by(Profile2, Study, CellType) %>% 
    summarise(medianMGP = median(MGP)) %>% data.frame
  
  GroupMedian <-  data %>% filter(!is.na(MGP)) %>% 
    group_by(CellType, Profile2) %>% 
    summarise(medianMGP = median(MGP)) %>% data.frame %>%
    mutate(x = rep(c(1:3)-0.25, length(cellGenes)), 
           xend=rep(c(1:3)+0.25, length(cellGenes)))

  plot <- ggplot(data, aes(Profile2,
                           MGP),color = "study") +
    theme_classic(base_size = txtSize) +
    theme(axis.text.x = element_text(),
          axis.text.y = element_text(size = rel(0.8)),
          panel.grid = element_blank(),
          legend.title = element_text(size=rel(0.8)),
          axis.line.x = element_line(color="black", size=0.5),
          axis.line.y = element_line(color="black", size=0.5),
          panel.background = element_blank()) +
    labs(x = "", y= ylab)+
    geom_violin( width=0.8, color="grey", size=0.1, fill= "gray", alpha=0.8) +
    geom_boxplot(outlier.shape = NA, width = 0.4, colour="black", alpha=0.8) +
    geom_segment(data = GroupMedian, mapping=aes(x = x, y = medianMGP,
                 xend=xend, yend = medianMGP),
                 color="black", size=2) +
    geom_jitter(data=studyMedians, aes(y = medianMGP, color = Study),size=3, width = 0.3, alpha = 0.8) +
    scale_colour_manual(values = c("gray35", "chartreuse4", "dodgerblue4", "brown4",
                                   "darksalmon", "hotpink", "red", "darkorchid4", "darkorange1"),
                        name = "Dataset") +
    geom_text(data = pvalAll,  aes(x, y, label = pVal), size=pvalTxtSize, inherit.aes = FALSE) +
    facet_wrap(~CellType)
  print(plot)
  return(plot)
}

PlotSingleStudy <- function(data, celltype, txtSize = 14, pvalSize = 10, lineCol = "red"){
  ContMedian <- data %>% filter(CellType == celltype, Profile == "Cont") %>% group_by(Study2) %>%
    summarise(median = median(MGP, na.rm = T)) %>% data.frame
  ggplot(data %>% filter(CellType == celltype), aes(Profile2, MGP)) +
    theme_classic(base_size = txtSize) +
    labs(x = "", y = paste0(celltype, " MGP")) + 
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(), 
          axis.text.x = element_text(size = rel(0.9)),
          axis.text.y = element_text(size = rel(0.9))) +
    ylim(0, 1.1) +
    #geom_violin(width = 0.8)+
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size=1)+
    facet_wrap(~Study2, scales = "free_x") +
    geom_text(data = EstimateDeltaCellMelt %>% filter(CellType == celltype),
              aes(x, y, label = Star), size=pvalSize, inherit.aes = FALSE) +
    geom_hline(data = ContMedian, aes(yintercept = median), color = lineCol, linetype="dashed")
}

GetLM <- function(data, celltype, group){
  data %<>% filter(Profile %in% c("Cont", group))
  data2 <- data[!is.na(data$pH),]
  lm0 <- lmer(as.formula(paste0(celltype, "~1 + (1|Study) + (1|SubjectID)")), data = data2, REML = FALSE)
  lm1a <- lmer(as.formula(paste0(celltype, "~Sex + (1|Study) + (1|SubjectID)")), data = data2 , REML = FALSE)
  lm1b <- lmer(as.formula(paste0(celltype, "~Age + (1|Study) + (1|SubjectID)")), data = data2 , REML = FALSE)
  lm1c <- lmer(as.formula(paste0(celltype, "~pH + (1|Study) + (1|SubjectID)")), data = data2 , REML = FALSE)
  lm1d <- lmer(as.formula(paste0(celltype, "~PMI + (1|Study) + (1|SubjectID)")), data = data2 , REML = FALSE)
  
  covarEff = "(1|Study) + (1|SubjectID)"
  if (anova(lm0, lm1a)$`Pr(>Chisq)`[2] < 0.05 ) {
    covarEff <- paste0("Sex + ", covarEff)
  }
  if (anova(lm0, lm1b)$`Pr(>Chisq)`[2] < 0.05) {
    covarEff <- paste0("Age +", covarEff)
  }
  if (anova(lm0, lm1c)$`Pr(>Chisq)`[2] < 0.05) {
    covarEff <- paste0("pH +", covarEff)
  }
  if (anova(lm0, lm1d)$`Pr(>Chisq)`[2] < 0.05) {
    covarEff <- paste0("PMI +", covarEff)
  }
  
  lm2 <- lmer(as.formula(paste0(celltype,"~", covarEff)), data = data, REML = FALSE)
  lm2ast <- lmer(as.formula(paste0(celltype,"~ Astrocyte_Genes + ", covarEff)), data = data, REML = FALSE)
  lm2PV <- lmer(as.formula(paste0(celltype,"~ GabaPV_Genes + ", covarEff)), data = data, REML = FALSE)
  lm2b <- lmer(as.formula(paste0(celltype, "~Sex + Age + pH + PMI + (1|Study) + (1|SubjectID)")), data = data , REML = FALSE)
  lm2bast <- lmer(as.formula(paste0(celltype, "~Sex + Age + pH + PMI + Astrocyte_Genes + (1|Study) + (1|SubjectID)")), data = data , REML = FALSE)
  lm2bPV <- lmer(as.formula(paste0(celltype, "~Sex + Age + pH + PMI + GabaPV_Genes + (1|Study) + (1|SubjectID)")), data = data , REML = FALSE)
  lm3 <- lmer(as.formula(paste0(celltype, "~Profile +", covarEff)), data = data, REML = FALSE)
  lm3ast <- lmer(as.formula(paste0(celltype, "~Profile + Astrocyte_Genes + ", covarEff)), data = data, REML = FALSE)
  lm3PV <- lmer(as.formula(paste0(celltype, "~Profile + GabaPV_Genes + ", covarEff)), data = data, REML = FALSE)
  lm3b <- lmer(as.formula(paste0(celltype, "~Profile + Sex + Age + pH + PMI +  (1|Study) + (1|SubjectID)")), data = data, REML = FALSE)
  lm3bast <- lmer(as.formula(paste0(celltype, "~Profile + Sex + Age + pH + PMI + Astrocyte_Genes + (1|Study) + (1|SubjectID)")), data = data, REML = FALSE)
  lm3bPV <- lmer(as.formula(paste0(celltype, "~Profile + Sex + Age + pH + PMI +GabaPV_Genes + (1|Study) + (1|SubjectID)")), data = data, REML = FALSE)
  return(list(lm0 = lm0,
              lm1a = lm1a,
              lm1b = lm1b,
              lm1c = lm1c,
              lm2 = lm2,
              lm2ast = lm2ast,
              lm2PV = lm2PV,
              lm2b = lm2b,
              lm2bast = lm2bast,
              lm2bPV = lm2bPV,
              lm3 = lm3,
              lm3ast = lm3ast,
              lm3PV = lm3PV,
              lm3b = lm3b,
              lm3bast=lm3bast,
              lm3bPV=lm3bPV))
  
}

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

GetCI <- function(lmResults, model, CellTypes){
  CIresults <- sapply(CellTypes, function(celltype){
    CellType = paste0(as.character(celltype), "_Genes")
    if(as.character(celltype) %in% CellTypeRM2){
      lmResult = NA
      CI = NA
    } else {
      lmResult <- lmResults[[CellType]][[model]]
      CI <- lmResult %>% confint.merMod(method = "boot")
    }
    list(lmResult = lmResult, CI = CI)
  }, simplify = FALSE)
  names(CIresults) <- unique(CellTypes)
  return(CIresults)
}

GetResultSummary <- function(CIdata, CellTypes, Effects, CellTypeRM2){
  ResultSummary <- sapply(CellTypes, function(celltype){
    CellType = as.character(celltype)
    temp <- sapply(Effects, function(effect){
      if(CellType %in% CellTypeRM2){
        DF <- data.frame(row.names = CellType, CellType, NA,NA,NA)
      } else {
        lmResult <- CIdata[[CellType]]$lmResult %>% summary
        CI <- CIdata[[CellType]]$CI
        if(!effect %in% rownames(CI)){
          DF = data.frame(row.names = CellType, CellType, NA,NA,NA)
        } else {
          DF = data.frame(row.names = CellType, CellType, lmResult$coefficients[effect,"Estimate"], CI[effect,1],CI[effect,2])
        }
        
      }
      effect = gsub("Profile", "", effect)
      names(DF) <- c("CellType", effect, paste0(effect, "min"), paste0(effect, "max"))
      DF
    }, simplify = FALSE)
    temp %>% Reduce(function(x, y) merge(x, y, all=TRUE, by="CellType"), .)
  }, simplify = FALSE) %>% do.call(rbind, .)
  ResultSummary$CellType <- factor(ResultSummary$CellType, levels = c("Astrocyte" ,"Endothelial","Oligo", "OligoPrecursors",
                                                                "Microglia", "Microglia_activation", "Microglia_deactivation",
                                                                "GabaPV", "GabaRelnCalb",  "GabaVIPReln",
                                                                "Layer4Pyra", "Layer6bPyra", "PyramidalAll"))
  return(ResultSummary)
}

plotCI <- function(data, title = "Mixed model effects", rowNum = 1, scales = "fixed",
                   color_x = TRUE, setYlimAuto = FALSE, ylimMin = NA, ylimMax = NA, addXlab = TRUE){
  if(color_x){
    colorX = c("darkred", rep("black", 6), "darkred", rep("black", 6))
  } else {
    colorX = "black"
  }
  
  if(setYlimAuto){
    ylimMin = -1.05*max(c(abs(data$Min), abs(data$Max)))
    ylimMax = 1.05*max(c(abs(data$Min), abs(data$Max)))
  }
  
  if(addXlab){
    xlabSize = 1.3
  } else {
    xlabSize = 0
  }
  
  Plot <- ggplot(data, aes(x=CellType, y=Coefficient)) +
    theme_light(base_size = 16) +
    labs(title = title, y = "Coefficient", x="") +
    theme(axis.text.y = element_text(size = rel(1.3)),
          axis.text.x = element_text(size = rel(xlabSize), colour = colorX,
                                     angle = 90, hjust = 1, vjust=0.5),
          strip.background = element_rect(fill="grey40"),
          panel.grid.major.y =  element_blank(),
          panel.grid.minor.y =  element_blank(),
          panel.grid.major.x =  element_blank(),
          panel.grid.minor.x =  element_blank()) +
    scale_y_continuous(limits = c(ylimMin, ylimMax)) +
    facet_wrap(~Effect, nrow = rowNum, scales = scales ) +
    geom_point() +
    geom_errorbar(aes(ymin = Min, ymax = Max), na.rm = TRUE) +
    geom_abline(intercept = 0, slope = 0, color = "red", linetype="dashed")
  return(Plot)
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

MetadataAllStanley$CharVec2 <- apply(MetadataAllStanley %>% select(Sex, BrainPH, PMIh, Profile), 1, function(x) paste0(x, collapse="_"))


load("AnalysisStudies.rda")

for(study in AnalysisStudies){
  load(paste0(GeneralResultsPath, "/studyFinal", study, ".rda"))
  assign(paste0("studyFinal", study), studyFinal)
}

studyFinalGSE35978$Cortex$Metadata$CharVec <- apply(studyFinalGSE35978$Cortex$Metadata %>% select(Sex, pH, PMI, Profile), 1, function(x) paste0(x, collapse="_"))

Matched <- data.frame(GSE35978 = studyFinalGSE35978$Cortex$Metadata$CharVec,
                      Stanley = NA)
Matched$Stanley <- MetadataAllStanley$CharVec2[match(studyFinalGSE35978$Cortex$Metadata$CharVec, MetadataAllStanley$CharVec2)]

studyFinalGSE35978$Cortex$Metadata <- cbind(studyFinalGSE35978$Cortex$Metadata,MetadataAllStanley[match(Matched$Stanley,
                                                    MetadataAllStanley$CharVec2),
                                              c("CharVec2","StanleyID", "CollectionType")])
rm(MetadataAllStanley)

studyFinalMcLeanCortex$Cortex$Metadata$CollectionType <- "McLean66" %>% factor

studyFinalGSE53987 <- lapply(studyFinalGSE53987, function(region){
  region$Metadata$CollectionType <- "Pittsburgh" %>% factor
  region
})

studyFinalGSE21138$Cortex$Metadata$CollectionType <- "VBBN" %>% factor

studyFinalGSE80655$DLPFC$Metadata$CollectionType <- "PNDRC" %>% factor

for(study in names(studyFinalStanleyArray)[-1]){
 assign(paste0("Metadata", study), studyFinalStanleyArray[[study]]$Metadata) 
}

for(study in names(studyFinalStanleyConsortium)[-1]){
  assign(paste0("Metadata", study), studyFinalStanleyConsortium[[study]]$Metadata) 
}

MetadataGSE53987 <- studyFinalGSE53987$Cortex$Metadata
MetadataGSE35978 <- studyFinalGSE35978$Cortex$Metadata
MetadataMcLean <- studyFinalMcLeanCortex$Cortex$Metadata
MetadataGSE21138 <- studyFinalGSE21138$Cortex$Metadata
MetadataGSE80655 <- studyFinalGSE80655$DLPFC$Metadata

for(meta in ls(pat = "Metadata")){
  metaData <- eval(as.name(meta))
  metaData$SubjectID <- if("StanleyID" %in% names(metaData)){
    metaData$StanleyID
  } else {
    metaData$Filename
  }
  assign(meta, metaData)
}

#rm(list = ls(pat = "studyFinal"))

metaCombined <- sapply(ls(pat = "Metadata"), function(studyMeta){
  meta <- eval(as.name(studyMeta)) %>% select(matches("CommonName|SubjectID|Profile|BioGender|age$|_Genes|ph$|CollectionType|pmi$|rin$", ignore.case = TRUE))
  names(meta)[grepl("BioGender", names(meta), ignore.case = TRUE)] <- "Sex"
  names(meta)[grepl("age", names(meta), ignore.case = TRUE)] <- "Age"
  names(meta)[grepl("ph$", names(meta), ignore.case = TRUE)] <- "pH"
  if(sum(grepl("pmi$", names(meta), ignore.case = TRUE)) == 1){
    names(meta)[grepl("pmi$", names(meta), ignore.case = TRUE)] <- "PMI"
  } else if(sum(grepl("pmi$", names(meta), ignore.case = TRUE)) > 1){
    warning("more than one column matches the regex for pmi")
  } else {
    meta$PMI <- NA
  }
  
  if(sum(grepl("rin$", names(meta), ignore.case = TRUE)) == 1){
    names(meta)[grepl("rin$", names(meta), ignore.case = TRUE)] <- "RIN"
  } else if(sum(grepl("rin$", names(meta), ignore.case = TRUE)) > 1){
    warning("more than one column matches the regex for rin")
  } else {
    meta$RIN <- NA
  }
  meta$Study <- gsub("Metadata", "", studyMeta) %>% factor
  meta %<>% select(CommonName, Profile, Age, pH, Sex, PMI, RIN, SubjectID, Study, CollectionType, matches("_Genes"))
  meta
}, simplify = FALSE) %>% rbindlist %>% data.frame %>% droplevels()

#Remove cell types for which MGP was never calculated
CellTypeRM <- names(metaCombined)[apply(metaCombined, 2, function(x) sum(is.na(x))) == nrow(metaCombined)]
metaCombined <- metaCombined[,!names(metaCombined) %in% CellTypeRM]
metaCombined$Profile <- factor(as.character(metaCombined$Profile), levels = c("Cont", "BP", "SCZ"))
StudySum <- group_by(metaCombined %>% filter(!is.na(.$Astrocyte_Genes)), Study, Profile) %>% summarise(n=n()) %>% data.frame
names(metaCombined) <- sapply(names(metaCombined), function(x) gsub(" ", "", x))
metaCombined$Profile2 <- apply(metaCombined, 1, function(sbj){
  temp <- StudySum %>% filter(Study == sbj["Study"], Profile == sbj["Profile"])
  paste0(temp$Profile, "\n(n=", temp$n, ")")
})
metaCombined$Profile2 <- factor(metaCombined$Profile2, levels = unique(metaCombined$Profile2)[c(grep("Cont", unique(metaCombined$Profile2)),
                                                                                                grep("BP", unique(metaCombined$Profile2)),
                                                                                                grep("SCZ", unique(metaCombined$Profile2)))])
metaCombined$Study2 <-sapply(metaCombined$Study, function(study){
  if(study == "GSE21138"){
    "VBBN"
  } else if(study == "GSE35978"){
    "SMRI_A+C"
  } else if(study == "McLean"){
    "McLean66"
  } else if(study == "GSE53987"){
    "ACOME"
  } else if (study == "GSE80655"){
    "PNDRC"
  } else if(grepl("AltarC|Chen", study)){
    "SMRI_C"
  } else {
    "SMRI_A"
  }
})

metaCombined$Study2 <- paste0(metaCombined$Study, " (", metaCombined$Study2, ")")
metaCombined$SubjStudy <- apply(metaCombined %>% select(CommonName, Study), 1, function(x) paste0(x, collapse = "."))
names(metaCombined) <- sapply(names(metaCombined), function(x) gsub("\\.", "", x))

ExpAll <- list()
for(study in ls(pat = "studyFinal.*")){
  if(grepl("Stanl", study)){
    temp = eval(as.name(study))
    CrtxDS = names(temp)[names(temp) %in% levels(metaCombined$Study)] 
    for(dataset in CrtxDS){
      ExpAll[[dataset]] <- temp[[dataset]]$aned_high
    }
  } else if(grepl("GSE80655", study)){
    temp = eval(as.name(study))
    ExpAll[[gsub("studyFinal", "", study)]] <-  temp$DLPFC$aned_high
    } else {
    temp = eval(as.name(study))
    ExpAll[[gsub("studyFinal", "", study)]] <-  temp$Cortex$aned_high
  }
}

save(ExpAll, metaCombined, file = paste0(GeneralResultsPath, "DataCombined.rda"))

lmResultsBP <- sapply(grep("_Genes", names(metaCombined), value = TRUE), function(celltype){
  GetLM(metaCombined, celltype, "BP")
}, simplify  = FALSE)

lmResultsSCZ <- sapply(grep("_Genes", names(metaCombined), value = TRUE), function(celltype){
  GetLM(metaCombined, celltype, "SCZ")
}, simplify  = FALSE)

lmResultsCombined <- sapply(grep("_Genes", names(metaCombined), value = TRUE), function(celltype){
  GetLM(metaCombined, celltype, c("SCZ", "BP"))
}, simplify  = FALSE)

metaCombinedMelt <- melt(metaCombined, id.vars = grep("_Genes", names(metaCombined), invert = TRUE, value = TRUE),
                         variable.name = "CellType",
                         value.name = "MGP")

metaCombinedMelt$CellType <- sapply(metaCombinedMelt$CellType, function(x) gsub("_Genes|\\.", "", x))
metaCombinedMelt$Study2 <- factor(metaCombinedMelt$Study2, levels = rev(c("Study2AltarC (SMRI_C)", "Study4Chen (SMRI_C)", "GSE35978 (SMRI_A+C)",
                                                                          "Study1AltarA (SMRI_A)", "Study3Bahn (SMRI_A)", "Study5Dobrin (SMRI_A)", "Study7Kato (SMRI_A)",
                                                                          "McLean (McLean66)", "GSE21138 (VBBN)", "GSE53987 (ACOME)", "GSE80655 (PNDRC)")))


EstimateDeltaCell <- sapply(unique(metaCombinedMelt$Study), function(study){
  studyData <- sapply(metaCombinedMelt %>% filter(Study == study) %>% .$CellType %>% unique, function(celltype){
    data = metaCombinedMelt %>% filter(Study == study, CellType == celltype) %>% droplevels()
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
               Study2 = data$Study2 %>% unique,
               CellType = celltype,
               BP = BP,
               SCZ = SCZ,
               BPpval = BPpval,
               SCZpval = SCZpval)
  }, simplify = FALSE) %>% do.call(rbind, .)
}, simplify = FALSE) %>% do.call(rbind, .)

EstimateDeltaCell$BPstar <- GetStar(EstimateDeltaCell$BPpval)

EstimateDeltaCell$SCZstar <- GetStar(EstimateDeltaCell$SCZpval)

EstimateDeltaCell$Study <- factor(EstimateDeltaCell$Study, levels = c("Study2AltarC", "Study4Chen", "GSE35978",
                                                                      "Study1AltarA", "Study3Bahn", "Study5Dobrin", "Study7Kato",
                                                                      "McLean", "GSE21138", "GSE53987", "GSE80655"))
EstimateDeltaCell$Study2 <- factor(EstimateDeltaCell$Study2, levels = c("Study2AltarC (SMRI_C)", "Study4Chen (SMRI_C)", "GSE35978 (SMRI_A+C)",
                                                                        "Study1AltarA (SMRI_A)", "Study3Bahn (SMRI_A)", "Study5Dobrin (SMRI_A)", "Study7Kato (SMRI_A)",
                                                                        "McLean (McLean66)", "GSE21138 (VBBN)", "GSE53987 (ACOME)", "GSE80655 (PNDRC)"))

EstimateDeltaCell$CellType <- factor(EstimateDeltaCell$CellType, levels = c("Astrocyte" ,"Endothelial","Oligo", "OligoPrecursors",
                                                                            "Microglia", "Microglia_activation", "Microglia_deactivation",
                                                                            "GabaPV", "GabaRelnCalb",  "GabaVIPReln",
                                                                            "Layer4Pyra", "Layer6bPyra", "PyramidalAll"))

EstimateDeltaCellMelt <- melt(data.table(EstimateDeltaCell), id.vars = c("Study", "Study2", "CellType"), measure.vars = patterns("BP$|SCZ$", "pval", "star", cols = names(EstimateDeltaCell)),
                              variable.name = "Profile",
                              value.name = c("Estimate", "pVal", "Star")) %>% data.frame
EstimateDeltaCellMelt$x <- EstimateDeltaCellMelt$Profile %>% as.numeric() %>% +1
EstimateDeltaCellMelt$x[EstimateDeltaCellMelt$Study == "GSE21138"] <- 2
EstimateDeltaCellMelt$y <- 1.03

EstimateDeltaCellMelt$Profile <- sapply(EstimateDeltaCellMelt$Profile, function(x) {
  switch(as.character(x), "1"="BP", "2"="SCZ")
}) %>% factor

# Plot heatmap of estimates for all cell types in all studies
EstimatePlot <- ggplot(EstimateDeltaCellMelt, aes(CellType, Study2)) +
  theme_light(base_size = 16) +
  theme(axis.text.y = element_text(size = rel(1.2)),
        axis.text.x = element_text(size = rel(1.3),
                                   colour = c("darkred", rep("black", 6), "darkred", rep("black", 6))),
        strip.background = element_rect(fill="grey40")) +
  labs(y="", x="") +
  facet_wrap(~Profile) +
  geom_tile(aes(fill = Estimate), colour = "white") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(low = "darkolivegreen", high = "darkorange3", mid = "grey90",na.value = "white",
                       midpoint = 0,
                       limit = c(-max(abs(c(EstimateDeltaCellMelt$Estimate)), na.rm=TRUE),
                                               max(abs(c(EstimateDeltaCellMelt$Estimate)), na.rm=TRUE)),
                       name = "Effect") +
  geom_text(aes(label=Star), color="black", size=7) 
ggsave("GeneralResults/EstimatePlotPMI.pdf", plot = EstimatePlot, width = 12, height = 6, units = "in",
       dpi=300)

#Plot violin plotsof individual studies
AstroAllStudyPlot <- PlotSingleStudy(data = metaCombinedMelt, celltype = "Astrocyte", pvalSize = 5)
ggsave("GeneralResults/AstroAllStudyPlot.pdf", plot = AstroAllStudyPlot, width = 10,  height = 6, units = "in",
       dpi=300)
GabaPVAllStudyPlot <- PlotSingleStudy(data = metaCombinedMelt, celltype = "GabaPV", pvalSize = 5)
ggsave("GeneralResults/GabaPVAllStudyPlot.pdf", plot = GabaPVAllStudyPlot, width = 10,  height = 6, units = "in",
       dpi=300)


MissingCell <- EstimateDeltaCellMelt %>% filter(is.na(Estimate)) %>% group_by(CellType) %>% summarise(n = n()) %>% data.frame()
CellTypeRM2 <- MissingCell %>% filter(n/2 > 0.5*length(levels(EstimateDeltaCellMelt$Study))) %>% .$CellType

#Calculate confidence intervals
AllresultsCI <- GetCI(lmResults = lmResultsCombined, model = "lm3", CellTypes = unique(EstimateDeltaCell$CellType))
AllresultsCIb <- GetCI(lmResults = lmResultsCombined, model = "lm3b", CellTypes = unique(EstimateDeltaCell$CellType))

AllresultsCIast <-  GetCI(lmResults = lmResultsCombined, model = "lm3ast", CellTypes = unique(EstimateDeltaCell$CellType))
AllresultsCIastb <-  GetCI(lmResults = lmResultsCombined, model = "lm3bast", CellTypes = unique(EstimateDeltaCell$CellType))

AllresultsCIpv <-  GetCI(lmResults = lmResultsCombined, model = "lm3PV", CellTypes = unique(EstimateDeltaCell$CellType))
AllresultsCIpvb <-  GetCI(lmResults = lmResultsCombined, model = "lm3bPV", CellTypes = unique(EstimateDeltaCell$CellType))

#Get summary of the mixed models results
ResultsSummary <- GetResultSummary(CIdata = AllresultsCI, CellTypes = unique(EstimateDeltaCell$CellType),
                               Effects = c("ProfileSCZ","ProfileBP"), CellTypeRM2 = CellTypeRM2)
                                
ResultsSummaryAst <-GetResultSummary(CIdata = AllresultsCIast, CellTypes = unique(EstimateDeltaCell$CellType),
                                 Effects = c("ProfileSCZ","ProfileBP"), CellTypeRM2 = c(CellTypeRM2, "Astrocyte"))

ResultsSummaryPV <-GetResultSummary(CIdata = AllresultsCIpv, CellTypes = unique(EstimateDeltaCell$CellType),
                                     Effects = c("ProfileSCZ","ProfileBP"), CellTypeRM2 = c(CellTypeRM2, "GabaPV"))

ResultsSummaryFull <- GetResultSummary(CIdata = AllresultsCIb, CellTypes = unique(EstimateDeltaCell$CellType),
                               Effects = c("pH", "PMI", "Age", "SexM",  "ProfileSCZ","ProfileBP"), CellTypeRM2 = CellTypeRM2)

ResultsSummaryAstFull <-GetResultSummary(CIdata = AllresultsCIastb, CellTypes = unique(EstimateDeltaCell$CellType),
                                 Effects = c("pH", "PMI", "Age", "SexM", "ProfileSCZ","ProfileBP"), CellTypeRM2 = c(CellTypeRM2, "Astrocyte"))
ResultsSummaryPVFull <-GetResultSummary(CIdata = AllresultsCIpvb, CellTypes = unique(EstimateDeltaCell$CellType),
                                         Effects = c("pH", "PMI", "Age", "SexM", "ProfileSCZ","ProfileBP"), CellTypeRM2 = c(CellTypeRM2, "GabaPV"))

#Modify tables for plotting
for(ResultSum in ls(pat = "^ResultsSummary")){
  temp <- melt(data.table(eval(as.name(ResultSum)) %>% droplevels()), id.vars = "CellType",
               measure.vars = patterns("BP$|SCZ$|^pH$|PMI$|Age$|SexM$", "min$", "max$",
                                       cols = as.character(names(eval(as.name(ResultSum))))),
               variable.name = "Effect", variable.factor = FALSE,
               value.name = c("Coefficient", "Min", "Max")) %>% data.frame
  measureNames <- grep("min$", names(eval(as.name(ResultSum))), value = T)
  measureNames <- sapply(measureNames, function(x) gsub("min", "", x))
  for(i in 1:length(measureNames)){
    temp$Effect[temp$Effect == i] <- measureNames[i]
  } %>% factor
  temp$Effect <- relevel(as.factor(temp$Effect), ref = "BP")
  assign(paste0(ResultSum, "Melt"),temp)
}

CombinedResultsSummary <- rbind(ResultsSummaryMelt, ResultsSummaryAstMelt, ResultsSummaryPVMelt)
CombinedResultsSummary$Adjustment <- c(rep("No adjustment", 26), rep("Astrocyte adjusted", 26), rep("Gaba PV adjusted", 26)) %>%
  factor(levels = c("No adjustment","Astrocyte adjusted","Gaba PV adjusted"))

# lmer CI plot for all cell types
EstimateCIPlot <- plotCI(data = ResultsSummaryMelt %>% filter(Effect %in% c("BP", "SCZ")), title = "Group effect (95%CI)")
ggsave("GeneralResults/EstimateCIPlot.pdf", plot = EstimateCIPlot, width = 10, height = 6, units = "in",
       dpi=300, useDingbats = FALSE)

# lmer CI plot for all cell types, adjusted for astrocyte MGP
EstimateCIastPlot <- plotCI(data = ResultsSummaryAstMelt %>% filter(Effect %in% c("BP", "SCZ")), title = "Group effect (95%CI) - Astrocyte adjusted")
ggsave("GeneralResults/EstimateCIastPlot.pdf", plot = EstimateCIastPlot, width = 10, height = 6, units = "in",
       dpi=300, useDingbats = FALSE)

# lmer CI plot for all cell types, adjusted for GabaPV MGP
EstimateCIpvPlot <- plotCI(data = ResultsSummaryPVMelt %>% filter(Effect %in% c("BP", "SCZ")), title = "Group effect (95%CI) - GabaPV adjusted")
ggsave("GeneralResults/EstimateCIpvPlot.pdf", plot = EstimateCIpvPlot, width = 10, height = 6, units = "in",
       dpi=300, useDingbats = FALSE)

# lmer CI plot for all cell types, covariate effects
EstimateCIcombinedPlot <- ggplot(CombinedResultsSummary, aes(x=CellType, y=Coefficient)) +
  theme_light(base_size = 12) +
  labs(title = "Group effect (95%CI)", y = "Coefficient", x="") +
  theme(axis.text.y = element_text(size = rel(1.2)),
        axis.text.x = element_text(size = rel(1.2),
                                   colour = c("darkred", rep("black", 6), "darkred", rep("black", 6))),
        strip.background = element_rect(fill="grey40"),
        panel.grid.major.y =  element_blank(),
        panel.grid.minor.y =  element_blank()) +
  facet_grid(Effect~Adjustment) +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust=0.5)) +
  geom_point() +
  geom_errorbar(aes(ymin = Min, ymax = Max), na.rm = TRUE) +
  geom_abline(intercept = 0, slope = 0, color = "red", linetype="dashed")
ggsave("GeneralResults/EstimateCIcombinedPlotWithPMI.pdf", plot = EstimateCIcombinedPlot, width = 10, height = 6, units = "in",
       dpi=300, useDingbats = FALSE)

EstimateCIFullPlot <- plotCI(data = ResultsSummaryFullMelt %>% filter(Effect %in% c("pH", "PMI", "Age", "SexM")), title = "Covariate effect (95%CI)",
                             rowNum = 2, scales = "free_y", color_x = FALSE)
ggsave("GeneralResults/EstimateCIFullPlotPMI.pdf", plot = EstimateCIFullPlot, width = 10, height = 10, units = "in",
       dpi=300, useDingbats = FALSE)

# lmer CI plot for Astrocyte MGP only
EstimateCIAstOnlyPlot <- ggplot(ResultsSummaryMelt %>% filter(CellType %in% c("Astrocyte")), aes(x=Effect, y=Coefficient)) +
  theme_classic(base_size = 16) +
  labs(title = "Mixed effects model Group effect (95% CI)", y = "Coefficient", x="") +
  theme(axis.text.y = element_text(size = rel(1.4)),
        axis.text.x = element_text(size = rel(1.4)),
        strip.background = element_rect(fill="grey40"),
        panel.grid.major.y =  element_blank(),
        panel.grid.minor.y =  element_blank(),
        panel.grid.major.x =  element_blank(),
        panel.grid.minor.x =  element_blank()) +
  scale_y_continuous(limits = c(-max(ResultsSummaryMelt %>% filter(CellType %in% c("Astrocyte")) %>% .$Max), NA)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Min, ymax = Max), na.rm = TRUE, width = 0.6) +
  geom_abline(intercept = 0, slope = 0, color = "red", linetype="dashed")
ggsave("GeneralResults/EstimateCIAstOnlyPlot.pdf", plot = EstimateCIAstOnlyPlot, width = 6, height = 3, units = "in",
       dpi=300, useDingbats = FALSE)

# lmer CI plot for GabaPV MGP only
EstimateCIGabaPVOnlyPlot <- ggplot(ResultsSummaryMelt %>% filter(CellType %in% c("GabaPV")), aes(x=Effect, y=Coefficient)) +
  theme_classic(base_size = 16) +
  labs(title = "Mixed effects model Group effect (95% CI)", y = "Coefficient", x="") +
  theme(axis.text.y = element_text(size = rel(1.4)),
        axis.text.x = element_text(size = rel(1.4)),
        strip.background = element_rect(fill="grey40"),
        panel.grid.major.y =  element_blank(),
        panel.grid.minor.y =  element_blank(),
        panel.grid.major.x =  element_blank(),
        panel.grid.minor.x =  element_blank()) +
  scale_y_continuous(limits = c(NA, -min(ResultsSummaryMelt %>% filter(CellType %in% c("GabaPV")) %>% .$Min))) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Min, ymax = Max), na.rm = TRUE, width = 0.6) +
  geom_abline(intercept = 0, slope = 0, color = "red", linetype="dashed")
ggsave("GeneralResults/EstimateCIGabaPVOnlyPlot.pdf", plot = EstimateCIGabaPVOnlyPlot, width = 6, height = 3, units = "in",
       dpi=300, useDingbats = FALSE)

ResultTableCovariates <- ResultsSummaryFull %>% select(matches("CellType|pH|PMI|Age|SexM")) %>% .[complete.cases(.),]
ResultTableCovariates[-1] <- apply(ResultTableCovariates[-1], c(1,2), function(x) signif(x,digits=1)) 
ResultTableCovariates %<>% mutate("pH (95%CI)" = paste0(pH, " (", pHmin, ",", pHmax, ")"),
                                  "Age (95%CI)" = paste0(Age, " (", Agemin, ",", Agemax, ")"),
                                  "SexM (95%CI)" = paste0(pH, " (", SexMmin, ",", SexMmax, ")"),
                                  "pH (95%CI)" = paste0(pH, " (", pHmin, ",", pHmax, ")"),
                                  "Age (95%CI)" = paste0(Age, " (", Agemin, ",", Agemax, ")"),
                                  "PMI (95%CI)" = paste0(PMI, " (", PMImin, ",", PMImax, ")"))

ResultTableCovariates %<>% select(matches("CellType|pH |PMI |Age |SexM ")) %>% droplevels
ResultTableCovariates <- ResultTableCovariates[match(levels(ResultTableCovariates$CellType),
                                                     ResultTableCovariates$CellType),] 
write.table(ResultTableCovariates, "GeneralResults/ResultTableCovariatesPMI.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

save.image("GeneralResults/CombinedStudies.Rdata")

#lm including RIN for GSE53987 for astrocyte and GabaPV cells
dataGSE53987 <- metaCombined %>% filter(Study == "GSE53987") %>% droplevels
lm1Astro <- lm(Astrocyte_Genes~Profile + pH + Age + Sex + PMI + RIN, data = dataGSE53987) %>% summary
lm2Astro <- lm(Astrocyte_Genes~Profile + pH + Age + Sex + PMI, data = dataGSE53987) %>% summary
lm1GabaPV <- lm(GabaPV_Genes~Profile + pH + Age + Sex + PMI + RIN, data = dataGSE53987) %>% summary
lm2GabaPV <- lm(GabaPV_Genes~Profile + pH + Age + Sex + PMI, data = dataGSE53987) %>% summary

write.table(lm1Astro$coefficients, file = paste0(GeneralResultsPath, "/lmGSE53987AstroMod1.tsv"), sep = "\t",
            row.names = TRUE, col.names = TRUE)
write.table(lm2Astro$coefficients, file = paste0(GeneralResultsPath, "/lmGSE53987AstroMod2.tsv"), sep = "\t",
            row.names = TRUE, col.names = TRUE)
write.table(lm1GabaPV$coefficients, file = paste0(GeneralResultsPath, "/lmGSE53987GabaPVMod1.tsv"), sep = "\t",
            row.names = TRUE, col.names = TRUE)
write.table(lm2GabaPV$coefficients, file = paste0(GeneralResultsPath, "/lmGSE53987GabaPVMod2.tsv"), sep = "\t",
            row.names = TRUE, col.names = TRUE)
