source("SetUp.R")
packageF("parallel")
packageF("ggdendro")
packageF("cowplot")
packageF("pheatmap")


load(paste0(GeneralResultsPath, "DataCombined.rda"))

metaCombined$SubjStudy %>% head()

AllGenes <- lapply(ExpAll, function(study){
  as.character(study$GeneSymbol)
}) %>% unlist

CommonGenes <- names(table(AllGenes))[table(AllGenes) >= 5]
CommonGenesList <- as.list(CommonGenes)
names(CommonGenesList) <- CommonGenes



GetMMgene <- function(geneList = CommonGenesList, coreNum){
  MMgenes <- mclapply(geneList, function(gene){
    temp <- MetaExp[!is.na(MetaExp[[gene]]),] %>% .$SubjectID %>% unlist
    if(length(unique(temp)) < length(temp)){
      mod1 = "Sex + pH + Age + PMI + (1|SubjectID) + (1|Study)"
      mod2 = "Sex + pH + Age + PMI + GabaPV_Genes + Astrocyte_Genes + Endothelial_Genes + GabaVIPReln_Genes + Microglia_Genes + Microglia_activation_Genes + Microglia_deactivation_Genes + Oligo_Genes + OligoPrecursors_Genes + PyramidalAll_Genes + (1|SubjectID) + (1|Study)"
    } else {
      mod1 = "Sex + pH + Age + PMI + (1|Study)"
      mod2 = "Sex + pH + Age + PMI + GabaPV_Genes + Astrocyte_Genes + Endothelial_Genes + GabaVIPReln_Genes + Microglia_Genes + Microglia_activation_Genes + Microglia_deactivation_Genes + Oligo_Genes + OligoPrecursors_Genes + PyramidalAll_Genes + (1|Study)"
    }
    gene2 = paste0("`", gene, "`")
    BPmod1null <- lmer(as.formula(paste0(gene2, "~", mod1))  , data = MetaExpBP, REML = FALSE)
    BPmod1Full <- lmer(as.formula(paste0(gene2, "~ Profile +" , mod1)), data = MetaExpBP, REML = FALSE)
    BPmod2null <- lmer(as.formula(paste0(gene2, "~", mod2))  , data = MetaExpBP, REML = FALSE)
    BPmod2Full <- lmer(as.formula(paste0(gene2, "~ Profile +" , mod2)), data = MetaExpBP, REML = FALSE)
    
    temp <- MetaExpSCZ[!is.na(MetaExpSCZ[[gene]]),] %>% .$SubjectID %>% unlist
    if(length(unique(temp)) < length(temp)){
      mod1 = "Sex + pH + Age + PMI + (1|SubjectID) + (1|Study)"
      mod2 = "Sex + pH + Age + PMI + GabaPV_Genes + Astrocyte_Genes + Endothelial_Genes + GabaVIPReln_Genes + Microglia_Genes + Microglia_activation_Genes + Microglia_deactivation_Genes + Oligo_Genes + OligoPrecursors_Genes + PyramidalAll_Genes + (1|SubjectID) + (1|Study)"
    } else {
      mod1 = "Sex + pH + Age + PMI + (1|Study)"
      mod2 = "Sex + pH + Age + PMI + GabaPV_Genes + Astrocyte_Genes + Endothelial_Genes + GabaVIPReln_Genes + Microglia_Genes + Microglia_activation_Genes + Microglia_deactivation_Genes + Oligo_Genes + OligoPrecursors_Genes + PyramidalAll_Genes + (1|Study)"
    }
    
    SCZmod1null <- lmer(as.formula(paste0(gene2, "~", mod1))  , data = MetaExpSCZ, REML = FALSE)
    SCZmod1Full <- lmer(as.formula(paste0(gene2, "~ Profile +" , mod1)), data = MetaExpSCZ, REML = FALSE)
    SCZmod2null <- lmer(as.formula(paste0(gene2, "~", mod2))  , data = MetaExpSCZ, REML = FALSE)
    SCZmod2Full <- lmer(as.formula(paste0(gene2, "~ Profile +" , mod2)), data = MetaExpSCZ, REML = FALSE)
    
    list(BPpvalMod1 = anova(BPmod1null, BPmod1Full), BPfullMod1 = BPmod1Full,
         BPpvalMod2 = anova(BPmod2null, BPmod2Full), BPfullMod2 = BPmod2Full,
         SCZpvalMod1 = anova(SCZmod1null, SCZmod1Full), SCZfullMod1 = SCZmod1Full,
         SCZpvalMod2 = anova(SCZmod2null, SCZmod2Full), SCZfullMod2 = SCZmod2Full)
  }, mc.cores = coreNum)
  return(MMgenes)
}

GetStatTable <- function(group, model){
  Stat <- mclapply(MMgenes, function(geneLM){
    temp <- geneLM[[paste0(group, "pval", model)]] %>% data.frame
    temp2 <- geneLM[[paste0(group,"full", model)]] %>% summary %>% .$coefficients
    c(temp2[2,1], temp$Pr[2])
  }, mc.cores = detectCores()/6) %>% do.call(rbind, .) %>% data.frame()
  names(Stat) <- c("FC", "pVal")
  StatAdj <- data.frame(GeneSymbol = rownames(Stat), FC = 2^(Stat$FC),
                        pVal  = Stat$pVal, pAdj = p.adjust(Stat$pVal, "BH"))
  StatAdj$pValDown <- apply(StatAdj %>% select(FC, pVal), 1, function(gene){
    if(gene[1] < 1){
      gene[2]/2
    } else {
      1-gene[2]/2
    }
  })
  StatAdj$pAdjDown <- p.adjust(StatAdj$pValDown, "BH")
  return(StatAdj)
}

GeneStudyPlot <- function(data, gene, title, ptSize = 1){
  ggplot(data, aes_string("Profile", gene, fill = "Profile")) +
    theme_classic() +
    labs(title = title) +
    geom_boxplot(outlier.size = NULL) + geom_jitter(width = 0.2, size = ptSize) +
    facet_wrap(~Study, scales = "free_y")
}

GetAdjustedExpr <- function(data, gene, adjVar){
  data <- data[c("CommonName", "SubjectID", "Study", "SubjStudy", "Profile", "Age", "Sex", "pH", "PMI",  gene, "Astrocyte_Genes", "GabaPV_Genes")]
  names(data)[names(data) == gene] <- "Gene"
  data <- data[!is.na(data$Gene),] %>% droplevels()
  LM <- sapply(unique(data$Study), function(study){
    subData = data[data$Study == study,]
    subData <- subData[!is.na(subData$pH),]
    Mod <- lm(as.formula(paste0("Gene~", paste0(c("Profile",adjVar), collapse = "+"))), data = subData)
    modFrame <- model.frame(Mod)
    modMatrix <- model.matrix.lm(Mod)
    Coef <- coef(Mod)
    for(Var in adjVar){
      if(Var %in% colnames(modMatrix)){
        modMatrix[,grepl(Var, colnames(modMatrix))] <- mean(modMatrix[,grepl(Var, colnames(modMatrix))])
      } else {
        Var = grep(Var, colnames(modMatrix), value = T)
        modMatrix[,grepl(Var, colnames(modMatrix))] <- rep(0, nrow(modMatrix))
      }
    }
    Resid <- resid(Mod)
    AdjExpr <- modMatrix %*% Coef + Resid
    modFrame <- cbind(AdjExpr, modFrame, study)
    colnames(modFrame)[1] <- "AdjExpr"
    colnames(modFrame)[grepl("study", colnames(modFrame))] <- "Study"
    modFrame$SubjectID <- subData$SubjectID
    modFrame$SubjStudy <- subData$SubjStudy
    modFrame
  }, simplify = FALSE)
  names(LM) <- unique(data$Study)
  return(LM)
}

GeneStudyPlotCombine <- function(data = MetaExpSCZ, gene, ptSize = 0.5){
  LM0 <- GetAdjustedExpr(data, gene = gene, adjVar = c("Age", "Sex", "pH")) %>% rbindlist() %>% data.frame()
  LM1 <- GetAdjustedExpr(data, gene = gene, adjVar = c("Age", "Sex", "pH", "Astrocyte_Genes", "GabaPV_Genes")) %>% rbindlist() %>% data.frame()
  p1 <- GeneStudyPlot(LM0, gene = "AdjExpr", title = paste0(gene, " - Sex, pH, Age adjusted"), ptSize = ptSize)
  p2 <- GeneStudyPlot(LM1, gene = "AdjExpr", title = paste0(gene, " - Sex, pH, Age, Astrocyte, GabaPV adjusted"), ptSize = ptSize)
  return(list(LM0plot = p1, LM1plot = p2))
}

ScaleAdjExpr <- function(genes, data = NULL, group = NULL, adjVar = c("Sex", "Age", "pH", "PMI")){
  if(is.null(data)){
    if(group == "SCZ"){
      groups = c("Cont", "SCZ")
      dataset = "Study4Chen"
    } else {
      groups = c("Cont", "BP")
      dataset = "GSE21138"
    }
    meta <- metaCombined %>% filter(!is.na(Astrocyte_Genes), Profile %in% groups, Study != dataset) %>% droplevels()
    GeneList = as.list(genes)
    names(GeneList) <- genes
    Expr <- mclapply(GeneList, function(gene){
      sapply(meta$SubjStudy, function(subj){
        study = meta[meta$SubjStudy == subj,] %>% .$Study %>% as.character()
        ID = meta[meta$SubjStudy == subj,] %>% .$CommonName %>% as.character()
        if(gene %in% as.character(ExpAll[[study]]$GeneSymbol)) {
          ExpAll[[study]] %>% .[.$GeneSymbol == gene,] %>% .[[ID]]
        } else {
          NA
        }
      })
    }, mc.cores = 0.8*detectCores()) %>% do.call(cbind, .)
    data <- cbind(meta, Expr)
  }
  temp <- sapply(genes, function(gene){
    AdjExp <- GetAdjustedExpr(data = data, gene = gene, adjVar = adjVar)
    AdjExp2 <- lapply(AdjExp, function(study){
      study$GeneSymbol <- as.factor(gene)
      study %>% mutate(AdjExpScaled = scale(AdjExpr),
                       NonAdjScaled = scale(Gene))
    })
    AdjExp2 %>% rbindlist() %>% data.frame()
  }, simplify = F) %>% rbindlist
  return(temp)
}

AddMGPdata <- function(AdjData, Meta, celltype){
  MGPdf <- data.frame(AdjExpr = Meta[[celltype]],
                      Gene = Meta[[celltype]])
  MGPdf <- cbind(MGPdf, Meta %>% select(Profile, Sex, Age, pH, Study, PMI, SubjectID, SubjStudy))
  MGPdf$GeneSymbol <- celltype
  MGPdf$AdjExpScaled <- scale(MGPdf$AdjExpr) %>% as.numeric
  MGPdf$NonAdjScaled <- scale(MGPdf$Gene) %>% as.numeric
  AdjData <- rbind(MGPdf, AdjData)
  AdjData$Type <- sapply(AdjData$GeneSymbol, function(x){
    if(grepl("Genes", x)){
      "MGP"
    } else {
      "Gene"
    } 
  })
  return(AdjData)
}

AdjScaleHeatmap <- function(AdjData, Meta, title, xlab = "Subjects", ylab = "Genes",
                             ExpCol = AdjExpScaled,
                             GrpCount = GroupNum, geneType = NULL, celltype = NULL,
                             Ysize = NULL, KnownMarkers = c("ALDH1L1", "SOX9", "AQP4", "SLC1A3"),
                             GeneCol = c("cornflowerblue", "gold3", "coral2", "darkolivegreen4", "grey20")){
  AdjData$Study <- factor(AdjData$Study, levels = names(GrpCount))
  GeneNum = sapply(levels(AdjData$Study), function(study){
    nrow(PCAresults[[study]][[1]]$All[[celltype]]$rotation)
  }) %>% unlist
  
  GrpcolBar <- data.frame(Study = names(GrpCount),
                          xmin1 = 0.5,
                          xmax1 = lapply(GrpCount, function(x) x$n[1]+0.5) %>% unlist,
                          xmin2 = lapply(GrpCount, function(x) x$n[1]+0.5) %>% unlist,
                          xmax2 = lapply(GrpCount, function(x) sum(x$n[1] + x$n[2])+0.5) %>% unlist,
                          ymin = -GeneNum/10,
                          ymax = 0)
  plots <- list()
  for(study in levels(AdjData$Study)){
    
    AdjDataSub <- AdjData %>% filter(Study == study) %>% droplevels()
    MetaSub <- Meta %>% filter(Study == study) %>% droplevels()
    #Filter and arrange genes and subjects
    #Subject order
    CellGenes <- names(PCAresults[[study]][[1]]$All[[celltype]]$rotation[,1] %>% sort)
    AdjDataSub <- AdjDataSub[AdjDataSub$GeneSymbol %in% c(celltype, CellGenes),] %>% droplevels()
    
    SubjOrder <- MetaSub %>% arrange_(.dots = c("Profile", geneType)) %>% .$SubjStudy %>% as.character()
    AdjDataSub$SubjStudy <- factor(AdjDataSub$SubjStudy, levels = unique(SubjOrder))
    
    #Gene order
    ScaleExpDF <- data.frame(SubjStudy = unique(AdjDataSub$SubjStudy), Study = study)
    for(gene in unique(AdjDataSub$GeneSymbol)){
      temp <- AdjDataSub[AdjDataSub$GeneSymbol == gene,] %>% select(SubjStudy, AdjExpScaled) %>% data.frame() %>% droplevels()
      names(temp)[2] <- gene
      ScaleExpDF <- merge(ScaleExpDF, temp, by = "SubjStudy")
    }

    AdjDataSub$GeneSymbol <- factor(AdjDataSub$GeneSymbol, levels = CellGenes)
    
    if(is.null(Ysize)){
      Ysize = 0.2
    }
    
    GrpcolBarSub <- GrpcolBar %>% filter(Study == study) %>% droplevels()
    
    #Creatte the boxplot for marker gene expression
    temp <- AdjDataSub %>% filter(Type == "Gene") %>% droplevels
    GrpEffect <- temp %>% group_by(GeneSymbol,Profile) %>%
      summarise(MeanExp = mean(NonAdjScaled), MedianExp = median(NonAdjScaled)) %>% data.frame() 
    
    GrpEffect$GeneType <- sapply(GrpEffect$GeneSymbol, function(gene){
      if(gene %in% KnownMarkers){
        as.character(gene)
      } else {
        "Other"
      }
    }) %>% factor
    GrpEffect$GeneType2 <- sapply(GrpEffect$GeneSymbol, function(gene){
      if(gene %in% KnownMarkers){
        "KnownMarker"
      } else {
        "Other"
      }
    }) %>% factor
    
    
    GeneColors <- GeneCol
    names(GeneColors) <- c(KnownMarkers, "Other")
    GrpEffect$GeneType <- factor(GrpEffect$GeneType, levels = c(KnownMarkers, "Other"))

    p1 <- ggplot(GrpEffect, aes(Profile, MedianExp)) + geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes(color = GeneType,
                      alpha = GeneType2,
                      size = GeneType2 )) +
      theme(axis.title.y = element_text(size=12),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            legend.position = c(0, 0.85), legend.justification=c(0,0),
            legend.text=element_text(size=6), legend.direction = "horizontal",
            plot.margin = unit(c(0,0,-0.6,0), "cm")) +
      labs(x = "", y = "Median scaled\nexpression") +
      scale_color_manual(values=GeneColors, name = "") +
      scale_alpha_discrete(range=c(1, 0.3), guide = "none") +
      scale_size_discrete(range=c(1.5,0.6),  guide = "none")
    
    p2 <- ggplot(AdjDataSub %>% filter(Type == "MGP") %>% droplevels() , aes(SubjStudy, GeneSymbol)) +
      theme_minimal()+
      geom_tile(aes_string(fill = ExpCol)) +
      scale_fill_gradient2(low = "darkslateblue", mid = "white", high = "gold1", midpoint = 0, limits = c(-2.5,2.5)) +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(angle = 0, size=12, vjust=0.5, hjust = 0.9),
            axis.text.y = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(-0.3,0,0,0), "cm")) +
      labs(title = "", x=xlab, y="MGP")
    
    p3 <- ggplot(AdjDataSub %>% filter(Type == "Gene") %>% droplevels() , aes(SubjStudy, GeneSymbol)) +
      theme_minimal()+
      geom_tile(aes_string(fill = ExpCol)) +
      scale_fill_gradient2(low = "darkslateblue", mid = "white", high = "gold1", midpoint = 0, limits = c(-2.5,2.5)) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = 12),
            axis.text.y = element_text(size=Ysize),
            legend.position = "none",
            plot.margin = unit(c(-0.3,0,0,0), "cm")) +
      labs(title="", x=xlab, y="Marker\ngenes") +
      geom_rect(data = GrpcolBarSub, aes(xmin = GrpcolBarSub$xmin2,
                                         xmax =  GrpcolBarSub$xmax2,
                                         ymin = GrpcolBarSub$ymin,
                                         ymax = GrpcolBarSub$ymax),
                fill = "brown2", inherit.aes = FALSE) +
      geom_rect(data = GrpcolBarSub, aes(xmin = GrpcolBarSub$xmin1,
                                         xmax =  GrpcolBarSub$xmax1,
                                         ymin = GrpcolBarSub$ymin,
                                         ymax = GrpcolBarSub$ymax),
                fill = "cyan4", inherit.aes = FALSE)
    plots[[study]] <- plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(4,1,5))
  }
  plot_grid(plotlist = plots, nrow = 3, labels = names(plots), vjust = 0.5, hjust = -1, scale = 0.85,
            label_size = 12)
}

  
if(length(list.files(path = GeneralResultsPath, pattern = "MetaExprCombined.rds")) == 0){
  meta <- metaCombined %>% filter(!is.na(Astrocyte_Genes))
  
  Expr <- mclapply(CommonGenesList, function(gene){
    sapply(meta$SubjStudy, function(subj){
      study = meta[meta$SubjStudy == subj,] %>% .$Study %>% as.character()
      ID = meta[meta$SubjStudy == subj,] %>% .$CommonName %>% as.character()
      if(gene %in% as.character(ExpAll[[study]]$GeneSymbol)) {
        ExpAll[[study]] %>% .[.$GeneSymbol == gene,] %>% .[[ID]]
      } else {
        NA
      }
    })
  }, mc.cores = 0.8*detectCores()) %>% do.call(cbind, .)
  
  MetaExp <- cbind(meta, Expr)
  MetaExp$Subj <-  paste(MetaExp$Profile, 1:nrow(MetaExp), sep="_")
  saveRDS(MetaExp, file = paste0(GeneralResultsPath, "MetaExprCombined.rds"))
} else {
  MetaExp <- readRDS( paste0(GeneralResultsPath, "MetaExprCombined.rds"))
}

MetaExpBP <- MetaExp %>% filter(Profile != "SCZ", Study != "GSE21138") %>% droplevels()
MetaExpSCZ <- MetaExp %>% filter(Profile != "BP", Study != "Study4Chen") %>% droplevels()


MMgenes <- GetMMgene(geneList = CommonGenesList, coreNum = 0.8*detectCores())

BPpval <- GetStatTable(group = "BP", model = "Mod1")
BPpvalMGP <- GetStatTable(group = "BP", model = "Mod2")
SCZpval <- GetStatTable(group = "SCZ", model = "Mod1")
SCZpvalMGP <- GetStatTable(group = "SCZ", model = "Mod2")

BPcombine <- merge(BPpval, BPpvalMGP, by="GeneSymbol", suffixes = c("_mod0", "_mod1"))
BPcombine %<>% mutate(pAdjLog_mod0 = -log10(pAdj_mod0),
                      pAdjLog_mod1 = -log10(pAdj_mod1))
SCZcombine <- merge(SCZpval, SCZpvalMGP, by="GeneSymbol", suffixes = c("_mod0", "_mod1"))
SCZcombine %<>% mutate(pAdjLog_mod0 = -log10(pAdj_mod0),
                       pAdjLog_mod1 = -log10(pAdj_mod1))

rm(MMgenes)
.rs.restartR()

#Comparing to marker genes and Mistry et al genes
source(paste0(GenScriptPath, "Cell_type_PCA.R"))
HumanMarkers <- list()
HumanMarkers$Astrocyte <- sapply(mouseMarkerGenesCombined$Cortex$Astrocyte, function(gene){
  mouse2human(gene)$humanGene
}) %>% unlist

HumanMarkers$GabaPV <- sapply(mouseMarkerGenesCombined$Cortex$GabaPV, function(gene){
  mouse2human(gene)$humanGene
}) %>% unlist

BPcombine$Marker <- sapply(BPcombine$GeneSymbol, function(gene){
  if(gene %in% HumanMarkers$Astrocyte){
    "Astrocyte"
  } else if(gene %in% HumanMarkers$GabaPV){
    "GabaPV"
  } else{
    NA
  }
})
save(BPcombine, file = paste(GeneralResultsPath, "MMgenesBP_allMGP.rda"))
write.table(BPcombine, file = "BPcombine_allMGP.txt", row.names = FALSE, sep = "\t")

SCZcombine$Marker <- sapply(SCZcombine$GeneSymbol, function(gene){
  if(gene %in% HumanMarkers$Astrocyte){
    "Astrocyte"
  } else if(gene %in% HumanMarkers$GabaPV){
    "GabaPV"
  } else{
    NA
  }
})

MistryUP <- read.table("MistryUp.txt", header = T, sep = "\t")
MistryDown <- read.table("MistryDown.txt", header = T, sep = "\t")

#Update Mistry et al gene names
Anno_file <- read.table("GPL570.gz", comment="#", header=T, quote='"', sep="\t")

MistryDown$GeneSymbolOrg <- MistryDown$GeneSymbol
MistryDown$GeneSymbol <- sapply(MistryDown$Probe, function(x) Anno_file %>% filter(ProbeName == as.character(x)) %>%
                                  .$GeneSymbol) %>% make.names
names(MistryDown) <- sapply(names(MistryDown), function(x) gsub("Down", "", x))

MistryUP$GeneSymbolOrg <- MistryUP$GeneSymbol
MistryUP$GeneSymbol <- sapply(MistryUP$Probe, function(x) Anno_file %>% filter(ProbeName == as.character(x)) %>%
                                .$GeneSymbol) %>% make.names
names(MistryUP) <- sapply(names(MistryUP), function(x) gsub("Up", "", x))

MistryCombined <- rbind(MistryDown, MistryUP)

SCZcombine$MistryFC <- MistryCombined$FoldChange[match(SCZcombine$GeneSymbol, MistryCombined$GeneSymbol)]
SCZcombine$MistryqVal <- MistryCombined$Qvalue[match(SCZcombine$GeneSymbol, MistryCombined$GeneSymbol)]
save(SCZcombine, file = paste(GeneralResultsPath, "MMgenesSCZ_allMGP.rda"))
write.table(SCZcombine, file = paste0(GeneralResultsPath, "SCZcombine_allMGP.txt"), row.names = FALSE, sep = "\t")
write.table(SCZcombine %>% filter(FC_mod0 < 1), file = paste0(GeneralResultsPath, "SCZcombineDownMod0_allMGP.txt"), row.names = FALSE, sep = "\t", col.names = TRUE)
write.table(SCZcombine %>% filter(FC_mod1 < 1), file = paste0(GeneralResultsPath, "SCZcombineDownMod1_allMGP.txt"), row.names = FALSE, sep = "\t", col.names = TRUE)

SCZmistry <- SCZcombine[!is.na(SCZcombine$MistryFC),]
SCZmistryMelt <- melt(data.table(SCZmistry),
                      id = grep("GeneSymbol|Mistry|Marker", names(SCZmistry)),
                      measure=patterns("FC_mod", "pVal_mod", "pAdj_", "pAdjLog", "pValDown", "pAdjDown"),
                      value.name = c("FC", "pVal", "pAdj", "pAdjLog", "pValDown", "pAdjDown"),
                      variable.name = "Model",
                      variable.factor = TRUE) %>% data.frame
levels(SCZmistryMelt$Model) <- c("no_MGPadj", "MGPadj")

corMod0 <- cor.test(~MistryFC+FC, data = SCZmistryMelt %>% filter(Model == "no_MGPadj"), method = "spearman")
corMod1 <- cor.test(~MistryFC+FC, data = SCZmistryMelt %>% filter(Model == "MGPadj"), method = "spearman")

labelmod0 <- paste0('no_MGPadj~r[s] == ', signif(corMod0$estimate, digits = 3))
labelmod1 <- paste0("MGPadjr[s] == ", signif(corMod1$estimate, digits = 3))

#Plotting correlation between FC in Mistry et all and current analysis with and without correcting for MGPs
ggplot(SCZmistryMelt, aes(MistryFC,FC, color = Model)) +
  theme_classic(base_size = 16) +
  scale_color_manual(values=c("darkorange", "darkgreen")) +
  labs(x = "Fold change Mistry et al", y = "Fold change") +
  geom_point(size=0.8) + stat_smooth(method=lm) +
  annotate("text", x = 0.885, y = 1.2,  label = labelmod0, parse=T, color = "darkorange", size = 5) +
  annotate("text", x = 1.02, y = 1.2 , label = labelmod1, parse=T, color = "darkgreen", size = 5)
ggsave("MistryUndCurrentFCcor.pdf", path = GeneralResultsPath, width = 8, height = 6, units = "in", dpi=300)


GeneStudyPlot(MetaExpSCZ, gene = "B4GAT1", title = "Not adjusted")

#Show adjusted expression for NECAB3 that was given as an example in Mistry et al
temp <- GeneStudyPlotCombine(gene = "NECAB3")
plot(temp$LM0plot)
plot(temp$LM1plot)

#Create Heatmaps
load(paste0(GeneralResultsPath, "GeneMGPcor.rda"))
load(paste0(GeneralResultsPath, "PCAresultsList.rda"))
names(PCAresults) <- sapply(names(PCAresults), function(x) {gsub("PCAresults|Cortex", "", x)})

downGenes <- SCZmistry %>% filter(MistryFC < 1) %>% .$GeneSymbol %>% as.character()

MistryDownAdjO <- ScaleAdjExpr(genes = downGenes, data = MetaExpSCZ, group = "SCZ") %>% data.frame()

AstroGenesAdj <- ScaleAdjExpr(genes = HumanMarkers$Astrocyte, group = "SCZ") %>% data.frame()
GabaPVGenesAdj <- ScaleAdjExpr(genes = HumanMarkers$GabaPV, group = "SCZ") %>% data.frame()

#Add the MGP to AdjData
AstroGenesAdj <- AddMGPdata(AdjData = AstroGenesAdj, Meta = MetaExpSCZ, celltype = "Astrocyte_Genes")
GabaPVGenesAdj <- AddMGPdata(AdjData = GabaPVGenesAdj, Meta = MetaExpSCZ, celltype = "GabaPV_Genes")

SCZGrpEffectAstro <- sapply(levels(AstroGenesAdj$Study), function(study){
  SubData <- AstroGenesAdj %>% filter(Study == study, Type == "Gene") %>% droplevels
  SubData %>% group_by(GeneSymbol,Profile) %>% summarise(MeanExp = mean(NonAdjScaled), MedianExp = median(NonAdjScaled)) %>%
    data.frame() %>% mutate(Study = study)
}, simplify = F) %>% rbindlist %>% data.frame

SCZGrpEffectAstro$GeneType <- sapply(SCZGrpEffectAstro$GeneSymbol, function(gene){
  if(gene %in% c("ALDH1L1", "SOX9", "AQP4", "SLC1A3")){
    gene
  } else {
    "Other"
  }
}) %>% factor
SCZGrpEffectAstro$GeneType2 <- sapply(SCZGrpEffectAstro$GeneSymbol, function(gene){
  if(gene %in% c("ALDH1L1", "SOX9", "AQP4", "SLC1A3")){
    "KnownMarker"
  } else {
    "Other"
  }
}) %>% factor

SCZGrpEffectAstro$GeneType <- factor(SCZGrpEffectAstro$GeneType, levels = c("ALDH1L1", "SOX9", "AQP4", "SLC1A3", "Other"))

GeneCol <- c("cornflowerblue", "gold3", "coral2", "darkolivegreen4", "grey")
ggplot(SCZGrpEffectAstro, aes(Profile, MedianExp)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = GeneType,
                  alpha = GeneType2,
                  size = GeneType2 )) +
  theme(legend.position = "top", legend.direction = "horizontal") +
  scale_color_manual(values=GeneCol, name = "") +
  scale_alpha_discrete(range=c(1, 0.3), guide = "none") +
  scale_size_discrete(range=c(1.5,0.6),  guide = "none") +
  facet_wrap(~Study, nrow = 3)

GroupNum <- sapply(unique(MetaExpSCZ$Study), function(study){
  meta <- MetaExpSCZ %>% select(Study, Profile) 
  group_by(meta %>% filter(Study == study), Profile) %>% summarise(n = n()) %>% data.frame()
}, simplify = F)
names(GroupNum) <- unique(MetaExpSCZ$Study)

#Plotting with gene boxplots and MGP
pdf(paste0(GeneralResultsPath, "/GabaPVGenesNotAdjusted2.pdf"), width = 12, height = 8,  pointsize = 12, useDingbats = FALSE)
AdjScaleHeatmap2(GabaPVGenesAdj, Meta = MetaExpSCZ, title = "GabaPVGenesAdj, not Adjusted",ExpCol = "NonAdjScaled",
                 geneType = "GabaPV_Genes", celltype = "GabaPV_Genes", KnownMarkers = c("TAC1", "PVALB"),
                 GeneCol = c("cornflowerblue", "gold3", "grey20"))
dev.off()

pdf(paste0(GeneralResultsPath, "/AstrocyteGenesNotAdjusted2.pdf"), width = 12, height = 8,  pointsize = 12, useDingbats = FALSE)
AdjScaleHeatmap(AstroGenesAdj, Meta = MetaExpSCZ, title = "AstrocyteGenesAdj, not Adjusted", ExpCol = "NonAdjScaled",
                 geneType = "Astrocyte_Genes", celltype = "Astrocyte_Genes",
                 GeneCol = c("cornflowerblue", "gold3", "coral2", "darkolivegreen4", "grey20"))
dev.off()

#Run ErmineJ
if(!"ermineR" %in% rownames(installed.packages())){
  install_github("oganm/ermineR", force = T)
}
library(ermineR)

rownames(SCZcombine) <- SCZcombine$GeneSymbol
rownames(BPcombine) <- BPcombine$GeneSymbol

Anno_fileGeneric <- read.table("Generic_human_noParents.an.txt", sep = "\t", header = TRUE, comment = "#", quote = "")

CommonGeneAnnoFile <- Anno_fileGeneric %>% filter(GeneSymbols %in% CommonGenes) %>%
  droplevels

write.table(CommonGeneAnnoFile, file ="AnnotationsCommonGenes.txt", row.names = FALSE, col.names = TRUE, sep = "\t")

temp <- gsr(scores = SCZcombine, scoreColumn = "pValDown_mod0", bigIsBetter = F, logTrans = T, annotation = CommonGeneAnnoFile, aspects = "B")
temp2 <- gsr(scores = SCZcombine, scoreColumn = "pValDown_mod1", bigIsBetter = F, logTrans = T, annotation = CommonGeneAnnoFile, aspects = "B")

temp3 <- gsr(scores = BPcombine, scoreColumn = "pValDown_mod0", bigIsBetter = F, logTrans = T, annotation = CommonGeneAnnoFile, aspects = "B")
temp4 <- gsr(scores = BPcombine, scoreColumn = "pValDown_mod1", bigIsBetter = F, logTrans = T, annotation = CommonGeneAnnoFile, aspects = "B")



MitoTerm <- temp2$results %>% data.frame() %>% select(-GeneMembers, -Same.as) %>% .[grep("mitochond|respiratory ele|ATP|oxidative phosph|nucleoside triphosphate",.$Name),] %>% .$Name

MitoRanksSCZ <- data.frame(No_MGPadj = which(temp$results$Name %in% MitoTerm),
                           pAdjNo_MGP = temp$results$CorrectedPvalue[which(temp$results$Name %in% MitoTerm)],
                           MGPadj = which(temp2$results$Name %in% MitoTerm),
                           pAdjMGP = temp2$results$CorrectedPvalue[which(temp2$results$Name %in% MitoTerm)],
                           Group = "SCZ") %>% data.table %>% melt(id.vars = "Group", measure=patterns("MGPadj", "pAdj"), variable.factor = TRUE,
                                                                  variable.name = "Model", value.name = c("Rank", "AdjPval"))
                        
MitoRanksBP <- data.frame(No_MGPadj = which(temp3$results$Name %in% MitoTerm),
                          pAdjNo_MGP = temp3$results$CorrectedPvalue[which(temp3$results$Name %in% MitoTerm)],
                          MGPadj = which(temp4$results$Name %in% MitoTerm),
                          pAdjMGP = temp4$results$CorrectedPvalue[which(temp4$results$Name %in% MitoTerm)],
                          Group = "BP") %>% data.table %>% melt(id.vars = "Group", measure=patterns("MGPadj", "pAdj"), variable.factor = TRUE,
                                                                variable.name = "Model", value.name = c("Rank", "AdjPval"))

MitoRanks <- rbind(MitoRanksSCZ, MitoRanksBP) %>% data.frame()
levels(MitoRanks$Model) <- c("No_MGPadj", "MGPadj")

MitoRanks$Signif <- sapply(MitoRanks$AdjPval, function(x){
  if(x < 0.1){
    "Yes"
  } else {
    "No"
  }
}) %>% factor(levels = c("Yes", "No"))

SCZnoMGPadjSig <- temp$results %>% data.frame()  %>% select(-GeneMembers, -Same.as, -matches("Multifunc|MFP"), -NumProbes) %>% filter(CorrectedPvalue < 0.1)
SCZMGPadjSig <- temp2$results %>% data.frame()  %>% select(-GeneMembers, -Same.as, -matches("Multifunc|MFP"), -NumProbes) %>% filter(CorrectedPvalue < 0.1)
sum(MitoTerm %in% SCZnoMGPadjSig$Name)
sum(MitoTerm %in% SCZMGPadjSig$Name)

ggplot(MitoRanks, aes(Model, Rank)) + 
  theme_classic(base_size = 12) +
  labs(title = "Ranks of mitochondria-related gene-sets (ErmineJ GSR method)") +
  facet_wrap(~Group) +
  geom_violin() +
  scale_y_reverse() +
  geom_boxplot(width = 0.08, outlier.shape = NA) +
  geom_jitter(width = 0.08, alpha = 0.4, size = 1, aes(color = Signif)) +
  scale_color_manual(values = c("darkred", "black"), guide = FALSE)
ggsave("MitoRanks.pdf", path = GeneralResultsPath, width = 8, height = 6, units = "in", dpi=300)

