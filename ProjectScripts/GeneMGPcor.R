source("SetUp.R")

if(!"markerGeneProfile" %in% rownames(installed.packages())){
  install_github("oganm/markerGeneProfile", force = T)
}
library(markerGeneProfile)
data("mouseMarkerGenesCombined")

plotTopCor <- function(cellMGP, markerInfo = NULL, corInfo = NULL, geneNum = 10,
                       txtSize = 16, ptSize = 0.4, title = "top correlated genes",
                       colors = c("darkorange", "firebrick1", "darkorchid1",
                                  "darkolivegreen4", "darkseagreen", "dodgerblue3")) {
  plot <- GetHumanExp(genes = allStudyCorMean %>% filter(CellType == cellMGP) %>% .$GeneSymbol %>% head(10),
                      markerInfo =  allStudyCorMean %>% filter(CellType == cellMGP) %>% select(GeneSymbol, GeneType) %>% head(geneNum),
                      corInfo = allStudyCorMean %>% filter(CellType == cellMGP) %>% select(GeneSymbol, meanCor) %>% head(geneNum),
                      CellType = strsplit(cellMGP, "_")[[1]][1],
                      txtSize = txtSize, ptSize = ptSize,
                      title = paste0(gsub("_", " ", cellMGP), " - ", title), color = colors)
  return(plot)
}

plotMistryMGP <- function(data, MGP, geneType, geneTypeVal, study = "All", groups = "All", title, txtSize = 14, ptSize = 1.5, alpha = 1,
                          bxptSize = 0.3){
  celltype = strsplit(MGP, "_")[[1]][1]
  if(groups %in% c("All")) {
    groups = unique(data$Profile)
  }
  if(study %in% c("All")){
    study = unique(data$Study)
  }
  data = data  %>% filter(Study %in% study,Profile %in% groups, CellType == MGP) %>%
    .[.[[geneType]] %in% geneTypeVal,] %>% droplevels
  data$Profile <- factor(as.character(data$Profile), levels = c("Cont", "BP", "SCZ"))
  p <- ggplot(data, aes(Profile, y = Cor))
  if(study %in% c("Combined")){
    Plot <- p +
      labs(x = "", y = paste0("Gene-", "MGP cor"), title = title)+
      scale_y_continuous(limits = c(-1, 1))+
      theme_classic(base_size = txtSize) +
      theme(axis.text.y = element_text(size = rel(1.1)),
            axis.text.x = element_text(size = rel(1.2)),
            legend.position = "none",
            panel.grid = element_blank()) +
      geom_boxplot(width = bxptSize, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = ptSize, alpha = alpha) +
      geom_abline(intercept = 0, slope = 0, color = "red", linetype="dashed")
  } else {
    Plot <- p +
      labs(x = "", y = paste0("Gene-", " MGP correlation"), title = title)+
      scale_y_continuous(limits = c(-1, 1))+
      theme_classic(base_size = txtSize) +
      theme(axis.text.y = element_text(size = rel(1.2)),
            axis.text.x = element_text(size = rel(1.2)),
            legend.position = "none",
            panel.grid = element_blank()) +
      geom_violin() +
      geom_boxplot(width = bxptSize, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = ptSize, alpha = alpha) +
      geom_abline(intercept = 0, slope = 0, color = "red", linetype="dashed") +
      facet_wrap(~Study, nrow = 3)
  }
  return(Plot)
}

PlotSingleGeneMGPCor <- function(ExpData, MGPdata, gene,
                                 method = "spearman", celltype, celltypeExtent = "_Genes",
                                 NamepColMGPData = "CommonName",
                                 base_size = 14){
  cellColumn = paste0(celltype, celltypeExtent) 
  temp <- ExpData %>%
    filter(GeneSymbol %in% gene)
  temp <- temp[sapply(temp[1,], function(x) is.numeric(x))] %>% unlist
  MGPdata[[gene]] <- temp[match(MGPdata[[NamepColMGPData]], names(temp))]
  MGPdataSub <- MGPdata %>% select_(.dots = c(gene, cellColumn, "Profile"))
  names(MGPdataSub)[1:2] <- c("Gene", "MGP")
  Cor <-sapply(levels(MGPdataSub$Profile), function(grp){
    cor.test(~Gene+MGP,
             data = MGPdataSub %>% filter(Profile == grp), method = method) %>%
      .$estimate %>% signif(digits = 2) %>% paste0("(rho = ", ., ")")
  }, simplify = FALSE) %>% do.call(rbind, .) %>% data.frame()
  names(Cor) <- "rho"
  levels(MGPdataSub$Profile) <- paste(levels(MGPdataSub$Profile), Cor$rho)
  ggplot(MGPdataSub, aes(Gene,MGP)) +
    labs(y = paste0(celltype," MGP"), x="Expression (log2)", title = gene) +
    theme_classic(base_size = base_size) +
    theme(axis.text.y = element_text(size = rel(1.1)),
          axis.text.x = element_text(size = rel(1.1)),
          panel.grid = element_blank()) +
    facet_wrap(~Profile) +
    geom_point()
  
}

MistryGenesUP <- read.table("MistryUp.txt", header = T, sep="\t") 
MistryGenesDown <- read.table("MistryDown.txt", header = T, sep="\t") 

#Update Mistry et al gene names
Anno_file <- read.table("GPL570.gz", comment="#", header=T, quote='"', sep="\t")
MistryGenesDown$GeneSymbolOrg <- MistryGenesDown$GeneSymbol
MistryGenesDown$GeneSymbol <- sapply(MistryGenesDown$Probe, function(x) Anno_file %>% filter(ProbeName == as.character(x)) %>%
                                       .$GeneSymbol) %>% make.names

MistryGenesUP$GeneSymbolOrg <- MistryGenesUP$GeneSymbol
MistryGenesUP$GeneSymbol <- sapply(MistryGenesUP$Probe, function(x) Anno_file %>% filter(ProbeName == as.character(x)) %>%
                                       .$GeneSymbol) %>% make.names

if(length(list.files(path = "GeneralResults", pattern = "GeneMGPcor.rda")) == 1){
  load("GeneralResults/GeneMGPcor.rda")
} else {
  load("AnalysisStudies.rda")
  StanleyStudies <- read.table("StanleyStudies.txt", sep="\t", header=T)
  for(name in AnalysisStudies){
    load(paste0(GeneralResultsPath, "MGPcorAllgenes_", name,".rda"))
    rm(list = ls(pat=paste0(name, "_Cortex$")))
    MGPcorAllgenes <- lapply(MGPcorAllgenes, function(std){
      lapply(std, function(cellType){
        lapply(cellType, function(grp){
          grp$GeneSymbol <- rownames(grp)
          grp %>% arrange(desc(Cor))
        })
      })
    })
    
    if(grepl("Stanley", name)){
      for(study in names(MGPcorAllgenes)){
        investigator = gsub("Study.[0-9]?", "", study)
        region = StanleyStudies %>%
          filter(Investigator == investigator) %>%
          select(NeuExpRegion) %>%
          unlist %>% as.character
        assign(paste0(study, "_", region), MGPcorAllgenes[[study]])
      }
    } else if(name == "GSE80655"){
      assign(paste0(name, "_Cortex"), MGPcorAllgenes$DLPFC)
    } else {
      assign(paste0(name, "_Cortex"), MGPcorAllgenes$Cortex)
    }
  }
  
  allStudyCor <- sapply(ls(pat = "_Cortex$"), function(study){
    data = eval(as.name(study))
    data2 <- sapply(names(data), function(celltype){
      cellData <- data[[celltype]]
      DF <- do.call(rbind, cellData)
      DF$CellType <- celltype
      DF$Profile <- sapply(rownames(DF), function(x) strsplit(x, "\\.")[[1]][1])
      DF$Study <- study
      DF
    }, simplify = FALSE)
    do.call(rbind, data2)
  }, simplify = FALSE) %>% do.call(rbind, .) %>% data.frame()
  allStudyCor %<>% mutate_each(funs(as.factor), -Cor)
}

allStudyGroupMean <- allStudyCor %>% group_by(CellType, Profile, GeneSymbol) %>%
  summarise(Cor = mean(Cor)) %>% arrange(CellType, Profile, desc(Cor)) %>% data.frame()

# Add information regarding Mistry et al gene signature
allStudyCor$MistryGenes <- sapply(allStudyCor$GeneSymbol, function(gene){
  if(gene %in% MistryGenesUP$GeneSymbol){
    "MistryUP"
  } else if (gene %in% MistryGenesDown$GeneSymbol) {
    "MistryDown"
  } else {
    NA
  }
})
allStudyCorMean$MistryGenes <- sapply(allStudyCorMean$GeneSymbol, function(gene){
  if(gene %in% MistryGenesUP$GeneSymbol){
    "MistryUP"
  } else if (gene %in% MistryGenesDown$GeneSymbol) {
    "MistryDown"
  } else {
    NA
  }
})
allStudyGroupMean$MistryGenes <- sapply(allStudyGroupMean$GeneSymbol, function(gene){
  if(gene %in% MistryGenesUP$GeneSymbol){
    "MistryUP"
  } else if (gene %in% MistryGenesDown$GeneSymbol) {
    "MistryDown"
  } else {
    NA
  }
})

allStudyCorMean$Study <- "Combined"
allStudyGroupMean$Study <- "Combined"

AstroMarker <- mouseMarkerGenesCombined$Cortex$Astrocyte
AstroMarkerHuman <- sapply(AstroMarker,
                           function(gene) mouse2human(gene)$humanGene) %>% unlist
GabaPVMarker <- mouseMarkerGenesCombined$Cortex$GabaPV
GabaPVMarkerHuman <- sapply(GabaPVMarker,
                           function(gene) mouse2human(gene)$humanGene) %>% unlist

#Remove PV markers which are not specific to neuros based on Darmanis data
GabaPVMarkerHuman <- GabaPVMarkerHuman[!GabaPVMarkerHuman %in% c("WIF1", "TMEM132C", "BTN2A2")] 

# Add information regarding whether the gene is a marker
allStudyCor$GeneType <- sapply(allStudyCor$GeneSymbol, function(gene){
  if(gene %in% AstroMarkerHuman){
    "Astrocyte"
  } else if (gene %in% GabaPVMarkerHuman) {
    "GabaPV"
  } else {
    NA
  }
})
allStudyCorMean$GeneType <- sapply(allStudyCorMean$GeneSymbol, function(gene){
  if(gene %in% AstroMarkerHuman){
    "Astrocyte"
  } else if (gene %in% GabaPVMarkerHuman) {
    "GabaPV"
  } else {
    NA
  }
})
allStudyCorMean$GeneType <- sapply(allStudyCorMean$GeneSymbol, function(gene){
  if(gene %in% AstroMarkerHuman){
    "Astrocyte"
  } else if (gene %in% GabaPVMarkerHuman) {
    "GabaPV"
  } else {
    NA
  }
})
allStudyGroupMean$GeneType <- sapply(allStudyGroupMean$GeneSymbol, function(gene){
  if(gene %in% AstroMarkerHuman){
    "Astrocyte"
  } else if (gene %in% GabaPVMarkerHuman) {
    "GabaPV"
  } else {
    NA
  }
})


GranthGWAS <- read.table("GranthSCZ.csv", header = TRUE, sep = "\t") %>% select(matches("^Sch"))
names(GranthGWAS) <- sapply(names(GranthGWAS), function(x) gsub("(\\.)*", "", x))

# Add information regarding whether the gene was implicated in SCZ GWAS study
allStudyCor$GWAS <- sapply(allStudyCor$GeneSymbol, function(gene){
  if(gene %in% GranthGWAS$SchizophreniaGWAS){
    "YES"
  } else {
    "NO"
  }
})
allStudyCorMean$GWAS <- sapply(allStudyCorMean$GeneSymbol, function(gene){
  if(gene %in% GranthGWAS$SchizophreniaGWAS){
    "YES"
  } else {
    "NO"
  }
})
allStudyCorMean$GWAS <- sapply(allStudyCorMean$GeneSymbol, function(gene){
  if(gene %in% GranthGWAS$SchizophreniaGWAS){
    "YES"
  } else {
    "NO"
  }
})
allStudyGroupMean$GWAS <- sapply(allStudyGroupMean$GeneSymbol, function(gene){
  if(gene %in% GranthGWAS$SchizophreniaGWAS){
    "YES"
  } else {
    "NO"
  }
})

ChngUP <- read.table("ASDup.csv", header = TRUE, sep = "\t")$Gene.Symbol
ChngDown <- read.table("ASDdown.csv", header = TRUE, sep = "\t")$Gene.Symbol

# Add information regarding whether the gene was changed in Ch'ng ASD meta-analysis
allStudyCor$ASDgenes <- sapply(allStudyCor$GeneSymbol, function(gene){
  if(gene %in% ChngUP){
    "UP"
  } else if (gene %in% ChngDown) {
    "Down"
  } else {
    NA
  }
})
allStudyCorMean$ASDgenes <- sapply(allStudyCorMean$GeneSymbol, function(gene){
  if(gene %in% ChngUP){
    "UP"
  } else if (gene %in% ChngDown) {
    "Down"
  } else {
    NA
  }
})
allStudyCorMean$ASDgenes <- sapply(allStudyCorMean$GeneSymbol, function(gene){
  if(gene %in% ChngUP){
    "UP"
  } else if (gene %in% ChngDown) {
    "Down"
  } else {
    NA
  }
})
allStudyGroupMean$ASDgenes <- sapply(allStudyGroupMean$GeneSymbol, function(gene){
  if(gene %in% ChngUP){
    "UP"
  } else if (gene %in% ChngDown) {
    "Down"
  } else {
    NA
  }
})


#Plot top correlated genes in human single cell data
source(paste0(ProjScriptPath, "HumanSingleCell.R"))

TopCorAstroPlot <- plotTopCor("Astrocyte_MGP", txtSize = 15)
ggsave("GeneralResults/TopCorAstroPlot.pdf", TopCorAstroPlot, width = 10, height = 6, units = "in",device = "pdf",
       dpi=300, useDingbats = FALSE)

TopCorGabaPVPlot <- plotTopCor("GabaPV_MGP", txtSize = 15)
ggsave("GeneralResults/TopCorGabaPVPlot.pdf", TopCorGabaPVPlot, width = 10, height = 6, units = "in",device = "pdf",
       dpi=300, useDingbats = FALSE)

#Plot individual studies
GabaPVdwnAllPlot <- plotMistryMGP(data = allStudyCor, MGP = "GabaPV_MGP", geneType = "MistryGenes", geneTypeVal = "MistryDown",
                                  title = "Mistry et al.- downregulated genes", ptSize = 0.6, alpha = 0.6)
  
ggsave("GeneralResults/GabaPVdwnAllPlot.pdf", GabaPVdwnAllPlot, width = 10, height = 6, units = "in",device = "pdf",
       dpi=300, useDingbats = FALSE)

AstroupAllPlot <- plotMistryMGP(data = allStudyCor, MGP = "Astrocyte_MGP", geneType = "MistryGenes", geneTypeVal = "MistryUP",
                                  title = "Mistry et al.- upregulated genes", ptSize = 1.2, alpha = 0.6)
ggsave("GeneralResults/AstroupAllPlot.pdf", AstroupAllPlot, width = 10, height = 6, units = "in",device = "pdf",
       dpi=300, useDingbats = FALSE)


GabaPVMGPmarker <- plotMistryMGP(data = allStudyCor, MGP = "GabaPV_MGP", geneType = "GeneType", geneTypeVal = "GabaPV",
                                 title = "GabaPV markers", ptSize = 1.2, alpha = 0.6)
ggsave("GeneralResults/GabaPVMGPmarker.pdf", GabaPVMGPmarker ,width = 10, height = 6, units = "in",device = "pdf",
       dpi=300, useDingbats = FALSE)

AstrocyteMGPmarker <- plotMistryMGP(data = allStudyCor, MGP = "Astrocyte_MGP", geneType = "GeneType", geneTypeVal = "Astrocyte",
                                 title = "Astrocyte markers", ptSize = 0.6, alpha = 0.6)
ggsave("GeneralResults/AstrocyteMGPmarker.pdf", AstrocyteMGPmarker ,width = 10, height = 6, units = "in",device = "pdf",
       dpi=300, useDingbats = FALSE)


#Plot mean correlation
GabaPVdwnMeanPlot <- plotMistryMGP(data = allStudyGroupMean, MGP = "GabaPV_MGP", study = "Combined", geneType = "MistryGenes", geneTypeVal = "MistryDown",
                                  title = "SCZ downregulated genes - GabaPV MGP correlation", alpha = 0.5, bxptSize = 0.6)
AstroupMeanPlot <- plotMistryMGP(data = allStudyGroupMean, MGP = "Astrocyte_MGP", study = "Combined", geneType = "MistryGenes", geneTypeVal = "MistryUP",
                                title = "SCZ upregulated genes - Astrocyte MGP correlation", alpha = 0.5, bxptSize = 0.6)
GabaPVMGPmarkerMean <- plotMistryMGP(data = allStudyGroupMean, MGP = "GabaPV_MGP", study = "Combined", geneType = "GeneType", geneTypeVal = "GabaPV",
                                 title = "GabaPV markers - GabaPV MGP correlation", alpha = 0.5, bxptSize = 0.6)
AstrocyteMGPmarkerMean <- plotMistryMGP(data = allStudyGroupMean, MGP = "Astrocyte_MGP", study = "Combined",geneType = "GeneType", geneTypeVal = "Astrocyte",
                                    title = "Astrocyte markers - Astrocyte MGP correlation", alpha = 0.5, bxptSize = 0.6)

#Plot single gene MGP correlation - this is just for a specific example
load("GeneralResults/studyFinalGSE35978.rda")
PvalbMGPcor <- PlotSingleGeneMGPCor(ExpData = studyFinal$Cortex$aned_high,
                     MGPdata = studyFinal$Cortex$Metadata,
                     gene = "PVALB",celltype = "GabaPV")
Fbxo9MGPor <- PlotSingleGeneMGPCor(ExpData = studyFinal$Cortex$aned_high,
                     MGPdata = studyFinal$Cortex$Metadata,
                     gene = "FBXO9",celltype = "GabaPV")

ggdraw() +
  draw_plot(PvalbMGPcor, 0, .66, .5, .33) +
  draw_plot(Fbxo9MGPor, .5, .66, .5, .33) +
  draw_plot(GabaPVMGPmarkerMean,0, .33, .5, 0.33) +
  draw_plot(GabaPVdwnMeanPlot, 0.5, .33, .5, .33) +
  draw_plot(AstrocyteMGPmarkerMean, 0, 0, .5, .33) +
  draw_plot(AstroupMeanPlot, 0.5, 0, .5, .33) +
  draw_plot_label(c("A", "B", "C", "D", "E", "F"),
                  c(0, 0.5, 0, 0.5, 0, 0.5),
                  c(1, 1,0.67, 0.67, 0.34, 0.34), size = 18)
ggsave("GeneralResults/FullMistryPlot.pdf" ,width = 12, height = 8, units = "in",device = "pdf",
       dpi=300, useDingbats = FALSE)


save.image("GeneralResults/GeneMGPcorWS.Rdata")
