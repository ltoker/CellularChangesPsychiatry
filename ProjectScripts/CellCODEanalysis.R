source("SetUp.R")
packageF("tidyverse")
packageF("purrr")

install_github("oganm/markerGeneProfile", force = T)
library(markerGeneProfile)
install_github("oganm/homologene", force = T)
library(homologene)
install_github("oganm/CellCODE", force = T)
library(CellCODE)

data('mouseMarkerGenesCombined')

# this function will turn gene lists into CellCODE inputs

cellCodeInput = function(genes){
  
  geneCount = genes %>% map(length)
  
  allGenes = genes %>% unlist %>% unique
  
  
  
  tags = genes %>% lapply(function(x){
    
    out = rep(0,length(allGenes))
    
    out[allGenes %in% x] = 1
    
    return(out)
    
  }) %>% as.data.frame()
  
  
  
  rownames(tags) = allGenes
  
  return(tags)
  
}

# this is marker gene input

geneList = mouseMarkerGenesCombined$Cortex %>% map(mouse2human) %>% map('humanGene')
geneList$GabaPV <- geneList$GabaPV[!geneList$GabaPV %in% c("WIF1", "TMEM132C", "BTN2A2")]

load(paste0(GeneralResultsPath, "DataCombined.rda"))
AnalysisStudies <- names(ExpAll)
  

CellCodeResults <- sapply(AnalysisStudies, function(study){
  data = ExpAll[[study]]
  meta = metaCombined %>% filter(Study == study) %>% droplevels()
  meta %<>% filter(!is.na(Astrocyte_Genes)) %>% droplevels
  
  ExpData <- data %>% select_(.dots = as.character(meta$CommonName)) %>% as.matrix()
  rownames(ExpData) <- data$GeneSymbol
  
  geneList %<>% map(function(x){
    
    x[x %in% rownames(ExpData)]
    
  })
  
  # code works with gene count > 2. Not sure why but it is a decent cutoff anyway
  
  geneList = geneList[sapply(geneList,length)>2]
  
  tagList = cellCodeInput(geneList)
  
  # dataset input is a MATRIX (not data frame) with genes as rownames
  
  groups = meta$Profile
  
  SPVs=getAllSPVs(data = ExpData, 
                  
                  grp = groups %>% as.factor,
                  
                  dataTag = tagList, 
                  
                  method = "mixed", 
                  
                  plot =  T) 
  
  SPVs$SPVs %<>% data.frame()
  SPVs$SPVs <- cbind(meta, SPVs$SPVs)
  SPVs
}, simplify = FALSE)


CellCodeResultsCombined <- lapply(CellCodeResults, function(study){
  study$SPVs %>% select(CommonName, Profile, Profile2, Study, Study2, CollectionType,
                        Age,  pH, Sex,  PMI, RIN, Astrocyte_Genes, Astrocyte, GabaPV_Genes, GabaPV) 
}) %>% do.call(rbind, .)

CellCodeResultsCombined$Profile <- factor(CellCodeResultsCombined$Profile, levels = c("Cont", "BP", "SCZ"))

ggplot(CellCodeResultsCombined, aes(Astrocyte_Genes, Astrocyte)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "MGP", y = "CellCODE") +
  geom_point(size = 0.6, alpha = 0.4) +
  facet_wrap(~Study, scales = "free_x") #+
  #geom_text(data = PC1_PVCor,  mapping = aes(x = x, y = y, label = Cor), inherit.aes = FALSE, color = "red")

ggsave("CellCODE_MGP_GabaPV.pdf", path = GeneralResultsPath, width = 8, height = 6, units = "in", dpi=300)

ggplot(CellCodeResultsCombined, aes(GabaPV_Genes, GabaPV)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "MGP", y = "CellCODE") +
  geom_point(size = 0.6, alpha = 0.4) +
  facet_wrap(~Study, scales = "free_x")
ggsave("CellCODE_MGP_Astro.pdf", path = GeneralResultsPath, width = 8, height = 6, units = "in", dpi=300)

PlotGnes <- function(study){
  temp <- ExpAll[[study]] %>% filter(GeneSymbol %in% CellCodeResults[[study]]$usedGenes$GabaPV)
  tempDF <- temp %>% select(matches("_")) %>% t
  tempDF <- apply(tempDF, 2, scale) %>% data.frame()
  names(tempDF) <- temp$GeneSymbol
  tempDF %<>% mutate(CommonName = grep("_", names(temp), value = T))
  tempDF$Profile <- sapply(tempDF$CommonName, function(x) strsplit(x, "_")[[1]][1])
  temp2 <- metaCombined %>% filter(Study == study, !is.na(Astrocyte_Genes)) %>% .$CommonName %>% as.character
  tempDF %<>% filter(CommonName %in% temp2)
  tempDFmelt <- melt(tempDF, id.vars = c("CommonName", "Profile"), variable.name = "Gene", value.name = "Expression")
  tempDFmelt$Profile <- relevel(as.factor(tempDFmelt$Profile), ref = "Cont")
  
  GrpMedian <- tempDFmelt %>% group_by(Gene,Profile) %>%
    summarise(Median = median(Expression)) %>% data.frame() %>%
    filter(Profile == "Cont") %>% droplevels
  
  ggplot(tempDFmelt, aes(Profile, Expression)) +
    labs(x = "", y = "Normalized expression",title = study) + 
    geom_boxplot(outlier.shape = NA, aes(fill = Profile), alpha = 0.6) +
    scale_fill_manual(values = c("white", "dodgerblue2", "darkolivegreen4"), name = "Group") +
    geom_jitter(size = 0.6) +
    facet_wrap(~Gene, scales = "free_y") +
    geom_hline(data = GrpMedian, aes(yintercept = Median), color = "red", linetype="dashed")
}

PlotGnes("GSE80655")
PlotGnes("GSE53987")
PlotGnes("GSE21138")

#Compare correlation to actual cell count data
Reynolds <- read.table("ReynoldsINCountsSMRI.txt.csv", header = TRUE, sep = "\t")
ConsrtData <- lapply(CellCodeResults[c("Study2AltarC", "Study4Chen")], function(study){
  study$SPVs %>% select(Profile, SubjectID, CommonName, GabaPV, GabaPV_Genes)
})

names(Reynolds) <- sapply(names(Reynolds), function(x){
  x = gsub("\\.", "", x)
  gsub("PV", "BA", x)
})

names(Reynolds)[grepl("[0-9]", names(Reynolds))] <- sapply(grep("[0-9]", names(Reynolds), value = T), function(x){
  paste0("Count_", x)
})

Reynolds <- merge(Reynolds, ConsrtData$Study2AltarC, by.x = "ID", by.y = "SubjectID", all = TRUE)
Reynolds <- merge(Reynolds, ConsrtData$Study4Chen, by.x = "ID", by.y = "SubjectID",
                  suffixes = c("_AltarC","_Chen"), all = TRUE)

names(Reynolds) <- sapply(names(Reynolds), function(x){
  x <- gsub("GabaPV_Genes", "MGP", x)
  gsub("GabaPV", "CellCODE", x)
})

CorExpCount <- cor(Reynolds %>% select(Count_BA46, Count_BA9, CellCODE_Chen, CellCODE_AltarC, MGP_Chen, MGP_AltarC),
                   method = "spearman", use = "pairwise.complete.obs")
CorExpCount <- apply(CorExpCount, c(1,2), function(x) round(x, digits = 2))

diag(CorExpCount) <- NA

CorExpCount2 <- CorExpCount %>%
  tbl_df() %>%
  rownames_to_column('Var1') %>%
  gather(Var2, value, -Var1) %>%
  mutate(Var1 = factor(rownames(CorExpCount)[as.numeric(Var1)], levels = rownames(CorExpCount)),
         Var2 = factor(Var2, levels = rev(rownames(CorExpCount)))) 

ggplot(CorExpCount2, aes(Var1, Var2)) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1 )) +
  labs(x = "", y = "") +
  geom_tile(aes(fill = value)) + 
  geom_text(aes(label = value)) +
  scale_fill_gradient(low = "white", high = "red", name = "Spearman\ncorrelation", na.value = "grey50")          
ggsave("CorExpCount.pdf", path = GeneralResultsPath, width = 8, height = 6, units = "in", dpi=300)
