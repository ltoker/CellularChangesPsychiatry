source("SetUp.R")

load(paste0(GeneralResultsPath, "studyFinalStanleyConsortium.rda"))

plotCountsMGPgroupResults <- function(data){
  ContMedian <- data %>% filter(Profile == "Cont") %>% group_by(MeasureType) %>%
    summarise(median = median(Value, na.rm = T)) %>% data.frame
  
  #Plot group level results for counts vs MGP
  ggplot(data, aes(Profile, Value)) +
    theme_bw(base_size = 16) +
    labs(x = "", y ="Cell count / MGP") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(), 
          axis.text.x = element_text(size = rel(0.9)),
          axis.text.y = element_text(size = rel(0.9))) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size=1)+
    facet_wrap(~MeasureType, scales = "free") +
    geom_hline(data = ContMedian, aes(yintercept = median), color = "red", linetype="dashed")
}

plotCountsMGPpairs <- function(data, labels, pch = 16, col = "grey", method = "spearman"){
  panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, method = "spearman", ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = "pairwise.complete.obs", method = method)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor)
  }
  pairs(data, labels,
        lower.panel = function(...) {
          panel.cor(method = method, cex.cor=1.5, prefix = "r = ", ...)
        },
        xaxt = "n", yaxt = "n", pch = 16, col = "grey")
}


MetadataAltarC <- studyFinal$Study2AltarC$Metadata
MetadataChen4 <- studyFinal$Study4Chen$Metadata


plotCountsMGPgroupResults <- function(data, pvalTxtSize = 6){
  ContMedian <- data %>% filter(Profile == "Cont") %>% group_by(MeasureType) %>%
    summarise(median = median(Value, na.rm = T)) %>% data.frame
  
  #Plot group level results for counts vs MGP
  ggplot(data, aes(Profile, Value)) +
    theme_bw(base_size = 16) +
    labs(x = "", y = "Cell count (cell/mm^2)   /   MGP") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(), 
          axis.text.x = element_text(size = rel(0.9)),
          axis.text.y = element_text(size = rel(0.9))) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size=1.5)+
    facet_wrap(~MeasureType, scales = "free", nrow = 2) +
    geom_text(data = pVal,  aes(x, y, label = pValue), size=pvalTxtSize, inherit.aes = FALSE) +
    geom_hline(data = ContMedian, aes(yintercept = median), color = "red", linetype="dashed")
}

plotCountsMGPpairs <- function(data, labels, pch = 16, col = "grey", method = "spearman"){
  panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, method = "spearman", ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor.test(x, y, use = "pairwise.complete.obs", method = method)
    
    txt <- format(c(r$estimate, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    txt2 <- paste0("p = ", scientific(r$p.value, digits = 2))
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor)
    text(0.5, 0.3, txt2, cex = cex.cor)
  }
  pairs(data, labels,
        lower.panel = function(...) {
          panel.cor(method = method, cex.cor=1, prefix = "r = ", ...)
        },
        xaxt = "n", yaxt = "n", pch = 16, col = "grey")
}



Reynolds <- read.table("ReynoldsINCountsSMRI.txt.csv", header = TRUE, sep = "\t")

Reynolds$GabaPV_MGPAltarC <- MetadataAltarC$GabaPV_Genes[match(Reynolds$ID, MetadataAltarC$StanleyID)]
Reynolds$GabaPV_MGPChen4 <- MetadataChen4$GabaPV_Genes[match(Reynolds$ID, MetadataChen4$StanleyID)]
Reynolds$CommonName <- MetadataAltarC$CommonName[match(Reynolds$ID, MetadataAltarC$StanleyID)]
Reynolds$Profile <- MetadataAltarC$Profile[match(Reynolds$ID, MetadataAltarC$StanleyID)]
Reynolds <- merge(Reynolds, MetadataAltarC %>% select(StanleyID, Age, Sex, pH, PMI), by.x = "ID", by.y = "StanleyID", all = TRUE)

Reynolds %<>% arrange(GabaPV_MGPChen4) 
Reynolds$ID <- factor(Reynolds$ID, levels = Reynolds$ID)
Reynolds$GabaPV_meanCount <- apply(Reynolds %>% select(PV.46, PV.9), 1, function(x) mean(x))

Reynolds %<>% mutate(PV.46norm = rescale(PV.46, c(0,1)),
                      PV.9norm = rescale(PV.9, c(0,1)))


ReynoldsMelt <- melt(Reynolds %>% select(matches("ID|CommonName|PV\\.46$|PV\\.9$|GabaPV_MGP|Profile|Sex|pH|PMI|Age")),
                     id.vars = c("ID", "CommonName",
                                 "Profile", "Sex", "Age", "pH", "PMI"),
                     variable.name = "MeasureType", value.name = "Value") %>% .[complete.cases(.),] %>% droplevels()

ReynoldsMelt$MeasureType <- factor(ReynoldsMelt$MeasureType, levels = c("GabaPV_MGPAltarC","GabaPV_MGPChen4","PV.9","PV.46"))
levels(ReynoldsMelt$MeasureType) <- c("GabaPV MGP (BA46/BA10)\nStudy2AltarC",
                                      "GabaPV MGP (BA6)\nStudy4Chen",
                                      "GabaPV cell count (BA9) \nBeasley 2002",
                                      "GabaPV cell count (BA46) \nBeasley 2002")
ReynoldsMelt$MeasureType2 <- sapply(ReynoldsMelt$MeasureType, function(x){
  if(grepl("MGP", x)){
    "MGP"
  } else {
    "cell/mm^2"
  }
})

#Plot group level results for counts vs MGP
pVal <- data.frame(Study = levels(ReynoldsMelt$MeasureType))
Vals <- sapply(pVal$Study, function(study){
  study = ReynoldsMelt %>% filter(MeasureType == study) %>% droplevels()
  BP <- wilcox.test(Value~Profile, data = study %>% filter(Profile != "SCZ"))
  if("SCZ" %in% study$Profile){
    SCZ <-  wilcox.test(Value~Profile, data = study %>% filter(Profile != "BP"))
  } else {
    SCZ <- list(p.value=NA)
  }
  c(scientific(BP$p.value), scientific(SCZ$p.value))
})

pVal <- cbind(pVal, t(Vals))
names(pVal) <- c("MeasureType", "BP", "SCZ")

pVal <- melt(pVal, id.vars = "MeasureType", variable.name = "Profile", value.name = "pValue") %>% mutate(pValue = as.numeric(pValue))
pVal <- pVal[!is.na(pVal$pValue),]
pVal %<>% mutate(pValue = paste0("p = ", scientific(pValue, digits = 2)))

pVal$x <- c(rep(2, 4), rep(3, 3))
pVal$y <- c(1,1, 40, 43, 1, 40, 43)

p <- plotCountsMGPgroupResults(ReynoldsMelt, pvalTxtSize = 2.5)
ggsave(filename = paste0(GeneralResultsPath, "/GabaPV_SMRIcountsMGP.pdf"), plot = p, width = 6, height = 6, units = "in",
              dpi=300)


#Plot correlations of cell counts vs MGPs
pdf(paste0(GeneralResultsPath, "/GabaPV_SMRIcountsMGPcor.pdf"), width = 6, height = 6, pointsize = 14, useDingbats = FALSE)
plotCountsMGPpairs(data  = Reynolds %>% select(PV.46, PV.9, GabaPV_MGPChen4, GabaPV_MGPAltarC),
                   labels = c("Cell count\nBA46\n\nBeasley 2002 " ,
                              "Cell count\nBA9\n\nBeasley 2002",
                              "MGP \nBA6\n\nStudy4Chen",
                              "MGP \nBA46/BA10\n\nStudy2AltarC"))
dev.off()

ReynoldsMelt2 <- melt(Reynolds %>% select(matches("ID|CommonName|PV\\.46$|PV\\.9$|GabaPV_MGP|Profile|Sex|pH|PMI|Age")),
                      id.vars = c("ID", "CommonName",
                                  "Profile", "PV.46", "PV.9", "Sex", "Age", "pH", "PMI"),
                      variable.name = "Study", value.name = "MGP") %>% .[complete.cases(.),] %>% droplevels()
ReynoldsMelt2$Study <- sapply(ReynoldsMelt2$Study, function(x){
  gsub("GabaPV_MGP", "", x)
}) %>% factor

lmm0 <- lmer(MGP~Age + Sex + pH + PMI + (1|Study) + (1|ID), REML = FALSE, data = ReynoldsMelt2)
lmm1 <- lmer(MGP~Age + Sex + pH + PMI + PV.9 + (1|Study) + (1|ID), REML = FALSE, data = ReynoldsMelt2)
lmm1b <- lmer(MGP~Age + Sex + pH + PMI + PV.46 + (1|Study) + (1|ID), REML = FALSE, data = ReynoldsMelt2)
lmm1c <- lmer(MGP~Profile + Age + Sex + pH + PMI + (1|Study) + (1|ID), REML = FALSE, data = ReynoldsMelt2)
lmm2 <- lmer(MGP~Age + Sex + pH + PMI + PV.46 + PV.9 + (1|Study) + (1|ID), REML = FALSE, data = ReynoldsMelt2)
lmm3 <- lmer(MGP~Profile + Age + Sex + pH + PMI + PV.9 + (1|Study) + (1|ID), REML = FALSE, data = ReynoldsMelt2)
lmm3b <- lmer(MGP~Profile + Age + Sex + pH + PMI + PV.9 + PV.46 + (1|Study) + (1|ID), REML = FALSE, data = ReynoldsMelt2)
lmm3c <- lmer(MGP~Profile + Age + Sex + pH + PMI + (1|Study) + (1|ID), REML = FALSE, data = ReynoldsMelt2)


ANOVAres <- anova(lmm0, lmm3c, lmm3b)
Models <- attr(ANOVAres, "heading")[3:5] %>% data.frame()
names(Models) <- "Model"

write.table(Models, paste0(GeneralResultsPath, "MGPcountModel.tsv"), append = FALSE, row.names = FALSE, sep = "\t")
write.table(ANOVAres, paste0(GeneralResultsPath, "MGPcountModel.tsv"), append = TRUE, row.names = TRUE, sep = "\t")

NoAdjMod <- summary(lmm3c)$coefficients
NoAdjMod <- rbind(NoAdjMod, NA)
rownames(NoAdjMod)[nrow(NoAdjMod)] <- "PV.9"
AdjMod <- summary(lmm3)$coefficients

EffectTable <- cbind(NoAdjMod, AdjMod)
EffectTable <- apply(EffectTable, c(1,2), function(x) round(x, digits  = 2))
write.table(EffectTable, paste0(GeneralResultsPath, "MGPcountEffectTable.tsv"), row.names = TRUE, sep = "\t")
