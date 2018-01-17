source("SetUp.R")

load(paste0(GeneralResultsPath, "studyFinalStanleyConsortium.rda"))

plotCountsMGPgroupResults <- function(data){
  ContMedian <- data %>% filter(Profile == "Cont") %>% group_by(MeasureType) %>%
    summarise(median = median(Value, na.rm = T)) %>% data.frame
  
  #Plot group level results for counts vs MGP
  ggplot(data, aes(Profile, Value)) +
    theme_bw(base_size = 16) +
    labs(x = "", y = "Normalized cell count / MGP") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(), 
          axis.text.x = element_text(size = rel(0.9)),
          axis.text.y = element_text(size = rel(0.9))) +
    ylim(0, 1.1) +
    #geom_violin(width = 0.8)+
    geom_boxplot(width = 0.6, outlier.size = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size=1)+
    facet_wrap(~MeasureType, scales = "free_x", nrow = 2) +
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

Uranova <- read.table("URANOVAraw.txt", header = FALSE, sep = "\t")
Uranova %<>% select(V1, V2, V3, V4, V5, V6, V7, V8) %>% filter(V1 %in% c("Code", "Area 9", "Layer VI"))
Uranova$V1 <- as.character(Uranova$V1)
Uranova$V1[seq(4, 40, 4)] <- "WM"
Uranova$V1[seq(3, 39, 4)] <- "GMse"
Uranova$V1[seq(2, 38, 4)] <- "GM"

Uranova <- data.frame(StanleyID = Uranova %>% filter(V1 == "Code") %>% select(-V1, -V2) %>% as.list() %>% unlist %>% as.character(),
                      GMCount = Uranova %>% filter(V1 == "GM") %>% select(-V1, -V2) %>% as.list() %>% unlist %>% as.character %>% as.numeric(),
                      GMse = Uranova %>% filter(V1 == "GMse") %>% select(-V1, -V2) %>% as.list() %>% unlist %>% as.character %>% as.numeric(),
                      WMCount = Uranova %>% filter(V1 == "WM") %>% select(-V1, -V2) %>% as.list() %>% unlist %>% as.character %>% as.numeric())
Uranova$Oligo_MGPAltarC <- MetadataAltarC$Oligo_Genes[match(Uranova$StanleyID, MetadataAltarC$StanleyID)]
Uranova$Oligo_MGPChen4 <- MetadataChen4$Oligo_Genes[match(Uranova$StanleyID, MetadataChen4$StanleyID)]
Uranova$CommonName <- MetadataAltarC$CommonName[match(Uranova$StanleyID, MetadataAltarC$StanleyID)]
Uranova$Profile <- MetadataAltarC$Profile[match(Uranova$StanleyID, MetadataAltarC$StanleyID)]

Uranova %<>% arrange(GMCount)
Uranova %<>% mutate(GMCountnorm = rescale(GMCount, c(0,1)))
#Uranova %<>% filter(!is.na(.$Oligo_MGPAltarC))
UranovaMelt <- melt(Uranova %>% select(matches("ID|CommonName|norm|MGP|Profile")), id.vars = c("StanleyID", "CommonName", "Profile"),
                    variable.name = "MeasureType", value.name = "Value") %>% .[complete.cases(.),] %>% droplevels()

levels(UranovaMelt$MeasureType) <- c("Oligo MGP (BA46/BA10)\nStudy2AltarC", "Oligo MGP (BA6)\nStudy4Chen", "Oligo cell count (BA9, L4)\nUranova 2004")

plotCountsMGPgroupResults <- function(data){
  ContMedian <- data %>% filter(Profile == "Cont") %>% group_by(MeasureType) %>%
    summarise(median = median(Value, na.rm = T)) %>% data.frame
  
  #Plot group level results for counts vs MGP
  ggplot(data, aes(Profile, Value)) +
    theme_bw(base_size = 16) +
    labs(x = "", y = "Normalized cell count / MGP") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(), 
          axis.text.x = element_text(size = rel(0.9)),
          axis.text.y = element_text(size = rel(0.9))) +
    ylim(0, 1.1) +
    #geom_violin(width = 0.8)+
    geom_boxplot(width = 0.6, outlier.size = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size=1.5)+
    facet_wrap(~MeasureType, scales = "free_x", nrow = 2) +
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

plotCountsMGPgroupResults(UranovaMelt)
ggsave("GeneralResults/Oligo_SMRIcountsMGP.pdf", width = 6, height = 6, units = "in",
       dpi=300, useDingbats = FALSE)


pdf("GeneralResults/Oligo_SMRIcountsMGPcor.pdf", width = 6, height = 6,  pointsize = 14, useDingbats = FALSE)
plotCountsMGPpairs(data  = Uranova %>% select(GMCount, Oligo_MGPAltarC),
                   labels = c("Cell count\nBA9, L4\n\nUranova 2004 " ,
                              #"MGP \nBA6\n\nStudy4Chen",
                              "MGP \nBA46/BA10\n\nStudy2AltarC"))
dev.off()


Reynolds <- read.table("ReynoldsINCountsSMRI.txt.csv", header = TRUE, sep = "\t")

Reynolds$GabaPV_MGPAltarC <- MetadataAltarC$GabaPV_Genes[match(Reynolds$ID, MetadataAltarC$StanleyID)]
Reynolds$GabaPV_MGPChen4 <- MetadataChen4$GabaPV_Genes[match(Reynolds$ID, MetadataChen4$StanleyID)]
Reynolds$CommonName <- MetadataAltarC$CommonName[match(Reynolds$ID, MetadataAltarC$StanleyID)]
Reynolds$Profile <- MetadataAltarC$Profile[match(Reynolds$ID, MetadataAltarC$StanleyID)]

Reynolds %<>% arrange(GabaPV_MGPChen4) 
Reynolds$ID <- factor(Reynolds$ID, levels = Reynolds$ID)
Reynolds$GabaPV_meanCount <- apply(Reynolds %>% select(PV.46, PV.9), 1, function(x) mean(x))

Reynolds %<>% mutate(PV.46norm = rescale(PV.46, c(0,1)),
                      PV.9norm = rescale(PV.9, c(0,1)))


ReynoldsMelt <- melt(Reynolds %>% select(matches("ID|CommonName|norm|GabaPV_MGP|Profile")), id.vars = c("ID", "CommonName", "Profile")
                     , variable.name = "MeasureType", value.name = "Value") %>% .[complete.cases(.),] %>% droplevels()

ReynoldsMelt$MeasureType <- factor(ReynoldsMelt$MeasureType, levels = c("GabaPV_MGPAltarC","GabaPV_MGPChen4","PV.9norm","PV.46norm"))
levels(ReynoldsMelt$MeasureType) <- c("GabaPV MGP (BA46/BA10)\nStudy2AltarC",
                                      "GabaPV MGP (BA6)\nStudy4Chen",
                                      "GabaPV cell count (BA9) \nBeasley 2002",
                                      "GabaPV cell count (BA46) \nBeasley 2002")

#Plot group level results for counts vs MGP
p <- plotCountsMGPgroupResults(ReynoldsMelt)
ggsave("GeneralResults/GabaPV_SMRIcountsMGP.pdf", p, width = 6, height = 6, units = "in",
              dpi=300, useDingbats = FALSE)


#Plot correlations of cell counts vs MGPs
pdf("GeneralResults/GabaPV_SMRIcountsMGPcor.pdf", width = 6, height = 6, pointsize = 14, useDingbats = FALSE)
plotCountsMGPpairs(data  = Reynolds %>% select(PV.46, PV.9, GabaPV_MGPChen4, GabaPV_MGPAltarC),
                   labels = c("Cell count\nBA46\n\nBeasley 2002 " ,
                              "Cell count\nBA9\n\nBeasley 2002",
                              "MGP \nBA6\n\nStudy4Chen",
                              "MGP \nBA46/BA10\n\nStudy2AltarC"))
dev.off()

Buttner <- read.table("ButtnerAstroCountsSMRI.txt.csv", header = TRUE, sep = "\t")
names(Buttner)[1] <- "ID"
Buttner$Astro_MGPAltarC <- MetadataAltarC$Astrocyte_Genes[match(Buttner$ID, MetadataAltarC$StanleyID)]
Buttner$Astro_MGPChen4 <- MetadataChen4$Astrocyte_Genes[match(Buttner$ID, MetadataChen4$StanleyID)]
Buttner$CommonName <- MetadataAltarC$CommonName[match(Buttner$ID, MetadataAltarC$StanleyID)]
Buttner$Profile <- MetadataAltarC$Profile[match(Buttner$ID, MetadataAltarC$StanleyID)]

Buttner %<>% arrange(Astro_MGPChen4) 
Buttner$ID <- factor(Buttner$ID, levels = Buttner$ID)

Buttner %<>% mutate(GFAPnorm = rescale(Gfap.gm, c(0,1)))
Buttner %<>% filter(!is.na(.$Astro_MGPAltarC), !is.na(GFAPnorm))
ButtnerMelt <- melt(Buttner %>% select(matches("ID|CommonName|norm|MGP|Profile")), id.vars = c("ID", "CommonName", "Profile")
                     , variable.name = "MeasureType", value.name = "Value") %>% .[complete.cases(.),] %>% droplevels()

levels(ButtnerMelt$MeasureType) <- c("Astro MGP (BA46/BA10)\nStudy2AltarC",
                                      "Astro MGP (BA6)\nStudy4Chen",
                                      "Astro cell count (Cingulate cortex) \nBuettner")

#Plot group level results for counts vs MGP
plotCountsMGPgroupResults(ButtnerMelt)

ggsave("GeneralResults/Astro_SMRIcountsMGP.pdf", width = 10, height = 10, units = "in",
       dpi=300, useDingbats = FALSE)


#Plot correlations of cell counts vs MGPs
plotCountsMGPpairs(data  = Buttner %>% select(Gfap.gm, Astro_MGPChen4, Astro_MGPAltarC),
                   labels = c("Cell count\nCingulate cortex\n\nBuettner " ,
                              "MGP \nBA6\n\nStudy4Chen",
                              "MGP \nBA46/BA10\n\nStudy2AltarC"))

ggsave("GeneralResults/Astro_SMRIcountsMGPcor.pdf", width = 10, height = 10, units = "in",
       dpi=300, useDingbats = FALSE)

load("GeneralResults/studyFinalMcLeanCortex.rda")
MetadataMcLean <- studyFinal$Cortex$Metadata
