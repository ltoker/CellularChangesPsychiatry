setwd("/home/ltoker/CellularProportionsPsychiatry")
GenScriptPath = "/home/ltoker/Rscripts/"
source(paste0(GenScriptPath,"general_functions.R"))
source("projectFunc.R")
packageF("sva")
packageF("gplots")
packageF("scales")

packageF("parallel")

OkatyMeta <- read.table("FSmaturation_GSE17806/1483_GSE17806_expdesign.data.txt", header=T, sep="\t", quote = "")
OkatyMeta$AgeNumeric <- sapply(OkatyMeta$Age, function(x) gsub("P", "", x) %>% as.numeric )
OkatyMeta %<>% arrange(AgeNumeric) 
OkatyMeta$Age <- factor(OkatyMeta$Age, levels = unique(OkatyMeta$Age))
OkatyMeta$Bioassay <- sapply(OkatyMeta$Bioassay, function(x) make.names(x))

OkatyExp <- read.table("FSmaturation_GSE17806/1483_GSE17806_expmat.data.txt", header=T, sep="\t", quote = "")
OkatyExp %<>% select(-Sequence, -GemmaId, -NCBIid)
OkatyExp <- cbind(OkatyExp[c(1:3)],OkatyExp[match(OkatyMeta$Bioassay, names(OkatyExp))])
rownames(OkatyExp) <- OkatyExp$Probe



OkatyExpMelt <- melt(OkatyExp, id.vars = c("Probe", "GeneSymbol"),
                     measure.vars = grep("GSE", names(OkatyExp), value=T),
                     variable.name = "SampleName", value.name = "Exp")
OkatyExpMelt$Age <- OkatyMeta$Age[match(OkatyExpMelt$SampleName, OkatyMeta$Bioassay)] %>% ordered()

temp <- group_by(OkatyExpMelt, Age, Probe) %>% summarise(Mean = mean(Exp)) %>% data.frame
temp <- group_by(temp, Probe) %>% summarize(Max = max(Mean))
ProbeRM <- temp %>% filter(temp$Max < 6) %>% data.frame

OkatyExpMelt %<>% filter(!Probe %in% ProbeRM$Probe)


lmResults <- mclapply(unique(OkatyExpMelt$Probe),function(probe){ 
  lm.mod <- lm(Exp ~ 1 + Age, data = OkatyExpMelt %>% filter(Probe == probe))
  anv <- anova(lm.mod)
  cbind(t(lm.mod$coefficients), anv["Age",])
},mc.cores = detectCores())

names(lmResults) <- OkatyExpMelt$Probe %>% unique()
lmResults <-  rbindlist(lmResults, use.names = TRUE, idcol="Probe")
lmResults$BH <- p.adjust(lmResults$`Pr(>F)`, method = "BH")
lmResults %<>% arrange(BH)

save(lmResults, file = paste0("DevelopmentLmResults.rda"))
