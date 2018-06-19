source("SetUp.R")
packageF("parallel")

packageF("biomaRt")
packageF("grid")
packageF("pheatmap")
packageF("RColorBrewer")

## Adjustment for pheatmap
draw_colnames_90 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}

assignInNamespace(x="draw_colnames", value="draw_colnames_90",
                  ns=asNamespace("pheatmap"))


Count2CPM <- function(countData){
  apply(countData, 2, function(smp){
    TotalCount = sum(smp)
    (10^6)*smp/TotalCount
  })
}

GetCountMatrix <- function(DataLocm, MitoGenes){
  countMatrix <- read.csv(paste0(DataLoc, "countMatrix.genes"), header=TRUE, sep = "\t", quote = "")
  CountSum <- apply(countMatrix[grepl("GSM", colnames(countMatrix))], 2, sum)
  MitoCountSum <- apply(countMatrix[countMatrix$genes %in% MitoGenes$ensembl_gene_id,
                                    grepl("GSM", colnames(countMatrix))], 2, sum)
  MitoCountFiltered <- countMatrix[!countMatrix$genes %in% MitoGenes$ensembl_gene_id,]
  CPMmatrix <- Count2CPM(MitoCountFiltered[-1])
  #log2 transformation
  log2CPMmatrix <- apply(CPMmatrix, c(1,2), function(x) {log2(x+1)})
  GeneSymbol <- geneNames$hgnc_symbol[match(MitoCountFiltered$genes, geneNames$ensembl_gene_id)]
  ExpData <- data.frame(GeneSymbol = GeneSymbol,
                        Probe = GeneSymbol,
                        ensemblID = MitoCountFiltered$genes)
  ExpDataCPM <- cbind(ExpData, CPMmatrix)
  ExpDataLog2CPM <- cbind(ExpData, log2CPMmatrix)
  return(list(ExpDataCPM = ExpDataCPM,
              ExpDataLog2CPM = ExpDataLog2CPM))
}

GetLenghStats <- function(expData, geneNameMap = geneNames, countThresh = 10){
  expData$Mean <- apply(expData %>% select(matches("GSM")), 1, mean)
  expData$SD <- apply(expData %>% select(matches("GSM")), 1, sd)
  GeneLength <- sapply(grep("GSM", names(expData), value = T), function(GSM){
    data <- read.csv(paste0(DataLoc, GSM, ".genes.results"), header=TRUE, sep = "\t", quote = "")
    data %>% select(matches("gene|length"))
    }, simplify = FALSE)

  #Get the ensembl gene length
  GeneLengthReal <- lapply(GeneLength, function(GSM){
    GSM$length
    }) %>% do.call(cbind, .) %>% data.frame()

  #Get the effective gene length
  GeneLengthEffect <- lapply(GeneLength, function(GSM){
    GSM$effective_length
    }) %>% do.call(cbind, .) %>% data.frame()

  #Add the ensembl and hugu gene name
  EnsemblGene <- as.character(GeneLength[[1]]$gene_id)
  GeneLengthEffect <- cbind(EnsemblGene, GeneLengthEffect)
  GeneLengthEffect$GeneSymbol <- geneNameMap$hgnc_symbol[match(GeneLengthEffect$EnsemblGene,
                                                               geneNameMap$ensembl_gene_id)]

  #Calculate mean and sd of the gene length
  GeneLengthEffect$lengthMean <- apply(GeneLengthEffect %>% select(matches("GSM")), 1, mean)
  GeneLengthEffect$lengthSD <- apply(GeneLengthEffect %>% select(matches("GSM")), 1, sd)

  #Add mean and sd of gene CPM
  GeneLengthEffect$countMean <- expData$Mean[match(GeneLengthEffect$EnsemblGene, expData$ensemblID)]
  GeneLengthEffect$countSD <- expData$SD[match(GeneLengthEffect$EnsemblGene, expData$ensemblID)]

  #Get the number of isoforms for each gene
  Isoforms <- read.csv(paste0(DataLoc, names(GeneLength)[1], ".isoforms.results"), header=TRUE, sep = "\t", quote = "")
  Isoforms %<>% filter(gene_id %in% GeneLengthEffect$EnsemblGene) %>% droplevels()
  IsoNum <- group_by(Isoforms, gene_id) %>% summarise(n = n()) %>% data.frame()
  
  GeneLengthEffect$TransNum <- IsoNum$n[match(GeneLengthEffect$EnsemblGene, IsoNum$gene_id)]
  GeneLengthEffect %<>% select_(.dots = c("EnsemblGene", "GeneSymbol", "TransNum", "lengthMean", "lengthSD", "countMean", "countSD",
                                          grep("GSM", names(GeneLengthReal), value = T)))

  #Remove genes with no Gene symbols
  GeneLengthEffect %<>% filter(!is.na(GeneSymbol)) %>% droplevels()
  GeneLengthEffect <- GeneLengthEffect[!GeneLengthEffect$GeneSymbol == "",]

  #Remove miRs and mitochondria-encoded genes
  GeneLengthEffect <- GeneLengthEffect[!grepl("^MIR|MT-|SNORD|^RN7", GeneLengthEffect$GeneSymbol),]

  #Remove genes with short length
  GeneLengthEffect %<>% filter(lengthMean > 100) %>% droplevels()

  #Log trandform cmp and length measures
  GeneLengthEffect %<>% mutate(length_mean = log10(lengthMean),
                               length_sd = log10(lengthSD),
                               CPM_mean = log2(countMean + 1),
                               CPM_sd = log2(countSD + 1))

  names(GeneLengthEffect)[grepl("length_", names(GeneLengthEffect))] <- paste0("log10(",grep("length_", names(GeneLengthEffect), value = T),")")

  names(GeneLengthEffect)[grepl("CPM", names(GeneLengthEffect))] <- paste0("log2(", grep("CPM", names(GeneLengthEffect), value = T),")")
  
  GeneLengthEffect %<>% mutate(lengthCV = lengthSD/lengthMean,
                               countCV = countSD/countMean)

  #Separate to genes with high and low expression, and with one isoform
  GeneLengthEffect_low <- GeneLengthEffect %>% filter(countMean <= countThresh, TransNum > 1)
  GeneLengthEffect_high <- GeneLengthEffect %>% filter(countMean > countThresh,  TransNum > 1)
  
  return(list(GeneLengthReal = GeneLengthReal,
              GeneLengthEffect = GeneLengthEffect,
              GeneLengthEffect_low = GeneLengthEffect_low,
              GeneLengthEffect_high = GeneLengthEffect_high))
}

GetGeneTsof <- function(EnsGene, expData){
  IsoPercnt <- sapply(grep("GSM", names(expData), value = T), function(GSM){
    data <- read.csv(paste0(DataLoc, GSM, ".isoforms.results"), header=TRUE, sep = "\t", quote = "") %>% .[.$gene_id == EnsGene,]
  }, simplify = FALSE)
  browser()
}

ensembl <- useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl")

geneNames <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), mart = ensembl)

MitoGenes <- geneNames[grepl("MT-", geneNames$hgnc_symbol),]

DataLoc = "/space/grp/Pipelines/rnaseq-pipeline/Pipelines/rsem/quantified/GSE80655/"

##############################################################33

MetaData <- read.table("GSE80655/Metadata.tsv", header = T, sep = "\t")

MetaDataDLPFC <- MetaData %>% filter(brain.region == "DLPFC") %>% droplevels()
MetaDataAnCg <- MetaData %>% filter(brain.region == "AnCg") %>% droplevels()
MetaDataNAcc <- MetaData %>% filter(brain.region == "nAcc") %>% droplevels()

ExpDataList <- GetCountMatrix(DataLoc = DataLoc, MitoGenes = MitoGenes)

ExpDataDLPFC <- ExpDataList$ExpDataCPM %>% select_(.dots = c("GeneSymbol", "Probe", "ensemblID", levels(MetaDataDLPFC$Series_sample_id)))
ExpDataAnCg <- ExpDataList$ExpDataCPM %>% select_(.dots = c("GeneSymbol", "Probe", "ensemblID", levels(MetaDataAnCg$Series_sample_id)))
ExpDataNAcc <- ExpDataList$ExpDataCPM %>% select_(.dots = c("GeneSymbol", "Probe", "ensemblID", levels(MetaDataNAcc$Series_sample_id)))

LengthSumListDLPFC <- GetLenghStats(expData = ExpDataDLPFC)
LengthSumListAnCg <- GetLenghStats(expData = ExpDataAnCg)
LengthSumListNAcc <- GetLenghStats(expData = ExpDataNAcc)


DataLoc = "/space/grp/Pipelines/rnaseq-pipeline/Pipelines/rsem/quantified/GSE65540/"
ExpDataGSE65540 <- GetCountMatrix(DataLoc = DataLoc, MitoGenes = MitoGenes)

LengthSumListGSE65540 <- GetLenghStats(expData = ExpDataGSE65540$ExpDataCPM)


#Create merged data frame for length SD an CV
mergedLengthSum <- merge(LengthSumListDLPFC$GeneLengthEffect_high %>% select(GeneSymbol, lengthSD, lengthCV, countSD, countCV, countMean, TransNum),
                         LengthSumListAnCg$GeneLengthEffect_high %>% select(GeneSymbol, lengthSD, lengthCV, countSD, countCV, countMean, TransNum),
                         by = "GeneSymbol", all = TRUE, suffixes = c("_GSE80655_DLPFC", "_GSE80655_AnCg"))

temp <-  merge(LengthSumListNAcc$GeneLengthEffect_high %>% select(GeneSymbol, lengthSD, lengthCV, countSD, countCV, countMean, TransNum),
               LengthSumListGSE65540$GeneLengthEffect_high %>% select(GeneSymbol, lengthSD, lengthCV, countSD, countCV, countMean, TransNum),
               by = "GeneSymbol", all = TRUE, suffixes = c("_GSE80655_NAcc", "_GSE65540"))

mergedLengthSum <- merge(mergedLengthSum, temp, by = "GeneSymbol", all = T)


CorDF <- cor(mergedLengthSum %>% select(-matches("GeneSymbol")), use = "na.or.complete", method = "spearman") %>% round(digits = 2)
diag(CorDF) <- NA

CorDFpart <- cor(mergedLengthSum %>% select(-matches("GeneSymbol|CV|countSD")), use = "na.or.complete", method = "spearman") %>% round(digits = 2)
diag(CorDFpart) <- NA

pheatmap(as.matrix(CorDF), filename = "/home/ltoker/LengthSDcorFull.png")
pheatmap(as.matrix(CorDFpart), filename = "/home/ltoker/LengthSDcor.png")

mergedLengthSumLow <- merge(LengthSumListDLPFC$GeneLengthEffect_low %>% select(GeneSymbol, lengthSD, lengthCV, countSD, countCV, countMean, TransNum),
                         LengthSumListAnCg$GeneLengthEffect_low %>% select(GeneSymbol, lengthSD, lengthCV, countSD, countCV, countMean, TransNum),
                         by = "GeneSymbol", all = TRUE, suffixes = c("_GSE80655_DLPFC", "_GSE80655_AnCg"))

temp <-  merge(LengthSumListNAcc$GeneLengthEffect_low %>% select(GeneSymbol, lengthSD, lengthCV, countSD, countCV, countMean, TransNum),
               LengthSumListGSE65540$GeneLengthEffect_low %>% select(GeneSymbol, lengthSD, lengthCV, countSD, countCV, countMean, TransNum),
               by = "GeneSymbol", all = TRUE, suffixes = c("_GSE80655_NAcc", "_GSE65540"))

mergedLengthSumLow <- merge(mergedLengthSumLow, temp, by = "GeneSymbol", all = T)


CorDFlow <- cor(mergedLengthSumLow %>% select(-matches("GeneSymbol")), use = "na.or.complete", method = "spearman") %>% round(digits = 2)
diag(CorDFlow) <- NA

CorDFpartlow <- cor(mergedLengthSumLow %>% select(-matches("GeneSymbol|CV|countSD")), use = "na.or.complete", method = "spearman") %>% round(digits = 2)
diag(CorDFpartlow) <- NA

pheatmap(as.matrix(CorDFlow), filename = "/home/ltoker/LengthSDcorFullLow.png")
pheatmap(as.matrix(CorDFpartlow), filename = "/home/ltoker/LengthSDcorLow.png")

DataLoc = "/space/grp/Pipelines/rnaseq-pipeline/Pipelines/rsem/quantified/GSE80655/"
IsoRYR2 <- GetGeneTsof(EnsGene = "ENSG00000198626", expData = ExpDataDLPFC)
IsoNEAT1 <- GetGeneTsof(EnsGene = "ENSG00000245532", expData = ExpDataDLPFC)
IsoFAM208B <- GetGeneTsof(EnsGene = "ENSG00000108021", expData = ExpDataDLPFC)
IsoCDH6 <- GetGeneTsof(EnsGene = "ENSG00000124177", expData = ExpDataDLPFC)

pheatmap(as.matrix(CorDF), color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList)

#Create annotation file to run ermineJ
AnnoFileGemma <- read.table("Generic_human_noParents.an.txt", sep = "\t", header = TRUE, comment = "#", quote = "")

AnnoFileData <- AnnoFileGemma[match(GeneLengthReal$GeneSymbol, AnnoFileGemma$GeneSymbols),]
AnnoFileData$ProbeName <- GeneLengthReal$EnsemblGene

rownames(GeneLengthEffect) <- GeneLengthEffect$EnsemblGene
rownames(GeneLengthReal) <- GeneLengthReal$EnsemblGene

GeneLengthEffectMelt <- melt(GeneLengthEffect[,!grepl("GSM", names(GeneLengthEffect))],
                             id.vars = c("EnsemblGene", "GeneSymbol", "lengthSD"), value.name = "Value", variable.name = "Measure")

ggplot(GeneLengthEffectMelt, aes(lengthSD, Value)) +
  theme_classic() +
  labs(y = "") +
  geom_point(alpha = 0.2, size = 0.3) +
  facet_wrap(~Measure, scales = "free_y")

#Run ErmineJ
if(!"ermineR" %in% rownames(installed.packages())){
  install_github("oganm/ermineR", force = T)
}
library(ermineR)

tempSDhigh <- gsr(scores = GeneLengthEffect, scoreColumn = "SD", bigIsBetter = T, logTrans = F, annotation = AnnoFileData, aspects = "B")
tempSDlow <- gsr(scores = GeneLengthEffect, scoreColumn = "SD", bigIsBetter = F, logTrans = F, annotation = AnnoFileData, aspects = "B")


tempSDhigh$results %>% data.frame() %>% select(matches("Name|Pvalue")) %>% arrange(CorrectedPvalue) %>% head(20)
tempSDlow$results %>% data.frame() %>% select(matches("Name|Pvalue")) %>% arrange(CorrectedPvalue) %>% head(20)

plot(density(log10(GeneLengthEffect$SD)))

GeneLengthEffect %>% arrange(SD) %>% .[1:4] %>% head(20)
GeneLengthEffect %>% arrange(desc(SD)) %>% .[1:4] %>% head(20)


