Metadata <- read.table("McLeanCortex/mclean66-design.txt", sep="\t", header=T)
Meta2 <- read.table("McLeanCortex/Metadata.csv", sep="\t", header=T)
Metadata_org <- Metadata
rownames(Metadata) <- Metadata$sample
Metadata$sample <- sapply(Metadata$sample, function(x) strsplit(as.character(x), "\\.")[[1]][2])
levels(Metadata$Diagnosis) <- c("Bipolar Disorder", "Control", "Schizoaffective", "Schizophrenia" )
levels(Metadata$DCode) <- c("BP", "Cont", "SCZ", "SCZA")
Metadata <- Metadata[order(Metadata$DCode),]

source("McLeanCortex/MaAtching samples.R")

Metadata$pH <- "NA"
Metadata$pH <- Meta2$pH[match(rownames(Metadata), sampleMatch[,6])]
names(Metadata)[names(Metadata)=="sample"] <- "Series_sample_id"
Metadata$region <- rep("Cortex", nrow(Metadata))
Metadata$Platform <- "GPL96"
write.table(Metadata, "McLeanCortex/Metadata.tsv", sep="\t", row.names=F)
rm(list=ls(pat="Meta2|match|DS|Match"))
