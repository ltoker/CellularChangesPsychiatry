source("SetUp.R")
packageF("biomaRt")
packageF("parallel")
name = "PsychEncode"

####### Preparing analysis directories #########
path = paste0(getwd(), "/")

if(length(ls(pat = "^name$")) == 0){
  print("no study name")
  browser()
}

if(!name %in% list.dirs(full.names = FALSE)){
  dir.create(name)
  dir.create(paste0(name, "/data"))
}

if(!paste0(name, "/data") %in% list.dirs(full.names = FALSE)){
  dir.create(paste0(name, "/data"))
}

resultsPath = paste0(path, name, "/")

subCol = "NeuExpRegion"

########## BioMart #########
ensembl <- useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl")

geneNames <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), mart = ensembl)
MitoGenes <- geneNames[grepl("MT-", geneNames$hgnc_symbol),]
###### Getting the metadata#########
Meta1 <- read.table("PEC_BrainGVEX_ATACseq_metadata.tsv", header = T, sep = "\t") %>% filter(!is.na(PMI)) %>% droplevels()
Meta1Stanley <- Meta1 %>% filter(BrainBank != "BSHRI") %>% droplevels() 

write.table(Meta1Stanley$RNASeq.BID, "PsychEncodeIDsStanley.tsv", sep = "\t", row.names = FALSE, col.names = FALSE)

Meta2 <- read.table("UIC-UChicago-U01MH103340_BrainGVEX_ClinicalMetadata_AgeCensored_August2016Release.tsv", header = T, sep = "\t")

Meta1Stanley$pH <- Meta2$pH[match(Meta1Stanley$RNASeq.BID, Meta2$Individual_ID..RNAseq.library.BID..Primary.ID.)]

RNAseqMeta <- read.table("UIC-UChicago-U01MH103340_BrainGVEX_RNAseqMetadata_August2016Release.txt", header = T, sep = "\t")
RNAseqMeta$Individual_ID..RNAseq.Library.BID..Synapse. <- sapply(RNAseqMeta$File_Name..on.Synapse..if.available, function(x){
  strsplit(as.character(x), "_")[[1]][7]
})

RNAseqMetaStanley <- RNAseqMeta %>%
  filter(Individual_ID..RNAseq.Library.BID..Synapse. %in% as.character(Meta1Stanley$RNASeq.BID)) %>%
  droplevels()

RNAseqMetaStanley %<>% mutate_at(grep("Reads", names(RNAseqMetaStanley), value = T), as.character)
RNAseqMetaStanley %<>% mutate_at(grep("Reads", names(RNAseqMetaStanley), value = T), as.numeric)
RNAseqMetaStanley$rRNARate <- sapply(RNAseqMetaStanley$rRNARate, function(x){
  x <- gsub("%", "", as.character(x))
  as.numeric(x)
})

Meta1Stanley$SequencingPlatform <- RNAseqMetaStanley$SequencingPlatform[match(Meta1Stanley$RNASeq.BID,
                                                                              RNAseqMetaStanley$Individual_ID..RNAseq.Library.BID..Synapse.)]
Meta1Stanley$rRNArate <- RNAseqMetaStanley$rRNARate[match(Meta1Stanley$RNASeq.BID,
                                                          RNAseqMetaStanley$Individual_ID..RNAseq.Library.BID..Synapse.)]

Metadata <- Meta1Stanley %>% select(-ATACSeq.BID)
Metadata %<>% mutate(NeuExpRegion = "Cortex")
Metadata$Profile <- sapply(as.character(Metadata$Diagnosis), function(x){
  if(x == "BP"){
    "BP"
  } else if (x == "SCZ"){
    "SCZ"
  } else if(x == "Control"){
    "Cont"
  } else {
    NA
  }
}) %>% factor(levels = c("Cont", "BP", "SCZ"))

SampleNum <- Metadata %>% group_by(Profile) %>% summarise(n = n()) %>% data.frame()

Metadata %<>% arrange(Profile) %<>% droplevels()

Metadata$CommonName <- apply(SampleNum, 1, function(x){
  paste0(x[1], "_", 1:x[2])
}) %>% unlist

Meta

######################
