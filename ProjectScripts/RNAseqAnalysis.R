source("SetUp.R")
source("parallel")

packageF("biomaRt")

Count2CPM <- function(countData){
  apply(countData, 2, function(smp){
    TotalCount = sum(smp)
    (10^6)*smp/TotalCount
    })
  }


ensembl <- useMart(biomart = "ensembl", dataset="hsapiens_gene_ensembl")

geneNames <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), mart = ensembl)

MitoGenes <- geneNames[grepl("MT-", geneNames$hgnc_symbol),]
countMatrix <- read.csv("GSE80655/data/countMatrix.genes", header=TRUE, sep = "\t", quote = "")
CountSum <- apply(countMatrix[grepl("GSM", colnames(countMatrix))], 2, sum)

MitoCountSum <- apply(countMatrix[countMatrix$genes %in% MitoGenes$ensembl_gene_id,
                                  grepl("GSM", colnames(countMatrix))], 2, sum)

MitoCountFiltered <- countMatrix[!countMatrix$genes %in% MitoGenes$ensembl_gene_id,]

CPMmatrix <- Count2CPM(MitoCountFiltered[-1])

#log2 transformation
CPMmatrix <- apply(CPMmatrix, c(1,2), function(x) {log2(x+1)})

GeneSymbol <- geneNames$hgnc_symbol[match(MitoCountFiltered$genes, geneNames$ensembl_gene_id)]

ExpData <- data.frame(GeneSymbol = GeneSymbol,
                      Probe = GeneSymbol,
                      ensemblID = MitoCountFiltered$genes)

ExpData <- cbind(ExpData, CPMmatrix) 












##############################################################33
IsoformGSM2132201 <- read.csv("GSE80655/data/GSM2132201.isoforms.results", header=TRUE, sep = "\t", quote = "")
IsoformGSM2132202 <- read.csv("GSE80655/data/GSM2132202.isoforms.results", header=TRUE, sep = "\t", quote = "")
GeneGSM2132201 <- read.csv("GSE80655/data/GSM2132201.genes.results", header=TRUE, sep = "\t", quote = "")
GeneGSM2132202 <- read.csv("GSE80655/data/GSM2132202.genes.results", header=TRUE, sep = "\t", quote = "")

ENSG00000000003_GSM2132201length <- apply(IsoformGSM2132201 %>% filter(gene_id == "ENSG00000000003"), 1, function(isoform){
  effLength = isoform[["effective_length"]] %>% as.numeric()
  Pct = isoform[["IsoPct"]] %>% as.numeric()
  effLength*Pct/100
})

ENSG00000000003_GSM2132202length <- apply(IsoformGSM2132202 %>% filter(gene_id == "ENSG00000000003"), 1, function(isoform){
  effLength = isoform[["effective_length"]] %>% as.numeric()
  Pct = isoform[["IsoPct"]] %>% as.numeric()
  effLength*Pct/100
})


GeneGSM2132201 %<>% mutate(RPK =expected_count/effective_length)
AllRPKGSM2132201 = sum(GeneGSM2132201$RPK, na.rm = T)
AllCountGSM2132201 = sum(GeneGSM2132201$expected_count)
