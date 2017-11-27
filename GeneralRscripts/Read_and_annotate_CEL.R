packageF("affy")
packageF("affyPLM")
packageF("arrayQualityMetrics")
packageF ("hgu133a.db")
packageF("hgu133acdf")
packageF ("hgu133a.db")
packageF("hgu133plus2.db")
packageF ("hgu133plus2.db")

add_symbols <- function(data, Anno_file){
  ColNum <- ncol(data)
  A <- c(ColNum +1) ; B <- A+2
  data <- as.data.frame(data)
  data$Probe <- rownames(data)
  data$GeneSymbol <- Anno_file$GeneSymbols[match(rownames(data), Anno_file$ProbeName)]
  data$Annotation <- Anno_file$GeneNames[match(rownames(data), Anno_file$ProbeName)]
  data <- data[,c(A:B, 1:ColNum)]
  data <- data[which(data$GeneSymbol != ""),]
  return(data)
}

download.file("http://www.chibi.ubc.ca/Gemma/arrays/downloadAnnotationFile.html?id=4&fileType=noParents", destfile="GPL570.gz")
Anno_file <- read.table("GPL570.gz", comment="#", header=T, quote='"', sep="\t")


# 1 Read in probe level data
# read the data. In order to read all the CEL files in the wd, use (), otherwise - indicate the cel_files

affydata <- ReadAffy(celfile.path=path)


# raw expression data
ed <- exprs(affydata)

samp <- sampleNames(affydata)
probes <- featureNames(affydata)


# 2 Normalizing Data   
#
# The Affy package has implementations of a number of normalization methods
# for single-channel arrays. This includes (among others):
#   - mas5() - Affymetrix's own normalization program
#   - rma() - 'Robust Multi-Chip' average
#   - gcrma() - A bias-corrected RMA
# GCRMA is good but takes significantly longer than RMA, so RMA is the
# most commonly used


nvals <- rma(affydata)

# normalised expression data
ned <- exprs(nvals)

nsamp <- sampleNames(nvals)
nprobes <- featureNames(nvals)

#Add Gene symbol and annotation
aned <- add_symbols(ned, Anno_file)
names(aned) <- sapply(names(aned), function(x) strsplit(x, "\\.")[[1]][1])

#Identify scanDate
scanDate <- protocolData(affydata)$ScanDate
scanDate <- sapply(scanDate, function(x) strsplit(x, " ")[[1]][1])
names(scanDate) <- rownames(pData(affydata))
names(scanDate) <- sapply(names(scanDate), function(x) strsplit(x, "\\.")[[1]][1])



                       

