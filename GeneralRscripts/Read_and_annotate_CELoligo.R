packageF("oligoClasses")

if(isPackageLoaded(pkg="affyPLM")){
  detach("package:affyPLM",force=T)
}

if(isPackageLoaded(pkg="affy")){
  detach("package:affy", force=T)
}

packageF("oligo")
packageF("pd.hugene.1.0.st.v1")
#packageF("arrayQualityMetrics")

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


# 1 Read in probe level data
# read the data. In order to read all the CEL files in the wd, use (), otherwise - indicate the cel_files

data <- read.celfiles(paste0(path, CelFiles))

#RMA - normalization core transcript level
geneCore <- rma(data, target="core")

#Retrieving NetAffx Biological Annotation
featureData(geneCore) <- getNetAffx(geneCore, "transcript")

#Extract the expression data
exp_value <- get("exprs", pos=assayData(geneCore))


                       

