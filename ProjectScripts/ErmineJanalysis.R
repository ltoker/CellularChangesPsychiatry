source("SetUp.R")
# Runing ermineJ for cell type profiles.

# First steps to run ErmineJ from command line:
# Download ErmineJ 'wget http://www.chibi.ubc.ca/ermineJ/distributions/ermineJ-3.0.2-generic-bundle.zip' 
# set in the ~/.bashrc file: ERMINEJ_HOME=ermineJ-3.0.2/bin ; PATH=$ERMINEJ_HOME:$PATH, JAVA_HOME = /usr/
# Documenatation: http://erminej.chibi.ubc.ca/help/tutorials/erminej-cli/


#Getting the NeuroExpresso cortical cell expression profiles
download.file(url = "https://github.com/oganm/neuroExpressoAnalysis/blob/master/data/n_expressoSamples.rda?raw=true", destfile = "n_expressoMeta.rda")
load("n_expressoMeta.rda")
#Removes samples not included in NeuroExpresso
n_expressoSamples %<>% filter(!is.na(.$PyramidalDeep))

download.file(url = "https://github.com/oganm/neuroExpressoAnalysis/blob/master/data/n_expressoExpr.rda?raw=true", destfile = "n_expressoExpr.rda")
load("n_expressoExpr.rda")



download.file(url = "https://github.com/oganm/neuroExpressoAnalysis/blob/master/data/regionHierarchy.rda?raw=true", destfile = "regionHierarchy.rda")
load("regionHierarchy.rda")

#Scale the expression for ErmineJ analysis
n_expressoExpr<- n_expressoExpr[!grepl("\\|", n_expressoExpr$Gene.Symbol),]
NeuroExpScaled <- scale(n_expressoExpr[names(n_expressoExpr) %in% n_expressoSamples$sampleName] %>% as.matrix %>% t) %>% t %>% data.frame()
rownames(NeuroExpScaled) <- n_expressoExpr$Probe

CortexCells <- n_expressoSamples[grepl("Cortex", n_expressoSamples$Region),] %>% .$ShinyNames %>% unique

NeuroExpCellsScaled <- sapply(CortexCells, function(celltype){
  samples <- n_expressoSamples %>% filter(ShinyNames == celltype) %>% .$sampleName
  expr <- NeuroExpScaled %>%  select_(.dots = samples)
  rowMeans(expr)
}) %>% data.frame() %>% mutate(Probe = rownames(.))

names(NeuroExpCellsScaled) <- sapply(names(NeuroExpCellsScaled), function(x) gsub("\\.", "", x))
NeuroExpCellsScaled %<>% select(Probe, Astrocyte,  Microglia, Oligodendrocyte,
                                FSBasketG42, MartinottiGIN,  VIPRelnG30, PyramidalGlt25d2, PyramidalS100a10, PyramidalCrtThalamic )


NeuroExpCellsScaled2 <- apply(NeuroExpCellsScaled[-1], c(1,2), function(x) -1*x) %>% data.frame()
NeuroExpCellsScaled2$Probe <- NeuroExpCellsScaled$Probe
NeuroExpCellsScaled2 %<>% select(Probe, Astrocyte,  Microglia, Oligodendrocyte,
                                FSBasketG42, MartinottiGIN,  VIPRelnG30, PyramidalGlt25d2, PyramidalS100a10, PyramidalCrtThalamic )

write.table(NeuroExpCellsScaled, "GeneralResults/NeuroExpCellsScaled.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(NeuroExpCellsScaled2, "GeneralResults/NeuroExpCellsScaled2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

download.file(paste0("http://chibi.ubc.ca/microannots/", "GPL339",  "_noParents.an.txt.gz"), destfile="GPL339.gz")
GPL339 <- read.table("GPL339.txt", header = TRUE, sep = "\t", quote = "")

write.table(GPL339 %>% filter(ProbeName %in% n_expressoExpr$Probe), "NeuExprProbes.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


ReadErmineJ <- function(file){
  temp = readLines(file, ok = TRUE)
  DataStart <- grep("RawScore", temp)
  Descript <- temp[1:DataStart]
  ErmineJ_Names <- strsplit(temp[DataStart], "\\\t")[[1]]
  ErmineJ_DF <- vector()
  DataStart = DataStart + 1
  for(i in DataStart:length(temp)){
    ErmineJ_DF <- rbind(ErmineJ_DF, strsplit(temp[i], "\\\t")[[1]])
  } 
  
  ErmineJ_DF %<>% data.frame() 
  names(ErmineJ_DF) <- ErmineJ_Names
  ErmineJ_DF <- ErmineJ_DF[-1]
  ErmineJ_DF <- ErmineJ_DF[-nrow(ErmineJ_DF),] %>% droplevels()
  for(col in names(ErmineJ_DF)[-c(1:2)]){
    ErmineJ_DF[col] <- as.numeric(as.character(ErmineJ_DF[col] %>% unlist))
  }
  return(list(Descript = Descript,
              Result = ErmineJ_DF)) 
  }

scriptLoc = "~/CellularProportionsPsychiatry/ErmineJscript.sh"
AnnoFile = "/home/ltoker/CellularProportionsPsychiatry/NeuExprProbes.txt"
OutDir = "/home/ltoker/CellularProportionsPsychiatry/GeneralResults/"

for(cell in names(NeuroExpCellsScaled)[-1]){
  #Most upregulated genes
  ScoreFile = "/home/ltoker/CellularProportionsPsychiatry/GeneralResults/NeuroExpCellsScaled.txt"
  ScoreCol = grep(cell, names(NeuroExpCellsScaled))
  OutFile = paste0(OutDir, "ErmineJ",cell, "PC_top.txt")
  RexExp = paste(scriptLoc, ScoreFile, ScoreCol, AnnoFile, OutFile)
  system(paste("source ~/.bashrc; ERMINEJ_HOME=$ERMINEJ_HOME; cd $HOME/$ERMINEJ_HOME;", RexExp), intern = T)
  
  #Most downregulated genes
  ScoreFile = "/home/ltoker/CellularProportionsPsychiatry/GeneralResults/NeuroExpCellsScaled2.txt"
  OutFile = paste0(OutDir, "ErmineJ",cell, "PC_bottom.txt")
  RexExp = paste(scriptLoc, ScoreFile, ScoreCol, AnnoFile, OutFile)
  system(paste("source ~/.bashrc; ERMINEJ_HOME=$ERMINEJ_HOME; cd $HOME/$ERMINEJ_HOME;", RexExp), intern = T)
}



for(celltype in list.files(path = "GeneralResults/", pattern = "^ErmineJ.*PC")){
  temp <- ReadErmineJ(paste0(ResultsPath,celltype))
  assign(gsub(".txt", "", celltype), temp)
}



for(celltype in ls(pat = "^ErmineJ.*top")){
  cellData <-  eval(as.name(celltype))
  temp <- cellData$Result %>% select(Name, ID, NumGenes, RawScore,
                                     Pval, CorrectedPvalue, MFPvalue,
                                     CorrectedMFPvalue) %>% filter(CorrectedMFPvalue <= 0.05) %>% arrange(CorrectedMFPvalue) %>% head(10)
  
  if(nrow(temp) == 0){
    temp2 <- rep("---", ncol(temp)) %>% t %>% data.frame()
    names(temp2) <- names(temp)
    temp <- temp2
  }
    for(col in names(temp)){
    if(is.numeric(temp[[col]])){
      temp[[col]] <- signif(temp[[col]], digits=2)
    }
  }
  
  write.table(data.frame(gsub("ErmineJ|PC_.*", "", celltype)), append = TRUE, file = "GeneralResults/ErmineJSummaryMA.tsv", sep = "\t", col.names = FALSE, row.names = FALSE)
  write.table(temp, append = TRUE, file = "GeneralResults/ErmineJSummaryMA.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)
}

#Repeat anaysis with Tasic data
#Getting the NeuroExpresso cortical cell expression profiles
download.file(url = "https://github.com/oganm/neuroExpressoAnalysis/blob/master/data/TasicMouseExp.rda?raw=true", destfile = "TasicMouseExp.rda")
load("TasicMouseExp.rda")
download.file(url = "https://github.com/oganm/neuroExpressoAnalysis/blob/master/data/TasicMouseMeta.rda?raw=true", destfile = "TasicMouseMeta.rda")
load("TasicMouseMeta.rda")
download.file(url = "https://github.com/oganm/neuroExpressoAnalysis/blob/master/data/meltedSingleCells.rda?raw=true", destfile = "TasicNEmatch.rda")
load("TasicNEmatch.rda")
TasicMouseMeta$ShinyNames <- meltedSingleCells$ShinyNames[match(TasicMouseMeta$primary_type, meltedSingleCells$sampleName)]

TasicExpScaled <- scale(TasicMouseExp %>% as.matrix %>% t) %>% t %>% data.frame()

TasicExpCellsScaled <- sapply(unique(TasicMouseMeta$ShinyNames), function(celltype){
  samples <- TasicMouseMeta %>% filter(ShinyNames == celltype) %>% .$sample_title %>% make.names()
  expr <- TasicExpScaled %>%  select_(.dots = samples)
  rowMeans(expr)
}) %>% data.frame() %>% mutate(Probe = rownames(.))

names(TasicExpCellsScaled) <- sapply(names(TasicExpCellsScaled), function(x) gsub("\\.| ", "", x))
TasicExpCellsScaled <- TasicExpCellsScaled[-2]
TasicExpCellsScaled %<>% select(ncol(TasicExpCellsScaled), 1:c(ncol(TasicExpCellsScaled)-1))

TasicExpCellsScaled2 <- apply(TasicExpCellsScaled[-1], c(1,2), function(x) -1*x) %>% data.frame()
TasicExpCellsScaled2$Probe <- TasicExpCellsScaled$Probe
TasicExpCellsScaled2 %<>% select(ncol(TasicExpCellsScaled2), 1:c(ncol(TasicExpCellsScaled2)-1))


write.table(TasicExpCellsScaled, "GeneralResults/TasicExpCellsScaled.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(TasicExpCellsScaled2, "GeneralResults/TasicExpCellsScaled2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



scriptLoc = "~/CellularProportionsPsychiatry/ErmineJscript.sh"
AnnoFile = "/home/ltoker/CellularProportionsPsychiatry/Generic_human_noParents.an.txt"
OutDir = "/home/ltoker/CellularProportionsPsychiatry/GeneralResults/"

for(cell in names(TasicExpCellsScaled)[-1]){
  #Most upregulated genes
  ScoreFile = "/home/ltoker/CellularProportionsPsychiatry/GeneralResults/TasicExpCellsScaled.txt"
  ScoreCol = grep(cell, names(TasicExpCellsScaled))
  OutFile = paste0(OutDir, "TasicErmineJ",cell, "PC_top.txt")
  RexExp = paste(scriptLoc, ScoreFile, ScoreCol, AnnoFile, OutFile)
  system(paste("source ~/.bashrc; ERMINEJ_HOME=$ERMINEJ_HOME; cd $HOME/$ERMINEJ_HOME;", RexExp), intern = T)
  
  #Most downregulated genes
  ScoreFile = "/home/ltoker/CellularProportionsPsychiatry/GeneralResults/TasicExpCellsScaled2.txt"
  OutFile = paste0(OutDir, "TasicErmineJ",cell, "PC_bottom.txt")
  RexExp = paste(scriptLoc, ScoreFile, ScoreCol, AnnoFile, OutFile)
  system(paste("source ~/.bashrc; ERMINEJ_HOME=$ERMINEJ_HOME; cd $HOME/$ERMINEJ_HOME;", RexExp), intern = T)
}



for(celltype in list.files(path = "GeneralResults/", pattern = "^TasicErmineJ")){
  temp <- ReadErmineJ(paste0(ResultsPath,celltype))
  assign(gsub(".txt", "", celltype), temp)
}

for(celltype in ls(pat = "^TasicErmineJ.*top")){
  cellData <-  eval(as.name(celltype))
  temp <- cellData$Result %>% select(Name, ID, NumGenes, RawScore,
                                     Pval, CorrectedPvalue, MFPvalue,
                                     CorrectedMFPvalue) %>% filter(CorrectedMFPvalue <= 0.05) %>% arrange(CorrectedMFPvalue) %>% head(10)
  
  if(nrow(temp) == 0){
    temp2 <- rep("---", ncol(temp)) %>% t %>% data.frame()
    names(temp2) <- names(temp)
    temp <- temp2
  }
  for(col in names(temp)){
    if(is.numeric(temp[[col]])){
      temp[[col]] <- signif(temp[[col]], digits=1)
    }
  }
  
  write.table(data.frame(gsub("TasicErmineJ|PC_.*", "", celltype)), append = TRUE, file = "GeneralResults/ErmineJSummaryRNAseq.tsv", sep = "\t", col.names = FALSE, row.names = FALSE)
  write.table(temp, append = TRUE, file = "GeneralResults/ErmineJSummaryRNAseq.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)
}
