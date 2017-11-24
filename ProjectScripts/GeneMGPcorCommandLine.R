setwd("/home/ltoker/CellularProportionsPsychiatry")
GenScriptPath = "/home/ltoker/Rscripts/"
ProjScriptPath = "ProjectScripts/"
source(paste0(GenScriptPath,"general_functions.R"))
source(paste0(ProjScriptPath,"projectFunc.R"))
packageF("sva")
packageF("gplots")
packageF("scales")
detach("package:dplyr", unload=TRUE)
packageF("dplyr")
packageF("doParallel")

MistryGenesUP <- read.table("MistryUp.txt", header = T, sep="\t") 
MistryGenesDown <- read.table("MistryDown.txt", header = T, sep="\t") 

registerDoSEQ()
registerDoParallel(detectCores())
MGPcorAllgenes <- lapply(studyFinal, function(std){
  data <- std$aned_high
  rownames(data) <- data$GeneSymbol
  rownames(data) <- sapply(rownames(data), function(x) make.names(x))
  meta <- std$Metadata
  meta <- meta[apply(meta, 2, function(x) sum(!is.na(x))) > 0]
  names(meta) <- sapply(names(meta), function(x) gsub(" ", "", x))
  groups = levels(meta$Profile)
  allMGP<- names(meta %>% select(matches("_Genes")))
  temp <- foreach(mgp = allMGP) %dopar% {
    dataMGP <- meta %>% select_(.dots = c("CommonName", "Profile", mgp))
    rownames(dataMGP) <- dataMGP$CommonName %>% as.character
    GeneExp <- t(data[sapply(names(data), function(x) is.numeric(data[[x]]))])
    temp <- merge(dataMGP, GeneExp, by = "row.names")
    GeneMGPcor <- sapply(groups, function(grp){
      tempCor <- sapply(names(temp)[names(temp) %in% rownames(data)], function(gene){
        subdata = temp %>% select_(.dots = c("Profile", mgp, gene)) %>% filter_(paste0("Profile","=='", grp,"'"))
        cor(subdata[[mgp]], subdata[[gene]], use="complete.obs", method = "spearman")
      }, simplify=TRUE) %>% data.frame
      names(tempCor) <- "Cor"
      tempCor
    }, simplify=FALSE) 
  }
  names(temp) <- allMGP
  temp
})
registerDoSEQ()
save(MGPcorAllgenes, file = paste0("GeneralResults/MGPcorAllgenes_", name, ".rda"))