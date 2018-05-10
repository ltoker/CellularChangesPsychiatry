source("SetUp.R")
source(paste0(GenScriptPath, "Cell_type_PCA.R"))

load("AnalysisStudies.rda")

PCAresults <- list()
for(study in AnalysisStudies){
  load(paste0(GeneralResultsPath, "PCAresults", study, ".Rda"))
  if(grepl("Stanley", study)){
    if(study == "StanleyArray"){
      PCA_results = PCA_results[!names(PCA_results) == "Study17Laeng"]
      } else if(study == "StanleyConsortium") {
        PCA_results = PCA_results[!names(PCA_results) == "Study16Kemether"]
        }
    for(study2 in names(PCA_results)){
      PCAresults[[paste0("PCAresults", study2)]] <-  PCA_results[[study2]]
      }
  } else if (study == "GSE80655"){
    PCAresults[[paste0("PCAresults", study)]] <- PCA_results$DLPFC
    } else {
      PCAresults[[paste0("PCAresults", study)]] <- PCA_results$Cortex
    }
}

rm(PCA_results, study, study2)
save(PCAresults, file = paste0(GeneralResultsPath, "PCAresultsList.rda"))

VarExplained <- lapply(PCAresults, function(study){
  temp <- sapply(1:length(study), function(x){
    study[[x]]$All$Astrocyte_Genes %>%
      summary() %>% .$importance %>% data.frame() %>%
      .[2,1:3] %>% unlist
  }) %>% t %>% data.frame()
  
  temp2 <- sapply(1:length(study), function(x){
    study[[x]]$All$GabaPV_Genes %>%
      summary() %>% .$importance %>% data.frame() %>%
      .[2,1:3] %>% unlist
  }) %>% t %>% data.frame()
  list(Astro = temp, GabaPV = temp2)
})

VarExplainedSum <- lapply(VarExplained, function(study){
  lapply(study, function(celltype){
    apply(celltype, 2, function(x) {percent(mean(x))})
  })
})


GabaPVvarDF <- lapply(VarExplainedSum, function(study){
  study$GabaPV %>% t %>% data.frame()
}) %>% rbindlist() %>% data.frame %>%
  mutate(Dataset = sapply(names(VarExplainedSum), function(x) gsub("PCAresults", "", x))) %>%
  select(4, 1:3)

write.table(GabaPVvarDF, paste0(GeneralResultsPath, "GabaPVvarDF.tsv"),
            sep = "\t", row.names = F, col.names = T )

AstrocytevarDF <- lapply(VarExplainedSum, function(study){
  study$Astro %>% t %>% data.frame()
}) %>% rbindlist() %>% data.frame %>%
  mutate(Dataset = sapply(names(VarExplainedSum), function(x) gsub("PCAresults", "", x))) %>%
  select(4, 1:3)
write.table(AstrocytevarDF, paste0(GeneralResultsPath, "AstrovarDF.tsv"),
            sep = "\t", row.names = F, col.names = T )
