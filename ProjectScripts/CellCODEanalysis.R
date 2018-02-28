packageF("dplyr")

packageF("purrr")

packageF("magrittr")

install_github("oganm/markerGeneProfile", force = T)
install_github("oganm/CellCODE", force = T)
library(CellCODE)

data('mouseMarkerGenesCombined')



# this function will turn gene lists into CellCODE inputs

cellCodeInput = function(genes){
  
  geneCount = genes %>% map(length)
  
  allGenes = genes %>% unlist %>% unique
  
  
  
  tags = genes %>% lapply(function(x){
    
    out = rep(0,length(allGenes))
    
    out[allGenes %in% x] = 1
    
    return(out)
    
  }) %>% as.data.frame()
  
  
  
  rownames(tags) = allGenes
  
  return(tags)
  
}



# this is marker gene input

geneList = mouseMarkerGenesCombined$Cortex %>% map(mouse2human) %>% map('humanGene')
geneList$GabaPV <- geneList$GabaPV[!geneList$GabaPV %in% c("WIF1", "TMEM132C", "BTN2A2")]

# remove genes that are not in the dataset

load("GeneralResults/studyFinalGSE35978.rda")
ExpData <- studyFinal$Cortex$aned_high %>%
  select_(.dots = as.character(studyFinal$Cortex$Metadata$CommonName)) %>% as.matrix()
rownames(ExpData) <- studyFinal$Cortex$aned_high$GeneSymbol

geneList %<>% map(function(x){
  
  x[x %in% rownames(ExpData)]
  
})

# code works with gene count > 2. Not sure why but it is a decent cutoff anyway

geneList = geneList[sapply(geneList,length)>2]

tagList = cellCodeInput(geneList)



# dataset input is a MATRIX (not data frame) with genes as rownames

groups = studyFinal$Cortex$Metadata$Profile



SPVs=getAllSPVs(data = ExpData, 
                
                grp = groups %>% as.factor,
                
                dataTag = tagList, 
                
                method = "mixed", 
                
                plot =  T) 
