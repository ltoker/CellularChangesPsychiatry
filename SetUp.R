#Set up the path for general scripts
GenScriptPath = "/home/ltoker/Rscripts/"

#Set up the path for project specific scripts
ProjScriptPath = paste0(getwd(),"/ProjectScripts/")

#Set up the path for results directory
if(!"GeneralResults" %in% list.dirs(full.names = FALSE)){
  dir.create("GeneralResults")
}
GeneralResultsPath = paste0(getwd(), "/GeneralResults/")

# Source required functions
source(paste0(GenScriptPath,"general_functions.R"))
source(paste0(ProjScriptPath,"projectFunc.R"))

packageF("sva")
packageF("gplots")
packageF("scales")
select = dplyr::select
filter = dplyr::filter
mutate = dplyr::mutate
