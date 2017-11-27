source(paste0(GenScriptPath,"graphic_functions.R")) #networks and heatmaps
source(paste0(GenScriptPath,"correlation_functions.R")) #corrlations and distributions functions

Metadata %<>% droplevels()
if(length(ls(pat="resultsPath", envir = .GlobalEnv)) == 0){
  print("No results path specified, using working directory")
  resultsPath = paste0(getwd(), "/")
}
if(length(ls(pat="study")) == 0){
  study = ""
}
if(length(ls(pat="region")) == 0){
  region = ""
}

aned_good <- aned[complete.cases(aned),]
#Get max signal for noise data based on expression of known non-expressed genes
GeneGender <- GeneSex(aned=aned_good,
                      Metadata=data.frame(CommonName = names(aned_good[sapply(aned_good,
                                                                              function(x) is.numeric(x))])))
Fgene <- grep("XIST", aned_good$GeneSymbol, value=TRUE) %>% unique
Mgene <- grep("KDM5D|RPS4Y1", aned_good$GeneSymbol, value=TRUE)

if(length(c(Fgene, Mgene)) > 0){
  Noise <- sapply(c(Fgene, Mgene), function(gene){
    if(gene %in% Fgene){
      gender = "M"
    } else if (gene %in% Mgene) {
      gender = "F"
    }
    
    temp <- aned_good %>% filter(GeneSymbol == gene) %>%
      select_(.dots=GeneGender %>% filter(BioGender == gender) %>% .$CommonName %>% as.character)  %>% unlist
    #quantile(temp, 0.9)
  })  %>% unlist
  
  MaxNoise <- quantile(Noise, 0.95)
} else {
  MaxNoise = 6
}
print(paste("Noise threshold:", MaxNoise))

#Define genes above and below noise threshold
ProbeSum <- apply(aned_good[,-c(1:3)], 1, function(x) quantile(x, 0.95) > MaxNoise)
aned_high <- aned_good[ProbeSum,]
aned_low <- aned_good[!ProbeSum,]


#Create correlation histograms for original high and low signals
#png(paste0(resultsPath, study, "_High and Low correlations before and after.png") ,width = 980, height = 700, units = "px",  bg = "white")
#split.screen(c(2,2))
#screen(1) ; hist(cor(t(aned_low[,-c(1:3)])), main="Low signals begining", freq=F)
#screen(2) ; hist(cor(t(aned_high[,-c(1:3)])), main="High signals begining", freq=F)

#Remove batch effect
source(paste0(GenScriptPath,"ComBat_script.R"), local=T)

#Remove probesets below noise threshold
ProbeSum <- apply(aned_good[,-c(1:3)], 1, function(x) quantile(x, 0.95) > MaxNoise)
aned_high <- aned_good[ProbeSum,]
aned_low <- aned_good[!ProbeSum,]

print("Getting the heatmaps")
#Create Z-score heatmaps for low and high signals
samples <- 4:ncol(aned_good)
Z_scores_LOW <- MOD_z_scores_ALL_genes(samples, aned_low , paste0(resultsPath, study,"_Low expression all samples"))
Z_scores_High <- MOD_z_scores_ALL_genes(samples, aned_high , paste0(resultsPath, study, "_High expression all samples"))


#Exclude outliers
if (length(labels(Z_scores_LOW$colDendrogram[[1]])) < 0.2*length(labels(Z_scores_LOW$colDendrogram))) {
  exclude_samples_low <- labels(Z_scores_LOW$colDendrogram[[1]])
} else if (length(labels(Z_scores_LOW$colDendrogram[[2]])) < 0.2*length(labels(Z_scores_LOW$colDendrogram))) {
  exclude_samples_low <- labels(Z_scores_LOW$colDendrogram[[2]])
} else {
  exclude_samples_low <- NULL
}

if (length(labels(Z_scores_High$colDendrogram[[1]])) < 0.2*length(labels(Z_scores_High$colDendrogram))) {
  exclude_samples_high <- labels(Z_scores_High$colDendrogram[[1]])
} else if (length(labels(Z_scores_High$colDendrogram[[2]])) < 0.2*length(labels(Z_scores_High$colDendrogram))) {
  exclude_samples_high <- labels(Z_scores_High$colDendrogram[[2]])
} else {
  exclude_samples_high <- NULL
}

#Criteria for sample exclusion: if both high and low genes have outliers - use the low ones (noise). 
if (!is.null(exclude_samples_low) & !is.null(exclude_samples_high)) {
  exclude_samples <- exclude_samples_low
} else if (is.null(exclude_samples_low) | is.null(exclude_samples_high)){
  exclude_samples <- union(exclude_samples_high, exclude_samples_low)
} else {
  exclude_samples <- "none"
}

if(ncol(aned_good) > 23){
  aned_good_org <- aned_good
  aned_good <- aned_good[,!names(aned_good) %in% exclude_samples]
  
  Metadata_org <- Metadata
  Metadata <- Metadata[!Metadata$CommonName %in% exclude_samples,]
} else {
  if(exclude_samples != "none"){
    warning(paste0(paste0(exclude_samples, collapse = ","), "are suspected outliers. Less than 20 samples, no outliers removed")
)  }
}

#Define genes above and below noise threshold
ProbeSum <- apply(aned_good[,-c(1:3)], 1, function(x) quantile(x, 0.95) > MaxNoise)
aned_high <- aned_good[ProbeSum,]
aned_low <- aned_good[!ProbeSum,]

#Multiple-probset treatment
aned_high <- aned_high[names(rev(sort(apply(aned_high[,-c(1:3)], 1, sd)))),]
aned_high <- aned_high[!duplicated(aned_high[,2]),]

#print("Working on the last histograms")
#screen(3) ; hist(cor(t(aned_low[,-c(1:3)])), main="Low signals end", freq=F)
#screen(4) ; hist(cor(t(aned_high[,-c(1:3)])), main="High signals end", freq=F)

#close.screen(all.screens=T)
closeDev()
