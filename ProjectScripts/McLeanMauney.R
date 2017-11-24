McLeanMauney <- read.table("McLeanCortex/McLean_Mauney2013.csv", sep = "\t", header=TRUE, na.strings = "NA")
McLeanMauney$Profile <- sapply(McLeanMauney$Profile, function(x) if(x == "Control"){
  "Cont"
  } else if(x == "Bipolar"){
    "BP"
    } else {
      "SCZ"
      }) %>% factor()

McLeanMauney$Profile <- relevel(McLeanMauney$Profile, ref = "Cont")

McLeanMetaFiltered <- MetadataMcLean %>% filter(CommonName %in% McLeanMauney$CommonName)
p <- ggplot(McLeanMetaFiltered, aes(Profile, GabaPV_Genes))
p + 
  geom_violin() +
  geom_boxplot(width = 0.2) +
  geom_jitter(width = 0.2)

temp <- PlotPCggplot(data = McLeanMetaFiltered, CellName =  "GabaPV", name = "McLeanMetaFiltered")
ggsave("McLeanCortex/Cortex/GabaPV_filtered.pdf",
       plot = temp, width = 12, height = 8, units = "in",
       dpi=300)
temp <- PlotPCggplot(data = McLeanMetaFiltered, CellName =  "Astrocyte", name = "McLeanMetaFiltered")
ggsave("McLeanCortex/Cortex/Astrocyte_filtered.pdf",
       plot = temp, width = 12, height = 8, units = "in",
       dpi=300)


MetaBitanihirwe <- read.table("McLeanBitanihirwe.csv", sep = "\t", header = TRUE)
MetaBitanihirwe$Psychosis <- sapply(MetaBitanihirwe$Diagnosis, function(x) {
  if(grepl("b$", as.character(x))){
    "Yes"
  } else {
    "No"
  }
}) %>% factor

MetaBitanihirwe$Case <- paste0("Case_", MetaBitanihirwe$Case) %>% factor
MetaBitanihirwe$Diagnosis <- sapply(as.character(MetaBitanihirwe$Diagnosis), function(x) gsub("b$", "", x)) %>% factor
levels(MetaBitanihirwe$Diagnosis) <- c("BP", "Cont")
MetaBitanihirwe$pH <- sapply(MetaBitanihirwe$pH, function(x) round(as.numeric(as.character(x)), digits = 1))
MetaBitanihirwe$CharVec <- apply(MetaBitanihirwe %>% select(Age, Sex, PMI, pH), 1, function(x) paste(x, collapse = "_"))

# First match on based on common metadata (exact values for Age, Sex, PMI and pH)
SampleMatch <- list()
SampleMatch$NoMatch <- vector()
SampleMatch$MatchMcLean <- vector()
SampleMatch$MatchBitanihirwe <- vector()

for(meta in unique(MetadataMcLean$CharVec)){
  MatchBitanihirwe <- MetaBitanihirwe %>% filter(CharVec %in% meta) %>% .$Case %>% as.character()
    if(length(MatchBitanihirwe) == 0){
    SampleMatch$NoMatch <- c(SampleMatch$NoMatch, MetadataMcLean %>% filter(CharVec == meta) %>% .$sample %>% as.character)
  } else if(length(McLeanMatch) == 1){
    SampleMatch$MatchMcLean <- c(SampleMatch$MatchMcLean, MetadataMcLean %>% filter(CharVec == meta) %>% .$sample %>% as.character)
    SampleMatch$MatchBitanihirwe <- c(SampleMatch$MatchBitanihirwe, MatchBitanihirwe)
  }
}

# For the non matched samples - match based on common metadata with rounded PMI and pH
subMetaMcLean <- MetadataMcLean %>% filter(!sample %in% SampleMatch$MatchMcLean)
subMetaMcLean$CharVec <- paste(subMetaMcLean$Profile, subMetaMcLean$Sex, round(subMetaMcLean$PMI, digits=0),
                               round(subMetaMcLean$pH, digits=1), sep="_")
subMetaBitanihirwe <- MetaBitanihirwe %>% filter(!Case %in% SampleMatch$MatchBitanihirwe)
subMetaBitanihirwe$CharVec <- paste(subMetaBitanihirwe$Diagnosis, subMetaBitanihirwe$Sex, round(subMetaBitanihirwe$PMI, digits=0),
                                    round(subMetaBitanihirwe$pH, digits=1), sep="_")
subSampleMatch <- list()
subSampleMatch$NoMatch <- vector()
subSampleMatch$MatchMcLean <- vector()
subSampleMatch$MatchBitanihirwe <- vector()

for(meta in unique(subMetaMcLean$CharVec)){
  MatchBitanihirwe <- subMetaBitanihirwe %>% filter(CharVec %in% meta) %>% .$Case %>% as.character()
  if(length(MatchBitanihirwe) == 0){
    subSampleMatch$NoMatch <- c(subSampleMatch$NoMatch, subMetaMcLean %>% filter(CharVec == meta) %>% .$sample %>% as.character)
  } else if(length(MatchBitanihirwe) == 1){
    subSampleMatch$MatchMcLean <- c(subSampleMatch$MatchMcLean, subMetaMcLean %>% filter(CharVec == meta) %>% .$sample %>% as.character)
    subSampleMatch$MatchBitanihirwe <- c(subSampleMatch$MatchBitanihirwe, MatchBitanihirwe)
  } 
}


# Look up the non matched samples in the Mauny metadata
McLeanCommon <- subMetaMcLean %>% filter(!sample %in% subSampleMatch$MatchMcLean, Profile != "SCZ") %>% select(CommonName) %>% unlist
McLeanMauney %>% filter(CommonName %in% McLeanCommon) %>% select(Case1, Profile, Sex, PMI, CauseOfDeath)
subMetaBitanihirwe %>% filter(!Case %in% subSampleMatch$MatchBitanihirwe) %>% arrange(CharVec) %>% select(Case, Diagnosis, Sex, PMI, Cause.of.death)

manualSampleMatch <- list()
manualSampleMatch$NoMatch <- c("1", "2", "5", "6", "11", "14", "19", "33", "34", "35", "37", "40", "41")
manualSampleMatch$MatchMcLean <- c("mclean66.1063_01")
manualSampleMatch$MatchBitanihirwe <- c("Case_29")

MatchedSamples <- data.frame(McLean  = c(SampleMatch$MatchMcLean, subSampleMatch$MatchMcLean, manualSampleMatch$MatchMcLean),
                             Bitanihirwe = c(SampleMatch$MatchBitanihirwe, subSampleMatch$MatchBitanihirwe, manualSampleMatch$MatchBitanihirwe))

p <- ggplot(MetadataMcLean %>% filter(sample %in% MatchedSamples$McLean), aes(Profile, GabaPV_Genes))
p + 
  geom_violin() +
  geom_boxplot(width = 0.2) +
  geom_jitter(width = 0.2)

p <- ggplot(MetaBitanihirwe %>% filter(Case %in% MatchedSamples$Case), aes(Profile, GabaPV_Genes))
p + 
  geom_violin() +
  geom_boxplot(width = 0.2) +
  geom_jitter(width = 0.2)