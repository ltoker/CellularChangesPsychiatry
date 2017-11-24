packageF("pROC")
packageF("cowplot")

GetMeta <- function(StudyName, subCol, groups=TRUE){
  if(grepl("Stanley", StudyName)){
    collection = gsub("Stanley", "", StudyName)
    StanleyStudies <- read.table("StanleyStudies.txt", sep="\t", header=T)
    MetaAll <- read.table(paste0(path, "MetadataStanley.csv"), header = T, sep="\t")
    names(MetaAll) <- sapply(names(MetaAll), function(x) gsub("\\.*", "", x))
    Metadata <- subset(MetaAll, MetaAll$CollectionType == collection)
    
    for(meta in list.files(path = StudyName, pattern = "Study.*?csv")){
      assign(paste0("Meta", sapply(meta, function(x) strsplit(x, "\\.")[[1]][1])),
             read.table(paste0(StudyName, "/", meta), header = T, sep="\t"))
      
    }
    #######################################################################################################
    
    #Create data_frame matching the subject ids from different studies
    samples <- rep(NA, length = nrow(Metadata))
    SampleIDs <- as.data.frame(sapply(ls(pat="MetaStudy"),
                                      function(x) {
                                        temp <- eval(as.name(x)) 
                                        names(temp) <- gsub("\\.", "", names(temp))
                                        samples[match(temp$StanleyID, Metadata$StanleyID)] <- as.character(temp$Filename)
                                        return(samples)
                                      }, simplify = "array" ))
    rm(samples)
    rownames(SampleIDs) <- Metadata$StanleyID
    names(SampleIDs) <- gsub("Meta", "", names(SampleIDs))
    Metadata <- cbind(SampleIDs, Metadata)
    Metadata <- Metadata[order(Metadata$Profile),]
    Metadata %<>% melt(measure.vars = grep("Study", names(.)), variable.name="Study", value.name="Filename", na.rm=T)
    Metadata$Filename <- sapply(Metadata$Filename, function(filename) list.files(paste0(name, "/data/"),
                                                                                 pattern=paste0("^",filename)))
    Metadata$Filename <- sapply(Metadata$Filename, function(x) gsub(".cel", "", x, ignore.case = T))
    Metadata$Series_sample_id <- Metadata$Filename
    Metadata$Series_sample_id <- sapply(Metadata$Series_sample_id, function(x) strsplit(x, "\\.|_")[[1]][1])
    names(Metadata)[grep("pmi", names(Metadata), ignore.case=TRUE)] <- "PMI"
    names(Metadata)[grep("gender", names(Metadata), ignore.case=TRUE)] <- "Sex"
    names(Metadata)[grep("ph", names(Metadata), ignore.case=TRUE)] <- "pH"
    
    #Set up the region
    Metadata$Region <- StanleyStudies$Region[match(sapply(Metadata$Study,
                                                          function(x) gsub("Study.[0-9]?", "", x)),
                                                   StanleyStudies$Investigator)]
    #Set up the platform
    Metadata$Platform <- StanleyStudies$ArrayType[match(sapply(Metadata$Study,
                                                               function(x) gsub("Study.[0-9]?", "", x)),
                                                        StanleyStudies$Investigator)]
    
    Metadata$Platform <- sapply(Metadata$Platform, function(x) {
      if(grepl("hgu133a", tolower(x))){
        "GPL570"
        } else if (grepl("hgu133p", tolower(x))){
          "GPL96"
          } else if (grepl("hgu95av2", tolower(x))){
            "GPL8300"
            } else {
              NA
              }
      }, simplify=TRUE)
    
    } else {
      if(length(list.files(name, pattern="Metadata.tsv")) == 0){
        if(length(list.files(path = paste0(StudyName,"/data"), pattern=".soft")) == 0){
          softDown(StudyName, paste0(StudyName, "/data/",paste0(StudyName, ".soft")))
        }
      Metadata <- ReadSoft(paste0(name,"/data/", list.files(paste0(name,"/data/"), pattern=".soft")))
      if(name == "GSE1138"){
        names(Metadata)[grepl("stage.of.illness", names(Metadata))] <- "Stage_of_illness"
        Metadata$Profile = sapply(Metadata$Stage_of_illness, function(x) strsplit(as.character(x), " ?- ?")[[1]][2])
        Metadata$Stage_of_illness <- sapply(Metadata$Stage_of_illness, function(x) strsplit(as.character(x), " ?- ?")[[1]][1]) %>% factor
      }
      write.table(Metadata, paste0(name,"/Metadata.tsv"), sep="\t", row.names=FALSE)
      }
      
    Metadata <- read.table(paste0(name,"/Metadata.tsv"), sep="\t", header=T)
    names(Metadata)[grep("pmi|mortem", names(Metadata), ignore.case=TRUE)] <- "PMI"
    names(Metadata)[grep("ph$|tissue.ph", names(Metadata), ignore.case=TRUE)] <- "pH"
    names(Metadata)[grep("^age$|age.*?death", names(Metadata), ignore.case=TRUE)] <- "Age"
    names(Metadata)[grep("disease|Diagnosis|profile|phenotype", names(Metadata), ignore.case=TRUE)] <- "Profile"
    names(Metadata)[grep("gender|sex", names(Metadata), ignore.case=TRUE)] <- "Sex"
    names(Metadata)[grep("platform", names(Metadata), ignore.case = TRUE)] <- "Platform" 
    Metadata$Filename = Metadata$Series_sample_id
    names(Metadata) <- sapply(names(Metadata), function(x) gsub(" ", "_", x))
    }
  if(groups & sum(grepl("disease|diagnosis|phenotype|profile", tolower(names(Metadata)))) == 0){
    print("No group information, modify Metadata file")
  } else if (groups & sum(grepl("disease|diagnosis|profile", tolower(names(Metadata)))) > 0){
    # Harmonize profile names
    Metadata$Profile <- sapply(Metadata$Profile, function(x) {
      if(grepl("bipolar", tolower(x))){
        "BP"
      } else if (grepl("control|unaffected|normal", tolower(x))){
        "Cont"
      } else if (grepl("schizo|scz", tolower(x))){
        "SCZ"
      } else if (grepl("depres", tolower(x))){
        "MD"
      } else {
        NA
      }
    }, simplify=TRUE) %>% factor
    
    
    #remove MD subjects
    Metadata %<>% filter(!is.na(Profile),  Profile != "MD") %>% droplevels()
    Metadata <- Metadata[order(Metadata$Profile),]
    
  } else {
    Metadata$Profile <- "Cont"
  }
  
  Metadata <- GetCommonName(Metadata,char=c("age", "sex", "PMI", "ph"))
  OrgRegion <- grep("(tissue|region|brain.?region)", names(Metadata),ignore.case = TRUE, value = T)
  if(length(OrgRegion) == 0){
    print("No brain region specified, modify Metadata file")
  }
  
  if(length(OrgRegion) > 1) {
    print("Multiple columns define region, modify Metadata file")
    OrgRegion <- grep("(tissue$|region|brain.?region)", names(Metadata),ignore.case = TRUE, value = T)
    if(length(OrgRegion == 0)){
      print("No columns matches tissue description")
    } else if(length(OrgRegion > 1)){
      print("more than one column matched tissue description, check out the Metadata file")
    }
  }
  Metadata$NeuExpRegion = Metadata %>%
    select_(.dots = OrgRegion) %>% unlist
  
  Metadata$NeuExpRegion <- sapply(Metadata$NeuExpRegion, function(x) {
    if (grepl("cerebe",tolower(x))){
      "Cerebellum"
    } else if (grepl("cortex|pfc|dlpfc|frontalba|^ba|gyrus",tolower(x))){
      "Cortex"
    } else if (grepl("hippocampus|hip|hpc",tolower(x))){
      "Hippocampus"
    } else if (grepl("thalamus",tolower(x))) {
      "Thalamus"
    } else if (grepl("str",tolower(x))){
      "Striatum"
    } else if (grepl("putamen",tolower(x))){
      "Putamen"
    } else if (grepl("nigra",tolower(x))){
      "SubstantiaNigra"
    } else if (grepl("brain",tolower(x))){
      "Cortex"
    } else {
      NA
    }
  }, simplify=TRUE) %>% factor
  
  if(Metadata[[OrgRegion]] %>% unique %in% c("brain", "Brain")){
    warning("No specific brain region specified, rigion was set to Cortex")
  }
  
  if(length(Metadata[[OrgRegion]] %>% unique) > length(Metadata$NeuExpRegion %>% unique)){
    subCol = OrgRegion
  }
  #If there are replicates, change their names accordngly
  repSamples <- sapply(unique(Metadata[[subCol]]), function(subdata){
    subMeta <- Metadata %>% filter_(.dots = paste0(subCol, " =='", subdata, "'"))
    RepSample <- subMeta$CommonName[!duplicated(subMeta$CommonName)]
    DF <- subMeta %>% filter(CommonName %in% RepSample) %>% select(Series_sample_id, CommonName)
    RepName <- sapply(as.character(unique(DF$CommonName)), function(SampName){
      names <- subMeta %>% filter(CommonName == SampName) %>% select(Series_sample_id, CommonName)
      cbind(as.character(names$Series_sample_id), paste0(unique(names$CommonName), "rep",seq(1:nrow(names))))
    }, simplify=TRUE) %>% do.call(rbind, .) %>% data.frame()
    DF$RepName <- RepName$X2[match(DF$Series_sample_id, RepName$X1)]
    DF
  }, simplify=FALSE)
  names(repSamples) <- unique(Metadata[[subCol]])
  Metadata$CommonName <- as.character(Metadata$CommonName)
  for(region in names(repSamples)){
    Metadata$CommonName[match(repSamples[[region]]$Series_sample_id,
                              Metadata$Series_sample_id)] <- as.character(repSamples[[region]]$RepName) 
  }
  Metadata$CommonName <- as.factor(Metadata$CommonName)
  Metadata[[name]] <- Metadata$Series_sample_id
  return(Metadata)
}


ReadCell <- function(path, QC = 0, platform, CelFiles=list.files(paste0(name, "/data"))){ # QC =1 / 0, QC = 1 runs arrayQualityMetrics. Takes some time.. 
  Anno_file <- GetAnnoFiles(platform)
  if(platform %in% c("GPL96", "GPL97", "GPL570", "GPL1352", "GPL1261", "GPL390")){
    source(paste0(GenScriptPath,"Read_and_annotate_CELaffy.R"), local = T)
    
    #Add Gene symbol and annotation
    aned <- add_symbols(ned, Anno_file)
    names(aned) <- sapply(names(aned), function(x) strsplit(x, "\\.")[[1]][1])
    #Identify scanDate
    scanDate <- protocolData(affydata)$ScanDate
    scanDate <- sapply(scanDate, function(x) strsplit(x, " |T")[[1]][1])
    names(scanDate) <- rownames(pData(affydata))
    names(scanDate) <- sapply(names(scanDate), function(x) strsplit(x, "\\.")[[1]][1])
  
    } else if (platform %in% c("GPL6244", "GPL5175", "GPL5188")) {
    source(paste0(GenScriptPath,"Read_and_annotate_CELoligo.R"), local = T)
    
    #Add Gene symbol and annotation
    aned <- add_symbols(exp_value, Anno_file)
    names(aned) <- sapply(names(aned), function(x) strsplit(x, "\\.|_")[[1]][1])
    
    #Identify scanDate
    scanDate <- pData(protocolData(data))$dates %>% as.character
    scanDate <- sapply(scanDate, function(x) strsplit(x, "T| ")[[1]][1])
    names(scanDate) <- rownames(pData(protocolData(data)))
    names(scanDate) <- sapply(names(scanDate), function(x) strsplit(x, "\\.|_")[[1]][1])
  }
  
  aned[,-c(1:3)] <- apply(aned[,-c(1:3)], 2,
                          function(x) as.numeric(as.character(x))) 
  study <-list("aned" = aned, "scanDate" = scanDate)
  return(study)
}

makeSet <- function(data, meta, name=NULL){
  expData <- data$aned
  #Add probset names as rownames
  if(is.null(rownames(expData))){
    rownames(expData) <- expData$Probe
  }
  
  #Add the scan date to the meta file
  meta$ScanDate <- data$scanDate[pmatch(meta$Series_sample_id, names(data$scanDate))]
  
  studyName <- name
  
  if(is.null(name)){
    studyName <- readline("What is study name?")  
  }
  index <- grep(studyName, names(meta), ignore.case=T) #get the relevan colomn for study sample names
  meta <- meta[!is.na(meta[index]),]
  names(meta) <- sapply(names(meta), function(x) strsplit(x, "\\.\\.")[[1]][1])
  rownames(meta) <- as.character(meta$CommonName) 
  expr <- as.matrix(expData[,sapply(expData[1,], function(x) is.numeric(x))])
  phenoData <- AnnotatedDataFrame(meta[match(colnames(expr), rownames(meta)),])
  featureData <-  AnnotatedDataFrame(expData[,c(1:3)])
  eSetObj <- ExpressionSet(assayData=expr,
                           phenoData = phenoData, featureData=featureData ) 
  return(eSetObj)                        
}

PreProcces <- function(eSet, study=NULL){
  aned <- as.data.frame(cbind(pData(featureData(eSet)),exprs(eSet)))
  Metadata <- pData(eSet)
  source(paste0(GenScriptPath,"pre-proccess.R"), local=T)
  output <- list(aned_high, aned_good, aned_low, MaxNoise,Z_scores_High,
                 exclude_samples_low, exclude_samples_high, exclude_samples, Metadata)
  names(output) <- c("aned_high", "aned_good", "aned_low", "NoiseThreshold", "Z_scores_High",
                     "exclude_samples_low", "exclude_samples_high", "exclude_samples", "Metadata")
  return(output)
}

datasGenerate <- function(genes){
  datas <- lapply(sapply(names(studyFinal), function(s) studyFinal[[s]]$aned_good,
                         simplify=FALSE),
                  function(x) subset(x, GeneSymbol %in% genes & Probe != "243712_at"))
    

  names(datas) <- names(studyFinal)
  return(datas)
}

HeatMapGen <- function(datas, Meta, path=NULL, save = 1){
  out <- sapply(names(datas), function(x) {
    if(nrow(datas[[x]]) > 1){
      study <- x
      x <- datas[[x]]
      dummy <- data.frame("CommonName"=names(x)[sapply(x[1,], function(y) is.numeric(y))])
      KmeanGen <- GeneSex(x, dummy )
      KmeanGen$BioGender[KmeanGen$BioGender=="F"] <- "grey"; KmeanGen$BioGender[KmeanGen$BioGender=="M"] <- "black"
      KmeanSex <- KmeanGen$BioGender
      data <- t(x[sapply(x[1,], function(y) is.numeric(y))]) 
      values <- data.frame(row.names=Meta$CommonName)
      values[,1:ncol(data)] <- NA 
      values[match(rownames(data), rownames(values)),] <- data
      MetaSex = as.character(Meta$Sex[match(rownames(data), Meta$CommonName)])
      MetaSex[grepl("F|female|^wom", MetaSex, ignore.case=T)] <- "deeppink"
      MetaSex[grepl("^M|^male|^man", MetaSex, ignore.case=T)] <- "darkslateblue"
      Row_col = as.character(x$GeneSymbol)
      AllSexGen <- data.frame(gene = unique(Row_col))
      AllSexGen$sex <- sapply(unique(Row_col), function(gene){
        if(grepl("KDM5D|RPS4Y1", gene)){
          "(Male)"
        } else if (grepl("XIST", gene)){
          "(Female)"
        } else{
          NA
        }
      })
      AllSexGen$color <- sapply(AllSexGen$gene, function(x){
        switch(as.character(x),
               "XIST" = "darkred",
               "KDM5D" = "blue",
               "RPS4Y1|RPS4Y2" = "darkblue")
      })
      Row_col <- sapply(Row_col, function(x){
        switch(x,
               "XIST" = "darkred",
               "KDM5D" = "blue",
               "RPS4Y1|RPS4Y2" = "darkblue")
      })
      
      Row_col <- t(as.matrix(Row_col))
      Col_col <- cbind(MetaSex, KmeanSex)
      colnames(Col_col) <- c("Metadata Gender", "Kmeans cluster")
      rownames(Row_col) <- "Gene"
      myPalette <- colorRampPalette(c("skyblue1", "#123766"))(99)
      if(save==1){
        pdf(paste0(path,study, "Gender.pdf"),width=20, height=15, useDingbats = FALSE)
        cex1 = 2
        cex2 = 2
        cexCol = 2
        cex.axis=3
      } else if (save==2) {
        cex1  = 0.8
        cex2 = 1
        cexCol=1
        cex.axis=1.5
      }
      heatmap.2b(t(scale(data)), 
                 dendrogram="none",Rowv="none", density.info="density",
                 ColSideColorsSize = 2, cexCol = cexCol,
                 col=myPalette, RowSideColors=Row_col,
                 ColSideColors = Col_col,
                 symbreaks=T, trace="none", symkey=T,
                 margins = c(10, 10), na.color="black",  
                 key=T,  keysize = 0.8,  key.title = "", KeyValueName="Normalized expression", 
                 na.rm=T, cexRow = 1.5 , cex.axis=cex.axis)
      legend("left",  y.intersp=2, x.intersp = 0.2, cex=cex1, bty="n",
             legend=apply(AllSexGen, 1,  function(x){
               paste(x[1:2], collapse="\n")
             }),
             fill=AllSexGen$color)
      legend("top", xjust=0.5,legend=c("Male\n(Metadata)", "Female\n(Metadata)", "","Male\n(cluster)", "Female\n(cluster)"),
             fill=c("darkslateblue", "deeppink", "white", "black", "grey"), bty="n",
             border=FALSE, x.intersp=0.8, cex=cex2, horiz=T)
      if (save==1){
        dev.off()
      }
    } else {
      print("Only one sex gene, better use a plot")
    }
  })
  names(out) <- names(datas)
}


HeatMapGen2 <- function(datas, Meta, missmatched, path=NULL, save = 1){
  if(length(names(datas)) >1){
    Combined <- lapply(datas, function(x) {
      data <- t(x[sapply(x[1,], function(y) is.numeric(y))]) 
      values <- data.frame(row.names=Meta$CommonName)
      values[,1:ncol(data)] <- NA 
      values[match(rownames(data), rownames(values)),] <- data
      colnames(values) <- colnames(data)
      values
    })
    Features <- lapply(datas, function(x) {
      Probes <- x[sapply(x[1,], function(y) !is.numeric(y))]
      Probes
    })
    

    AllData <- do.call(cbind, args=Combined)
    
    #Add region/study to colnames in case it is not there already
    if(!grepl("\\.", names(AllData)[1])){
      names(AllData) <- paste(names(datas), names(AllData), sep=".")
    }
    
    
    #Remove samples that were excluded in all of the datasets
    SampleRM <- apply(AllData, 1, complete.cases) %>% t %>% rowSums
    SampleRM <- names(SampleRM)[SampleRM == 0]
    
    AllData <- subset(AllData, !rownames(AllData) %in% SampleRM)
    Meta <- subset(Meta, !CommonName %in% SampleRM)
    
    AllData <- AllData[order(Meta$Sex),]
    AllFeatures <- unique(do.call(rbind, args=Features))
    MetaSex = as.character(Meta$Sex[match(rownames(AllData), Meta$CommonName)])
    MetaSex[MetaSex == "F"] <- "deeppink" ; MetaSex[MetaSex == "M"] <- "darkslateblue"
    Probe_col <-  sapply(names(AllData), function(x) strsplit(x, "\\.")[[1]][2])
    Probe_col <- as.character(AllFeatures$GeneSymbol[match(Probe_col, AllFeatures$Probe)])
    AllSexGen <- data.frame(gene = unique(Probe_col))
    AllSexGen$sex <- sapply(unique(Probe_col), function(gene){
      if(grepl("KDM5D|RPS4Y1", gene)){
        "(Male)"
      } else if (grepl("XIST", gene)){
        "(Female)"
      } else{
        NA
      }
    })
    AllSexGen$color <- sapply(AllSexGen$gene, function(x){
      switch(as.character(x),
             "XIST" = "darkred",
             "KDM5D" = "blue",
             "RPS4Y1|RPS4Y2" = "darkblue")
    })
    
    Probe_col <- sapply(Probe_col, function(x){
      switch(x,
             "XIST" = "darkred",
             "KDM5D" = "blue",
             "RPS4Y1|RPS4Y2" = "darkblue")
    })
    
    Study <- sapply(colnames(AllData), function(x)
      strsplit(x, "\\.")[[1]][1])
    Study_col <- c("palegreen4", "black", "palevioletred4", "steelblue", "aquamarine3")[as.factor(Study)]
    Row_col <- rbind(Study_col, Probe_col)
    Col_col <- as.matrix(MetaSex)
    names(AllData) <- sapply(names(AllData), function(x) strsplit(x, "\\.")[[1]][2])
    colnames(Col_col) <- "Metadata Gender"
    rownames(Row_col) <- c("Study", "Gene")
    myPalette <- colorRampPalette(c("skyblue1", "#123766"))(99)
    SampleCol <- rep("grey", nrow(AllData))
    StudyNum <- unique(tolower(sapply(Study, function(x) strsplit(x, "\\.")[[1]][1])))
    sapply(StudyNum, function(x){
      MM <- sapply(missmatched[[grep(x, names(missmatched), ignore.case=T)]],
                   function(subject) grep(paste0(subject,"$"), rownames(AllData)), simplify=T)
      if(length(MM)>0){
        studies <- unique(tolower(Study))
        SampleCol[MM] <<- unique(Study_col)[grep(x, studies)] 
      }
    })
    
    if(save==1){
      pdf(paste0(path,"CombinedGenderHeatmap.pdf"),width=20, height=15, useDingbats = FALSE)
      cex1 = 2
      cex2 = 2
      cexCol = 2
      cex.axis=3
    } else if (save==2) {
      cex1  = 0.8
      cex2 = 1
      cexCol=1
      cex.axis=1.5
    }
    heatmap.2b(t(scale(AllData)), density.info="density",
               dendrogram="none",Rowv="none", Colv="none",
               colCol=SampleCol, colRow = "black", 
               RowSideColorsSize = 1.5,
               cexCol=cexCol, cexRow=2.5,
               col=myPalette, RowSideColors=Row_col,
               ColSideColors = Col_col,
               symbreaks=T, trace="none", 
               margins = c(12, 12), na.color="grey80",  
               key=T,  keysize = 0.8,  key.title = "", KeyValueName="Normalized expression", 
               na.rm=F, cex.axis=cex.axis)
    legend("left",  y.intersp=2, x.intersp = 0.2, cex=cex1, bty="n",
           legend=apply(AllSexGen, 1,  function(x){
             paste(x[1:2], collapse="\n")
           }),
           fill=AllSexGen$color)
    legend("bottomleft",  y.intersp=2, x.intersp = 0.2, cex=cex1, bty="n", 
           legend=unique(sapply(Study,
                                function(x) {
                                  name <- strsplit(x, "\\.")[[1]][1]
                                  name <- strsplit(name, "_")[[1]]
                                  if(length(name) > 1){
                                    paste0(name[1],"\n",name[2])
                                  } else{
                                    name
                                  }
                                })),
           fill=unique(Study_col))
    legend("top", xjust=0.5, legend=c("Female", "Male"),
           fill=c("deeppink", "darkslateblue"), bty="n",
           border=FALSE, x.intersp=0.8, cex=cex2, horiz=T)
    
    
    if (save==1){
      dev.off()
    }
  } else {
    print("Only one dataset")
  }
}

PlotAllStudies <- function(Cells=CellType_genes, data, Meta=MetaConsort, remove = "Marker", main=NULL) { 
  Cells[[remove]] <- NULL #remove the non-relevant list elements
  Data <- lapply(data, function(x) x$modified)
  PC_combined <- list()
  for(i in names(Cells)){
    print(i)
    PC1 <- lapply(Data, function(dat){
      pc_1 <- data.frame(row.names=Meta$CommonName)
      pc_1[,1] <- NA
      pc_1[match(rownames(dat[[i]]$x), as.character(Meta$CommonName)),1] <- dat[[i]]$x[,1]
      pc_1
    }) 
                                                 
    all_PCA <- as.data.frame(PC1)
    colnames(all_PCA) <- c("study1", "study3", "study5", "study7")
    exclude <- apply(all_PCA, 1, function(x) sum(is.na(x)))
    all_PCA <- all_PCA[!exclude %in% c(3,4),]
    if(nrow(all_PCA) > 0){
      all_PCA <- all_PCA[order(apply(all_PCA,1, function(x) median(x, na.rm=T))),]
      stripchart(as.data.frame(t(all_PCA)),
                 main=paste(strsplit(i, "_Genes")[[1]][1], main),
                 vertical=T,
                 col="white",
                 xaxt="none",
                 las=2)
      boxplot(t(all_PCA[grep("BP", rownames(all_PCA)),]),
              col="brown4",
              at=grep("BP", rownames(all_PCA)),
              add=T, xaxt="none",
              yaxt="none")
      boxplot(t(all_PCA[grep("Cont", rownames(all_PCA)),]),
              col="burlywood",
              add=T,
              at=grep("Cont", rownames(all_PCA)),
              xaxt="none",
              yaxt="none")
      boxplot(t(all_PCA[grep("SCZ", rownames(all_PCA)),]),
              col="chartreuse4",
              at=grep("SCZ", rownames(all_PCA)),
              xaxt="none",
              yaxt="none",
              add=T)
      stripchart(as.data.frame(t(all_PCA)),
                 cex=0.7,
                 pch=16,
                 vertical=T,
                 xaxt="none",
                 yaxt="none",
                 add=T)
      axis(side=1, at=grep("SCZ", rownames(all_PCA)),
           labels=rownames(all_PCA)[grep("SCZ", rownames(all_PCA))],
           cex.axis=0.6, las=2, col.axis="chartreuse4")
      axis(side=1, at=grep("BP", rownames(all_PCA)),
           labels=rownames(all_PCA)[grep("BP", rownames(all_PCA))],
           cex.axis=0.6, las=2, col.axis="brown4")
      axis(side=1, at=which(rownames(all_PCA) %in% BP_II),
           labels=paste("***", rownames(all_PCA)[which(rownames(all_PCA) %in% BP_II)]),
           cex.axis=0.6, las=2, col.axis="brown4")
      axis(side=1, at=grep("Cont", rownames(all_PCA)),
           labels=rownames(all_PCA)[grep("Cont", rownames(all_PCA))],
           cex.axis=0.6, las=2, col.axis="burlywood")
      
      legend("bottomright", legend=c("Control","Bipolar", "Schizophrenia"), fill=c("burlywood", "brown4", "chartreuse4"), cex=1.5)
    
      PC_combined[[i]] <- all_PCA
    }
  }
  return(PC_combined)
}

FindRep <- function(Metadata, char){
  names(Metadata) <- tolower(names(Metadata))
  char <- tolower(char)
  char <- char[char %in% names(Metadata)]
  AllChar <- Metadata %>% select_(.dots=char)
  Metadata$CharVec <- apply(AllChar, 1, function(x) paste0(x, collapse="_"))
  RepTable <- table(Metadata$CharVec)
  Replicates <- Metadata %>% filter(CharVec %in% names(RepTable[RepTable > 1])) %>% arrange(CharVec)
  Unique <- Metadata %>% filter(!CharVec %in% names(RepTable[RepTable > 1]))
  return(list(Replic = Replicates,
              Unique = Unique))
}

GetCommonName <- function(Metadata, char=char){
  SameSub <-  FindRep(Metadata=Metadata, char=char)
  UniqueName <- rbind(SameSub$Replic[!duplicated(SameSub$Replic$CharVec),],
                      SameSub$Unique) %>% arrange(profile)
  UniqueName$CommonName <- UniqueName$profile
  sampleSum <- table(UniqueName$CommonName)
  UniqueName$CommonName <- melt(sapply(names(sampleSum), function(x) {
    paste(x, seq(1:sum(UniqueName$CommonName == x)), sep="_")
    
  }))$value
  Metadata$CharVec <- sapply(Metadata$Filename,
                             function(x) rbind(SameSub$Replic, SameSub$Unique) %>%
                               filter(filename == x) %>%
                               .$CharVec) %>% unlist
  Metadata$CommonName <- UniqueName$CommonName[match(Metadata$CharVec, UniqueName$CharVec)]
  return(Metadata)
}

GeneMGPcombined <- function(dataGenes, metaGenes,  dataMGP, NameVarMGP="CommonName", GeneVar = "GeneSymbol", GeneList){
  rownames(dataMGP) <- dataMGP[[NameVarMGP]] %>% as.character
  GeneExp <- dataGenes %>% filter_(paste0(GeneVar,  " %in% ", paste0("c(",paste0("'",GeneList,"'", collapse=","), ")")))
  rownames(GeneExp) <- dataGenes[[GeneVar]][dataGenes[[GeneVar]] %in% GeneList] %>% as.character
  rownames(GeneExp) <- paste0(rownames(GeneExp), "_gene")
  GeneExp <- t(GeneExp[sapply(names(GeneExp), function(x) is.numeric(GeneExp[[x]]))])
  dataCombined <- merge(dataMGP, GeneExp, by = "row.names")
  return(dataCombined)
} 

plotGeneMGPcor <- function(dataGenes, dataMGP, GeneList,
                           ListName = NULL,
                           grpVar = "Profile",CellVar=NULL, 
                           CellExtend = "_Genes", CellName = NULL,
                           grpRef = "Cont",
                           groups=c("Cont", "BP", "SCZ")){
  temp <- GeneMGPcombined(dataGenes = dataGenes, dataMGP = dataMGP , GeneList = GeneList)
  GeneMGPcor <- sapply(groups, function(grp){
    sapply(names(temp)[grepl("_gene", names(temp))], function(gene){
      cor.test(formula(paste0("~", CellVar, "+", gene)), data=temp %>% filter_(paste0(grpVar,"=='", grp,"'")))$estimate
    })
  }) %>% data.frame
  
  GeneMGPcor$GeneSymbol <- sapply(rownames(GeneMGPcor), function(x) strsplit(x, "_gene")[[1]][1])
  GeneMGPcor %<>% arrange(desc(.[[grpRef]]))
  GeneMGPcor <- melt(GeneMGPcor, id.vars="GeneSymbol", variable.name = grpVar, value.name="Cor")
  GeneMGPcor$GeneSymbol <- factor(GeneMGPcor$GeneSymbol, levels = unique(GeneMGPcor$GeneSymbol))
  
  grpColors = c("burlywood3", "cornflowerblue", "indianred4")
  ggplot(GeneMGPcor, aes(x = GeneSymbol, y=Cor))+
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_text(size = rel(0.8)),
          axis.text.x = element_text(size = rel(0.8), angle=90),
          panel.grid = element_blank()) +
    labs(title = CellName, x = ListName, y=paste0("Correlation to ",CellName, " MGP"))+
    scale_color_manual(values = grpColors, name="Group") +
    geom_point(aes_string(color=grpVar), size=3, shape=16)
}

getCerebellumAstro <- function(data){
  AstroGenes <- neuroExpressoAnalysis::mouseMarkerGenes$Cortex$Astrocyte
  AstroHuman <- mouse2human(AstroGenes)$humanGene
  BregmannGenes <- neuroExpressoAnalysis::mouseMarkerGenes$Cerebellum$Bergmann
  BregmannHuman <- mouse2human(BregmannGenes)$humanGene
  dataAstro <- data %>% filter(GeneSymbol %in% AstroHuman)
  dataAstro$Mean <- apply(dataAstro[sapply(names(data), function(x) is.numeric(data[[x]]))], 1, mean)
  dataAstro$Cell <- "Astrocyte"
  dataBregmann <- data %>% filter(GeneSymbol %in% BregmannHuman)
  dataBregmann$Mean <- apply(dataBregmann[sapply(names(data), function(x) is.numeric(data[[x]]))], 1, mean)
  dataBregmann$Cell <- "BregmannGlia"
  dataBoth <- rbind(dataAstro %>% select(GeneSymbol, Mean, Cell),
                    dataBregmann %>% select(GeneSymbol, Mean, Cell))
  return(dataBoth)
}

getROC <- function(MGPcorData){
  data <- sapply(names(MGPcorData), function(cellName){
    celltype = MGPcorData[[cellName]]
    cellName <- gsub("_Genes", "MGPcor", cellName)
    sapply(names(celltype), function(grpName){
      grp = celltype[[grpName]]
      grp$GeneType <-  sapply(grp$GeneSymbol, function(gene){
        if(gene %in% MistryGenesUP$GeneSymbol){
          "MistryUP"
        } else if(gene %in% MistryGenesDown$GeneSymbol){
          "MistryDown"
        } else {
          "NonMistry"
        }
      }, simplify = TRUE) %>% factor
      ROCup <- roc(GeneType~Cor, data = grp, ci=TRUE, levels = c("NonMistry", "MistryUP"),
                   plot=FALSE, direction="<",
                   legacy.axes=TRUE, print.auc=TRUE,
                   xlab="FPR", ylab = "TPR",
                   main = paste0(cellName, "(", grpName, "), ", "MistryUP"))
      ROCdown <- roc(GeneType~Cor, data = grp, ci=TRUE, levels = c("NonMistry", "MistryDown"),
                     plot=FALSE, direction="<",
                     legacy.axes=TRUE, print.auc=TRUE,
                     xlab="FPR", ylab = "TPR",
                     main = paste0(cellName, "(", grpName, "), ", "MistryDown"))
      list(ROCup = ROCup, ROCdown = ROCdown)
    }, simplify=FALSE)
  }, simplify=FALSE)
return(data)
}

plotROC <- function(rocObj, title, xlab="FPR", ylab = "TPR", auc = TRUE){
  
}

PlotAllStudyOneGeneMGPcor <- function(exclGRP = "MD", gene, MGPname = "GabaPV_Genes"){
  #get correlations
  corStat <- sapply(ls(pat="^Cortex", .GlobalEnv), function(std){
    study = eval(as.name(std))
    exp <- study$aned_high[,!grepl(exclGRP , names(study$aned_high))]
    meta <- study$Metadata %>% filter(Profile != exclGRP) %>% droplevels()
    meta <- meta[match(names(exp)[-c(1:3)], meta$CommonName),] 
    exp <- exp %>% filter(GeneSymbol == gene) %>% unlist %>% .[-c(1:3)] %>% as.numeric
    corStat = cor(exp, meta[[MGPname]], method="spearman", use="complete.obs") %>% round(digits=2)
    paste0("rho = ", corStat)
  }, simplify=FALSE)
  
  #combine all studies
  allStudyCor <- sapply(ls(pat="^Cortex", .GlobalEnv), function(std){
    study = eval(as.name(std))
    exp <- study$aned_high[,!grepl(exclGRP, names(study$aned_high))]
    meta <- study$Metadata %>% filter(Profile != exclGRP) %>% droplevels() %>%
      select_(.dots = c("CommonName", "Profile", MGPname))
    exp <- exp %>% filter(GeneSymbol == gene) %>% .[,-c(1:3)] %>% t
    colnames(exp) <- gene
    data <- merge(meta, exp, by.x="CommonName", by.y="row.names")
  }, simplify=FALSE) %>% do.call(rbind, .)
  
  allStudyCor$Study <- gsub("Cortex|\\..*", "", rownames(allStudyCor)) %>% as.factor
  
  #plot
  p <- ggplot(allStudyCor, aes_string(gene, MGPname, color = "Profile"))
  MGPname2 = gsub("_Genes", "", MGPname)
  txtSize = 12
  plot <- p + labs(x = paste0(gene, "expression(log2)"), y=paste0(MGPname2, "relative MGP"), fill="Profile")+
    theme_grey(base_size = txtSize) +
    theme(axis.text.y = element_text(size = rel(0.8)),
          axis.text.x = element_text(size = rel(1.2)),
          panel.grid = element_blank()) +
    scale_color_manual(values = c("burlywood3", "cornflowerblue", "indianred4")) +
    geom_point(pch=16) + 
    facet_wrap(~Study) +
    annotate("text", label = unlist(corStat), size = 0.25*txtSize, x = 6.5, y = 0.9)
  ggsave(paste0(gene, "_", MGPname2, "\\.pdf"), width = 8, height = 6, units = "in", plot = plot, dpi=300)
}

CreateMGPcellTypeDF <- function(cellData, cellMeta, MGPcorGenes, sampleRegex = "GSM",
                                cellTypeMetaVar = "PyramidalDeep",
                                idVars = c("Probe", "Gene.Symbol", "GeneNames"),title = NULL){
  dataMGPcorGenes <- sapply(names(MGPcorGenes), function(cellType){
    cellGenes <- na.omit(MGPcorGenes[[cellType]]) %>% as.vector
    data <- cellData %>% filter(Gene.Symbol %in% cellGenes)
    data
  }, simplify =  FALSE) %>% rbindlist(use.names = TRUE, idcol = "GeneType")
  
  #Add information regarding whether a gene is a marker
  dataMGPcorGenes$Marker <- NA
  for(i in 1:nrow(dataMGPcorGenes)){
    gene = dataMGPcorGenes$Gene.Symbol[i]
    mgp = dataMGPcorGenes[i,] %>% .$GeneType %>% gsub("_MGP", "", .)
    if(gene %in% markerGenes[[mgp]]){
      dataMGPcorGenes$Marker[i] <- "YES"
    } else {
      dataMGPcorGenes$Marker[i] <- "NO"
    }
  }
  
  #Normalize signals 0-1
  dataMGPcorGenesNorm <- apply(dataMGPcorGenes %>% select(matches(sampleRegex)), 1, function(gene){
    rescale(gene, c(0,1))
  }) %>% t
  
  dataMGPcorGenesNorm <- cbind(dataMGPcorGenes %>% select(-matches(sampleRegex)), dataMGPcorGenesNorm)
  
  dataMGPcorGenesNorm$Gene.Symbol  <- factor(dataMGPcorGenes$Gene.Symbol,
                                               levels = unique(dataMGPcorGenes$Gene.Symbol))
  dataMGPcorGenesMelt <- melt(dataMGPcorGenesNorm, id.vars = c("GeneType", idVars, "Marker"),
                                variable.name = "SampleName", value = "Exp")
  dataMGPcorGenesMelt$CellType <- factor(cellMeta[[cellTypeMetaVar]][match(dataMGPcorGenesMelt$SampleName, cellMeta$sampleNameCol)],
                                           levels = c("Astrocyte", "Microglia", "Oligo", "GabaPV", "GabaVIPReln",
                                                      "Pyramidal_S100a10", "PyramidalCorticoThalam"))
  dataMGPcorGenesMelt %<>% filter(!CellType %in% c("GabaRelnCalb", "Pyramidal_Glt_25d2", "Pyramidal_Thy1",
                                                     "Microglia_activation_MGP", "Microglia_deactivation_MGP")) %>% droplevels()
  dataMGPcorGenesMelt$GeneType <- factor(dataMGPcorGenesMelt$GeneType, levels = dataMGPcorGenesMelt$GeneType %>% unique())
  return(dataMGPcorGenesMelt)
}


PlotMGPcellTypes <- function(data, title=NULL, ylab="Normalized expression",
                             size=12, width=0.1, save = FALSE, path = NULL,
                             CellColors = c("goldenrod1", "grey", "darkolivegreen4",
                                                  "firebrick1", "firebrick", "dodgerblue3",
                                                  "deepskyblue2")){
  GeneNum <- group_by(data %>% filter(SampleName==levels(data$SampleName)[1]), GeneType) %>% summarise(n = n()) %>% data.frame
  data$GeneType2 <- sapply(data$GeneType, function(genetype){
    num = GeneNum %>% filter(GeneType == genetype) %>% .$n
    paste0("Top correlated genes - ", genetype, " (", num, ")")
  }) %>% factor(levels = unique(.)) 
  
  plot <- ggplot(data, aes(CellType, Exp)) +
    labs(title = title, y=ylab) +
    theme_bw(base_size = size) +
    theme(axis.text.x = element_text(angle = 40,hjust = 1),
          panel.grid = element_blank()) +
    scale_fill_manual(values = CellColors) +
    geom_violin(alpha=0.8, aes(fill = CellType)) +
    geom_boxplot(width=width, outlier.size = 0) +
    facet_wrap(~GeneType2, nrow = length(levels(data$CellType)))
  if(!save){
    print(plot)
  } else {
    ggsave(ggsave(paste0(path, "/MGPcellTypes_", title, ".pdf"),
                  plot = plot, width = 12, height = 8, units = "in",
                  dpi=300))
  }
}