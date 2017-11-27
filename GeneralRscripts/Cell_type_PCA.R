packageF("devtools")
#packageF("listwiewer")
packageF("affy")
packageF("limma")

if(!"homologene" %in% rownames(installed.packages())){
  install_github("oganm/homologene", force = T)
}
library(homologene)

if(!"markerGeneProfile" %in% rownames(installed.packages())){
  install_github("oganm/markerGeneProfile", force = T)
}
library(markerGeneProfile)
data("mouseMarkerGenesCombined")

correct_sign <- function(Celltype, object){
  sign = sum(object$rotation[,1])  
  if(sign < 0) {
    object$rotation[,1] <- -1*object$rotation[,1]
    object$x[,1] <- -1*object$x[,1]
  }
  
  return(object)
}

PCA_genes <- function(dataset_id, dataset, CellType_genes, extreme_genes, extreme_marker, metadata){ #Atp5g is marker for the extreme genes
  
  controls <- grep("Cont", colnames(dataset) )
  bipolar <- grep("BP", colnames(dataset))
  scz <- grep("SCZ", colnames(dataset))
  md <- grep("MD", colnames(dataset))
  
  PCA_all <- list()
  PCA_all[[1]] <- dataset_id
  for(i in 1:c(length(names(CellType_genes))-1)){
    if(length(which(rownames(dataset) %in% CellType_genes[[i]])) > 2){
      
      PCA_all[[c(i+1)]] <- prcomp(t(dataset[which(rownames(dataset) %in% CellType_genes[[i]]),c(bipolar,controls,scz)]), scale=T)
      
      if(length(PCA_all[[c(i+1)]]$rotation[,1]) > 2){
        PCA_all[[c(i+1)]] <- correct_sign(CellType_genes$Marker[i], PCA_all[[c(i+1)]])
      }else{
        print(paste("Looki here!!!! Less than three genes for", names(CellType_genes)[i], "ignore PCA"))
      }
      
    } else {
      print(paste("Looki here!!! No genes for", names(CellType_genes)[i]))
      PCA_all[[c(i+1)]] <- "NA"
    }               
  }
  names(PCA_all) <- c("dataset_id", names(CellType_genes)[-(length(names(CellType_genes)))])
  
  PCA_all$Extreme_genes <- prcomp(t(dataset[which(rownames(dataset) %in% extreme_genes),c(bipolar,controls, scz)]), scale=T)
  PCA_all$Extreme_genes <- correct_sign(extreme_marker, PCA_all$Extreme_genes)
  
  return(PCA_all)
}

PCA_genes_All_based <- function(dataset_id, dataset, CellType_genes, NoiseThershold, contName = "Cont"){ 
  print("########################################")
  print("##### Scales based on all samples ####")
  print("########################################")
  groups <- sapply(grep("_|GSM", names(dataset), value=T),
                   function(x) gsub("_.*", "", x)) %>% table %>% names
  for(grp in groups){
    assign(grp, grep(grp, names(dataset)))
  }
  PCAresults <- list()
  PCAresults$dataset_id <- dataset_id
  PCAresults$ControlOnly <- list()
  PCAresults$All <- list()
  PCAresults$modified <- list()
  for(i in 1:c(length(names(CellType_genes)))){
    print(names(CellType_genes)[i])
    if(sum(dataset$GeneSymbol %in% CellType_genes[[i]]) > 2){
      data <- dataset %>% filter(GeneSymbol %in% CellType_genes[[i]])
      #Remove genes with expression level bellow the noise threshold in 1/3 of control samples. Ths step is done to ensure that the marker genes can be detected in human bulk tissue
      geneRM <- sapply(data$GeneSymbol, function(gene){
        if(data %>% filter(GeneSymbol == gene) %>% select(matches(contName)) %>% unlist %>% quantile( 0.67) < NoiseThershold){
          TRUE
        } else {
          FALSE
        }
      })
      data <- data[!geneRM,]
      PCAresults$ControlOnly[[i]] <- data %>%
        select(matches(contName)) %>% t %>% prcomp(scale=T) #The assumption here is that there is a control group
      rownames(PCAresults$ControlOnly[[i]]$rotation) <- data$GeneSymbol
      PCAresults$ControlOnly[[i]] <- correct_sign(names(CellType_genes)[i], PCAresults$ControlOnly[[i]])
      PCAresults$All[[i]] <- data %>% select(matches("_|GSM")) %>% t %>% prcomp(scale=TRUE) 
      PCAresults$All[[i]] <- correct_sign(names(CellType_genes)[i], PCAresults$All[[i]])
      rownames(PCAresults$All[[i]]$rotation) <- data$GeneSymbol
      while (sum(PCAresults$All[[i]]$rotation[,1] > 0) < nrow(PCAresults$All[[i]]$rotation)){
        minorGene <- rownames(PCAresults$All[[i]]$rotation)[PCAresults$All[[i]]$rotation[,1] < 0]
        data %<>% filter(!GeneSymbol %in% minorGene)
        rownames(data) <- data$GeneSymbol
        PCAresults$All[[i]] <- data %>% select(matches("_|GSM")) %>% t %>% prcomp(scale=TRUE) 
        PCAresults$All[[i]] <- correct_sign(names(CellType_genes)[i], PCAresults$All[[i]])
      }
      
      PCAresults$modified[[i]] <- PCAresults$All[[i]]$x %>% list
      PCAresults$modified[[i]] <- apply(PCAresults$modified[[i]][[1]],2,
                                        function(x) rescale(x,c(0,1))) %>% list
      names(PCAresults$modified[[i]]) <- "x"
      
    } else {
      print(paste("Looki here!!! No genes for", names(CellType_genes)[i]))
      x <- matrix(nrow = c(sum(grepl("_|GSM", names(data)))), ncol=1)
      rownames(x) <- names(dataset)[sapply(groups, function(grp){
        eval(as.name(grp))
      }) %>% unlist] ;colnames(x)="x"
      PCAresults$ControlOnly[[i]] <- list(x)
      names(PCAresults$ControlOnly[[i]]) <- "x"
      PCAresults$All[[i]] <- list(x)
      names(PCAresults$All[[i]]) <- "x"
      PCAresults$modified[[i]] <- list(x)
      names(PCAresults$modified[[i]]) <- "x"
    }
    
  }
  for(i in 2:4){
    names(PCAresults[[i]]) <- names(CellType_genes)
  }
  return(PCAresults)
}


PlotPCggplot <- function(data, grpVar = "Profile", CellVar=NULL, 
                         CellExtend = "_Genes", CellName = NULL,
                         grpRef = "Cont",
                         name="studyName", stat="wilcox",
                         txtSize=12, ptSize=3, pValSize=NULL,
                         grpColors = c("burlywood3", "cornflowerblue", "indianred4")){
  if("study" %in% names(data)){
    data$study <- factor(data$study)
  } else {
    data$study <- factor(name)
  }
  
  #If the name of the cell type is not specified, set it based on the variable name
  if(is.null(CellName)){
    CellName = gsub(CellExtend, "", CellVar)
  }
  
  #If the variable corresponding to the cell type is not specified, set it based on the cell name + extension
  if(is.null(CellVar)){
    CellVar = paste0(CellName,CellExtend)
  }
  
  pvalAll <- sapply(levels(data$study), function(std){
    sapply(levels(data[[grpVar]])[-1], function(grp){
      data = data %>% filter_(paste0(grpVar, " %in% c('",grpRef,"','", grp, "')"))
      WlCx <- wilcox.test(formula(paste0(CellVar, "~", grpVar)), data=data)
      if(WlCx$p.value > 0.05){
        pval=paste0("p = ", round(WlCx$p.value, digits=2))
      } else{
        pval = paste0("p = ", scientific(WlCx$p.value))
      }
      char=GetSigChar(WlCx$p.value)
      paste0(char, pval)
    })
  })
  
  if(is.null(pValSize)){
    pValSize = 0.25*txtSize
  }
  
  rmSamples <- sapply(data$Filename, function(smpl){
    temp1 = data %>% filter(Filename == smpl) %>% select(matches("_Genes")) %>% is.na %>% sum
    temp2 = grepl("_Genes", names(data)) %>% sum 
    if(temp2 == temp1){
      as.character(smpl)
    }
  }) %>% unlist
  data %<>% filter(!Filename %in% rmSamples) %>% droplevels()
  
  data$Profile2 <- sapply(unique(data$study), function(std){
    data2 <- data %>% filter(study == std)
    GroupSum <- group_by_(data2, .dots=grpVar) %>% summarise(n=n())
    data2$Profile2 <- sapply(data2[[grpVar]], function(grp){
      paste0(grp, "\n(n=", GroupSum %>% filter_(paste0(grpVar, "=='", grp, "'")) %>% .$n, ")")
    }) %>% factor %>% relevel(ref = paste0(grpRef,"\n(n=", GroupSum %>% filter_(paste0(grpVar,"=='",grpRef, "'")) %>% .$n, ")"))
    
  }, simplify=FALSE) %>% unlist
  title = paste(CellName, "marker gene profile")
  
  plot <- ggplot(data, aes_string("Profile2", CellVar)) +
    labs(title = title, x = "", y="Relative profile")+
    theme_grey(base_size = txtSize) +
    theme(axis.text.y = element_text(size = rel(0.8)),
          axis.text.x = element_text(size = rel(1.2)),
          panel.grid = element_blank()) +
    scale_fill_manual(values = grpColors, name="Group")+
    geom_violin(aes_string(fill=grpVar), width=0.9, size=0.1, alpha=0.6) +
    geom_boxplot(outlier.shape = NA, width = 0.25, colour="black", alpha=0.8) +
    geom_jitter(colour="black", size=ptSize,width=0.1, shape=16, alpha=0.5)+
    facet_wrap(~study, scales = "free_x")+
    geom_text(data = data.frame(x=c(2:length(levels(data[[grpVar]]))), y=1.1,
                                label=pvalAll,
                                study=levels(data$study)),
              aes(x,y,label=pvalAll),size=pValSize, inherit.aes = FALSE)
  return(plot)
}


GetMarkers <- function(region){
  #path = paste0(path, region)
  region=as.character(region)
  CellType_genes <- mouseMarkerGenesCombined[[region]]
  #just for the cerebellum, add Astrocyte genes from to cortex to compare with Bergmann glia results
  #for hippocampus - add cortex Astrocyte from cortex just to make sure the robustness of the results
  
  if(region != "Cortex"){
    CellType_genes$AstrocyteCortex <-  mouseMarkerGenesCombined$Cortex$Astrocyte
  }
  
  for(i in 1:length(CellType_genes)){
    
    CellType_genes[[i]] <- as.vector(mouse2human(CellType_genes[[i]])$humanGene)
    
  }
  
  names(CellType_genes) <- sapply(names(CellType_genes), function(x) paste(x, "Genes", sep="_"))
  
  if(region == "Cortex"){
    #Exclude specific genes: 5HTR3A - also expressed in VIP negative 5HTR3A cells.
  CellType_genes$GabaVIPReln_Genes <- CellType_genes$GabaVIPReln_Genes[!CellType_genes$GabaVIPReln_Genes %in% "HTR3A"]
  }
  
  return(CellType_genes)
} 
