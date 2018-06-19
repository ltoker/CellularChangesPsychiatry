Metadata %<>% droplevels
print("Batch correcting... tata - taaaam")
packageF("sva")
if("Batch" %in%  names(Metadata)){
  batches <- table(Metadata$Batch)
  
  bad_batch <- names(batches)[batches == 1]
  
  ProfileBatch <- data.frame(row.names = unique(Metadata$Batch))
  ProfileBatch$nSample <- sapply(rownames(ProfileBatch), function(batch){
    Metadata %>% filter(Batch == batch) %>% .$Profile %>% as.character() %>% length
  })
  
  ProfileBatch$nProfile <- sapply(rownames(ProfileBatch), function(batch){
    Metadata %>% filter(Batch == batch) %>% .$Profile %>% as.character() %>% unique %>% length
  })
  
  ProfileBatch$Confound <- apply(ProfileBatch, 1, function(batch){
    if(length(Metadata$Profile %>% unique) > 1 & batch[1] != 1 & batch[2] == 1 ){
      "Yes"
    } else {
      "No"
    }
  })
  
  #Run combat only if group factor is not confounded with batch
  if(sum(ProfileBatch$Confound == "Yes") == 0){
    #Excluse the samples that were scaned alone
    single_sample <- Metadata$CommonName[Metadata$Batch == bad_batch]
    
    aned_good <- aned_good[,!names(aned_good) %in% single_sample]
    if(!"CommonName" %in% names(Metadata)){
      Metadata$CommonName = Metadata$Series_sample_id
    }
    Metadata <- Metadata[!Metadata$CommonName %in% single_sample,]
    if(length(single_sample) > 0){
      print(paste0(single_sample, " removed since was scanned alone"))
    }
    
    batch = Metadata$Batch
    
    if(length(unique(batch)) > 1) {
      #Create the files for ComBat
      exp_file <- aned_good[,-c(1:3)]
      rownames(exp_file) <- aned_good$Probe
      
      
      if(length(levels(as.factor(Metadata$Profile))) > 1){
        mod=model.matrix(~as.factor(Profile),
                         Metadata[match(names(aned)[-c(1:3)], Metadata$CommonName),])
      } else {
        mod=matrix(rep(1, ncol(exp_file)), dimnames=list(names(exp_file),
                                                         "(Intercept)"))
      }
      
      #ComBat function
      aned_combat<- ComBat(dat=exp_file, mod=mod, batch=batch, par.prior=T, prior.plots=F)
      aned_good <- cbind(aned_good[,1:3],aned_combat)
    } else {
      print("All samples analysed during the same day")
    }
    rm(mod, aned_combat, exp_file, bad_batch, batch, batches, single_sample, list=ls(pat="temp"))
  } else {
    warning("Group factor confounded with batch - no batch correction")
  }
} else {
  batches <- table(Metadata$ScanDate)
  if(length(batches) < nrow(Metadata)/4) {
    bad_batch <- names(batches)[batches == 1]
  } else {
    print(paste(length(batches), "batches, cannot use scan day, group by month"))
    if(sum(grepl("-", Metadata$ScanDate)) < nrow(Metadata) & sum(grepl("/", Metadata$ScanDate)) < nrow(Metadata)){
      warning("Different time formats for ScanDate")
    }
    Format <- sapply(c("-", "/"), function(x){
      grep(x, Metadata$ScanDate)
    }, simplify = FALSE)
    temp3 <- sapply(names(Format), function(x){
      if(length(Format[[x]]) == 0){
        ""
      } else {
        temp <- sapply(Metadata$ScanDate[Format[[x]]], function(Date){
          strsplit(Date, x)[[1]]
        },simplify=FALSE) %>% do.call(rbind, .)
        temp2 <- temp %>% apply(c(1,2), as.numeric) %>% apply(2, max)
        if(max(nchar(temp2)) == 4){ #For 4 digit year formats
          Y_temp <- which(temp2 > 1000)
          M_temp <- which(temp2 <= 12)
          D_temp <- which(temp2 < 32 & temp2 > 12)
          temp %>% apply(1, function(x) paste(x[Y_temp], x[M_temp], sep="_"))
        } else { #For 2 digit year format                                                        
          Y_temp <- which(temp2 > 90) # Just in case there are arrays from the 90th
          M_temp <- which(temp2 <= 12)
          D_temp <- which(temp2 < 32 & temp2 > 12)
          temp %>% apply(1, function(x) paste(paste0("19",x[Y_temp]), x[M_temp], sep="_")) 
          if(length(Y_temp) == 0){ #This would be the case for most arrays
            Values <- apply(temp, 2, table) %>% lapply(length) %>% unlist %>% order
            Y_temp = Values[1] #The year is the least variable among the the three
            M_temp <- M_temp[M_temp != Y_temp]
            D_temp <- D_temp[D_temp != Y_temp]
            temp %>% apply(1, function(x) paste(paste0("20",x[Y_temp]), x[M_temp], sep="_")) 
          }
        }
      }
    }, simplify = FALSE) %>% unlist()
    names(temp3) <- sapply(names(temp3), function(x) strsplit(x, "\\.")[[1]][2])
    temp3 <- temp3[match(Metadata$ScanDate, names(temp3))]
    
    batches <- table(temp3)
    Metadata$ScanDate2 <- Metadata$ScanDate
    Metadata$ScanDate <- temp3
    bad_batch <- names(batches)[batches == 1]
  }
  
  ProfileBatch <- data.frame(row.names = unique(Metadata$ScanDate))
  ProfileBatch$nSample <- sapply(rownames(ProfileBatch), function(batch){
    Metadata %>% filter(ScanDate == batch) %>% .$Profile %>% as.character() %>% length
  })
  
  ProfileBatch$nProfile <- sapply(rownames(ProfileBatch), function(batch){
    Metadata %>% filter(ScanDate == batch) %>% .$Profile %>% as.character() %>% unique %>% length
  })
  
  ProfileBatch$Confound <- apply(ProfileBatch, 1, function(batch){
    if(length(Metadata$Profile %>% unique) > 1 & batch[1] != 1 & batch[2] == 1 ){
      "Yes"
    } else {
      "No"
    }
  })
  
  #Run combat only if group factor is not confounded with batch
  
  if(sum(ProfileBatch$Confound == "Yes") == 0){
    #Excluse the samples that were scaned alone
    single_sample <- Metadata$CommonName[Metadata$ScanDate == bad_batch]
    
    aned_good <- aned_good[,!names(aned_good) %in% single_sample]
    if(!"CommonName" %in% names(Metadata)){
      Metadata$CommonName = Metadata$Series_sample_id
    }
    Metadata <- Metadata[!Metadata$CommonName %in% single_sample,]
    if(length(single_sample) > 0){
      print(paste0(single_sample, " removed since was scanned alone"))
    }
    
    batch = Metadata$ScanDate
    
    if(length(unique(batch)) > 1) {
      #Create the files for ComBat
      exp_file <- aned_good[,-c(1:3)]
      rownames(exp_file) <- aned_good$Probe
      
      
      if(length(levels(as.factor(Metadata$Profile))) > 1){
        mod=model.matrix(~as.factor(Profile),
                         Metadata[match(names(aned)[-c(1:3)], Metadata$CommonName),])
      } else {
        mod=matrix(rep(1, ncol(exp_file)), dimnames=list(names(exp_file),
                                                         "(Intercept)"))
      }
      
      #ComBat function
      aned_combat<- ComBat(dat=exp_file, mod=mod, batch=batch, par.prior=TRUE, prior.plots=F)
      aned_good <- cbind(aned_good[,1:3],aned_combat)
    } else {
      print("All samples analysed during the same day")
    }
    rm(mod, aned_combat, exp_file, bad_batch, batch, batches, single_sample, list=ls(pat="temp"))
  } else {
    warning("Group factor confounded with batch - no batch correction")
  }
}
