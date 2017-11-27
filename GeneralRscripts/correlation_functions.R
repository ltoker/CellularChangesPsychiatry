data_set<-rev(as.vector(strsplit(getwd(),"/"))[[1]])[1]
require("corpcor")
#FUNCTIONS

Connections_Num <-function(samples, group){
  #Create correlation matrix of all the genes on the array with signal above the threshold
  data.frame(do.call("rbind", by(aned_high[,samples], aned_high$GeneSymbol, colMeans)))-> try_genes
  m_genes2 <- t(try_genes)
  cor<-cor(m_genes2)
  cor_all_genes<-cor
  
  #Creatre adjastecy matrix for all the array genes
  cor_all_vec<-sm2vec(cor_all_genes)
  
  #Correlations histogram
  png(paste(group,"correlation_histogram.png", sep="_"))
  hist(cor_all_vec, main=paste(data_set, group, sep="_"))
  dev.off()
  
  treshold = quantile(abs(cor_all_vec), 0.95)
  #treshold = mean(abs(cor_all_vec))+2*sd(abs(cor_all_vec))
  adj_cor<-cor
  adj_cor[which(abs(adj_cor)<treshold)]<-0
  diag(adj_cor)<-0
  adj_cor[which(abs(adj_cor)>=treshold)]<-1
  
  
  #Create matrix with the connection number for each gene on the array (from aned_high)
  connections_num <- rev(sort(rowSums(adj_cor)))
  max_hist=(round(max(connections_num)/100)+1)*100
  png(paste(group,"corr_num_all.png", sep="_"))
  hist_all_genes <- rev(hist(connections_num, breaks=seq(0, max_hist,100), main=paste(data_set,group, sep="_"), labels=T, )[[2]])
  dev.off()
  
  #Create list of genes organised by their number of connections (breaks of 100)
  
  cor_hist_list <- list()
  
  
  cor_hist_list$con_0_99 <- connections_num[which(connections_num >= 1 &connections_num < 100)]
  cor_hist_list$con_100_199 <- connections_num[which(connections_num >= 100 & connections_num < 200)]
  cor_hist_list$con_200_299 <- connections_num[which(connections_num >= 200 & connections_num < 300)]
  cor_hist_list$con_300_399 <- connections_num[which(connections_num >= 300 & connections_num < 400)]
  cor_hist_list$con_400_499 <- connections_num[which(connections_num >= 400 & connections_num < 500)]
  cor_hist_list$con_500_599 <- connections_num[which(connections_num >= 500 & connections_num < 600)]
  cor_hist_list$con_600_699 <- connections_num[which(connections_num >= 600 & connections_num < 700)]
  cor_hist_list$con_700_799 <- connections_num[which(connections_num >= 700 & connections_num < 800)]
  cor_hist_list$con_800_899 <- connections_num[which(connections_num >= 800 & connections_num < 900)]
  cor_hist_list$con_900_999 <- connections_num[which(connections_num >= 900 & connections_num < 1000)]
  cor_hist_list$con_1000_1099 <- connections_num[which(connections_num >= 1000 & connections_num < 1100)]
  cor_hist_list$con_1100_1199 <- connections_num[which(connections_num >= 1100 & connections_num < 1200)]
  cor_hist_list$con_1200_1299 <- connections_num[which(connections_num >= 1200 & connections_num < 1300)]
  cor_hist_list$con_1300_1399 <- connections_num[which(connections_num >= 1300 & connections_num < 1400)]
  cor_hist_list$con_1400_1499 <- connections_num[which(connections_num >= 1400 & connections_num < 1500)]
  cor_hist_list$con_1500_1599 <- connections_num[which(connections_num >= 1500 & connections_num < 1600)]
  cor_hist_list$con_1600_1699 <- connections_num[which(connections_num >= 1600 & connections_num < 1700)]
  cor_hist_list$con_1700_1799 <- connections_num[which(connections_num >= 1700 & connections_num < 1800)]
  cor_hist_list$con_1800_1899 <- connections_num[which(connections_num >= 1800 & connections_num < 1900)]
  cor_hist_list$con_1900_1999 <- connections_num[which(connections_num >= 1900 & connections_num < 2000)]
  cor_hist_list$con_2000_2099 <- connections_num[which(connections_num >= 2000 & connections_num < 2100)]
  cor_hist_list$con_2100_2199 <- connections_num[which(connections_num >= 2100 & connections_num < 2200)]
  cor_hist_list$con_2200_2299 <- connections_num[which(connections_num >= 2200 & connections_num < 2300)]
  cor_hist_list$con_2300_2399 <- connections_num[which(connections_num >= 2300 & connections_num < 2400)]
  cor_hist_list$con_2400_2499 <- connections_num[which(connections_num >= 2400 & connections_num < 2500)]
  cor_hist_list$con_2500_2600 <- connections_num[which(connections_num >= 2500)]
  cor_hist_list$con_1400_max <- connections_num[which(connections_num >= 1600)]
  
  
  
  #Find the number of connection for the genes in "genes" vector
  genes_conn_num <- connections_num[which(names(connections_num) %in% genes)]
  
  max_hist_genes=(round(max(genes_conn_num)/100)+1)*100
  png(paste(group,"cor_num_genes.png", sep="_"))
  hist_genes <- hist(genes_conn_num, breaks=seq(0, max_hist_genes,100), main=paste(data_set,group, sep="_"))[[2]]
  dev.off()
  print(hist_genes)
  
  Connections_Num_list<-list(cor_all_genes, treshold, adj_cor, connections_num, cor_hist_list, genes_conn_num, hist_genes)
  names(Connections_Num_list)<-c("cor_all_genes", "treshold", "adj_cor", "connections_num", "cor_hist_list", "genes_conn_num", "hist_genes")
  
  return(Connections_Num_list)
}


#drow function
drow_g<-function(genes){
  aned_high[which(aned_high[,2] %in% genes),] -> m_genes
  return(m_genes)
}

#signal_average function
signal_average<-function(m_genes){
  data.frame(do.call("rbind", by(m_genes[,samples], m_genes$GeneSymbol, colMeans)))-> try_genes
  return(try_genes)
}


#correlation function
cor_func<-function(try_genes) {
  
  #Create correlation matrix
  t(try_genes)->m_genes2
  cor(m_genes2)->cor
  
  #Fishers's transformaion of the r values
  0.5*log((1+cor)/(1-cor))-> cor
  
  #Create correlation vector
  sm2vec(cor)->correlations
  
  cor_list<-list(cor,correlations)
  return(cor_list)
}

#Connectivity function
network_connectivity <- function(test_genes){
  adj_cor <- Connections_Num_list$adj_cor
  connections_num <- Connections_Num_list$connections_num
  genes_conn_num <- connections_num[which(names(connections_num) %in% test_genes)]
  genes_conn_num <- genes_conn_num[sort(names(genes_conn_num))]
  network_connections <- adj_cor[which(rownames(adj_cor) %in% test_genes), which(colnames(adj_cor) %in% test_genes)]
  network_conn_num <- apply(network_connections, 1, sum)
  network_conn_num <- network_conn_num[sort(names(network_conn_num))]
  normalized_connections <- network_conn_num/genes_conn_num
  return(normalized_connections)
}


#Find genes with extreme values in bipolar based on z-scores
extreme_genes <- function(z_scores, bipolar, controls, Z, BP_num) {
  compared_z <- vector(mode="numeric")
  for(i in 1:nrow(z_scores)){
    if(max(z_scores[i, c(bipolar-3)]) > Z){
      compared_z[i] <- length(which(z_scores[i, c(bipolar-3)] > max(z_scores[i,c(controls-3)]) &
                                    z_scores[i, c(bipolar-3)] > Z))
      } else if (min(z_scores[i, c(bipolar-3)]) < -Z){
        compared_z[i] <- length(which(z_scores[i, c(bipolar-3)] < min(z_scores[i,c(controls-3)]) &
                                    z_scores[i, c(bipolar-3)] < -Z))
        } else {
          compared_z[i] <- 0
        }
  }
  names(compared_z)<-rownames(z_scores)
  extreme_z <- names(compared_z)[which(compared_z > BP_num)]
  return(extreme_z)
}

extreme_genes_mod_z <- function(mod_z_scores, bipolar, controls, Z, BP_num) {
  compared_mod_z <- vector(mode="numeric")
  for(i in 1:nrow(mod_z_scores)){
    if(max(mod_z_scores[i, c(bipolar-3)]) > Z){
      compared_mod_z[i] <- length(which(mod_z_scores[i, c(bipolar-3)] > max(mod_z_scores[i,c(controls-3)]) &
                                        mod_z_scores[i, c(bipolar-3)] > Z))
      } else if (min(mod_z_scores[i, c(bipolar-3)]) < -Z){
        compared_mod_z[i] <- length(which(mod_z_scores[i, c(bipolar-3)] < min(mod_z_scores[i,c(controls-3)]) &
                                        mod_z_scores[i, c(bipolar-3)] < -Z))
        } else {
          compared_mod_z[i] <- 0
        }
  }
  names(compared_mod_z)<-rownames(mod_z_scores)
  extreme_mod_z <- names(compared_mod_z)[which(compared_mod_z > BP_num)]
  return(extreme_mod_z)
}


rand_Z<-function(samples,group){
  
  source("Coexpression_connections.R")
  
  list(samples,cor,correlations,rand_mean,Z)->par_list
  
  print(paste(group, "Z is", as.name(Z), sep=" "))
  
  return(par_list)
}

rand_Z2<-function(samples,group){
  
  source("Coexpression_connections2.R")
  
  par_list <- list(samples, normalized_connections,rand_median, Z, rand_mean)
  
  print(paste(group, "Z is", as.name(Z), sep=" "))
  
  return(par_list)
}
#####################################################################################################