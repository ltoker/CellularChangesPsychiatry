
#All kinf of graphich functions
###########################################################################################
require("igraph")
require("gplots")
require("sm")
#Create i_grapj from adjustenct matrix
network_g<-function(samples, treshold, group, data_set){
aned_high[which(aned_high[,2] %in% genes),] -> m_genes
data.frame(do.call("rbind", by(m_genes[,samples], m_genes$GeneSymbol, colMeans)))-> try_genes
cor(t(try_genes))-> cor
cor[which(abs(cor) < treshold)]<-0
diag(cor)<-0
connections_num <- rev(sort(rowSums(cor)))
graph.adjacency(cor, mode="upper", weighted=TRUE)-> ig_group
png(paste(group,"Network.png", sep="_"))
plot(ig_group, main=paste(data_set, group, sep="_"))
dev.off()

return(ig_group)
}


#Calculate Z-score
z_scores_genes <- function(samples){
  m_genes <- aned_high[which(aned_high[,2] %in% genes),] 
  try_genes <- data.frame(do.call("rbind", by(m_genes[,samples], m_genes$GeneSymbol, colMeans)))
  mean_expression <- apply(try_genes, 1, mean)
  apply(try_genes, 1, sd) -> sd_expression
  try_genes_z <- matrix(nrow=nrow(try_genes),ncol=ncol(try_genes)) 
  for(i in 1:nrow(try_genes)){
    try_genes_z[i,] <- sapply(try_genes[i,], function(x) as.vector(x-mean_expression[i])/sd_expression[i]) 
  }
  colnames(try_genes_z)<-colnames(try_genes)
  rownames(try_genes_z)<-rownames(try_genes)
  my_palette <- colorRampPalette(c("#123766", "#FFF5D4","#7D0D19"))(n = 999)
  heatmap.2(try_genes_z, Rowv=T, Colv=T, dendrogram="both", trace="none", density.info="none", 
            col=my_palette, keysize=0.8,main= paste(data_set, "Z scores genes", data, sep=" "), cex.main = 0.7) -> Z_scores_genes
  return(try_genes_z)
}

z_scores_ALL_genes <- function(samples, signals, what_data){
  try_genes <- data.frame(do.call("rbind", by(signals[,samples], signals$GeneSymbol, colMeans)))
  mean_expression <- apply(try_genes, 1, mean)
  apply(try_genes, 1, sd) -> sd_expression
  try_genes_z <- matrix(nrow=nrow(try_genes),ncol=ncol(try_genes)) 
  for(i in 1:nrow(try_genes)){
    try_genes_z[i,] <- sapply(try_genes[i,], function(x) as.vector(x-mean_expression[i])/sd_expression[i]) 
  }
  colnames(try_genes_z)<-colnames(try_genes)
  rownames(try_genes_z)<-rownames(try_genes)
  my_palette <- colorRampPalette(c("#123766", "#FFF5D4","#7D0D19"))(n = 999)
  png(paste(what_data, "png", sep="."), width = 980, height = 700, units = "px",  bg = "white")
  heatmap.2(try_genes_z, Rowv=T, Colv=T, dendrogram="both", trace="none", density.info="none", 
            col=my_palette, keysize=0.8,main= paste(data_set, "Z scores", what_data, sep=" ")) -> Z_scores_ALL_genes
  dev.off()
  return(Z_scores_ALL_genes)
}

MOD_z_scores_ALL_genes <- function(samples, signals, what_data){
  try_genes <- signals[samples]
  median_expression <- apply(try_genes, 1, median)
  MAD <- apply(try_genes, 1, function(x) mad(x, center = median(x), constant = 1.4826))
  try_genes_MAD <- matrix(nrow=nrow(try_genes),ncol=ncol(try_genes)) 
  for(i in 1:nrow(try_genes)){
    try_genes_MAD[i,] <- sapply(try_genes[i,], function(x) as.vector(x-median_expression[i])/MAD[i])
  }
  colnames(try_genes_MAD)<-colnames(try_genes)
  rownames(try_genes_MAD)<-rownames(try_genes)
  #remove genes with median of 0
  try_genes_MAD <- try_genes_MAD[!MAD == 0, ]
  my_palette <- colorRampPalette(c("#123766", "#FFF5D4","#7D0D19"))(n = 999)
  png(paste(what_data, "png", sep="."), width = 980, height = 700, units = "px",  bg = "white")
  heatmap.2(try_genes_MAD, Rowv=T, Colv=T, dendrogram="both", hclustfun= function(x) hclust(x, method="ward.D2"), trace="none", density.info="none", 
            col=my_palette, keysize=0.8,main= paste(data_set, "MOD Z scores", what_data, sep=" ")) -> MOD_Z_scores_ALL_genes
  dev.off()
  return(MOD_Z_scores_ALL_genes)
}  

#Heatmap samples correlations
samples_cor <- function(){
  all_samples_cor <- cor(aned_high[,4:ncol(aned_high)])
  my_palette <- colorRampPalette(c("#123766", "#FFF5D4","#7D0D19"))(n = 999)
  unique_scan <- unique(Samples[3,])
  scan_colors<-rich.colors(length(unique_scan))
  unique_pH <- unique(sort(round(as.numeric(Samples[9,]), digit=1)))
  pH_colors <-  colorRampPalette(c("white", "grey", "black"))(length(unique_pH))
  scan_date_col <- character()
  for(i in 1:length(unique(Samples[3,]))){
    scan_date_col[which(Samples[3,] %in% unique_scan[i])]<-scan_colors[i]
  }
  
  pH_col <- character()
  for(i in 1:length(unique(Samples[9,]))){
    pH_col[which(round(as.numeric(Samples[9,]), digit=1) %in% unique_pH[i])]<-pH_colors[i]
  }
  png("Samples_correlation_both_combat_sd.png", width = 980, height = 700, units = "px",  bg = "white")
  heatmap.2(all_samples_cor, Rowv=T, Colv=T, symm=T, dendrogram="none", trace="none", density.info="none", col=my_palette,
  keysize=0.6, RowSideColors=scan_date_col, ColSideColors=pH_col, cexCol=1.2, cexRow=1.2, main= paste(data_set, "Both excl + ComBat + sd", sep=" "))
  legend("left",legend=unique_scan, fill=scan_colors, border=FALSE, bty="n", y.intersp = 1, x.intersp = 0.3, cex=0.9, title = "Batch")
  legend("top",legend=unique_pH, fill=pH_colors, border=FALSE, bty="n", x.intersp = 1, y.intersp=0.8, cex=0.9, title = "pH", horiz=T)
  dev.off()
}


#Heatmap - correlations
heatmap_cor <- function(samples, group){
  m_genes <- aned_high[which(aned_high$GeneSymbol %in% genes),] 
  try_genes <- data.frame(do.call("rbind", by(m_genes[,samples], m_genes$GeneSymbol, colMeans)))
  cor(t(try_genes))-> try_group_quantile
  my_palette <- colorRampPalette(c("#123766", "#FFF5D4","#7D0D19"))(n = 999)
  png(paste(group,"correlations heatmap.png", sep="_") , width = 980, height = 700, units = "px",  bg = "white")
  heatmap.2(try_group_quantile, Rowv=T, dendrogram="row", symm=T, trace="none", density.info="none", col=my_palette,
  keysize=1, breaks=c(seq(-1,-0.5, length=250), seq(-0.49,0.5, length=500), seq(0.51,1,length=250)), symbreaks=T,
  main=paste(data_set, group))
  dev.off()
}


#Heatmap - correlation ranks
heatmap_cor_ranks <- function(cor_all_genes){
  #Creating correlations ranks for all gene pairs
  m_genes <- aned_high[which(aned_high[,2] %in% genes),] 
  try_genes <- data.frame(do.call("rbind", by(m_genes[,samples], m_genes$GeneSymbol, colMeans)))
  try_group_quantile <- cor(t(try_genes))
  cor_ranks_group <- rank(cor_all_genes) #might kill R...
  ranks_m_group <- matrix(nrow=nrow(cor_all_genes),ncol=ncol(cor_all_genes))
  colnames(ranks_m_group)<-colnames(cor_all_genes)
  rownames(ranks_m_group)<-rownames(cor_all_genes)
  ranks_m_group[1:length(cor_ranks_group)]<-cor_ranks_group
  gene_cor_ranks_group <- ranks_m_group[which(rownames(ranks_m_group) %in% rownames(try_group_quantile)),
                                  which(colnames(ranks_m_group) %in% colnames(try_group_quantile))]
  
  #Creating the heatmap
  my_palette <- colorRampPalette(c("#123766", "#FFF5D4","#7D0D19"))(n = 999)
  cor_ranks_vec_group <- sm2vec(cor_ranks_group)
  q_25 <-quantile(cor_ranks_vec_group, 0.25)
  q_75 <-quantile(cor_ranks_vec_group, 0.75)
  min <- min(cor_ranks_vec_group)
  max <- max(cor_ranks_vec_group)
  png(paste(group,"cor_ranks_heatmap.png", sep="_"), width = 980, height = 700, units = "px",  bg = "white")
  heatmap.2(gene_cor_ranks_group, Rowv=T, dendrogram="row", symm=T, trace="none", density.info="none", col=my_palette,
  keysize=1, breaks=c(seq(min, q_25, length=250), seq(c(q_25+1),q_75, length=500), 
  seq(c(q_75+1), max, length=250)), symbreaks=F, main=paste(data_set, group, "Correlation ranks", sep=" "))
  dev.off()
  return(ranks_m_group)
  
}