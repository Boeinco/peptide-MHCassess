peps <- read.csv("8mers_uniq.csv",  header=FALSE)
peps_features <- read.csv("8mers_features_all_uiq.csv")
set.seed(1)
pca <-prcomp(peps_features, scale=TRUE)
pcs <- as.data.frame(pca$x[,1:2])
#require(ggfortify)
pc_table <- cbind(peps, pcs)
colnames(pc_table)[1] <- "Peptide"
rm(pca)
rm(pcs)
rm(peps_features)
rm(peps)

setwd("/Users/nguyenau/Documents/Coding/PeptideAlleleBinding/")
merged_peps_umap8 <- read.csv("merged_peps_umap8mer.csv")
merged_peps_umap8 <- merged_peps_umap8[,c(1,6)]


merged_pc <- merge(pc_table, merged_peps_umap8, by="Peptide", all=TRUE)
merged_pc <- distinct(merged_pc)

i <- 1
pca_matrix_bins_8mers <- tibble(x=character(), y=character(), diff=numeric())
names <- sort(unique(merged_pc$Source))
library(plyr)
merged_pc$Source <- mapvalues(merged_pc$Source, matrix_names[,1], matrix_names[,2])
ab <- matrix( c(min(merged_pc$PC1),min(merged_pc$PC2),max(merged_pc$PC1),max(merged_pc$PC2)), 2, 2)      
library(ash)
#library(ggplot2)
for (i in 1:length(names))
{
  for (j in i:length(names))
  {
    x <- as.matrix(merged_pc[merged_pc$Source == names[i],2:3])  
    nbin <- c( 40, 40)          
    bins <- bin2(x, ab, nbin)
    bin_sol <- bins$nc
    bin_sol1 <- bin_sol/(dim(x)[1])
    
    y <- as.matrix(merged_pc[merged_pc$Source == names[j],2:3])  
    #x <- matrix( rnorm(200), 100 , 2)
    bins <- bin2(y, ab, nbin)
    bin_sol <- bins$nc
    bin_sol2 <- bin_sol/(dim(y)[1])
    diff_bins <- bin_sol1 - bin_sol2
    sum_bins <- sum(abs(diff_bins))
    pca_matrix_bins_8mers <- pca_matrix_bins_8mers %>% add_row(x=names[i], y=names[j], diff=sum_bins)
  }
}

pca_matrix_bins_8mers$diff <- pca_matrix_bins_8mers$diff * 100/2
#umap_bins_8mers <- read.csv("matrix_bins_8mers.csv")
merged_bins <- merge(pca_matrix_bins_8mers, umap_bins_8mers, by=c("x","y"), all=TRUE)
write.csv(merged_bins, "pca_vs_umap_8mers_same.csv", row.names=FALSE, quote=FALSE)

ggplot(data = pca_matrix_bins_8mers, aes(y, x, fill = diff))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(0,100), 
                       name="Percent Peptides Difference") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ggtitle("8mer Differential Peptides PCA")
