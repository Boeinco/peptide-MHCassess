peps <- read.csv("8mers_uniq.csv",  header=FALSE)
peps_features <- read.csv("8mers_features_all_uiq.csv")
library(uwot)

umaps <- umap(peps_features)
merged_peps_umap8 <- cbind(peps, umaps)

peps <- read.csv("9mers_uniq.csv",  header=FALSE)
peps_features <- read.csv("9mers_features_all_uiq.csv")

umaps <- umap(peps_features)
merged_peps_umap9 <- cbind(peps, umaps)

peps <- read.csv("10mers_uniq.csv",  header=FALSE)
peps_features <- read.csv("10mers_features_all_uiq.csv")

umaps <- umap(peps_features)
merged_peps_umap10 <- cbind(peps, umaps)

peps <- read.csv("11mers_uniq.csv",  header=FALSE)
peps_features <- read.csv("11mers_features_all_uiq.csv")

umaps <- umap(peps_features)
merged_peps_umap11 <- cbind(peps, umaps)

