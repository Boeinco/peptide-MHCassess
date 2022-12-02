library(dplyr)
library(reshape2)
library(ggplot2)
matrix_names <- matrix(c("bk", "covid", "betaherpes_6a", "human", "random",
                         "herpes_1_17", 
                         "herpes_2_hg52", 
                         "ebv", 
                         "cmv","BK", "Covid", "HHV-6","Human", "Random", "HSV-1","HSV-2","HSV-4", "HHV-5"),
                       nrow=9, ncol=2)

big_allele_list <- c("HLA-A02:02","HLA-A02:06","HLA-A68:01","HLA-A02:03","HLA-A02:01","H-2-Kb","HLA-A11:01","HLA-A30:02","HLA-A31:01","HLA-B15:01","HLA-B35:01","HLA-A30:01","HLA-B07:02","HLA-A03:01","HLA-A29:02","HLA-A68:02","HLA-B08:01","HLA-B58:01","HLA-A23:01","HLA-A24:02","H-2-Db","HLA-A33:01","HLA-B40:01","HLA-B27:05","HLA-B44:02","HLA-B57:01","HLA-A01:01","HLA-B51:01","HLA-A26:01","HLA-A69:01","HLA-B18:01")
bighlathena_list <- gsub("HLA-", "", big_allele_list)
bighlathena_list <- gsub(":", "", bighlathena_list)

##Figure 1

#creating merged random for all tools
mhcnuggets_random <- read.csv("mhcnuggets_random.csv", header = FALSE)
mhcnuggets_random$V2 <- 1-log(mhcnuggets_random$V2)/log(50000)
mhcflurry_random <- read.csv("mhcflurry_random_filt.csv")
mhcflurry_random$mhcflurry_affinity <- 1-log(mhcflurry_random$mhcflurry_affinity)/log(50000)
netmhcpan_random <- read.csv("netmhcpan_joint_alleles_random.csv")
netmhcpan_random$Allele <- gsub("HLA-", "", netmhcpan_random$Allele)
netmhcpan_random$Allele <- gsub(":", "", netmhcpan_random$Allele)
colnames(mhcnuggets_random) <- c("Peptide", "mhcnuggets_affinity", "Allele")
mhcnuggets_random$Allele <- gsub("HLA-", "", mhcnuggets_random$Allele)
mhcnuggets_random$Allele <- gsub(":", "", mhcnuggets_random$Allele)
mhcnuggets_random$Allele <- gsub("_out.txt", "", mhcnuggets_random$Allele)
colnames(mhcflurry_random) <- c("Peptide", "Allele", "mhcflurry_affinity")
netmhcpan_random_alleles <- unique(netmhcpan_random$Allele)
mhcflurry_random <- mhcflurry_random[mhcflurry_random$Allele %in% netmhcpan_random_alleles,]
mhcnuggets_random <- mhcnuggets_random[mhcnuggets_random$Allele %in% netmhcpan_random_alleles,]
merged_peps_random_all <- merge(mhcflurry_random, mhcnuggets_random, by=c("Peptide", "Allele"))
merged_peps_random_all <- merge(merged_peps_random_all, netmhcpan_random, by=c("Peptide", "Allele"))
colnames(merged_peps_random_all)[5] <- "netmhcpan_affinity"
hlathena_random <- read.csv("hlathena_joint_alleles_random.csv")
colnames(hlathena_random) <- c("Peptide", "Allele", "hlathena_affinity")
merged_peps_random_all <- merge(merged_peps_random_all, hlathena_random, by=c("Peptide", "Allele"))
write.csv(merged_peps_random_all, "merged_peps_random_all.csv", row.names = FALSE)
merged_peps_random_all <- read.csv("merged_peps_random_all.csv")

#Fig 1a
###ICC per peptide rather than per allele
kappa_stat_pep <- data.frame(peptide=character(),value=numeric())
i <- 1
list_peps <- unique(merged_peps_random_all$Peptide)
list_peps <- as.data.frame(list_peps)
for (i in 1:dim(list_peps)[1])
{
  bigtemp <- merged_peps_random_all[merged_peps_random_all$Peptide == list_peps[i,1],3:6]
  kappa_stat_pep[nrow(kappa_stat_pep) + 1,] = c(list_peps[i,1], icc(bigtemp, model="twoway", type="agreement")[7])
}
ggplot(kappa_stat_pep, aes(x=ICC)) + geom_histogram(bins=100, color="black", fill="white") + 
  geom_vline(aes(xintercept=mean(ICC)), color="blue", linetype="dashed", size=0.5) +
  theme_bw()+ylab("Peptides (#)")

#Fig 1b
iota_stat  <- data.frame(allele=character(),value=numeric())
library(irr)
list_alleles <- as.data.frame(table(merged_peps_random_all$Allele))
for (i in 1:dim(list_alleles)[1])
{
  bigtemp <- merged_peps_random_all[merged_peps_random_all$Allele == list_alleles[i,1],3:6]
  iota_stat[nrow(iota_stat) + 1,] = c(as.character(list_alleles[i,1]), iota(list(bigtemp),scaledata="quantitative", standardize=FALSE)[5])
}
write.table(iota_stat, "iota_stat.csv", row.names = FALSE, quote=FALSE, sep=",")
iota_stat <- read.csv("iota_stat.csv")
iota_stat$Class <- ifelse(grepl("A", iota_stat$allele), 'A', ifelse(grepl("B", iota_stat$allele), 'B', 'C'))
colnames(iota_stat) <- c("Allele", "ICC", "Class")
library(ggpubr)
gghistogram(iota_stat, x = "ICC",
            add = "mean", rug = TRUE,
            color = "Class", fill = "Class",
            bins=10) + ylab("Alleles (#)")

#Fig 1c
merged_peps_random_all_filt <- merged_peps_random_all[merged_peps_random_all$mhcflurry_affinity >= 0.5,]
merged_peps_random_all_filt <- merged_peps_random_all_filt[merged_peps_random_all_filt$mhcnuggets_affinity >= 0.5,]
merged_peps_random_all_filt <- merged_peps_random_all_filt[merged_peps_random_all_filt$netmhcpan_affinity >= 0.5,]
merged_peps_random_all_filt <- merged_peps_random_all_filt[merged_peps_random_all_filt$hlathena_affinity >= 0.5,]

merged_peps_random_all_filt1 <- merged_peps_random_all[merged_peps_random_all$mhcflurry_affinity >= 0.5,]
merged_peps_random_all_filt2 <- merged_peps_random_all[merged_peps_random_all$mhcnuggets_affinity >= 0.5,]
merged_peps_random_all_filt3 <- merged_peps_random_all[merged_peps_random_all$netmhcpan_affinity >= 0.5,]
merged_peps_random_all_filt4 <- merged_peps_random_all[merged_peps_random_all$hlathena_affinity >= 0.5,]

temp <- rbind(merged_peps_random_all_filt1, merged_peps_random_all_filt2, merged_peps_random_all_filt3, merged_peps_random_all_filt4)
any_bind_table <- table(temp$Allele)
all_bind_table <- table(merged_peps_random_all_filt$Allele)
comb_table <- merge(as.data.frame(any_bind_table), as.data.frame(all_bind_table), by="Var1",all=TRUE )
comb_table[is.na(comb_table)] <- 0

ggplot(comb_table, aes(x=Freq.x, y=Freq.y)) + 
  geom_point() + geom_smooth(method="lm", formula=y~x+0) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Unique Peptide Binding >=1 Tool") + ylab("Unique Peptide Binding 4 Tools") + ggtitle("HLAthena Correlation vs Peptide Difference")

##upset
temp[temp$mhcflurry_affinity >= 0.5,]$mhcflurry_affinity <- 1
temp[temp$mhcflurry_affinity < 0.5,]$mhcflurry_affinity <- 0
temp[temp$mhcnuggets_affinity >= 0.5,]$mhcnuggets_affinity <- 1
temp[temp$mhcnuggets_affinity < 0.5,]$mhcnuggets_affinity <- 0
temp[temp$netmhcpan_affinity >= 0.5,]$netmhcpan_affinity <- 1
temp[temp$netmhcpan_affinity < 0.5,]$netmhcpan_affinity <- 0
temp[temp$hlathena_affinity >= 0.5,]$hlathena_affinity <- 1
temp[temp$hlathena_affinity < 0.5,]$hlathena_affinity <- 0

require(ggupset)
library(dplyr)
temp$Peptide <- paste(temp$Peptide, temp$Allele, sep="-")
temp <- distinct(temp)
temp$mhcflurry_affinity <- as.logical(temp$mhcflurry_affinity)
temp$mhcnuggets_affinity <- as.logical(temp$mhcnuggets_affinity)
temp$netmhcpan_affinity <- as.logical(temp$netmhcpan_affinity)
temp$hlathena_affinity <- as.logical(temp$hlathena_affinity)

temp2 <- temp[,-1]
rownames(temp2) <- temp[,1]
temp2 <- temp2[,c(2:5)]
#temp2 <- head(temp2,100)
temp2 <- t(temp2)
temp_tibble <- temp2 %>%
  as_tibble(rownames = "Tool") %>%
  gather(Peptide, Member, -Tool) %>%
  filter(Member) %>%
  select(-Member)
#temp_tibble %>%
#  group_by(Peptide) %>%
#  summarize(Tool = list(Tool))

options(scipen=999)
temp_tibble$Tool <- gsub("_affinity", "", temp_tibble$Tool)
temp_tibble$Tool <- gsub("mhcflurry", "MHCflurry", temp_tibble$Tool)
temp_tibble$Tool <- gsub("netmhcpan", "netMHCpan", temp_tibble$Tool)
temp_tibble$Tool <- gsub("mhcnuggets", "MHCnuggets", temp_tibble$Tool)
temp_tibble$Tool <- gsub("hlathena", "HLAthena", temp_tibble$Tool)

write.csv(temp_tibble, "upset_tibble.csv", quote=FALSE, row.names=FALSE)
temp_tibble <- read.csv("upset_tibble.csv")
options(scipen=999)
temp_tibble %>%
  group_by(Peptide) %>%
  summarize(Tool = list(Tool)) %>%
  ggplot(aes(x = Tool)) + xlab("")+
  geom_bar() + ylab("Predicted Peptide Bound (#)") + 
  scale_x_upset(order_by="degree")



##Figure 2


hlathena_alleles <- as.data.frame(table(merged_peps_random_all[merged_peps_random_all$hlathena_affinity > 0.5,]$Allele))
netmhcpan_alleles <- as.data.frame(table(merged_peps_random_all[merged_peps_random_all$netmhcpan_affinity > 0.5,]$Allele))
mhcflurry_alleles <- as.data.frame(table(merged_peps_random_all[merged_peps_random_all$mhcflurry_affinity > 0.5,]$Allele))
mhcnuggets_alleles <- as.data.frame(table(merged_peps_random_all[merged_peps_random_all$mhcnuggets_affinity > 0.5,]$Allele))

#hlathena_alleles <- as.data.frame(table(merged_peps_random_all[merged_peps_random_all$hlathena_affinity > 0.6,]$Allele))
#netmhcpan_alleles <- as.data.frame(table(merged_peps_random_all[merged_peps_random_all$netmhcpan_affinity > 0.6,]$Allele))
#mhcflurry_alleles <- as.data.frame(table(merged_peps_random_all[merged_peps_random_all$mhcflurry_affinity > 0.6,]$Allele))
#mhcnuggets_alleles <- as.data.frame(table(merged_peps_random_all[merged_peps_random_all$mhcnuggets_affinity > 0.6,]$Allele))

#hlathena_alleles <- as.data.frame(table(merged_peps_random_all[merged_peps_random_all$hlathena_affinity > 0.7,]$Allele))
#netmhcpan_alleles <- as.data.frame(table(merged_peps_random_all[merged_peps_random_all$netmhcpan_affinity > 0.7,]$Allele))
#mhcflurry_alleles <- as.data.frame(table(merged_peps_random_all[merged_peps_random_all$mhcflurry_affinity > 0.7,]$Allele))
#mhcnuggets_alleles <- as.data.frame(table(merged_peps_random_all[merged_peps_random_all$mhcnuggets_affinity > 0.7,]$Allele))

merge_alleles <- merge(hlathena_alleles, netmhcpan_alleles, by="Var1", all=TRUE)
merge_alleles <- merge(merge_alleles, mhcflurry_alleles, by="Var1", all=TRUE)
merge_alleles <- merge(merge_alleles, mhcnuggets_alleles, by="Var1", all=TRUE)
colnames(merge_alleles) <- c("Allele", "HLAthena", "netMHCpan", "MHCflurry", "MHCnuggets")
write.csv(merge_alleles, "merge_alleles_random_0.5.csv", quote=FALSE, row.names=FALSE)
#write.csv(merge_alleles, "merge_alleles_random_0.6.csv", quote=FALSE, row.names=FALSE)
#write.csv(merge_alleles, "merge_alleles_random_0.7.csv", quote=FALSE, row.names=FALSE)

merge_alleles <- read.csv("merge_alleles_random_0.5.csv")
merge_alleles[is.na(merge_alleles)] <- 0

#make sure to source the mod_pairs and mod_pairs_panel.R files
pairs.panels(merge_alleles[,-1], 
             method = "spearman", # correlation method
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)


#Figure 3
###HLAthena training data vs metrics
path <- "NIHMS1541512-supplement-Sup_Data1.xlsx"
library(readxl)
hlathena_data <- lapply(excel_sheets(path), read_excel, path = path)
hlathena_data <- hlathena_data[-1]

hlathena_seqs <- matrix(data=NA, nrow=0, ncol=2)
for (i in 1:length(hlathena_data))
{
  test <- as.data.frame(hlathena_data[i])
  tempdf <- cbind(test$sequence,substr(test$directory,1,5))
  hlathena_seqs <- rbind(hlathena_seqs, tempdf)
}
hlathena_seqs <- as.data.frame(hlathena_seqs)
hlathena_training_table <- table(hlathena_seqs$V2)
hlathena_training_table <- as.data.frame(hlathena_training_table)
colnames(hlathena_training_table)[1] <- "Var1"
hlathena_random_table <- table(merged_peps_random_all[merged_peps_random_all$hlathena_affinity >= 0.5,]$Allele)
hlathena_random_train_merge <- merge(hlathena_training_table, hlathena_random_table, by="Var1")


###MHCflurry training data vs random
mhcflurry_train <- read.csv("allele_train_list2.csv")
mhcflurry_train_table <- as.data.frame(table(mhcflurry_train))
rm(mhcflurry_train)
mhcflurry_train_table$mhcflurry_train <- gsub("HLA-", "", mhcflurry_train_table$mhcflurry_train)
mhcflurry_train_table$mhcflurry_train <- gsub("[*]", "", mhcflurry_train_table$mhcflurry_train)
mhcflurry_train_table$mhcflurry_train <- gsub(":", "", mhcflurry_train_table$mhcflurry_train)
colnames(mhcflurry_train_table)[1] <- "Var1"

mhcflurry_random_table <- table(merged_peps_random_all[merged_peps_random_all$mhcflurry_affinity >= 0.5,]$Allele)
mhcflurry_random_train_merge <- merge(mhcflurry_random_table, mhcflurry_train_table, by="Var1")

##MHCnuggets training data vs random
mhcnuggets_train <- read.csv("curated_training_data.csv")
mhcnuggets_train_table <- as.data.frame(table(mhcnuggets_train$mhc))
rm(mhcnuggets_train)
mhcnuggets_train_table <- mhcnuggets_train_table[grep("HLA", mhcnuggets_train_table$Var1),]
mhcnuggets_train_table$Var1 <- gsub("HLA-", "", mhcnuggets_train_table$Var1)
mhcnuggets_train_table$Var1 <- gsub("[*]", "", mhcnuggets_train_table$Var1)
mhcnuggets_train_table$Var1 <- gsub(":", "", mhcnuggets_train_table$Var1)

mhcnuggets_random_table <- table(merged_peps_random_all[merged_peps_random_all$mhcnuggets_affinity >= 0.5,]$Allele)
head(mhcnuggets_random_table)
mhcnuggets_random_table$Var1 <- gsub(":", "", mhcnuggets_random_table$Var1)
mhcnuggets_random_train_merge <- merge(mhcnuggets_train_table, mhcnuggets_random_table, by="Var1")

##netmhcpan training
netmhcpan_train <- read.table("netmhcpan_training_data", header=FALSE, sep=",")
netmhcpan_train$V1 <- gsub("HLA-", "", netmhcpan_train$V1)
netmhcpan_train$V1 <- gsub(":", "", netmhcpan_train$V1)

netmhcpan_random_bind <- as.data.frame(table(netmhcpan_random_all_alleles_filt$Allele))
rm(netmhcpan_random_all_alleles_filt)
colnames(netmhcpan_train)[1] <- "Var1"

netmhcpan_random_table <- table(merged_peps_random_all[merged_peps_random_all$netmhcpan_affinity >= 0.5,]$Allele)

netmhcpan_random_train_merge <- merge(netmhcpan_train, netmhcpan_random_table, by="Var1")


###All training data
hlathena_random_train_merge$Tool <- "HLAthena"
colnames(hlathena_random_train_merge) <- c("Allele", "Train", "Random Bind", "Tool")
netmhcpan_random_train_merge$Tool <- "netMHCpan"
colnames(netmhcpan_random_train_merge) <- c("Allele", "Train", "Random Bind", "Tool")
mhcflurry_random_train_merge$Tool <- "MHCflurry"
colnames(mhcflurry_random_train_merge) <- c("Allele", "Train", "Random Bind", "Tool")
mhcnuggets_random_train_merge$Tool <- "MHCnuggets"
colnames(mhcnuggets_random_train_merge) <- c("Allele", "Train", "Random Bind", "Tool")


all_random_train_merge <- rbind(hlathena_random_train_merge, netmhcpan_random_train_merge, 
                                mhcflurry_random_train_merge, mhcnuggets_random_train_merge)

library(ggplot2)

merge <- merge(all_train, iota_stat[,1:2], by="Allele")
ggplot(data=merge, aes(x=mean, y=`ICC`)) + geom_point() + xlab("Training Data (Mean # Peptides)") + ylab("ICC")


all_random_train_merge$`Random Bind` <- all_random_train_merge$`Random Bind`*100/800000
ggplot(data=all_random_train_merge, aes(x=Train, y=`Random Bind`, color=Tool, shape=Tool)) +
  geom_point() + xlab("Training Data (# Peptides)") + ylab("Predicted Peptides Bound (%)") 

#Figure 4/Sup 4-5

##MHCnuggets
mhcnuggets_betaherpes_6a <- read.csv("mhcnuggets_betaherpes_6a.csv", header=FALSE)
colnames(mhcnuggets_betaherpes_6a) <- c("Peptide", "mhcnuggets_affinity", "allele")
mhcnuggets_betaherpes_6a <- distinct(mhcnuggets_betaherpes_6a)
mhcnuggets_betaherpes_6a$mhcnuggets_affinity <- 1-log(mhcnuggets_betaherpes_6a$mhcnuggets_affinity)/log(50000)
mhcnuggets_betaherpes_6a$allele <- gsub("HLA-", "", mhcnuggets_betaherpes_6a$allele)
mhcnuggets_betaherpes_6a$allele <- gsub("_out.txt", "", mhcnuggets_betaherpes_6a$allele)

mhcnuggets_bk <- read.csv("mhcnuggets_bk.csv", header=FALSE)
colnames(mhcnuggets_bk) <- c("Peptide", "mhcnuggets_affinity", "allele")
mhcnuggets_bk <- distinct(mhcnuggets_bk)
mhcnuggets_bk$mhcnuggets_affinity <- 1-log(mhcnuggets_bk$mhcnuggets_affinity)/log(50000)
mhcnuggets_bk$allele <- gsub("HLA-", "", mhcnuggets_bk$allele)
mhcnuggets_bk$allele <- gsub("_out.txt", "", mhcnuggets_bk$allele)

mhcnuggets_cmv <- read.csv("mhcnuggets_cmv.csv", header=FALSE)
colnames(mhcnuggets_cmv) <- c("Peptide", "mhcnuggets_affinity", "allele")
mhcnuggets_cmv <- distinct(mhcnuggets_cmv)
mhcnuggets_cmv$mhcnuggets_affinity <- 1-log(mhcnuggets_cmv$mhcnuggets_affinity)/log(50000)
mhcnuggets_cmv$allele <- gsub("HLA-", "", mhcnuggets_cmv$allele)
mhcnuggets_cmv$allele <- gsub("_out.txt", "", mhcnuggets_cmv$allele)

mhcnuggets_covid <- read.csv("mhcnuggets_covid.csv", header=FALSE)
colnames(mhcnuggets_covid) <- c("Peptide", "mhcnuggets_affinity", "allele")
mhcnuggets_covid <- distinct(mhcnuggets_covid)
mhcnuggets_covid$mhcnuggets_affinity <- 1-log(mhcnuggets_covid$mhcnuggets_affinity)/log(50000)
mhcnuggets_covid$allele <- gsub("HLA-", "", mhcnuggets_covid$allele)
mhcnuggets_covid$allele <- gsub("_out.txt", "", mhcnuggets_covid$allele)

mhcnuggets_ebv <- read.csv("mhcnuggets_ebv.csv", header=FALSE)
colnames(mhcnuggets_ebv) <- c("Peptide", "mhcnuggets_affinity", "allele")
mhcnuggets_ebv <- distinct(mhcnuggets_ebv)
mhcnuggets_ebv$mhcnuggets_affinity <- 1-log(mhcnuggets_ebv$mhcnuggets_affinity)/log(50000)
mhcnuggets_ebv$allele <- gsub("HLA-", "", mhcnuggets_ebv$allele)
mhcnuggets_ebv$allele <- gsub("_out.txt", "", mhcnuggets_ebv$allele)

mhcnuggets_herpes_1_17 <- read.csv("mhcnuggets_herpes_1_17.csv", header=FALSE)
colnames(mhcnuggets_herpes_1_17) <- c("Peptide", "mhcnuggets_affinity", "allele")
mhcnuggets_herpes_1_17 <- distinct(mhcnuggets_herpes_1_17)
mhcnuggets_herpes_1_17$mhcnuggets_affinity <- 1-log(mhcnuggets_herpes_1_17$mhcnuggets_affinity)/log(50000)
mhcnuggets_herpes_1_17$allele <- gsub("HLA-", "", mhcnuggets_herpes_1_17$allele)
mhcnuggets_herpes_1_17$allele <- gsub("_out.txt", "", mhcnuggets_herpes_1_17$allele)

mhcnuggets_herpes_2_hg52 <- read.csv("mhcnuggets_herpes_2_hg52.csv", header=FALSE)
colnames(mhcnuggets_herpes_2_hg52) <- c("Peptide", "mhcnuggets_affinity", "allele")
mhcnuggets_herpes_2_hg52 <- distinct(mhcnuggets_herpes_2_hg52)
mhcnuggets_herpes_2_hg52$mhcnuggets_affinity <- 1-log(mhcnuggets_herpes_2_hg52$mhcnuggets_affinity)/log(50000)
mhcnuggets_herpes_2_hg52$allele <- gsub("HLA-", "", mhcnuggets_herpes_2_hg52$allele)
mhcnuggets_herpes_2_hg52$allele <- gsub("_out.txt", "", mhcnuggets_herpes_2_hg52$allele)

mhcnuggets_human <- read.csv("mhcnuggets_human_filt.csv", header=FALSE)
colnames(mhcnuggets_human) <- c("Peptide", "mhcnuggets_affinity", "allele")
mhcnuggets_human <- distinct(mhcnuggets_human)
mhcnuggets_human$mhcnuggets_affinity <- 1-log(mhcnuggets_human$mhcnuggets_affinity)/log(50000)
mhcnuggets_human$allele <- gsub("HLA-", "", mhcnuggets_human$allele)
mhcnuggets_human$allele <- gsub("_out.txt", "", mhcnuggets_human$allele)

mhcnuggets_random <- read.csv("mhcnuggets_random.csv", header=FALSE)
colnames(mhcnuggets_random) <- c("Peptide", "mhcnuggets_affinity", "allele")
mhcnuggets_random <- distinct(mhcnuggets_random)
mhcnuggets_random$mhcnuggets_affinity <- 1-log(mhcnuggets_random$mhcnuggets_affinity)/log(50000)
mhcnuggets_random$allele <- gsub("HLA-", "", mhcnuggets_random$allele)
mhcnuggets_random$allele <- gsub("_out.txt", "", mhcnuggets_random$allele)


mhcnuggets_betaherpes_6a_table <- as.data.frame(table(mhcnuggets_betaherpes_6a[mhcnuggets_betaherpes_6a$mhcnuggets_affinity >= 0.5,]$allele))
mhcnuggets_bk_table <- as.data.frame(table(mhcnuggets_bk[mhcnuggets_bk$mhcnuggets_affinity >= 0.5,]$allele))
mhcnuggets_cmv_table <- as.data.frame(table(mhcnuggets_cmv[mhcnuggets_cmv$mhcnuggets_affinity >= 0.5,]$allele))
mhcnuggets_covid_table <- as.data.frame(table(mhcnuggets_covid[mhcnuggets_covid$mhcnuggets_affinity >= 0.5,]$allele))
mhcnuggets_ebv_table <- as.data.frame(table(mhcnuggets_ebv[mhcnuggets_ebv$mhcnuggets_affinity >= 0.5,]$allele))
mhcnuggets_herpes_1_17_table <- as.data.frame(table(mhcnuggets_herpes_1_17[mhcnuggets_herpes_1_17$mhcnuggets_affinity >= 0.5,]$allele))
mhcnuggets_herpes_2_hg52_table <- as.data.frame(table(mhcnuggets_herpes_2_hg52[mhcnuggets_herpes_2_hg52$mhcnuggets_affinity >= 0.5,]$allele))
mhcnuggets_human_table <- as.data.frame(table(mhcnuggets_human[mhcnuggets_human$mhcnuggets_affinity >= 0.5,]$allele))
mhcnuggets_random_table <- as.data.frame(table(mhcnuggets_random[mhcnuggets_random$mhcnuggets_affinity >= 0.5,]$allele))

rm(mhcnuggets_betaherpes_6a)
rm(mhcnuggets_bk)
rm(mhcnuggets_cmv)
rm(mhcnuggets_covid)
rm(mhcnuggets_ebv)
rm(mhcnuggets_herpes_1_17)
rm(mhcnuggets_herpes_2_hg52)
rm(mhcnuggets_human)
rm(mhcnuggets_random)

merge_mhcnuggets <- merge(mhcnuggets_betaherpes_6a_table, mhcnuggets_bk_table, by="Var1", all=TRUE)
merge_mhcnuggets <- merge(merge_mhcnuggets, mhcnuggets_cmv_table, by="Var1", all=TRUE)
merge_mhcnuggets <- merge(merge_mhcnuggets, mhcnuggets_covid_table, by="Var1", all=TRUE)
merge_mhcnuggets <- merge(merge_mhcnuggets, mhcnuggets_ebv_table, by="Var1", all=TRUE)
merge_mhcnuggets <- merge(merge_mhcnuggets, mhcnuggets_herpes_1_17_table, by="Var1", all=TRUE)
merge_mhcnuggets <- merge(merge_mhcnuggets, mhcnuggets_herpes_2_hg52_table, by="Var1", all=TRUE)
merge_mhcnuggets <- merge(merge_mhcnuggets, mhcnuggets_human_table, by="Var1", all=TRUE)
merge_mhcnuggets <- merge(merge_mhcnuggets, mhcnuggets_random_table, by="Var1", all=TRUE)

#merged_mhcnuggets_pep
#big_table <- as.data.frame(table(merge_mhcnuggets_peps$Allele, merge_mhcnuggets_peps$Source))
#big_table <- big_table %>% spread(Var2, Freq)
#big_table <- big_table[big_table$Var1 %in% list_alleles[,1],]
#big_table_cor <- cor(big_table[,2:10], method="spearman")


head(merge_mhcnuggets)
colnames(merge_mhcnuggets) <- c("Allele", "Betaherpes 6a", "BK", "CMV", "Covid", "EBV", "Herpes 1", "Herpes 2", "Human", "Random")
merge_mhcnuggets[is.na(merge_mhcnuggets)] <- 0
merge_mhcnuggets$Viral <- rowSums(merge_mhcnuggets[,2:8])
#list_alleles$Allele <- gsub('^(.{3})(.*)$', '\\1:\\2', list_alleles$Allele)

write.csv(merge_mhcnuggets, "merge_mhcnuggets.csv", quote=FALSE, row.names=FALSE)
merge_mhcnuggets <- read.csv("merge_mhcnuggets.csv")
merge_mhcnuggets$Allele <- gsub(":", "", merge_mhcnuggets$Allele)
merge_mhcnuggets <- merge_mhcnuggets[merge_mhcnuggets$Allele %in% list_alleles$Allele,]

#cor(merge_mhcnuggets[,9:11])
colnames(merge_mhcnuggets) <- mapvalues(colnames(merge_mhcnuggets), map2$Name, map2$Map)
merge_mhcnuggets <- merge_mhcnuggets %>% select(order(colnames(merge_mhcnuggets)))

cormat <- cor(merge_mhcnuggets[2:10], method="spearman", use="complete.obs")
temp <- merge_mhcnuggets[,c(9,11,10)]
#list_all_2000 <- gsub('^(.{3})(.*)$', '\\1:\\2', list_all_2000)
cormat2 <- cor(temp,  method="spearman", use="complete.obs")


get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
upper_tri

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ geom_text(aes(label = round(value,2)), color = "white", size = 4) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ggtitle("MHCnuggets")


#MHCflurry
mhcflurry_betaherpes_6a <- read.csv("mhcflurry_betaherpes_6a_preds.csv")
mhcflurry_betaherpes_6a$mhcflurry_affinity <- 1-log(mhcflurry_betaherpes_6a$mhcflurry_affinity)/log(50000)
mhcflurry_betaherpes_6a <- distinct(mhcflurry_betaherpes_6a)

mhcflurry_bk <- read.csv("mhcflurry_bk_preds.csv")
mhcflurry_bk$mhcflurry_affinity <- 1-log(mhcflurry_bk$mhcflurry_affinity)/log(50000)
mhcflurry_bk <- distinct(mhcflurry_bk)

mhcflurry_cmv <- read.csv("mhcflurry_cmv_filt.csv")
mhcflurry_cmv$mhcflurry_affinity <- 1-log(mhcflurry_cmv$mhcflurry_affinity)/log(50000)
mhcflurry_cmv <- distinct(mhcflurry_cmv)

mhcflurry_covid <- read.csv("mhcflurry_covid_preds.csv")
mhcflurry_covid$mhcflurry_affinity <- 1-log(mhcflurry_covid$mhcflurry_affinity)/log(50000)
mhcflurry_covid <- distinct(mhcflurry_covid)

mhcflurry_ebv <- read.csv("mhcflurry_ebv_preds.csv")
mhcflurry_ebv$mhcflurry_affinity <- 1-log(mhcflurry_ebv$mhcflurry_affinity)/log(50000)
mhcflurry_ebv <- distinct(mhcflurry_ebv)

mhcflurry_herpes_1_17 <- read.csv("mhcflurry_herpes_1_17_preds.csv")
mhcflurry_herpes_1_17$mhcflurry_affinity <- 1-log(mhcflurry_herpes_1_17$mhcflurry_affinity)/log(50000)
mhcflurry_herpes_1_17 <- distinct(mhcflurry_herpes_1_17)

mhcflurry_herpes_2_hg52 <- read.csv("mhcflurry_herpes_2_hg52_preds.csv")
mhcflurry_herpes_2_hg52$mhcflurry_affinity <- 1-log(mhcflurry_herpes_2_hg52$mhcflurry_affinity)/log(50000)
mhcflurry_herpes_2_hg52 <- distinct(mhcflurry_herpes_2_hg52)


mhcflurry_human <- read.csv("mhcflurry_human_filt.csv", header=FALSE)
colnames(mhcflurry_human) <- c("peptide", "allele", "mhcflurry_affinity")
mhcflurry_human$mhcflurry_affinity <- 1-log(as.numeric(mhcflurry_human$mhcflurry_affinity))/log(50000)
mhcflurry_human <- distinct(mhcflurry_human)


mhcflurry_random <- read.csv("mhcflurry_random_filt.csv")
mhcflurry_random$mhcflurry_affinity <- 1-log(mhcflurry_random$mhcflurry_affinity)/log(50000)
mhcflurry_random <- distinct(mhcflurry_random)

mhcflurry_betaherpes_6a_table <- as.data.frame(table(mhcflurry_betaherpes_6a[mhcflurry_betaherpes_6a$mhcflurry_affinity >= 0.5,]$allele))
mhcflurry_bk_table <- as.data.frame(table(mhcflurry_bk[mhcflurry_bk$mhcflurry_affinity >= 0.5,]$allele))
mhcflurry_cmv_table <- as.data.frame(table(mhcflurry_cmv[mhcflurry_cmv$mhcflurry_affinity >= 0.5,]$allele))
mhcflurry_covid_table <- as.data.frame(table(mhcflurry_covid[mhcflurry_covid$mhcflurry_affinity >= 0.5,]$allele))
mhcflurry_ebv_table <- as.data.frame(table(mhcflurry_ebv[mhcflurry_ebv$mhcflurry_affinity >= 0.5,]$allele))
mhcflurry_herpes_1_17_table <- as.data.frame(table(mhcflurry_herpes_1_17[mhcflurry_herpes_1_17$mhcflurry_affinity >= 0.5,]$allele))
mhcflurry_herpes_2_hg52_table <- as.data.frame(table(mhcflurry_herpes_2_hg52[mhcflurry_herpes_2_hg52$mhcflurry_affinity >= 0.5,]$allele))
mhcflurry_human_table <- as.data.frame(table(mhcflurry_human[mhcflurry_human$mhcflurry_affinity >= 0.5,]$allele))
mhcflurry_random_table <- as.data.frame(table(mhcflurry_random[mhcflurry_random$mhcflurry_affinity >= 0.5,]$allele))

rm(mhcflurry_betaherpes_6a)
rm(mhcflurry_bk)
rm(mhcflurry_cmv)
rm(mhcflurry_covid)
rm(mhcflurry_ebv)
rm(mhcflurry_herpes_1_17)
rm(mhcflurry_herpes_2_hg52)
rm(mhcflurry_human)
rm(mhcflurry_random)

merge_mhcflurry <- merge(mhcflurry_betaherpes_6a_table, mhcflurry_bk_table, by="Var1", all=TRUE)
merge_mhcflurry <- merge(merge_mhcflurry, mhcflurry_cmv_table, by="Var1", all=TRUE)
merge_mhcflurry <- merge(merge_mhcflurry, mhcflurry_covid_table, by="Var1", all=TRUE)
merge_mhcflurry <- merge(merge_mhcflurry, mhcflurry_ebv_table, by="Var1", all=TRUE)
merge_mhcflurry <- merge(merge_mhcflurry, mhcflurry_herpes_1_17_table, by="Var1", all=TRUE)
merge_mhcflurry <- merge(merge_mhcflurry, mhcflurry_herpes_2_hg52_table, by="Var1", all=TRUE)
merge_mhcflurry <- merge(merge_mhcflurry, mhcflurry_human_table, by="Var1", all=TRUE)
merge_mhcflurry <- merge(merge_mhcflurry, mhcflurry_random_table, by="Var1", all=TRUE)

head(merge_mhcflurry)
colnames(merge_mhcflurry) <- c("Allele", "Betaherpes 6a", "BK", "CMV", "Covid", "EBV", "Herpes 1", "Herpes 2", "Human", "Random")
merge_mhcflurry[is.na(merge_mhcflurry)] <- 0
merge_mhcflurry$Viral <- rowSums(merge_mhcflurry[,2:8])
cor(merge_mhcflurry[,9:11])

merge_mhcflurry <- merge_mhcflurry[merge_mhcflurry$Allele %in% list_alleles[,1],]
write.csv(merge_mhcflurry, "merge_mhcflurry.csv", quote=FALSE, row.names=FALSE)
merge_mhcflurry <- read.csv("merge_mhcflurry.csv")

map2 <- as.data.frame(colnames(merge_mhcflurry))
map2$Map <- c("Allele", "HHV-6", "BK", "HHV-5", "Covid", "HSV-4", "HSV-1", "HSV-2", "Human", "Random", "Viral")
colnames(map2)[1] <- "Name"

colnames(merge_mhcflurry) <- mapvalues(colnames(merge_mhcflurry), map2$Name, map2$Map)
merge_mhcflurry <- merge_mhcflurry %>% select(order(colnames(merge_mhcflurry)))

temp <- merge_mhcflurry[,c(9,11,10)]
#merge_mhcflurry_filt <- merge_mhcflurry[merge_mhcflurry$Allele %in% list_all_2000,]
cormat <- cor(merge_mhcflurry[2:10], method="spearman", use="complete.obs")
cormat2 <- cor(temp, method="spearman")
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat2)
upper_tri

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ geom_text(aes(label = round(value,2)), color = "white", size = 4) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ggtitle("MHCflurry")

merge_mhcflurry <- merge_mhcflurry[merge_mhcflurry$Allele %in% bighlathena_list,]
cormat <- cor(merge_mhcflurry[2:10], method="spearman", use="complete.obs")
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
upper_tri

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 5, hjust = 1))+
  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())


#netMHCpan
###Netmhcpan all peptides filtered all alleles
betaherpes_6a_nmc_filt <- read.csv("netmhcpan_betaherpes_6a_filt.csv")
betaherpes_6a_nmc_filt <- distinct(betaherpes_6a_nmc_filt)

bk_nmc_filt <- read.csv("netmhcpan_bk_filt.csv")
bk_nmc_filt <- distinct(bk_nmc_filt)

cmv_nmc_filt <- read.csv("netmhcpan_cmv_filt.csv")
cmv_nmc_filt <- distinct(cmv_nmc_filt)

covid_nmc_filt <- read.csv("netmhcpan_covid_filt.csv")
covid_nmc_filt <- distinct(covid_nmc_filt)

ebv_nmc_filt <- read.csv("netmhcpan_ebv_filt.csv")
ebv_nmc_filt <- distinct(ebv_nmc_filt)

herpes_1_17_nmc_filt <- read.csv("netmhcpan_herpes_1_17_filt.csv")
herpes_1_17_nmc_filt <- distinct(herpes_1_17_nmc_filt)

herpes_2_hg52_nmc_filt <- read.csv("netmhcpan_herpes_2_hg52_filt.csv")
herpes_2_hg52_nmc_filt <- distinct(herpes_2_hg52_nmc_filt)

netmhcpan_betaherpes_6a_table <- as.data.frame(table(betaherpes_6a_nmc_filt$Allele))
netmhcpan_bk_table <- as.data.frame(table(bk_nmc_filt$Allele))
netmhcpan_cmv_table <- as.data.frame(table(cmv_nmc_filt$Allele))
netmhcpan_covid_table <- as.data.frame(table(covid_nmc_filt$Allele))
netmhcpan_ebv_table <- as.data.frame(table(ebv_nmc_filt$Allele))
netmhcpan_herpes1_table <- as.data.frame(table(herpes_1_17_nmc_filt$Allele))
netmhcpan_herpes2_table <- as.data.frame(table(herpes_2_hg52_nmc_filt$Allele))
netmhcpan_human_table

nmc_alleles <- intersect(gsub("HLA-", "", gsub(":", "", netmhcpan_human_table$Var1)), list_alleles$Allele)

netmhcpan_group_table <- merge(netmhcpan_betaherpes_6a_table, netmhcpan_bk_table, by="Var1", all=TRUE)
netmhcpan_group_table <- merge(netmhcpan_group_table, netmhcpan_cmv_table, by="Var1", all=TRUE)
netmhcpan_group_table <- merge(netmhcpan_group_table, netmhcpan_covid_table, by="Var1", all=TRUE)
netmhcpan_group_table <- merge(netmhcpan_group_table, netmhcpan_ebv_table, by="Var1", all=TRUE)
netmhcpan_group_table <- merge(netmhcpan_group_table, netmhcpan_herpes1_table, by="Var1", all=TRUE)
netmhcpan_group_table <- merge(netmhcpan_group_table, netmhcpan_herpes2_table, by="Var1", all=TRUE)
netmhcpan_group_table <- merge(netmhcpan_group_table, netmhcpan_human_table, by="Var1", all=TRUE)

netmhcpan_group_table$Var1 <- gsub("HLA-", "", gsub(":", "", netmhcpan_group_table$Var1))
netmhcpan_group_table <- netmhcpan_group_table[netmhcpan_group_table$Var1 %in% list_alleles$Allele,]
colnames(netmhcpan_group_table) <- c("Allele", "Betaherpes 6a", "BK", "CMV", "Covid", "EBV", "Herpes 1", "Herpes 2", "Human")


nmc_alleles <- unique(ebv_nmc_filt$Allele)
nmc_alleles <- gsub("HLA-", "", nmc_alleles)
nmc_alleles <- gsub(":", "", nmc_alleles)
length(intersect(nmc_alleles, list_alleles$Allele))


netmhcpan_random_table <- table(netmhcpan_random[netmhcpan_random$Binding_affinity >= 0.5,]$Allele)
netmhcpan_random_table <- as.data.frame(netmhcpan_random_table)
netmhcpan_random_table <- netmhcpan_random_table[netmhcpan_random_table$Var1 %in% list_alleles$Allele,]
colnames(netmhcpan_random_table) <- c("Allele", "Random")
merge_netmhcpan <- merge(netmhcpan_group_table, netmhcpan_random_table, by="Allele")
colnames(merge_netmhcpan) <- mapvalues(colnames(merge_netmhcpan), map2$Name, map2$Map)
merge_netmhcpan <- merge_netmhcpan %>% select(order(colnames(merge_netmhcpan)))
write.csv(merge_netmhcpan, "merge_netmhcpan.csv", quote=FALSE, row.names=FALSE)
merge_netmhcpan <- read.csv("merge_netmhcpan.csv")


#temp <- merge_netmhcpan[,c(9,11,10)]
#colnames(temp) <- c("Human", "Viral", "Random")
cormat <- cor(merge_netmhcpan[merge_netmhcpan$Allele %in% list_all_2000,2:10], method="spearman", use="complete.obs")
#cormat <- cor(temp, method="spearman", use="complete.obs")
#


get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
upper_tri

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+   geom_text(aes(label = round(value,2)), color = "white", size = 4) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + ggtitle("netMHCpan")



#HLAthena
hlathena_betaherpes_6a <- read.csv("betaherpes_6a_hlathena_out.txt")
hlathena_betaherpes_6a_table <- as.data.frame(table(hlathena_betaherpes_6a[hlathena_betaherpes_6a$score >=0.5,]$allele))
rm(hlathena_betaherpes_6a)
gc()

hlathena_bk <- read.csv("bk_hlathena_out.txt")
hlathena_bk_table <- as.data.frame(table(hlathena_bk[hlathena_bk$score >=0.5,]$allele))
rm(hlathena_bk)
gc()

hlathena_cmv <- read.csv("cmv_hlathena_out.txt")
hlathena_cmv_table <- as.data.frame(table(hlathena_cmv[hlathena_cmv$score >=0.5,]$allele))
rm(hlathena_cmv)
gc()

hlathena_covid <- read.csv("covid_hlathena_out.txt")
hlathena_covid_table <- as.data.frame(table(hlathena_covid[hlathena_covid$score >=0.5,]$allele))
rm(hlathena_covid)
gc()

hlathena_ebv <- read.csv("ebv_hlathena_out.txt")
hlathena_ebv_table <- as.data.frame(table(hlathena_ebv[hlathena_ebv$score >=0.5,]$allele))
rm(hlathena_ebv)
gc()


hlathena_herpes1 <- read.csv("herpes_1_17_hlathena_out.txt")
hlathena_herpes1_table <- as.data.frame(table(hlathena_herpes1[hlathena_herpes1$score >=0.5,]$allele))
rm(hlathena_herpes1)
gc()

hlathena_herpes2 <- read.csv("herpes_2_hg52_hlathena_out.txt")
hlathena_herpes2_table <- as.data.frame(table(hlathena_herpes2[hlathena_herpes2$score >=0.5,]$allele))
rm(hlathena_herpes2)
gc()

hlathena_human <- read.csv("human_hlathena_52_filt.txt")
hlathena_human_table <- as.data.frame(table(hlathena_human[hlathena_human$score >=0.5,]$allele))
rm(hlathena_human)
gc()

library(dplyr)
library(tidyr)
hlathena_group_table <- merge(hlathena_betaherpes_6a_table, hlathena_bk_table, by="Var1")
hlathena_group_table <- merge(hlathena_group_table, hlathena_cmv_table, by="Var1")
hlathena_group_table <- merge(hlathena_group_table, hlathena_covid_table, by="Var1")
hlathena_group_table <- merge(hlathena_group_table, hlathena_ebv_table, by="Var1")
hlathena_group_table <- merge(hlathena_group_table, hlathena_herpes1_table, by="Var1")
hlathena_group_table <- merge(hlathena_group_table, hlathena_herpes2_table, by="Var1")
hlathena_group_table <- merge(hlathena_group_table, hlathena_human_table, by="Var1", all = TRUE)
hlathena_group_table <- hlathena_group_table[hlathena_group_table$Var1 %in% list_alleles$Allele,]
colnames(hlathena_group_table) <- c("Allele", "Betaherpes 6a", "BK", "CMV", "Covid", "EBV", "Herpes 1", "Herpes 2", "Human")

merge_hlathena <- merge(hlathena_group_table, hlathena_random_table, by="Allele")
merge_hlathena <- merge_hlathena %>% select(order(colnames(merge_hlathena)))

#merged_table<- merged_table[rownames(merged_table) %in% list_all_2000,]

cormat <- cor(merge_hlathena[2:10], method="spearman", use = "complete.obs")
#temp <- merge_hlathena[,c(8,10,9)]
#cormat <- cor(temp, method="spearman", use = "complete.obs")

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
upper_tri

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
melted_cormat_hlathena <- melted_cormat

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                       midpoint = 0.4, limit = c(-0.2,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal() + geom_text(aes(label = round(value,2)), color = "white", size = 4) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + ggtitle("HLAthena")


##For 2k+
merge_hlathena <- read.csv("merge_hlathena.csv")
merge_netmhcpan <- read.csv("merge_netmhcpan.csv")
merge_mhcflurry <- read.csv("merge_mhcflurry.csv")
merge_mhcnuggets <- read.csv("merge_mhcnuggets.csv")

##Get binned difference metrics
###Binned peps graphs
merged_peps_umap8 <- read.csv("merged_peps_umap8mer.csv")
merged_peps_umap9 <- read.csv("merged_peps_umap9mer.csv")
merged_peps_umap10 <- read.csv("merged_peps_umap10mer.csv")
merged_peps_umap11 <- read.csv("merged_peps_umap11mer.csv")

library(dplyr)
merged_peps_umap8 <- distinct(merged_peps_umap8)
merged_peps_umap9 <- distinct(merged_peps_umap9)
merged_peps_umap10 <- distinct(merged_peps_umap10)
merged_peps_umap11 <- distinct(merged_peps_umap11)

###Binned vs Cor graphs
##Get 8mer cors
library(tidyr)
merged_peps_umap8_filt <- merged_peps_umap8[merged_peps_umap8$Score >= 0.5,]
hlathena_8mers_counts <- as.data.frame(table(merged_peps_umap8_filt$Allele, merged_peps_umap8_filt$Source))
hlathena_8mers_counts <- hlathena_8mers_counts %>% spread(Var2, Freq)
head(hlathena_8mers_counts)
hlathena_8mer_cor <- melt(cor(hlathena_8mers_counts[,2:10]))

##8mer diffs
ab <- matrix( c(min(merged_peps_umap8$X),min(merged_peps_umap8$Y),max(merged_peps_umap8$X),max(merged_peps_umap8$Y)), 2, 2)      
merged_peps_umap8$Source <- mapvalues(merged_peps_umap8$Source, matrix_names[,1], matrix_names[,2])
names <- sort(unique(merged_peps_umap8$Source))
#merged_peps_umap8 <- merged_peps_umap8[,-c(4,5)]
#merged_peps_umap8 <- distinct(merged_peps_umap8)
matrix_bins_8mers <- tibble(x=character(), y=character(), diff=numeric())
library(ash)
for (i in 1:9)
{
  for (j in i:9)
  {
    x <- as.matrix(merged_peps_umap8[merged_peps_umap8$Source == names[i],2:3])  
    nbin <- c( 40, 40)          
    bins <- bin2(x, ab, nbin)
    bin_sol <- bins$nc
    bin_sol1 <- bin_sol/(dim(x)[1])
    
    y <- as.matrix(merged_peps_umap8[merged_peps_umap8$Source == names[j],2:3])  
    #x <- matrix( rnorm(200), 100 , 2)
    bins <- bin2(y, ab, nbin)
    bin_sol <- bins$nc
    bin_sol2 <- bin_sol/(dim(y)[1])
    diff_bins <- bin_sol1 - bin_sol2
    sum_bins <- sum(abs(diff_bins))
    matrix_bins_8mers <- matrix_bins_8mers %>% add_row(x=names[i], y=names[j], diff=sum_bins)
  }
}


matrix_bins_8mers$diff <- matrix_bins_8mers$diff * 100/2
write.csv(matrix_bins_8mers, "matrix_bins_8mers.csv", row.names=FALSE, quote=FALSE)
library(ggplot2)
matrix_bins_8mers <- as.data.frame(matrix_bins_8mers)
ggplot(data = matrix_bins_8mers, aes(y, x, fill = diff))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(0,100), 
                       name="Percent Peptides Difference") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ggtitle("8mer Differential Peptides")

merged_hlathena_8mer <- merge(hlathena_8mer_cor, matrix_bins_8mers, by=c("Var1", "Var2"))
merged_hlathena_8mer <- merged_hlathena_8mer[!merged_hlathena_8mer$Diff==0,]
merged_hlathena_8mer$length <- "8mer"
rm(merged_peps_umap8)
rm(merged_peps_umap8_filt)

merged_peps_umap9 <- read.csv("merged_peps_umap9mer.csv")
###Binned vs Cor graphs
##Get 9mer cors
merged_peps_umap9_filt <- merged_peps_umap9[merged_peps_umap9$Score >= 0.5,]
hlathena_9mers_counts <- as.data.frame(table(merged_peps_umap9_filt$Allele, merged_peps_umap9_filt$Source))
hlathena_9mers_counts <- hlathena_9mers_counts %>% spread(Var2, Freq)
head(hlathena_9mers_counts)
hlathena_9mer_cor <- melt(cor(hlathena_9mers_counts[,2:10]))

##9mer diffs
ab <- matrix( c(min(merged_peps_umap9$X),min(merged_peps_umap9$Y),max(merged_peps_umap9$X),max(merged_peps_umap9$Y)), 2, 2)      
merged_peps_umap9$Source <- mapvalues(merged_peps_umap9$Source, matrix_names[,1], matrix_names[,2])
names <- sort(unique(merged_peps_umap9$Source))
#merged_peps_umap9 <- merged_peps_umap9[,-c(4,5)]
#merged_peps_umap9 <- distinct(merged_peps_umap9)
matrix_bins_9mers <- tibble(x=character(), y=character(), diff=numeric())
library(ash)
for (i in 1:9)
{
  for (j in i:9)
  {
    x <- as.matrix(merged_peps_umap9[merged_peps_umap9$Source == names[i],2:3])  
    nbin <- c( 40, 40)          
    bins <- bin2(x, ab, nbin)
    bin_sol <- bins$nc
    bin_sol1 <- bin_sol/(dim(x)[1])
    
    y <- as.matrix(merged_peps_umap9[merged_peps_umap9$Source == names[j],2:3])  
    #x <- matrix( rnorm(200), 100 , 2)
    bins <- bin2(y, ab, nbin)
    bin_sol <- bins$nc
    bin_sol2 <- bin_sol/(dim(y)[1])
    diff_bins <- bin_sol1 - bin_sol2
    sum_bins <- sum(abs(diff_bins))
    matrix_bins_9mers <- matrix_bins_9mers %>% add_row(x=names[i], y=names[j], diff=sum_bins)
  }
}

matrix_bins_9mers$diff <- matrix_bins_9mers$diff * 100/2
write.csv(matrix_bins_9mers, "matrix_bins_9mers.csv", row.names=FALSE, quote=FALSE)
library(ggplot2)
matrix_bins_9mers <- as.data.frame(matrix_bins_9mers)
ggplot(data = matrix_bins_9mers, aes(y, x, fill = diff))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(0,100), 
                       name="Percent Peptides Difference") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ggtitle("9mer Differential Peptides")


head(hlathena_9mer_cor)
matrix_bins_9mers <- as.data.frame(matrix_bins_9mers)
colnames(matrix_bins_9mers) <- c("Var1", "Var2", "Diff")
head(matrix_bins_9mers)

merged_hlathena_9mer <- merge(hlathena_9mer_cor, matrix_bins_9mers, by=c("Var1", "Var2"))
merged_hlathena_9mer <- merged_hlathena_9mer[!merged_hlathena_9mer$Diff==0,]
merged_hlathena_9mer$length <- "9mer"
rm(merged_peps_umap9)
rm(merged_peps_umap9_filt)

merged_peps_umap10 <- read.csv("merged_peps_umap10mer.csv")
###Binned vs Cor graphs
##Get 10mer cors
merged_peps_umap10_filt <- merged_peps_umap10[merged_peps_umap10$Score >= 0.5,]
hlathena_10mers_counts <- as.data.frame(table(merged_peps_umap10_filt$Allele, merged_peps_umap10_filt$Source))
hlathena_10mers_counts <- hlathena_10mers_counts %>% spread(Var2, Freq)
head(hlathena_10mers_counts)
hlathena_10mer_cor <- melt(cor(hlathena_10mers_counts[,2:10]))

##10mer diffs
merged_peps_umap10$Source <- mapvalues(merged_peps_umap10$Source, matrix_names[,1], matrix_names[,2])
names <- sort(unique(merged_peps_umap9$Source))

ab <- matrix( c(min(merged_peps_umap10$X),min(merged_peps_umap10$Y),max(merged_peps_umap10$X),max(merged_peps_umap10$Y)), 2, 2)      
#merged_peps_umap10 <- merged_peps_umap10[,-c(4,5)]
#merged_peps_umap10 <- distinct(merged_peps_umap10)
matrix_bins_10mers <- tibble(x=character(), y=character(), diff=numeric())
library(ash)
for (i in 1:9)
{
  for (j in i:9)
  {
    x <- as.matrix(merged_peps_umap10[merged_peps_umap10$Source == names[i],2:3])  
    nbin <- c( 40, 40)          
    bins <- bin2(x, ab, nbin)
    bin_sol <- bins$nc
    bin_sol1 <- bin_sol/(dim(x)[1])
    
    y <- as.matrix(merged_peps_umap10[merged_peps_umap10$Source == names[j],2:3])  
    #x <- matrix( rnorm(200), 100 , 2)
    bins <- bin2(y, ab, nbin)
    bin_sol <- bins$nc
    bin_sol2 <- bin_sol/(dim(y)[1])
    diff_bins <- bin_sol1 - bin_sol2
    sum_bins <- sum(abs(diff_bins))
    matrix_bins_10mers <- matrix_bins_10mers %>% add_row(x=names[i], y=names[j], diff=sum_bins)
  }
}


matrix_bins_10mers$diff <- matrix_bins_10mers$diff * 100/2
write.csv(matrix_bins_10mers, "matrix_bins_10mers.csv", row.names=FALSE, quote=FALSE)
library(ggplot2)
matrix_bins_10mers <- as.data.frame(matrix_bins_10mers)
ggplot(data = matrix_bins_10mers, aes(y, x, fill = diff))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(0,100), 
                       name="Percent Peptides Difference") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ggtitle("10mer Differential Peptides")


head(hlathena_10mer_cor)
matrix_bins_10mers <- as.data.frame(matrix_bins_10mers)
colnames(matrix_bins_10mers) <- c("Var1", "Var2", "Diff")
head(matrix_bins_10mers)

merged_hlathena_10mer <- merge(hlathena_10mer_cor, matrix_bins_10mers, by=c("Var1", "Var2"))
merged_hlathena_10mer <- merged_hlathena_10mer[!merged_hlathena_10mer$Diff==0,]
merged_hlathena_10mer$length <- "10mer"
rm(merged_peps_umap10)
rm(merged_peps_umap10_filt)


merged_peps_umap11 <- read.csv("merged_peps_umap11mer.csv")
merged_peps_umap11$Source <- mapvalues(merged_peps_umap11$Source, matrix_names[,1], matrix_names[,2])

###Binned vs Cor graphs
##Get 11mer cors
merged_peps_umap11_filt <- merged_peps_umap11[merged_peps_umap11$Score >= 0.5,]
hlathena_11mers_counts <- as.data.frame(table(merged_peps_umap11_filt$Allele, merged_peps_umap11_filt$Source))
hlathena_11mers_counts <- hlathena_11mers_counts %>% spread(Var2, Freq)
head(hlathena_11mers_counts)
hlathena_11mer_cor <- melt(cor(hlathena_11mers_counts[,2:10]))

##11mer diffs
ab <- matrix( c(min(merged_peps_umap11$X),min(merged_peps_umap11$Y),max(merged_peps_umap11$X),max(merged_peps_umap11$Y)), 2, 2)      
#names <- unique(merged_peps_umap11$Source)
#merged_peps_umap11 <- merged_peps_umap11[,-c(4,5)]
#merged_peps_umap11 <- distinct(merged_peps_umap11)
matrix_bins_11mers <- tibble(x=character(), y=character(), diff=numeric())
library(ash)
for (i in 1:9)
{
  for (j in i:9)
  {
    x <- as.matrix(merged_peps_umap11[merged_peps_umap11$Source == names[i],2:3])  
    nbin <- c( 40, 40)          
    bins <- bin2(x, ab, nbin)
    bin_sol <- bins$nc
    bin_sol1 <- bin_sol/(dim(x)[1])
    
    y <- as.matrix(merged_peps_umap11[merged_peps_umap11$Source == names[j],2:3])  
    #x <- matrix( rnorm(200), 100 , 2)
    bins <- bin2(y, ab, nbin)
    bin_sol <- bins$nc
    bin_sol2 <- bin_sol/(dim(y)[1])
    diff_bins <- bin_sol1 - bin_sol2
    sum_bins <- sum(abs(diff_bins))
    matrix_bins_11mers <- matrix_bins_11mers %>% add_row(x=names[i], y=names[j], diff=sum_bins)
  }
}


matrix_bins_11mers$diff <- matrix_bins_11mers$diff * 100/2
write.csv(matrix_bins_11mers, "matrix_bins_11mers.csv", row.names=FALSE, quote=FALSE)
library(ggplot2)
matrix_bins_11mers <- as.data.frame(matrix_bins_11mers)
ggplot(data = matrix_bins_11mers, aes(y, x, fill = diff))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(0,100), 
                       name="Percent Peptides Difference") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ggtitle("11mer Differential Peptides")


head(hlathena_11mer_cor)
matrix_bins_11mers <- as.data.frame(matrix_bins_11mers)
colnames(matrix_bins_11mers) <- c("Var1", "Var2", "Diff")
head(matrix_bins_11mers)

merged_hlathena_11mer <- merge(hlathena_11mer_cor, matrix_bins_11mers, by=c("Var1", "Var2"))
merged_hlathena_11mer <- merged_hlathena_11mer[!merged_hlathena_11mer$Diff==0,]
merged_hlathena_11mer$length <- "11mer"
rm(merged_peps_umap11)
rm(merged_peps_umap11_filt)

merged_hlathena_all <- rbind(merged_hlathena_8mer,merged_hlathena_9mer,merged_hlathena_10mer,
                             merged_hlathena_11mer)

write.csv(merged_hlathena_all, "merged_hlathena_all_cor_vs_diff.csv")
merged_hlathena_all$length <- factor(merged_hlathena_all$length, levels=c("8mer", "9mer", "10mer", "11mer"))
ggplot(merged_hlathena_all, aes(x=Diff, y=value, color=length)) + 
  geom_point() + geom_smooth(method="lm", formula=y~x, aes(color=length), se=F) + ylim(0,1) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Matrix Binned Difference") + ylab("Correlation Coefficient") + ggtitle("HLAthena Correlation vs Peptide Difference")

write.csv(matrix_bins_8mers, "matrix_bins_8mers.csv", row.names = FALSE)
write.csv(matrix_bins_9mers, "matrix_bins_9mers.csv", row.names = FALSE)
write.csv(matrix_bins_10mers, "matrix_bins_10mers.csv", row.names = FALSE)
write.csv(matrix_bins_11mers, "matrix_bins_11mers.csv", row.names = FALSE)

#Figure 6
matrix_bins_8mers <- read.csv("matrix_bins_8mers.csv")
matrix_bins_9mers <- read.csv("matrix_bins_9mers.csv")
matrix_bins_10mers <- read.csv("matrix_bins_10mers.csv")
matrix_bins_11mers <- read.csv("matrix_bins_11mers.csv")


###See master_merge_code_below
netmhcpan_merged_8mer_counts <- as.data.frame(table(netmhcpan_merged_8mer$Allele, netmhcpan_merged_8mer$Source))
netmhcpan_merged_8mer_counts <- netmhcpan_merged_8mer_counts %>% spread(Var2, Freq)
netmhcpan_merged_8mer_counts <- netmhcpan_merged_8mer_counts[netmhcpan_merged_8mer_counts$Var1 %in% big_allele_list,]
netmhcpan_8mer_cor <- melt(cor(netmhcpan_merged_8mer_counts[,2:10]))
merged_netmhcpan_8mer <- merge(netmhcpan_8mer_cor, matrix_bins_8mers, by=c("Var1", "Var2"))
merged_netmhcpan_8mer <- merged_netmhcpan_8mer[!merged_netmhcpan_8mer$Diff==0,]
merged_netmhcpan_8mer$length <- "8mer"


netmhcpan_merged_9mer_counts <- as.data.frame(table(netmhcpan_merged_9mer$Allele, netmhcpan_merged_9mer$Source))
netmhcpan_merged_9mer_counts <- netmhcpan_merged_9mer_counts %>% spread(Var2, Freq)
netmhcpan_merged_9mer_counts <- netmhcpan_merged_9mer_counts[netmhcpan_merged_9mer_counts$Var1 %in% big_allele_list,]
netmhcpan_9mer_cor <- melt(cor(netmhcpan_merged_9mer_counts[,2:10]))
merged_netmhcpan_9mer <- merge(netmhcpan_9mer_cor, matrix_bins_9mers, by=c("Var1", "Var2"))
merged_netmhcpan_9mer <- merged_netmhcpan_9mer[!merged_netmhcpan_9mer$Diff==0,]
merged_netmhcpan_9mer$length <- "9mer"


netmhcpan_merged_10mer_counts <- as.data.frame(table(netmhcpan_merged_10mer$Allele, netmhcpan_merged_10mer$Source))
netmhcpan_merged_10mer_counts <- netmhcpan_merged_10mer_counts %>% spread(Var2, Freq)
netmhcpan_merged_10mer_counts <- netmhcpan_merged_10mer_counts[netmhcpan_merged_10mer_counts$Var1 %in% big_allele_list,]
netmhcpan_10mer_cor <- melt(cor(netmhcpan_merged_10mer_counts[,2:10]))
merged_netmhcpan_10mer <- merge(netmhcpan_10mer_cor, matrix_bins_10mers, by=c("Var1", "Var2"))
merged_netmhcpan_10mer <- merged_netmhcpan_10mer[!merged_netmhcpan_10mer$Diff==0,]
merged_netmhcpan_10mer$length <- "10mer"


netmhcpan_merged_11mer_counts <- as.data.frame(table(netmhcpan_merged_11mer$Allele, netmhcpan_merged_11mer$Source))
netmhcpan_merged_11mer_counts <- netmhcpan_merged_11mer_counts %>% spread(Var2, Freq)
netmhcpan_merged_11mer_counts <- netmhcpan_merged_11mer_counts[netmhcpan_merged_11mer_counts$Var1 %in% big_allele_list,]
netmhcpan_11mer_cor <- melt(cor(netmhcpan_merged_11mer_counts[,2:10]))
merged_netmhcpan_11mer <- merge(netmhcpan_11mer_cor, matrix_bins_11mers, by=c("Var1", "Var2"))
merged_netmhcpan_11mer <- merged_netmhcpan_11mer[!merged_netmhcpan_11mer$Diff==0,]
merged_netmhcpan_11mer$length <- "11mer"

merged_netmhcpan_all <- rbind(merged_netmhcpan_8mer,merged_netmhcpan_9mer,merged_netmhcpan_10mer,
                              merged_netmhcpan_11mer)
write.csv(merged_netmhcpan_all, "merged_netmhcpan_all.csv")
merged_netmhcpan_all$length <- factor(merged_netmhcpan_all$length, levels=c("8mer", "9mer", "10mer", "11mer"))
library(ggplot2)
library(tidyverse)
ggplot(merged_netmhcpan_all, aes(x=Diff, y=value, color=length)) + 
  geom_point() + geom_smooth(method="lm", formula=y~x, aes(color=length), se=F) + ylim(0,1) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Matrix Binned Difference") + ylab("Correlation Coefficient") + ggtitle("netMHCpan Correlation vs Peptide Difference")


####MHCFLURRY Binned diff vs COR
library(dplyr)
mhcflurry_betaherpes_6a$Source <- "betaherpes_6a"
mhcflurry_bk$Source <- "bk"
mhcflurry_cmv$Source <- "cmv"
mhcflurry_covid$Source <- "covid"
mhcflurry_ebv$Source <- "ebv"
mhcflurry_herpes_1_17$Source <- "herpes_1_17"
mhcflurry_herpes_2_hg52$Source <- "herpes_2_hg52"
mhcflurry_human$Source <- "human"
mhcflurry_random$Source <- "random"

mhcflurry_betaherpes_6a <- mhcflurry_betaherpes_6a[mhcflurry_betaherpes_6a$allele %in% list_alleles[,1] & mhcflurry_betaherpes_6a$mhcflurry_affinity >= 0.5,]
mhcflurry_bk <- mhcflurry_bk[mhcflurry_bk$allele %in% list_alleles[,1] & mhcflurry_bk$mhcflurry_affinity >= 0.5,]
mhcflurry_cmv <- mhcflurry_cmv[mhcflurry_cmv$allele %in% list_alleles[,1] & mhcflurry_cmv$mhcflurry_affinity >= 0.5,]
mhcflurry_covid <- mhcflurry_covid[mhcflurry_covid$allele %in% list_alleles[,1] & mhcflurry_covid$mhcflurry_affinity >= 0.5,]
mhcflurry_ebv <- mhcflurry_ebv[mhcflurry_ebv$allele %in% list_alleles[,1] & mhcflurry_ebv$mhcflurry_affinity >= 0.5,]
mhcflurry_herpes_1_17 <- mhcflurry_herpes_1_17[mhcflurry_herpes_1_17$allele %in% list_alleles[,1] & mhcflurry_herpes_1_17$mhcflurry_affinity >= 0.5,]
mhcflurry_herpes_2_hg52 <- mhcflurry_herpes_2_hg52[mhcflurry_herpes_2_hg52$allele %in% list_alleles[,1] & mhcflurry_herpes_2_hg52$mhcflurry_affinity >= 0.5,]
mhcflurry_human <- mhcflurry_human[mhcflurry_human$allele %in% list_alleles[,1] & mhcflurry_human$mhcflurry_affinity >= 0.5,]
mhcflurry_random <- mhcflurry_random[mhcflurry_random$allele %in% list_alleles[,1] & mhcflurry_random$mhcflurry_affinity >= 0.5,]

merge_mhcflurry_peps <- rbind(mhcflurry_betaherpes_6a,mhcflurry_bk )
merge_mhcflurry_peps <- rbind(merge_mhcflurry_peps,mhcflurry_cmv )
merge_mhcflurry_peps <- rbind(merge_mhcflurry_peps,mhcflurry_covid )
merge_mhcflurry_peps <- rbind(merge_mhcflurry_peps,mhcflurry_ebv)
merge_mhcflurry_peps <- rbind(merge_mhcflurry_peps,mhcflurry_herpes_1_17)
merge_mhcflurry_peps <- rbind(merge_mhcflurry_peps,mhcflurry_herpes_2_hg52 )
merge_mhcflurry_peps <- rbind(merge_mhcflurry_peps,mhcflurry_human )
merge_mhcflurry_peps <- rbind(merge_mhcflurry_peps,mhcflurry_random )

rm(mhcflurry_betaherpes_6a)
rm(mhcflurry_bk)
rm(mhcflurry_cmv)
rm(mhcflurry_covid)
rm(mhcflurry_ebv)
rm(mhcflurry_herpes_1_17)
rm(mhcflurry_herpes_2_hg52)
rm(mhcflurry_human)
rm(mhcflurry_random)
gc()

colnames(merge_mhcflurry_peps)[2] <- "Allele"
colnames(merge_mhcflurry_peps)[1] <- "Peptide"
merge_mhcflurry_peps <- merge_mhcflurry_peps[merge_mhcflurry_peps$mhcflurry_affinity >= 0.5,]
mhcflurry_merged_8mer <- merge_mhcflurry_peps[nchar(merge_mhcflurry_peps$Peptide) == 8,]
mhcflurry_merged_9mer <- merge_mhcflurry_peps[nchar(merge_mhcflurry_peps$Peptide) == 9,]
mhcflurry_merged_10mer <- merge_mhcflurry_peps[nchar(merge_mhcflurry_peps$Peptide) == 10,]
mhcflurry_merged_11mer <- merge_mhcflurry_peps[nchar(merge_mhcflurry_peps$Peptide) == 11,]
#rm(merge_mhcflurry_peps)

library(tidyr)
library(reshape2)
mhcflurry_merged_8mer_counts <- as.data.frame(table(mhcflurry_merged_8mer$Allele, mhcflurry_merged_8mer$Source))
mhcflurry_merged_8mer_counts <- mhcflurry_merged_8mer_counts %>% spread(Var2, Freq)
mhcflurry_merged_8mer_counts <- mhcflurry_merged_8mer_counts[mhcflurry_merged_8mer_counts$Var1 %in% bighlathena_list,]
mhcflurry_8mer_cor <- melt(cor(mhcflurry_merged_8mer_counts[,2:10]))
merged_mhcflurry_8mer <- merge(mhcflurry_8mer_cor, matrix_bins_8mers, by=c("Var1", "Var2"))
merged_mhcflurry_8mer <- merged_mhcflurry_8mer[!merged_mhcflurry_8mer$Diff==0,]
merged_mhcflurry_8mer$length <- "8mer"


mhcflurry_merged_9mer_counts <- as.data.frame(table(mhcflurry_merged_9mer$Allele, mhcflurry_merged_9mer$Source))
mhcflurry_merged_9mer_counts <- mhcflurry_merged_9mer_counts %>% spread(Var2, Freq)
mhcflurry_merged_9mer_counts <- mhcflurry_merged_9mer_counts[mhcflurry_merged_9mer_counts$Var1 %in% bighlathena_list,]
mhcflurry_9mer_cor <- melt(cor(mhcflurry_merged_9mer_counts[,2:10]))
merged_mhcflurry_9mer <- merge(mhcflurry_9mer_cor, matrix_bins_9mers, by=c("Var1", "Var2"))
merged_mhcflurry_9mer <- merged_mhcflurry_9mer[!merged_mhcflurry_9mer$Diff==0,]
merged_mhcflurry_9mer$length <- "9mer"


mhcflurry_merged_10mer_counts <- as.data.frame(table(mhcflurry_merged_10mer$Allele, mhcflurry_merged_10mer$Source))
mhcflurry_merged_10mer_counts <- mhcflurry_merged_10mer_counts %>% spread(Var2, Freq)
mhcflurry_merged_10mer_counts <- mhcflurry_merged_10mer_counts[mhcflurry_merged_10mer_counts$Var1 %in% bighlathena_list,]
mhcflurry_10mer_cor <- melt(cor(mhcflurry_merged_10mer_counts[,2:10]))
merged_mhcflurry_10mer <- merge(mhcflurry_10mer_cor, matrix_bins_10mers, by=c("Var1", "Var2"))
merged_mhcflurry_10mer <- merged_mhcflurry_10mer[!merged_mhcflurry_10mer$Diff==0,]
merged_mhcflurry_10mer$length <- "10mer"


mhcflurry_merged_11mer_counts <- as.data.frame(table(mhcflurry_merged_11mer$Allele, mhcflurry_merged_11mer$Source))
mhcflurry_merged_11mer_counts <- mhcflurry_merged_11mer_counts %>% spread(Var2, Freq)
mhcflurry_merged_11mer_counts <- mhcflurry_merged_11mer_counts[mhcflurry_merged_11mer_counts$Var1 %in% bighlathena_list,]
mhcflurry_11mer_cor <- melt(cor(mhcflurry_merged_11mer_counts[,2:10]))
merged_mhcflurry_11mer <- merge(mhcflurry_11mer_cor, matrix_bins_11mers, by=c("Var1", "Var2"))
merged_mhcflurry_11mer <- merged_mhcflurry_11mer[!merged_mhcflurry_11mer$Diff==0,]
merged_mhcflurry_11mer$length <- "11mer"

merged_mhcflurry_all <- rbind(merged_mhcflurry_8mer,merged_mhcflurry_9mer,merged_mhcflurry_10mer,
                              merged_mhcflurry_11mer, all=TRUE)

merged_mhcflurry_all[is.na(merged_mhcflurry_all)] <- 0
library(ggplot2)
library(tidyverse)
merged_mhcflurry_all$length <- factor(merged_mhcflurry_all$length, levels=c("8mer", "9mer", "10mer", "11mer"))
merged_mhcflurry_all <- head(merged_mhcflurry_all,-1)
write.csv(merged_mhcflurry_all, "merged_mhcflurry_all.csv")
ggplot(merged_mhcflurry_all, aes(x=Diff, y=value, color=length)) + 
  geom_point() + geom_smooth(method="lm", formula=y~x, aes(color=length), se=F) + ylim(0,1) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Matrix Binned Difference") + ylab("Correlation Coefficient") + ggtitle("MHCflurry Correlation vs Peptide Difference")





####MHCNUGGETS BINNED DIFF VS COR
library(dplyr)
mhcnuggets_betaherpes_6a$Source <- "betaherpes_6a"
mhcnuggets_bk$Source <- "bk"
mhcnuggets_cmv$Source <- "cmv"
mhcnuggets_covid$Source <- "covid"
mhcnuggets_ebv$Source <- "ebv"
mhcnuggets_herpes_1_17$Source <- "herpes_1_17"
mhcnuggets_herpes_2_hg52$Source <- "herpes_2_hg52"
mhcnuggets_human$Source <- "human"
mhcnuggets_random$Source <- "random"

list_alleles$Allele <- gsub('^(.{3})(.*)$', '\\1:\\2', list_alleles$Allele)
mhcnuggets_betaherpes_6a <- mhcnuggets_betaherpes_6a[mhcnuggets_betaherpes_6a$allele %in% list_alleles[,1] & mhcnuggets_betaherpes_6a$mhcnuggets_affinity >= 0.5,]
mhcnuggets_bk <- mhcnuggets_bk[mhcnuggets_bk$allele %in% list_alleles[,1] & mhcnuggets_bk$mhcnuggets_affinity >= 0.5,]
mhcnuggets_cmv <- mhcnuggets_cmv[mhcnuggets_cmv$allele %in% list_alleles[,1] & mhcnuggets_cmv$mhcnuggets_affinity >= 0.5,]
mhcnuggets_covid <- mhcnuggets_covid[mhcnuggets_covid$allele %in% list_alleles[,1] & mhcnuggets_covid$mhcnuggets_affinity >= 0.5,]
mhcnuggets_ebv <- mhcnuggets_ebv[mhcnuggets_ebv$allele %in% list_alleles[,1] & mhcnuggets_ebv$mhcnuggets_affinity >= 0.5,]
mhcnuggets_herpes_1_17 <- mhcnuggets_herpes_1_17[mhcnuggets_herpes_1_17$allele %in% list_alleles[,1] & mhcnuggets_herpes_1_17$mhcnuggets_affinity >= 0.5,]
mhcnuggets_herpes_2_hg52 <- mhcnuggets_herpes_2_hg52[mhcnuggets_herpes_2_hg52$allele %in% list_alleles[,1] & mhcnuggets_herpes_2_hg52$mhcnuggets_affinity >= 0.5,]
mhcnuggets_human <- mhcnuggets_human[mhcnuggets_human$allele %in% list_alleles[,1] & mhcnuggets_human$mhcnuggets_affinity >= 0.5,]
mhcnuggets_random <- mhcnuggets_random[mhcnuggets_random$allele %in% list_alleles[,1] & mhcnuggets_random$mhcnuggets_affinity >= 0.5,]

merge_mhcnuggets_peps <- rbind(mhcnuggets_betaherpes_6a,mhcnuggets_bk )
merge_mhcnuggets_peps <- rbind(merge_mhcnuggets_peps,mhcnuggets_cmv )
merge_mhcnuggets_peps <- rbind(merge_mhcnuggets_peps,mhcnuggets_covid )
merge_mhcnuggets_peps <- rbind(merge_mhcnuggets_peps,mhcnuggets_ebv)
merge_mhcnuggets_peps <- rbind(merge_mhcnuggets_peps,mhcnuggets_herpes_1_17)
merge_mhcnuggets_peps <- rbind(merge_mhcnuggets_peps,mhcnuggets_herpes_2_hg52 )
merge_mhcnuggets_peps <- rbind(merge_mhcnuggets_peps,mhcnuggets_human )
merge_mhcnuggets_peps <- rbind(merge_mhcnuggets_peps,mhcnuggets_random )

rm(mhcnuggets_betaherpes_6a)
rm(mhcnuggets_bk)
rm(mhcnuggets_cmv)
rm(mhcnuggets_covid)
rm(mhcnuggets_ebv)
rm(mhcnuggets_herpes_1_17)
rm(mhcnuggets_herpes_2_hg52)
rm(mhcnuggets_human)
rm(mhcnuggets_random)
gc()

merge_mhcnuggets_peps$allele <- gsub(":", "", merge_mhcnuggets_peps$allele)
colnames(merge_mhcnuggets_peps)[3] <- "Allele"
mhcnuggets_merged_8mer <- merge_mhcnuggets_peps[nchar(merge_mhcnuggets_peps$Peptide) == 8,]
mhcnuggets_merged_9mer <- merge_mhcnuggets_peps[nchar(merge_mhcnuggets_peps$Peptide) == 9,]
mhcnuggets_merged_10mer <- merge_mhcnuggets_peps[nchar(merge_mhcnuggets_peps$Peptide) == 10,]
mhcnuggets_merged_11mer <- merge_mhcnuggets_peps[nchar(merge_mhcnuggets_peps$Peptide) == 11,]
#rm(merge_mhcnuggets_peps)


library(tidyr)
library(reshape2)
mhcnuggets_merged_8mer_counts <- as.data.frame(table(mhcnuggets_merged_8mer$Allele, mhcnuggets_merged_8mer$Source))
mhcnuggets_merged_8mer_counts <- mhcnuggets_merged_8mer_counts %>% spread(Var2, Freq)
mhcnuggets_8mer_cor <- melt(cor(mhcnuggets_merged_8mer_counts[,2:10]))
merged_mhcnuggets_8mer <- merge(mhcnuggets_8mer_cor, matrix_bins_8mers, by=c("Var1", "Var2"), all=TRUE)
#merged_mhcnuggets_8mer <- merged_mhcnuggets_8mer[!merged_mhcnuggets_8mer$Diff==0,]
merged_mhcnuggets_8mer$length <- "8mer"


mhcnuggets_merged_9mer_counts <- as.data.frame(table(mhcnuggets_merged_9mer$Allele, mhcnuggets_merged_9mer$Source))
mhcnuggets_merged_9mer_counts <- mhcnuggets_merged_9mer_counts %>% spread(Var2, Freq)
mhcnuggets_9mer_cor <- melt(cor(mhcnuggets_merged_9mer_counts[,2:10]))
merged_mhcnuggets_9mer <- merge(mhcnuggets_9mer_cor, matrix_bins_9mers, by=c("Var1", "Var2"), all=TRUE)
#merged_mhcnuggets_9mer <- merged_mhcnuggets_9mer[!merged_mhcnuggets_9mer$Diff==0,]
merged_mhcnuggets_9mer$length <- "9mer"


mhcnuggets_merged_10mer_counts <- as.data.frame(table(mhcnuggets_merged_10mer$Allele, mhcnuggets_merged_10mer$Source))
mhcnuggets_merged_10mer_counts <- mhcnuggets_merged_10mer_counts %>% spread(Var2, Freq)
mhcnuggets_10mer_cor <- melt(cor(mhcnuggets_merged_10mer_counts[,2:10]))
merged_mhcnuggets_10mer <- merge(mhcnuggets_10mer_cor, matrix_bins_10mers, by=c("Var1", "Var2"), all=TRUE)
#merged_mhcnuggets_10mer <- merged_mhcnuggets_10mer[!merged_mhcnuggets_10mer$Diff==0,]
merged_mhcnuggets_10mer$length <- "10mer"


mhcnuggets_merged_11mer_counts <- as.data.frame(table(mhcnuggets_merged_11mer$Allele, mhcnuggets_merged_11mer$Source))
mhcnuggets_merged_11mer_counts <- mhcnuggets_merged_11mer_counts %>% spread(Var2, Freq)
mhcnuggets_11mer_cor <- melt(cor(mhcnuggets_merged_11mer_counts[,2:10]))
merged_mhcnuggets_11mer <- merge(mhcnuggets_11mer_cor, matrix_bins_11mers, by=c("Var1", "Var2"), all=TRUE)
#merged_mhcnuggets_11mer <- merged_mhcnuggets_11mer[!merged_mhcnuggets_11mer$Diff==0,]
merged_mhcnuggets_11mer$length <- "11mer"

merged_mhcnuggets_all <- rbind(merged_mhcnuggets_8mer,merged_mhcnuggets_9mer,merged_mhcnuggets_10mer,
                               merged_mhcnuggets_11mer)
head(merged_mhcnuggets_all)
merged_mhcnuggets_all[is.na(merged_mhcnuggets_all)] <- 0
merged_mhcnuggets_all <- merged_mhcnuggets_all[!merged_mhcnuggets_all$Diff==0,]

library(ggplot2)
library(tidyverse)
merged_mhcnuggets_all$length <- factor(merged_mhcnuggets_all$length, levels=c("8mer", "9mer", "10mer", "11mer"))

ggplot(merged_mhcnuggets_all, aes(x=Diff, y=value, color=length)) + 
  geom_point() + geom_smooth(method="lm", formula=y~x, aes(color=length), se=F) + ylim(0,1) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Matrix Binned Difference") + ylab("Correlation Coefficient") + ggtitle("MHCnuggets Correlation vs Peptide Difference")





#Figure 7
merged_peps_random_all <- read.csv("merged_peps_random_all.csv")
#8mers
merged_peps_random_8mers <- merge(merged_peps_random_all, merged_peps_umap8, by="Peptide")
rm(merged_peps_random_all)
melt_merged_peps_random_8mers <- melt(merged_peps_random_8mers, id=c("Peptide", "Allele", "X", "Y"))
rm(merged_peps_random_8mers)
gc()

i <- 1
allele_diff_bins_master_8mers <- tibble(X=character(), Y=character(), Diff=numeric(), Length=character(), Allele=character())
colnames(melt_merged_peps_random_8mers)[5:6] <- c("Tool", "Score")

#fix kmer length as var
ab <- matrix( c(min(melt_merged_peps_random_8mers$X),min(melt_merged_peps_random_8mers$Y),max(melt_merged_peps_random_8mers$X),max(melt_merged_peps_random_8mers$Y)), 2, 2)      
library(ash)
#library(ggplot2)
for (i in 1:length(list_alleles[,1]))
{
  
  x <- as.matrix(unique(melt_merged_peps_random_8mers[melt_merged_peps_random_8mers$Allele == list_alleles[i,1] & melt_merged_peps_random_8mers$Score >= 0.5,3:4]))
  nbin <- c( 40, 40)                      # 1600 bins #compare binning to umap -> more bins(?)
  bins <- bin2(x, ab, nbin)
  bin_sol <- bins$nc
  sum(bin_sol)
  dim(x)[1]
  bin_sol1 <- bin_sol/sum(bin_sol)
  
  y <- as.matrix(unique(melt_merged_peps_random_8mers[melt_merged_peps_random_8mers$Allele == list_alleles[i,1] & melt_merged_peps_random_8mers$Score < 0.5,3:4]))
  bins <- bin2(y, ab, nbin)
  bin_sol <- bins$nc
  sum(bin_sol)
  dim(y)[1]
  bin_sol2 <- bin_sol/sum(bin_sol)
  diff_bins <- bin_sol1 - bin_sol2
  colnames(diff_bins) <- 1:40
  rownames(diff_bins) <- 1:40
  
  melted_bins <- melt(diff_bins, na.rm = TRUE)
  melted_bins$value <- melted_bins$value*100
  melted_bins <- tibble(melted_bins)
  melted_bins$Length <- "8mer"
  melted_bins$Allele <- list_alleles[i,1]
  colnames(melted_bins)[1:3] <- c("X", "Y", "Diff")
  allele_diff_bins_master_8mers <- rbind(allele_diff_bins_master_8mers, melted_bins)
}
gc()
allele_diff_bins_master_8mers_filt <- allele_diff_bins_master_8mers
rm(allele_diff_bins_master_8mers)
allele_diff_bins_master_8mers_filt[allele_diff_bins_master_8mers_filt == 0] <- NA
allele_diff_bins_master_8mers_filt<-allele_diff_bins_master_8mers_filt[complete.cases(allele_diff_bins_master_8mers_filt),]

allele_diff_bins_master_8mers_filt <- merge(allele_diff_bins_master_8mers_filt, list_alleles, by="Allele")
temp <- allele_diff_bins_master_8mers_filt
temp$Label <- paste(temp$X, temp$Y, sep="-")
library(tidyr)
temp <- temp %>%
  group_by(Label) %>%
  slice(which.max(Diff))

library(ggplot2)
temp <- dplyr::filter(temp, Diff>0.2)
temp <- temp[order(-temp$Diff),]
plot1 <- ggplot(temp, aes(x=Y, y=X, alpha=Diff, fill=Supertype)) + 
  geom_tile()  +theme_minimal() +guides(alpha="none")+expand_limits(x = 0, y = 0)+
  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ggtitle("Binders Enrichment")

allele_diff_bins_master_8mers_filt[which.max(allele_diff_bins_master_8mers_filt$Diff),]
##Example B0801 highest 8mer

x <- as.matrix(unique(melt_merged_peps_random_8mers[melt_merged_peps_random_8mers$Allele == "B0801" & melt_merged_peps_random_8mers$Score >= 0.5,3:4]))
ab <- matrix( c(min(melt_merged_peps_random_8mers$X),min(melt_merged_peps_random_8mers$Y),max(melt_merged_peps_random_8mers$X),max(melt_merged_peps_random_8mers$Y)), 2, 2)      
nbin <- c( 40, 40)                      # 1600 bins #compare binning to umap -> more bins(?)
bins <- bin2(x, ab, nbin)
bin_sol <- bins$nc
sum(bin_sol)
dim(x)[1]
bin_sol1 <- bin_sol/sum(bin_sol)

y <- as.matrix(unique(melt_merged_peps_random_8mers[melt_merged_peps_random_8mers$Allele == "B0801" & melt_merged_peps_random_8mers$Score < 0.5,3:4]))
bins <- bin2(y, ab, nbin)
bin_sol <- bins$nc
sum(bin_sol)
dim(y)[1]
bin_sol2 <- bin_sol/sum(bin_sol)
diff_bins <- bin_sol1 - bin_sol2
colnames(diff_bins) <- 1:40
rownames(diff_bins) <- 1:40

melted_bins <- melt(diff_bins, na.rm = TRUE)
melted_bins$value <- melted_bins$value*100
title <- "B0801 Binder vs Nonbinder"

ggplot(data = melted_bins, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Enrichment") +
  theme_minimal()+ 
  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ggtitle(title)

###

#9mers
list_alleles <- read.csv("list_alleles_merged.csv", header=FALSE)
colnames(list_alleles) <- c("Allele", "Supertype")

library(reshape2)
library(ggplot2)

merged_peps_random_all <- read.csv("merged_peps_random_all.csv")
merged_peps_umap9 <- read.csv("merged_peps_umap9mer.csv")
library(dplyr)
merged_peps_umap9 <- distinct(merged_peps_umap9[,1:3])

merged_peps_random_9mers <- merge(merged_peps_random_all, merged_peps_umap9, by="Peptide")
rm(merged_peps_umap9)
rm(merged_peps_random_all)
gc()
melt_merged_peps_random_9mers <- melt(merged_peps_random_9mers, id=c("Peptide", "Allele", "X", "Y"))
rm(merged_peps_random_9mers)
gc()
i <- 1
colnames(melt_merged_peps_random_9mers)[5:6] <- c("Tool", "Score")

allele_diff_bins_master_9mers <- tibble(X=character(), Y=character(), Diff=numeric(), Length=character(), Allele=character())

#fix kmer length as var
ab <- matrix( c(min(melt_merged_peps_random_9mers$X),min(melt_merged_peps_random_9mers$Y),max(melt_merged_peps_random_9mers$X),max(melt_merged_peps_random_9mers$Y)), 2, 2)      
library(ash)
#library(ggplot2)
for (i in 1:dim(list_alleles)[1])
{
  
  x <- as.matrix(unique(melt_merged_peps_random_9mers[melt_merged_peps_random_9mers$Allele == list_alleles[i,1] & melt_merged_peps_random_9mers$Score >= 0.5,3:4]))
  nbin <- c( 40, 40)                      # 1600 bins #compare binning to umap -> more bins(?)
  bins <- bin2(x, ab, nbin)
  bin_sol <- bins$nc
  sum(bin_sol)
  dim(x)[1]
  bin_sol1 <- bin_sol/sum(bin_sol)
  
  y <- as.matrix(unique(melt_merged_peps_random_9mers[melt_merged_peps_random_9mers$Allele == list_alleles[i,1] & melt_merged_peps_random_9mers$Score < 0.5,3:4]))
  bins <- bin2(y, ab, nbin)
  bin_sol <- bins$nc
  sum(bin_sol)
  dim(y)[1]
  bin_sol2 <- bin_sol/sum(bin_sol)
  diff_bins <- bin_sol1 - bin_sol2
  colnames(diff_bins) <- 1:40
  rownames(diff_bins) <- 1:40
  #master_diff_bins <- master_diff_bins+diff_bins
  #sum_bins <- sum(abs(diff_bins))
  #allele_bind_vs_nonbind <- allele_bind_vs_nonbind %>% add_row(Source="9mer", Allele=list_alleles[i], diff=sum_bins)
  
  melted_bins <- melt(diff_bins, na.rm = TRUE)
  melted_bins$value <- melted_bins$value*100
  melted_bins <- tibble(melted_bins)
  melted_bins$Length <- "9mer"
  melted_bins$Allele <- list_alleles[i,1]
  colnames(melted_bins)[1:3] <- c("X", "Y", "Diff")
  allele_diff_bins_master_9mers <- rbind(allele_diff_bins_master_9mers, melted_bins)
}

gc()
allele_diff_bins_master_9mers_filt <- allele_diff_bins_master_9mers

rm(allele_diff_bins_master_9mers)
allele_diff_bins_master_9mers_filt[allele_diff_bins_master_9mers_filt == 0] <- NA
allele_diff_bins_master_9mers_filt<-allele_diff_bins_master_9mers_filt[complete.cases(allele_diff_bins_master_9mers_filt),]

allele_diff_bins_master_9mers_filt <- merge(allele_diff_bins_master_9mers_filt, list_alleles, by="Allele")
head(allele_diff_bins_master_9mers_filt)
temp <- allele_diff_bins_master_9mers_filt
temp$Label <- paste(temp$X, temp$Y, sep="-")
library(tidyr)
temp <- temp %>%
  group_by(Label) %>%
  slice(which.max(Diff))

library(ggplot2)
#options(repr.plot.width =9, repr.plot.height =9)
#ggplot(dplyr::filter(allele_diff_bins_master_9mers_filt, Diff>0.2), aes(x=X, y=Y, alpha=Diff, color=Supertype)) + 
#  geom_point()  +theme_minimal() +guides(alpha="none")+expand_limits(x = 0, y = 0)+
#  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ggtitle("9mers Binders Enrichment")
temp <- dplyr::filter(temp, Diff>0.2)
temp <- temp[order(-temp$Diff),]
plot1 <- ggplot(temp, aes(x=Y, y=X, alpha=Diff, fill=Supertype)) + 
  geom_tile()  +theme_minimal() +guides(alpha="none")+expand_limits(x = 0, y = 0)+
  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ggtitle("Binders Enrichment")

#+ theme_minimal()+guides(alpha="none")+
#  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ggtitle("9mers Binders Diff>0.2%")

allele_diff_bins_master_9mers_filt[which.max(allele_diff_bins_master_9mers_filt$Diff),]
allele_diff_bins_master_9mers_filt[which.min(allele_diff_bins_master_9mers_filt$Diff),]
colnames(melt_merged_peps_random_9mers)[5:6] <- c("Tool", "Score")

##Example A0101 highest 9mer
#27    30 0.868 9mer   A0101 

x <- as.matrix(unique(melt_merged_peps_random_9mers[melt_merged_peps_random_9mers$Allele == "A0101" & melt_merged_peps_random_9mers$Score >= 0.5,3:4]))
ab <- matrix( c(min(melt_merged_peps_random_9mers$X),min(melt_merged_peps_random_9mers$Y),max(melt_merged_peps_random_9mers$X),max(melt_merged_peps_random_9mers$Y)), 2, 2)      
nbin <- c( 40, 40)                      # 1600 bins #compare binning to umap -> more bins(?)
bins <- bin2(x, ab, nbin)
bin_sol <- bins$nc
sum(bin_sol)
dim(x)[1]
bin_sol1 <- bin_sol/sum(bin_sol)

y <- as.matrix(unique(melt_merged_peps_random_9mers[melt_merged_peps_random_9mers$Allele == "A0101" & melt_merged_peps_random_9mers$Score < 0.5,3:4]))
bins <- bin2(y, ab, nbin)
bin_sol <- bins$nc
sum(bin_sol)
dim(y)[1]
bin_sol2 <- bin_sol/sum(bin_sol)
diff_bins <- bin_sol1 - bin_sol2
colnames(diff_bins) <- 1:40
rownames(diff_bins) <- 1:40

melted_bins <- melt(diff_bins, na.rm = TRUE)
melted_bins$value <- melted_bins$value*100
title <- "A0101 Binder vs Nonbinder"

plot2 <- ggplot(data = melted_bins, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Enrichment") +
  theme_minimal()+ 
  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ggtitle(title)


###MIN


x <- as.matrix(unique(melt_merged_peps_random_9mers[melt_merged_peps_random_9mers$Allele == "B4501" & melt_merged_peps_random_9mers$Score >= 0.5,3:4]))
ab <- matrix( c(min(melt_merged_peps_random_9mers$X),min(melt_merged_peps_random_9mers$Y),max(melt_merged_peps_random_9mers$X),max(melt_merged_peps_random_9mers$Y)), 2, 2)      
nbin <- c( 40, 40)                      # 1600 bins #compare binning to umap -> more bins(?)
bins <- bin2(x, ab, nbin)
bin_sol <- bins$nc
sum(bin_sol)
dim(x)[1]
bin_sol1 <- bin_sol/sum(bin_sol)

y <- as.matrix(unique(melt_merged_peps_random_9mers[melt_merged_peps_random_9mers$Allele == "B4501" & melt_merged_peps_random_9mers$Score < 0.5,3:4]))
bins <- bin2(y, ab, nbin)
bin_sol <- bins$nc
sum(bin_sol)
dim(y)[1]
bin_sol2 <- bin_sol/sum(bin_sol)
diff_bins <- bin_sol1 - bin_sol2
colnames(diff_bins) <- 1:40
rownames(diff_bins) <- 1:40

melted_bins <- melt(diff_bins, na.rm = TRUE)
melted_bins$value <- melted_bins$value*100
title <- "B4501 Binder vs Nonbinder"

melted_bins[melted_bins == 0] <- NA
melted_bins<-melted_bins[complete.cases(melted_bins),]
plot3 <- ggplot(data = melted_bins, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Enrichment") +
  theme_minimal()+ 
  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ggtitle(title)

#rm(melt_merged_peps_random_9mers)
#gc()
###
#B0801

x <- as.matrix(unique(melt_merged_peps_random_9mers[melt_merged_peps_random_9mers$Allele == "B0801" & melt_merged_peps_random_9mers$Score >= 0.5,3:4]))
ab <- matrix( c(min(melt_merged_peps_random_9mers$X),min(melt_merged_peps_random_9mers$Y),max(melt_merged_peps_random_9mers$X),max(melt_merged_peps_random_9mers$Y)), 2, 2)      
nbin <- c( 40, 40)                      # 1600 bins #compare binning to umap -> more bins(?)
bins <- bin2(x, ab, nbin)
bin_sol <- bins$nc
sum(bin_sol)
dim(x)[1]
bin_sol1 <- bin_sol/sum(bin_sol)

y <- as.matrix(unique(melt_merged_peps_random_9mers[melt_merged_peps_random_9mers$Allele == "B0801" & melt_merged_peps_random_9mers$Score < 0.5,3:4]))
bins <- bin2(y, ab, nbin)
bin_sol <- bins$nc
sum(bin_sol)
dim(y)[1]
bin_sol2 <- bin_sol/sum(bin_sol)
diff_bins <- bin_sol1 - bin_sol2
colnames(diff_bins) <- 1:40
rownames(diff_bins) <- 1:40

melted_bins <- melt(diff_bins, na.rm = TRUE)
melted_bins$value <- melted_bins$value*100
title <- "B0801 Binder vs Nonbinder"

melted_bins[melted_bins == 0] <- NA
melted_bins<-melted_bins[complete.cases(melted_bins),]
plot4 <- ggplot(data = melted_bins, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Enrichment") +
  theme_minimal()+ 
  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ggtitle(title)



##

g1 = ggplotGrob(plot1)
g2 = ggplotGrob(plot2)
g4 = ggplotGrob(plot3)
g3 = ggplotGrob(plot4)

g1$widths <- g2$widths
g3$widths <- g2$widths
g4$widths <- g2$widths
library(gridExtra)
grid.arrange(g1,g2,g3,g4, nrow=2)


#10mers
merged_peps_random_10mers <- merge(merged_peps_random_all, merged_peps_umap10, by="Peptide")
melt_merged_peps_random_10mers <- melt(merged_peps_random_10mers, id=c("Peptide", "Allele", "X", "Y"))
rm(merged_peps_random_10mers)
gc()

i <- 1
allele_diff_bins_master_10mers <- tibble(X=character(), Y=character(), Diff=numeric(), Length=character(), Allele=character())
colnames(melt_merged_peps_random_10mers)[5:6] <- c("Tool", "Score")

#fix kmer length as var
ab <- matrix( c(min(melt_merged_peps_random_10mers$X),min(melt_merged_peps_random_10mers$Y),max(melt_merged_peps_random_10mers$X),max(melt_merged_peps_random_10mers$Y)), 2, 2)      
library(ash)
#library(ggplot2)
for (i in 1:length(list_alleles[,1]))
{
  
  x <- as.matrix(unique(melt_merged_peps_random_10mers[melt_merged_peps_random_10mers$Allele == list_alleles[i,1] & melt_merged_peps_random_10mers$Score >= 0.5,3:4]))
  nbin <- c( 40, 40)                      # 1600 bins #compare binning to umap -> more bins(?)
  bins <- bin2(x, ab, nbin)
  bin_sol <- bins$nc
  sum(bin_sol)
  dim(x)[1]
  bin_sol1 <- bin_sol/sum(bin_sol)
  
  y <- as.matrix(unique(melt_merged_peps_random_10mers[melt_merged_peps_random_10mers$Allele == list_alleles[i,1] & melt_merged_peps_random_10mers$Score < 0.5,3:4]))
  bins <- bin2(y, ab, nbin)
  bin_sol <- bins$nc
  sum(bin_sol)
  dim(y)[1]
  bin_sol2 <- bin_sol/sum(bin_sol)
  diff_bins <- bin_sol1 - bin_sol2
  colnames(diff_bins) <- 1:40
  rownames(diff_bins) <- 1:40
  
  melted_bins <- melt(diff_bins, na.rm = TRUE)
  melted_bins$value <- melted_bins$value*100
  melted_bins <- tibble(melted_bins)
  melted_bins$Length <- "10mer"
  melted_bins$Allele <- list_alleles[i,1]
  colnames(melted_bins)[1:3] <- c("X", "Y", "Diff")
  allele_diff_bins_master_10mers <- rbind(allele_diff_bins_master_10mers, melted_bins)
}
gc()
allele_diff_bins_master_10mers_filt <- allele_diff_bins_master_10mers
rm(allele_diff_bins_master_10mers)
allele_diff_bins_master_10mers_filt[allele_diff_bins_master_10mers_filt == 0] <- NA
allele_diff_bins_master_10mers_filt<-allele_diff_bins_master_10mers_filt[complete.cases(allele_diff_bins_master_10mers_filt),]

allele_diff_bins_master_10mers_filt <- merge(allele_diff_bins_master_10mers_filt, list_alleles, by="Allele")
temp <- allele_diff_bins_master_10mers_filt
temp$Label <- paste(temp$X, temp$Y, sep="-")
library(tidyr)
temp <- temp %>%
  group_by(Label) %>%
  slice(which.max(Diff))

library(ggplot2)
temp <- dplyr::filter(temp, Diff>0.2)
temp <- temp[order(-temp$Diff),]
ggplot(temp, aes(x=Y, y=X, alpha=Diff, fill=Supertype)) + 
  geom_tile()  +theme_minimal() +guides(alpha="none")+expand_limits(x = 0, y = 0)+
  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ggtitle("Binders Enrichment")

allele_diff_bins_master_10mers_filt[which.max(allele_diff_bins_master_10mers_filt$Diff),]
##Example B4402 highest 10mer

x <- as.matrix(unique(melt_merged_peps_random_10mers[melt_merged_peps_random_10mers$Allele == "B4402" & melt_merged_peps_random_10mers$Score >= 0.5,3:4]))
ab <- matrix( c(min(melt_merged_peps_random_10mers$X),min(melt_merged_peps_random_10mers$Y),max(melt_merged_peps_random_10mers$X),max(melt_merged_peps_random_10mers$Y)), 2, 2)      
nbin <- c( 40, 40)                      # 1600 bins #compare binning to umap -> more bins(?)
bins <- bin2(x, ab, nbin)
bin_sol <- bins$nc
sum(bin_sol)
dim(x)[1]
bin_sol1 <- bin_sol/sum(bin_sol)

y <- as.matrix(unique(melt_merged_peps_random_10mers[melt_merged_peps_random_10mers$Allele == "B4402" & melt_merged_peps_random_10mers$Score < 0.5,3:4]))
bins <- bin2(y, ab, nbin)
bin_sol <- bins$nc
sum(bin_sol)
dim(y)[1]
bin_sol2 <- bin_sol/sum(bin_sol)
diff_bins <- bin_sol1 - bin_sol2
colnames(diff_bins) <- 1:40
rownames(diff_bins) <- 1:40

melted_bins <- melt(diff_bins, na.rm = TRUE)
melted_bins$value <- melted_bins$value*100
title <- "B4402 Binder vs Nonbinder"

ggplot(data = melted_bins, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Enrichment") +
  theme_minimal()+ 
  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ggtitle(title)

rm(melt_merged_peps_random_10mers)
gc()
###

#11mers
merged_peps_random_11mers <- merge(merged_peps_random_all, merged_peps_umap11, by="Peptide")
rm(merged_peps_random_all)
melt_merged_peps_random_11mers <- melt(merged_peps_random_11mers, id=c("Peptide", "Allele", "X", "Y"))
rm(merged_peps_random_11mers)
gc()

i <- 1
allele_diff_bins_master_11mers <- tibble(X=character(), Y=character(), Diff=numeric(), Length=character(), Allele=character())
colnames(melt_merged_peps_random_11mers)[5:6] <- c("Tool", "Score")

#fix kmer length as var
ab <- matrix( c(min(melt_merged_peps_random_11mers$X),min(melt_merged_peps_random_11mers$Y),max(melt_merged_peps_random_11mers$X),max(melt_merged_peps_random_11mers$Y)), 2, 2)      
library(ash)
#library(ggplot2)
for (i in 1:length(list_alleles[,1]))
{
  
  x <- as.matrix(unique(melt_merged_peps_random_11mers[melt_merged_peps_random_11mers$Allele == list_alleles[i,1] & melt_merged_peps_random_11mers$Score >= 0.5,3:4]))
  nbin <- c( 40, 40)                      # 1600 bins #compare binning to umap -> more bins(?)
  bins <- bin2(x, ab, nbin)
  bin_sol <- bins$nc
  sum(bin_sol)
  dim(x)[1]
  bin_sol1 <- bin_sol/sum(bin_sol)
  
  y <- as.matrix(unique(melt_merged_peps_random_11mers[melt_merged_peps_random_11mers$Allele == list_alleles[i,1] & melt_merged_peps_random_11mers$Score < 0.5,3:4]))
  bins <- bin2(y, ab, nbin)
  bin_sol <- bins$nc
  sum(bin_sol)
  dim(y)[1]
  bin_sol2 <- bin_sol/sum(bin_sol)
  diff_bins <- bin_sol1 - bin_sol2
  colnames(diff_bins) <- 1:40
  rownames(diff_bins) <- 1:40
  
  melted_bins <- melt(diff_bins, na.rm = TRUE)
  melted_bins$value <- melted_bins$value*100
  melted_bins <- tibble(melted_bins)
  melted_bins$Length <- "11mer"
  melted_bins$Allele <- list_alleles[i,1]
  colnames(melted_bins)[1:3] <- c("X", "Y", "Diff")
  allele_diff_bins_master_11mers <- rbind(allele_diff_bins_master_11mers, melted_bins)
}
gc()
allele_diff_bins_master_11mers_filt <- allele_diff_bins_master_11mers
rm(allele_diff_bins_master_11mers)
allele_diff_bins_master_11mers_filt[allele_diff_bins_master_11mers_filt == 0] <- NA
allele_diff_bins_master_11mers_filt<-allele_diff_bins_master_11mers_filt[complete.cases(allele_diff_bins_master_11mers_filt),]

allele_diff_bins_master_11mers_filt <- merge(allele_diff_bins_master_11mers_filt, list_alleles, by="Allele")
temp <- allele_diff_bins_master_11mers_filt
temp$Label <- paste(temp$X, temp$Y, sep="-")
library(tidyr)
temp <- temp %>%
  group_by(Label) %>%
  slice(which.max(Diff))

library(ggplot2)
temp <- dplyr::filter(temp, Diff>0.2)
temp <- temp[order(-temp$Diff),]
ggplot(temp, aes(x=Y, y=X, alpha=Diff, fill=Supertype)) + 
  geom_tile()  +theme_minimal() +guides(alpha="none")+expand_limits(x = 0, y = 0)+
  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ggtitle("Binders Enrichment")

allele_diff_bins_master_11mers_filt[which.max(allele_diff_bins_master_11mers_filt$Diff),]
##Example B01101 highest 11mer

x <- as.matrix(unique(melt_merged_peps_random_11mers[melt_merged_peps_random_11mers$Allele == "A0101" & melt_merged_peps_random_11mers$Score >= 0.5,3:4]))
ab <- matrix( c(min(melt_merged_peps_random_11mers$X),min(melt_merged_peps_random_11mers$Y),max(melt_merged_peps_random_11mers$X),max(melt_merged_peps_random_11mers$Y)), 2, 2)      
nbin <- c( 40, 40)                      # 1600 bins #compare binning to umap -> more bins(?)
bins <- bin2(x, ab, nbin)
bin_sol <- bins$nc
sum(bin_sol)
dim(x)[1]
bin_sol1 <- bin_sol/sum(bin_sol)

y <- as.matrix(unique(melt_merged_peps_random_11mers[melt_merged_peps_random_11mers$Allele == "A0101" & melt_merged_peps_random_11mers$Score < 0.5,3:4]))
bins <- bin2(y, ab, nbin)
bin_sol <- bins$nc
sum(bin_sol)
dim(y)[1]
bin_sol2 <- bin_sol/sum(bin_sol)
diff_bins <- bin_sol1 - bin_sol2
colnames(diff_bins) <- 1:40
rownames(diff_bins) <- 1:40

melted_bins <- melt(diff_bins, na.rm = TRUE)
melted_bins$value <- melted_bins$value*100
title <- "A0101 Binder vs Nonbinder"

ggplot(data = melted_bins, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Enrichment") +
  theme_minimal()+ 
  coord_fixed() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ggtitle(title)


rm(melt_merged_peps_random_11mers)
gc()
###
