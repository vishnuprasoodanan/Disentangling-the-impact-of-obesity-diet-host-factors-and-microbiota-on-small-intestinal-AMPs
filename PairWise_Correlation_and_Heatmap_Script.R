library("ccrepe")
library("reshape")
library("reshape2")
library("dplyr")
library("knitr")
library("foreach")
library("NMF")
#vignette("nested")
library("psych")

Rel_abun_asv <- read.table(file="FP001_AMP_exp.txt", sep = "\t", header = T, row.names = 1)

df <- as.data.frame(t(Rel_abun_asv)) 
names(df) <- gsub("^X", "", names(df))
names(df) <- gsub("\\.", "-", names(df))
# calculate the relative abundance
df <- as.data.frame(t(df))
df_relative <- apply(df, 1, function(x) x/sum(x))

Rel_abun_asv <- as.matrix(t(df_relative))

Rel_abun_asv2 <- read.table(file="FP001_Genus_abundance.txt", sep = "\t", header = T, row.names = 1)
#Rel_abun_asv2 <- as.matrix(Rel_abun_asv2)
df2 <- as.data.frame(Rel_abun_asv2) 
names(df2) <- gsub("^X", "", names(df2))
names(df2) <- gsub("\\.", "-", names(df2))
names(df2) <- gsub("MAB", "C", names(df2))
names(df2) <- gsub("_.*$", "", colnames(df2))
df2 <- as.data.frame(t(df2))
df2 <- subset(df2, rownames(df2) != "53-Si8-C")
Rel_abun_asv2 <- df2 %>% select(where(~ sum(. > 0) >= 3))
Rel_abun_asv2 <- as.matrix(Rel_abun_asv2)

out <- corr.test(Rel_abun_asv, Rel_abun_asv2, use = "pairwise",method="spearman",adjust="BH", alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
R_value_all <- out$r
P_value_all <- out$p
P_adj_all <- out$p.adj
write.table(R_value_all, file= "Genera_AMP_Corr_All.txt", sep='\t', quote=F)
write.table(P_value_all, file= "Genera_AMP_Pval_All.txt", sep='\t', quote=F)
write.table(P_adj_all, file= "Genera_AMP_AdjPval_All.txt", sep='\t', quote=F)
# Iterate through rows and columns to find values less than 0.01
file_conn <- file("significant_correlations.txt", "w")

for (i in 1:nrow(P_value_all)) {
  for (j in 1:ncol(P_value_all)) {
    if (P_value_all[i, j] < 0.01) {
      out = cat(row.names(P_value_all)[i], "\t", colnames(P_value_all)[j], "\t", R_value_all[row.names(P_value_all)[i], colnames(P_value_all)[j]],"\t", P_value_all[row.names(P_value_all)[i], colnames(P_value_all)[j]], "\t", P_adj_all[row.names(P_value_all)[i], colnames(P_value_all)[j]],"\n", file = file_conn)
    }
  }
}
close(file_conn)

#----------------- Heat-map using correlation values
library(dplyr)
library(NMF)
library(RColorBrewer)
result_matrix <- R_value_all
aheatmap(result_matrix, color = "-BrBG:50", breaks = 0, cellwidth = 10, cellheight =10, border_color = "white",
         main = "Taxa-AMP expression",
         distfun = "spearman",
         hclustfun = "complete",
         fontsize=10,
         filename="Taxa_AMP_expression_heatmap_spearman.pdf")

#----------------- Another Method
foreach(i = 1:ncol(Rel_abun_asv),.combine = rbind)%do%{
  foreach(j = 1:ncol(Rel_abun_asv2),.combine = cbind)%do%{
    cor.temp = cor.test(Rel_abun_asv[,i],Rel_abun_asv2[,j],method = "spearman")
    spearmanR = cor.temp$estimate
    spearmanR
  }
}-> R_value
rownames(R_value) <- colnames(Rel_abun_asv)
colnames(R_value) <- colnames(Rel_abun_asv2)
R_value_matrix <- as.matrix(R_value)

foreach(i = 1:ncol(Rel_abun_asv),.combine = rbind)%do%{
  foreach(j = 1:ncol(Rel_abun_asv2),.combine = cbind)%do%{
    cor.temp = cor.test(Rel_abun_asv[,i],Rel_abun_asv2[,j],method = "spearman")
    spearmanP = cor.temp$p.value
    spearmanP
  }
}-> P_value
rownames(P_value) <- colnames(Rel_abun_asv)
colnames(P_value) <- colnames(Rel_abun_asv2)
P_value_matrix <- as.matrix(P_value)

write.table(R_value_matrix, file= "Genera_AMP_association.txt", sep='\t', quote=F)
write.table(R_value_matrix, file= "Genera_AMP_association_p-value.txt", sep='\t', quote=F)
