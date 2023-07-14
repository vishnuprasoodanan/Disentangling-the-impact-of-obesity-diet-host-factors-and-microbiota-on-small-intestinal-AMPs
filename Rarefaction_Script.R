library("GUniFrac")
OTU_table <- read.table(file = "selected_feature-table.txt", sep = "\t", header = T, row.names = 1) #ASV-IDs in rows and Samples in columns
OTU_table <- as.data.frame(t(OTU_table)) #Samples in rows and ASV-IDs in columns
otu.tab.rff <- Rarefy(OTU_table, depth = min(rowSums(OTU_table)))$otu.tab.rff # Samples in Row and OTUs in column
OTU_table_rar <- data.frame(otu.tab.rff)

write.table(OTU_table_rar, file = "Selected_FeatureTable_rarefied.txt", quote = FALSE, sep = '\t')
