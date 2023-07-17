#This Script is written by Vishnu Prasoodanan P K
#The input for this script will be the ASV-count table (feature-table) and the metadatafile

library(tibble)
library(tidyverse)
library(plyr)
library(dplyr)
library(reshape2)
library(phyloseq)
library (ape)
library (ggplot2)
library(vegan)

#--------------------------------------------------------
Colors <- c("darkolivegreen4", "red") #colors used for different diets
sex_Colors <- c("lightblue", "gray") #colors used for gender

Colors1 <- c("darkolivegreen","red4") #colors used for different diets
sex_Colors1 <- c("midnightblue", "darkgray") #colors used for gender
#-------------------------------------------------------

df <- read.table(file="selected_feature-table.txt", sep = "\t", header = T, row.names = 1)
names(df) <- gsub("^X", "", names(df))
names(df) <- gsub("\\.", "-", names(df))
# calculate the relative abundance
df <- as.data.frame(t(df))
df_relative <- t(apply(df, 1, function(x) x/sum(x)))

# transpose the result back to its original orientation
df_relative <- as.data.frame(t(df_relative))
write.table(df_relative, file = "sel_feature_table_abundance.txt", sep = '\t', quote = FALSE, row.names = TRUE)

#-------------------------

#otu_table_in <- read.csv("sel_feature_table_abundance.txt", sep = "\t")
otu_table_in <- df_relative
otu_table_t <- as.data.frame(t(otu_table_in ))

#otu_table_t <- setNames(data.frame(t(otu_table_in[,-1])), otu_table_in[,1])
otu_table_t <- tibble::rownames_to_column(otu_table_t)

names(otu_table_t)[names(otu_table_t) == "rowname"] <- "SampleID"

otu_table_t$SampleID <- gsub('^[X]', '', otu_table_t$SampleID) 
otu_table_t$SampleID <- gsub('[.]', '-', otu_table_t$SampleID)

metadata_all <- read.table("selected_Metadata.tsv", sep="\t", row.names = 1, header=T)
metadata_all <- tibble::rownames_to_column(metadata_all)
names(metadata_all)[names(metadata_all) == "rowname"] <- "SampleID"

merge_otu_table_t <- merge(otu_table_t, metadata_all, by.x = "SampleID")
merge_otu_table_t$experiment <- as.factor(merge_otu_table_t$experiment)
merge_otu_table_t$source <- as.factor(merge_otu_table_t$source)
# Get the factor column name
factor_column <- "experiment"

# Iterate through each factor variable
factor_levels <- levels(merge_otu_table_t[[factor_column]])
for (level in factor_levels) {
  # Filter the data frame for the current factor level
  Si8_FP001 <- merge_otu_table_t[merge_otu_table_t[[factor_column]] == level, ]
  numeric_columns <- Si8_FP001[, sapply(Si8_FP001, is.numeric)]
  numeric_columns <- numeric_columns[, colSums(numeric_columns) > 0]
  
  non_numeric_columns <- Si8_FP001[, !sapply(Si8_FP001, is.numeric)]
  
  # Create a new data frame
  filtered_df <- data.frame(numeric_columns, non_numeric_columns)
  row.names(filtered_df) <- Si8_FP001$SampleID
  
  Si8_FP001_Final1 <- as.data.frame(t(filtered_df[,1:(ncol(filtered_df)-8)]))
  Si8_FP001_Metadata1 <- filtered_df[,(ncol(filtered_df)-7):ncol(filtered_df)]
  Si8_FP001_Metadata <- rownames_to_column(Si8_FP001_Metadata1, var = "RowID")
  colnames(Si8_FP001_Metadata)[1] <- "SampleID"
  Si8_FP001_Final <- rownames_to_column(Si8_FP001_Final1, var = "RowID")
  colnames(Si8_FP001_Final)[1] <- "SampleID"
  file_name1 <- paste(level, "_feature-table.txt", sep = "_")
  file_name2 <- paste(level, "_Metadata.txt", sep = "_")
  write.table(Si8_FP001_Final, file = file_name1, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(Si8_FP001_Metadata, file = file_name2, sep = "\t", quote = FALSE, row.names = FALSE)
  filtered_df$source <- as.factor(filtered_df$source)
  factor_column2 <- "source"
  factor_levels2 <- levels(filtered_df[[factor_column2]])
  for (level2 in factor_levels2) {
    gut <- filtered_df[filtered_df[[factor_column2]] == level2, ]
    numeric_columns2 <- gut[, sapply(gut, is.numeric)]
    numeric_columns2 <- numeric_columns2[, colSums(numeric_columns2) > 0]
    
    non_numeric_columns2 <- gut[, !sapply(gut, is.numeric)]
    # Create a new data frame
    filtered_df2 <- data.frame(numeric_columns2, non_numeric_columns2)
    row.names(filtered_df2) <- gut$SampleID
    gut_Final1 <- as.data.frame(t(filtered_df2[,1:(ncol(filtered_df2)-8)]))
    gut_Metadata1 <- filtered_df2[,(ncol(filtered_df2)-7):ncol(filtered_df2)]
    gut_Metadata <- rownames_to_column(gut_Metadata1, var = "RowID")
    colnames(gut_Metadata)[1] <- "SampleID"
    gut_Final <- rownames_to_column(gut_Final1, var = "RowID")
    colnames(gut_Final)[1] <- "SampleID"
    file_name3 <- paste(level, level2, "_feature-table.txt", sep = "_")
    file_name4 <- paste(level, level2, "_Metadata.txt", sep = "_")
    write.table(gut_Final, file = file_name3, sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(gut_Metadata, file = file_name4, sep = "\t", quote = FALSE, row.names = FALSE)
    
    gut_Final_v1 <- as.data.frame(t(gut_Final1))
    gut_Final_v1 <- tibble::rownames_to_column( gut_Final_v1)
    names(gut_Final_v1)[names(gut_Final_v1) == "rowname"] <- "SampleID"
    gut_Metadata_v1 <- gut_Metadata
    gut_Metadata_v2 <- gut_Metadata_v1[, -c(1, 2, 4:7)]
    gut_Final_v2 <- merge(gut_Final_v1, gut_Metadata_v2, by.x = "SampleID")
    
    gut_Final_diet <- gut_Final_v2 %>% select(SampleID, diet, everything())%>% select(-sex)
    row.names(gut_Final_diet) <- gut_Final_diet$SampleID
    gut_Final_diet <- gut_Final_diet[, -1]
    gut_Final_diet_t <- as.data.frame(t(gut_Final_diet))
    
    gut_Final_diet_numeric <- gut_Final_diet
    colnames(gut_Final_diet_numeric) <- NULL
    rownames(gut_Final_diet_numeric) <- NULL
    gut_Final_diet_numeric <- gut_Final_diet_numeric[, -1]
    gut_Final_diet_numeric_t <- as.data.frame(t(gut_Final_diet_numeric))
    
    gut_Final_sex <- gut_Final_v2 %>% select(SampleID, sex, everything())%>% select(-diet)
    row.names(gut_Final_sex) <- gut_Final_sex$SampleID
    gut_Final_sex <- gut_Final_sex[, -1]
    gut_Final_sex_t <- as.data.frame(t(gut_Final_sex))
    
    gut_Final_sex_numeric <- gut_Final_sex
    colnames(gut_Final_sex_numeric) <- NULL
    rownames(gut_Final_sex_numeric) <- NULL
    gut_Final_sex_numeric <- gut_Final_sex_numeric[, -1]
    gut_Final_sex_numeric_t <- as.data.frame(t(gut_Final_sex_numeric))
    
    KO <- apply(gut_Final_diet_numeric_t, 2, as.numeric)
    #KO <- read.table(file = "pcoa_numeric_data.txt", sep = '\t', colClasses = "numeric") # data without row and column names
    
    #KO_test <- read.table(file = "pcoa_whole_data.txt", sep = '\t', header = TRUE, row.names = 1)
    KO_test <- gut_Final_diet_t
    KO_test_T <- as.data.frame(t(KO_test))
    KO_test_1 <- as.data.frame(KO_test_T[2:ncol(KO_test_T)])
    rownames(KO) <- colnames(KO_test_1)
    colnames(KO) <- row.names(KO_test_1)
    
    KO_proportions <- KO
    KO_proportions1 <- as.data.frame(t(KO_proportions))
    
    class <- KO_test_T$diet
    Bray_pcoa <-pcoa(vegdist(KO_proportions1, "bray"))
    Bray_pcoa$values[1:2,]
    mds.var.per = round(Bray_pcoa$values$Eigenvalues/sum(Bray_pcoa$values$Eigenvalues)*100, 1)
    Bray_PCoA_MATRIX <- Bray_pcoa$vectors[,1:2]
    Bray_PCoA_MATRIX <- as.data.frame(Bray_PCoA_MATRIX)
    
    Bray_distances <-vegdist(KO_proportions1, "bray")
    adonis2(Bray_distances ~ KO_test_T$diet)
    
    Bray_PCoA_MATRIX_New <- cbind(Bray_PCoA_MATRIX, class)
    write.table(Bray_PCoA_MATRIX_New, file = "PCA_data", quote = FALSE, sep = '\t')
    
    pc <- c(1,2)
    file_name6 <- paste(level, level2, "_PCOA_rel_abundance_diet.jpg", sep = "_")
    jpeg(file_name6, height = 10, width = 10, units = 'in', res = 600)
    
    plot(Bray_pcoa$vectors[,1:2], bg=Colors1[as.factor(Bray_PCoA_MATRIX_New$class)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
    ordiellipse(Bray_pcoa$vectors[,1:2], Bray_PCoA_MATRIX_New$class, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 70, col = Colors1)
    ordispider(Bray_pcoa$vectors[,1:2], Bray_PCoA_MATRIX_New$class, lty=3, spider ="centroid", lwd=1, col="black")
    legend("topright", legend = c("Chow", "WSD"), col = Colors1,lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
    text(Bray_pcoa$vectors[,1:2], labels=as.factor(rownames(KO_proportions1)), cex=0.6, font=1, pos=1)
    abline(h=0, v=0, col = "gray60")
    dev.off ()
    
    KO <- apply(gut_Final_sex_numeric_t, 2, as.numeric)
    KO_test <- gut_Final_sex_t
    KO_test_T <- as.data.frame(t(KO_test))
    KO_test_1 <- as.data.frame(KO_test_T[2:ncol(KO_test_T)])
    rownames(KO) <- colnames(KO_test_1)
    colnames(KO) <- row.names(KO_test_1)
    
    KO_proportions <- KO
    KO_proportions1 <- as.data.frame(t(KO_proportions))
    
    class <- KO_test_T$sex
    Bray_pcoa <-pcoa(vegdist(KO_proportions1, "bray"))
    Bray_pcoa$values[1:2,]
    mds.var.per = round(Bray_pcoa$values$Eigenvalues/sum(Bray_pcoa$values$Eigenvalues)*100, 1)
    Bray_PCoA_MATRIX <- Bray_pcoa$vectors[,1:2]
    Bray_PCoA_MATRIX <- as.data.frame(Bray_PCoA_MATRIX)
    
    Bray_distances <-vegdist(KO_proportions1, "bray")
    adonis2(Bray_distances ~ KO_test_T$sex)
    
    Bray_PCoA_MATRIX_New <- cbind(Bray_PCoA_MATRIX, class)
    write.table(Bray_PCoA_MATRIX_New, file = "PCA_data", quote = FALSE, sep = '\t')
    
    pc <- c(1,2)
    file_name6 <- paste(level, level2, "_PCOA_rel_abundance_sext.jpg", sep = "_")
    jpeg(file_name6, height = 10, width = 10, units = 'in', res = 600)
    
    plot(Bray_pcoa$vectors[,1:2], bg=sex_Colors1[as.factor(Bray_PCoA_MATRIX_New$class)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
    ordiellipse(Bray_pcoa$vectors[,1:2], Bray_PCoA_MATRIX_New$class, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 70, col = sex_Colors1)
    ordispider(Bray_pcoa$vectors[,1:2], Bray_PCoA_MATRIX_New$class, lty=3, spider ="centroid", lwd=1, col="black")
    legend("topright", legend = c("Female", "Male"), col = sex_Colors1,lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
    text(Bray_pcoa$vectors[,1:2], labels=as.factor(rownames(KO_proportions1)), cex=0.6, font=1, pos=1)
    abline(h=0, v=0, col = "gray60")
    dev.off ()
    
    
    }
}
