# Alpha Beta Diversity analysis by diet
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
#Colors <- c("lightblue", "gray") #colors used for gender

Colors1 <- c("darkolivegreen","red4") #colors used for different diets
#Colors1 <- c("midnightblue", "darkgray") #colors used for gender
#-------------------------------------------------------

otu_table_in <- read.csv("selected_feature-table.txt", sep = "\t")
otu_table_t <- setNames(data.frame(t(otu_table_in[,-1])), otu_table_in[,1])
otu_table_t <- tibble::rownames_to_column(otu_table_t)

names(otu_table_t)[names(otu_table_t) == "rowname"] <- "SampleID"

otu_table_t$SampleID <- gsub('^[X]', '', otu_table_t$SampleID) 
otu_table_t$SampleID <- gsub('[.]', '-', otu_table_t$SampleID)

metadata_all <- read.table("Metadata_corrected_analysis.txt", sep="\t", row.names = 1, header=T)
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
  
  Si8_FP001_Final1 <- as.data.frame(t(filtered_df[,1:(ncol(filtered_df)-9)]))
  Si8_FP001_Metadata1 <- filtered_df[,(ncol(filtered_df)-8):ncol(filtered_df)]
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
    gut_Final1 <- as.data.frame(t(filtered_df2[,1:(ncol(filtered_df2)-9)]))
    gut_Metadata1 <- filtered_df2[,(ncol(filtered_df2)-8):ncol(filtered_df2)]
    gut_Metadata <- rownames_to_column(gut_Metadata1, var = "RowID")
    colnames(gut_Metadata)[1] <- "SampleID"
    gut_Final <- rownames_to_column(gut_Final1, var = "RowID")
    colnames(gut_Final)[1] <- "SampleID"
    file_name3 <- paste(level, level2, "_feature-table.txt", sep = "_")
    file_name4 <- paste(level, level2, "_Metadata.txt", sep = "_")
    write.table(gut_Final, file = file_name3, sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(gut_Metadata, file = file_name4, sep = "\t", quote = FALSE, row.names = FALSE)
    
    otu_table_final <- gut_Final
    row.names(otu_table_final) <- otu_table_final$SampleID
    otu_table_final <- otu_table_final[, -1]
    metadata <- gut_Metadata
    row.names(metadata) <- metadata$SampleID
    metadata <- metadata[, -1]
    taxonomy_all <- read.csv("selected_taxonomy.tsv", sep = "\t", row.names = 1)
    taxonomy <- taxonomy_all[rownames(taxonomy_all)  %in% gsub('^[X]', '', row.names(otu_table_final)), ]
    taxonomy <- as.matrix(taxonomy)
    # Sanity checks for consistent OTU names
    counts  <- otu_table_final  # Abundance table (e.g. ASV data; to assay data)
    tax     <- taxonomy     # Taxonomy table (to rowData)
    samples <- metadata  # collate data (to colData)
    # Same sample names
    
    names(counts) <- gsub("^X", "", names(counts))
    names(counts) <- gsub("\\.", "-", names(counts))
    row.names(counts) <- gsub("^X", "", row.names(counts))
    row.names(counts) <- gsub("\\.", "-", row.names(counts))
    counts[, -1] <- apply(counts[, -1], 2, function(x) as.numeric(x))
    #counts <- apply(counts, 2, function(x) as.numeric(x))
    counts <- as.matrix(counts)  
    se <- SummarizedExperiment(assays = list(counts = counts),
                               colData = samples,
                               rowData = tax)
    tse <- as(se, "TreeSummarizedExperiment")
    tse <- transformAssay(se, method = "clr", pseudocount = 1)
    clr_assay <- assays(tse)$clr
    clr_assay <- t(clr_assay)
    
    euclidean_dist <- vegan::vegdist(clr_assay, method = "euclidean")
    
    #UniFrac_distances <- UniFrac(Bushman2, weighted=TRUE)
    UniFrac_distances <- as.matrix(euclidean_dist)
    UniFrac_dist_column <- melt(UniFrac_distances)
    write.table (UniFrac_dist_column, file = "UniFrac_distances.txt", sep = "\t")
    
    Uni_pcoa <- pcoa(UniFrac_distances)
    Uni_pcoa$values[1:2,]
    mds.var.per = round(Uni_pcoa$values$Eigenvalues/sum(Uni_pcoa$values$Eigenvalues)*100, 1)
    pc <- c(1,2)
    file_name1 <- paste(level, level2, "Eucledian.jpg", sep = "_")
    jpeg(file_name1, height = 10, width = 10, units = 'in', res = 600)
    plot(Uni_pcoa$vectors[,1:2], bg= c("darkgreen", "darkolivegreen4", "red", "salmon4")[as.factor(as.data.frame(samples)$diet_sex)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
    #text(Uni_pcoa$vectors[,1:2], labels=rownames(Bushman2@sam_data), cex=0.3, font=1, pos=1)
    ordiellipse(Uni_pcoa$vectors[,1:2], as.factor(as.data.frame(samples)$diet_sex), kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = c("darkgreen", "darkolivegreen4", "red", "salmon4"))
    #ordispider(Uni_pcoa$vectors[,1:2], as.data.frame(samples)$Status, lty=3, spider ="centroid", lwd=1, col="black")
    legend("bottomleft", legend = levels(as.factor(as.data.frame(samples)$diet_sex)), col = c("darkgreen", "darkolivegreen4", "red", "salmon4"), lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
    abline(h=0, v=0, col = "gray60")
    dev.off ()
    adonis2(UniFrac_distances ~ samples$diet_sex)
  }
}
