# correlation between clr-genus abundance and amp expression
library(tibble)
library(tidyverse)
library(plyr)
library(dplyr)
library(reshape2)
library(phyloseq)
library(ape)
library(ggplot2)
library(vegan)

#-------------
clean_names <- function(data) {
  names(data) <- gsub("^X", "", names(data))
  names(data) <- gsub("\\.", "-", names(data))
  row.names(data) <- gsub("^X", "", row.names(data))
  row.names(data) <- gsub("\\.", "-", row.names(data))
  if ("SampleID" %in% colnames(data)) {
    data$SampleID <- gsub('^[X]', '', data$SampleID) 
    data$SampleID <- gsub('[.]', '-', data$SampleID) 
  }
  return(data)
}
#---------------

otu_table_in <- read.csv("feature-table.txt", sep = "\t")
otu_table_t <- setNames(data.frame(t(otu_table_in[,-1])), otu_table_in[,1])
otu_table_t <- tibble::rownames_to_column(otu_table_t)

names(otu_table_t)[names(otu_table_t) == "rowname"] <- "SampleID"

otu_table_t <- clean_names(otu_table_t)

# Split the values in the 'SampleID' column based on '_'
split_values <- strsplit(as.character(otu_table_t$SampleID), "_")

# Keep only the first part of the split values
otu_table_t$SampleID <- sapply(split_values, function(x) x[1])
otu_table_t$SampleID <- gsub("-content", "-C", otu_table_t$SampleID)
#---------------------------------------------------------------------------------------------------
taxonomy_in <- read.csv("taxonomy.txt", sep = "\t", row.names = 1, header=T)
#---------------------------------------------------------------------------------------------------
metadata_all <- read.table("Metadata_corrected_analysis.txt", sep="\t", row.names = 1, header=T)
metadata_all <- tibble::rownames_to_column(metadata_all)
names(metadata_all)[names(metadata_all) == "rowname"] <- "SampleID"

# Split the values in the 'SampleID' column based on '_'
split_values <- strsplit(as.character(metadata_all$SampleID), "_")

# Keep only the first part of the split values
metadata_all$SampleID <- sapply(split_values, function(x) x[1])
metadata_all$SampleID <- gsub("-content", "-C", metadata_all$SampleID)
#---------------------------------------------------------------------------------------------------
amp_exp <- read.table(file="AMP_exp_data.txt", sep = "\t", header = T, row.names = 1)
amp_exp_df <- as.data.frame(t(amp_exp)) 

amp_exp_df <- clean_names(amp_exp_df)

# calculate the relative abundance
amp_exp_df <- as.data.frame(t(amp_exp_df))
amp_exp_df_relative <- apply(amp_exp_df, 1, function(x) x/sum(x))

amp_exp_df_rl <- as.matrix(t(amp_exp_df_relative))
# Replace '-MAB' with '-C' in the row names
rownames(amp_exp_df_rl) <- gsub("-MAB", "-C", rownames(amp_exp_df_rl))
#---------------------------------------------------------------------------------------------------
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
    numeric_columns2 <- numeric_columns2[, colSums(numeric_columns2) >= 10]
    # Identify columns with more than one value greater than 0
    selected_cols <- sapply(numeric_columns2, function(col) sum(col > 0) > 1)
    
    # Subset the dataframe to keep only selected columns
    numeric_columns2 <- numeric_columns2[, selected_cols]
    
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
#----------------------------------------------------------------------------------------------------    
    otu_table_final <- gut_Final
    
    row.names(otu_table_final) <- otu_table_final$SampleID
    otu_table_final <- otu_table_final[, -1]
    metadata <- gut_Metadata
    row.names(metadata) <- metadata$SampleID
    metadata <- metadata[, -1]
    
    taxonomy_check <- taxonomy_in[rownames(taxonomy_in) %in% rownames(otu_table_final), , drop = FALSE]
    taxonomy <- as.matrix(taxonomy_check)
    
    colnames(otu_table_final) <- gsub("-MAB", "-C", colnames(otu_table_final))
    rownames(metadata) <- gsub("-MAB", "-C", rownames(metadata))
    # Sanity checks for consistent OTU names
    counts  <- otu_table_final  # Abundance table (e.g. ASV data; to assay data)
    tax     <- taxonomy     # Taxonomy table (to rowData)
    samples <- metadata  # collate data (to colData)
    # Same sample names
    counts <- clean_names(counts)
    counts[, -1] <- apply(counts[, -1], 2, function(x) as.numeric(x))
    #counts <- apply(counts, 2, function(x) as.numeric(x))
    counts <- as.matrix(counts)  
    se <- SummarizedExperiment(assays = list(counts = counts),
                               colData = samples,
                               rowData = tax)
    tse <- as(se, "TreeSummarizedExperiment")
    tse <- transformAssay(se, method = "clr", pseudocount = 1)
    clr_assay <- assays(tse)$clr
    #clr_assay <- t(clr_assay)
    clr_values <- as.data.frame(t(clr_assay))
    
    # Find common row names
    common_rows_final <- intersect(rownames(clr_values), rownames(amp_exp_df_rl))
    
    # Subset both dataframes using common row names
    subset_clr_values_check <- as.matrix(t(clr_values[common_rows_final, , drop = FALSE]))
    # Make sure rownames are turned into a column in both dataframes
    subset_clr_values_check <- data.frame(row.names = row.names(subset_clr_values_check), subset_clr_values_check)
    taxonomy_check <- data.frame(row.names = row.names(taxonomy_check), taxonomy_check)
    
    # Merge the dataframes based on row names
    merged_df_check <- merge(subset_clr_values_check, taxonomy_check, by = "row.names", all = TRUE)
    
    # Rename the first column to "Row.names" (optional)
    colnames(merged_df_check)[1] <- "Row.names"
    # Remove the redundant column containing row names
    merged_df_check <- merged_df_check[, -1]
    
    # Set the 'Taxon' column as row names
    rownames(merged_df_check) <- merged_df_check$Taxon
    
    # Remove the 'Taxon' column from the dataframe
    merged_df_check <- merged_df_check[, -which(names(merged_df_check) == 'Taxon')]
    merged_df_check <- clean_names(merged_df_check)
    
    
    subset_clr_values <- as.matrix(t(merged_df_check))
    #subset_clr_values <- as.matrix(clr_values[common_rows_final, , drop = FALSE])
    subset_amp_exp_df_rl <- as.matrix(amp_exp_df_rl[common_rows_final, , drop = FALSE])
    
    
    out <- corr.test(subset_amp_exp_df_rl, subset_clr_values, use = "pairwise",method="spearman",adjust="BH", alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
    R_value_all <- out$r
    P_value_all <- out$p
    P_adj_all <- out$p.adj
    file_name5 <- paste(level, level2, "_Genera_AMP_Corr.txt", sep = "_")
    write.table(R_value_all, file= file_name5, sep='\t', quote=F)
    
    file_name6 <- paste(level, level2, "_Genera_AMP_Pval.txt", sep = "_")
    write.table(P_value_all, file= file_name6, sep='\t', quote=F)
    
    file_name7 <- paste(level, level2, "_Genera_AMP_AdjPval.txt", sep = "_")
    write.table(P_adj_all, file= file_name7, sep='\t', quote=F)
    # Iterate through rows and columns to find values less than 0.01
    file_conn <- file(paste(level, level2, "_Genera_AMP_signif_corr.txt", sep = "_"), "w")
    
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
             filename=paste(level, level2, "_Taxa_AMP_expression_heatmap_spearman.pdf", sep = "_"))

  }
}
