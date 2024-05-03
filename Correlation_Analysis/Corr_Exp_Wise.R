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
library(NMF)
library(RColorBrewer)
library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(mia)
library(psych)
library(ggpubr)
################################################################################################################
#Experiment
#------------------------------------------------------------------------------------------------------------------
# Load Functions
#Function-1------------------------------------------------------------------------------------------------------------------
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
#Function-2-------------------------------------------------------------------------------------------------------------------
generate_heatmap <- function(result_matrix, filename) {
  aheatmap(result_matrix, color = "-BrBG:50", breaks = 0, cellwidth = 10, cellheight = 10, border_color = "white",
           main = "Taxa-AMP expression",
           distfun = "spearman",
           hclustfun = "complete",
           fontsize = 10,
           filename = filename)
}
#Function-3-------------------------------------------------------------------------------------------------------------------
write_significant_results <- function(p_val_matrix, r_val_matrix, adj_p_val_matrix, file_conn) {
  for (i in 1:nrow(p_val_matrix)) {
    for (j in 1:ncol(p_val_matrix)) {
      if (p_val_matrix[i, j] < 0.01) {
        cat(row.names(p_val_matrix)[i], "\t", colnames(p_val_matrix)[j], "\t", 
            r_val_matrix[row.names(p_val_matrix)[i], colnames(p_val_matrix)[j]], "\t", 
            p_val_matrix[row.names(p_val_matrix)[i], colnames(p_val_matrix)[j]], "\t", 
            adj_p_val_matrix[row.names(p_val_matrix)[i], colnames(p_val_matrix)[j]], "\n", 
            file = file_conn)
      }
    }
  }
}
#Function-5------------------------------------------------------------------------------------------------------------
rel_abundance <- function(data) {
  data_relative <- apply(data, 1, function(x) x/sum(x))
  data_relative_matrix <- as.matrix(t(data_relative))
  return(data_relative_matrix)
}
#Function-6------------------------------------------------------------------------------------------------------------
generate_scatter_plots <- function(amp_data, genus_clr_data, out_file_name) {
  pdf(out_file_name)
  par(mfrow = c(3, 3))
  
  for (col_name_df1 in colnames(amp_data)) {
    # Extract the column from df1
    col_df1 <- amp_data[, col_name_df1]
    
    # Create a new dataframe with the column values and row names from df1
    new_df1 <- data.frame(col_df1, row.names = rownames(amp_data))
    colnames(new_df1) <- col_name_df1
    
    # Loop through the columns of df2
    for (col_name_df2 in colnames(genus_clr_data)) {
      # Extract the column from df2
      col_df2 <- genus_clr_data[, col_name_df2]
      
      # Create a new dataframe with the column values and row names from df2
      new_df2 <- data.frame(col_df2, row.names = rownames(genus_clr_data))
      colnames(new_df2) <- col_name_df2
      
      # Merge new_df1 and new_df2
      merged_lm_df <- cbind(new_df1, new_df2)
      corr_results <- corr.test(merged_lm_df[,1], merged_lm_df[,2], use = "pairwise", method = "spearman", adjust = "BH", alpha = .05, ci = TRUE, minlength = 5, normal = TRUE)
      corr_val <- round(corr_results$r, 3)
      corr_p_val <- round(corr_results$p, 3)
      #y_title <- split_into_n_lines(colnames(merged_lm_df)[2], 3, 70)
      p <- ggplot(merged_lm_df, aes(x = merged_lm_df[,1], y = merged_lm_df[,2])) + 
        geom_point(colour = 'red', size = 4) +
        geom_smooth(method = lm) + 
        labs(x = colnames(merged_lm_df)[1], y =colnames(merged_lm_df)[2]) + 
        annotate("text", x = Inf, y = max(merged_lm_df[,2]), hjust = 1, vjust = 0, size = 4, 
                 label = paste("cor-value = ", corr_val, ", p-value = ", corr_p_val)) +
        theme_minimal() + 
        theme(axis.text = element_text(size = 8), axis.title = element_text(size = 8, face = "bold"))
      
      print(p)
    }
  }
  
  dev.off()
}
#Function7-------------------------------------------------------------------------------------------------------------
subset_matrices <- function(P_value_all, R_value_all, threshold = 0.05) {
  # Find row indices where any value is less than the threshold
  rows_to_keep <- apply(P_value_all, 1, function(row) any(row < threshold))
  
  # Find column indices where any value is less than the threshold
  cols_to_keep <- apply(P_value_all, 2, function(col) any(col < threshold))
  
  # Subset the matrix based on the identified rows and columns
  subset_matrix <- P_value_all[rows_to_keep, cols_to_keep]
  
  # Extract row and column names from the subset matrix
  rows_to_keep <- rownames(subset_matrix)
  cols_to_keep <- colnames(subset_matrix)
  
  # Subset R_value_all using the extracted row and column names
  subset_R <- R_value_all[rows_to_keep, cols_to_keep]
  
  return(list(P_value = subset_matrix, R_value = subset_R))
}
#----------------------------------------------------------------------------------------------------------------------
#code starts here
#----------------------------------------------------------------------------------------------------------------------
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

amp_exp_df <- as.data.frame(t(clean_names(amp_exp_df)))
# Replace '-MAB' with '-C' in the row names
rownames(amp_exp_df) <- gsub("-MAB", "-C", rownames(amp_exp_df))

# calculate the relative abundance
amp_exp_df_rl <- rel_abundance(amp_exp_df)
# Replace '-MAB' with '-C' in the row names
rownames(amp_exp_df_rl) <- gsub("-MAB", "-C", rownames(amp_exp_df_rl))
#---------------------------------------------------------------------------------------------------
merge_otu_table_t <- merge(otu_table_t, metadata_all, by.x = "SampleID")
merge_otu_table_t$experiment <- as.factor(merge_otu_table_t$experiment)
merge_otu_table_t$source <- as.factor(merge_otu_table_t$source)

# Get the factor column name
factor_column <- "experiment"

# Iterate through each factor variable
factor_levels <- levels(as.factor(merge_otu_table_t[[factor_column]]))
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
    #clr-transformation
    se <- SummarizedExperiment(assays = list(counts = counts),
                               colData = samples,
                               rowData = tax)
    
    # Convert to TreeSummarizedExperiment
    tse <- as(se, "TreeSummarizedExperiment")
    
    # Transform assay to clr
    tse <- transformAssay(tse, method = "clr", pseudocount = 1)
    
    # Extract clr assay values
    clr_assay <- assays(tse)$clr
    #clr_assay <- assays(tse)[[1]]
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
    
    # Rename the first column to "Row.names"
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
    amp_exp_df_raw <- as.matrix(amp_exp_df[common_rows_final, , drop = FALSE])
    
    file_name13 <- paste(level, level2, "_relAMP_scatterplot.pdf", sep = "_")
    file_name14 <- paste(level, level2, "_rawAMP_scatterplot.pdf", sep = "_")
    generate_scatter_plots(subset_amp_exp_df_rl, subset_clr_values, file_name13)
    generate_scatter_plots(amp_exp_df_raw, subset_clr_values, file_name14)
    #ggsave(file_name13, p, width = 20, height = 10, limitsize = FALSE)
    #---------------- evaluate correlation using relative expression of AMP
    out <- corr.test(subset_amp_exp_df_rl, subset_clr_values, use = "pairwise",method="spearman",adjust="BH", alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
    R_value_all <- out$r
    P_value_all <- out$p
    P_adj_all <- out$p.adj
    #---------------- evaluate correlation using raw expression of AMP
    out_raw <- corr.test(amp_exp_df_raw, subset_clr_values, use = "pairwise",method="spearman",adjust="BH", alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
    raw_R_value_all <- out_raw$r
    raw_P_value_all <- out_raw$p
    raw_P_adj_all <- out_raw$p.adj
    #---------------- Write Correlation matrices (using rel. expression) in text files 
    file_name5 <- paste(level, level2, "_Genera_REL_AMP_Corr.txt", sep = "_")
    write.table(R_value_all, file= file_name5, sep='\t', quote=F)
    
    file_name6 <- paste(level, level2, "_Genera_REL_AMP_Pval.txt", sep = "_")
    write.table(P_value_all, file= file_name6, sep='\t', quote=F)
    
    file_name7 <- paste(level, level2, "_Genera_REL_AMP_AdjPval.txt", sep = "_")
    write.table(P_adj_all, file= file_name7, sep='\t', quote=F)
    
    #---------------- Write Correlation matrices (using raw expression) in text files 
    file_name8 <- paste(level, level2, "_Genera_RAW_AMP_Corr.txt", sep = "_")
    write.table(raw_R_value_all, file= file_name8, sep='\t', quote=F)
    
    file_name9 <- paste(level, level2, "_Genera_RAW_AMP_Pval.txt", sep = "_")
    write.table(raw_P_value_all, file= file_name9, sep='\t', quote=F)
    
    file_name10 <- paste(level, level2, "_Genera_RAW_AMP_AdjPval.txt", sep = "_")
    write.table(raw_P_adj_all, file= file_name10, sep='\t', quote=F)
    
    #----------------- Iterate through rows and columns to find values less than 0.01
    file_conn <- file(paste(level, level2, "_Genera_REL_AMP_signif_corr.txt", sep = "_"), "w")
    write_significant_results(P_value_all, R_value_all, P_adj_all, file_conn)
    close(file_conn)
    #----------------- Iterate through rows and columns to find values less than 0.01
    file_conn1 <- file(paste(level, level2, "_Genera_RAW_AMP_signif_corr.txt", sep = "_"), "w")
    write_significant_results(raw_P_value_all, raw_R_value_all, raw_P_adj_all, file_conn)
    close(file_conn1)
    #----------------- Heat-map using correlation values
    result_matrix <- R_value_all
    file_name11 <- paste(level, level2, "_Taxa_REL_AMP_expression_heatmap_spearman.pdf", sep = "_")
    generate_heatmap(result_matrix, file_name11)
    
    result_matrix1 <- raw_R_value_all
    file_name12 <- paste(level, level2, "_Taxa_RAW_AMP_expression_heatmap_spearman.pdf", sep = "_")
    generate_heatmap(result_matrix1, file_name12)
    
    
    # Call the function with your matrices
    subset_result <- subset_matrices(raw_P_value_all, raw_R_value_all)
    subset_result1 <- subset_matrices(P_value_all, R_value_all)
    # Extract the subset matrices
    subset_P_value <- subset_result$P_value
    subset_R_value <- subset_result$R_value
    
    # Extract the subset matrices
    subset_P_value1 <- subset_result1$P_value
    subset_R_value1 <- subset_result1$R_value
    
    # Check if result_matrix2 is not empty
    if (nrow(subset_R_value) != 0 && ncol(subset_R_value) != 0) {
      # If result_matrix2 is not empty, proceed with the code
      result_matrix2 <- subset_R_value
      file_name12 <- paste(level, level2, "_SIGNIFICANT_RAW_HEATMAP.pdf", sep = "_")
      generate_heatmap(result_matrix2, file_name12)
    }
    else{
      print("EMPTY MATRIX!!!")
    }
    # Check if result_matrix2 is not empty
    if (nrow(subset_R_value1) != 0 && ncol(subset_R_value1) != 0) {
      # If result_matrix2 is not empty, proceed with the code
      result_matrix3 <- subset_R_value1
      file_name13 <- paste(level, level2, "_SIGNIFICANT_REL_HEATMAP.pdf", sep = "_")
      generate_heatmap(result_matrix3, file_name13)
    }
    else{
      print("EMPTY MATRIX!!!")
    }
  }
}
rm(list = ls())
#####################################################################################################
#Diet
#------------------------------------------------------------------------------------------------------------------
# Load Functions
#Function-1------------------------------------------------------------------------------------------------------------------
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
#Function-2-------------------------------------------------------------------------------------------------------------------
generate_heatmap <- function(result_matrix, filename) {
  aheatmap(result_matrix, color = "-BrBG:50", breaks = 0, cellwidth = 10, cellheight = 10, border_color = "white",
           main = "Taxa-AMP expression",
           distfun = "spearman",
           hclustfun = "complete",
           fontsize = 10,
           filename = filename)
}
#Function-3-------------------------------------------------------------------------------------------------------------------
write_significant_results <- function(p_val_matrix, r_val_matrix, adj_p_val_matrix, file_conn) {
  for (i in 1:nrow(p_val_matrix)) {
    for (j in 1:ncol(p_val_matrix)) {
      if (p_val_matrix[i, j] < 0.01) {
        cat(row.names(p_val_matrix)[i], "\t", colnames(p_val_matrix)[j], "\t",
            r_val_matrix[row.names(p_val_matrix)[i], colnames(p_val_matrix)[j]], "\t",
            p_val_matrix[row.names(p_val_matrix)[i], colnames(p_val_matrix)[j]], "\t",
            adj_p_val_matrix[row.names(p_val_matrix)[i], colnames(p_val_matrix)[j]], "\n",
            file = file_conn)
      }
    }
  }
}
#Function-5------------------------------------------------------------------------------------------------------------
rel_abundance <- function(data) {
  data_relative <- apply(data, 1, function(x) x/sum(x))
  data_relative_matrix <- as.matrix(t(data_relative))
  return(data_relative_matrix)
}
#Function-6------------------------------------------------------------------------------------------------------------
generate_scatter_plots <- function(amp_data, genus_clr_data, out_file_name) {
  pdf(out_file_name)
  par(mfrow = c(3, 3))

  for (col_name_df1 in colnames(amp_data)) {
    # Extract the column from df1
    col_df1 <- amp_data[, col_name_df1]

    # Create a new dataframe with the column values and row names from df1
    new_df1 <- data.frame(col_df1, row.names = rownames(amp_data))
    colnames(new_df1) <- col_name_df1

    # Loop through the columns of df2
    for (col_name_df2 in colnames(genus_clr_data)) {
      # Extract the column from df2
      col_df2 <- genus_clr_data[, col_name_df2]

      # Create a new dataframe with the column values and row names from df2
      new_df2 <- data.frame(col_df2, row.names = rownames(genus_clr_data))
      colnames(new_df2) <- col_name_df2

      # Merge new_df1 and new_df2
      merged_lm_df <- cbind(new_df1, new_df2)
      corr_results <- corr.test(merged_lm_df[,1], merged_lm_df[,2], use = "pairwise", method = "spearman", adjust = "BH", alpha = .05, ci = TRUE, minlength = 5, normal = TRUE)
      corr_val <- round(corr_results$r, 3)
      corr_p_val <- round(corr_results$p, 3)
      #y_title <- split_into_n_lines(colnames(merged_lm_df)[2], 3, 70)
      p <- ggplot(merged_lm_df, aes(x = merged_lm_df[,1], y = merged_lm_df[,2])) +
        geom_point(colour = 'red', size = 4) +
        geom_smooth(method = lm) +
        labs(x = colnames(merged_lm_df)[1], y =colnames(merged_lm_df)[2]) +
        annotate("text", x = Inf, y = max(merged_lm_df[,2]), hjust = 1, vjust = 0, size = 4,
                 label = paste("cor-value = ", corr_val, ", p-value = ", corr_p_val)) +
        theme_minimal() +
        theme(axis.text = element_text(size = 8), axis.title = element_text(size = 8, face = "bold"))

      print(p)
    }
  }

  dev.off()
}
#Function7-------------------------------------------------------------------------------------------------------------
subset_matrices <- function(P_value_all, R_value_all, threshold = 0.05) {
  # Find row indices where any value is less than the threshold
  rows_to_keep <- apply(P_value_all, 1, function(row) any(row < threshold))

  # Find column indices where any value is less than the threshold
  cols_to_keep <- apply(P_value_all, 2, function(col) any(col < threshold))

  # Subset the matrix based on the identified rows and columns
  subset_matrix <- P_value_all[rows_to_keep, cols_to_keep]

  # Extract row and column names from the subset matrix
  rows_to_keep <- rownames(subset_matrix)
  cols_to_keep <- colnames(subset_matrix)

  # Subset R_value_all using the extracted row and column names
  subset_R <- R_value_all[rows_to_keep, cols_to_keep]

  return(list(P_value = subset_matrix, R_value = subset_R))
}
#----------------------------------------------------------------------------------------------------------------------
#code starts here
#----------------------------------------------------------------------------------------------------------------------
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

amp_exp_df <- as.data.frame(t(clean_names(amp_exp_df)))
# Replace '-MAB' with '-C' in the row names
rownames(amp_exp_df) <- gsub("-MAB", "-C", rownames(amp_exp_df))

# calculate the relative abundance
amp_exp_df_rl <- rel_abundance(amp_exp_df)
# Replace '-MAB' with '-C' in the row names
rownames(amp_exp_df_rl) <- gsub("-MAB", "-C", rownames(amp_exp_df_rl))
#---------------------------------------------------------------------------------------------------
merge_otu_table_t <- merge(otu_table_t, metadata_all, by.x = "SampleID")
merge_otu_table_t$experiment <- as.factor(merge_otu_table_t$experiment)
merge_otu_table_t$source <- as.factor(merge_otu_table_t$source)

# Get the factor column name
factor_column <- "diet"

# Iterate through each factor variable
factor_levels <- levels(as.factor(merge_otu_table_t[[factor_column]]))
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
    #clr-transformation
    se <- SummarizedExperiment(assays = list(counts = counts),
                               colData = samples,
                               rowData = tax)

    # Convert to TreeSummarizedExperiment
    tse <- as(se, "TreeSummarizedExperiment")

    # Transform assay to clr
    tse <- transformAssay(tse, method = "clr", pseudocount = 1)

    # Extract clr assay values
    clr_assay <- assays(tse)$clr
    #clr_assay <- assays(tse)[[1]]
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

    # Rename the first column to "Row.names"
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
    amp_exp_df_raw <- as.matrix(amp_exp_df[common_rows_final, , drop = FALSE])

    file_name13 <- paste(level, level2, "_relAMP_scatterplot.pdf", sep = "_")
    file_name14 <- paste(level, level2, "_rawAMP_scatterplot.pdf", sep = "_")
    generate_scatter_plots(subset_amp_exp_df_rl, subset_clr_values, file_name13)
    generate_scatter_plots(amp_exp_df_raw, subset_clr_values, file_name14)
    #ggsave(file_name13, p, width = 20, height = 10, limitsize = FALSE)
    #---------------- evaluate correlation using relative expression of AMP
    out <- corr.test(subset_amp_exp_df_rl, subset_clr_values, use = "pairwise",method="spearman",adjust="BH", alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
    R_value_all <- out$r
    P_value_all <- out$p
    P_adj_all <- out$p.adj
    #---------------- evaluate correlation using raw expression of AMP
    out_raw <- corr.test(amp_exp_df_raw, subset_clr_values, use = "pairwise",method="spearman",adjust="BH", alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
    raw_R_value_all <- out_raw$r
    raw_P_value_all <- out_raw$p
    raw_P_adj_all <- out_raw$p.adj
    #---------------- Write Correlation matrices (using rel. expression) in text files
    file_name5 <- paste(level, level2, "_Genera_REL_AMP_Corr.txt", sep = "_")
    write.table(R_value_all, file= file_name5, sep='\t', quote=F)

    file_name6 <- paste(level, level2, "_Genera_REL_AMP_Pval.txt", sep = "_")
    write.table(P_value_all, file= file_name6, sep='\t', quote=F)

    file_name7 <- paste(level, level2, "_Genera_REL_AMP_AdjPval.txt", sep = "_")
    write.table(P_adj_all, file= file_name7, sep='\t', quote=F)

    #---------------- Write Correlation matrices (using raw expression) in text files
    file_name8 <- paste(level, level2, "_Genera_RAW_AMP_Corr.txt", sep = "_")
    write.table(raw_R_value_all, file= file_name8, sep='\t', quote=F)

    file_name9 <- paste(level, level2, "_Genera_RAW_AMP_Pval.txt", sep = "_")
    write.table(raw_P_value_all, file= file_name9, sep='\t', quote=F)

    file_name10 <- paste(level, level2, "_Genera_RAW_AMP_AdjPval.txt", sep = "_")
    write.table(raw_P_adj_all, file= file_name10, sep='\t', quote=F)

    #----------------- Iterate through rows and columns to find values less than 0.01
    file_conn <- file(paste(level, level2, "_Genera_REL_AMP_signif_corr.txt", sep = "_"), "w")
    write_significant_results(P_value_all, R_value_all, P_adj_all, file_conn)
    close(file_conn)
    #----------------- Iterate through rows and columns to find values less than 0.01
    file_conn1 <- file(paste(level, level2, "_Genera_RAW_AMP_signif_corr.txt", sep = "_"), "w")
    write_significant_results(raw_P_value_all, raw_R_value_all, raw_P_adj_all, file_conn)
    close(file_conn1)
    #----------------- Heat-map using correlation values
    result_matrix <- R_value_all
    file_name11 <- paste(level, level2, "_Taxa_REL_AMP_expression_heatmap_spearman.pdf", sep = "_")
    generate_heatmap(result_matrix, file_name11)

    result_matrix1 <- raw_R_value_all
    file_name12 <- paste(level, level2, "_Taxa_RAW_AMP_expression_heatmap_spearman.pdf", sep = "_")
    generate_heatmap(result_matrix1, file_name12)

    # Call the function with your matrices
    subset_result <- subset_matrices(raw_P_value_all, raw_R_value_all)
    subset_result1 <- subset_matrices(P_value_all, R_value_all)
    # Extract the subset matrices
    subset_P_value <- subset_result$P_value
    subset_R_value <- subset_result$R_value
    
    subset_P_value1 <- subset_result1$P_value
    subset_R_value1 <- subset_result1$R_value

    # Check if result_matrix2 is not empty
    if (nrow(subset_R_value) != 0 && ncol(subset_R_value) != 0) {
      # If result_matrix2 is not empty, proceed with the code
      result_matrix2 <- subset_R_value
      file_name12 <- paste(level, level2, "_SIGNIFICANT_RAW_HEATMAP.pdf", sep = "_")
      generate_heatmap(result_matrix2, file_name12)
    }
    else{
      print("EMPTY MATRIX!!!")
    }
    # Check if result_matrix2 is not empty
    if (nrow(subset_R_value1) != 0 && ncol(subset_R_value1) != 0) {
      # If result_matrix2 is not empty, proceed with the code
      result_matrix3 <- subset_R_value1
      file_name13 <- paste(level, level2, "_SIGNIFICANT_REL_HEATMAP.pdf", sep = "_")
      generate_heatmap(result_matrix3, file_name13)
    }
    else{
      print("EMPTY MATRIX!!!")
    }
  }
}
rm(list = ls())

#####################################################################################################
#Source
#------------------------------------------------------------------------------------------------------------------
# Load Functions
#Function-1------------------------------------------------------------------------------------------------------------------
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
#Function-2-------------------------------------------------------------------------------------------------------------------
generate_heatmap <- function(result_matrix, filename) {
  aheatmap(result_matrix, color = "-BrBG:50", breaks = 0, cellwidth = 10, cellheight = 10, border_color = "white",
           main = "Taxa-AMP expression",
           distfun = "spearman",
           hclustfun = "complete",
           fontsize = 10,
           filename = filename)
}
#Function-3-------------------------------------------------------------------------------------------------------------------
write_significant_results <- function(p_val_matrix, r_val_matrix, adj_p_val_matrix, file_conn) {
  for (i in 1:nrow(p_val_matrix)) {
    for (j in 1:ncol(p_val_matrix)) {
      if (p_val_matrix[i, j] < 0.01) {
        cat(row.names(p_val_matrix)[i], "\t", colnames(p_val_matrix)[j], "\t",
            r_val_matrix[row.names(p_val_matrix)[i], colnames(p_val_matrix)[j]], "\t",
            p_val_matrix[row.names(p_val_matrix)[i], colnames(p_val_matrix)[j]], "\t",
            adj_p_val_matrix[row.names(p_val_matrix)[i], colnames(p_val_matrix)[j]], "\n",
            file = file_conn)
      }
    }
  }
}
#Function-5------------------------------------------------------------------------------------------------------------
rel_abundance <- function(data) {
  data_relative <- apply(data, 1, function(x) x/sum(x))
  data_relative_matrix <- as.matrix(t(data_relative))
  return(data_relative_matrix)
}
#Function-6------------------------------------------------------------------------------------------------------------
generate_scatter_plots <- function(amp_data, genus_clr_data, out_file_name) {
  pdf(out_file_name)
  par(mfrow = c(3, 3))

  for (col_name_df1 in colnames(amp_data)) {
    # Extract the column from df1
    col_df1 <- amp_data[, col_name_df1]

    # Create a new dataframe with the column values and row names from df1
    new_df1 <- data.frame(col_df1, row.names = rownames(amp_data))
    colnames(new_df1) <- col_name_df1

    # Loop through the columns of df2
    for (col_name_df2 in colnames(genus_clr_data)) {
      # Extract the column from df2
      col_df2 <- genus_clr_data[, col_name_df2]

      # Create a new dataframe with the column values and row names from df2
      new_df2 <- data.frame(col_df2, row.names = rownames(genus_clr_data))
      colnames(new_df2) <- col_name_df2

      # Merge new_df1 and new_df2
      merged_lm_df <- cbind(new_df1, new_df2)
      corr_results <- corr.test(merged_lm_df[,1], merged_lm_df[,2], use = "pairwise", method = "spearman", adjust = "BH", alpha = .05, ci = TRUE, minlength = 5, normal = TRUE)
      corr_val <- round(corr_results$r, 3)
      corr_p_val <- round(corr_results$p, 3)
      #y_title <- split_into_n_lines(colnames(merged_lm_df)[2], 3, 70)
      p <- ggplot(merged_lm_df, aes(x = merged_lm_df[,1], y = merged_lm_df[,2])) +
        geom_point(colour = 'red', size = 4) +
        geom_smooth(method = lm) +
        labs(x = colnames(merged_lm_df)[1], y =colnames(merged_lm_df)[2]) +
        annotate("text", x = Inf, y = max(merged_lm_df[,2]), hjust = 1, vjust = 0, size = 4,
                 label = paste("cor-value = ", corr_val, ", p-value = ", corr_p_val)) +
        theme_minimal() +
        theme(axis.text = element_text(size = 8), axis.title = element_text(size = 8, face = "bold"))

      print(p)
    }
  }

  dev.off()
}
#Function7-------------------------------------------------------------------------------------------------------------
subset_matrices <- function(P_value_all, R_value_all, threshold = 0.05) {
  # Find row indices where any value is less than the threshold
  rows_to_keep <- apply(P_value_all, 1, function(row) any(row < threshold))

  # Find column indices where any value is less than the threshold
  cols_to_keep <- apply(P_value_all, 2, function(col) any(col < threshold))

  # Subset the matrix based on the identified rows and columns
  subset_matrix <- P_value_all[rows_to_keep, cols_to_keep]

  # Extract row and column names from the subset matrix
  rows_to_keep <- rownames(subset_matrix)
  cols_to_keep <- colnames(subset_matrix)

  # Subset R_value_all using the extracted row and column names
  subset_R <- R_value_all[rows_to_keep, cols_to_keep]

  return(list(P_value = subset_matrix, R_value = subset_R))
}
#----------------------------------------------------------------------------------------------------------------------
#code starts here
#----------------------------------------------------------------------------------------------------------------------
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

amp_exp_df <- as.data.frame(t(clean_names(amp_exp_df)))
# Replace '-MAB' with '-C' in the row names
rownames(amp_exp_df) <- gsub("-MAB", "-C", rownames(amp_exp_df))

# calculate the relative abundance
amp_exp_df_rl <- rel_abundance(amp_exp_df)
# Replace '-MAB' with '-C' in the row names
rownames(amp_exp_df_rl) <- gsub("-MAB", "-C", rownames(amp_exp_df_rl))
#---------------------------------------------------------------------------------------------------
merge_otu_table_t <- merge(otu_table_t, metadata_all, by.x = "SampleID")
merge_otu_table_t$experiment <- as.factor(merge_otu_table_t$experiment)
merge_otu_table_t$source <- as.factor(merge_otu_table_t$source)

# Get the factor column name
factor_column <- "source"

# Iterate through each factor variable
factor_levels <- levels(as.factor(merge_otu_table_t[[factor_column]]))
for (level in factor_levels) {
  # Filter the data frame for the current factor level
  Si8_FP001 <- merge_otu_table_t[merge_otu_table_t[[factor_column]] == level, ]
  numeric_columns <- Si8_FP001[, sapply(Si8_FP001, is.numeric)]
  numeric_columns <- numeric_columns[, colSums(numeric_columns) >= 10]

  numeric_columns <- numeric_columns[, colSums(numeric_columns) > 1]

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
  otu_table_final <- Si8_FP001_Final

  row.names(otu_table_final) <- otu_table_final$SampleID
  otu_table_final <- otu_table_final[, -1]
  metadata <- Si8_FP001_Metadata
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
  #clr-transformation
  se <- SummarizedExperiment(assays = list(counts = counts),
                             colData = samples,
                             rowData = tax)

  # Convert to TreeSummarizedExperiment
  tse <- as(se, "TreeSummarizedExperiment")

  # Transform assay to clr
  tse <- transformAssay(tse, method = "clr", pseudocount = 1)

  # Extract clr assay values
  clr_assay <- assays(tse)$clr
  #clr_assay <- assays(tse)[[1]]
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

  # Rename the first column to "Row.names"
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
  amp_exp_df_raw <- as.matrix(amp_exp_df[common_rows_final, , drop = FALSE])

  file_name13 <- paste(level,  "_relAMP_scatterplot.pdf", sep = "_")
  file_name14 <- paste(level,  "_rawAMP_scatterplot.pdf", sep = "_")
  generate_scatter_plots(subset_amp_exp_df_rl, subset_clr_values, file_name13)
  generate_scatter_plots(amp_exp_df_raw, subset_clr_values, file_name14)
  #ggsave(file_name13, p, width = 20, height = 10, limitsize = FALSE)
  #---------------- evaluate correlation using relative expression of AMP
  out <- corr.test(subset_amp_exp_df_rl, subset_clr_values, use = "pairwise", method="spearman",adjust="BH", alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
  R_value_all <- out$r
  P_value_all <- out$p
  P_adj_all <- out$p.adj
  #---------------- evaluate correlation using raw expression of AMP
  out_raw <- corr.test(amp_exp_df_raw, subset_clr_values, use = "pairwise", method="spearman",adjust="BH", alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
  raw_R_value_all <- out_raw$r
  raw_P_value_all <- out_raw$p
  raw_P_adj_all <- out_raw$p.adj
  #---------------- Write Correlation matrices (using rel. expression) in text files
  file_name5 <- paste(level,  "_Genera_REL_AMP_Corr.txt", sep = "_")
  write.table(R_value_all, file= file_name5, sep='\t', quote=F)

  file_name6 <- paste(level,  "_Genera_REL_AMP_Pval.txt", sep = "_")
  write.table(P_value_all, file= file_name6, sep='\t', quote=F)

  file_name7 <- paste(level,  "_Genera_REL_AMP_AdjPval.txt", sep = "_")
  write.table(P_adj_all, file= file_name7, sep='\t', quote=F)

  #---------------- Write Correlation matrices (using raw expression) in text files
  file_name8 <- paste(level,  "_Genera_RAW_AMP_Corr.txt", sep = "_")
  write.table(raw_R_value_all, file= file_name8, sep='\t', quote=F)

  file_name9 <- paste(level,  "_Genera_RAW_AMP_Pval.txt", sep = "_")
  write.table(raw_P_value_all, file= file_name9, sep='\t', quote=F)

  file_name10 <- paste(level,  "_Genera_RAW_AMP_AdjPval.txt", sep = "_")
  write.table(raw_P_adj_all, file= file_name10, sep='\t', quote=F)

  #----------------- Iterate through rows and columns to find values less than 0.01
  file_conn <- file(paste(level,  "_Genera_REL_AMP_signif_corr.txt", sep = "_"), "w")
  write_significant_results(P_value_all, R_value_all, P_adj_all, file_conn)
  close(file_conn)
  #----------------- Iterate through rows and columns to find values less than 0.01
  file_conn1 <- file(paste(level,  "_Genera_RAW_AMP_signif_corr.txt", sep = "_"), "w")
  write_significant_results(raw_P_value_all, raw_R_value_all, raw_P_adj_all, file_conn)
  close(file_conn1)
  #----------------- Heat-map using correlation values
  result_matrix <- R_value_all
  file_name11 <- paste(level,  "_Taxa_REL_AMP_expression_heatmap_spearman.pdf", sep = "_")
  generate_heatmap(result_matrix, file_name11)

  result_matrix1 <- raw_R_value_all
  file_name12 <- paste(level,  "_Taxa_RAW_AMP_expression_heatmap_spearman.pdf", sep = "_")
  generate_heatmap(result_matrix1, file_name12)

    # Call the function with your matrices
    subset_result <- subset_matrices(raw_P_value_all, raw_R_value_all)
    subset_result1 <- subset_matrices(P_value_all, R_value_all)
    # Extract the subset matrices
    subset_P_value <- subset_result$P_value
    subset_R_value <- subset_result$R_value
    
    subset_P_value1 <- subset_result1$P_value
    subset_R_value1 <- subset_result1$R_value
    
    # Check if result_matrix2 is not empty
    if (nrow(subset_R_value) != 0 && ncol(subset_R_value) != 0) {
      # If result_matrix2 is not empty, proceed with the code
      result_matrix2 <- subset_R_value
      file_name12 <- paste(level, "_SIGNIFICANT_RAW_HEATMAP.pdf", sep = "_")
      generate_heatmap(result_matrix2, file_name12)
    }
    else{
      print("EMPTY MATRIX!!!")
    }
    # Check if result_matrix2 is not empty
    if (nrow(subset_R_value1) != 0 && ncol(subset_R_value1) != 0) {
      # If result_matrix2 is not empty, proceed with the code
      result_matrix3 <- subset_R_value1
      file_name13 <- paste(level, "_SIGNIFICANT_REL_HEATMAP.pdf", sep = "_")
      generate_heatmap(result_matrix3, file_name13)
    }
    else{
      print("EMPTY MATRIX!!!")
    }
}
rm(list = ls())
