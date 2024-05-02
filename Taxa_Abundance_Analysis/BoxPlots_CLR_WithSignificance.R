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
library(pgirmess)
library(gridExtra)
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
    
    #colnames(otu_table_final) <- gsub("-MAB", "-C", colnames(otu_table_final))
    #rownames(metadata) <- gsub("-MAB", "-C", rownames(metadata))
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
    genus_rel <- clr_values
    # Create a new column with row names
    genus_rel$SampleID <- rownames(genus_rel)
    
    # Move the 'SampleID' column to the first position
    genus_rel <- genus_rel[, c("SampleID", setdiff(names(genus_rel), "SampleID"))]
    rownames(gut_Metadata) <- gut_Metadata$SampleID
    
    #genus_rel <- genus_rel[, c(ncol(genus_rel), 1:(ncol(genus_rel)-1))]
    # Compare row names
    if(identical(rownames(genus_rel), rownames(gut_Metadata))) {
      print("Matching")
    } else {
      print("Not matching")
    }
    boxplot_input <- merge(genus_rel, gut_Metadata[, "diet", drop = FALSE], by = "row.names")
    boxplot_input$Row.names <- NULL
    boxplot_input$diet <- factor(boxplot_input$diet)
    rownames(boxplot_input) <- boxplot_input$SampleID
    
    # Move the diet column to the second position
    boxplot_input <- boxplot_input %>%
      select(SampleID, diet, everything())
    boxplot_input[, -c(1, 2)] <- apply(boxplot_input[, -c(1, 2)], 2, as.numeric)
    
    # Print the modified dataframe
    #print(boxplot_input)
    # Extract levels of the factor column diet
    lev_group <- levels(boxplot_input$diet)
    #print(kruskalmc(boxplot_input$d__Bacteria_p__Firmicutes_c__Clostridia_o__Lachnospirales_f__Lachnospiraceae_g__Roseburia ~ boxplot_input$diet, data = boxplot_input, probs = 0.05, alpha = 0.05))
    ####boxplot script
    #print(lev_group)
    file_name2 <- paste(level, level2, "diet","_boxplot.pdf", sep = "_")
    pdf(file = file_name2)
    colnames(boxplot_input[,1:ncol(boxplot_input)]) -> Species_name
    
   # Colors <- c("darkgreen", "darkolivegreen4", "red", "red4","#1F77B4", "#FF7F0E","#2CA02C", "#D62728", "deepskyblue", "goldenrod1", "darkolivegreen", "salmon", "mediumpurple3", "paleturquoise4")
   # Colors1 <- c("darkolivegreen4", "darkgreen", "red4","red", "midnightblue", "orange", "green", "salmon4", "purple", "blue")
    
    Colors <- c("darkgreen", "red", "red4","#1F77B4", "#FF7F0E","#2CA02C", "#D62728", "deepskyblue", "goldenrod1", "darkolivegreen", "salmon", "mediumpurple3", "paleturquoise4")
    Colors1 <- c("darkolivegreen4", "red4","red", "midnightblue", "orange", "green", "salmon4", "purple", "blue")
    
    plot_lst <- vector("list", length = ncol(boxplot_input))
    for (i in 3:ncol(boxplot_input)) {    
      
      species = Species_name[i]
      
      # Search variable in rownames of df1
      variable <- species
      index <- match(variable, rownames(taxonomy))
      
      # Print corresponding value of first column of df2
      if (!is.na(index)) {
        corresponding_value <- taxonomy[index, 1]
        print(corresponding_value)
      } else {
        print("Variable not found in rownames of df1.")
      }
      species <- corresponding_value
      data1 = boxplot_input[,c(1:2,i)]
      colnames(data1) <- c("ID", "Status", "name")
      
      print(species)
      species <- gsub("_", " ", species)
      
      #CHECK SPECIFIC TAXA
      # if (species == "d Bacteria p Firmicutes c Clostridia o Peptostreptococcales Tissierellales f Peptostreptococcaceae g Romboutsia" && level == "Si8_FP001" && level2 == "Content")
      # {
      #   print(kruskalmc(name ~ Status, data = data1, probs = 0.05, alpha = 0.05))
      # }
      
      # Check the number of unique values in the 'Status' column
      unique_values <- unique(data1$Status)
      num_unique <- length(unique_values)
      
      # Initialize a vector to store significance stars
      significance <- rep("", length(data1$Status))
      significant_pairs <- NULL
      # Perform the appropriate test based on the number of unique values
      if (num_unique == 2) {
        # Wilcoxon test
        p_val <- wilcox.test(name ~ Status, data = data1, paired = FALSE, p.adjust.method = "fdr", alternative = "two.sided", conf.level = 0.95)$p.value
        # Check significance
        if (p_val < 0.05) {
          significant_levels <- levels(data1$Status)
          # Combine levels into pairs
          significant_pairs <- combn(significant_levels, 2, simplify = FALSE)
        }
      } else if (num_unique > 2) {
        # Kruskal-Wallis test
        kruskal_test_df <- as.data.frame(kruskalmc(name ~ Status, data = data1, probs = 0.05, alpha = 0.05))
        kruskal_test_df$Comparisons <- rownames(kruskal_test_df)
        
        
        # Move the new column to the first position
        kruskal_test_df <- kruskal_test_df[, c(ncol(kruskal_test_df), 1:(ncol(kruskal_test_df)-1))]
        significant_pairs <- kruskal_test_df %>%
          filter(dif.com.stat.signif) %>%
          select(Comparisons)
        
        # Extract significant pairs
        significant_pairs <- strsplit(significant_pairs$Comparisons, "-")
        significant_pairs <- lapply(significant_pairs, trimws)
        #print(significant_pairs)
      }
      p <- ggplot(data1, aes(Status, name, fill = Status)) +
        ggtitle(species) +
        labs(y = "CLR-transformed") +
        geom_boxplot(outlier.shape = NA) + 
        #geom_jitter(width = 0.2, size = 4, alpha = 0.6) +
        #geom_point(aes(colour = factor(Status)), size = 4, alpha = 0.6) +
        geom_point(aes(colour = factor(Status)), position = position_jitterdodge(jitter.width = 0.5),  size = 4, alpha = 0.6) +
        scale_color_manual(values = Colors1) +
        scale_fill_manual(values = Colors) +
        theme(plot.title = element_text(size = 8, face = "bold")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
      # Add stars between significant pairs
      for (pair in significant_pairs) {
        #print(pair)
        p <- p + geom_signif(comparisons = list(pair), map_signif_level = TRUE, textsize = 6)
      }
      plot_lst[[i]] <- p
    }
    ml <- marrangeGrob(plot_lst, nrow = 1, ncol = 1)
    print(ml)
    dev.off()
  }
}
