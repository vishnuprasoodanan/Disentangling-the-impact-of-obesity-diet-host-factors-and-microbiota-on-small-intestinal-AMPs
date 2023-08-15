library(tibble)
library(tidyverse)
library(plyr)
library(dplyr)
library(reshape2)
library(phyloseq)
library (ape)
library (ggplot2)
library(vegan)

otu_table_in <- read.csv("selected_feature-table.txt", sep = "\t")
otu_table_t <- setNames(data.frame(t(otu_table_in[,-1])), otu_table_in[,1])
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
adonis_results <- data.frame()
# Get the factor column name
factor_column <- "experiment"
Wilcox_df_final <- data.frame()
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
    
    otu_table_final <- gut_Final
    row.names(otu_table_final) <- otu_table_final$SampleID
    otu_table_final <- otu_table_final[, -1]
    metadata <- gut_Metadata
    row.names(metadata) <- metadata$SampleID
    metadata <- metadata[, -1]
    taxonomy <- read.csv("taxonomy.tsv", sep = "\t", row.names = 1)
    taxonomy <- as.matrix(taxonomy)
    # Read in metadata
    # Read in tree
    phy_tree <- read_tree("tree.nwk")
    # Import all as phyloseq objects
    OTU <- otu_table(otu_table_final, taxa_are_rows = TRUE)
    TAX <- tax_table(taxonomy)
    META <- sample_data(metadata)
    # Sanity checks for consistent OTU names
    head(taxa_names(TAX))
    head(taxa_names(OTU))
    head(taxa_names(phy_tree))
    # Same sample names
    sample_names(OTU) <- gsub('[.]', '-', sample_names(OTU))
    sample_names(OTU) <- gsub('^[X]', '', sample_names(OTU)) 
    head(sample_names(OTU))
    head(sample_names(META))
    
    # Finally merge to create Phyloseq object!
    ps <- phyloseq(OTU, TAX, META, phy_tree)
    
    covar <- as.vector(colnames(metadata))
    # Elements to remove
    elements_to_remove <- c("mouse.number", "sample", "SampleID", "Mouse_provider", "experiment", "source")
    
    # Remove specific elements
    covar_upd <- covar[!(covar %in% elements_to_remove)]
    for (i in 1:length(covar_upd)) {
      element <- covar_upd[i]
      
      #---------------------------PCoA Analysis
      
      Bushman2  = transform_sample_counts(ps, function(x) x / sum(x) )
      
      #---------------------- weighted UniFrac PCoA
      UniFrac_distances <- UniFrac(Bushman2, weighted=TRUE)
      UniFrac_distances <- as.matrix(UniFrac_distances)
      
      #column_to_check <- "diet"
      unique_values_count <- Bushman2@sam_data %>%
        pull({{ element }}) %>%
        n_distinct()
      
      if (unique_values_count > 1){
        W_UniFrac_adonis_results <- data.frame()
        out <- as.data.frame(adonis2(UniFrac_distances ~ Bushman2@sam_data[[element]]))
        out1 <- bind_cols(out, status = element)
        out1 <- bind_cols(out1, Dist = "Weighted-Unifrac", Group = paste(level, level2, sep = "_"))
        W_UniFrac_adonis_results  <- bind_rows(out1, W_UniFrac_adonis_results )
        adonis_results <- bind_rows(adonis_results, W_UniFrac_adonis_results)
        # Print the combined data frame
        print(W_UniFrac_adonis_results )
      }
      
      #---------------------- Unweighted UniFrac PCoA
      UWUniFrac_distances <- UniFrac(Bushman2, weighted=FALSE)
      UWUniFrac_distances <- as.matrix(UWUniFrac_distances)
      
      #column_to_check <- "diet"
      unique_values_count <- Bushman2@sam_data %>%
        pull({{ element }}) %>%
        n_distinct()
      
      if (unique_values_count > 1){
        UW_UniFrac_adonis_results <- data.frame()
        out <- as.data.frame(adonis2(UWUniFrac_distances ~ Bushman2@sam_data[[element]]))
        out1 <- bind_cols(out, status = element)
        out1 <- bind_cols(out1, Dist = "UnWeighted-Unifrac", Group = paste(level, level2, sep = "_"))
        UW_UniFrac_adonis_results  <- bind_rows(out1, UW_UniFrac_adonis_results )
        adonis_results <- bind_rows(adonis_results, UW_UniFrac_adonis_results)
        # Print the combined data frame
        print(UW_UniFrac_adonis_results )
      }
      
      #---------------------- Bray-Curtis PCoA
      BC_distances <- distance(Bushman2, method="bray")
      BC_distances <- as.matrix(BC_distances)
      
      #column_to_check <- "diet"
      unique_values_count <- Bushman2@sam_data %>%
        pull({{ element }}) %>%
        n_distinct()
      
      if (unique_values_count > 1){
        BC_adonis_results <- data.frame()
        out <- as.data.frame(adonis2(BC_distances ~ Bushman2@sam_data[[element]]))
        out1 <- bind_cols(out, status = element)
        out1 <- bind_cols(out1, Dist = "Bray-Curtis", Group = paste(level, level2, sep = "_"))
        BC_adonis_results  <- bind_rows(out1, BC_adonis_results )
        adonis_results <- bind_rows(adonis_results, BC_adonis_results)
        # Print the combined data frame
        print(BC_adonis_results )
      }
    }
  }
}
write.table (adonis_results, file = "Experiment_Wise_PERMANOVA_Results.txt", sep = "\t", row.names = TRUE, col.names = NA)
