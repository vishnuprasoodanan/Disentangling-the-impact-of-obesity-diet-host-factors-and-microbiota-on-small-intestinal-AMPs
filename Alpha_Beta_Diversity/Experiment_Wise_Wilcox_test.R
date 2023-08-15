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
    Wilcox_df <- data.frame(
      Observed = numeric(2),
      Chao1 = numeric(2),
      Shannon = numeric(2),
      Group = character(2),
      row.names = c("Diet", "Sex")
    )
    erich <- estimate_richness(ps, measures = c("Observed", "Chao1", "Shannon"))
    erich_Status <- cbind(erich, metadata)
    write.table (erich_Status, file = "Richness_total.txt", sep = "\t")
    column_to_check <- "diet"
    unique_values_count <- erich_Status %>%
    pull({{ column_to_check }}) %>%
    n_distinct()
      
    if (unique_values_count == 2){
      diet_Obs <- wilcox.test(erich_Status$Observed ~ erich_Status$diet, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)$p.value
      diet_Chao1 <- wilcox.test(erich_Status$Chao1 ~ erich_Status$diet, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)$p.value
      diet_Shannon <- wilcox.test(erich_Status$Shannon ~ erich_Status$diet, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)$p.value
      Wilcox_df[1, ] <- c(diet_Obs, diet_Chao1, diet_Shannon, paste(level, level2, sep = "_"))
      } else {
      Wilcox_df[1, ] <- c(0, 0, 0, paste(level, level2, sep = "_"))
      }
      # Check if a column has more than 2 unique values
    column_to_check <- "sex"
    unique_values_count <- erich_Status %>%
    pull({{ column_to_check }}) %>%
    n_distinct()
      
    if (unique_values_count == 2){
      sex_Obs <- wilcox.test(erich_Status$Observed ~ erich_Status$sex, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)$p.value
      sex_Chao1 <- wilcox.test(erich_Status$Chao1 ~ erich_Status$sex, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)$p.value
      sex_Shannon <- wilcox.test(erich_Status$Shannon ~ erich_Status$sex, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)$p.value
      Wilcox_df[2, ] <- c(sex_Obs, sex_Chao1, sex_Shannon, paste(level, level2, sep = "_"))
      }else{
      Wilcox_df[2, ] <- c(0, 0, 0, paste(level, level2, sep = "_"))
      }
      #file_name <- paste(level, level2, "_Diet_Sex_Wilcox_Test.txt", sep = "_")
      #write.table (Wilcox_df, file = file_name, sep = "\t", row.names = TRUE, col.names = NA)
      Wilcox_df_final  <- bind_rows(Wilcox_df, Wilcox_df_final)
    }
}

write.table (Wilcox_df_final, file = "Experiment_Wise_Wilcox_test.txt", sep = "\t", row.names = TRUE, col.names = NA)

