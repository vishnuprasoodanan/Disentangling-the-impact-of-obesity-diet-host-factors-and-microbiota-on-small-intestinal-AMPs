library(phyloseq)
library (ape)
library (ggplot2)

#------------------------- Alpha diversity analysis on complete data (Total-Data)
#-- Make a phyloseq object after filtering of OTUs
# Read in OTU table
otu_table_in <- read.csv("selected_feature-table.txt", sep = "\t", row.names = 1)
total_asv <- colSums(otu_table_in)
range (total_asv)
sd (total_asv)
mean (total_asv)

# Read in taxonomy
# Separated by kingdom, phylum, class, order, family, genus, species
taxonomy <- read.csv("taxonomy.tsv", sep = "\t", row.names = 1)
taxonomy <- as.matrix(taxonomy)
# Read in metadata
metadata <- read.table("selected_Metadata.tsv", sep="\t", row.names = 1, header=T)
# Read in tree
phy_tree <- read_tree("tree.nwk")
# Import all as phyloseq objects
OTU <- otu_table(otu_table_in, taxa_are_rows = TRUE)
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
adonis_results <- data.frame()
# Finally merge to create Phyloseq object!
ps <- phyloseq(OTU, TAX, META, phy_tree)

covar <- as.vector(colnames(metadata))
# Elements to remove
elements_to_remove <- c("mouse.number", "sample")

# Remove specific elements
covar_upd <- covar[!(covar %in% elements_to_remove)]
for (i in 1:length(covar_upd)) {
  element <- covar_upd[i]

  #---------------------------PCoA Analysis
  
  Bushman2  = transform_sample_counts(ps, function(x) x / sum(x) )
  
  #---------------------- weighted UniFrac PCoA
  UniFrac_distances <- UniFrac(Bushman2, weighted=TRUE)
  UniFrac_distances <- as.matrix(UniFrac_distances)
  W_UniFrac_adonis_results <- data.frame()
  out <- as.data.frame(adonis2(UniFrac_distances ~ Bushman2@sam_data[[element]]))
  out1 <- bind_cols(out, status = element)
  out1 <- bind_cols(out1, Dist = "Weighted-Unifrac")
  W_UniFrac_adonis_results  <- bind_rows(out1, W_UniFrac_adonis_results )
  adonis_results <- bind_rows(adonis_results, W_UniFrac_adonis_results)
  # Print the combined data frame
  print(W_UniFrac_adonis_results )
  #--------------------- unweighted UniFrac PCoA
  
  UWUniFrac_distances <- UniFrac(Bushman2, weighted=FALSE)
  UWUniFrac_distances <- as.matrix(UWUniFrac_distances)
  UW_Unifrac_adonis_results <- data.frame()
  out <- as.data.frame(adonis2(UWUniFrac_distances ~ Bushman2@sam_data[[element]]))
  out1 <- bind_cols(out, status = element)
  out1 <- bind_cols(out1, Dist = "UnWeighted-Unifrac")
  UW_Unifrac_adonis_results <- bind_rows(out1, UW_Unifrac_adonis_results)
  # Print the combined data frame
  print(UW_Unifrac_adonis_results)
  adonis_results <- bind_rows(adonis_results, UW_Unifrac_adonis_results)

  #-------------------------- Bray-Curtis PCoA
  BC_distances <- distance(Bushman2, method="bray")
  BC_distances <- as.matrix(BC_distances)

  BC_adonis_results <- data.frame()
  element <- covar_upd[i]
  out <- as.data.frame(adonis2(BC_distances ~ Bushman2@sam_data[[element]]))
  out1 <- bind_cols(out, status = element)
  out1 <- bind_cols(out1, Dist = "Bray-Curtis")
  BC_adonis_results <- bind_rows(out1, BC_adonis_results)
  # Print the combined data frame
  print(BC_adonis_results)
  adonis_results <- bind_rows(adonis_results, BC_adonis_results)

}
write.table (adonis_results, file = "Adonis_Test_Results.txt", sep = "\t", col.names = NA, row.names = TRUE)
