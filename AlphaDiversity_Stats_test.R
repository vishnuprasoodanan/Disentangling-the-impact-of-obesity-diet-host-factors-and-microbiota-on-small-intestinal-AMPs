library(phyloseq)
library(ape)
library(ggplot2)
library(vegan)
library(pgirmess)
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

# Finally merge to create Phyloseq object!
ps <- phyloseq(OTU, TAX, META, phy_tree)

wilcox_results <- data.frame()
KW_results <- data.frame()
KW_results_exp <- data.frame()

Wilcox_df <- data.frame(
  Observed = numeric(3),
  Chao1 = numeric(3),
  Shannon = numeric(3),
  row.names = c("MouseProvider", "Sex", "Source")
)

erich <- estimate_richness(ps, measures = c("Observed", "Chao1", "Shannon"))
erich_Status <- cbind(erich, metadata)
write.table (erich_Status, file = "Richness_total_v4.txt", sep = "\t")

covar <- as.vector(colnames(metadata))
elements_to_remove <- c("mouse.number", "sample")
covar_upd <- covar[!(covar %in% elements_to_remove)]
for (i in 1:length(covar_upd)) {
  element <- covar_upd[i]
  if (element == "diet") {
    test1 <- as.data.frame(kruskalmc(erich_Status$Observed ~ erich_Status$diet, probs=0.05))
    KW_results  <- bind_rows(test1, KW_results)
    test2 <- as.data.frame(kruskalmc(erich_Status$Chao1 ~ erich_Status$diet, probs=0.05))
    KW_results  <- bind_rows(test2, KW_results)
    test3 <- as.data.frame(kruskalmc(erich_Status$Shannon ~ erich_Status$diet, probs=0.05))
    KW_results  <- bind_rows(test3, KW_results)
    KW_results <- bind_cols(KW_results, AlphaMeasure = c(rep("Shannon", times = 3), rep("Chao1", times = 3), rep("Observed", times = 3)))
  }else if (element == "experiment"){
    testA1 <- as.data.frame(kruskalmc(erich_Status$Observed ~ erich_Status$experiment, probs=0.05))
    KW_results_exp  <- bind_rows(testA1, KW_results_exp)
    testA2 <- as.data.frame(kruskalmc(erich_Status$Chao1 ~ erich_Status$experiment, probs=0.05))
    KW_results_exp  <- bind_rows(testA2, KW_results_exp)
    testA3 <- as.data.frame(kruskalmc(erich_Status$Shannon ~ erich_Status$experiment, probs=0.05))
    KW_results_exp  <- bind_rows(testA3, KW_results_exp)
    KW_results_exp <- bind_cols(KW_results_exp, AlphaMeasure = c(rep("Shannon", times = 3), rep("Chao1", times = 3), rep("Observed", times = 3)))
  }else if (element == "Mouse_provider"){
    Mouse_provider_Obs <- wilcox.test(erich_Status$Observed ~ erich_Status$Mouse_provider, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)$p.value
    Mouse_provider_Chao1 <- wilcox.test(erich_Status$Chao1 ~ erich_Status$Mouse_provider, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)$p.value
    Mouse_provider_Shannon <- wilcox.test(erich_Status$Shannon ~ erich_Status$Mouse_provider, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)$p.value
    Wilcox_df[1, ] <- c(Mouse_provider_Obs, Mouse_provider_Chao1, Mouse_provider_Shannon)
  }else if (element == "sex"){
    sex_Obs <- wilcox.test(erich_Status$Observed ~ erich_Status$sex, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)$p.value
    sex_Chao1 <- wilcox.test(erich_Status$Chao1 ~ erich_Status$sex, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)$p.value
    sex_Shannon <- wilcox.test(erich_Status$Shannon ~ erich_Status$sex, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)$p.value
    Wilcox_df[2, ] <- c(sex_Obs, sex_Chao1, sex_Shannon)
  }else if (element == "source"){
    source_Obs <- wilcox.test(erich_Status$Observed ~ erich_Status$source, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)$p.value
    source_Chao1 <- wilcox.test(erich_Status$Chao1 ~ erich_Status$source, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)$p.value
    source_Shannon <- wilcox.test(erich_Status$Shannon ~ erich_Status$source, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)$p.value
    Wilcox_df[3, ] <- c(source_Obs, source_Chao1, source_Shannon)  
  }
}
write.table (KW_results, file = "Kruskal_Wallis_Diet.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table (KW_results_exp, file = "Kruskal_Wallis_Experiment.txt", sep = "\t", row.names = TRUE, col.names = NA)
write.table (Wilcox_df, file = "Wilcox_Test_MouseProvider.txt", sep = "\t", row.names = TRUE, col.names = NA)
