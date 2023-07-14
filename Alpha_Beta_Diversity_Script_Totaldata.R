#This script is used to carry out analysis using the ASV table obtained after QIIME2 analysis

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

# Finally merge to create Phyloseq object!
ps <- phyloseq(OTU, TAX, META, phy_tree)

covar <- as.vector(colnames(metadata))
# Elements to remove
elements_to_remove <- c("mouse.number", "sample")

# Remove specific elements
covar_upd <- covar[!(covar %in% elements_to_remove)]
for (i in 1:length(covar_upd)) {
  element <- covar_upd[i]
  if (element == "diet") {
    Colors <- c("darkolivegreen4", "red", "orange")
    Colors1 <- c("darkolivegreen","red4", "orangered")
  } else if (element == "experiment"){
    Colors <- c("mediumpurple3", "snow4", "aquamarine")
    Colors1 <- c("violet", "gray30", "deepskyblue")
  }else if (element == "Mouse_provider"){
    Colors <- c("goldenrod1", "salmon4")
    Colors1 <- c("goldenrod4", "salmon")
  }else if (element == "sex"){
    Colors <- c("lightblue", "gray")
    Colors1 <- c("midnightblue", "darkgray")
  }else if (element == "sex"){
    source_Colors <- c("darkseagreen", "khaki")
    source_Colors1 <- c("cadelblue4", "khaki4")
  }
  #f_col1 <- paste(element, "Colors1", sep = "_")
  #f_col <- paste(element, "Colors", sep = "_")
  p <- plot_richness(ps, element, measures = c("Observed", "Chao1", "Shannon"), color = NULL, shape = NULL)
  p <- p + geom_boxplot(data = p$data, aes_string(x= element, y = p$data$value, fill = element)) + 
    labs(x="",y="Alpha Diversity Measure") + 
    theme_bw() +
    scale_color_manual(values = Colors1)+
    scale_fill_manual(values = Colors)+
    geom_point(aes_string(colour = element), position=position_jitterdodge(jitter.width = 0.5), size = 3) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 24))
  file_name1 <- paste(element, "phyloseq_analysis-richness_estimates.pdf", sep = "_")
  file_name2 <- paste(element, "phyloseq_analysis-richness_estimates.jpg", sep = "_")
  ggsave(file_name1, p, width = 10, height = 10, limitsize = FALSE)
  jpeg(file_name2, height = 10, width = 15, units = 'in', res = 600)
  p
  dev.off ()
  #---------------------------PCoA Analysis
  
  Bushman2  = transform_sample_counts(ps, function(x) x / sum(x) )

  #---------------------- weighted UniFrac PCoA
  file_name3 <- paste(element, "W-UniFrac_distances.txt", sep = "_")
  UniFrac_distances <- UniFrac(Bushman2, weighted=TRUE)
  UniFrac_distances <- as.matrix(UniFrac_distances)
  UniFrac_dist_column <- melt(UniFrac_distances)
  write.table (UniFrac_dist_column, file = file_name3, sep = "\t")
  
  Uni_pcoa <- pcoa(UniFrac_distances)
  Uni_pcoa$values[1:2,]
  mds.var.per = round(Uni_pcoa$values$Eigenvalues/sum(Uni_pcoa$values$Eigenvalues)*100, 1)
  pc <- c(1,2)
  file_name4 <- paste(element, "W-Unifrac_PCOA.jpg", sep = "_")
  jpeg(file_name4, height = 10, width = 10, units = 'in', res = 600)
  plot(Uni_pcoa$vectors[,1:2], bg= Colors1[as.factor(Bushman2@sam_data[[element]])], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
  text(Uni_pcoa$vectors[,1:2], labels=rownames(Bushman2@sam_data), cex=0.3, font=1, pos=1)
  ordiellipse(Uni_pcoa$vectors[,1:2], Bushman2@sam_data[[element]], kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = Colors1)
  ordispider(Uni_pcoa$vectors[,1:2], Bushman2@sam_data[[element]], lty=3, spider ="centroid", lwd=1, col="black")
  legend("bottomleft", legend = sort(unique(as.factor(Bushman2@sam_data[[element]]))), col = Colors1,lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
  abline(h=0, v=0, col = "gray60")
  dev.off ()
  #adonis2(UniFrac_distances ~ Bushman2@sam_data[[element]])
  
  #--------------------- unweighted UniFrac PCoA
  
  UWUniFrac_distances <- UniFrac(Bushman2, weighted=FALSE)
  UWUniFrac_distances <- as.matrix(UWUniFrac_distances)
  UniFrac_dist_column <- melt(UWUniFrac_distances)
  file_name5 <- paste(element, "unweighted_UWUniFrac_distances.txt", sep = "_")
  write.table (UniFrac_dist_column, file = file_name5, sep = "\t")
  Uni_pcoa <- pcoa(UWUniFrac_distances)
  Uni_pcoa$values[1:2,]
  mds.var.per = round(Uni_pcoa$values$Eigenvalues/sum(Uni_pcoa$values$Eigenvalues)*100, 1)
  pc <- c(1,2)
  
  file_name6 <- paste(element, "unweighted_UWUniFrac_distances.jpg", sep = "_")
  jpeg(file_name6, height = 10, width = 10, units = 'in', res = 600)
  plot(Uni_pcoa$vectors[,1:2], bg= Colors1[as.factor(Bushman2@sam_data[[element]])], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
  text(Uni_pcoa$vectors[,1:2], labels=rownames(Bushman2@sam_data), cex=0.3, font=1, pos=1)
  ordiellipse(Uni_pcoa$vectors[,1:2], Bushman2@sam_data[[element]], kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = Colors1)
  ordispider(Uni_pcoa$vectors[,1:2], Bushman2@sam_data[[element]], lty=3, spider ="centroid", lwd=1, col="black")
  legend("bottomleft", legend = sort(unique(as.factor(Bushman2@sam_data[[element]]))), col = Colors1,lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
  abline(h=0, v=0, col = "gray60")
  dev.off ()
  #adonis2(UWUniFrac_distances ~ Bushman2@sam_data$diet)
  
  #-------------------------- Bray-Curtis PCoA
  BC_distances <- distance(Bushman2, method="bray")
  BC_distances <- as.matrix(BC_distances)
  UniFrac_dist_column <- melt(BC_distances)
  
  file_name7 <- paste(element, "BC_distances.txt", sep = "_")
  write.table (UniFrac_dist_column, file = file_name7, sep = "\t")
  Uni_pcoa <- pcoa(BC_distances)
  Uni_pcoa$values[1:2,]
  mds.var.per = round(Uni_pcoa$values$Eigenvalues/sum(Uni_pcoa$values$Eigenvalues)*100, 1)
  pc <- c(1,2)
  
  file_name8 <- paste(element, "BC_distances.jpg", sep = "_")
  jpeg(file_name8, height = 10, width = 10, units = 'in', res = 600)
  plot(Uni_pcoa$vectors[,1:2], bg= Colors1[as.factor(Bushman2@sam_data[[element]])], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
  text(Uni_pcoa$vectors[,1:2], labels=rownames(Bushman2@sam_data), cex=0.3, font=1, pos=1)
  ordiellipse(Uni_pcoa$vectors[,1:2], Bushman2@sam_data[[element]], kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = Colors1)
  ordispider(Uni_pcoa$vectors[,1:2], Bushman2@sam_data[[element]], lty=3, spider ="centroid", lwd=1, col="black")
  legend("topleft", legend = sort(unique(as.factor(Bushman2@sam_data[[element]]))), col = Colors1,lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
  abline(h=0, v=0, col = "gray60")
  dev.off ()
  #adonis2(BC_distances ~ Bushman2@sam_data$diet)
}
erich <- estimate_richness(ps, measures = c("Observed", "Chao1", "Shannon"))
erich_all <- cbind(erich, metadata)
write.table (erich_all, file = "Richness_total.txt", sep = "\t")



adonis2(UniFrac_distances ~ Bushman2@sam_data[[element]])
