
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
taxonomy <- read.csv("selected_taxonomy.tsv", sep = "\t", row.names = 1)
taxonomy <- as.matrix(taxonomy)
# Read in metadata
metadata <- read.table("Metadata_corrected_analysis.txt", sep="\t", row.names = 1, header=T)
# Read in tree
phy_tree <- read_tree("tree.nwk")
# Import all as phyloseq objects

counts  <- otu_table_in  # Abundance table (e.g. ASV data; to assay data)
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
  }else if (element == "diet_sex"){
    source_Colors <- c("darkgreen", "darkolivegreen4", "red", "salmon4")
    source_Colors1 <- c("darkolivegreen4", "darkgreen", "salmon4", "red")
  }
  #---------------------------PCoA Analysis
  euclidean_dist <- vegan::vegdist(clr_assay, method = "euclidean")
  #UniFrac_distances <- UniFrac(Bushman2, weighted=TRUE)
  UniFrac_distances <- as.matrix(euclidean_dist)
  UniFrac_dist_column <- melt(UniFrac_distances)
  write.table (UniFrac_dist_column, file = "UniFrac_distances.txt", sep = "\t")
  
  Uni_pcoa <- pcoa(UniFrac_distances)
  Uni_pcoa$values[1:2,]
  mds.var.per = round(Uni_pcoa$values$Eigenvalues/sum(Uni_pcoa$values$Eigenvalues)*100, 1)
  pc <- c(1,2)
  file_name1 <- paste(element, "Eucledian.jpg", sep = "_")
  jpeg(file_name1, height = 10, width = 10, units = 'in', res = 600)
  plot(Uni_pcoa$vectors[,1:2], bg= Colors1[as.factor(samples[[element]])], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
  #text(Uni_pcoa$vectors[,1:2], labels=rownames(Bushman2@sam_data), cex=0.3, font=1, pos=1)
  ordiellipse(Uni_pcoa$vectors[,1:2], as.factor(samples[[element]]), kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = Colors1)
  #ordispider(Uni_pcoa$vectors[,1:2], as.data.frame(samples)$Status, lty=3, spider ="centroid", lwd=1, col="black")
  legend("bottomleft", legend = sort(unique(as.factor(samples[[element]]))), col = Colors1, lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
  abline(h=0, v=0, col = "gray60")
  dev.off ()
  adonis2(UniFrac_distances ~ samples$diet_sex)

}
