
# ---------------------------
# Load required packages
# ---------------------------
pkgs <- c("vegan", "ggplot2", "ggrepel", "compositions", "ade4", "gridExtra")
to_install <- pkgs[!pkgs %in% installed.packages()]
if(length(to_install)) install.packages(to_install, repos="http://cran.us.r-project.org")
lapply(pkgs, library, character.only = TRUE)

# ---------------------------
# Load data
# ---------------------------
mb_data <- read.delim("mb_data_content.txt", header=TRUE, row.names=1, sep="\t", check.names=FALSE)
amp_data <- read.delim("amp_data_content.txt", header=TRUE, row.names=1, sep="\t", check.names=FALSE)
metadata <- read.delim("metadata_all.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# ---------------------------
# Preprocessing
# ---------------------------
mb_data <- as.matrix(mb_data)
storage.mode(mb_data) <- "numeric"
amp_data <- as.matrix(amp_data)
storage.mode(amp_data) <- "numeric"

common_samples <- Reduce(intersect, list(colnames(mb_data), colnames(amp_data), metadata$sample_id))
mb_data <- mb_data[, common_samples]
amp_data <- amp_data[, common_samples]
metadata <- metadata[match(common_samples, metadata$sample_id), ]

# ---------------------------
# AMP data type check
# ---------------------------
amp_range <- range(amp_data, na.rm=TRUE)
if(amp_range[2] > 1000 & all(amp_data >= 0)) {
  amp_type <- "raw_count"
} else if(amp_range[2] < 30 & all(amp_data >= 0)) {
  amp_type <- "log_transformed"
} else {
  stop("AMP data does not appear to be clearly raw or log-transformed.")
}
cat("AMP data type:", amp_type, "\n")

# ---------------------------
# AMP data transformations
# ---------------------------
transform_list <- list()
if(amp_type == "raw_count") {
  transform_list[["raw"]] <- amp_data
  transform_list[["log"]] <- log2(amp_data + 1)
} else if(amp_type == "log_transformed") {
  transform_list[["log"]] <- amp_data
  transform_list[["exp"]] <- exp(amp_data)
}

# ---------------------------
# Procrustes function with visualizations
# ---------------------------
perform_procrustes <- function(mb_sub, amp_sub, group_name, out_dir) {
  clr_mb <- t(mb_sub)
  #clr_mb <- clr(t(mb_sub + 1e-6))
  dist_mb <- dist(clr_mb)
  pcoa_mb <- cmdscale(dist_mb, k = 2, eig = TRUE)
  colnames(pcoa_mb$points) <- c("PC1_mb", "PC2_mb")

  amp_scaled <- scale(t(amp_sub))
  dist_amp <- dist(amp_scaled)
  pca_amp <- prcomp(amp_scaled)
  Y <- pca_amp$x[, 1:2]
  colnames(Y) <- c("PC1_amp", "PC2_amp")

  res <- protest(pcoa_mb$points, Y, permutations = 999, symmetric = TRUE)

  df_X <- as.data.frame(pcoa_mb$points)
  df_Y <- as.data.frame(Y)
  df_X$sample <- rownames(df_X)
  df_Y$sample <- rownames(df_Y)
  merged <- merge(df_X, df_Y, by = "sample")

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Procrustes plot
  p1 <- ggplot() +
    geom_segment(data = merged, aes(x = PC1_mb, y = PC2_mb, xend = PC1_amp, yend = PC2_amp),
                 arrow = arrow(length = unit(0.1, "cm")), alpha = 0.5) +
    geom_point(data = merged, aes(x = PC1_mb, y = PC2_mb), color = "blue", size = 2) +
    geom_point(data = merged, aes(x = PC1_amp, y = PC2_amp), color = "red", size = 2) +
    geom_text_repel(data = merged, aes(x = PC1_amp, y = PC2_amp, label = sample), size = 2.5, max.overlaps = Inf) +
    labs(title = paste0("Procrustes MDS (CLR) - ", group_name),
         subtitle = paste("m12² =", round(res$ss, 3), "p =", res$signif)) +
    theme_minimal()
  ggsave(file.path(out_dir, paste0("Procrustes_CLR_", group_name, ".pdf")), p1, width = 7, height = 5)

  # Residual plot
  resids <- sqrt(residuals(res))
  pdf(file.path(out_dir, paste0("Residuals_", group_name, ".pdf")))
  barplot(resids, las = 2, col = "steelblue", main = paste("Procrustes Residuals -", group_name))
  dev.off()

  # Scree plots
  pdf(file.path(out_dir, paste0("Scree_", group_name, ".pdf")))
  par(mfrow = c(1, 2))
  barplot(pcoa_mb$eig / sum(pcoa_mb$eig), main = "Scree - Microbiome (PCoA)")
  barplot(pca_amp$sdev^2 / sum(pca_amp$sdev^2), main = "Scree - AMP (PCA)")
  dev.off()

  # Mantel test scatterplot
  mantel_res <- mantel(dist_mb, dist_amp)
  pdf(file.path(out_dir, paste0("Mantel_", group_name, ".pdf")))
  plot(as.vector(dist_mb), as.vector(dist_amp),
       xlab = "Microbiome Distance", ylab = "AMP Distance",
       main = paste("Mantel Test (r =", round(mantel_res$statistic, 2),
                    ", p =", mantel_res$signif, ")"))
  abline(lm(as.vector(dist_amp) ~ as.vector(dist_mb)), col = "red")
  dev.off()

  # Co-Inertia Analysis
  pca1 <- dudi.pca(clr_mb, scannf = FALSE, nf = 2)
  pca2 <- dudi.pca(amp_scaled, scannf = FALSE, nf = 2)
  coi <- coinertia(pca1, pca2, scannf = FALSE)
  pdf(file.path(out_dir, paste0("CoInertia_", group_name, ".pdf")))
  plot(coi, main = paste("Co-Inertia Analysis -", group_name))
  dev.off()

  return(res)
}

# ---------------------------
# Run analyses for each grouping column
# ---------------------------
for(grouping_col in colnames(metadata)[colnames(metadata) != "sample_id"]) {
  group_vector <- metadata[[grouping_col]]
  names(group_vector) <- metadata$sample_id

  for(trans_name in names(transform_list)) {
    amp_trans <- transform_list[[trans_name]]
    out_root <- file.path("Results", grouping_col, paste0(trans_name, "_amp"))
    dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
    sink(file = file.path(out_root, "procrustes_summary_CLR.txt"))

    cat(">>> Procrustes Analysis (CLR) using", trans_name, "AMP data - Grouping:", grouping_col, "<<<\n\n")

    res_all <- perform_procrustes(mb_data, amp_trans, "All", out_root)
    cat("All Samples:\n")
    cat("  m12² =", res_all$ss, "p =", res_all$signif, "\n\n")

    for(status in unique(group_vector)) {
      samp_grp <- names(group_vector)[group_vector == status]
      if(length(samp_grp) < 3) {
        cat("Group", status, ": Too few samples\n")
        next
      }
      mb_grp <- mb_data[, samp_grp]
      amp_grp <- amp_trans[, samp_grp]
      res_grp <- perform_procrustes(mb_grp, amp_grp, status, out_root)
      cat("Group", status, ":\n")
      cat("  m12² =", res_grp$ss, "p =", res_grp$signif, "\n\n")
    }

    sink()
  }
}
