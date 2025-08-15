# monocyte-redox-scRNA

Data and outputs for the manuscript on redox-associated programs in human peripheral blood mononuclear cells (PBMCs).

- **Dataset:** 10x Genomics **“PBMCs from Citrate-Treated Cell Preparation Tubes (3′ v3.1 Chemistry)”**  
  File used: `3p_Citrate_CPT_molecule_info.h5` (gene-expression only; no ADT).  
- **Cells used:** **3,940 pass-filter cells**  
- **HVGs:** top **2,000** highly variable genes  
- **PCA:** **25 components**, cumulative explained variance **0.668**  
- **Normalization:** library-size scaling to **10,000 counts/cell** → **natural log** transform (`log1p`)  
- **PCA method:** sparse SVD via **`irlba`** on the log-normalized matrix

---

## Repository layout

## Make the figures from the CSVs

```r
# (paste your R code here exactly)
pkgs <- c("ggplot2","ggrepel","patchwork","readr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

root <- "."  # repo root

# --- reload (if needed)
scores   <- readr::read_csv(file.path(root, "data/PCA_scores_topHVGs.csv"), show_col_types = FALSE)
loadings <- readr::read_csv(file.path(root, "data/PCA_loadings_topHVGs.csv"), show_col_types = FALSE)
ev       <- readr::read_csv(file.path(root, "data/PCA_explained_variance.csv"), show_col_types = FALSE)

cell_ids <- scores[[1]]; scores <- as.data.frame(scores[,-1]); rownames(scores) <- cell_ids
gene_ids <- loadings[[1]]; loadings <- as.data.frame(loadings[,-1]); rownames(loadings) <- gene_ids

pc_names <- colnames(scores); pc1 <- pc_names[1]; pc2 <- pc_names[2]; pc3 <- if (length(pc_names) >= 3) pc_names[3] else NA
ev$PC_num <- seq_len(nrow(ev)); cum_last <- round(tail(ev$cumulative_explained_variance, 1), 3)

# --- Scree
p_scree <- ggplot(ev, aes(x = PC_num, y = explained_variance_ratio)) +
  geom_col() +
  labs(x = "Principal Component", y = "Explained variance ratio",
       title = "Scree plot (per-PC variance)") +
  theme_minimal(base_size = 12)

# --- Cumulative
p_cum <- ggplot(ev, aes(x = PC_num, y = cumulative_explained_variance)) +
  geom_line() + geom_point(size = 1.8) +
  geom_hline(yintercept = 0.85, linetype = 2) +
  annotate("text", x = max(ev$PC_num), y = 0.85, vjust = -0.5, hjust = 1,
           label = "85% target", size = 3) +
  labs(x = "Principal Component", y = "Cumulative explained variance",
       title = paste0("Cumulative variance (last = ", cum_last, ")")) +
  theme_minimal(base_size = 12)

# --- PC1 vs PC2 scatter (color by PC3 if available)
p_scatter <- if (!is.na(pc3)) {
  ggplot(scores, aes_string(x = pc1, y = pc2, color = pc3)) +
    geom_point(alpha = 0.6, size = 0.8) +
    labs(title = paste("PC scatter:", pc1, "vs", pc2),
         x = pc1, y = pc2, color = pc3) +
    theme_minimal(base_size = 12)
} else {
  ggplot(scores, aes_string(x = pc1, y = pc2)) +
    geom_point(alpha = 0.6, size = 0.8) +
    labs(title = paste("PC scatter:", pc1, "vs", pc2), x = pc1, y = pc2) +
    theme_minimal(base_size = 12)
}

# --- Biplot with top-loading genes
topN <- 15
ld <- loadings[, c(pc1, pc2), drop = FALSE]
ld$gene <- rownames(loadings)
ld$vec_len <- sqrt(ld[[pc1]]^2 + ld[[pc2]]^2)
top_ld <- ld[order(-ld$vec_len), ][1:topN, ]

rng_x <- range(scores[[pc1]], na.rm = TRUE)
rng_y <- range(scores[[pc2]], na.rm = TRUE)
scale_factor <- 0.5 * min(diff(rng_x), diff(rng_y)) / max(top_ld$vec_len)
top_ld$xend <- top_ld[[pc1]] * scale_factor
top_ld$yend <- top_ld[[pc2]] * scale_factor

p_biplot <- ggplot() +
  geom_point(data = scores, aes_string(x = pc1, y = pc2), alpha = 0.3, size = 0.6) +
  geom_segment(data = top_ld, aes(x = 0, y = 0, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.015, "npc")), linewidth = 0.3) +
  ggrepel::geom_text_repel(data = top_ld, aes(x = xend, y = yend, label = gene),
                           size = 3, max.overlaps = 100) +
  labs(title = paste0("Biplot (", pc1, " vs ", pc2, ") with top ", topN, " gene loadings"),
       x = pc1, y = pc2) +
  theme_minimal(base_size = 12)

# --- Top-20 loadings bars for PC1 & PC2
plot_top <- function(loadings, pc, N = 20, out = NULL) {
  vals <- loadings[[pc]]; genes <- rownames(loadings)
  ord  <- order(abs(vals), decreasing = TRUE)[1:N]
  df   <- data.frame(gene = genes[ord], loading = vals[ord],
                     direction = ifelse(vals[ord] > 0, "Positive", "Negative"))
  p <- ggplot(df, aes(x = reorder(gene, abs(loading)), y = loading, fill = direction)) +
    geom_col() + coord_flip() +
    labs(title = paste("Top", N, "gene loadings –", pc), x = NULL, y = "Loading") +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = 9)) +
    scale_fill_manual(values = c("Positive"="#3d7dbd","Negative"="#bd3d3d"))
  if (!is.null(out)) ggsave(out, p, width = 6.5, height = 6.5, dpi = 300)
  p
}

# --- Save PNGs
dir.create(file.path(root, "figures"), showWarnings = FALSE)
ggsave(file.path(root, "figures/PCA_scree.png"),           p_scree,   width = 7, height = 4.5, dpi = 300)
ggsave(file.path(root, "figures/PCA_cumulative.png"),      p_cum,     width = 7, height = 4.5, dpi = 300)
ggsave(file.path(root, "figures/PCA_PC1_PC2_scatter.png"), p_scatter, width = 6, height = 5, dpi = 300)
ggsave(file.path(root, "figures/PCA_biplot_topgenes.png"), p_biplot,  width = 6.5, height = 6, dpi = 300)
plot_top(loadings, pc1, 20, file.path(root, "figures/PCA_PC1_loadings_top20.png"))
plot_top(loadings, pc2, 20, file.path(root, "figures/PCA_PC2_loadings_top20.png"))

# Optional: show scree + cumulative side-by-side in RStudio
(p_scree | p_cum)
