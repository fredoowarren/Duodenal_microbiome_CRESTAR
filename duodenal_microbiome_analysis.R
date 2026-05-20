################################################################################
#
# Custom code for:
#   Warren, F.J. et al. "Postprandial profiling of the duodenal microbiome
#   reveals the impact of food structure and association with luminal
#   metabolite and gut hormone responses."
#
# This script reproduces all statistical analyses and figures presented in
# the main text (Figures 1-6) and the corresponding supplementary figures.
#
# Sections:
#   1. Setup, libraries, plot theme
#   2. Input data loading and metadata harmonisation
#   3. Figure 1   - Baseline (fasted) microbiome composition and PCoA
#   4. Figure 2   - Shannon alpha diversity over time (LME)
#   5. Figure 3   - Compositional area plot and PCoA following meal
#                   consumption + PERMANOVA (full and pairwise)
#   6. Figure 4   - Differential abundance (arcsine sqrt + LME, ANOVA Type III)
#   7. Figure 5   - Microbe-metabolite-hormone network
#                   (CLR + Spearman + FDR + Louvain modularity + FR layout)
#   8. Figure 6   - D/L amino acid ratio analysis (paired and independent
#                   t-tests)
#   9. Supplementary - per-metabolite postprandial trajectories
#                       (Kruskal-Wallis between groups), time-stratified
#                       PERMANOVA
#
# Required inputs (paths assume the script is run from the project root):
#   data/physeq_hormones.rds      - phyloseq object (MetaPhlAn4 species-level
#                                   relative abundances, sample metadata,
#                                   taxonomy table)
#   data/metabolites_hormones.xlsx - quantified NMR metabolites, bile acids
#                                    and gut hormones (GIP, GLP-1), one row
#                                    per sample, matched by PPT/SampleID/Time
#   data/metabolites_plot.xlsx    - per-metabolite long-format file used for
#                                    Kruskal-Wallis trajectory plots
#   data/amino_acids.xlsx         - chiral-LC-MS D- and L-amino acid
#                                    concentrations at 0 and 60 min
#
# Tested under R 4.5.1.
#
################################################################################

################################################################################
# 1. SETUP
################################################################################

suppressPackageStartupMessages({
  library(phyloseq)
  library(tidyverse)
  library(readxl)
  library(vegan)
  library(lme4)
  library(lmerTest)
  library(compositions)
  library(Hmisc)
  library(igraph)
  library(ggraph)
  library(ggforce)
  library(ggrepel)
  library(ggtext)
  library(patchwork)
  library(RColorBrewer)
  library(viridis)
  library(scales)
  library(pheatmap)
})

set.seed(123)                   # global seed for reproducibility

# ---- Nature-style plotting theme ----
theme_nature <- function(base_size = 8, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title         = element_text(size = 9, face = "bold", hjust = 0),
      axis.title         = element_text(size = 8, face = "bold"),
      axis.text          = element_text(size = 7, colour = "black"),
      legend.title       = element_text(size = 8, face = "bold"),
      legend.text        = element_text(size = 7),
      strip.text         = element_text(size = 8, face = "bold"),
      legend.key.size    = unit(0.4, "cm"),
      legend.background  = element_rect(fill = "white", colour = NA),
      panel.grid.major   = element_line(colour = "grey90", size = 0.25),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(colour = "black", fill = NA, size = 0.5),
      strip.background   = element_rect(fill = "grey95", colour = "black", size = 0.5),
      plot.margin        = margin(5, 5, 5, 5, "pt")
    )
}
theme_set(theme_nature())

# ---- Colour palettes (colour-blind friendly) ----
meal_colors      <- c("rr_flour" = "#1B9E77", "RR_flour" = "#D95F02",
                      "rr_peas"  = "#7570B3", "RR_peas"  = "#E7298A")
genotype_colors  <- c("RR" = "#1B9E77", "rr" = "#D95F02")
food_colors      <- c("Peas" = "#7570B3", "Flour" = "#E7298A")

# Output directory
out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# 2. INPUT DATA
################################################################################

# ---- Phyloseq object (species-level MetaPhlAn4 relative abundances) ----
ps <- readRDS("data/physeq_hormones.rds")

# Drop samples missing SampleID, PPT, or Time
ps <- subset_samples(ps, !is.na(SampleID) & !is.na(PPT) & !is.na(Time))

# Encode meal factors:
#   A = RR peas, B = rr peas, C = RR flour, D = rr flour
sample_data(ps)$Meal_Full <- with(sample_data(ps), case_when(
  SampleID == "A" ~ "RR_peas",
  SampleID == "B" ~ "rr_peas",
  SampleID == "C" ~ "RR_flour",
  SampleID == "D" ~ "rr_flour"
))
sample_data(ps)$Genotype  <- factor(ifelse(sample_data(ps)$SampleID %in% c("A", "C"), "RR", "rr"))
sample_data(ps)$Food_Type <- factor(ifelse(sample_data(ps)$SampleID %in% c("A", "B"), "Peas", "Flour"))
sample_data(ps)$PPT       <- factor(sample_data(ps)$PPT)
sample_data(ps)$SampleID  <- factor(sample_data(ps)$SampleID)
sample_data(ps)$Time      <- as.numeric(as.character(sample_data(ps)$Time))

# Drop empty samples / taxa
ps <- prune_samples(sample_sums(ps) >= 1, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# Helper: best (lowest available) taxonomic name for an OTU
best_tax_name <- function(tax_row) {
  for (lvl in c("Species", "Genus", "Family", "Order", "Class", "Phylum")) {
    nm <- tax_row[[lvl]]
    if (!is.null(nm) && !is.na(nm) && nm != "" && nm != "NA" &&
        !grepl("uncultured|unidentified|unknown", nm, ignore.case = TRUE)) {
      return(as.character(nm))
    }
  }
  NA_character_
}

################################################################################
# 3. FIGURE 1 - BASELINE (FASTED) MICROBIOME COMPOSITION AND PCoA
################################################################################

ps_baseline <- subset_samples(ps, Time == 0)

# ---- Figure 1A: top-30 species stacked bar per participant ----
top30_baseline <- names(sort(taxa_sums(ps_baseline), decreasing = TRUE))[1:30]
ps_baseline_top <- prune_taxa(top30_baseline, ps_baseline)

tax_df_base <- as.data.frame(tax_table(ps_baseline_top))
new_names   <- vapply(seq_len(nrow(tax_df_base)),
                     function(i) best_tax_name(tax_df_base[i, , drop = FALSE]),
                     character(1))
new_names[is.na(new_names)] <- rownames(tax_df_base)[is.na(new_names)]
taxa_names(ps_baseline_top) <- new_names

baseline_long <- psmelt(ps_baseline_top)

fig1a <- ggplot(baseline_long,
                aes(x = factor(PPT), y = Abundance, fill = OTU)) +
  geom_col(position = "fill", width = 0.9) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(30)) +
  facet_wrap(~ Meal_Full, nrow = 1, scales = "free_x") +
  labs(x = "Participant", y = "Relative abundance", fill = "Taxon",
       title = "Figure 1A: Baseline duodenal microbiome composition") +
  theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1),
        legend.text = element_text(size = 6, face = "italic"),
        legend.key.size = unit(0.3, "cm"))

# ---- Figure 1B: PCoA (Bray-Curtis) of baseline samples ----
bc_baseline   <- phyloseq::distance(ps_baseline, method = "bray")
pcoa_baseline <- ordinate(ps_baseline, method = "PCoA", distance = bc_baseline)

# PERMANOVA for inter-individual variability at baseline
permanova_baseline <- adonis2(
  bc_baseline ~ PPT,
  data        = data.frame(sample_data(ps_baseline)),
  permutations = 999
)
cat("\nFigure 1: baseline PERMANOVA (Participant effect):\n")
print(permanova_baseline)

fig1b <- plot_ordination(ps_baseline, pcoa_baseline, color = "PPT") +
  geom_point(size = 1.5, alpha = 0.7) +
  stat_ellipse(aes(group = PPT), type = "t", level = 0.95, size = 0.4) +
  labs(title = "Figure 1B: Baseline PCoA (Bray-Curtis)",
       x = sprintf("PC1 (%.1f%%)", pcoa_baseline$values$Relative_eig[1] * 100),
       y = sprintf("PC2 (%.1f%%)", pcoa_baseline$values$Relative_eig[2] * 100),
       color = "Participant")

ggsave(file.path(out_dir, "Figure1A_baseline_composition.pdf"),
       fig1a, width = 183, height = 100, units = "mm", device = cairo_pdf)
ggsave(file.path(out_dir, "Figure1B_baseline_pcoa.pdf"),
       fig1b, width = 89, height = 89, units = "mm", device = cairo_pdf)

################################################################################
# 4. FIGURE 2 - SHANNON ALPHA DIVERSITY (LINEAR MIXED EFFECTS)
################################################################################

# Shannon index from relative abundance (percentage scale -> proportion)
otu_mat <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)
prop_mat <- otu_mat / 100

shannon_vec <- apply(prop_mat, 2, function(x) {
  x <- x[x > 0]; -sum(x * log(x))
})

alpha_div <- data.frame(
  Shannon   = shannon_vec,
  PPT       = sample_data(ps)$PPT,
  Time      = sample_data(ps)$Time,
  Meal_Full = sample_data(ps)$Meal_Full,
  Genotype  = sample_data(ps)$Genotype,
  Food_Type = sample_data(ps)$Food_Type
)

# Linear mixed-effects model:
#   Shannon ~ Time * Genotype * Food_Type, random intercept for Participant
shannon_lmer <- lmer(
  Shannon ~ Time * Genotype * Food_Type + (1 | PPT),
  data = alpha_div
)
cat("\nFigure 2: Shannon diversity - Type III ANOVA (Satterthwaite):\n")
print(anova(shannon_lmer))

fig2 <- ggplot(alpha_div,
               aes(x = Time, y = Shannon,
                   colour = Meal_Full,
                   group  = interaction(PPT, Meal_Full))) +
  geom_line(alpha = 0.3, size = 0.3) +
  geom_smooth(aes(group = Meal_Full, fill = Meal_Full),
              method = "loess", se = TRUE, size = 0.8, alpha = 0.2) +
  facet_wrap(~ Meal_Full, ncol = 2) +
  scale_x_continuous(breaks = seq(0, 180, 60)) +
  scale_colour_manual(values = meal_colors) +
  scale_fill_manual(values   = meal_colors) +
  labs(title = "Figure 2: Shannon diversity over time",
       x = "Time (min)", y = "Shannon diversity") +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "Figure2_alpha_diversity.pdf"),
       fig2, width = 183, height = 120, units = "mm", device = cairo_pdf)

################################################################################
# 5. FIGURE 3 - COMMUNITY COMPOSITION AND BETA DIVERSITY OVER TIME
################################################################################

# ---- Bray-Curtis distance on full dataset ----
bc_dist  <- phyloseq::distance(ps, method = "bray")
metadata <- data.frame(sample_data(ps))

# ---- Full PERMANOVA model (Methods text: R2 = 0.365, F = 9.97, p = 0.001) ----
permanova_full <- adonis2(
  bc_dist ~ Genotype * Food_Type + PPT,
  data         = metadata,
  permutations = 999
)
cat("\nFigure 3: full PERMANOVA model (Genotype * Food_Type + Participant):\n")
print(permanova_full)

# ---- Pairwise PERMANOVA between meal types with BH correction ----
pairwise_permanova <- function(dist_obj, group_var, meta, perms = 999) {
  groups <- unique(meta[[group_var]])
  comp   <- combn(groups, 2, simplify = FALSE)
  res    <- map_dfr(comp, function(pair) {
    idx <- which(meta[[group_var]] %in% pair)
    d_sub <- as.dist(as.matrix(dist_obj)[idx, idx])
    set.seed(123)
    a <- adonis2(d_sub ~ meta[[group_var]][idx], permutations = perms)
    tibble(Comparison = paste(pair, collapse = " vs "),
           R2 = a$R2[1], F = a$F[1], P = a$`Pr(>F)`[1])
  })
  res$P_adj <- p.adjust(res$P, method = "BH")
  res
}

pairwise_meals <- pairwise_permanova(bc_dist, "Meal_Full", metadata)
cat("\nFigure 3: pairwise PERMANOVA between meals (BH-corrected):\n")
print(pairwise_meals)
write.csv(pairwise_meals,
          file.path(out_dir, "pairwise_meal_permanova.csv"),
          row.names = FALSE)

# ---- Figure 3A: top-20 species area plot, averaged per meal x time ----
top20 <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:20]
ps_top20 <- prune_taxa(top20, ps)

tax_df20  <- as.data.frame(tax_table(ps_top20))
nm20      <- vapply(seq_len(nrow(tax_df20)),
                    function(i) best_tax_name(tax_df20[i, , drop = FALSE]),
                    character(1))
nm20[is.na(nm20)] <- rownames(tax_df20)[is.na(nm20)]
taxa_names(ps_top20) <- nm20

ps_top20_long <- psmelt(ps_top20) %>%
  mutate(Time = as.numeric(as.character(Time))) %>%
  group_by(Meal_Full, Time, OTU) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

fig3a <- ggplot(ps_top20_long,
                aes(x = Time, y = Abundance, fill = OTU)) +
  geom_area(position = "fill") +
  facet_wrap(~ Meal_Full, ncol = 2) +
  scale_x_continuous(breaks = seq(0, 180, 60)) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(20)) +
  labs(title = "Figure 3A: Top-20 species over time",
       x = "Time (min)", y = "Relative abundance", fill = "Taxon") +
  theme(legend.text = element_text(size = 6, face = "italic"),
        legend.key.size = unit(0.3, "cm"))

# ---- Figure 3B: PCoA coloured by meal ----
pcoa <- ordinate(ps, method = "PCoA", distance = bc_dist)
fig3b <- plot_ordination(ps, pcoa, color = "Meal_Full") +
  geom_point(size = 1.5, alpha = 0.7) +
  stat_ellipse(aes(group = Meal_Full), type = "t", level = 0.95, size = 0.4) +
  scale_colour_manual(values = meal_colors) +
  labs(title = "Figure 3B: PCoA by meal type",
       x = sprintf("PC1 (%.1f%%)", pcoa$values$Relative_eig[1] * 100),
       y = sprintf("PC2 (%.1f%%)", pcoa$values$Relative_eig[2] * 100),
       colour = "Meal")

ggsave(file.path(out_dir, "Figure3A_composition_areaplot.pdf"),
       fig3a, width = 183, height = 140, units = "mm", device = cairo_pdf)
ggsave(file.path(out_dir, "Figure3B_pcoa_by_meal.pdf"),
       fig3b, width = 89, height = 89, units = "mm", device = cairo_pdf)

################################################################################
# 6. FIGURE 4 - DIFFERENTIAL ABUNDANCE (TOP-50 TAXA, ARCSINE SQRT + LME)
################################################################################

# Test the 50 most abundant species. Relative abundances are transformed with
# the arcsine square-root transform to stabilise variance, then fitted with:
#   asin(sqrt(p)) ~ Time * Genotype * Food_Type + (1 | PPT)
# and Type III ANOVA (Satterthwaite df) is read out.

otu_df <- as.data.frame(as(otu_table(ps), "matrix"))
if (!taxa_are_rows(ps)) otu_df <- as.data.frame(t(otu_df))

top50 <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:50]

# Map OTU IDs to readable taxonomic names
tax_full <- as.data.frame(tax_table(ps))
otu_label <- sapply(top50, function(otu) {
  if (otu %in% rownames(tax_full)) {
    nm <- best_tax_name(tax_full[otu, , drop = FALSE])
    if (is.na(nm)) otu else nm
  } else otu
})

diff_res <- list()
for (otu in top50) {
  y <- as.numeric(otu_df[otu, ])
  if (all(is.na(y))) next
  td <- data.frame(
    Abundance = y,
    PPT       = sample_data(ps)$PPT,
    Time      = sample_data(ps)$Time,
    Genotype  = sample_data(ps)$Genotype,
    Food_Type = sample_data(ps)$Food_Type
  )
  td$y_trans <- asin(sqrt(td$Abundance / 100))
  m <- tryCatch(
    lmer(y_trans ~ Time * Genotype * Food_Type + (1 | PPT),
         data = td,
         control = lmerControl(optimizer = "bobyqa",
                               optCtrl   = list(maxfun = 20000))),
    error = function(e) NULL
  )
  if (is.null(m)) next
  ao <- anova(m)
  diff_res[[otu]] <- data.frame(
    Taxon   = otu_label[otu],
    OTU_ID  = otu,
    Effect  = rownames(ao),
    F_value = ao[, "F value"],
    P_value = ao[, "Pr(>F)"]
  )
}
diff_res <- bind_rows(diff_res) %>%
  filter(!is.na(P_value)) %>%
  arrange(P_value)

write.csv(diff_res, file.path(out_dir, "differential_abundance.csv"),
          row.names = FALSE)

# ---- Figure 4A: dot plot of top-15 significant effects ----
fig4a_data <- diff_res %>%
  filter(P_value < 0.05) %>%
  head(15) %>%
  mutate(neg_log_p     = -log10(P_value),
         Effect_clean  = gsub(":", " x ", Effect))

fig4a <- ggplot(fig4a_data,
                aes(x = neg_log_p,
                    y = reorder(Taxon, neg_log_p),
                    colour = Effect_clean, size = F_value)) +
  geom_point(alpha = 0.85) +
  geom_vline(xintercept = -log10(0.05),
             linetype = "dashed", colour = "grey50", size = 0.3) +
  scale_colour_manual(values = brewer.pal(8, "Dark2")) +
  scale_size_continuous(range = c(2, 5)) +
  labs(title = "Figure 4A: Differentially abundant taxa",
       x = expression(-log[10](italic(P))), y = NULL,
       colour = "Effect", size = expression(italic(F))) +
  theme(axis.text.y = element_text(face = "italic", size = 6))

# ---- Figure 4B: trajectories of the 6 most significant taxa ----
top6_otu <- unique(diff_res$OTU_ID)[1:6]
top6_long <- map_dfr(top6_otu, function(otu) {
  data.frame(
    Taxon     = otu_label[otu],
    Abundance = as.numeric(otu_df[otu, ]),
    Time      = sample_data(ps)$Time,
    Meal_Full = sample_data(ps)$Meal_Full,
    PPT       = sample_data(ps)$PPT
  )
})

fig4b <- ggplot(top6_long,
                aes(x = Time, y = Abundance,
                    colour = Meal_Full,
                    group  = interaction(PPT, Meal_Full))) +
  geom_line(alpha = 0.3, size = 0.3) +
  geom_smooth(aes(group = Meal_Full, fill = Meal_Full),
              method = "loess", se = TRUE, size = 0.6, alpha = 0.2) +
  facet_wrap(~ Taxon, ncol = 3, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 180, 60)) +
  scale_colour_manual(values = meal_colors) +
  scale_fill_manual(values   = meal_colors) +
  labs(title = "Figure 4B: Trajectories of top differential taxa",
       x = "Time (min)", y = "Relative abundance (%)") +
  theme(strip.text = element_text(size = 6, face = "italic"))

ggsave(file.path(out_dir, "Figure4A_diff_abundance_dotplot.pdf"),
       fig4a, width = 120, height = 100, units = "mm", device = cairo_pdf)
ggsave(file.path(out_dir, "Figure4B_diff_abundance_trajectories.pdf"),
       fig4b, width = 183, height = 120, units = "mm", device = cairo_pdf)

################################################################################
# 7. FIGURE 5 - MICROBIOME - METABOLITE - HORMONE NETWORK
#    CLR-transformed microbial abundances + z-scored metabolites/hormones,
#    pairwise Spearman correlations (Hmisc::rcorr), Benjamini-Hochberg FDR,
#    |rho| > 0.3 and q < 0.05, Louvain community detection on |rho|,
#    Fruchterman-Reingold layout, community hulls via ggforce.
################################################################################

# ---- Load metabolite / hormone data and match to microbiome samples ----
metabolites <- read_excel("data/metabolites_hormones.xlsx") %>%
  tibble::column_to_rownames("sample")

# Microbial abundances at species level, CLR-transformed
microbe_raw  <- t(otu_table(ps))
microbe_raw  <- microbe_raw[, colSums(microbe_raw) > 0]
microbe_clr  <- as.data.frame(compositions::clr(microbe_raw + 1))

tax_lookup   <- as.data.frame(tax_table(ps))
species_map  <- tax_lookup$Species
names(species_map) <- rownames(tax_lookup)
colnames(microbe_clr) <- species_map[colnames(microbe_clr)]

# z-score metabolite / hormone concentrations
metabolite_z <- as.data.frame(scale(metabolites))

# Align on shared samples
common_samples <- intersect(rownames(microbe_clr), rownames(metabolite_z))
microbe_clr    <- microbe_clr[common_samples, , drop = FALSE]
metabolite_z   <- metabolite_z[common_samples, , drop = FALSE]
combined       <- cbind(microbe_clr, metabolite_z)

# ---- Pairwise Spearman correlations ----
cor_obj  <- Hmisc::rcorr(as.matrix(combined), type = "spearman")
cor_mat  <- cor_obj$r
p_mat    <- cor_obj$P

microbe_names    <- colnames(microbe_clr)
metabolite_names <- colnames(metabolite_z)
cor_sub <- cor_mat[microbe_names, metabolite_names, drop = FALSE]
p_sub   <- p_mat  [microbe_names, metabolite_names, drop = FALSE]

# Benjamini-Hochberg FDR on the microbe x metabolite submatrix
p_adj   <- matrix(p.adjust(as.vector(p_sub), method = "fdr"),
                  nrow = nrow(p_sub), dimnames = dimnames(p_sub))

sig_threshold <- 0.05
cor_threshold <- 0.3
sig_idx <- which(abs(cor_sub) > cor_threshold & p_adj < sig_threshold,
                 arr.ind = TRUE)

hormone_names <- metabolite_names[grepl("GIP|GLP", metabolite_names, ignore.case = TRUE)]

if (nrow(sig_idx) == 0) {
  message("Figure 5: no significant correlations met the thresholds; skipping.")
} else {

  edges <- data.frame(
    from        = rownames(cor_sub)[sig_idx[, 1]],
    to          = colnames(cor_sub)[sig_idx[, 2]],
    correlation = cor_sub[sig_idx],
    p_adj       = p_adj[sig_idx],
    stringsAsFactors = FALSE
  )

  cat(sprintf("\nFigure 5: %d significant microbe-metabolite/hormone edges\n",
              nrow(edges)))

  # Metabolite class lookup (from Methods)
  scfa_names       <- c("Acetate", "Propionate", "Formate")
  amino_acid_names <- c("Valine", "Alanine", "Methionine", "Lysine", "Tyrosine",
                        "Phenylalanine", "Tryptophan", "Threonine")
  bile_acid_names  <- c("BA1", "BA2", "BA3", "BA4", "BA6", "BA7")

  g <- igraph::graph_from_data_frame(edges, directed = FALSE)

  V(g)$type <- with(igraph::as_data_frame(g, what = "vertices"),
                    ifelse(name %in% microbe_names,    "Microbe",
                    ifelse(name %in% hormone_names,    "Hormone",
                    ifelse(name %in% scfa_names,       "SCFA",
                    ifelse(name %in% amino_acid_names, "Amino Acid",
                    ifelse(name %in% bile_acid_names,  "Bile Acid",
                                                       "Other Metabolite"))))))
  V(g)$degree <- igraph::degree(g)

  # Louvain modularity on absolute correlations
  set.seed(123)
  louvain         <- igraph::cluster_louvain(g, weights = abs(E(g)$correlation))
  V(g)$community  <- louvain$membership
  modularity_val  <- modularity(g, louvain$membership,
                                weights = abs(E(g)$correlation))

  # Per-node table for supplementary export
  node_table <- data.frame(
    node        = V(g)$name,
    type        = V(g)$type,
    community   = V(g)$community,
    degree      = V(g)$degree,
    betweenness = igraph::betweenness(g),
    closeness   = igraph::closeness(g)
  ) %>%
    arrange(community, desc(degree))

  write.csv(edges,      file.path(out_dir, "network_edges.csv"),      row.names = FALSE)
  write.csv(node_table, file.path(out_dir, "network_nodes.csv"),      row.names = FALSE)

  # ---- Fruchterman-Reingold layout ----
  set.seed(123)
  layout    <- ggraph::create_layout(g, layout = "fr")
  layout_df <- as.data.frame(layout)

  node_colors <- c("Microbe"          = "#5A6C7D",
                   "SCFA"             = "#F39C12",
                   "Amino Acid"       = "#27AE60",
                   "Bile Acid"        = "#8E44AD",
                   "Hormone"          = "#F1C40F",
                   "Other Metabolite" = "#7F8C8D")

  cluster_colors <- viridis::viridis(length(unique(V(g)$community)),
                                     option = "plasma")

  fig5 <- ggraph(layout) +
    ggforce::geom_mark_hull(
      data = layout_df,
      aes(x = x, y = y, group = community, fill = factor(community)),
      alpha = 0.12, colour = NA, expand = unit(4, "mm")
    ) +
    geom_edge_link(aes(colour = correlation, width = abs(correlation)),
                   alpha = 0.6, show.legend = TRUE) +
    geom_node_point(aes(colour = type, size = degree), alpha = 0.95) +
    ggrepel::geom_text_repel(
      data = layout_df,
      aes(x = x, y = y, label = name, colour = type),
      size = 3, box.padding = 0.3, point.padding = 0.3,
      segment.colour = "grey70", max.overlaps = Inf, show.legend = FALSE
    ) +
    scale_edge_colour_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
                                midpoint = 0, name = expression(rho)) +
    scale_colour_manual(values = node_colors, name = "Node type") +
    scale_fill_manual(values   = cluster_colors, guide = "none") +
    scale_size(range = c(3, 10), guide = "none") +
    scale_edge_width(range = c(0.2, 1.5), guide = "none") +
    theme_void(base_family = "Helvetica") +
    theme(legend.position = "right",
          plot.background  = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA))

  ggsave(file.path(out_dir, "Figure5_network.pdf"),
         fig5, width = 183, height = 150, units = "mm", device = cairo_pdf)

  cat(sprintf("Figure 5: %d nodes, %d edges, %d Louvain communities (modularity = %.3f)\n",
              vcount(g), ecount(g),
              length(unique(V(g)$community)),
              modularity_val))
}

################################################################################
# 8. FIGURE 6 - D/L AMINO ACID RATIOS
#    Paired t-tests for the time effect within each meal arm, independent
#    t-tests for the meal effect at each timepoint (0 and 60 min).
################################################################################

aa_data <- read_excel("data/amino_acids.xlsx")

d_cols     <- colnames(aa_data)[grepl("^D ", colnames(aa_data))]
amino_acids <- sub("^D ", "", d_cols)

dl_ratios <- aa_data %>%
  select(Participant, Timepoint, Randomization, all_of(d_cols)) %>%
  mutate(across(all_of(d_cols), as.numeric))

for (aa in amino_acids) {
  dl_ratios[[paste0("DL_", aa)]] <- aa_data[[paste0("D ", aa)]] /
                                    aa_data[[paste0("L ", aa)]]
}

ratio_cols <- grep("^DL_", colnames(dl_ratios), value = TRUE)

dl_long <- dl_ratios %>%
  select(Participant, Timepoint, Randomization, all_of(ratio_cols)) %>%
  pivot_longer(cols = all_of(ratio_cols),
               names_to  = "Amino_Acid",
               values_to = "DL_Ratio") %>%
  mutate(Amino_Acid = sub("^DL_", "", Amino_Acid),
         Treatment  = ifelse(Randomization == "A", "RR peas", "rr peas"),
         Timepoint  = factor(Timepoint, levels = c("0", "60"))) %>%
  filter(is.finite(DL_Ratio) & !is.na(DL_Ratio))

# Helper for paired t-test on (T0, T60) within an arm
paired_t <- function(d, treatment) {
  w <- d %>%
    filter(Treatment == treatment) %>%
    group_by(Participant, Timepoint) %>%
    summarise(DL_Ratio = mean(DL_Ratio, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Timepoint, values_from = DL_Ratio,
                names_prefix = "T") %>%
    filter(!is.na(T0) & !is.na(T60))
  if (nrow(w) < 2) return(NA_real_)
  tryCatch(t.test(w$T0, w$T60, paired = TRUE)$p.value, error = function(e) NA_real_)
}

# Helper for independent t-test (RR vs rr) at a single timepoint
indep_t <- function(d, tp) {
  s <- d %>%
    filter(Timepoint == tp) %>%
    group_by(Participant, Treatment) %>%
    summarise(DL_Ratio = mean(DL_Ratio, na.rm = TRUE), .groups = "drop")
  if (nrow(s) < 4) return(NA_real_)
  tryCatch(t.test(DL_Ratio ~ Treatment, data = s)$p.value,
           error = function(e) NA_real_)
}

dl_stats <- dl_long %>%
  group_by(Amino_Acid) %>%
  group_modify(~ tibble(
    time_RR  = paired_t(.x, "RR peas"),
    time_rr  = paired_t(.x, "rr peas"),
    trt_T0   = indep_t(.x, "0"),
    trt_T60  = indep_t(.x, "60")
  )) %>%
  ungroup()

write.csv(dl_stats, file.path(out_dir, "DL_ratio_statistics.csv"),
          row.names = FALSE)

p_to_stars <- function(p) {
  case_when(is.na(p) ~ "",
            p < 0.001 ~ "***",
            p < 0.01  ~ "**",
            p < 0.05  ~ "*",
            TRUE      ~ "ns")
}

dl_summary <- dl_long %>%
  group_by(Amino_Acid, Treatment, Timepoint) %>%
  summarise(mean_ratio = mean(DL_Ratio, na.rm = TRUE),
            se         = sd(DL_Ratio, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

fig6a <- ggplot(dl_summary,
                aes(x = Timepoint, y = mean_ratio,
                    colour = Treatment, group = Treatment)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_ratio - se, ymax = mean_ratio + se),
                width = 0.2, linewidth = 0.6) +
  facet_wrap(~ Amino_Acid, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = c("RR peas" = "#E78AC3",
                                 "rr peas" = "#8DA0CB")) +
  labs(title = "Figure 6A: D/L amino acid ratios over time",
       x = "Time (min)", y = "D/L ratio")

p_matrix <- dl_stats %>%
  pivot_longer(cols = -Amino_Acid,
               names_to  = "Comparison",
               values_to = "p_value") %>%
  mutate(sig = p_to_stars(p_value),
         Comparison = recode(Comparison,
                             "time_RR" = "Time (RR)",
                             "time_rr" = "Time (rr)",
                             "trt_T0"  = "Trt (0 min)",
                             "trt_T60" = "Trt (60 min)"))

fig6b <- ggplot(p_matrix,
                aes(x = Comparison, y = Amino_Acid,
                    fill = -log10(p_value))) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sig), size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "white", mid = "lightyellow", high = "firebrick",
                       midpoint = 1, na.value = "grey90",
                       name = expression(-log[10](italic(p)))) +
  labs(title = "Figure 6B: p-value summary", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig6 <- fig6a / fig6b + plot_layout(heights = c(3, 1))

ggsave(file.path(out_dir, "Figure6_DL_amino_acid_ratios.pdf"),
       fig6, width = 183, height = 240, units = "mm", device = cairo_pdf)

################################################################################
# 9. SUPPLEMENTARY ANALYSES
################################################################################

# ---- 9A. Per-metabolite postprandial trajectories with Kruskal-Wallis ----
#         between-group testing at each timepoint (for Supplementary Figures
#         showing per-metabolite curves).

metab_plot_data <- read_excel("data/metabolites_plot.xlsx") %>%
  pivot_longer(cols      = -(Participant:Group),
               names_to  = "Metabolite",
               values_to = "Value") %>%
  mutate(Time = as.numeric(Time))

plot_metabolite <- function(metabolite_name, data = metab_plot_data) {

  d_sub <- filter(data, Metabolite == metabolite_name)

  d_sum <- d_sub %>%
    group_by(Group, Time) %>%
    summarise(mean_value = mean(Value, na.rm = TRUE),
              sem        = sd(Value, na.rm = TRUE) / sqrt(n()),
              .groups = "drop")

  # Per-timepoint Kruskal-Wallis between meal groups
  pvals <- d_sub %>%
    group_by(Time) %>%
    summarise(p = tryCatch({
      if (length(unique(Group[!is.na(Value)])) < 2) NA_real_
      else kruskal.test(Value ~ Group)$p.value
    }, error = function(e) NA_real_),
    .groups = "drop") %>%
    mutate(stars = p_to_stars(p))

  y_max <- max(d_sub$Value, na.rm = TRUE)
  pvals$ypos <- y_max + y_max * 0.15

  ggplot() +
    geom_line(data = d_sub,
              aes(x = Time, y = Value,
                  group = Participant, colour = Group),
              alpha = 0.25, size = 0.5, show.legend = FALSE) +
    geom_line(data = d_sum,
              aes(x = Time, y = mean_value,
                  colour = Group, group = Group),
              size = 1.2) +
    geom_errorbar(data = d_sum,
                  aes(x = Time,
                      ymin = mean_value - sem,
                      ymax = mean_value + sem,
                      colour = Group),
                  width = max(1, diff(range(d_sum$Time, na.rm = TRUE)) / 50),
                  size = 0.6) +
    geom_point(data = d_sum,
               aes(x = Time, y = mean_value, colour = Group),
               size = 2.5) +
    geom_text(data = filter(pvals, stars != "ns" & stars != ""),
              aes(x = Time, y = ypos, label = stars),
              inherit.aes = FALSE, size = 4) +
    scale_colour_brewer(palette = "Dark2") +
    scale_x_continuous(breaks = sort(unique(d_sub$Time))) +
    labs(title = metabolite_name,
         x = "Time (min)", y = "Concentration (a.u.)",
         colour = "Group")
}

# Example: render one panel per metabolite
all_metabs <- unique(metab_plot_data$Metabolite)
for (m in all_metabs) {
  ggsave(file.path(out_dir, paste0("Sup_metabolite_", m, ".pdf")),
         plot_metabolite(m),
         width = 100, height = 80, units = "mm", device = cairo_pdf)
}

# ---- 9B. Time-stratified PERMANOVA (when do groups separate?) ----
time_strat <- map_dfr(sort(unique(metadata$Time)), function(tp) {
  idx <- which(metadata$Time == tp)
  if (length(idx) < 10) return(NULL)
  d_sub <- as.dist(as.matrix(bc_dist)[idx, idx])
  m_sub <- metadata[idx, ]
  set.seed(123)
  pg <- adonis2(d_sub ~ Genotype,  data = m_sub, permutations = 999)
  set.seed(123)
  pf <- adonis2(d_sub ~ Food_Type, data = m_sub, permutations = 999)
  bind_rows(
    tibble(Time = tp, Comparison = "Genotype",
           R2 = pg$R2[1], P = pg$`Pr(>F)`[1]),
    tibble(Time = tp, Comparison = "Food_Type",
           R2 = pf$R2[1], P = pf$`Pr(>F)`[1])
  )
})
time_strat$P_adj <- p.adjust(time_strat$P, method = "BH")
write.csv(time_strat, file.path(out_dir, "time_stratified_permanova.csv"),
          row.names = FALSE)

cat("\nAll analyses complete. Outputs are in:", normalizePath(out_dir), "\n")

################################################################################
# End of script
################################################################################
