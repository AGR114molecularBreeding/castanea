#########################################################################################################
# 1_Comparative_analysis.R
# Comparative analysis of FANTASIA vs Blast2GO annotations in C. sativa
#
# Produces: Figure 1, Figure 2 (A semantic similarity vs null model; B specificity),
#           Supplementary Figure S1 (Venn: proteins and GO vocabulary),
#           and all quantitative results reported in the manuscript
#           (coverage, vocabulary overlap, richness, exact agreement, semantic
#            similarity with null model, specificity, effect-size CIs and
#            sensitivity analysis).
#
# Run this BEFORE the Python taxonomic classification step: at the end it exports
# results/go_terms_to_classify.tsv, the input for 2_taxonomic_classification.py
#
# Postigo-Luque et al.
#########################################################################################################

# Set the working directory to the repository folder before running, e.g.:
# setwd("/path/to/Functional annotation comparative analysis FANTASIA vs Blast2Go")

source("0_Load_data.R")
library(eulerr)
library(ontologyIndex)
library(ontologySimilarity)
library(boot)
library(patchwork)

dir.create("results", showWarnings = FALSE)

# ============================================================
# A.1  Annotation coverage
# ============================================================
cat("\n=== A.1 Annotation coverage ===\n")
cat("Total proteins:       ", n_total, "\n")
cat("FANTASIA annotated:   ", n_fantasia, "(", round(n_fantasia/n_total*100, 1), "%)\n")
cat("Blast2GO annotated:   ", n_homologia, "(", round(n_homologia/n_total*100, 1), "%)\n")
cat("Both methods:         ", n_both, "\n")
cat("Only FANTASIA:        ", length(only_fantasia), "\n")
cat("Only Blast2GO:        ", length(only_homologia), "\n")

# ------------------------------------------------------------
# Supplementary Figure S1A: Venn (proteins)
# ------------------------------------------------------------
venn_prot <- euler(c(
  "FANTASIA"            = length(only_fantasia),
  "Blast2GO"            = length(only_homologia),
  "FANTASIA&Blast2GO"   = n_both
))
png("results/Figure_S1A_venn_proteins.png", width = 9, height = 9,
    units = "cm", res = 300)
print(plot(venn_prot, quantities = TRUE,
           fills = c("#4BACC6", "#F79646"),
           main = "Annotated proteins"))
dev.off()
cat("Figure S1A (Venn proteins) saved.\n")

# ============================================================
# A.2  GO vocabulary diversity
# ============================================================
cat("\n=== A.2 GO vocabulary diversity ===\n")

fan_bp <- fantasia_go_largo[dominio_map[fantasia_go_largo] == "BP" & !is.na(dominio_map[fantasia_go_largo])]
fan_mf <- fantasia_go_largo[dominio_map[fantasia_go_largo] == "MF" & !is.na(dominio_map[fantasia_go_largo])]
fan_cc <- fantasia_go_largo[dominio_map[fantasia_go_largo] == "CC" & !is.na(dominio_map[fantasia_go_largo])]
hom_bp <- homologia_go_largo[dominio_map[homologia_go_largo] == "BP" & !is.na(dominio_map[homologia_go_largo])]
hom_mf <- homologia_go_largo[dominio_map[homologia_go_largo] == "MF" & !is.na(dominio_map[homologia_go_largo])]
hom_cc <- homologia_go_largo[dominio_map[homologia_go_largo] == "CC" & !is.na(dominio_map[homologia_go_largo])]

calc_vocab <- function(set_fan, set_hom, label) {
  n_fan   <- length(set_fan)
  n_hom   <- length(set_hom)
  n_comun <- length(intersect(set_fan, set_hom))
  n_union <- length(union(set_fan, set_hom))
  data.frame(
    Domain             = label,
    FANTASIA_unique    = n_fan,
    Blast2Go_unique    = n_hom,
    Shared             = n_comun,
    Exclusive_FANTASIA = length(setdiff(set_fan, set_hom)),
    Exclusive_Blast2Go = length(setdiff(set_hom, set_fan)),
    Jaccard            = round(n_comun / n_union, 4),
    Dice               = round(2 * n_comun / (n_fan + n_hom), 4),
    Overlap            = round(n_comun / min(n_fan, n_hom), 4)
  )
}

vocab_table <- rbind(
  calc_vocab(fantasia_go_largo, homologia_go_largo, "Total"),
  calc_vocab(fan_bp, hom_bp, "BP"),
  calc_vocab(fan_mf, hom_mf, "MF"),
  calc_vocab(fan_cc, hom_cc, "CC")
)
print(vocab_table)
write.table(vocab_table, "results/A2_vocabulary_comparison.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ------------------------------------------------------------
# Supplementary Figure S1B: Venn (GO vocabulary, Total)
# ------------------------------------------------------------
venn_go <- euler(c(
  "FANTASIA"          = vocab_table$Exclusive_FANTASIA[1],
  "Blast2GO"          = vocab_table$Exclusive_Blast2Go[1],
  "FANTASIA&Blast2GO" = vocab_table$Shared[1]
))
png("results/Figure_S1B_venn_GO.png", width = 9, height = 9,
    units = "cm", res = 300)
print(plot(venn_go, quantities = TRUE,
           fills = c("#4BACC6", "#F79646"),
           main = "GO vocabulary"))
dev.off()
cat("Figure S1B (Venn GO) saved.\n")

# ============================================================
# A.3  Annotation richness per protein  (Figure 1)
# ============================================================
cat("\n=== A.3 Annotation richness ===\n")

go_counts_fantasia <- data.frame(
  ProteinID = fantasia$ProteinID,
  n_GO      = sapply(strsplit(fantasia$GOterms, ","),
                     function(x) length(unique(trimws(x)))),
  Method    = "FANTASIA"
)
go_counts_homologia <- aggregate(Annotation.GO.ID ~ Sequence.Name,
                                 data = homologia,
                                 FUN = function(x) length(unique(x)))
colnames(go_counts_homologia) <- c("ProteinID", "n_GO")
go_counts_homologia$Method <- "Blast2GO"

go_counts        <- rbind(go_counts_fantasia[, c("ProteinID","n_GO","Method")],
                          go_counts_homologia[, c("ProteinID","n_GO","Method")])
go_counts$n_GO   <- as.numeric(go_counts$n_GO)
go_counts$Method <- factor(go_counts$Method, levels = c("FANTASIA", "Blast2GO"))

# Mann-Whitney U test + rank-biserial effect size (absolute value)
fantasia_go  <- go_counts$n_GO[go_counts$Method == "FANTASIA"]
homologia_go <- go_counts$n_GO[go_counts$Method == "Blast2GO"]
test_mw      <- wilcox.test(fantasia_go, homologia_go, exact = FALSE)
n1 <- length(fantasia_go)
n2 <- length(homologia_go)
r_biserial <- abs(1 - (2 * as.numeric(test_mw$statistic)) / (as.numeric(n1) * as.numeric(n2)))

cat("  FANTASIA median:", median(fantasia_go), "\n")
cat("  Blast2GO median:", median(homologia_go), "\n")
cat("  Mann-Whitney p =", test_mw$p.value, "\n")
cat("  Rank-biserial |r| =", round(r_biserial, 4), "\n")

# Figure 1: Violin plot
medians <- aggregate(n_GO ~ Method, data = go_counts, FUN = median)

p_fig1 <- ggplot(go_counts, aes(x = Method, y = n_GO, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA, coef = Inf, color = "grey30") +
  geom_text(data = medians,
            aes(x = Method, y = n_GO, label = paste0("Md=", n_GO)),
            vjust = -0.8, size = 3, fontface = "bold", color = "grey20") +
  scale_fill_manual(values = c("FANTASIA" = "#4BACC6", "Blast2GO" = "#F79646")) +
  scale_y_log10(breaks = c(1, 2, 5, 10, 25, 50, 100, 200)) +
  annotate("text", x = 1.5, y = 100,
           label = paste0("Mann-Whitney U\np < 2.2e-16\n|r| = ",
                          round(r_biserial, 2)),
           size = 2.8, hjust = 0.5, color = "grey20") +
  labs(x = NULL, y = "Number of GO terms per protein (log scale)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        axis.title = element_text(size = 9),
        axis.text  = element_text(color = "grey20", size = 8))
ggsave("results/Figure_1_violin_GO_per_protein.jpg", p_fig1,
       width = 9, height = 9, units = "cm", dpi = 300, bg = "white")
cat("Figure 1 saved.\n")

# ============================================================
# B.1.1  Exact agreement between methods
# ============================================================
cat("\n=== B.1.1 Exact GO term agreement ===\n")

acuerdo_exacto <- sapply(proteinas_comunes, function(prot) {
  length(intersect(fantasia_comun[[prot]], homologia_comun[[prot]])) > 0
})
n_acuerdo    <- sum(acuerdo_exacto)
perc_acuerdo <- n_acuerdo / length(proteinas_comunes) * 100

cat("  Common proteins:            ", length(proteinas_comunes), "\n")
cat("  Sharing at least 1 GO term: ", n_acuerdo, "\n")
cat("  Percentage:                 ", round(perc_acuerdo, 1), "%\n")

# ============================================================
# B.1.2  Semantic similarity BMA/Resnik  (observed)
# ============================================================
cat("\n=== B.1.2 Semantic similarity BMA/Resnik (observed) ===\n")

go_obo <- get_ontology(GO_OBO_FILE, extract_tags = "everything")
ic     <- descendants_IC(go_obo)
namespace_map <- setNames(go_obo$namespace, go_obo$id)

resnik_sim <- function(t1, t2, go_obo, ic) {
  if (is.na(ic[t1]) || is.na(ic[t2])) return(0)
  anc1   <- c(t1, get_ancestors(go_obo, t1))
  anc2   <- c(t2, get_ancestors(go_obo, t2))
  common <- intersect(anc1, anc2)
  if (length(common) == 0) return(0)
  max(ic[common], na.rm = TRUE)
}

calc_bma <- function(set1, set2, go_obo, ic) {
  set1 <- set1[set1 %in% names(ic) & !is.na(ic[set1])]
  set2 <- set2[set2 %in% names(ic) & !is.na(ic[set2])]
  if (length(set1) == 0 || length(set2) == 0) return(NA)
  best_1to2 <- sapply(set1, function(t1)
    max(sapply(set2, function(t2) resnik_sim(t1, t2, go_obo, ic)), na.rm = TRUE))
  best_2to1 <- sapply(set2, function(t2)
    max(sapply(set1, function(t1) resnik_sim(t1, t2, go_obo, ic)), na.rm = TRUE))
  (mean(best_1to2) + mean(best_2to1)) / 2
}

filter_by_domain <- function(go_terms, domain) {
  go_terms[go_terms %in% names(namespace_map) & namespace_map[go_terms] == domain]
}

sim_results <- data.frame(
  ProteinID = proteinas_comunes,
  sim_BP = NA_real_, sim_MF = NA_real_, sim_CC = NA_real_
)

cat("Calculating observed BMA/Resnik for", length(proteinas_comunes),
    "proteins (may take several hours)...\n")
start_time <- Sys.time()
for (i in seq_along(proteinas_comunes)) {
  prot <- proteinas_comunes[i]
  fan  <- fantasia_comun[[prot]]
  hom  <- homologia_comun[[prot]]
  sim_results$sim_BP[i] <- calc_bma(filter_by_domain(fan, "biological_process"),
                                    filter_by_domain(hom, "biological_process"),
                                    go_obo, ic)
  sim_results$sim_MF[i] <- calc_bma(filter_by_domain(fan, "molecular_function"),
                                    filter_by_domain(hom, "molecular_function"),
                                    go_obo, ic)
  sim_results$sim_CC[i] <- calc_bma(filter_by_domain(fan, "cellular_component"),
                                    filter_by_domain(hom, "cellular_component"),
                                    go_obo, ic)
  if (i %% 2000 == 0) {
    elapsed   <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    remaining <- (length(proteinas_comunes) - i) / (i / elapsed)
    cat("  Observed:", i, "/", length(proteinas_comunes),
        "| Elapsed:", round(elapsed, 1), "min",
        "| Remaining:", round(remaining, 1), "min\n")
  }
}
cat("Observed calculation complete.\n")

for (dom in c("BP", "MF", "CC")) {
  vals <- sim_results[[paste0("sim_", dom)]]
  vals <- vals[!is.na(vals)]
  cat("  ", dom, "- Median:", round(median(vals), 2),
      "| Mean:", round(mean(vals), 2), "\n")
}
write.table(sim_results, "results/B1_semantic_similarity_BMA_Resnik.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ============================================================
# B.1.3  Null model for semantic similarity  (Reviewer 1, comment 6)
# Shuffle the pairing between FANTASIA and Blast2GO annotations across
# proteins, preserving the GO graph structure and term frequencies.
# ============================================================
cat("\n=== B.1.3 Null model for semantic similarity ===\n")

set.seed(123)
shuffled_order <- sample(seq_along(proteinas_comunes))

sim_null <- data.frame(
  ProteinID = proteinas_comunes,
  sim_BP = NA_real_, sim_MF = NA_real_, sim_CC = NA_real_
)
cat("Calculating null BMA/Resnik (similar runtime to the observed analysis)...\n")
start_time <- Sys.time()
for (i in seq_along(proteinas_comunes)) {
  fan <- fantasia_comun[[proteinas_comunes[i]]]               # FANTASIA real
  hom <- homologia_comun[[proteinas_comunes[shuffled_order[i]]]]  # Blast2GO shuffled
  sim_null$sim_BP[i] <- calc_bma(filter_by_domain(fan, "biological_process"),
                                 filter_by_domain(hom, "biological_process"),
                                 go_obo, ic)
  sim_null$sim_MF[i] <- calc_bma(filter_by_domain(fan, "molecular_function"),
                                 filter_by_domain(hom, "molecular_function"),
                                 go_obo, ic)
  sim_null$sim_CC[i] <- calc_bma(filter_by_domain(fan, "cellular_component"),
                                 filter_by_domain(hom, "cellular_component"),
                                 go_obo, ic)
  if (i %% 2000 == 0) {
    elapsed   <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    remaining <- (length(proteinas_comunes) - i) / (i / elapsed)
    cat("  Null:", i, "/", length(proteinas_comunes),
        "| Elapsed:", round(elapsed, 1), "min",
        "| Remaining:", round(remaining, 1), "min\n")
  }
}
cat("Null model complete.\n")
write.table(sim_null, "results/B1_semantic_similarity_NULL.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Observed vs Null semantic similarity ===\n")
for (dom in c("BP", "MF", "CC")) {
  real <- sim_results[[paste0("sim_", dom)]]; real <- real[!is.na(real)]
  null <- sim_null[[paste0("sim_", dom)]];    null <- null[!is.na(null)]
  wt <- wilcox.test(real, null)
  cat(" ", dom,
      "| Observed median:", round(median(real), 2),
      "| Null median:", round(median(null), 2),
      "| Wilcoxon p =", format(wt$p.value, digits = 3), "\n")
}

# ============================================================
# B.2  Annotation specificity
# ============================================================
cat("\n=== B.2 Annotation specificity ===\n")

# B.2.1 Corpus-based IC
corpus_fantasia  <- trimws(unlist(strsplit(fantasia$GOterms, ",")))
corpus_homologia <- homologia$Annotation.GO.ID
corpus_total     <- c(corpus_fantasia, corpus_homologia)
term_counts      <- table(corpus_total)
ic_corpus        <- -log2(as.numeric(term_counts) / length(corpus_total))
names(ic_corpus) <- names(term_counts)

# B.2.2 GO term depth
all_go_terms <- union(
  unique(trimws(unlist(strsplit(fantasia$GOterms, ",")))),
  unique(homologia$Annotation.GO.ID)
)
depth_map_full <- sapply(all_go_terms, function(go_id) {
  tryCatch({ length(get_ancestors(go_obo, go_id)) }, error = function(e) NA)
})
names(depth_map_full) <- all_go_terms

# B.2.3 Max IC and depth per protein (common proteins)
max_ic_corpus <- data.frame(ProteinID = proteinas_comunes,
                            max_IC_fan = NA_real_, max_IC_hom = NA_real_)
max_ic_topo   <- data.frame(ProteinID = proteinas_comunes,
                            max_IC_fan = NA_real_, max_IC_hom = NA_real_)
max_depth     <- data.frame(ProteinID = proteinas_comunes,
                            max_depth_fan = NA_real_, max_depth_hom = NA_real_)

for (i in seq_along(proteinas_comunes)) {
  prot    <- proteinas_comunes[i]
  fan_gos <- fantasia_comun[[prot]]
  hom_gos <- homologia_comun[[prot]]

  fan_ic_c <- ic_corpus[fan_gos[fan_gos %in% names(ic_corpus)]]
  hom_ic_c <- ic_corpus[hom_gos[hom_gos %in% names(ic_corpus)]]
  max_ic_corpus$max_IC_fan[i] <- ifelse(length(fan_ic_c) > 0, max(fan_ic_c, na.rm = TRUE), NA)
  max_ic_corpus$max_IC_hom[i] <- ifelse(length(hom_ic_c) > 0, max(hom_ic_c, na.rm = TRUE), NA)

  fan_ic_t <- ic[fan_gos[fan_gos %in% names(ic)]]
  hom_ic_t <- ic[hom_gos[hom_gos %in% names(ic)]]
  max_ic_topo$max_IC_fan[i] <- ifelse(length(fan_ic_t) > 0, max(fan_ic_t, na.rm = TRUE), NA)
  max_ic_topo$max_IC_hom[i] <- ifelse(length(hom_ic_t) > 0, max(hom_ic_t, na.rm = TRUE), NA)

  fan_d <- depth_map_full[fan_gos[fan_gos %in% names(depth_map_full)]]
  hom_d <- depth_map_full[hom_gos[hom_gos %in% names(depth_map_full)]]
  max_depth$max_depth_fan[i] <- ifelse(length(fan_d) > 0, max(fan_d, na.rm = TRUE), NA)
  max_depth$max_depth_hom[i] <- ifelse(length(hom_d) > 0, max(hom_d, na.rm = TRUE), NA)
}

# B.2.4 Paired Wilcoxon tests (effect size as absolute rank-biserial)
run_wilcoxon_paired <- function(x, y, label) {
  valid <- !is.na(x) & !is.na(y)
  x_val <- x[valid]; y_val <- y[valid]
  test  <- wilcox.test(x_val, y_val, paired = TRUE, exact = FALSE)
  n     <- sum(valid)
  r     <- abs(1 - (2 * as.numeric(test$statistic)) / (as.numeric(n) * (as.numeric(n) + 1) / 2))
  cat("\n ---", label, "---\n")
  cat("  FANTASIA median:", round(median(x_val), 2),
      "| Blast2GO median:", round(median(y_val), 2), "\n")
  cat("  FANTASIA mean:  ", round(mean(x_val), 2),
      "| Blast2GO mean:  ", round(mean(y_val), 2), "\n")
  cat("  Wilcoxon p =", test$p.value, "| |r| =", round(r, 4), "\n")
  invisible(list(test = test, r = r, n = n,
                 mean_fan = mean(x_val), mean_hom = mean(y_val),
                 sd_fan = sd(x_val), sd_hom = sd(y_val)))
}

res_corpus <- run_wilcoxon_paired(max_ic_corpus$max_IC_fan, max_ic_corpus$max_IC_hom, "Corpus IC")
res_topo   <- run_wilcoxon_paired(max_ic_topo$max_IC_fan,   max_ic_topo$max_IC_hom,   "Topological IC")
res_depth  <- run_wilcoxon_paired(max_depth$max_depth_fan,  max_depth$max_depth_hom,  "GO term depth")

# B.2.5 Spearman correlation between IC methods
common_ic_terms <- intersect(names(ic_corpus), names(ic))
common_ic_terms <- common_ic_terms[is.finite(ic_corpus[common_ic_terms]) &
                                     is.finite(as.numeric(ic[common_ic_terms]))]
spearman_cor <- cor(as.numeric(ic_corpus[common_ic_terms]),
                    as.numeric(ic[common_ic_terms]), method = "spearman")
cat("\nSpearman correlation between IC methods: rho =", round(spearman_cor, 4), "\n")

# ============================================================
# B.3  Effect-size confidence intervals and sensitivity analysis
#      (Reviewer 1, comment 5)
# ============================================================
cat("\n=== B.3 Effect-size CIs and sensitivity analysis ===\n")
set.seed(123)

# --- Rank-biserial statistics ---
rb_paired_stat <- function(x, y) {
  valid <- !is.na(x) & !is.na(y)
  x <- x[valid]; y <- y[valid]
  test <- suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE))
  n <- length(x)
  abs(1 - (2 * as.numeric(test$statistic)) / (n * (n + 1) / 2))
}
rb_unpaired_fast <- function(x, y, n_sample = 2e5) {
  xs <- sample(x, n_sample, replace = TRUE)
  ys <- sample(y, n_sample, replace = TRUE)
  mean(xs > ys) - mean(xs < ys)
}

# --- 5a. Bootstrap 95% CI ---
boot_rb_paired <- function(x, y, R = 1000) {
  valid <- !is.na(x) & !is.na(y)
  x <- x[valid]; y <- y[valid]; n <- length(x)
  b  <- boot(seq_len(n), function(d, i) rb_paired_stat(x[i], y[i]), R = R)
  ci <- boot.ci(b, type = "perc")$percent[4:5]
  c(r = round(b$t0, 3), lower = round(ci[1], 3), upper = round(ci[2], 3))
}

cat("Bootstrap 95% CI (richness, unpaired):\n")
r0_rich <- rb_unpaired_fast(fantasia_go, homologia_go, n_sample = 1e6)
rs_rich <- replicate(1000, {
  xi <- sample(fantasia_go,  length(fantasia_go),  replace = TRUE)
  yi <- sample(homologia_go, length(homologia_go), replace = TRUE)
  rb_unpaired_fast(xi, yi, n_sample = 2e5)
})
ci_rich <- c(r = round(r0_rich, 3),
             lower = round(quantile(rs_rich, 0.025), 3),
             upper = round(quantile(rs_rich, 0.975), 3))
print(ci_rich)

cat("Bootstrap 95% CI (corpus IC, paired):\n");        ci_corpus <- boot_rb_paired(max_ic_corpus$max_IC_fan, max_ic_corpus$max_IC_hom); print(ci_corpus)
cat("Bootstrap 95% CI (topological IC, paired):\n");   ci_topo   <- boot_rb_paired(max_ic_topo$max_IC_fan,   max_ic_topo$max_IC_hom);   print(ci_topo)
cat("Bootstrap 95% CI (GO depth, paired):\n");         ci_depth  <- boot_rb_paired(max_depth$max_depth_fan,  max_depth$max_depth_hom);  print(ci_depth)

ci_table <- rbind(Richness       = ci_rich[c("r","lower","upper")],
                  Corpus_IC      = ci_corpus[c("r","lower","upper")],
                  Topological_IC = ci_topo[c("r","lower","upper")],
                  GO_depth       = ci_depth[c("r","lower","upper")])
write.table(ci_table, "results/effect_size_CI.tsv",
            sep = "\t", quote = FALSE, col.names = NA)

# --- 5b. Sensitivity analysis by subsampling (1000 subsamples of 5000 proteins) ---
n_sub <- 5000; n_rep <- 1000
subsample_paired <- function(fan_vec, hom_vec) {
  N <- sum(!is.na(fan_vec) & !is.na(hom_vec))
  idx_all <- which(!is.na(fan_vec) & !is.na(hom_vec))
  rs <- replicate(n_rep, {
    idx <- sample(idx_all, min(n_sub, N))
    rb_paired_stat(fan_vec[idx], hom_vec[idx])
  })
  c(mean_r = round(mean(rs), 3), sd_r = round(sd(rs), 3))
}
rs_rich_sub <- replicate(n_rep, {
  xi <- sample(fantasia_go,  n_sub)
  yi <- sample(homologia_go, n_sub)
  rb_unpaired_fast(xi, yi, n_sample = 1e5)
})
sub_rich   <- c(mean_r = round(mean(rs_rich_sub), 3), sd_r = round(sd(rs_rich_sub), 3))
sub_corpus <- subsample_paired(max_ic_corpus$max_IC_fan, max_ic_corpus$max_IC_hom)
sub_topo   <- subsample_paired(max_ic_topo$max_IC_fan,   max_ic_topo$max_IC_hom)
sub_depth  <- subsample_paired(max_depth$max_depth_fan,  max_depth$max_depth_hom)

sens_table <- rbind(Richness = sub_rich, Corpus_IC = sub_corpus,
                    Topological_IC = sub_topo, GO_depth = sub_depth)
cat("\nSensitivity analysis (subsampling n =", n_sub, ",", n_rep, "reps):\n")
print(sens_table)
write.table(sens_table, "results/sensitivity_subsampling.tsv",
            sep = "\t", quote = FALSE, col.names = NA)

# ============================================================
# Figure 2  (A: semantic similarity observed vs null; B: specificity)
# ============================================================
cat("\n=== Building Figure 2 ===\n")

# --- Panel A: semantic similarity (Observed vs Null) ---
build_long <- function(df, type_label) {
  data.frame(
    Domain     = rep(c("Biological Process","Molecular Function","Cellular Component"),
                     each = nrow(df)),
    Similarity = c(df$sim_BP, df$sim_MF, df$sim_CC),
    Type       = type_label
  )
}
real_long <- build_long(sim_results, "Observed")
null_long <- build_long(sim_null,    "Null model")
sim_comb  <- rbind(real_long, null_long)
sim_comb  <- sim_comb[!is.na(sim_comb$Similarity) & sim_comb$Similarity > 0, ]
sim_comb$Domain <- factor(sim_comb$Domain,
                          levels = c("Biological Process","Molecular Function","Cellular Component"))
sim_comb$Type   <- factor(sim_comb$Type, levels = c("Observed","Null model"))
medians_comb    <- aggregate(Similarity ~ Domain + Type, data = sim_comb, FUN = median)

panelA <- ggplot(sim_comb, aes(x = Domain, y = Similarity, fill = Type)) +
  geom_violin(trim = TRUE, alpha = 0.6, color = NA, scale = "width", width = 0.85,
              position = position_dodge(width = 0.9)) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.9, color = "grey30",
               position = position_dodge(width = 0.9)) +
  geom_text(data = medians_comb,
            aes(label = paste0("Md=", round(Similarity, 2)), group = Type),
            position = position_dodge(width = 0.9), vjust = -0.6,
            size = 2.2, fontface = "bold", color = "grey20") +
  scale_fill_manual(values = c("Observed" = "#4BACC6", "Null model" = "#BFBFBF")) +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, 12, 2),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Semantic similarity (BMA/Resnik IC)", fill = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = "top", legend.text = element_text(size = 8),
        axis.title = element_text(size = 8.5),
        axis.text = element_text(color = "grey20", size = 7),
        axis.text.x = element_text(angle = 20, hjust = 1))

# --- Panel B: specificity (mean max IC / depth, with |r|) ---
n_paired     <- res_corpus$n
barplot_data <- data.frame(
  Method = factor(rep(c("FANTASIA", "Blast2GO"), 3), levels = c("FANTASIA", "Blast2GO")),
  Metric = factor(rep(c("Corpus IC", "Topological IC", "Depth"), each = 2),
                  levels = c("Corpus IC", "Topological IC", "Depth")),
  Value  = c(res_corpus$mean_fan, res_corpus$mean_hom,
             res_topo$mean_fan,   res_topo$mean_hom,
             res_depth$mean_fan,  res_depth$mean_hom),
  SE     = c(res_corpus$sd_fan / sqrt(n_paired), res_corpus$sd_hom / sqrt(n_paired),
             res_topo$sd_fan   / sqrt(n_paired), res_topo$sd_hom   / sqrt(n_paired),
             res_depth$sd_fan  / sqrt(n_paired), res_depth$sd_hom  / sqrt(n_paired))
)
r_labels <- data.frame(
  Metric = factor(c("Corpus IC", "Topological IC", "Depth"),
                  levels = c("Corpus IC", "Topological IC", "Depth")),
  r_text = c(paste0("|r| = ", round(res_corpus$r, 2)),
             paste0("|r| = ", round(res_topo$r, 2)),
             paste0("|r| = ", round(res_depth$r, 2)))
)
r_labels$y <- sapply(c("Corpus IC", "Topological IC", "Depth"), function(m)
  max(barplot_data$Value[barplot_data$Metric == m]) + 1.8)

panelB <- ggplot(barplot_data, aes(x = Metric, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = Value - SE, ymax = Value + SE),
                position = position_dodge(width = 0.7), width = 0.2, color = "grey30") +
  annotate("text", x = 1, y = max(barplot_data$Value[barplot_data$Metric == "Corpus IC"]) + 0.8,      label = "***", size = 3.5, color = "grey20") +
  annotate("text", x = 2, y = max(barplot_data$Value[barplot_data$Metric == "Topological IC"]) + 0.8, label = "***", size = 3.5, color = "grey20") +
  annotate("text", x = 3, y = max(barplot_data$Value[barplot_data$Metric == "Depth"]) + 0.8,          label = "***", size = 3.5, color = "grey20") +
  geom_text(data = r_labels, aes(x = Metric, y = y, label = r_text), inherit.aes = FALSE,
            size = 2.5, fontface = "italic", color = "grey20") +
  scale_fill_manual(values = c("FANTASIA" = "#4BACC6", "Blast2GO" = "#F79646")) +
  labs(x = NULL, y = "Mean max IC / depth (SE)", fill = "Method") +
  theme_classic(base_size = 9) +
  theme(legend.position = "top", legend.text = element_text(size = 8),
        axis.title = element_text(size = 8.5),
        axis.text = element_text(color = "grey20", size = 7))

fig2 <- panelA + panelB + plot_annotation(tag_levels = "A")
ggsave("results/Figure_2_consistency_specificity.jpg", fig2,
       width = 18, height = 10, units = "cm", dpi = 300, bg = "white")
cat("Figure 2 saved.\n")

# ============================================================
# Export GO terms for taxonomic classification (Python step)
# ============================================================
cat("\n=== Exporting GO terms for Python taxonomy step ===\n")
todos_go_export <- data.frame(
  GO_term = unique(c(
    trimws(unlist(strsplit(fantasia$GOterms, ","))),
    homologia$Annotation.GO.ID
  ))
)
todos_go_export <- todos_go_export[todos_go_export$GO_term != "", , drop = FALSE]
write.table(todos_go_export, "results/go_terms_to_classify.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Exported", nrow(todos_go_export), "unique GO terms\n")

cat("\n============================================================\n")
cat("1_Comparative_analysis.R complete.\n")
cat("Next step: run 2_taxonomic_classification.py\n")
cat("  Input:  results/go_terms_to_classify.tsv\n")
cat("  Output: results/taxon_map_v2.tsv\n")
cat("============================================================\n")
