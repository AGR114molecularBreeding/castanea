#########################################################################################################
# 1_ Comparative analysis
# Comparative analysis of FANTASIA vs Blast2Go annotations in C.sativa
#
# Produces:Figures 1, 2A, 2B, 3 and all quantitative reslts reported in the manuscript 
#
# Run before the python taxonomic classification step, and at the end this script will exports
# go_terms_to_classify.tsv which is the imput for the script "02_taxonomic_classification.py"
#
# Postigo-Luque et al. 
#########################################################################################################
setwd("F:/Tesis/Mis_artículos/Artículo_FANTASIA/Github")
source("0_Load_data.R")
library(eulerr)
library(ontologyIndex)
library(ontologySimilarity)



# Output directory
dir.create("results", showWarnings = FALSE)

# ============================================================
# A.1  Annotation coverage
# ============================================================
cat("\n=== A.1 Annotation coverage ===\n")
cat("Total proteins:       ", n_total, "\n")
cat("FANTASIA annotated:   ", n_fantasia, "(", round(n_fantasia/n_total*100, 1), "%)\n")
cat("Blast2Go annotated:   ", n_homologia, "(", round(n_homologia/n_total*100, 1), "%)\n")
cat("Both methods:         ", n_both, "\n")
cat("Only FANTASIA:        ", length(only_fantasia), "\n")
cat("Only Blast2Go:        ", length(only_homologia), "\n")

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
    Domain                  = label,
    FANTASIA_unique         = n_fan,
    Blast2Go_unique         = n_hom,
    Shared                  = n_comun,
    Exclusive_FANTASIA      = length(setdiff(set_fan, set_hom)),
    Exclusive_Blast2Go      = length(setdiff(set_hom, set_fan)),
    Jaccard                 = round(n_comun / n_union, 4),
    Dice                    = round(2 * n_comun / (n_fan + n_hom), 4),
    Overlap                 = round(n_comun / min(n_fan, n_hom), 4)
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
go_counts_homologia$Method <- "Blast2Go"

go_counts        <- rbind(go_counts_fantasia[, c("ProteinID","n_GO","Method")],
                          go_counts_homologia[, c("ProteinID","n_GO","Method")])
go_counts$n_GO   <- as.numeric(go_counts$n_GO)
go_counts$Method <- factor(go_counts$Method, levels = c("FANTASIA", "Blast2Go"))

# Mann-Whitney U test
fantasia_go  <- go_counts$n_GO[go_counts$Method == "FANTASIA"]
homologia_go <- go_counts$n_GO[go_counts$Method == "Blast2Go"]
test_mw      <- wilcox.test(fantasia_go, homologia_go, exact = FALSE)
n1 <- length(fantasia_go)
n2 <- length(homologia_go)
r_biserial <- 1 - (2 * as.numeric(test_mw$statistic)) / (as.numeric(n1) * as.numeric(n2))

cat("  FANTASIA median:", median(fantasia_go), "\n")
cat("  Blast2Go median:", median(homologia_go), "\n")
cat("  Mann-Whitney p =", test_mw$p.value, "\n")
cat("  Rank-biserial r =", round(r_biserial, 4), "\n")

# Figure 1: Violin plot
medians <- aggregate(n_GO ~ Method, data = go_counts, FUN = median)

p_fig1 <- ggplot(go_counts, aes(x = Method, y = n_GO, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA, coef = Inf, color = "grey30") +
  geom_text(data = medians,
            aes(x = Method, y = n_GO, label = paste0("Md=", n_GO)),
            vjust = -0.8, size = 3.2, fontface = "bold", color = "grey20") +
  scale_fill_manual(values = c("FANTASIA" = "#4BACC6", "Blast2Go" = "#F79646")) +
  scale_y_log10(breaks = c(1, 2, 5, 10, 25, 50, 100, 200)) +
  annotate("text", x = 1.5, y = 100,
           label = paste0("Mann-Whitney U\np < 2.2e-16\nr = ",
                          round(r_biserial, 3)),
           size = 3, hjust = 0.5, color = "grey20") +
  labs(x = NULL, y = "Number of GO terms per protein (log scale)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(color = "grey20"))
ggsave("results/Figure_1_violin_GO_per_protein.png", p_fig1,
       width = 6, height = 5, dpi = 300, bg = "white")
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

cat("  Common proteins:                  ", length(proteinas_comunes), "\n")
cat("  Sharing at least 1 GO term:       ", n_acuerdo, "\n")
cat("  Percentage:                        ", round(perc_acuerdo, 1), "%\n")

# ============================================================
# B.1.2  Semantic similarity BMA/Resnik  (Figure 2A)
# ============================================================
cat("\n=== B.1.2 Semantic similarity BMA/Resnik ===\n")

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

cat("Calculating BMA/Resnik for", length(proteinas_comunes),
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
    cat("  Progress:", i, "/", length(proteinas_comunes),
        "| Elapsed:", round(elapsed, 1), "min",
        "| Remaining:", round(remaining, 1), "min\n")
  }
}
cat("Calculation complete.\n")

for (dom in c("BP", "MF", "CC")) {
  vals <- sim_results[[paste0("sim_", dom)]]
  vals <- vals[!is.na(vals)]
  cat("  ", dom, "- Median:", round(median(vals), 2),
      "| Mean:", round(mean(vals), 2), "\n")
}

write.table(sim_results, "results/B1_semantic_similarity_BMA_Resnik.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Figure 2A: Semantic similarity violin
sim_long <- data.frame(
  Domain     = rep(c("Biological Process","Molecular Function","Cellular Component"),
                   each = nrow(sim_results)),
  Similarity = c(sim_results$sim_BP, sim_results$sim_MF, sim_results$sim_CC)
)
sim_long <- sim_long[!is.na(sim_long$Similarity) & sim_long$Similarity > 0, ]
sim_long$Domain <- factor(sim_long$Domain,
                          levels = c("Biological Process",
                                     "Molecular Function",
                                     "Cellular Component"))

medians_sim <- aggregate(Similarity ~ Domain, data = sim_long, FUN = median)

p_fig2a <- ggplot(sim_long, aes(x = Domain, y = Similarity, fill = Domain)) +
  geom_violin(trim = TRUE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA,
               color = "grey30", fill = "white", alpha = 0.8) +
  geom_text(data = medians_sim,
            aes(x = Domain, y = Similarity,
                label = paste0("Md=", round(Similarity, 2))),
            vjust = -0.8, size = 3, fontface = "bold", color = "grey20") +
  scale_fill_manual(values = c("Biological Process" = "#4BACC6",
                               "Molecular Function" = "#F79646",
                               "Cellular Component" = "#9DC3E6")) +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, 12, 2),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(x = "GO domain", y = "Semantic similarity (BMA/Resnik IC)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", axis.text = element_text(color = "grey20"))
ggsave("results/Figure_2A_semantic_similarity.png", p_fig2a,
       width = 7, height = 5, dpi = 300, bg = "white")
cat("Figure 2A saved.\n")

# ============================================================
# B.2  Annotation specificity  (Figure 2B)
# ============================================================
cat("\n=== B.2 Annotation specificity ===\n")

# B.2.1 Corpus-based IC
corpus_fantasia  <- trimws(unlist(strsplit(fantasia$GOterms, ",")))
corpus_homologia <- homologia$Annotation.GO.ID
corpus_total     <- c(corpus_fantasia, corpus_homologia)
term_counts      <- table(corpus_total)
ic_corpus        <- -log2(as.numeric(term_counts) / length(corpus_total))
names(ic_corpus) <- names(term_counts)

# B.2.3 GO term depth
all_go_terms <- union(
  unique(trimws(unlist(strsplit(fantasia$GOterms, ",")))),
  unique(homologia$Annotation.GO.ID)
)
depth_map_full <- sapply(all_go_terms, function(go_id) {
  tryCatch({ length(get_ancestors(go_obo, go_id)) }, error = function(e) NA)
})
names(depth_map_full) <- all_go_terms

# B.2.4 Max IC and depth per protein (common proteins)
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

# B.2.5 Paired Wilcoxon tests
run_wilcoxon_paired <- function(x, y, label) {
  valid <- !is.na(x) & !is.na(y)
  x_val <- x[valid]; y_val <- y[valid]
  test  <- wilcox.test(x_val, y_val, paired = TRUE, exact = FALSE)
  n     <- sum(valid)
  r     <- 1 - (2 * as.numeric(test$statistic)) / (as.numeric(n) * (as.numeric(n) + 1) / 2)
  cat("\n ---", label, "---\n")
  cat("  FANTASIA median:", round(median(x_val), 2),
      "| Blast2Go median:", round(median(y_val), 2), "\n")
  cat("  FANTASIA mean:  ", round(mean(x_val), 2),
      "| Blast2Go mean:  ", round(mean(y_val), 2), "\n")
  cat("  Wilcoxon p =", test$p.value, "| r =", round(r, 4), "\n")
  invisible(list(test = test, r = r, n = n,
                 mean_fan = mean(x_val), mean_hom = mean(y_val),
                 sd_fan = sd(x_val), sd_hom = sd(y_val)))
}

res_corpus <- run_wilcoxon_paired(max_ic_corpus$max_IC_fan, max_ic_corpus$max_IC_hom,
                                  "Corpus IC")
res_topo   <- run_wilcoxon_paired(max_ic_topo$max_IC_fan, max_ic_topo$max_IC_hom,
                                  "Topological IC")
res_depth  <- run_wilcoxon_paired(max_depth$max_depth_fan, max_depth$max_depth_hom,
                                  "GO term depth")

# B.2.10 Spearman correlation between IC methods
common_ic_terms <- intersect(names(ic_corpus), names(ic))
common_ic_terms <- common_ic_terms[is.finite(ic_corpus[common_ic_terms]) &
                                     is.finite(as.numeric(ic[common_ic_terms]))]
spearman_cor <- cor(as.numeric(ic_corpus[common_ic_terms]),
                    as.numeric(ic[common_ic_terms]), method = "spearman")
cat("\nSpearman correlation between IC methods: rho =", round(spearman_cor, 4), "\n")

# Figure 2B: Barplot of mean max IC / depth
n_paired     <- res_corpus$n
barplot_data <- data.frame(
  Method = rep(c("FANTASIA", "Blast2Go"), 3),
  Metric = rep(c("Corpus IC", "Topological IC", "Depth"), each = 2),
  Value  = c(res_corpus$mean_fan, res_corpus$mean_hom,
             res_topo$mean_fan,   res_topo$mean_hom,
             res_depth$mean_fan,  res_depth$mean_hom),
  SE     = c(res_corpus$sd_fan / sqrt(n_paired), res_corpus$sd_hom / sqrt(n_paired),
             res_topo$sd_fan   / sqrt(n_paired), res_topo$sd_hom   / sqrt(n_paired),
             res_depth$sd_fan  / sqrt(n_paired), res_depth$sd_hom  / sqrt(n_paired))
)
barplot_data$Method <- factor(barplot_data$Method, levels = c("FANTASIA", "Blast2Go"))
barplot_data$Metric <- factor(barplot_data$Metric,
                              levels = c("Corpus IC", "Topological IC", "Depth"))

p_fig2b <- ggplot(barplot_data, aes(x = Metric, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = Value - SE, ymax = Value + SE),
                position = position_dodge(width = 0.7),
                width = 0.2, color = "grey30") +
  annotate("text", x = 1, y = max(barplot_data$Value[barplot_data$Metric == "Corpus IC"]) + 0.8,
           label = "***", size = 5, color = "grey20") +
  annotate("text", x = 2, y = max(barplot_data$Value[barplot_data$Metric == "Topological IC"]) + 0.8,
           label = "***", size = 5, color = "grey20") +
  annotate("text", x = 3, y = max(barplot_data$Value[barplot_data$Metric == "Depth"]) + 0.8,
           label = "***", size = 5, color = "grey20") +
  scale_fill_manual(values = c("FANTASIA" = "#4BACC6", "Blast2Go" = "#F79646")) +
  labs(x = NULL, y = "Mean max IC / depth (SE)",
       fill = "Method") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right", axis.text = element_text(color = "grey20"))
ggsave("results/Figure_2B_specificity_barplot.png", p_fig2b,
       width = 8, height = 5, dpi = 300, bg = "white")
cat("Figure 2B saved.\n")

# ============================================================
# B.4  Functional enrichment analysis  (Figure 3)
# ============================================================
cat("\n=== B.4 Functional enrichment analysis ===\n")

pop             <- union(fantasia_annotated, homologia_annotated)
study_fantasia  <- only_fantasia
study_homologia <- only_homologia
assoc_fantasia  <- fantasia_go_list[names(fantasia_go_list) %in% pop]
assoc_homologia <- homologia_go_list[names(homologia_go_list) %in% pop]

run_goea_fisher <- function(study_prots, pop_prots, assoc, fdr_threshold = 0.05) {
  all_terms <- unique(unlist(assoc))
  if (length(all_terms) == 0) return(data.frame())
  n_pop   <- length(pop_prots)
  n_study <- length(study_prots)
  results <- lapply(all_terms, function(go_term) {
    pop_with     <- sum(sapply(assoc, function(gos) go_term %in% gos))
    study_with   <- sum(sapply(assoc[names(assoc) %in% study_prots],
                               function(gos) go_term %in% gos))
    if (study_with == 0) return(NULL)
    mat <- matrix(c(study_with, n_study - study_with,
                    pop_with - study_with,
                    n_pop - pop_with - n_study + study_with), nrow = 2)
    ft <- fisher.test(mat, alternative = "greater")
    data.frame(
      GO_ID       = go_term,
      GO_name     = unname(tryCatch(Term(go_term), error = function(e) go_term)),
      Ontology    = unname(tryCatch(Ontology(go_term), error = function(e) NA_character_)),
      study_count = study_with,
      study_total = n_study,
      pop_count   = pop_with,
      pop_total   = n_pop,
      p_value     = ft$p.value,
      stringsAsFactors = FALSE
    )
  })
  results_df <- do.call(rbind, Filter(Negate(is.null), results))
  if (nrow(results_df) == 0) return(data.frame())
  results_df$p_adjusted <- p.adjust(results_df$p_value, method = "BH")
  results_df <- results_df[results_df$p_adjusted < fdr_threshold, ]
  results_df[order(results_df$p_adjusted), ]
}

cat("Running GOEA for FANTASIA-exclusive proteins...\n")
goea_fantasia <- run_goea_fisher(study_fantasia, pop, assoc_fantasia)
cat("  Significant terms:", nrow(goea_fantasia), "\n")

cat("Running GOEA for Blast2Go-exclusive proteins...\n")
goea_homologia <- run_goea_fisher(study_homologia, pop, assoc_homologia)
cat("  Significant terms:", nrow(goea_homologia), "\n")

write.table(goea_fantasia, "results/B4_GOEA_FANTASIA_exclusive.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(goea_homologia, "results/B4_GOEA_Blast2Go_exclusive.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Figure 3: Top 15 enriched BP terms per method
plot_goea_barplot <- function(goea_df, method_label, fill_color, top_n = 15) {
  bp <- goea_df[!is.na(goea_df$Ontology) & goea_df$Ontology == "BP", ]
  if (nrow(bp) == 0) return(NULL)
  bp <- head(bp, top_n)
  bp$GO_short <- ifelse(nchar(bp$GO_name) > 50,
                        paste0(substr(bp$GO_name, 1, 47), "..."),
                        bp$GO_name)
  bp$GO_short <- factor(bp$GO_short, levels = rev(bp$GO_short))
  ggplot(bp, aes(x = -log10(p_adjusted), y = GO_short)) +
    geom_bar(stat = "identity", fill = fill_color, alpha = 0.85) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    labs(x = expression(-log[10](p[adjusted])), y = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_text(color = "grey20", size = 9),
          axis.text.x = element_text(color = "grey20"))
}

p_fig3a <- plot_goea_barplot(goea_fantasia,  "FANTASIA", "#4BACC6")
p_fig3b <- plot_goea_barplot(goea_homologia, "Blast2Go", "#F79646")

if (!is.null(p_fig3a)) {
  ggsave("results/Figure_3A_GOEA_FANTASIA.png", p_fig3a,
         width = 9, height = 6, dpi = 300, bg = "white")
  cat("Figure 3A saved.\n")
}
if (!is.null(p_fig3b)) {
  ggsave("results/Figure_3B_GOEA_Blast2Go.png", p_fig3b,
         width = 9, height = 6, dpi = 300, bg = "white")
  cat("Figure 3B saved.\n")
}

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
cat("01_comparative_analysis.R complete.\n")
cat("Next step: run 02_taxonomic_classification.py\n")
cat("  Input:  results/go_terms_to_classify.tsv\n")
cat("  Output: results/taxon_map_v2.tsv\n")
cat("============================================================\n")