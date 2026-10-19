#########################################################################################################
# 6_AHRD_comparison.R
# Supplementary comparison: FANTASIA vs AHRD in C. sativa (Reviewer 1, comment 1)
#
# Parallel analysis to the main FANTASIA vs Blast2GO comparison, using AHRD as a
# second, independent homology-based method. Self-contained: it loads its own data.
#
# Required input files (place in ./data/):
#   data/Casat.proteins.funct_ahrd.tsv          AHRD annotation
#   data/Csath1_FANTASIA_TopGO.txt              FANTASIA annotation
#   data/Cast.1_0.hap1.clean.fa                 C. sativa hap1 proteome (filtered)
#   data/Cast.1_0.hap1.clean_cdhit100.pep.clstr CD-HIT 100%-identity clusters
#   data/go.obo                                 Gene Ontology (go.obo)
#
# Produces: Figure S2 (Venn: proteins and GO vocabulary), Figure S3 (richness),
#           Figure S4 (semantic similarity vs null + specificity), and result tables.
#
# Postigo-Luque et al.
#########################################################################################################

# Set the working directory to the repository folder before running, e.g.:
# setwd("/path/to/Functional annotation comparative analysis FANTASIA vs Blast2Go")

library(Biostrings)
library(ggplot2)
library(GO.db)
library(ontologyIndex)
library(ontologySimilarity)
library(eulerr)
library(boot)
library(patchwork)

dir.create("results", showWarnings = FALSE)

#### File paths ####
AHRD_FILE     <- "data/Casat.proteins.funct_ahrd.tsv"
FANTASIA_FILE <- "data/Csath1_FANTASIA_TopGO.txt"
FASTA_FILE    <- "data/Cast.1_0.hap1.clean.fa"
CLSTR_FILE    <- "data/Cast.1_0.hap1.clean_cdhit100.pep.clstr"
GO_OBO_FILE   <- "data/go.obo"

# ============================================================
# Load annotations (FANTASIA + AHRD)
# ============================================================
ahrd_raw <- read.table(AHRD_FILE, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                       quote = "", fill = TRUE, comment.char = "#")
# AHRD GO column comes as "GO:0000166, GO:0004672, ..." (treated like Blast2GO)
ahrd_raw$GO_clean <- gsub(",\\s*$", "", trimws(ahrd_raw$Gene.Ontology.Term))
ahrd <- ahrd_raw[ahrd_raw$GO_clean != "" & !is.na(ahrd_raw$GO_clean), ]

fantasia <- read.table(FANTASIA_FILE, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(fantasia) <- c("ProteinID", "GOterms")

#### Harmonise protein IDs ####
ahrd$ProteinID <- ahrd$Protein.Accession
ahrd$ProteinID <- gsub("-mRNA-\\d+$", "", ahrd$ProteinID)
ahrd$ProteinID <- gsub("\\.t\\d+$",   "", ahrd$ProteinID)
ahrd$ProteinID <- gsub("\\..*$",      "", ahrd$ProteinID)
fantasia$ProteinID <- gsub("\\..*$",  "", fantasia$ProteinID)

#### Propagate FANTASIA annotations to identical sequences (CD-HIT) ####
clstr <- readLines(CLSTR_FILE)
.cluster_id <- NA
.clusters   <- list()
for (line in clstr) {
  if (grepl("^>Cluster", line)) {
    .cluster_id <- line
  } else {
    id <- gsub(".*>(CSH1-[^.]+).*", "\\1", line)
    .clusters[[.cluster_id]] <- c(.clusters[[.cluster_id]], id)
  }
}
.multi <- .clusters[sapply(.clusters, length) > 1]
.fantasia_ids <- unique(fantasia$ProteinID)
.rescatable   <- data.frame(removed = character(), representative = character(),
                            stringsAsFactors = FALSE)
for (cl in .multi) {
  en_fan <- cl[cl %in% .fantasia_ids]
  no_fan <- cl[!cl %in% .fantasia_ids]
  if (length(en_fan) > 0 & length(no_fan) > 0) {
    for (nf in no_fan) {
      .rescatable <- rbind(.rescatable, data.frame(removed = nf, representative = en_fan[1]))
    }
  }
}
.prop_rows <- data.frame(
  ProteinID = .rescatable$removed,
  GOterms   = fantasia$GOterms[match(.rescatable$representative, fantasia$ProteinID)],
  stringsAsFactors = FALSE)
.prop_rows <- .prop_rows[!is.na(.prop_rows$GOterms), ]
fantasia   <- rbind(fantasia, .prop_rows)
cat("Propagated FANTASIA annotations to", nrow(.prop_rows), "identical sequences\n")

#### Define annotation sets ####
fasta_seqs <- readAAStringSet(FASTA_FILE)
n_total    <- length(names(fasta_seqs))

ahrd_annotated     <- unique(ahrd$ProteinID)
fantasia_annotated <- unique(fantasia$ProteinID)
both_annotated <- intersect(fantasia_annotated, ahrd_annotated)
only_fantasia  <- setdiff(fantasia_annotated, ahrd_annotated)
only_ahrd      <- setdiff(ahrd_annotated, fantasia_annotated)
n_fantasia <- length(fantasia_annotated); n_ahrd <- length(ahrd_annotated); n_both <- length(both_annotated)

#### GO term lists per protein ####
fantasia_go_list <- strsplit(fantasia$GOterms, ",")
fantasia_go_list <- lapply(fantasia_go_list, function(x) unique(trimws(x)))
names(fantasia_go_list) <- fantasia$ProteinID
ahrd_go_list <- strsplit(ahrd$GO_clean, ",")
ahrd_go_list <- lapply(ahrd_go_list, function(x) unique(trimws(x)))
names(ahrd_go_list) <- ahrd$ProteinID

proteinas_comunes <- intersect(names(fantasia_go_list), names(ahrd_go_list))
fantasia_comun    <- fantasia_go_list[proteinas_comunes]
ahrd_comun        <- ahrd_go_list[proteinas_comunes]

#### GO domain mapping ####
fantasia_go_largo <- unique(trimws(unlist(strsplit(fantasia$GOterms, ","))))
ahrd_go_largo     <- unique(trimws(unlist(strsplit(ahrd$GO_clean, ","))))
todos_go          <- union(fantasia_go_largo, ahrd_go_largo)
dominio_df  <- AnnotationDbi::select(GO.db, keys = todos_go, columns = "ONTOLOGY", keytype = "GOID")
dominio_map <- setNames(dominio_df$ONTOLOGY, dominio_df$GOID)

cat("Data loaded:", n_total, "proteins |", n_fantasia, "FANTASIA |",
    n_ahrd, "AHRD |", n_both, "both\n")

# ============================================================
# A.1 Coverage
# ============================================================
cat("\n=== A.1 Coverage ===\n")
cat("FANTASIA:", n_fantasia, "(", round(n_fantasia/n_total*100, 2), "%) | ",
    "AHRD:", n_ahrd, "(", round(n_ahrd/n_total*100, 2), "%) | Both:", n_both, "\n")
cat("Only FANTASIA:", length(only_fantasia), "| Only AHRD:", length(only_ahrd), "\n")

# ============================================================
# A.2 GO vocabulary diversity
# ============================================================
calc_vocab <- function(set_fan, set_ah, label) {
  n_fan <- length(set_fan); n_ah <- length(set_ah)
  n_comun <- length(intersect(set_fan, set_ah)); n_union <- length(union(set_fan, set_ah))
  data.frame(Domain = label, FANTASIA_unique = n_fan, AHRD_unique = n_ah, Shared = n_comun,
             Exclusive_FANTASIA = length(setdiff(set_fan, set_ah)),
             Exclusive_AHRD = length(setdiff(set_ah, set_fan)),
             Jaccard = round(n_comun/n_union, 4), Dice = round(2*n_comun/(n_fan+n_ah), 4),
             Overlap = round(n_comun/min(n_fan, n_ah), 4))
}
fan_bp <- fantasia_go_largo[dominio_map[fantasia_go_largo] == "BP" & !is.na(dominio_map[fantasia_go_largo])]
fan_mf <- fantasia_go_largo[dominio_map[fantasia_go_largo] == "MF" & !is.na(dominio_map[fantasia_go_largo])]
fan_cc <- fantasia_go_largo[dominio_map[fantasia_go_largo] == "CC" & !is.na(dominio_map[fantasia_go_largo])]
ah_bp  <- ahrd_go_largo[dominio_map[ahrd_go_largo] == "BP" & !is.na(dominio_map[ahrd_go_largo])]
ah_mf  <- ahrd_go_largo[dominio_map[ahrd_go_largo] == "MF" & !is.na(dominio_map[ahrd_go_largo])]
ah_cc  <- ahrd_go_largo[dominio_map[ahrd_go_largo] == "CC" & !is.na(dominio_map[ahrd_go_largo])]
vocab_table <- rbind(
  calc_vocab(fantasia_go_largo, ahrd_go_largo, "Total"),
  calc_vocab(fan_bp, ah_bp, "BP"), calc_vocab(fan_mf, ah_mf, "MF"), calc_vocab(fan_cc, ah_cc, "CC"))
print(vocab_table)
write.table(vocab_table, "results/A2_vocabulary_comparison_AHRD.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# ------------------------------------------------------------
# Supplementary Figure S2: Venn (proteins + GO vocabulary)
# ------------------------------------------------------------
png("results/Figure_S2A_venn_proteins_AHRD.png", width = 9, height = 9, units = "cm", res = 300)
print(plot(euler(c("FANTASIA" = length(only_fantasia), "AHRD" = length(only_ahrd),
                   "FANTASIA&AHRD" = n_both)),
           quantities = TRUE, fills = c("#4BACC6", "#F79646"), main = "Annotated proteins"))
dev.off()
png("results/Figure_S2B_venn_GO_AHRD.png", width = 9, height = 9, units = "cm", res = 300)
print(plot(euler(c("FANTASIA" = vocab_table$Exclusive_FANTASIA[1], "AHRD" = vocab_table$Exclusive_AHRD[1],
                   "FANTASIA&AHRD" = vocab_table$Shared[1])),
           quantities = TRUE, fills = c("#4BACC6", "#F79646"), main = "GO vocabulary"))
dev.off()
cat("Figure S2 (Venn FANTASIA vs AHRD) saved.\n")

# ============================================================
# B.1 Exact agreement
# ============================================================
acuerdo <- sapply(proteinas_comunes, function(prot)
  length(intersect(fantasia_comun[[prot]], ahrd_comun[[prot]])) > 0)
cat("\nCommon proteins:", length(proteinas_comunes),
    "| Sharing >=1 GO:", sum(acuerdo),
    "(", round(sum(acuerdo)/length(proteinas_comunes)*100, 1), "%)\n")

# ============================================================
# Figure S3: Annotation richness (GO terms per protein)
# ============================================================
go_counts <- rbind(
  data.frame(n_GO = sapply(strsplit(fantasia$GOterms, ","), function(x) length(unique(trimws(x)))), Method = "FANTASIA"),
  data.frame(n_GO = sapply(strsplit(ahrd$GO_clean, ","),   function(x) length(unique(trimws(x)))), Method = "AHRD"))
go_counts$n_GO   <- as.numeric(go_counts$n_GO)
go_counts$Method <- factor(go_counts$Method, levels = c("FANTASIA", "AHRD"))
fan_go <- go_counts$n_GO[go_counts$Method == "FANTASIA"]; ah_go <- go_counts$n_GO[go_counts$Method == "AHRD"]
tmw  <- wilcox.test(fan_go, ah_go, exact = FALSE)
r_b  <- abs(1 - (2*as.numeric(tmw$statistic))/(length(fan_go)*length(ah_go)))
cat("Richness | FAN median:", median(fan_go), "| AHRD median:", median(ah_go),
    "| p =", format(tmw$p.value, digits = 3), "| |r| =", round(r_b, 3), "\n")
medians <- aggregate(n_GO ~ Method, data = go_counts, FUN = median)

p_s3 <- ggplot(go_counts, aes(x = Method, y = n_GO, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA, coef = Inf, color = "grey30") +
  geom_text(data = medians, aes(x = Method, y = n_GO, label = paste0("Md=", n_GO)),
            vjust = -0.8, size = 3, fontface = "bold", color = "grey20") +
  scale_fill_manual(values = c("FANTASIA" = "#4BACC6", "AHRD" = "#F79646")) +
  scale_y_log10(breaks = c(1, 2, 5, 10, 25, 50, 100, 200)) +
  annotate("text", x = 1.5, y = 100,
           label = paste0("Mann-Whitney U\np < 2.2e-16\n|r| = ", round(r_b, 2)),
           size = 2.8, hjust = 0.5, color = "grey20") +
  labs(x = NULL, y = "Number of GO terms per protein (log scale)") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none", axis.title = element_text(size = 9),
        axis.text = element_text(color = "grey20", size = 8))
ggsave("results/Figure_S3_AHRD_richness.jpg", p_s3, width = 9, height = 9, units = "cm", dpi = 300, bg = "white")
cat("Figure S3 (richness) saved.\n")

# ============================================================
# Semantic similarity BMA/Resnik: observed and null model
# ============================================================
cat("\n=== Semantic similarity FANTASIA vs AHRD ===\n")
go_obo <- get_ontology(GO_OBO_FILE, extract_tags = "everything")
ic     <- descendants_IC(go_obo)
namespace_map <- setNames(go_obo$namespace, go_obo$id)

resnik_sim <- function(t1, t2, go_obo, ic) {
  if (is.na(ic[t1]) || is.na(ic[t2])) return(0)
  common <- intersect(c(t1, get_ancestors(go_obo, t1)), c(t2, get_ancestors(go_obo, t2)))
  if (length(common) == 0) return(0)
  max(ic[common], na.rm = TRUE)
}
calc_bma <- function(set1, set2, go_obo, ic) {
  set1 <- set1[set1 %in% names(ic) & !is.na(ic[set1])]
  set2 <- set2[set2 %in% names(ic) & !is.na(ic[set2])]
  if (length(set1) == 0 || length(set2) == 0) return(NA)
  best_1to2 <- sapply(set1, function(t1) max(sapply(set2, function(t2) resnik_sim(t1, t2, go_obo, ic)), na.rm = TRUE))
  best_2to1 <- sapply(set2, function(t2) max(sapply(set1, function(t1) resnik_sim(t1, t2, go_obo, ic)), na.rm = TRUE))
  (mean(best_1to2) + mean(best_2to1)) / 2
}
filter_by_domain <- function(go_terms, domain)
  go_terms[go_terms %in% names(namespace_map) & namespace_map[go_terms] == domain]

# --- Observed ---
sim_results <- data.frame(ProteinID = proteinas_comunes, sim_BP = NA_real_, sim_MF = NA_real_, sim_CC = NA_real_)
cat("Calculating observed similarity for", length(proteinas_comunes), "proteins...\n")
for (i in seq_along(proteinas_comunes)) {
  fan <- fantasia_comun[[proteinas_comunes[i]]]; ah <- ahrd_comun[[proteinas_comunes[i]]]
  sim_results$sim_BP[i] <- calc_bma(filter_by_domain(fan, "biological_process"), filter_by_domain(ah, "biological_process"), go_obo, ic)
  sim_results$sim_MF[i] <- calc_bma(filter_by_domain(fan, "molecular_function"), filter_by_domain(ah, "molecular_function"), go_obo, ic)
  sim_results$sim_CC[i] <- calc_bma(filter_by_domain(fan, "cellular_component"), filter_by_domain(ah, "cellular_component"), go_obo, ic)
  if (i %% 2000 == 0) cat("  Observed:", i, "/", length(proteinas_comunes), "\n")
}
write.table(sim_results, "results/B1_semantic_similarity_AHRD_observed.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# --- Null model ---
set.seed(123)
shuffled <- sample(seq_along(proteinas_comunes))
sim_null <- data.frame(ProteinID = proteinas_comunes, sim_BP = NA_real_, sim_MF = NA_real_, sim_CC = NA_real_)
cat("Calculating null similarity...\n")
for (i in seq_along(proteinas_comunes)) {
  fan <- fantasia_comun[[proteinas_comunes[i]]]; ah <- ahrd_comun[[proteinas_comunes[shuffled[i]]]]
  sim_null$sim_BP[i] <- calc_bma(filter_by_domain(fan, "biological_process"), filter_by_domain(ah, "biological_process"), go_obo, ic)
  sim_null$sim_MF[i] <- calc_bma(filter_by_domain(fan, "molecular_function"), filter_by_domain(ah, "molecular_function"), go_obo, ic)
  sim_null$sim_CC[i] <- calc_bma(filter_by_domain(fan, "cellular_component"), filter_by_domain(ah, "cellular_component"), go_obo, ic)
  if (i %% 2000 == 0) cat("  Null:", i, "/", length(proteinas_comunes), "\n")
}
write.table(sim_null, "results/B1_semantic_similarity_AHRD_null.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Observed vs Null (FANTASIA vs AHRD) ===\n")
for (dom in c("BP","MF","CC")) {
  real <- sim_results[[paste0("sim_",dom)]]; real <- real[!is.na(real)]
  null <- sim_null[[paste0("sim_",dom)]];    null <- null[!is.na(null)]
  wt <- wilcox.test(real, null)
  cat(" ", dom, "| Observed median:", round(median(real),2),
      "| Null median:", round(median(null),2), "| p =", format(wt$p.value, digits = 3), "\n")
}

# ============================================================
# Annotation specificity (+ effect size CIs and sensitivity)
# ============================================================
corpus_total <- c(trimws(unlist(strsplit(fantasia$GOterms, ","))), trimws(unlist(strsplit(ahrd$GO_clean, ","))))
term_counts  <- table(corpus_total)
ic_corpus    <- -log2(as.numeric(term_counts) / length(corpus_total))
names(ic_corpus) <- names(term_counts)
all_go_terms <- union(unique(fantasia_go_largo), unique(ahrd_go_largo))
depth_map_full <- sapply(all_go_terms, function(go_id) tryCatch(length(get_ancestors(go_obo, go_id)), error = function(e) NA))
names(depth_map_full) <- all_go_terms

max_ic_corpus <- data.frame(fan = NA_real_, ah = NA_real_)[rep(1, length(proteinas_comunes)), ]
max_ic_topo   <- max_ic_corpus; max_depth <- max_ic_corpus
for (i in seq_along(proteinas_comunes)) {
  fg <- fantasia_comun[[proteinas_comunes[i]]]; ag <- ahrd_comun[[proteinas_comunes[i]]]
  fc <- ic_corpus[fg[fg %in% names(ic_corpus)]]; ac <- ic_corpus[ag[ag %in% names(ic_corpus)]]
  max_ic_corpus$fan[i] <- ifelse(length(fc) > 0, max(fc, na.rm = TRUE), NA)
  max_ic_corpus$ah[i]  <- ifelse(length(ac) > 0, max(ac, na.rm = TRUE), NA)
  ft <- ic[fg[fg %in% names(ic)]]; at <- ic[ag[ag %in% names(ic)]]
  max_ic_topo$fan[i] <- ifelse(length(ft) > 0, max(ft, na.rm = TRUE), NA)
  max_ic_topo$ah[i]  <- ifelse(length(at) > 0, max(at, na.rm = TRUE), NA)
  fd <- depth_map_full[fg[fg %in% names(depth_map_full)]]; ad <- depth_map_full[ag[ag %in% names(depth_map_full)]]
  max_depth$fan[i] <- ifelse(length(fd) > 0, max(fd, na.rm = TRUE), NA)
  max_depth$ah[i]  <- ifelse(length(ad) > 0, max(ad, na.rm = TRUE), NA)
}

run_wilcox <- function(x, y, label) {
  v <- !is.na(x) & !is.na(y); x <- x[v]; y <- y[v]
  test <- wilcox.test(x, y, paired = TRUE, exact = FALSE); n <- length(x)
  r <- abs(1 - (2*as.numeric(test$statistic))/(n*(n+1)/2))
  cat("\n", label, "| FAN median:", round(median(x),2), "| AHRD median:", round(median(y),2),
      "| p =", format(test$p.value, digits = 3), "| |r| =", round(r, 3), "\n")
  invisible(list(r = r, n = n, mean_fan = mean(x), mean_ah = mean(y), sd_fan = sd(x), sd_ah = sd(y)))
}
res_corpus <- run_wilcox(max_ic_corpus$fan, max_ic_corpus$ah, "Corpus IC")
res_topo   <- run_wilcox(max_ic_topo$fan,  max_ic_topo$ah,  "Topological IC")
res_depth  <- run_wilcox(max_depth$fan,    max_depth$ah,    "GO depth")

# Bootstrap CIs + sensitivity subsampling
set.seed(123)
rb_stat <- function(x, y) {
  test <- suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE)); n <- length(x)
  abs(1 - (2*as.numeric(test$statistic))/(n*(n+1)/2))
}
boot_ci <- function(x, y, R = 1000) {
  v <- !is.na(x) & !is.na(y); x <- x[v]; y <- y[v]; n <- length(x)
  b <- boot(seq_len(n), function(d, i) rb_stat(x[i], y[i]), R = R)
  ci <- boot.ci(b, type = "perc")$percent[4:5]
  c(r = round(b$t0, 3), lower = round(ci[1], 3), upper = round(ci[2], 3))
}
subsample <- function(x, y, n_sub = 5000, n_rep = 1000) {
  v <- !is.na(x) & !is.na(y); x <- x[v]; y <- y[v]; N <- length(x)
  rs <- replicate(n_rep, { idx <- sample.int(N, min(n_sub, N)); rb_stat(x[idx], y[idx]) })
  c(mean_r = round(mean(rs), 3), sd_r = round(sd(rs), 3))
}
cat("\n=== CI and sensitivity (FANTASIA vs AHRD) ===\n")
ci_tab <- rbind(Corpus_IC = boot_ci(max_ic_corpus$fan, max_ic_corpus$ah),
                Topological_IC = boot_ci(max_ic_topo$fan, max_ic_topo$ah),
                GO_depth = boot_ci(max_depth$fan, max_depth$ah))
sens_tab <- rbind(Corpus_IC = subsample(max_ic_corpus$fan, max_ic_corpus$ah),
                  Topological_IC = subsample(max_ic_topo$fan, max_ic_topo$ah),
                  GO_depth = subsample(max_depth$fan, max_depth$ah))
print(ci_tab); print(sens_tab)
write.table(ci_tab,   "results/effect_size_CI_AHRD.tsv",         sep = "\t", quote = FALSE, col.names = NA)
write.table(sens_tab, "results/sensitivity_subsampling_AHRD.tsv", sep = "\t", quote = FALSE, col.names = NA)

# ============================================================
# Figure S4 (A: semantic similarity observed vs null; B: specificity)
# ============================================================
build_long <- function(df, type_label)
  data.frame(Domain = rep(c("Biological Process","Molecular Function","Cellular Component"), each = nrow(df)),
             Similarity = c(df$sim_BP, df$sim_MF, df$sim_CC), Type = type_label)
sim_comb <- rbind(build_long(sim_results, "Observed"), build_long(sim_null, "Null model"))
sim_comb <- sim_comb[!is.na(sim_comb$Similarity) & sim_comb$Similarity > 0, ]
sim_comb$Domain <- factor(sim_comb$Domain, levels = c("Biological Process","Molecular Function","Cellular Component"))
sim_comb$Type   <- factor(sim_comb$Type, levels = c("Observed","Null model"))
medians_comb    <- aggregate(Similarity ~ Domain + Type, data = sim_comb, FUN = median)

panelA <- ggplot(sim_comb, aes(x = Domain, y = Similarity, fill = Type)) +
  geom_violin(trim = TRUE, alpha = 0.6, color = NA, scale = "width", width = 0.85, position = position_dodge(width = 0.9)) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.9, color = "grey30", position = position_dodge(width = 0.9)) +
  geom_text(data = medians_comb, aes(label = paste0("Md=", round(Similarity, 2)), group = Type),
            position = position_dodge(width = 0.9), vjust = -0.6, size = 2.2, fontface = "bold", color = "grey20") +
  scale_fill_manual(values = c("Observed" = "#4BACC6", "Null model" = "#BFBFBF")) +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, 12, 2), expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Semantic similarity (BMA/Resnik IC)", fill = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = "top", legend.text = element_text(size = 8),
        axis.title = element_text(size = 8.5), axis.text = element_text(color = "grey20", size = 7),
        axis.text.x = element_text(angle = 20, hjust = 1))

barplot_data <- data.frame(
  Method = factor(rep(c("FANTASIA","AHRD"), 3), levels = c("FANTASIA","AHRD")),
  Metric = factor(rep(c("Corpus IC","Topological IC","Depth"), each = 2), levels = c("Corpus IC","Topological IC","Depth")),
  Value = c(res_corpus$mean_fan, res_corpus$mean_ah, res_topo$mean_fan, res_topo$mean_ah, res_depth$mean_fan, res_depth$mean_ah),
  SE = c(res_corpus$sd_fan/sqrt(res_corpus$n), res_corpus$sd_ah/sqrt(res_corpus$n),
         res_topo$sd_fan/sqrt(res_topo$n), res_topo$sd_ah/sqrt(res_topo$n),
         res_depth$sd_fan/sqrt(res_depth$n), res_depth$sd_ah/sqrt(res_depth$n)))
r_labels <- data.frame(
  Metric = factor(c("Corpus IC","Topological IC","Depth"), levels = c("Corpus IC","Topological IC","Depth")),
  r_text = c(paste0("|r| = ", round(res_corpus$r, 2)), paste0("|r| = ", round(res_topo$r, 2)), paste0("|r| = ", round(res_depth$r, 2))))
r_labels$y <- sapply(c("Corpus IC","Topological IC","Depth"), function(m) max(barplot_data$Value[barplot_data$Metric == m]) + 1.8)

panelB <- ggplot(barplot_data, aes(x = Metric, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = Value - SE, ymax = Value + SE), position = position_dodge(width = 0.7), width = 0.2, color = "grey30") +
  annotate("text", x = 1, y = max(barplot_data$Value[barplot_data$Metric == "Corpus IC"]) + 0.8,      label = "***", size = 3.5, color = "grey20") +
  annotate("text", x = 2, y = max(barplot_data$Value[barplot_data$Metric == "Topological IC"]) + 0.8, label = "***", size = 3.5, color = "grey20") +
  annotate("text", x = 3, y = max(barplot_data$Value[barplot_data$Metric == "Depth"]) + 0.8,          label = "***", size = 3.5, color = "grey20") +
  geom_text(data = r_labels, aes(x = Metric, y = y, label = r_text), inherit.aes = FALSE, size = 2.5, fontface = "italic", color = "grey20") +
  scale_fill_manual(values = c("FANTASIA" = "#4BACC6", "AHRD" = "#F79646")) +
  labs(x = NULL, y = "Mean max IC / depth (SE)", fill = "Method") +
  theme_classic(base_size = 9) +
  theme(legend.position = "top", legend.text = element_text(size = 8),
        axis.title = element_text(size = 8.5), axis.text = element_text(color = "grey20", size = 7))

figS4 <- panelA + panelB + plot_annotation(tag_levels = "A")
ggsave("results/Figure_S4_AHRD_similarity_specificity.jpg", figS4, width = 18, height = 10, units = "cm", dpi = 300, bg = "white")
cat("Figure S4 (A + B) saved.\n")

cat("\n============================================================\n")
cat("6_AHRD_comparison.R complete. Figures S2, S3, S4 saved to results/.\n")
cat("============================================================\n")
