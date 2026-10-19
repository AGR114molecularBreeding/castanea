#########################################################################################################
# 4_anomalous_classification.R
# Classification of taxonomically anomalous GO terms (Metazoa / Fungi) into three
# categories by biological plausibility in plants (Supplementary Table S1)
#
# Reviewer 1 - comment 4 / Reviewer 2 - minor comment 2
#
# Run AFTER 0_Load_data.R and 3_taxonomic_reliability.R, which generate the object
# 'all_gos' (per-method, per-category taxonomic classification of GO terms).
#
# Produces:
#   results/metazoa_FANTASIA_classified.tsv   per-term classification (Metazoa)
#   results/fungi_FANTASIA_classified.tsv     per-term classification (Fungi)
#   results/breakdown_anomalous.tsv           counts per method/category
#   results/Table_S1_anomalous_classification.tsv   summary table (Table S1)
#
# Postigo-Luque et al.
#########################################################################################################

# Set the working directory to the repository folder before running, e.g.:
# setwd("/path/to/Functional annotation comparative analysis FANTASIA vs Blast2Go")

source("0_Load_data.R")
source("3_taxonomic_reliability.R")   # generates 'all_gos'

library(GO.db)

dir.create("results", showWarnings = FALSE)

# ============================================================
# Keyword-based classifiers
# ============================================================
# Metazoan-associated terms
classify_metazoa <- function(go_name) {
  if (is.na(go_name)) return("Questionable")
  n <- tolower(go_name)
  impossible_kw <- c("neuron","neural","synap","dendrit","axon","perikaryon",
                     "muscle","cardiac","heart","skeletal","bone","cartilage",
                     "behavior","memory","learning","cognition",
                     "sensory","perception","vision","visual","auditory","smell",
                     "epithelial","endothelial","epidermis",
                     "insulin","hormone","dopamine","serotonin","estrogen","melanin",
                     "liver","kidney","lung","blood","immune","inflammatory",
                     "microglia","lymphocyte","leukocyte","macrophage","t cell",
                     "embryo","gastrulation","endoderm","ectoderm","mesoderm",
                     "neurotransmitter","action potential","depolarization",
                     "tumor","cancer","melanocyte","myocyte")
  bias_kw <- c("phosphoryl","kinase","binding","transport","metabolic",
               "biosynthetic","catabolic","vesicle","membrane","mitochond",
               "cytochrome","apopto","calcium","atp","dna","rna","protein",
               "proteolysis","ubiquitin","autophag","oxidation","redox",
               "translation","transcription","replication","cell cycle",
               "cell division","cytoskeleton","microtubule","actin")
  if (any(sapply(impossible_kw, function(k) grepl(k, n, fixed=TRUE)))) return("Impossible")
  if (any(sapply(bias_kw,       function(k) grepl(k, n, fixed=TRUE)))) return("Database_bias")
  return("Questionable")
}

# Fungi-associated terms
classify_fungi <- function(go_name) {
  if (is.na(go_name)) return("Questionable")
  n <- tolower(go_name)
  impossible_kw <- c("hyphal","hypha","mycelium","mycelial","sporulation",
                     "conidium","conidia","ascospore","basidium","fruiting body",
                     "chitin","ergosterol","appressorium","mating type",
                     "pseudohyphal","budding")
  bias_kw <- c("phosphoryl","kinase","binding","transport","metabolic",
               "biosynthetic","catabolic","vesicle","membrane","mitochond",
               "cytochrome","apopto","calcium","atp","dna","rna","protein",
               "proteolysis","ubiquitin","autophag","oxidation","redox",
               "translation","transcription","replication","cell cycle",
               "cell division","cytoskeleton","microtubule","actin")
  if (any(sapply(impossible_kw, function(k) grepl(k, n, fixed=TRUE)))) return("Impossible")
  if (any(sapply(bias_kw,       function(k) grepl(k, n, fixed=TRUE)))) return("Database_bias")
  return("Questionable")
}

# ============================================================
# Breakdown of anomalous assignments by method and category
# ============================================================
breakdown <- function(method, taxon_cat, classifier) {
  freq <- table(all_gos$GO_term[all_gos$Method == method & all_gos$category == taxon_cat])
  if (length(freq) == 0) return(NULL)
  names_go <- sapply(names(freq), function(x) tryCatch(Term(x), error = function(e) NA))
  cls <- sapply(names_go, classifier)
  agg <- tapply(as.integer(freq), cls, sum); agg[is.na(agg)] <- 0
  total <- sum(as.integer(freq))
  list(counts = agg, total = total, perc = round(agg / total * 100, 2))
}

cat("=== METAZOA ===\n")
cat("FANTASIA:\n"); print(breakdown("FANTASIA", "Hallucination_Metazoa", classify_metazoa))
cat("Blast2Go:\n"); print(breakdown("Blast2Go", "Hallucination_Metazoa", classify_metazoa))
cat("\n=== FUNGI ===\n")
cat("FANTASIA:\n"); print(breakdown("FANTASIA", "Hallucination_Fungi", classify_fungi))
cat("Blast2Go:\n"); print(breakdown("Blast2Go", "Hallucination_Fungi", classify_fungi))

# ============================================================
# Breakdown table (counts per method / category)
# ============================================================
make_row <- function(method, taxon_cat, classifier, label) {
  b <- breakdown(method, taxon_cat, classifier)
  if (is.null(b)) return(NULL)
  cats <- c("Database_bias", "Questionable", "Impossible")
  n <- sapply(cats, function(c) ifelse(c %in% names(b$counts), b$counts[c], 0))
  p <- sapply(cats, function(c) ifelse(c %in% names(b$perc),   b$perc[c],   0))
  data.frame(Method = method, Taxon = label,
             Database_bias_n = n[1], Database_bias_pct = p[1],
             Questionable_n  = n[2], Questionable_pct  = p[2],
             Impossible_n    = n[3], Impossible_pct    = p[3],
             Total = b$total)
}
breakdown_table <- rbind(
  make_row("FANTASIA", "Hallucination_Metazoa", classify_metazoa, "Metazoa"),
  make_row("Blast2Go", "Hallucination_Metazoa", classify_metazoa, "Metazoa"),
  make_row("FANTASIA", "Hallucination_Fungi",   classify_fungi,   "Fungi"),
  make_row("Blast2Go", "Hallucination_Fungi",   classify_fungi,   "Fungi")
)
print(breakdown_table, row.names = FALSE)
write.table(breakdown_table, "results/breakdown_anomalous.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ============================================================
# Per-term classification tables (for traceability)
# ============================================================
make_classified <- function(taxon_cat, classifier, file) {
  freq <- sort(table(all_gos$GO_term[all_gos$Method == "FANTASIA" &
                                       all_gos$category == taxon_cat]), decreasing = TRUE)
  df <- data.frame(
    GO_term = names(freq),
    GO_name = sapply(names(freq), function(x) tryCatch(Term(x), error = function(e) NA)),
    Freq    = as.integer(freq), stringsAsFactors = FALSE)
  df$Class <- sapply(df$GO_name, classifier)
  write.table(df, file.path("results", file), sep = "\t", row.names = FALSE, quote = FALSE)
  df
}
meta_classified  <- make_classified("Hallucination_Metazoa", classify_metazoa, "metazoa_FANTASIA_classified.tsv")
fungi_classified <- make_classified("Hallucination_Fungi",   classify_fungi,   "fungi_FANTASIA_classified.tsv")

# ============================================================
# Supplementary Table S1 (summary with representative examples)
# ------------------------------------------------------------
# NOTE: the counts/percentages below are the values obtained from the breakdown
# above; the "Examples" column is curated. Keep it in sync with the final
# Supplementary Table S1 in the manuscript.
# ============================================================
tabla_s1 <- data.frame(
  Taxon = c("Metazoa","Metazoa","Metazoa","Metazoa",
            "Fungi","Fungi","Fungi","Fungi"),
  Category = c("Biologically implausible","Questionable","Database bias","Subtotal",
               "Biologically implausible","Questionable","Database bias","Subtotal"),
  FANTASIA_n = c(5075, 3109, 2715, 10899, 319, 172, 297, 788),
  FANTASIA_pct = c(46.6, 28.5, 24.9, NA, 40.5, 21.8, 37.7, NA),
  Blast2GO_n = c(3, 118, 439, 560, 0, 0, 117, 117),
  Blast2GO_pct = c(0.5, 21.1, 78.4, NA, 0, 0, 100, NA),
  Examples = c(
    "learning (GO:0007612); negative regulation of insulin secretion (GO:0046676)",
    "peptidyl-L-cysteine S-palmitoylation (GO:0018230); negative regulation of tubulin deacetylation (GO:1904428)",
    "peptidyl-tyrosine phosphorylation (GO:0018108); cytochrome c release from mitochondria (GO:0090200)",
    "",
    "ascospore formation (GO:0030437); regulation of ergosterol biosynthetic process (GO:0032443)",
    "acyl-CoA ceramide synthase complex (GO:0061576); peptidyl-histidine methylation, to form tele-methylhistidine (GO:0042038)",
    "terpenoid biosynthetic process, mevalonate-dependent (GO:0051485); negative regulation of anaphase-promoting complex-dependent catabolic process (GO:1905785)",
    ""
  ),
  stringsAsFactors = FALSE
)
write.table(tabla_s1, "results/Table_S1_anomalous_classification.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nTable S1 saved to: results/Table_S1_anomalous_classification.tsv\n")
print(tabla_s1, row.names = FALSE)

cat("\n============================================================\n")
cat("4_anomalous_classification.R complete.\n")
cat("============================================================\n")
