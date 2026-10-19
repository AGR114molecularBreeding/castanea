#########################################################################################################
# 5_enrichment_RNAseq.R
# Functional enrichment use case: C. sativa response to Phytophthora cinnamomi
# FANTASIA vs Blast2GO vs AHRD on the same set of DEGs and the same background.
#
# Data: Fernandes et al. 2024 (Front. Plant Sci.), BioProject PRJNA1108966.
#       C. sativa roots, control vs inoculated with P. cinnamomi (oomycete), 2 hai.
#       3 control + 3 infected (Illumina NovaSeq, paired-end).
# Upstream pipeline (not in this script): SRA -> fastp -> HISAT2 (hap1 genome)
#       -> featureCounts. The featureCounts table is the input below.
#
# Input (place in ./data/):
#   data/counts.txt                  featureCounts gene-level counts
#   data/Csath1_FANTASIA_TopGO.txt   FANTASIA annotation
#   data/GO_Hap1_JVD.txt             Blast2GO annotation
#   data/Casat.proteins.funct_ahrd.tsv   AHRD annotation
#
# Produces: Figure 3 (FANTASIA vs Blast2GO), Figure S5 (FANTASIA vs AHRD),
#           and the enrichment / DEG result tables in results/.
#
# Postigo-Luque et al.
#########################################################################################################

# Set the working directory to the repository folder before running, e.g.:
# setwd("/path/to/Functional annotation comparative analysis FANTASIA vs Blast2Go")

library(DESeq2)
library(clusterProfiler)
library(GO.db)
library(ggplot2)

dir.create("results", showWarnings = FALSE)

COUNTS_FILE   <- "data/counts.txt"
FANTASIA_FILE <- "data/Csath1_FANTASIA_TopGO.txt"
BLAST2GO_FILE <- "data/GO_Hap1_JVD.txt"
AHRD_FILE     <- "data/Casat.proteins.funct_ahrd.tsv"

# ============================================================
# 1. DESeq2 - differential expression (control vs infected, 2 hai)
# ============================================================
counts_raw <- read.table(COUNTS_FILE, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, comment.char = "#")

# featureCounts: col 1 = Geneid, cols 2-6 = Chr/Start/End/Strand/Length, cols 7-12 = samples
counts <- counts_raw[, 7:12]
rownames(counts) <- counts_raw$Geneid
colnames(counts) <- c("control_1","control_2","control_3",
                      "infected_1","infected_2","infected_3")
cat("Genes:", nrow(counts), "| Samples:", ncol(counts), "\n")

condition <- factor(c("control","control","control",
                      "infected","infected","infected"),
                    levels = c("control","infected"))
coldata <- data.frame(row.names = colnames(counts), condition = condition)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]              # low-expression filter
cat("Genes after filtering:", nrow(dds), "\n")
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition","infected","control"))
summary(res)

# ============================================================
# 2. Extract DEGs (padj < 0.05 and |log2FC| >= 1)
# ============================================================
res_df <- as.data.frame(res); res_df$gene <- rownames(res_df)
res_df <- res_df[!is.na(res_df$padj), ]
degs   <- res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1, ]
degs_up   <- degs[degs$log2FoldChange >=  1, ]
degs_down <- degs[degs$log2FoldChange <= -1, ]

cat("=== DEGs (padj<0.05, |log2FC|>=1) ===\n")
cat("Total DEGs:", nrow(degs), "| Up:", nrow(degs_up), "| Down:", nrow(degs_down), "\n")

write.table(res_df, "results/DESeq2_all_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(degs,   "results/DEGs_significativos.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(degs_up$gene,   "results/DEGs_up_ids.txt",   quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(degs_down$gene, "results/DEGs_down_ids.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

deg_genes      <- degs$gene
universe_genes <- rownames(dds)      # background = all expressed genes

# ============================================================
# 3. Build TERM2GENE tables for the three annotation methods
# ============================================================
# --- FANTASIA ---
fan <- read.table(FANTASIA_FILE, sep = "\t", stringsAsFactors = FALSE,
                  quote = "", header = FALSE)
colnames(fan) <- c("gene","GOs")
fan_list <- strsplit(fan$GOs, ",")
fan_t2g  <- data.frame(GO   = trimws(unlist(fan_list)),
                       gene = rep(fan$gene, sapply(fan_list, length)),
                       stringsAsFactors = FALSE)
fan_t2g  <- fan_t2g[fan_t2g$GO != "", ]
cat("FANTASIA:", length(unique(fan_t2g$gene)), "genes,",
    length(unique(fan_t2g$GO)), "GO terms\n")

# --- Blast2GO ---
b2g <- read.table(BLAST2GO_FILE, sep = "\t", header = TRUE,
                  stringsAsFactors = FALSE, quote = "", fill = TRUE)
colnames(b2g) <- c("gene","GO","term","category")
b2g_t2g <- data.frame(GO = b2g$GO, gene = b2g$gene, stringsAsFactors = FALSE)
b2g_t2g <- b2g_t2g[b2g_t2g$GO != "" & !is.na(b2g_t2g$GO), ]
cat("Blast2GO:", length(unique(b2g_t2g$gene)), "genes,",
    length(unique(b2g_t2g$GO)), "GO terms\n")

# --- AHRD ---
ahrd <- read.table(AHRD_FILE, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                   quote = "", fill = TRUE, comment.char = "#")
# Harmonise IDs to gene level (strip isoform suffixes) so they match the DEGs
ahrd$gene_clean <- gsub("-mRNA-\\d+$", "", ahrd$Protein.Accession)
ahrd$gene_clean <- gsub("\\.t\\d+$",   "", ahrd$gene_clean)
ahrd_go   <- ahrd[ahrd$Gene.Ontology.Term != "" & !is.na(ahrd$Gene.Ontology.Term), ]
ahrd_list <- strsplit(ahrd_go$Gene.Ontology.Term, ",")
ahrd_t2g  <- data.frame(GO   = trimws(unlist(ahrd_list)),
                        gene = rep(ahrd_go$gene_clean, sapply(ahrd_list, length)),
                        stringsAsFactors = FALSE)
ahrd_t2g  <- ahrd_t2g[ahrd_t2g$GO != "", ]
cat("AHRD:", length(unique(ahrd_t2g$gene)), "genes,",
    length(unique(ahrd_t2g$GO)), "GO terms\n")

# ============================================================
# 4. GO enrichment of the DEGs with each annotation
#    (same genes, same background; only the annotation changes)
# ============================================================
run_enrich <- function(t2g) {
  enricher(gene = deg_genes, TERM2GENE = t2g[, c("GO","gene")],
           universe = universe_genes, pvalueCutoff = 0.05,
           pAdjustMethod = "BH", qvalueCutoff = 0.2)
}
enrich_fan  <- run_enrich(fan_t2g)
enrich_b2g  <- run_enrich(b2g_t2g)
enrich_ahrd <- run_enrich(ahrd_t2g)

fan_res  <- as.data.frame(enrich_fan)
b2g_res  <- as.data.frame(enrich_b2g)
ahrd_res <- as.data.frame(enrich_ahrd)

cat("=== Enriched GO terms (FDR < 0.05) ===\n")
cat("FANTASIA:", nrow(fan_res), "| Blast2GO:", nrow(b2g_res), "| AHRD:", nrow(ahrd_res), "\n")

# Add GO names and save
add_names <- function(df) {
  if (nrow(df) == 0) return(df)
  df$name <- sapply(df$ID, function(x){ n <- tryCatch(Term(x), error = function(e) NA)
                                        ifelse(is.null(n) || is.na(n), x, n) })
  df
}
fan_n  <- add_names(fan_res)
b2g_n  <- add_names(b2g_res)
ahrd_n <- add_names(ahrd_res)
write.table(fan_n,  "results/enrichment_FANTASIA.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(b2g_n,  "results/enrichment_Blast2GO.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ahrd_n, "results/enrichment_AHRD.tsv",     sep = "\t", quote = FALSE, row.names = FALSE)

# ============================================================
# 5. Summary table
# ============================================================
go_hom       <- union(b2g_res$ID, ahrd_res$ID)
solo_fan     <- setdiff(fan_res$ID, go_hom)
genes_fan    <- unique(fan_t2g$gene); genes_b2g <- unique(b2g_t2g$gene); genes_ahrd <- unique(ahrd_t2g$gene)
deg_solo_fan <- deg_genes[deg_genes %in% genes_fan & !(deg_genes %in% genes_b2g) & !(deg_genes %in% genes_ahrd)]

resumen <- data.frame(
  Metric = c("DEGs total","DEGs up","DEGs down",
             "Enriched GO FANTASIA","Enriched GO Blast2GO","Enriched GO AHRD",
             "Only FANTASIA (vs Blast2GO+AHRD)","DEGs annotated only by FANTASIA"),
  Value  = c(nrow(degs), nrow(degs_up), nrow(degs_down),
             nrow(fan_res), nrow(b2g_res), nrow(ahrd_res),
             length(solo_fan), length(deg_solo_fan))
)
write.table(resumen, "results/summary_RNAseq.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(resumen, row.names = FALSE)

# ============================================================
# 6. Enrichment figures (defense-related GO terms)
# ============================================================
defensa_kw <- c("defense","defence","biotic","immune","pathogen","oomycete",
                "fungus","fungal","hypersensitive","jasmon","salicyl","ethylene",
                "brassinosteroid","chitin","glucan","lignin","pectin","callose",
                "hydrogen peroxide","oxidative","peroxidase","peptidase inhibitor",
                "endopeptidase inhibitor","pathogenesis","thaumatin","chitinase",
                "wounding","phytotoxin","microbial")
es_def  <- function(nm) any(sapply(defensa_kw, function(k) grepl(k, tolower(nm), fixed = TRUE)))
excluir <- c("GO:0042742","GO:0051607")   # defense response to bacterium / virus (not an oomycete)

mk_dot <- function(df, titulo, archivo, metodos, w = 16) {
  df$defensa <- sapply(df$name, es_def)
  ids <- setdiff(unique(df$ID[df$defensa]), excluir)
  df  <- df[df$ID %in% ids, ]
  ord <- aggregate(p.adjust ~ name, df, min); ord <- ord[order(ord$p.adjust), ]
  df$name   <- factor(df$name, levels = rev(ord$name))
  df$Method <- factor(df$Method, levels = metodos)
  p <- ggplot(df, aes(x = Method, y = name, size = Count, color = -log10(p.adjust))) +
    geom_point() +
    scale_color_gradient(low = "#F79646", high = "#C0392B", name = "-log10(FDR)") +
    scale_size(range = c(2, 7), name = "Gene count") +
    labs(x = NULL, y = NULL, title = titulo) +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 8),
          plot.title  = element_text(size = 10, face = "bold"))
  ggsave(file.path("results", archivo), p, width = w, height = 13,
         units = "cm", dpi = 300, bg = "white")
}

# Figure 3: FANTASIA vs Blast2GO
d_fb <- rbind(data.frame(fan_n[, c("ID","name","Count","p.adjust")], Method = "FANTASIA"),
              data.frame(b2g_n[, c("ID","name","Count","p.adjust")], Method = "Blast2GO"))
mk_dot(d_fb, "Enriched defense-related GO terms (FANTASIA vs Blast2GO)",
       "Figure_3_enrichment_Blast2GO.jpg", c("FANTASIA","Blast2GO"))

# Figure S5: FANTASIA vs AHRD
d_fa <- rbind(data.frame(fan_n[, c("ID","name","Count","p.adjust")], Method = "FANTASIA"),
              data.frame(ahrd_n[, c("ID","name","Count","p.adjust")], Method = "AHRD"))
mk_dot(d_fa, "Enriched defense-related GO terms (FANTASIA vs AHRD)",
       "Figure_S5_enrichment_AHRD.jpg", c("FANTASIA","AHRD"))

cat("\n============================================================\n")
cat("5_enrichment_RNAseq.R complete. Figures 3 and S5 saved to results/.\n")
cat("============================================================\n")
