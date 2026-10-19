#########################################################################################################
# 0_Load_data.R
# Data loading and preprocessing for the FANTASIA vs Blast2GO comparison in C. sativa
#
# Loads the FANTASIA and Blast2GO annotations, harmonises protein IDs, propagates
# FANTASIA predictions to identical sequences collapsed by CD-HIT, and defines the
# annotation sets and per-protein GO lists used by the downstream scripts
# (1_Comparative_analysis.R, 3_taxonomic_reliability.R, 4_anomalous_classification.R).
#
# This script is sourced by the others; it is not meant to be run on its own.
#
# Required input files (place them in a ./data/ folder, or edit the paths below):
#   data/GO_Hap1_JVD.txt                       Blast2GO annotations (OmicsBox v3.2.9)
#   data/Csath1_FANTASIA_TopGO.txt             FANTASIA annotations (v1.0, ProtT5)
#   data/Cast.1_0.hap1.clean.fa                C. sativa hap1 proteome (filtered)
#   data/Cast.1_0.hap1.clean_cdhit100.pep.clstr   CD-HIT 100%-identity clusters
#   data/go.obo                                Gene Ontology (go.obo)
#
# Postigo-Luque et al.
#########################################################################################################

library(Biostrings)
library(ggplot2)
library(GO.db)

# Set the working directory to the repository folder before sourcing, e.g.:
# setwd("/path/to/Functional annotation comparative analysis FANTASIA vs Blast2Go")

#### File paths ####
BLAST2GO_FILE <- "data/GO_Hap1_JVD.txt"
FANTASIA_FILE <- "data/Csath1_FANTASIA_TopGO.txt"
FASTA_FILE    <- "data/Cast.1_0.hap1.clean.fa"
CLSTR_FILE    <- "data/Cast.1_0.hap1.clean_cdhit100.pep.clstr"
GO_OBO_FILE   <- "data/go.obo"

#### Load Blast2GO annotations ####
homologia_raw <- read.table(
  BLAST2GO_FILE, sep = "\t", header = TRUE,
  stringsAsFactors = FALSE, quote = "", fill = TRUE)

colnames(homologia_raw) <- c("Sequence.Name", "Annotation.GO.ID",
                             "Annotation.GO.Term", "Annotation.GO.Category")

homologia <- homologia_raw[
  !is.na(homologia_raw$Annotation.GO.ID) &
    homologia_raw$Annotation.GO.ID != "", ]

#### Load FANTASIA annotations ####
fantasia <- read.table(
  FANTASIA_FILE, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(fantasia) <- c("ProteinID", "GOterms")

#### Harmonise protein IDs ####
homologia$Sequence.Name <- gsub("-mRNA-\\d+$", "", homologia$Sequence.Name)
homologia$Sequence.Name <- gsub("\\.t\\d+$",    "", homologia$Sequence.Name)
fantasia$ProteinID      <- gsub("\\..*$",        "", fantasia$ProteinID)

#### Propagate FANTASIA annotations to identical sequences ####
# CD-HIT (100% identity) collapsed identical sequences before FANTASIA.
# Since identical sequences yield identical embeddings and GO predictions,
# we propagate each representative's FANTASIA annotation to its identical copies,
# restoring symmetry with Blast2GO (run on the full proteome).
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
      .rescatable <- rbind(.rescatable,
                           data.frame(removed = nf, representative = en_fan[1]))
    }
  }
}

# Build propagated rows: each removed copy gets its representative's GO terms
.prop_rows <- data.frame(
  ProteinID = .rescatable$removed,
  GOterms   = fantasia$GOterms[match(.rescatable$representative, fantasia$ProteinID)],
  stringsAsFactors = FALSE
)
.prop_rows <- .prop_rows[!is.na(.prop_rows$GOterms), ]

fantasia <- rbind(fantasia, .prop_rows)
cat("Propagated FANTASIA annotations to", nrow(.prop_rows), "identical sequences\n")

#### Load proteome and define annotation sets ####
fasta_seqs <- readAAStringSet(FASTA_FILE)
all_ids    <- names(fasta_seqs)
n_total    <- length(all_ids)

homologia_annotated <- unique(homologia$Sequence.Name)
fantasia_annotated  <- unique(fantasia$ProteinID)

both_annotated <- intersect(fantasia_annotated, homologia_annotated)
only_fantasia  <- setdiff(fantasia_annotated,  homologia_annotated)
only_homologia <- setdiff(homologia_annotated, fantasia_annotated)

n_fantasia  <- length(fantasia_annotated)
n_homologia <- length(homologia_annotated)
n_both      <- length(both_annotated)

#### GO term lists per protein ####
fantasia_go_list <- strsplit(fantasia$GOterms, ",")
fantasia_go_list <- lapply(fantasia_go_list, function(x) unique(trimws(x)))
names(fantasia_go_list) <- fantasia$ProteinID

homologia_go_list <- split(homologia$Annotation.GO.ID, homologia$Sequence.Name)
homologia_go_list <- lapply(homologia_go_list, unique)

#### Common proteins ####
proteinas_comunes <- intersect(names(fantasia_go_list), names(homologia_go_list))
fantasia_comun    <- fantasia_go_list[proteinas_comunes]
homologia_comun   <- homologia_go_list[proteinas_comunes]

#### GO domain mapping ####
fantasia_go_largo  <- unique(trimws(unlist(strsplit(fantasia$GOterms, ","))))
homologia_go_largo <- unique(homologia$Annotation.GO.ID)
todos_go           <- union(fantasia_go_largo, homologia_go_largo)

dominio_df  <- AnnotationDbi::select(GO.db, keys = todos_go,
                                     columns = "ONTOLOGY", keytype = "GOID")
dominio_map <- setNames(dominio_df$ONTOLOGY, dominio_df$GOID)

cat("Data loaded:", n_total, "proteins |",
    n_fantasia, "FANTASIA |", n_homologia, "Blast2GO |",
    n_both, "both\n")
