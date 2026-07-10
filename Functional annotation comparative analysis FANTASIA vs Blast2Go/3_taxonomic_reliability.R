# ============================================================
# 3_taxonomic_reliability.R
# Taxonomic reliability analysis of GO annotations (Table 1)
# Postigo-Luque et al. — FANTASIA vs Blast2Go in C. sativa
#
# Run AFTER 02_taxonomic_classification.py has produced
# results/taxon_map_v2.tsv
# ============================================================

source("0_Load_data.R")
library(ontologyIndex)

# ============================================================
# Load taxon map from Python output
# ============================================================
taxon_map <- read.table(
  "results/taxon_map_v2.tsv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE,
  fill = TRUE, quote = ""
)
cat("Taxon map loaded:", nrow(taxon_map), "GO terms\n")
taxon_lookup <- setNames(taxon_map$taxon, taxon_map$GO_term)

# ============================================================
# Classification function
# ============================================================
general_noise_taxa <- c("root", "cellular organisms", "Eukaryota", "Opisthokonta")

valid_plant_keywords <- c("viridiplantae", "streptophyta", "streptophytina",
                          "embryophyta", "tracheophyta", "euphyllophyta",
                          "spermatophyta", "magnoliopsida", "mesangiospermae",
                          "eudicotyledons", "gunneridae", "pentapetalae",
                          "rosids", "asterids", "malvids", "brassicales",
                          "brassicaceae", "camelineae", "arabidopsis",
                          "poaceae", "fabaceae", "fagaceae", "castanea",
                          "solanaceae", "solanoideae", "solanum", "nicotiana",
                          "glycine", "gossypium", "oryza", "vitis", "zea",
                          "sorghum", "petrosaviidae", "myrtales",
                          "npaaa clade", "aparaglossata")

metazoa_keywords <- c("metazoa", "eumetazoa", "bilateria", "deuterostomia",
                      "chordata", "craniata", "vertebrata", "gnathostomata",
                      "teleostomi", "euteleostomi", "sarcopterygii",
                      "dipnotetrapodomorpha", "tetrapoda", "amniota",
                      "mammalia", "theria", "eutheria", "boreoeutheria",
                      "euarchontoglires", "primates", "murinae",
                      "homo sapiens", "mus musculus", "rattus norvegicus",
                      "danio rerio", "gallus gallus",
                      "drosophila", "diptera", "schizophora", "endopterygota",
                      "neoptera", "dicondylia", "pterygota",
                      "allotriocarida", "arthropoda",
                      "caenorhabditis", "nematoda", "rhabditida",
                      "protostomia", "ecdysozoa")

fungi_keywords <- c("fungi", "dikarya", "ascomycota", "basidiomycota",
                    "saccharomycotina", "saccharomyceta", "saccharomycetes",
                    "saccharomycetaceae", "saccharomyces cerevisiae",
                    "candida", "schizosaccharomyces", "taphrinomycotina",
                    "leotiomyceta", "aspergillus", "aspergillaceae")

bacteria_keywords <- c("bacteria", "pseudomonadota", "pseudomonadati",
                       "firmicutes", "bacillaceae", "actinomycetota",
                       "escherichia", "mycobacterium", "mycoplasmoides",
                       "chlorobiaceae", "chlorobiales")

archaea_keywords <- c("archaea", "euryarchaeota", "crenarchaeota", "thermoprotei")

virus_keywords <- c("viruses", "caudoviricetes", "orthornavirae", "riboviria",
                    "duplornaviricota", "polyploviricotina",
                    "adenoviridae", "mastadenovirus",
                    "alphaherpesvirinae", "orthoherpesviridae",
                    "autographa californica", "tequatrovirus",
                    "rotavirus", "moloney murine leukemia virus",
                    "human herpesvirus")

classify_taxon <- function(taxon) {
  if (is.na(taxon) || taxon == "" || taxon == "Unknown") return("Unknown")
  taxon_lower <- tolower(taxon)
  if (taxon %in% general_noise_taxa)                                        return("General_Noise")
  if (any(sapply(valid_plant_keywords, function(k) grepl(k, taxon_lower)))) return("Valid_Plant")
  if (taxon == "Non taxon-specific")                                        return("Valid_Non_Specific")
  if (any(sapply(metazoa_keywords,     function(k) grepl(k, taxon_lower)))) return("Hallucination_Metazoa")
  if (any(sapply(fungi_keywords,       function(k) grepl(k, taxon_lower)))) return("Hallucination_Fungi")
  if (any(sapply(bacteria_keywords,    function(k) grepl(k, taxon_lower)))) return("Contamination_Bacteria")
  if (any(sapply(archaea_keywords,     function(k) grepl(k, taxon_lower)))) return("Contamination_Archaea")
  if (any(sapply(virus_keywords,       function(k) grepl(k, taxon_lower)))) return("Contamination_Viruses")
  return("Other_Anomaly")
}

# ============================================================
# Classify all GO term assignments by method
# ============================================================
fantasia_all_gos <- data.frame(
  GO_term = trimws(unlist(strsplit(fantasia$GOterms, ","))),
  Method  = "FANTASIA"
)
homologia_all_gos <- data.frame(
  GO_term = homologia$Annotation.GO.ID,
  Method  = "Blast2Go"
)
all_gos <- rbind(fantasia_all_gos, homologia_all_gos)
all_gos <- all_gos[all_gos$GO_term != "", ]
all_gos$taxon    <- taxon_lookup[all_gos$GO_term]
all_gos$category <- sapply(all_gos$taxon, classify_taxon)

cat("\nClassification by method:\n")
print(table(all_gos$Method, all_gos$category))

# ============================================================
# Summary table by method
# ============================================================
all_categories <- c("General_Noise", "Valid_Non_Specific", "Valid_Plant",
                    "Hallucination_Metazoa", "Hallucination_Fungi",
                    "Contamination_Bacteria", "Contamination_Archaea",
                    "Contamination_Viruses", "Other_Anomaly", "Unknown")

summary_taxon <- do.call(rbind, lapply(c("FANTASIA", "Blast2Go"), function(method) {
  df      <- all_gos[all_gos$Method == method, ]
  n_total <- nrow(df)
  counts  <- table(factor(df$category, levels = all_categories))
  data.frame(
    Method                 = method,
    Total_assignments      = n_total,
    General_Noise          = as.integer(counts["General_Noise"]),
    Valid_Plant            = as.integer(counts["Valid_Plant"]),
    Hallucination_Metazoa  = as.integer(counts["Hallucination_Metazoa"]),
    Hallucination_Fungi    = as.integer(counts["Hallucination_Fungi"]),
    Contamination_Bacteria = as.integer(counts["Contamination_Bacteria"]),
    Contamination_Archaea  = as.integer(counts["Contamination_Archaea"]),
    Contamination_Viruses  = as.integer(counts["Contamination_Viruses"]),
    Other_Anomaly          = as.integer(counts["Other_Anomaly"]),
    Unknown                = as.integer(counts["Unknown"]),
    Perc_General_Noise     = round(counts["General_Noise"]         / n_total * 100, 2),
    Perc_Valid_Plant       = round(counts["Valid_Plant"]            / n_total * 100, 2),
    Perc_Hall_Metazoa      = round(counts["Hallucination_Metazoa"] / n_total * 100, 2),
    Perc_Hall_Fungi        = round(counts["Hallucination_Fungi"]   / n_total * 100, 2),
    Perc_Contamination     = round((counts["Contamination_Bacteria"] +
                                      counts["Contamination_Archaea"] +
                                      counts["Contamination_Viruses"]) / n_total * 100, 2),
    Perc_Other_Anomaly     = round(counts["Other_Anomaly"]         / n_total * 100, 2),
    Perc_Unknown           = round(counts["Unknown"]               / n_total * 100, 2)
  )
}))
summary_taxon[is.na(summary_taxon)] <- 0

# ============================================================
# Unknown breakdown (Obsolete / Active / Not in OBO)
# ============================================================
cat("\n=== Unknown GO terms analysis ===\n")

go_obo <- get_ontology(GO_OBO_FILE, extract_tags = "everything")

unknown_terms <- taxon_map$GO_term[taxon_map$taxon == "Unknown"]
cat("Unknown GO terms:", length(unknown_terms), "\n")

clasificacion_unknown <- data.frame(
  GO_term  = unknown_terms,
  en_obo   = unknown_terms %in% go_obo$id,
  obsoleto = NA,
  stringsAsFactors = FALSE
)
for (i in seq_along(unknown_terms)) {
  term <- unknown_terms[i]
  if (term %in% go_obo$id) {
    clasificacion_unknown$obsoleto[i] <- isTRUE(go_obo$obsolete[term])
  }
}

cat("  Obsolete:          ", sum(clasificacion_unknown$obsoleto == TRUE,  na.rm = TRUE), "\n")
cat("  Active (no UniProt):", sum(clasificacion_unknown$obsoleto == FALSE, na.rm = TRUE), "\n")
cat("  Not in OBO:        ", sum(!clasificacion_unknown$en_obo), "\n")

# Map Unknown type per GO term
unknown_type_map <- setNames(
  ifelse(!clasificacion_unknown$en_obo, "Unknown - Not in OBO",
         ifelse(clasificacion_unknown$obsoleto == TRUE, "Unknown - Obsolete",
                "Unknown - Active")),
  clasificacion_unknown$GO_term
)

# Count Unknown assignments per method and subtype
all_gos_unknown <- all_gos[all_gos$category == "Unknown", ]
all_gos_unknown$unknown_type <- unknown_type_map[all_gos_unknown$GO_term]
all_gos_unknown$unknown_type[is.na(all_gos_unknown$unknown_type)] <- "Unknown - Not in OBO"

n_fan <- nrow(all_gos[all_gos$Method == "FANTASIA", ])
n_hom <- nrow(all_gos[all_gos$Method == "Blast2Go", ])

get_unknown_counts <- function(method) {
  df <- all_gos_unknown[all_gos_unknown$Method == method, ]
  c(
    Obsolete   = sum(df$unknown_type == "Unknown - Obsolete",  na.rm = TRUE),
    Active     = sum(df$unknown_type == "Unknown - Active",    na.rm = TRUE),
    Not_in_OBO = sum(df$unknown_type == "Unknown - Not in OBO", na.rm = TRUE)
  )
}
unk_fan <- get_unknown_counts("FANTASIA")
unk_hom <- get_unknown_counts("Blast2Go")

# ============================================================
# Table 1: Consolidated table for paper
# ============================================================
paper_table <- data.frame(
  Category = c(
    "General Noise",
    "Valid - Viridiplantae",
    "Hallucination - Metazoa",
    "Hallucination - Fungi",
    "Contamination",
    "Other anomaly",
    "Unknown - Obsolete",
    "Unknown - Active",
    "Unknown - Not in OBO"
  ),
  FANTASIA_n = c(
    summary_taxon$General_Noise[1],
    summary_taxon$Valid_Plant[1],
    summary_taxon$Hallucination_Metazoa[1],
    summary_taxon$Hallucination_Fungi[1],
    summary_taxon$Contamination_Bacteria[1] +
      summary_taxon$Contamination_Archaea[1] +
      summary_taxon$Contamination_Viruses[1],
    summary_taxon$Other_Anomaly[1],
    unk_fan["Obsolete"],
    unk_fan["Active"],
    unk_fan["Not_in_OBO"]
  ),
  FANTASIA_perc = c(
    summary_taxon$Perc_General_Noise[1],
    summary_taxon$Perc_Valid_Plant[1],
    summary_taxon$Perc_Hall_Metazoa[1],
    summary_taxon$Perc_Hall_Fungi[1],
    summary_taxon$Perc_Contamination[1],
    summary_taxon$Perc_Other_Anomaly[1],
    round(unk_fan["Obsolete"]   / n_fan * 100, 2),
    round(unk_fan["Active"]     / n_fan * 100, 2),
    round(unk_fan["Not_in_OBO"] / n_fan * 100, 2)
  ),
  Blast2Go_n = c(
    summary_taxon$General_Noise[2],
    summary_taxon$Valid_Plant[2],
    summary_taxon$Hallucination_Metazoa[2],
    summary_taxon$Hallucination_Fungi[2],
    summary_taxon$Contamination_Bacteria[2] +
      summary_taxon$Contamination_Archaea[2] +
      summary_taxon$Contamination_Viruses[2],
    summary_taxon$Other_Anomaly[2],
    unk_hom["Obsolete"],
    unk_hom["Active"],
    unk_hom["Not_in_OBO"]
  ),
  Blast2Go_perc = c(
    summary_taxon$Perc_General_Noise[2],
    summary_taxon$Perc_Valid_Plant[2],
    summary_taxon$Perc_Hall_Metazoa[2],
    summary_taxon$Perc_Hall_Fungi[2],
    summary_taxon$Perc_Contamination[2],
    summary_taxon$Perc_Other_Anomaly[2],
    round(unk_hom["Obsolete"]   / n_hom * 100, 2),
    round(unk_hom["Active"]     / n_hom * 100, 2),
    round(unk_hom["Not_in_OBO"] / n_hom * 100, 2)
  )
)

cat("\n=== TABLE 1 (for paper) ===\n")
print(paper_table)
write.table(paper_table, "results/Table_1_taxonomic_classification.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nTable 1 saved to: results/Table_1_taxonomic_classification.tsv\n")

cat("\n============================================================\n")
cat("03_taxonomic_reliability.R complete.\n")
cat("============================================================\n")