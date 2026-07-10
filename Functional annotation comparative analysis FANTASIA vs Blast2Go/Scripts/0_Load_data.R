#########################################################################################################
# 0_ Load data
# Data loading and preprocessing
# Postigo-Luque et al. 
#########################################################################################################

library(Biostrings)
library(ggplot2)
library(GO.db)

#### File paths (modify based on the location of the data)
BLAST2GO_FILE <-"data/GO_Hap1_JVD.txt"
FANTASIA_FILE <-"data/Csath1_FANTASIA_TopGO.txt"
FASTA_FILE <- "data/Cast.1_0.hap1.clean.fa"
GO_OBO_FILE <- "data/go.obo"


                               #### Load Blast2GO annotations ####
homologia_raw <- read.table(
  BLAST2GO_FILE, sep = "\t", header = TRUE,
  stringsAsFactors = FALSE, quote = "", fill = TRUE)

colnames(homologia_raw) <- c("Sequence.Name", "Annotation.GO.ID", "Annotation.GO.Term", "Annotation.GO.Category")

homologia <- homologia_raw[
  !is.na(homologia_raw$Annotation.GO.ID) &
    homologia_raw$Annotation.GO.ID != "",]


                              #### Load FANTASIA annotations ####
fantasia <- read.table(
  FANTASIA_FILE, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(fantasia) <- c("ProteinID", "GOterms")


                             #### Harmonise protein IDs ####
homologia$Sequence.Name <- gsub("-mRNA-\\d+$","", homologia$Sequence.Name)
homologia$Sequence.Name <-gsub("\\.t\\d+$",    "", homologia$Sequence.Name)
fantasia$ProteinID      <- gsub("\\..*$",        "", fantasia$ProteinID)


                      #### Load proteome and define annotation sets #####
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
    n_fantasia, "FANTASIA |", n_homologia, "Blast2Go |",
    n_both, "both\n")



