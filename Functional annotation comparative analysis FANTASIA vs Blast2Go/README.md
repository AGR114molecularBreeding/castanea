# FANTASIA vs Blast2Go: Protein language model versus homology-based functional annotation in non-model forest trees: a case study in *Castanea sativa* Mill.

Companion repository for:

> **Postigo-Luque, R. et al.** Protein language model versus homology-based functional annotation in non-model forest trees: a case study in *Castanea sativa* Mill.

This repository contains the functional annotations, analysis scripts, and output files used in the study comparing [FANTASIA](https://github.com/FANTASIA-Annotation/FANTASIA) (protein language model-based) and Blast2Go (homology-based) functional annotation of the *Castanea sativa* (European chestnut) haplotype 1 proteome.

## Repository structure

```
├── annotations/
│   ├── Csath1_FANTASIA_TopGO.txt       # FANTASIA GO annotations (v1.0, ProtT5)
│   └── GO_Hap1_JVD.txt                 # Blast2Go GO annotations (OmicsBox v3.2.9)
├── scripts/
│   ├── 0_Load_data.R                   # Data loading and preprocessing
│   ├── 01_comparative_analysis.R       # Coverage, vocabulary, richness, similarity,
│   │                                   #   specificity, enrichment (Figures 1–3)
│   ├── 02_taxonomic_classification.py  # LCA-based taxonomic origin of GO terms
│   └── 03_taxonomic_reliability.R      # Taxonomic reliability analysis (Table 1)
└── results/
    ├── go_terms_to_classify.tsv        # Unique GO terms exported for Python step
    ├── taxon_map_v2.tsv                # Taxonomic classification of each GO term
    ├── A2_vocabulary_comparison.tsv    # GO vocabulary overlap metrics
    ├── B1_semantic_similarity_BMA_Resnik.tsv
    ├── B4_GOEA_FANTASIA_exclusive.tsv  # Enriched GO terms (FANTASIA-only proteins)
    ├── B4_GOEA_Blast2Go_exclusive.tsv  # Enriched GO terms (Blast2Go-only proteins)
    ├── Table_1_taxonomic_classification.tsv
    ├── Figure_1_violin_GO_per_protein.png
    ├── Figure_2A_semantic_similarity.png
    ├── Figure_2B_specificity_barplot.png
    ├── Figure_3A_GOEA_FANTASIA.png
    └── Figure_3B_GOEA_Blast2Go.png
```

## How to reproduce the analysis

### Prerequisites

**R (v4.4.2)** with packages:
- Biostrings, GO.db (Bioconductor)
- ggplot2, eulerr, ontologyIndex, ontologySimilarity (CRAN)

**Python (v3.12)** with packages:
- pandas.
- ete3.
- tqdm.

### External data (not included due to size)

The following files are required but not included in this repository. Download them and place them in a `data/` folder:

| File | Source | Size |
|------|--------|------|
| `Cast.1_0.hap1.clean.fa` | [TreeGenes](https://treegenesdb.org/FTP/Genomes/.Cast/v1.0/)
| `go.obo` | [Gene Ontology](http://purl.obolibrary.org/obo/go.obo)
| `goa_uniprot_all.gaf.gz` | [UniProt-GOA FTP](https://ftp.ebi.ac.uk/pub/databases/GO/goa/) (downloaded 22/05/2026) 

### Execution order

```
1. Rscript scripts/01_comparative_analysis.R    # ~6-8 hours (BMA/Resnik is slow)
2. python  scripts/02_taxonomic_classification.py  # ~30-60 min (GAF parsing + LCA)
3. Rscript scripts/03_taxonomic_reliability.R    # ~1 min
```

Step 1 exports `results/go_terms_to_classify.tsv`, which is the input for step 2. Step 2 produces `results/taxon_map_v2.tsv`, which is the input for step 3. All output files in `results/` are provided so that each step can also be run independently.

## Funding

This work is funded by project ProyExcel_00351 (Proyectos de investigación de Excelencia, Junta de Andalucía). JVD is a Ramón y Cajal postdoctoral fellow supported by MCIN/AEI/10.13039/501100011033 (Ref. RYC2019-028188-I).
