# Protein language model versus homology-based functional annotation in non-model forest trees: a case study in *Castanea sativa* Mill.

Companion repository for:

> **Postigo-Luque, R. et al.** Protein language model versus homology-based functional annotation in non-model forest trees: a case study in *Castanea sativa* Mill.

This repository contains the functional annotations and analysis scripts used in the study comparing [FANTASIA](https://github.com/FANTASIA-Annotation/FANTASIA) (protein language model-based) with two homology-based methods, **Blast2GO** (OmicsBox v3.2.9) and **AHRD** (Automated Assignment of Human Readable Descriptions), on the *Castanea sativa* (European chestnut) haplotype 1 proteome.

## Repository structure

```
├── Annotations/
│   ├── Csath1_FANTASIA_TopGO.txt          # FANTASIA GO annotations (v1.0, ProtT5)
│   ├── GO_Hap1_JVD.txt                    # Blast2GO GO annotations (OmicsBox v3.2.9)
│   └── Casat.proteins.funct_ahrd.tsv      # AHRD GO annotations (DIAMOND + Swiss-Prot/TrEMBL)
├── data/
│   ├── counts.txt                         # featureCounts gene-level counts (RNA-seq use case)
│   └── Cast.1_0.hap1.clean_cdhit100.pep.clstr   # CD-HIT 100%-identity clusters
├── Scripts/
│   ├── 0_Load_data.R                      # Data loading and preprocessing (sourced by 1, 3, 4)
│   ├── 1_Comparative_analysis.R           # FANTASIA vs Blast2GO: coverage, vocabulary, richness,
│   │                                      #   semantic similarity + null model, specificity,
│   │                                      #   effect-size CIs and sensitivity  (Fig 1, 2, S1)
│   ├── 2_taxonomic_classification.py      # LCA-based taxonomic origin of GO terms
│   ├── 3_taxonomic_reliability.R          # FANTASIA vs Blast2GO taxonomic reliability (Table 1)
│   ├── 4_anomalous_classification.R       # Classification of anomalous GO terms (Table S1)
│   ├── 5_enrichment_RNAseq.R              # RNA-seq use case: DESeq2 + enrichment, 3 methods (Fig 3, S5)
│   └── 6_AHRD_comparison.R                # FANTASIA vs AHRD parallel comparison (Fig S2, S3, S4)
└── results/                              # created when running the scripts (see below)
```

## Scripts and the manuscript outputs they produce

| Script | Comparison | Outputs |
|--------|-----------|---------|
| `1_Comparative_analysis.R` | FANTASIA vs Blast2GO | Figure 1, Figure 2, Figure S1 |
| `3_taxonomic_reliability.R` | FANTASIA vs Blast2GO | Table 1 |
| `4_anomalous_classification.R` | FANTASIA vs Blast2GO | Table S1 |
| `5_enrichment_RNAseq.R` | FANTASIA vs Blast2GO vs AHRD | Figure 3, Figure S5 |
| `6_AHRD_comparison.R` | FANTASIA vs AHRD | Figure S2, Figure S3, Figure S4 |

## How to reproduce the analysis

### Prerequisites

**R (v4.4.2)** with packages:
- Biostrings, GO.db, DESeq2, clusterProfiler (Bioconductor)
- ggplot2, eulerr, ontologyIndex, ontologySimilarity, boot, patchwork (CRAN)

**Python (v3.12)** with packages:
- pandas
- ete3
- tqdm

Before running, set the working directory to this repository folder, and place the external data files below in a `data/` folder.

### External data (not included due to size)

| File | Source | Used by |
|------|--------|---------|
| `Cast.1_0.hap1.clean.fa` | [TreeGenes](https://treegenesdb.org/FTP/Genomes/.Cast/v1.0/) (filtered proteome; see Methods) | 0, 6 |
| `go.obo` | [Gene Ontology](http://purl.obolibrary.org/obo/go.obo) | 0, 3, 6 |
| `goa_uniprot_all.gaf.gz` | [UniProt-GOA FTP](https://ftp.ebi.ac.uk/pub/databases/GO/goa/) (release 2025_01, downloaded 22/05/2026) | 2 |

### Execution order

**Main comparison (FANTASIA vs Blast2GO) and taxonomic reliability**
```
1. Rscript Scripts/1_Comparative_analysis.R      # sources 0_Load_data.R   (~6-8 h; BMA/Resnik + null model)
2. python  Scripts/2_taxonomic_classification.py # ~30-60 min (GAF parsing + LCA)
3. Rscript Scripts/3_taxonomic_reliability.R     # Table 1  (~1 min)
4. Rscript Scripts/4_anomalous_classification.R  # Table S1 (sources 0 and 3)
```
Step 1 exports `results/go_terms_to_classify.tsv` (input for step 2); step 2 produces
`results/taxon_map_v2.tsv` (input for steps 3 and 4).

**RNA-seq use case (enrichment, three methods)**
```
5. Rscript Scripts/5_enrichment_RNAseq.R         # Figures 3 and S5
```

**Supplementary comparison (FANTASIA vs AHRD)**
```
6. Rscript Scripts/6_AHRD_comparison.R           # Figures S2, S3 and S4
```

`0_Load_data.R` is sourced by scripts 1, 3 and 4; it is not run on its own.
All output files are written to `results/` when the scripts are run.

### RNA-seq use case: upstream pipeline

The enrichment use case (script 5) uses a public RNA-seq experiment of *C. sativa* roots,
control vs inoculated with the oomycete *Phytophthora cinnamomi* (2 h post-inoculation;
3 + 3 replicates) — Fernandes et al. 2024 (*Front. Plant Sci.*), **BioProject PRJNA1108966**.
The provided `data/counts.txt` was generated as:

```
SRA (prefetch / fasterq-dump)
  -> fastp (v1.4.0)                  quality trimming
  -> HISAT2 (v2.2.1)                 alignment to the C. sativa hap1 genome
  -> featureCounts (Subread v2.0.6)  gene-level counts  ->  data/counts.txt
```

Script 5 then runs DESeq2 (DEGs at padj < 0.05 and |log2FC| >= 1) and GO enrichment
(`clusterProfiler::enricher`, hypergeometric test, Benjamini-Hochberg correction,
FDR < 0.05) with the same set of DEGs and the same background (all expressed genes)
for the three annotation methods, so that the annotation is the only variable.

## Data availability

The functional annotations generated by FANTASIA, Blast2GO and AHRD, together with the
analysis scripts (R and Python), are provided in this repository. The *Castanea sativa*
haplotype 1 genome assembly and protein sequences are available at
[TreeGenes](https://treegenesdb.org/FTP/Genomes/.Cast/v1.0/).

## Funding

This work is funded by project ProyExcel_00351 (Proyectos de investigación de Excelencia,
Junta de Andalucía). JVD is a Ramón y Cajal postdoctoral fellow supported by
MCIN/AEI/10.13039/501100011033 (Ref. RYC2019-028188-I).
