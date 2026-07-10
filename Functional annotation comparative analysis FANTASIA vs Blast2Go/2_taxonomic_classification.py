# ============================================================
# 02_taxonomic_classification.py
# Taxonomic classification of GO terms via LCA (Lowest Common Ancestor)
# Postigo-Luque et al. — FANTASIA vs Blast2Go in C. sativa
#
# For each GO term, identifies all organisms annotated with it
# in UniProt-GOA, computes the LCA using NCBI Taxonomy, and
# outputs the taxonomic origin of each term.
#
# Input:  results/go_terms_to_classify.tsv  (from 01_comparative_analysis.R)
#         data/goa_uniprot_all.gaf.gz       (UniProt-GOA, downloaded 22/05/2026)
# Output: results/taxon_map_v2.tsv          (input for 03_taxonomic_reliability.R)
# ============================================================

import pandas as pd
import gzip
from ete3 import NCBITaxa
from tqdm import tqdm
from collections import defaultdict

# --- File paths (modify these to match your local setup) ---
GAF_FILE    = "data/goa_uniprot_all.gaf.gz"
GO_TERMS    = "results/go_terms_to_classify.tsv"
OUTPUT_FILE = "results/taxon_map_v2.tsv"

# ============================================================
# Step 1: Load GO terms to classify
# ============================================================
print("Loading GO terms to classify...")
go_terms_df  = pd.read_csv(GO_TERMS, sep="\t")
go_terms_set = set(go_terms_df["GO_term"].str.strip())
print(f"GO terms to classify: {len(go_terms_set)}")

# ============================================================
# Step 2: Parse GAF file and extract TaxIDs per GO term
# ============================================================
print("Reading GAF file (this may take 20-30 min)...")
go_to_taxids = defaultdict(set)

with gzip.open(GAF_FILE, "rt", encoding="utf-8") as f:
    for line in tqdm(f, desc="Processing GAF"):
        if line.startswith("!"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 13:
            continue
        go_term = parts[4].strip()
        taxon   = parts[12].strip()
        if go_term not in go_terms_set:
            continue
        if taxon.startswith("taxon:"):
            try:
                taxid = int(taxon.replace("taxon:", "").split("|")[0])
                go_to_taxids[go_term].add(taxid)
            except ValueError:
                continue

print(f"GO terms with TaxIDs found: {len(go_to_taxids)}")

# ============================================================
# Step 3: Compute LCA for each GO term
# ============================================================
print("Loading NCBI Taxonomy...")
ncbi = NCBITaxa()

print("Computing LCA for each GO term...")
results = []

for go_term in tqdm(go_terms_set, desc="Computing LCA"):
    taxids = go_to_taxids.get(go_term, set())

    if not taxids:
        results.append({"GO_term": go_term, "taxon": "Unknown", "taxon_id": None})
        continue

    try:
        taxids_list = list(taxids)

        if len(taxids_list) == 1:
            lca_id = taxids_list[0]
        else:
            tree   = ncbi.get_topology(taxids_list, intermediate_nodes=False)
            lca_id = int(tree.name)

        names    = ncbi.get_taxid_translator([lca_id])
        lca_name = names.get(lca_id, "Unknown")

        lineage   = ncbi.get_lineage(lca_id)
        lin_names = ncbi.get_taxid_translator(lineage)
        lin_list  = [lin_names.get(t, "") for t in lineage]

        results.append({
            "GO_term"  : go_term,
            "taxon"    : lca_name,
            "taxon_id" : lca_id,
            "lineage"  : " | ".join(lin_list)
        })

    except Exception as e:
        results.append({"GO_term": go_term, "taxon": "Unknown",
                        "taxon_id": None, "lineage": str(e)})

# ============================================================
# Step 4: Save results
# ============================================================
print(f"Saving {OUTPUT_FILE}...")
df_out = pd.DataFrame(results)
df_out.to_csv(OUTPUT_FILE, sep="\t", index=False)
print(f"Saved: {OUTPUT_FILE}")

print(f"\nTop 30 unique taxa:")
print(df_out["taxon"].value_counts().head(30))
print("COMPLETE")
