#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Create a small eukaryote test dataset for noveltree pipeline debugging.
#
# Strategy: sample ~50 gene families from a previous OrthoFinder run that
# contain >=80% of our 6 target Opisthokont species, then extract only
# sequences belonging to those species. This preserves real orthology
# signal so OrthoFinder can recover meaningful gene families.
#
# Prerequisites:
#   - Previous noveltree results in results-noveltree-model-euks/
#   - Docker (for build_reference_chronogram)
#   - Python 3
#
# Usage:
#   ./euk_test_data/create_test_data.sh --email <ncbi_email>
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# --- CLI parsing ---
EMAIL=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --email) EMAIL="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ -z "$EMAIL" ]]; then
    echo "Usage: $0 --email <ncbi_email>"
    exit 1
fi

# --- Configuration ---
GENE_COUNT_TSV="${PROJECT_DIR}/results-noveltree-model-euks/orthofinder/complete_dataset/Results_Inflation_1.5/Orthogroups/Orthogroups.GeneCount.tsv"
OG_SEQ_DIR="${PROJECT_DIR}/results-noveltree-model-euks/orthofinder/complete_dataset/Results_Inflation_1.5/Orthogroup_Sequences"
OUT_DIR="${SCRIPT_DIR}"
PROTEOME_DIR="${OUT_DIR}/proteomes"
NUM_OGS=50
MIN_SPECIES_FRAC=0.80
SEED=42

# 6 target Opisthokont species (4 animals + 2 fungi)
SPECIES=(
    "Homo-sapiens"
    "Mus-musculus"
    "Danio-rerio"
    "Drosophila-melanogaster"
    "Saccharomyces-cerevisiae"
    "Neurospora-crassa"
)
NUM_SPECIES=${#SPECIES[@]}
MIN_SPECIES=$(python3 -c "import math; print(math.ceil(${NUM_SPECIES} * ${MIN_SPECIES_FRAC}))")

echo "=== Noveltree Eukaryote Test Data Creator ==="
echo "Target species: ${NUM_SPECIES}"
echo "Min species per OG: ${MIN_SPECIES} (${MIN_SPECIES_FRAC})"
echo "OGs to sample: ${NUM_OGS}"
echo ""

# --- Verify inputs exist ---
if [[ ! -f "$GENE_COUNT_TSV" ]]; then
    echo "ERROR: Gene count file not found: $GENE_COUNT_TSV"
    echo "Make sure results-noveltree-model-euks/ is present."
    exit 1
fi
if [[ ! -d "$OG_SEQ_DIR" ]]; then
    echo "ERROR: Orthogroup sequences dir not found: $OG_SEQ_DIR"
    exit 1
fi

# --- Step 1: Select orthogroups with >=80% of target species ---
echo "Step 1: Selecting orthogroups with >= ${MIN_SPECIES}/${NUM_SPECIES} target species..."

mkdir -p "$PROTEOME_DIR"

python3 - "$GENE_COUNT_TSV" "$OG_SEQ_DIR" "$PROTEOME_DIR" "$NUM_OGS" "$MIN_SPECIES" "$SEED" "${SPECIES[@]}" <<'PYTHON_SCRIPT'
import csv
import os
import random
import sys
from collections import defaultdict

gene_count_tsv = sys.argv[1]
og_seq_dir = sys.argv[2]
proteome_dir = sys.argv[3]
num_ogs = int(sys.argv[4])
min_species = int(sys.argv[5])
seed = int(sys.argv[6])
target_species = set(sys.argv[7:])

random.seed(seed)

# --- Read gene counts and find eligible OGs ---
print(f"  Reading {gene_count_tsv}...")
eligible_ogs = []

with open(gene_count_tsv) as f:
    reader = csv.DictReader(f, delimiter='\t')
    # Verify our target species are in the header
    header_species = set(reader.fieldnames) - {'Orthogroup', 'Total'}
    missing = target_species - header_species
    if missing:
        print(f"  WARNING: Species not in gene counts: {missing}")
        target_species -= missing

    for row in reader:
        og = row['Orthogroup']
        # Count how many target species have >= 1 gene in this OG
        spp_present = sum(
            1 for sp in target_species if int(row.get(sp, 0)) > 0
        )
        if spp_present >= min_species:
            # Also compute total genes across target species (prefer smaller OGs)
            total_genes = sum(int(row.get(sp, 0)) for sp in target_species)
            eligible_ogs.append((og, spp_present, total_genes))

print(f"  Found {len(eligible_ogs)} OGs with >= {min_species} target species")

if len(eligible_ogs) < num_ogs:
    print(f"  WARNING: Only {len(eligible_ogs)} eligible OGs (wanted {num_ogs})")
    num_ogs = len(eligible_ogs)

# Sort by species coverage (desc), then by total gene count (asc, prefer smaller families)
# to get well-represented but not huge OGs, then sample randomly from top candidates
eligible_ogs.sort(key=lambda x: (-x[1], x[2]))

# Filter to OGs with reasonable size (< 50 total genes across target species)
# to keep the test fast
reasonable_ogs = [og for og in eligible_ogs if og[2] <= 50]
if len(reasonable_ogs) < num_ogs:
    # Relax if needed
    reasonable_ogs = eligible_ogs

# Sample from the reasonable set
selected = random.sample(reasonable_ogs[:max(num_ogs * 5, len(reasonable_ogs))], num_ogs)
selected_og_names = [og[0] for og in selected]

print(f"  Selected {len(selected_og_names)} OGs")
for og, nspp, ngenes in sorted(selected, key=lambda x: x[0]):
    print(f"    {og}: {nspp} species, {ngenes} genes")

# --- Step 2: Extract sequences per species from selected OGs ---
print(f"\n  Extracting sequences for target species from selected OGs...")
species_sequences = defaultdict(list)  # species -> [(header, seq)]

for og_name in selected_og_names:
    og_fasta = os.path.join(og_seq_dir, f"{og_name}.fa")
    if not os.path.exists(og_fasta):
        print(f"  WARNING: {og_fasta} not found, skipping")
        continue

    # Parse FASTA
    with open(og_fasta) as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    seq = ''.join(seq_lines)
                    # Extract species from header: >Species-name_proteinID
                    sp = header.split('_', 1)[0].lstrip('>')
                    if sp in target_species:
                        species_sequences[sp].append((header, seq))
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        # Last record
        if header is not None:
            seq = ''.join(seq_lines)
            sp = header.split('_', 1)[0].lstrip('>')
            if sp in target_species:
                species_sequences[sp].append((header, seq))

# --- Step 3: Write per-species FASTAs ---
print(f"\n  Writing per-species FASTA files...")
for sp in sorted(target_species):
    seqs = species_sequences.get(sp, [])
    outpath = os.path.join(proteome_dir, f"{sp}.fasta")
    with open(outpath, 'w') as f:
        for header, seq in seqs:
            f.write(f"{header}\n")
            # Write sequence in 60-char lines
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')
    print(f"    {sp}: {len(seqs)} sequences -> {outpath}")

# Write selected OG names for reference
og_list_path = os.path.join(os.path.dirname(proteome_dir), "selected_orthogroups.txt")
with open(og_list_path, 'w') as f:
    for og in sorted(selected_og_names):
        f.write(og + '\n')
print(f"\n  OG list written to {og_list_path}")
PYTHON_SCRIPT

echo ""
echo "Step 2: Building reference chronogram..."

# Write species names (underscore format) for build_reference_chronogram.py
SPECIES_FILE=$(mktemp)
for sp in "${SPECIES[@]}"; do
    echo "${sp//-/_}" >> "$SPECIES_FILE"
done

# Run build_reference_chronogram.py via Docker
docker run --rm \
    -v "${PROJECT_DIR}/bin/build_reference_chronogram.py:/app/build_reference_chronogram.py:ro" \
    -v "${SPECIES_FILE}:/app/species_names.txt:ro" \
    -v "${OUT_DIR}:/app/output" \
    arcadiascience/build_reference_chronogram:1.0.0 \
    python3 /app/build_reference_chronogram.py \
        --species-names /app/species_names.txt \
        --output /app/output/reference_chronogram.nwk \
        --email "$EMAIL"

rm -f "$SPECIES_FILE"

if [[ -f "${OUT_DIR}/reference_chronogram.nwk" ]]; then
    echo "  Reference chronogram written to ${OUT_DIR}/reference_chronogram.nwk"
else
    echo "  WARNING: Reference chronogram was not created."
    echo "  You can create it manually later."
fi

echo ""
echo "Step 3: Generating samplesheet..."

# Samplesheet columns: species,file,taxonomy,shallow_db,broad_db,mode,uniprot,mcl_test
# taxonomy, shallow_db, broad_db use BUSCO lineages
# mode: proteins for all
# uniprot: true for all
# mcl_test: true for animals, false for fungi

cat > "${OUT_DIR}/samplesheet.csv" <<SAMPLESHEET
species,file,taxonomy,shallow_db,broad_db,mode,uniprot,mcl_test
Homo-sapiens,euk_test_data/proteomes/Homo-sapiens.fasta,Opisthokonta,primates_odb10,eukaryota_odb10,proteins,true,true
Mus-musculus,euk_test_data/proteomes/Mus-musculus.fasta,Opisthokonta,glires_odb10,eukaryota_odb10,proteins,true,true
Danio-rerio,euk_test_data/proteomes/Danio-rerio.fasta,Opisthokonta,actinopterygii_odb10,eukaryota_odb10,proteins,true,true
Drosophila-melanogaster,euk_test_data/proteomes/Drosophila-melanogaster.fasta,Opisthokonta,diptera_odb10,eukaryota_odb10,proteins,true,true
Saccharomyces-cerevisiae,euk_test_data/proteomes/Saccharomyces-cerevisiae.fasta,Opisthokonta,saccharomycetes_odb10,eukaryota_odb10,proteins,true,false
Neurospora-crassa,euk_test_data/proteomes/Neurospora-crassa.fasta,Opisthokonta,sordariomycetes_odb10,eukaryota_odb10,proteins,true,false
SAMPLESHEET

echo "  Samplesheet written to ${OUT_DIR}/samplesheet.csv"

echo ""
echo "=== Done ==="
echo ""
echo "Output:"
echo "  Proteomes:            ${PROTEOME_DIR}/"
echo "  Samplesheet:          ${OUT_DIR}/samplesheet.csv"
echo "  Reference chronogram: ${OUT_DIR}/reference_chronogram.nwk"
echo "  Selected OGs:         ${OUT_DIR}/selected_orthogroups.txt"
echo ""
echo "--- Example run commands ---"
echo ""
echo "Standard mode:"
echo "  nextflow run . -profile docker \\"
echo "    --input euk_test_data/samplesheet.csv \\"
echo "    --outdir euk_test_results \\"
echo "    --ncbi_email ${EMAIL}"
echo ""
echo "Zoogle mode:"
echo "  nextflow run . -profile zoogle,docker \\"
echo "    --input euk_test_data/samplesheet.csv \\"
echo "    --outdir euk_test_results \\"
echo "    --reference_time_tree euk_test_data/reference_chronogram.nwk \\"
echo "    --ncbi_email ${EMAIL}"
echo ""
echo "Dry-run (preview):"
echo "  nextflow run . -profile zoogle,docker \\"
echo "    --input euk_test_data/samplesheet.csv \\"
echo "    --outdir euk_test_results \\"
echo "    --reference_time_tree euk_test_data/reference_chronogram.nwk \\"
echo "    --ncbi_email ${EMAIL} \\"
echo "    -preview"
