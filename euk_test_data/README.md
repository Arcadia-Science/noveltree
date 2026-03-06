# Eukaryote Test Data

This directory contains test samplesheets and small proteome files for testing the NovelTree pipeline.

## Samplesheets

### `samplesheet.csv` — Default test (fast, no internet required)

Uses pre-subsetted local proteome files (~100 sequences each). Suitable for CI and quick local testing. Does not exercise preprocessing.

```bash
nextflow run . -profile test,docker --outdir tests/results
```

| Species | Source | Notes |
|---|---|---|
| Homo sapiens | Local | uniprot=true, mcl_test=true |
| Mus musculus | Local | uniprot=true, mcl_test=true |
| Danio rerio | Local | uniprot=true, mcl_test=true |
| Drosophila melanogaster | Local | uniprot=true, mcl_test=true |
| Saccharomyces cerevisiae | Local | uniprot=true, mcl_test=false |
| Neurospora crassa | Local | uniprot=true, mcl_test=false |

### `samplesheet_preprocess.csv` — Preprocessing test (requires internet)

Exercises all preprocessing code paths: TransDecoder, isoform filtering, CD-HIT logic, URL downloads, and gzipped input handling. Downloads full proteomes from UniProt and NCBI, so runs are slower and require network access.

```bash
nextflow run . -profile test,docker --preprocess true \
    --input euk_test_data/samplesheet_preprocess.csv \
    --outdir tests/results_preprocess
```

| Species | Source | transdecoder | isoform | reference | Preprocessing path |
|---|---|---|---|---|---|
| Homo sapiens | Local file | no | no | yes | Skip CD-HIT (UniProt reference) |
| Mus musculus | Local file | no | no | yes | Skip CD-HIT (UniProt reference) |
| Danio rerio | Local file | no | yes | no | Isoform filter → CD-HIT 100% (exact dedup) |
| S. pombe | NCBI FTP URL (.fna.gz) | yes | yes | no | TransDecoder → isoform filter → CD-HIT 97% |
| S. cerevisiae | UniProt REST URL | no | yes | yes | Isoform filter → skip CD-HIT |
| N. crassa | UniProt REST URL | no | no | yes | URL download → skip CD-HIT |

#### What each species tests

- **Homo sapiens, Mus musculus**: Local files with `reference=yes` — verifies CD-HIT is skipped for UniProt reference proteomes.
- **Danio rerio**: Local file with `isoform=yes` — verifies isoform filtering runs, then CD-HIT at 100% (exact duplicate removal only).
- **S. pombe**: NCBI FTP download of gzipped CDS nucleotide file (`.fna.gz`) with `transdecoder=yes` — verifies gzip decompression, TransDecoder ORF prediction, isoform filtering (auto-triggered by transdecoder), and CD-HIT at 97% (collapse assembly/prediction artifacts).
- **S. cerevisiae**: UniProt REST API download with `isoform=yes` + `reference=yes` — verifies URL download, isoform filtering, and CD-HIT skip.
- **N. crassa**: UniProt REST API download with `reference=yes` — verifies URL download and CD-HIT skip with no other preprocessing.

All six species also go through stop codon removal, rare amino acid handling (U→C, J/B/Z→X), and minimum length filtering (default 50 aa).
