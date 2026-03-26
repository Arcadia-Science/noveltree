# Eukaryote Test Data

This directory contains test samplesheets and small proteome files for testing the NovelTree pipeline.

## Samplesheet `file` Column Formats

The `file` column in the samplesheet accepts four formats:

| Format | Example | Description |
|---|---|---|
| Local path | `data/Species.fasta` | Path to a local FASTA file |
| URL | `https://rest.uniprot.org/uniprotkb/stream?query=(proteome:UP000005640)&format=fasta` | Direct download URL (http/https/ftp/s3) |
| UniProt proteome ID | `UP000005640` | Auto-downloaded via UniProt REST API |
| NCBI RefSeq accession | `GCF_000001405.40` | Auto-downloaded via NCBI `datasets` CLI |

### Auto-set flags by source type

| Source type | `isoform` | `reference` | `transdecoder` |
|---|---|---|---|
| UniProt reference (`UP*`) | User-set | `yes` | `no` |
| UniProt representative (`UP*`) | User-set | `no` | `no` |
| NCBI RefSeq (`GCF_*`) | **Auto-set to `yes`** | `no` | `no` |
| URL | User-set | User-set | User-set |
| Local file | User-set | User-set | User-set |

NCBI RefSeq proteomes always include isoforms, so `isoform` is automatically overridden to `yes` regardless of what is specified in the samplesheet.

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

Exercises all preprocessing code paths: TransDecoder, isoform filtering, CD-HIT logic, URL downloads, accession-based downloads, and gzipped input handling. Downloads full proteomes from UniProt and NCBI, so runs are slower and require network access.

```bash
nextflow run . -profile test,docker --preprocess true \
    --input euk_test_data/samplesheet_preprocess.csv \
    --outdir tests/results_preprocess
```

| Species | Source | transdecoder | isoform | reference | Preprocessing path |
|---|---|---|---|---|---|
| Homo sapiens | Local file | no | no | yes | Skip CD-HIT (UniProt reference) |
| Mus musculus | Local file | no | no | yes | Skip CD-HIT (UniProt reference) |
| Danio rerio | Local file | no | yes | no | Isoform filter -> CD-HIT 100% (exact dedup) |
| S. pombe | NCBI RefSeq accession | yes | **auto: yes** | no | NCBI download -> TransDecoder -> isoform filter -> CD-HIT 97% |
| S. cerevisiae | UniProt proteome ID | no | yes | yes | UniProt download -> isoform filter -> skip CD-HIT |
| N. crassa | UniProt REST URL | no | no | yes | URL download -> skip CD-HIT |

#### What each species tests

- **Homo sapiens, Mus musculus**: Local files with `reference=yes` — verifies CD-HIT is skipped for UniProt reference proteomes.
- **Danio rerio**: Local file with `isoform=yes` — verifies isoform filtering runs, then CD-HIT at 100% (exact duplicate removal only).
- **S. pombe**: NCBI RefSeq accession (`GCF_000002945.2`) with `transdecoder=yes` — verifies NCBI datasets download, TransDecoder ORF prediction, auto-isoform filtering (triggered by both RefSeq auto-override and transdecoder), and CD-HIT at 97%.
- **S. cerevisiae**: UniProt proteome ID (`UP000002311`) with `isoform=yes` + `reference=yes` — verifies UniProt REST download, isoform filtering, and CD-HIT skip.
- **N. crassa**: UniProt REST API URL with `reference=yes` — verifies URL download and CD-HIT skip with no other preprocessing.

All six species also go through stop codon removal, rare amino acid handling (U->C, J/B/Z->X), and minimum length filtering (default 50 aa).

## Opt-In Features

### MCL Inflation Testing

By default, the pipeline uses a fixed MCL inflation value of `1.5` without testing alternatives. To enable inflation optimization (tests multiple values using COGEQC + InterPro annotations), pass `--test_run_mcl true` along with multiple inflation values:

```bash
nextflow run . -profile test,docker --outdir tests/results \
    --test_run_mcl true --mcl_inflation '1.1,1.3,1.5,2.0,3.0'
```

Note: `--test_run_mcl true` requires at least two `--mcl_inflation` values; the pipeline will exit with an error otherwise.

### BUSCO

BUSCO quality assessment is off by default. Enable it with `--busco true`:

```bash
nextflow run . -profile test,docker --outdir tests/results --busco true
```

This runs both shallow and broad taxonomic scale BUSCO analyses. BUSCO results are not used by downstream modules.
