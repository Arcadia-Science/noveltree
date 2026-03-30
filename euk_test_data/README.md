# Eukaryote Test Data

This directory contains test samplesheets and small proteome files for testing the NovelTree pipeline.

## Samplesheet `input_data` Column Formats

The `input_data` column in the samplesheet accepts four formats:

| Format | Example | Description |
|---|---|---|
| Local path | `data/Species.fasta` | Path to a local FASTA file |
| URL | `https://rest.uniprot.org/uniprotkb/stream?query=(proteome:UP000005640)&format=fasta` | Direct download URL (http/https/ftp/s3) |
| UniProt proteome ID | `UP000005640` | Auto-downloaded via UniProt REST API |
| NCBI RefSeq accession | `GCF_000001405.40` | Auto-downloaded via NCBI `datasets` CLI |

### Auto-set flags by source type

| Source type | `filter_isoforms` | `reference_proteome` | `transdecoder` |
|---|---|---|---|
| UniProt reference (`UP*`) | User-set | `yes` | `no` |
| UniProt representative (`UP*`) | User-set | `no` | `no` |
| NCBI RefSeq (`GCF_*`) | **Auto-set to `yes`** | `no` | `no` |
| URL | User-set | User-set | User-set |
| Local file | User-set | User-set | User-set |

NCBI RefSeq proteomes always include isoforms, so `filter_isoforms` is automatically overridden to `yes` regardless of what is specified in the samplesheet.

## Samplesheets

### `samplesheet.csv` — Default test (fast, no internet required)

Uses pre-subsetted local proteome files (~100 sequences each). Suitable for CI and quick local testing. Does not exercise preprocessing.

```bash
nextflow run . -profile test,docker --outdir tests/results
```

| Species | Source | Notes |
|---|---|---|
| Homo sapiens | Local | has_uniprot_ids=yes, include_in_mcl_test=yes |
| Mus musculus | Local | has_uniprot_ids=yes, include_in_mcl_test=yes |
| Danio rerio | Local | has_uniprot_ids=yes, include_in_mcl_test=yes |
| Drosophila melanogaster | Local | has_uniprot_ids=yes, include_in_mcl_test=yes |
| Saccharomyces cerevisiae | Local | has_uniprot_ids=yes, include_in_mcl_test=no |
| Neurospora crassa | Local | has_uniprot_ids=yes, include_in_mcl_test=no |

### `samplesheet_preprocess.csv` — Preprocessing test (requires internet)

Exercises all preprocessing code paths: TransDecoder, isoform filtering, CD-HIT logic, URL downloads, accession-based downloads, and gzipped input handling. Downloads full proteomes from UniProt and NCBI, so runs are slower and require network access.

```bash
nextflow run . -profile test,docker --preprocess true \
    --input euk_test_data/samplesheet_preprocess.csv \
    --outdir tests/results_preprocess
```

| Species | Source | transdecoder | filter_isoforms | reference_proteome | Preprocessing path |
|---|---|---|---|---|---|
| Homo sapiens | Local file | no | no | yes | Skip CD-HIT (UniProt reference) |
| Danio rerio | Local file | no | yes | no | Isoform filter -> CD-HIT 100% (exact dedup) |
| Nannochloropsis sp. | NCBI RefSeq accession | no | **auto: yes** | no | NCBI download -> isoform filter -> CD-HIT 100% |
| S. pombe | URL (nucleotide) | yes | **auto: yes** | no | URL download -> TransDecoder -> isoform filter -> CD-HIT 97% |
| S. cerevisiae | UniProt proteome ID | no | yes | yes | UniProt download -> isoform filter -> skip CD-HIT |
| N. crassa | UniProt REST URL | no | no | yes | URL download -> skip CD-HIT |

#### What each species tests

- **Homo sapiens**: Local file with `reference_proteome=yes` — verifies CD-HIT is skipped for UniProt reference proteomes.
- **Danio rerio**: Local file with `filter_isoforms=yes` — verifies isoform filtering runs, then CD-HIT at 100% (exact duplicate removal only).
- **Nannochloropsis sp.**: NCBI RefSeq accession (`GCF_000240725.1`) — verifies NCBI datasets download and auto-isoform filtering.
- **S. pombe**: URL to nucleotide FASTA with `transdecoder=yes` — verifies URL download, TransDecoder ORF prediction, auto-isoform filtering, and CD-HIT at 97%.
- **S. cerevisiae**: UniProt proteome ID (`UP000002311`) with `filter_isoforms=yes` + `reference_proteome=yes` — verifies UniProt REST download, isoform filtering, and CD-HIT skip.
- **N. crassa**: UniProt REST API URL with `reference_proteome=yes` — verifies URL download and CD-HIT skip with no other preprocessing.

All species also go through stop codon removal, rare amino acid handling (U->C, J/B/Z->X), and minimum length filtering (default 50 aa).

## Opt-In Features

### MCL Inflation Testing

By default, the pipeline uses a fixed MCL inflation value of `1.5` without testing alternatives. To enable inflation optimization (tests multiple values using COGEQC + InterPro annotations), pass `--test_run_mcl true` along with multiple inflation values:

```bash
nextflow run . -profile test,docker --outdir tests/results \
    --test_run_mcl true --mcl_inflation '1.1,1.3,1.5,2.0,3.0'
```

Note: `--test_run_mcl true` requires at least two `--mcl_inflation` values; the pipeline will exit with an error otherwise.

### BUSCO

BUSCO quality assessment is off by default. Enable it with `--busco true` and include `busco_shallow`/`busco_broad` columns in your samplesheet:

```bash
nextflow run . -profile test,docker --outdir tests/results --busco true
```

This runs both shallow and broad taxonomic scale BUSCO analyses. BUSCO results are not used by downstream modules.
