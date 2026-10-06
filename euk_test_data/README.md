# Eukaryote test data

This directory contains two NovelTree samplesheets and six small, committed protein FASTAs. The local FASTAs contain 59–91 proteins per species drawn from 50 orthogroups; `selected_orthogroups.txt` records the selected families. The FASTAs are ready to use. No fixture-generation or external preprocessing script is required.

Run the commands below from the repository root. The `test` profile limits resources, supplies the local samplesheet by default, and processes only gene families found in every input species. See [the usage guide](../docs/usage.md#preparation) for the full samplesheet format and preprocessing options.

## Small local dataset

`samplesheet.csv` lists the six committed protein FASTAs. To run the quick test without optional preprocessing:

```bash
nextflow run . -profile test,docker --outdir tests/results
```

To run the same files through NovelTree's built-in preprocessing:

```bash
nextflow run . -profile test,docker \
    --input euk_test_data/samplesheet.csv \
    --preprocess true \
    --outdir tests/results_preprocess_local
```

The local samplesheet supplies no `filter_isoforms`, `transdecoder`, or `reference_proteome` flags, so they default to `no`. With `--preprocess true`, all six protein files undergo stop-codon and ambiguous-residue cleanup, the 50-aa minimum-length filter, and CD-HIT at 100% identity. This test does not exercise isoform filtering or TransDecoder. The later FASTA normalization converts U and O to X on every run, including when optional preprocessing is off.

## Mixed-source preprocessing test

`samplesheet_preprocess.csv` exercises additional paths and downloads full proteomes or transcripts from UniProt and NCBI. It requires internet access and is substantially slower than the local test:

```bash
nextflow run . -profile test,docker \
    --input euk_test_data/samplesheet_preprocess.csv \
    --preprocess true \
    --outdir tests/results_preprocess
```

| Species | Input | Path exercised with `--preprocess true` |
| --- | --- | --- |
| Homo sapiens | Local protein FASTA | Cleanup and length filter; `reference_proteome=yes` skips CD-HIT. |
| Danio rerio | Local protein FASTA | Header-based isoform filter, then CD-HIT at 100%. All 91 fixture headers lack a recognized gene group, so this checks that the filter runs; it does not test removal of an isoform. |
| Nannochloropsis sp. | NCBI accession `GCF_000240725.1` | Paired protein and GFF3 download; longest annotated isoform per gene; CD-HIT skipped so distinct genes can retain identical sequences. |
| Schizosaccharomyces pombe | Gzipped nucleotide URL | TransDecoder, header-based isoform filter, then CD-HIT at 97%. |
| Saccharomyces cerevisiae | UniProt proteome ID `UP000002311` | Proteome download and isoform filter; `reference_proteome=yes` skips CD-HIT. |
| Neurospora crassa | UniProt FASTA URL | Download, cleanup, and length filter; `reference_proteome=yes` skips CD-HIT. |

`--preprocess true` is required: per-species samplesheet flags do not turn preprocessing on by themselves. For an NCBI genome accession, NovelTree sets `filter_isoforms=yes` even if the samplesheet says `no`. The `reference_proteome` flag is set by the user, including for UniProt IDs. This samplesheet does not exercise a user-supplied `annotation_gff3`; see [the usage guide](../docs/usage.md#preprocessing-behavior-when---preprocess-true) for that route.

Both samplesheets include BUSCO lineage columns, but BUSCO runs only with `--busco true`. MCL inflation optimization is also optional; see [the usage guide](../docs/usage.md) for its parameters.
