# Singularity Support for noveltree

This pipeline now supports Singularity/Apptainer for HPC environments!

## What was configured

1. **Singularity profile** in [nextflow.config](nextflow.config#L135-L144) with:
   - Auto-conversion from Docker images (`singularity.pullDockerContainer = true`)
   - Automatic cache directory for Singularity images
   - Auto-mounting of filesystems

2. **Lima VM setup** for local testing (since Singularity doesn't run on macOS natively)
   - VM: `apptainer-x86` with Apptainer 1.4.4
   - Nextflow 25.10.2 installed
   - Java 21 runtime

## Testing Locally (macOS)

### Quick Test
Run the test script:
```bash
./test_singularity.sh
```

This will run the pipeline's test suite with Singularity in the Lima VM.

## Running on HPC

On your HPC system, simply use the singularity profile:

```bash
nextflow run Arcadia-Science/noveltree \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results
```

### HPC-Specific Configurations

You may want to customize the Singularity cache location for your HPC:

```bash
nextflow run Arcadia-Science/noveltree \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results \
    -c custom.config
```

Where `custom.config` contains:
```groovy
singularity {
    cacheDir = '/path/to/shared/singularity/cache'
}
```

## How It Works

- **Docker → Singularity Auto-conversion**: When you use `-profile singularity`, Nextflow automatically:
  1. Pulls Docker images from the specified registries
  2. Converts them to Singularity `.sif` format
  3. Caches them for reuse (in `singularity_cache/` by default)
