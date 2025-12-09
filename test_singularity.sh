#!/bin/bash
# Script to test noveltree pipeline with Singularity in Lima VM
# Usage: ./test_singularity.sh

set -e

echo "=== Testing noveltree with Singularity ==="
echo ""
echo "This script will:"
echo "1. Start the Lima VM (apptainer-x86) if not already running"
echo "2. Run the noveltree test pipeline using Singularity/Apptainer"
echo "3. Docker images will be automatically converted to Singularity format"
echo ""

# Start Lima VM if not running
echo "Checking Lima VM status..."
if ! limactl list | grep -q "noveltree-singularity.*Running"; then
    echo "Starting Lima VM..."
    limactl start noveltree-singularity
else
    echo "Lima VM already running"
fi

# Get the project directory
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo ""
echo "Running Nextflow test with Singularity profile..."
echo "Project directory: $PROJECT_DIR"
echo ""

# Run the pipeline in Lima
limactl shell noveltree-singularity bash <<EOF
set -e

# Copy project to a writable location in the VM
WORK_DIR="\$HOME/noveltree_test"
echo "Copying project to VM working directory: \$WORK_DIR"

# Clean and create working directory
rm -rf "\$WORK_DIR"
mkdir -p "\$WORK_DIR"

# Copy project files (excluding work and results directories)
rsync -av --exclude='work/' --exclude='results/' --exclude='tests/results*/' --exclude='.nextflow*' \\
    "$PROJECT_DIR/" "\$WORK_DIR/"

cd "\$WORK_DIR"

# Clean previous test results if they exist
if [ -d "tests/results_singularity" ]; then
    echo "Cleaning previous test results..."
    rm -rf tests/results_singularity
fi

# Run the pipeline with test and singularity profiles
echo ""
echo "Running pipeline..."
nextflow run main.nf \\
    -profile test,singularity \\
    --outdir tests/results_singularity \\
    -resume

echo ""
echo "=== Test completed successfully! ==="
echo ""
echo "Copying results back to host..."
rsync -av tests/results_singularity/ "$PROJECT_DIR/tests/results_singularity/"

echo ""
echo "Results are in: $PROJECT_DIR/tests/results_singularity/"
echo "Singularity images cached in: $PROJECT_DIR/tests/results_singularity/singularity_cache/"
EOF

echo ""
echo "=== All done! ==="
echo ""
echo "To run your full pipeline with Singularity, use:"
echo "  limactl shell noveltree-singularity"
echo "  Copy your project to VM: rsync -av $PROJECT_DIR/ ~/noveltree/ --exclude='work/' --exclude='results/' --exclude='.nextflow*'"
echo "  cd ~/noveltree"
echo "  nextflow run main.nf -profile singularity --input <your_input> --outdir <your_outdir>"
