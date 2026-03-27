process BUILD_REFERENCE_CHRONOGRAM {
    tag "Build reference chronogram from TimeTree.org"
    label 'process_single'

    container 'arcadiascience/build_reference_chronogram:1.0.0'

    storeDir "${params.outdir}/species_trees/reference_chronogram"

    input:
    path species_names    // one species per line, Genus_species format
    val ncbi_email

    output:
    path "ref_chronogram.nwk"    , emit: chronogram
    path "taxid_cache.json"      , emit: taxid_cache
    path "timetree_cache.json"   , emit: timetree_cache
    path "dropped_species.txt"   , emit: dropped_species

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Ensure HOME is writable (Singularity mounts host HOME read-only)
    export HOME=\$PWD

    build_reference_chronogram.py \\
        --species-names ${species_names} \\
        --output ref_chronogram.nwk \\
        --email ${ncbi_email} \\
        --workers 4 \\
        --cache-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
