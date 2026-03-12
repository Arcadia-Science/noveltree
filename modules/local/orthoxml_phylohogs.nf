process ORTHOXML_PHYLOHOGS {
    tag "$meta.og"

    cpus 1
    time { 4.h * task.attempt }
    memory {
        def n = (meta?.n_seq ?: 100) as long
        // NHX→OrthoXML DOM is O(n); pair export buffers scale ~O(n²)
        def base_gb = Math.max(2L, (long)(n / 50L) + 2L)
        def capped_gb = (int) Math.min(base_gb, 64L)
        def requested = capped_gb.GB * task.attempt
        def max_mem = params.max_memory as nextflow.util.MemoryUnit
        requested.compareTo(max_mem) > 0 ? max_mem : requested
    }

    container 'arcadiascience/orthoxml_phylohogs:1.0.0'

    input:
    tuple val(meta), path(nhx_file)

    output:
    tuple val(meta), path("${meta.og}.orthoxml")     , emit: orthoxml
    path "${meta.og}_orthologs.tsv"                   , emit: orthologs
    path "${meta.og}_paralogs.tsv"                    , emit: paralogs
    path "${meta.og}_xenologs.tsv"                    , emit: xenologs
    path "${meta.og}_hog_membership.tsv"              , emit: hog_membership
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def og = "${meta.og}"
    """
    # Convert GeneRax NHX reconciliation to OrthoXML
    orthoxml-tools from-nhx --species-encode nhx --infile ${nhx_file} --outfile ${og}.orthoxml

    # Export ortholog pairs (no header in output, two tab-separated protId columns)
    orthoxml-tools export-pairs ortho --id protId --infile ${og}.orthoxml --outfile ${og}_orthologs_raw.tsv
    echo -e "gene1\\tgene2\\tog" > ${og}_orthologs.tsv
    awk -v og="${og}" '{print \$0 "\\t" og}' ${og}_orthologs_raw.tsv >> ${og}_orthologs.tsv

    # Export paralog pairs
    orthoxml-tools export-pairs para --id protId --infile ${og}.orthoxml --outfile ${og}_paralogs_raw.tsv
    echo -e "gene1\\tgene2\\tog" > ${og}_paralogs.tsv
    awk -v og="${og}" '{print \$0 "\\t" og}' ${og}_paralogs_raw.tsv >> ${og}_paralogs.tsv

    # Extract xenolog pairs from NHX (H=Y transfer nodes)
    extract_xenologs.py ${nhx_file} ${og}  > ${og}_xenologs.tsv

    # Extract HOG membership from OrthoXML
    extract_hog_membership.py ${og}.orthoxml ${og} > ${og}_hog_membership.tsv

    # Clean up intermediate files
    rm -f ${og}_orthologs_raw.tsv ${og}_paralogs_raw.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orthoxml-tools: \$(orthoxml-tools --version 2>&1 | head -1 || echo "1.3.0")
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
