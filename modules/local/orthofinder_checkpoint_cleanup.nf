process CLEAN_ORTHOFINDER_CHECKPOINT {
    tag "Remove completed OrthoFinder MCL checkpoint"
    label 'process_single'
    label 'error_ignore'

    container 'arcadiascience/orthofinder_2.5.4:1.0.0'

    input:
    path checkpoint_receipt

    output:
    path "orthofinder_checkpoint_cleanup.tsv", emit: cleanup_record, optional: true

    script:
    """
    orthofinder_checkpoint.py cleanup --receipt ${checkpoint_receipt}
    printf 'status\tcheckpoint_receipt\nremoved\t%s\n' '${checkpoint_receipt}' \
        > orthofinder_checkpoint_cleanup.tsv
    """
}
