process LIST_CHROMS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::awk=5.1.0 coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.chroms.txt"), emit: chroms

    script:
    """
    awk '
    BEGIN{OFS="\\t"}
    /^>/{
        if (name != "") print name, len
        name = substr(\$1, 2)
        len = 0
        next
    }
    {
        gsub(/[[:space:]]/, "", \$0)
        len += length(\$0)
    }
    END{
        if (name != "") print name, len
    }' ${fasta} > ${meta.id}.chroms.txt
    """
}
