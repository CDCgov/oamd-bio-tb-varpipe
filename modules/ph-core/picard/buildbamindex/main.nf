process PICARD_BUILDBAMINDEX {

    tag "${meta.id}"
    label 'process_medium'

    container 'oamd-bio-biocontainers-picard:3.3.0_47ba639_v2'

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.bai"), emit: bai
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
        picard BuildBamIndex \\
            INPUT=${input} \\
            VALIDATION_STRINGENCY=${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CreateSequenceDictionary --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CreateSequenceDictionary --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
