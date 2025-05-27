process VARPIPE_SNPEFF {
    tag "$meta.id"
    label 'process_medium'

    container 'oamd-bio-snpeff:4.3r_a5da4a7_v1'

    input:
    tuple val(meta), path(vcf)
    val   db
    tuple val(meta2), path(cache)

    output:
    tuple val(meta), path("*.ann.vcf"),       emit: vcf
    tuple val(meta), path("*.csv"),       emit: report
    tuple val(meta), path("*.html"),      emit: summary_html
    tuple val(meta), path("*.genes.txt"), emit: genes_txt
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = 6144
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cache_command = cache ? "-dataDir \${PWD}/${cache}" : ""
    """
    zcat ${vcf} | sed 's/NC_000962.3/NC_000962/g' | gzip > tmp.vcf.gz

    java -jar /usr/src/snpEff/snpEff/snpEff.jar \\
        $db \\
        $args \\
        -csvStats ${prefix}.csv \\
        $cache_command \\
        tmp.vcf.gz \\
        > ${prefix}.ann.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(java -jar /usr/src/snpEff/snpEff/snpEff.jar -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ann.vcf
    touch ${prefix}.csv
    touch ${prefix}.html
    touch ${prefix}.genes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(java -jar /usr/src/snpEff/snpEff/snpEff.jar -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """

}
