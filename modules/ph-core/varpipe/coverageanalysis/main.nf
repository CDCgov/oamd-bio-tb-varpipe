process VARPIPE_COVERAGEANALYSIS {

    tag "${meta.id}"
    label 'process_single'

    container 'oamd-bio-biocontainers-python:3.12_80c558c_v3'

    input:
    tuple val(meta), path(coverage), path(unmapped), path(mapped)
    path(amp) // gene regions (single BED file)
    path(features) // feature regions (multiple BED files)

    output:
    tuple val(meta), path("*_stats.txt"), emit: stats
    tuple val(meta), path("*_target_region_coverage.txt"), emit: target_region_coverage, optional: true
    tuple val(meta), path("*_genome_region_coverage.txt"), emit: genome_region_coverage, optional: true
    tuple val(meta), path("*_qc_log.txt"), emit: qc_failed, optional: true
    tuple val(meta), path(".qc_passed.txt"), emit: qc_passed, optional: true
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'coverage_stats.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_stats.txt
    touch ${prefix}_target_region_coverage.txt
    touch ${prefix}_genome_region_coverage.txt
    touch ${prefix}_structural_variants.txt
    touch .qc_passed.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
