//
// varpipe core subworkflow
//

include { TRIMMOMATIC                                                           } from '../../../modules/ph-core/trimmomatic/main'
include { VARPIPE_PROCESSBAM                                                    } from '../varpipe_processbam/main'
include { VARPIPE_COVERAGEANALYSIS                                              } from '../../../modules/ph-core/varpipe/coverageanalysis/main'
include { VARPIPE_STRUCTURALVARIANTS                                            } from '../../../modules/ph-core/varpipe/structuralvariants/main'
include { VARPIPE_TARGETEDVARIANTS                                              } from '../varpipe_targetedvariants/main'
include { VARPIPE_ALLVARIANTS                                                   } from '../varpipe_allvariants/main'
include { VARPIPE_CREATEANNOTATIONS as VARPIPE_CREATEANNOTATIONS_TARGETED       } from '../../../modules/ph-core/varpipe/createannotations/main'
include { VARPIPE_CREATEANNOTATIONS as VARPIPE_CREATEANNOTATIONS_ALL            } from '../../../modules/ph-core/varpipe/createannotations/main'
include { VARPIPE_PARSEANNOTATIONS as VARPIPE_PARSEANNOTATIONS_TARGETED         } from '../../../modules/ph-core/varpipe/parseannotations/main'
include { VARPIPE_PARSEANNOTATIONS as VARPIPE_PARSEANNOTATIONS_ALL              } from '../../../modules/ph-core/varpipe/parseannotations/main'
include { VARPIPE_LINEAGE                                                       } from '../../../modules/ph-core/varpipe/lineage/main'
include { VARPIPE_REPORT                                                        } from '../../../modules/ph-core/varpipe/report/main'
include { VARPIPE_INTERPRETATION                                                } from '../../../modules/ph-core/varpipe/interpretation/main'
include { VARPIPE_PDF                                                           } from '../../../modules/ph-core/varpipe/pdf/main'
include { VARPIPE_TAR                                                           } from '../../../modules/ph-core/varpipe/tar/main'

workflow VARPIPE_CORE {
    take:
    ch_cleaned_reads // channel (mandatory): [ val(meta), path(reads) ]
    datestamp        // string (mandatory): timestamp for output identification
    ref_fasta        // channel (mandatory): [ path(params.ref_fasta) ]
    bedlist_amp      // channel (mandatory): [ path(params.bedlist_amp) ]
    bedlist_one      // channel (mandatory): [ path(params.bedlist_one) ]
    bedlist_two      // channel (mandatory): [ path(params.bedlist_two) ]
    bedstruct        // channel (mandatory): [ path(params.bedstruct) ]
    mutationloci     // channel (mandatory): [ path(params.mutationloci) ]
    intervals        // channel (mandatory): [ path(params.intervals) ]
    reported         // channel (mandatory): [ path(params.reported) ]
    lineages         // channel (mandatory): [ path(params.lineages) ]
    snpeff_db        // channel (mandatory): [ path(params.snpeff_db) ]

    main:
    ch_versions = Channel.empty()

    // QC TRIMMOMATIC
    TRIMMOMATIC(
        ch_cleaned_reads
    )
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)

    ch_ref_fasta = Channel.fromPath(ref_fasta)
        .map { fasta -> [ [ id: 'reference' ], fasta] }
        .collect()

    VARPIPE_PROCESSBAM (
        TRIMMOMATIC.out.trimmed_reads,
        ch_ref_fasta
    )
    ch_versions = ch_versions.mix(VARPIPE_PROCESSBAM.out.versions)

    ch_coverageanalysis_input = VARPIPE_PROCESSBAM.out.coverage
        .join(VARPIPE_PROCESSBAM.out.unmapped)
        .join(VARPIPE_PROCESSBAM.out.mapped)

    VARPIPE_COVERAGEANALYSIS(
        ch_coverageanalysis_input,
        bedlist_amp,
        [bedlist_one, bedlist_two],
    )
    ch_versions = ch_versions.mix(VARPIPE_COVERAGEANALYSIS.out.versions)

    VARPIPE_STRUCTURALVARIANTS(
        VARPIPE_COVERAGEANALYSIS.out.qc_passed.join(ch_coverageanalysis_input),
        VARPIPE_COVERAGEANALYSIS.out.genome_region_coverage,
        bedstruct
    )
    ch_versions = ch_versions.mix(VARPIPE_STRUCTURALVARIANTS.out.versions)

    // CALLERS
    ch_gatk_targeted_input = VARPIPE_COVERAGEANALYSIS.out.qc_passed
        .map { meta1, _qc -> [meta1] }
        .join(VARPIPE_PROCESSBAM.out.bam)
        .join(VARPIPE_PROCESSBAM.out.bai)
        .map { meta1, bam, bai ->
            return [meta1, bam, bai, intervals]
        }

    VARPIPE_TARGETEDVARIANTS(
        ch_gatk_targeted_input,
        VARPIPE_PROCESSBAM.out.fai.collect(),
        VARPIPE_PROCESSBAM.out.reference_dict.collect(),
        ch_ref_fasta,
        snpeff_db
    )
    ch_versions = ch_versions.mix(VARPIPE_TARGETEDVARIANTS.out.versions)

    ch_gatk_whole_genome_input = VARPIPE_COVERAGEANALYSIS.out.qc_passed
        .map { meta1, _qc -> [meta1] }
        .join(VARPIPE_PROCESSBAM.out.bam)
        .join(VARPIPE_PROCESSBAM.out.bai)
        .map { meta1, bam, bai ->
            return [meta1, bam, bai, []]
        }

    VARPIPE_ALLVARIANTS(
        ch_gatk_whole_genome_input,
        VARPIPE_PROCESSBAM.out.fai.collect(),
        VARPIPE_PROCESSBAM.out.reference_dict.collect(),
        ch_ref_fasta,
        snpeff_db
    )
    ch_versions = ch_versions.mix(VARPIPE_ALLVARIANTS.out.versions)

    VARPIPE_CREATEANNOTATIONS_TARGETED(
        VARPIPE_TARGETEDVARIANTS.out.ann
    )
    ch_versions = ch_versions.mix(VARPIPE_CREATEANNOTATIONS_TARGETED.out.versions.first())

    VARPIPE_CREATEANNOTATIONS_ALL(
        VARPIPE_ALLVARIANTS.out.ann
    )
    ch_versions = ch_versions.mix(VARPIPE_CREATEANNOTATIONS_ALL.out.versions.first())


    VARPIPE_PARSEANNOTATIONS_TARGETED(
        VARPIPE_CREATEANNOTATIONS_TARGETED.out.drLociAnnotation,
        mutationloci
    )
    ch_versions = ch_versions.mix(VARPIPE_PARSEANNOTATIONS_TARGETED.out.versions.first())

    VARPIPE_PARSEANNOTATIONS_ALL(
        VARPIPE_CREATEANNOTATIONS_ALL.out.fullAnnotation,
        mutationloci
    )
    ch_versions = ch_versions.mix(VARPIPE_PARSEANNOTATIONS_ALL.out.versions.first())

    // RUN VARPIPE_LINEAGE ANALYSIS
    VARPIPE_LINEAGE(
        VARPIPE_PARSEANNOTATIONS_ALL.out.fullFinalAnnotation,
        lineages
    )
    ch_versions = ch_versions.mix(VARPIPE_LINEAGE.out.versions.first())

    VARPIPE_COVERAGEANALYSIS.out.stats
        .join(VARPIPE_COVERAGEANALYSIS.out.target_region_coverage)
        .join(VARPIPE_PARSEANNOTATIONS_TARGETED.out.drLociFinalAnnotation)
        .set { ch_varpipe_report }

    // GENERATE ANALYSIS REPORT, INTERPRETATION & PDF REPORT

    VARPIPE_REPORT(
        ch_varpipe_report
    )
    ch_versions = ch_versions.mix(VARPIPE_REPORT.out.versions.first())

    VARPIPE_REPORT.out.summary
        .join(VARPIPE_STRUCTURALVARIANTS.out.structural_variants)
        .join(VARPIPE_CREATEANNOTATIONS_TARGETED.out.drLociAnnotation)
        .join(VARPIPE_COVERAGEANALYSIS.out.target_region_coverage)
        .set { ch_varpipe_interpretation }

    VARPIPE_INTERPRETATION(
        ch_varpipe_interpretation,
        reported
    )
    ch_versions = ch_versions.mix(VARPIPE_INTERPRETATION.out.versions.first())

    VARPIPE_PDF(
        VARPIPE_INTERPRETATION.out.summary
    )
    ch_versions = ch_versions.mix(VARPIPE_PDF.out.versions.first())

    ch_varpipe_outputs = Channel.empty()
        .mix(
            VARPIPE_COVERAGEANALYSIS.out.stats,
            VARPIPE_COVERAGEANALYSIS.out.qc_failed,
            VARPIPE_COVERAGEANALYSIS.out.target_region_coverage,
            VARPIPE_COVERAGEANALYSIS.out.genome_region_coverage,
            VARPIPE_STRUCTURALVARIANTS.out.structural_variants,
            VARPIPE_LINEAGE.out.lineage,
            VARPIPE_LINEAGE.out.lineage_report,
            VARPIPE_INTERPRETATION.out.summary,
            VARPIPE_INTERPRETATION.out.interpretation,
            VARPIPE_PDF.out.report,
            VARPIPE_PARSEANNOTATIONS_TARGETED.out.drLociFinalAnnotation,
            VARPIPE_CREATEANNOTATIONS_TARGETED.out.drLociAnnotation,
            VARPIPE_TARGETEDVARIANTS.out.ann,
            VARPIPE_ALLVARIANTS.out.ann,
            VARPIPE_PARSEANNOTATIONS_ALL.out.fullFinalAnnotation,
            VARPIPE_CREATEANNOTATIONS_ALL.out.fullAnnotation
        )
        .groupTuple()
        .multiMap { meta, files ->
            meta: [ id: "${meta.id}", filenames: files.collect { it.getName() } ]
            files: files
        }

    VARPIPE_TAR(
        ch_varpipe_outputs.meta.collect(),
        ch_varpipe_outputs.files.collect(),
        datestamp
    )
    ch_versions = ch_versions.mix(VARPIPE_TAR.out.versions)

    emit:
    lineage = VARPIPE_LINEAGE.out.lineage_report // channel (mandatory): [ val(meta), path(lineage_report) ]
    versions = ch_versions // channel (mandatory): [ versions.yml ]
}
