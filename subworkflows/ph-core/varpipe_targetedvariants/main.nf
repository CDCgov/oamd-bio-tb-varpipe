/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GATK4_MUTECT2                       } from "../../../modules/ph-core/gatk4/mutect2/main"
include { GATK4_LEFTALIGNANDTRIMVARIANTS      } from "../../../modules/ph-core/gatk4/leftalignandtrimvariants/main"
include { GATK4_FILTERMUTECTCALLS             } from "../../../modules/ph-core/gatk4/filtermutectcalls/main"
include { VARPIPE_SNPEFF                       } from "../../../modules/ph-core/varpipe/snpeff/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARPIPE_TARGETEDVARIANTS {
    take:
    ch_input        // channel (mandatory): [ val(meta), path(bam), path(bai), [] ]
    ch_fai          // channel (mandatory): [ val(meta), path ("*.fai") ]
    ch_ref_dict     // channel (mandatory): [ val(meta), path(reference_dict) ]
    ch_ref_fasta    // channel (mandatory): [ val(meta),path(reference_dict) ]
    snpeff_db       // string (mandatory): [ path(snpeff_db) ]

    main:
    ch_versions = Channel.empty()

    GATK4_MUTECT2 (
        ch_input,
        ch_ref_fasta,
        ch_fai,
        ch_ref_dict,
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions)

    ch_leftalignandtrimvariants = GATK4_MUTECT2.out.vcf
        .join(GATK4_MUTECT2.out.tbi)
        .map { meta, vcf, tbi ->
            return [ meta, vcf, tbi, [] ]
        }

    GATK4_LEFTALIGNANDTRIMVARIANTS (
        ch_leftalignandtrimvariants,
        ch_ref_fasta.map{ meta, fasta -> [ fasta ] },
        ch_fai.map{ meta, fai -> [ fai ] },
        ch_ref_dict.map{ meta, dict -> [ dict ] }
    )
    ch_versions = ch_versions.mix(GATK4_LEFTALIGNANDTRIMVARIANTS.out.versions)

    ch_filtermutectcalls = GATK4_LEFTALIGNANDTRIMVARIANTS.out.vcf
        .join(GATK4_LEFTALIGNANDTRIMVARIANTS.out.tbi)
        .join(GATK4_MUTECT2.out.stats)
        .map { meta, vcf, tbi, stats ->
            return tuple(
                meta,
                vcf,
                tbi,
                stats,
                [],
                [],
                [],
                []
            )
        }

    GATK4_FILTERMUTECTCALLS (
        ch_filtermutectcalls,
        ch_ref_fasta,
        ch_fai,
        ch_ref_dict
    )
    ch_versions = ch_versions.mix(GATK4_FILTERMUTECTCALLS.out.versions)

    VARPIPE_SNPEFF (
        GATK4_FILTERMUTECTCALLS.out.vcf,
        snpeff_db,
        [[], []]
    )
    ch_versions = ch_versions.mix(VARPIPE_SNPEFF.out.versions)

    emit:
    ann      =  VARPIPE_SNPEFF.out.vcf // channel (mandatory): [ val(meta), path(vcf) ]
    versions = ch_versions // channel (mandatory): [ versions.yml ]
}
