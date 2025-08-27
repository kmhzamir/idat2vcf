//
// A subworkflow to annotate snvs in the genome
//

include { CSVTK_CONCAT                          } from '../modules/nf-core/bcftools/concat_tab/main'
include { BCFTOOLS_NORM                         } from '../modules/nf-core/bcftools/norm/main'
include { RHOCALL_ANNOTATE                      } from '../modules/nf-core/rhocall/annotate/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_SNV      } from '../modules/nf-core/ensemblvep/vep/main'
include { TABIX_BGZIPTABIX as ZIP_TABIX_ROHCALL } from '../modules/tabix/bgziptabix/main'
include { TABIX_TABIX as TABIX_VEP              } from '../modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_BCFTOOLS_CONCAT  } from '../modules/tabix/tabix/main'
include { MERGE_VCF_COLUMNS                     } from '../modules/nf-core/mergevcfcolumns/main'
include { POSTPROCESS_VEP                       } from '../modules/nf-core/postprocess_vep/main'
include { ANNOTATE_CLINVAR                      } from '../modules/nf-core/annotate_clinvar/main'
include { ANNOTATE_GT                           } from '../modules/nf-core/annotate_gt/main'
include { ANNOTATE_MICROARRAY                   } from '../modules/nf-core/annotate_microarray/main'
include { ANNOTATE_PHENOTYPE                    } from '../modules/nf-core/annotate_phenotype/main'
include { GENEBE_ANNOTATE                       } from '../modules/nf-core/genebe/main'



workflow ANNOTATE_GENOME_SNVS {

    take:
        ch_vcf                // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ]
        val_vep_genome        // string: [mandatory] GRCh37 or GRCh38
        val_vep_cache_version // string: [mandatory] default: 107
        ch_vep_cache          // channel: [mandatory] [ path(cache) ]
        ch_genome_fasta       // channel: [mandatory] [ val(meta), path(fasta) ]
        //ch_gnomad_af          // channel: [optional] [ path(tab), path(tbi) ]
        ch_vep_extra_files    // channel: [mandatory] [ path(files) ]
        clinvar_cache
        phenotype_cache

    main:
        ch_cadd_vcf       = Channel.empty()
        ch_versions       = Channel.empty()
        ch_vcf_scatter_in = Channel.empty()
        ch_vep_in         = Channel.empty()

        BCFTOOLS_NORM (ch_vcf, ch_genome_fasta)

        ch_genome_vcf       = BCFTOOLS_NORM.out.vcf
        ch_genome_tabix     = BCFTOOLS_NORM.out.tbi
        ch_genome_vcf_tabix = ch_genome_vcf.join(ch_genome_tabix, failOnMismatch:true, failOnDuplicate:true)

        //GENEBE_ANNOTATE (ch_genome_vcf_tabix)

        // Annotating with ensembl Vep
        ENSEMBLVEP_SNV(
            ch_genome_vcf_tabix,
            val_vep_genome,
            "homo_sapiens",
            val_vep_cache_version,
            ch_vep_cache,
            ch_genome_fasta,
            ch_vep_extra_files
        )

        POSTPROCESS_VEP(ENSEMBLVEP_SNV.out.tab)

        MERGE_VCF_COLUMNS (POSTPROCESS_VEP.out.cleaned_output,ch_genome_vcf_tabix)

        ANNOTATE_CLINVAR(MERGE_VCF_COLUMNS.out.cleaned_output,clinvar_cache)

        ANNOTATE_PHENOTYPE(ANNOTATE_CLINVAR.out.annotated_clinvar,phenotype_cache)

        ANNOTATE_PHENOTYPE.out.annotated_clinvar_phenotype
            .map { meta, tab -> [meta - meta.subMap('scatterid'), tab] }
            .set { ch_vep_out }


        ANNOTATE_GT(ANNOTATE_PHENOTYPE.out.annotated_clinvar_phenotype)

        ANNOTATE_MICROARRAY(ch_genome_vcf_tabix, ANNOTATE_GT.out.gt_output)
        
        //TABIX_BCFTOOLS_CONCAT (ANNOTATE_GT.out.gt_output)

        ch_vep_ann   = ANNOTATE_GT.out.gt_output
        //ch_vep_index = TABIX_BCFTOOLS_CONCAT.out.tbi

        //ch_versions = ch_versions.mix(BCFTOOLS_ROH.out.versions)
        //ch_versions = ch_versions.mix(RHOCALL_ANNOTATE.out.versions)
        //ch_versions = ch_versions.mix(ZIP_TABIX_ROHCALL.out.versions)
        //ch_versions = ch_versions.mix(VCFANNO.out.versions)
        //ch_versions = ch_versions.mix(UPD_SITES.out.versions)
        //ch_versions = ch_versions.mix(UPD_REGIONS.out.versions)
        //ch_versions = ch_versions.mix(CHROMOGRAPH_SITES.out.versions)
        //ch_versions = ch_versions.mix(CHROMOGRAPH_REGIONS.out.versions)
        //ch_versions = ch_versions.mix(ZIP_TABIX_VCFANNO.out.versions)
        //ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)
        //ch_versions = ch_versions.mix(TABIX_BCFTOOLS_VIEW.out.versions)
        //ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions.first())
        ch_versions = ch_versions.mix(ENSEMBLVEP_SNV.out.versions.first())
        //ch_versions = ch_versions.mix(TABIX_VEP.out.versions.first())
        //ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions)
        //ch_versions = ch_versions.mix(TABIX_BCFTOOLS_CONCAT.out.versions)
        //ch_versions = ch_versions.mix(ANNOTATE_RHOCALLVIZ.out.versions)

    emit:
        vcf_ann  = ch_vep_ann   // channel: [ val(meta), path(vcf) ]
        //tbi      = ch_vep_index // channel: [ val(meta), path(tbi) ]
        versions = ch_versions  // channel: [ path(versions.yml) ]
}
