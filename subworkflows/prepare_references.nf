//
// Prepare reference files
//

include { TABIX_TABIX as TABIX_GNOMAD_AF                     } from '../modules/tabix/tabix/main'
include { UNTAR as UNTAR_VEP_CACHE                           } from '../modules/untar/main'

workflow PREPARE_REFERENCES {
    take:
        //ch_gnomad_af_tab   // channel: [optional; used in for snv annotation] [ val(meta), path(tab) ]
        ch_vep_cache       // channel: [mandatory for annotation] [ path(cache) ]



    main:
        ch_versions    = Channel.empty()

        //TABIX_GNOMAD_AF(ch_gnomad_af_tab)
        UNTAR_VEP_CACHE (ch_vep_cache)

        // Gather versions
       //ch_versions = ch_versions.mix(TABIX_GNOMAD_AF.out.versions)
        ch_versions = ch_versions.mix(UNTAR_VEP_CACHE.out.versions)


    emit:
        //gnomad_af_idx         = TABIX_GNOMAD_AF.out.tbi.collect()                                // channel: [ val(meta), path(fasta) ]
        vep_resources         = UNTAR_VEP_CACHE.out.untar.map{meta, files -> [files]}.collect()  // channel: [ path(cache) ]
        versions              = ch_versions                                                      // channel: [ path(versions.yml) ]

}
