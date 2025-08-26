/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowRaredisease.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHECK MANDATORY PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def mandatoryParams = [
    "manifest_bpm",
    "manifest_csv",
    "clusterfile",
    "fasta",
    "input",

]
def missingParamsCount = 0

for (param in mandatoryParams.unique()) {
    if (params[param] == null) {
        println("params." + param + " not set.")
        missingParamsCount += 1
    }
}

if (missingParamsCount>0) {
    error("\nSet missing parameters and restart the run. For more information please check usage documentation on github.")
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES AND SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// SUBWORKFLOWS

include { DRAGENA_IDAT2GTC                                   } from '../subworkflows/dragena_idat2gtc'
include { DRAGENA_GTC2VCF                                    } from '../subworkflows/dragena_gtc2vcf'
include { ANNOTATE_GENOME_SNVS                               } from '../subworkflows/annotate_genome_snvs'
include { PREPARE_REFERENCES                                 } from '../subworkflows/prepare_references'
include { CREATE_HGNCIDS_FILE                                } from '../modules/create_hgncids_file'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow IDAT2VCF {

    ch_versions = Channel.empty()

    // Initialize read, sample, and case_info channels
    ch_input = Channel.fromPath(params.input)
    Channel.fromSamplesheet("input")
        .tap { ch_original_input }
        .map { meta, idat1, idat2, txt -> 
            [meta, idat1, idat2, txt]
        }
        .reduce([:]) { counts, entry -> // Count each sample
            def (meta, idat1, idat2, txt) = entry
            counts[meta.sample] = (counts[meta.sample] ?: 0) + 1
            counts
        }
        .combine(ch_original_input)
        .map { counts, meta, idat1, idat2, txt ->
            def new_meta = meta
            return [new_meta, [idat1, idat2], txt]
        }
        .tap { ch_input_counts }
        .map { meta, idats, txt -> [idats, txt] }
        .reduce([:]) { counts, entry -> // Create line numbers
            def (idats, txt) = entry
            counts[[idats, txt]] = counts.size() + 1
            return counts
        }
        .combine(ch_input_counts)
        .map { lineno, meta, idats, txt ->
            def new_meta = meta + [id: meta.sample]
            return [new_meta, idats, txt]
        }
        .set { ch_reads }




    // Initialize file channels for PREPARE_REFERENCES subworkflow
    ch_genome_fasta             = Channel.fromPath(params.fasta).map { it -> [[id:it[0].simpleName], it] }.collect()
    ch_genome_fai               = params.fai                        ? Channel.fromPath(params.fai).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                    : Channel.empty()                                                  
    manifest_bpm                = params.manifest_bpm               ? Channel.fromPath(params.manifest_bpm).collect()
                                                                    : Channel.value([])
    manifest_csv                = params.manifest_csv               ? Channel.fromPath(params.manifest_csv).collect()
                                                                    : Channel.value([]) 
    clusterfile                 = params.clusterfile                ? Channel.fromPath(params.clusterfile).collect()
                                                                    : Channel.value([])
    ch_gnomad_af_tab            = params.gnomad_af                  ? Channel.fromPath(params.gnomad_af).map{ it -> [[id:it[0].simpleName], it] }.collect()
                                                                    : Channel.value([[],[]])
    ch_vep_cache_unprocessed    = params.vep_cache                  ? Channel.fromPath(params.vep_cache).map { it -> [[id:'vep_cache'], it] }.collect()
                                                                    : Channel.value([[],[]])
    ch_vep_filters_std_fmt      = params.vep_filters                ? Channel.fromPath(params.vep_filters).map { it -> [[id:'standard'],it]}.collect()
                                                                    : Channel.empty()
    ch_vep_filters              = params.vep_filters                ? Channel.fromPath(params.vep_filters).collect()
                                                                    : Channel.value([])                                                                                                                                                                                

    // Prepare references and indices.
    PREPARE_REFERENCES (
        ch_gnomad_af_tab,
        ch_vep_cache_unprocessed
    )
    .set { ch_references }                                                                                                                   

    // Gather built indices or get them from the params
    phenotype_cache             = params.phenotype_cache            ? Channel.fromPath(params.phenotype_cache).collect()
                                                                    : Channel.value([])
    clinvar_cache               = params.clinvar_cache              ? Channel.fromPath(params.clinvar_cache).collect()
                                                                    : Channel.value([])
    ch_vep_cache                = ( params.vep_cache && params.vep_cache.endsWith("tar.gz") )  ? ch_references.vep_resources
                                                                    : ( params.vep_cache    ? Channel.fromPath(params.vep_cache).collect() : Channel.value([]) )
    ch_gnomad_afidx             = params.gnomad_af_idx              ? Channel.fromPath(params.gnomad_af_idx).collect()
                                                                    : ch_references.gnomad_af_idx
    ch_gnomad_af                = params.gnomad_af                  ? ch_gnomad_af_tab.join(ch_gnomad_afidx).map {meta, tab, idx -> [tab,idx]}.collect()
                                                                    : Channel.empty()

    ch_vep_extra_files_unsplit  = params.vep_plugin_files           ? Channel.fromPath(params.vep_plugin_files).collect()
                                                                    : Channel.value([])
    ch_vep_filters_scout_fmt    = params.vep_filters_scout_fmt      ? Channel.fromPath(params.vep_filters_scout_fmt).map { it -> [[id:'scout'],it]}.collect()
                                                                    : Channel.empty()                                                                     

    // Read and store paths in the vep_plugin_files file
    if (params.vep_plugin_files) {
        ch_vep_extra_files_unsplit.splitCsv ( header:true )
            .map { row ->
                f = file(row.vep_files[0])
                if(f.isFile() || f.isDirectory()){
                    return [f]
                } else {
                    error("\nVep database file ${f} does not exist.")
                }
            }
            .collect()
            .set {ch_vep_extra_files}
    }

    // Read and store hgnc ids in a channel
    ch_vep_filters_scout_fmt
        .mix (ch_vep_filters_std_fmt)
        .set {ch_vep_filters}

    CREATE_HGNCIDS_FILE(ch_vep_filters)
        .txt
        .set {ch_hgnc_ids}


    //
    // DRAGENA_IDAT2GTC.
    //
    DRAGENA_IDAT2GTC (
        ch_reads,
        manifest_bpm,
        clusterfile,
    )
    ch_versions   = ch_versions.mix(DRAGENA_IDAT2GTC.out.versions)

    DRAGENA_GTC2VCF (
        DRAGENA_IDAT2GTC.out.gtc,
        manifest_bpm,
        manifest_csv,
        ch_genome_fasta,
        ch_genome_fai
    )
    ch_versions   = ch_versions.mix(DRAGENA_GTC2VCF.out.versions)

    ch_genome_vcf       = DRAGENA_GTC2VCF.out.vcf
    ch_genome_tabix     = DRAGENA_GTC2VCF.out.tbi
    ch_genome_vcf_tabix = ch_genome_vcf.join(ch_genome_tabix, failOnMismatch:true, failOnDuplicate:true)

    //
    // ANNOTATE GENOME SNVs
    //
    if (!params.skip_snv_annotation) {
        ANNOTATE_GENOME_SNVS (
            ch_genome_vcf_tabix,
            params.genome,
            params.vep_cache_version,
            ch_vep_cache,
            ch_genome_fasta,
            ch_gnomad_af,
            ch_vep_extra_files,
            clinvar_cache,
            phenotype_cache
        ).set { ch_snv_annotate }
        ch_versions = ch_versions.mix(ch_snv_annotate.versions)
    }
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/