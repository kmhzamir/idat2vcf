process DRAGENA_GTC2VCF {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/kmhzamir/dragena:v1.2"

    input:
    tuple val(meta), path(gtc), path(txt)  // meta as metadata, gtc as a file path
    path manifest               // path for the manifest file
    path manifest_csv           // path for the manifest CSV file
    tuple val(meta2), path (genome)                 // path for the genome file
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("${meta.id}*.vcf.gz"), path(txt), emit: vcf
    tuple val(meta), path("${meta.id}*.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions

    script:
    def gtcFolder = "gtc_folder_${meta.id}"

    """
    mkdir -p ${gtcFolder}
    mv ${gtc} ${gtcFolder}/
    
    dragena genotype gtc-to-vcf \
    --bpm-manifest ${manifest} \
    --csv-manifest ${manifest_csv} \
    --genome-fasta-file ${genome} \
    --gtc-folder ${gtcFolder} \
    --output-folder .

    mv *.vcf.gz ${meta.id}.vcf.gz
    mv *.vcf.gz.tbi ${meta.id}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragena: v1.2
    END_VERSIONS
    """
}
