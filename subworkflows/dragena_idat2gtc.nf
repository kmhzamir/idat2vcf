process DRAGENA_IDAT2GTC {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/kmhzamir/dragena:v1.2"

    input:
    tuple val(meta), path(idats), path(txt)
    path manifest
    path clusterfile

    output:
    tuple val(meta), path("*.gtc"), path(txt), emit: gtc
    tuple val(meta), path("*.csv"), emit: csv
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml", emit: versions

    script:
    def idatFolder = "idat_folder_${meta.id}"

    """
    mkdir -p ${idatFolder}
    mv ${idats} ${idatFolder}/
    
    dragena genotype call \
    --bpm-manifest ${manifest} \
    --cluster-file ${clusterfile} \
    --idat-folder ${idatFolder} \
    --output-folder .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragena: v1.2
    END_VERSIONS
    """
}
