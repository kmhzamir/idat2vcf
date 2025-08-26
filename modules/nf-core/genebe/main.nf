process GENEBE_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/genebe/pygenebe:0.1.15"
    
    input:
    tuple val(meta), path(vcf), path(tbi)                      

    output:
    tuple val(meta), path("${meta.id}_annotated_genebe.vcf")   , emit: genebe

    script:
    """
    genebe annotate \\
    --input "${vcf}" \\
    --genome hg19 \\
    --output "${meta.id}_annotated_genebe.vcf" \\
    --username zamir@synapselaboratory.com \\
    --api_key ak-AtrExCAXp4Zp1UgcQsu3uk0jI
    """
}
