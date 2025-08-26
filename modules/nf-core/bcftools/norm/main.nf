process BCFTOOLS_NORM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.18--h8b25389_0':
        'biocontainers/bcftools:1.18--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(txt), path(tbi)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.vcf.gz"), path(txt), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--output-type z'
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Step 1: Normalize the VCF with bcftools
    bcftools norm -m -both \\
        -f ${fasta} \\
        -o raw_${prefix}.vcf \\
        --threads ${task.cpus} \\
        ${vcf}

    # Step 2: Clean RSIDs
    awk 'BEGIN{FS=OFS="\\t"} {
        if (\$3 ~ /^rs[0-9]/) {
            n = split(\$3, ids, ",");
            base_rs = ids[1];
            sub(/\\.[0-9]+\$/, "", base_rs);
            \$3 = base_rs;
        }
        print
        }' raw_${prefix}.vcf > ${prefix}_cleaned.vcf

    # Step 3: Compress and index
    bgzip -f ${prefix}_cleaned.vcf
    tabix -f -p vcf ${prefix}_cleaned.vcf.gz

    # Step 4: Version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
