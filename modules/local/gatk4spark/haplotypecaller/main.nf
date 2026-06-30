process GATK4SPARK_HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.6.0.0--py310hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.6.0.0--py310hdfd78af_0' }"

    input:
    tuple val(meta), path(input), path(input_index), path(intervals), path(dragstr_model)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    tuple val(meta5), path(dbsnp)
    tuple val(meta6), path(dbsnp_tbi)

    output:
    tuple val(meta), path("*.vcf.gz"),                emit: vcf
    tuple val(meta), path("*.tbi"),                   emit: tbi
    tuple val(meta), path("*.vcf.gz"), path("*.tbi"), emit: vcf_tbi
    path  "versions.yml",                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def intvls_a  = intervals     ? "--intervals ${intervals}" : ''
    def dragstr_a = dragstr_model ? "--dragstr-params-path ${dragstr_model}" : ''
    def dbsnp_a   = dbsnp         ? "--dbsnp ${dbsnp}" : ''
    def avail_mem = task.memory ? "-Xmx${task.memory.giga}g" : ''
    """
    gatk --java-options "${avail_mem}" HaplotypeCallerSpark \\
        --input ${input} \\
        --output ${prefix}.g.vcf.gz \\
        --reference ${fasta} \\
        --emit-ref-confidence GVCF \\
        --spark-master local[${task.cpus}] \\
        ${intvls_a} \\
        ${dragstr_a} \\
        ${dbsnp_a} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed 's/^.*GATK v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.g.vcf.gz
    touch ${prefix}.g.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.6.0.0
    END_VERSIONS
    """
}
